! *********************************************************************
! * Rocstar Simulation Suite                                          *
! * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
! *                                                                   *
! * Illinois Rocstar LLC                                              *
! * Champaign, IL                                                     *
! * www.illinoisrocstar.com                                           *
! * sales@illinoisrocstar.com                                         *
! *                                                                   *
! * License: See LICENSE file in top level of distribution package or *
! * http://opensource.org/licenses/NCSA                               *
! *********************************************************************
! *********************************************************************
! * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
! * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
! * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
! * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
! * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
! * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
! * Arising FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
! * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
! *********************************************************************
!******************************************************************************
!
! Purpose: Initialisation of TURB solution
!
! Description: Besides initialisation of TURB variables, in case non-uniform
!              filter is selected, non-uniform filter coefficients and
!              averaging coefficients are computed by calling the corresponding
!              routines. If wlm is active, define wlm mapping coefficients and
!              copy or interpolate wlm patch vals from level 1 to other levels.
!
! Input: region = data of current region
!
! Output: mueT, tCoT and cDyn get initial values. Wlm mapping coefficients
!         initiated and wlm patch vals interpolated to other levels.
!
! Notes: cDyn is set to zero every stage in ViscousFluxes.
!        Non-uniform quantities are recomputed in ViscousFluxes everytime
!        the grid moves.
!
!******************************************************************************
!
! $Id: TURB_InitSolution.F90,v 1.20 2009/08/26 12:28:52 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_InitSolution( region ) ! PUBLIC

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGrid, ONLY       : t_grid
  USE ModTurbulence, ONLY : t_turb
  USE ModGlobal, ONLY     : t_global
  USE TURB_ModInterfaces, ONLY : TURB_RansSAGetEddyVis, TURB_WlmInitia
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset
  USE TURB_ModInterfaces, ONLY : TURB_FloFaceVolume, &
                                 TURB_FloFaceWidth,  TURB_FloLesGenCoCC, &
                                 TURB_FloLesGenCoFF, TURB_FloWlmMetric 
#include "Indexing.h"
#endif
#ifdef RFLU
  USE ModInterfaces, ONLY : MixtureProperties
  USE TURB_ModInterfaces, ONLY : TURB_FluFaceVolume
#endif
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
#ifdef RFLO
  TYPE(t_region) :: region
#endif
#ifdef RFLU
  TYPE(t_region), TARGET :: region
#endif

! ... loop variables
  INTEGER         :: iPatch, iN, l, ifl

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_grid) , POINTER  :: grid
  TYPE(t_turb), POINTER   :: turb
  TYPE(t_patch), POINTER  :: patch1, patch

  LOGICAL               :: doWlm
  INTEGER               :: turbModel
  REAL(RFREAL)          :: tripLoc(XCOORD:ZCOORD), treshold
  REAL(RFREAL), POINTER :: tv(:,:), tcv(:,:), tcvOld(:,:), tdv(:,:), vort(:,:)
  REAL(RFREAL), POINTER :: dsterm(:,:)
#ifdef RFLO
  INTEGER :: iLev, iNOff, ijNOff, ibn, ien
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend, errorFlag
#endif

  CHARACTER(2*CHRLEN+17) :: fname
  LOGICAL                :: fileExists

#ifdef RFLU
  INTEGER :: ibc, iec, ibn, ien, ic, ifg, ifgBeg
#endif
!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_InitSolution.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_InitSolution',&
  'TURB_InitSolution.F90' )

! pre procedures -----------------------------------------------------------

#ifdef RFLU
  ibc = 1
  iec = region%grid%nCellsTot
  CALL MixtureProperties( region,ibc,iec,.true. )
#endif

! get dimensions -----------------------------------------------------------

#ifdef RFLO
  iLev =  region%currLevel
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)
#endif
#ifdef RFLU
  ibn = 1
  ien = region%grid%nFaces
#endif

! get parameters and pointers ----------------------------------------------

  turbModel =  region%mixtInput%turbModel

#ifdef RFLO
  grid => region%levels(iLev)%grid
  turb => region%levels(iLev)%turb
  tv   => region%levels(iLev)%mixt%tv

  IF (region%turbInput%nDv > 0) THEN
    tdv => turb%dv
  ENDIF
  IF (region%turbInput%nCv > 0) THEN
    tcv    => turb%cv
    tcvOld => turb%cvOld
    dsterm => turb%dsterm
  ENDIF
#endif
#ifdef RFLU
  grid => region%grid
  turb => region%turb
  tv   => region%mixt%tv

  IF (region%turbInput%nDv > 0) THEN
    tdv => region%turb%dv
  ENDIF
  IF (region%turbInput%nCv > 0) THEN
    tcv    => region%turb%cv
    tcvOld => region%turb%cvOld
    dsterm => region%turb%dsterm
  ENDIF
#endif

! initialize turbulent viscosity and thermal conductivity at cells (tv of
! NS system), turbulence variables (turb.cvOld, dv, vorticities, etc),
! and metrics pertinent to turbulence

! general ---------------------------------------------------------
! first check if it's a restart from a laminar run
      IF (global%solutFormat == FORMAT_ASCII .OR. global%solutFormat == FORMAT_HDF) THEN
        WRITE(fname,'(A,1PE11.5)') TRIM(global%inDir)//TRIM(global%casename)//'.turba_', &
                                   global%timeStamp
        INQUIRE(FILE=fname,EXIST=fileExists)
      ELSE IF (global%solutFormat == FORMAT_BINARY) THEN
        WRITE(fname,'(A,I6.6)') TRIM(global%inDir)//TRIM(global%casename)//'.turbb_', &
                                global%currentIter
        INQUIRE(FILE=fname,EXIST=fileExists)
      ELSE
        CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
      ENDIF

  IF ((global%flowType == FLOW_UNSTEADY .AND. &
       global%timeStamp <= 0._RFREAL)    .OR. &
      (global%flowType == FLOW_STEADY   .AND. &
       global%currentIter <= 0)         .OR.  &
       (global%flowType == FLOW_UNSTEADY .AND. &
       global%timeStamp > 0._RFREAL   .AND. &
       fileExists .EQV. .false.)    .OR. &
       (global%flowType == FLOW_STEADY   .AND. &
       global%currentIter > 0     .AND. &
       fileExists .EQV. .false.)) THEN
    tv(TV_MIXT_MUET,:)  = 0._RFREAL  
    global%esg1Sum      = 0._RFREAL                ! pertinent to
    global%esg4Sum      = 0._RFREAL                ! LES energy model 
    IF (region%turbInput%nCv > 0) THEN
#ifdef RFLO
      tcv(CV_SA_NUTIL,:)= tv(TV_MIXT_MUEL,:)       ! pertinent to RaNS 
#endif
#ifdef RFLU
      tcv(CV_SA_NUTIL,:)= region%mixtInput%refVisc ! mixtProp, dv unknown yet
#endif
      tcvOld(:,:)       = 0._RFREAL          
    ENDIF
  ENDIF
  tv(TV_MIXT_TCOT,:)  = 0._RFREAL

  IF (region%turbInput%nCv > 0) THEN
    dsterm(:,:) = 0._RFREAL
  ENDIF

  IF (region%turbInput%nDv > 0) THEN
    tdv(:,:) = 0._RFREAL
  ENDIF

  IF (ASSOCIATED( turb%vort )) THEN
    turb%vort = 0._RFREAL
  ENDIF

! model dependent ---------------------------------------

  IF ((region%mixtInput%turbModel == TURB_MODEL_SA) .OR. &
      (region%mixtInput%turbModel == TURB_MODEL_DESSA) .OR. &
      (region%mixtInput%turbModel == TURB_MODEL_HDESSA)) THEN
    CALL TURB_RansSAGetEddyVis( region ) ! nutilde to muet
  ENDIF

! metrics related ----------------------------------------

! compute face volumes needed for lesMij
  IF ((turbModel==TURB_MODEL_FIXSMAG) .OR. &
      (turbModel==TURB_MODEL_DYNSMAG) .OR. &
      (turbModel==TURB_MODEL_DYNMIXD)) THEN
#ifdef RFLO
    CALL TURB_FloFaceVolume( region,DIRI )
    CALL TURB_FloFaceVolume( region,DIRJ )
    CALL TURB_FloFaceVolume( region,DIRK )
#endif
#ifdef RFLU
    CALL TURB_FluFaceVolume( region )
#endif
  ENDIF

#ifdef RFLO
! compute filter coefficients

  IF (((turbModel==TURB_MODEL_SCALSIM)  .OR. &
       (turbModel==TURB_MODEL_DYNSMAG)  .OR. &
       (turbModel==TURB_MODEL_DYNMIXD)) .AND. &
       (region%turbInput%filterType == FILTYPE_NONUNIF)) THEN

    ALLOCATE( turb%workI(2,ibn:ien),stat=errorFlag )
    ALLOCATE( turb%workJ(2,ibn:ien),stat=errorFlag )
    ALLOCATE( turb%workK(2,ibn:ien),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    CALL TURB_FloFaceWidth( region )
    CALL TURB_FloLesGenCoCC( region )
    CALL TURB_FloLesGenCoFF( region )
    DEALLOCATE( turb%workI,turb%workJ,turb%workK )
  ENDIF
#endif

! initialize zero-one fields

  IF (region%turbInput%nZof > 0) THEN
#ifdef RFLO
    turb%zofi = 1._RFREAL
    turb%zofj = 1._RFREAL
    turb%zofk = 1._RFREAL
#endif
#ifdef RFLU
    turb%zofi  = 1._RFREAL
    turb%bZofi = 1._RFREAL
#endif
  ENDIF

! create zero-one fields

  IF (turbModel==TURB_MODEL_FIXSMAG .OR. &
      turbModel==TURB_MODEL_DYNSMAG .OR. &
      turbModel==TURB_MODEL_DYNMIXD) THEN
    tripLoc(:) =  region%turbInput%xyzSmag(:)
    treshold   = -HUGE( 1.0_RFREAL )/1000._RFREAL

    IF ((tripLoc(XCOORD) > treshold) .OR. &
        (tripLoc(YCOORD) > treshold) .OR. &
        (tripLoc(ZCOORD) > treshold)) THEN

#ifdef RFLO
      DO iN = ibn,ien
        DO l = XCOORD,ZCOORD
          IF (grid%cfcI(l,iN) < tripLoc(l)) THEN
            turb%zofi(l,ZOF_LES_EDDYVIS,iN) = 0._RFREAL
          ELSE
            turb%zofi(l,ZOF_LES_EDDYVIS,iN) = 1._RFREAL
          ENDIF
        ENDDO ! l
        DO l = XCOORD,ZCOORD
          IF (grid%cfcJ(l,iN) < tripLoc(l)) THEN
            turb%zofj(l,ZOF_LES_EDDYVIS,iN) = 0._RFREAL
          ELSE
            turb%zofj(l,ZOF_LES_EDDYVIS,iN) = 1._RFREAL
          ENDIF
        ENDDO ! l
        DO l = XCOORD,ZCOORD
          IF (grid%cfcK(l,iN) < tripLoc(l)) THEN
            turb%zofk(l,ZOF_LES_EDDYVIS,iN) = 0._RFREAL
          ELSE
            turb%zofk(l,ZOF_LES_EDDYVIS,iN) = 1._RFREAL
          ENDIF
        ENDDO ! l
      ENDDO   ! iN
#endif
#ifdef RFLU
      DO iN = ibn,ien
        DO l = XCOORD,ZCOORD
          IF (grid%fc(l,iN) < tripLoc(l)) THEN
            turb%zofi(l,ZOF_LES_EDDYVIS,iN) = 0._RFREAL
          ELSE
            turb%zofi(l,ZOF_LES_EDDYVIS,iN) = 1._RFREAL
          ENDIF
        ENDDO ! l
      ENDDO   ! iN
      DO iPatch = 1,grid%nPatches
        patch  => region%patches(iPatch)
! TEMPORARY : removing usage of bf2bg from everywhere
!        ifgBeg =  patch%bf2bg(BF2BG_BEG)

        DO ifl = 1,patch%nBFaces
          ic  = patch%bf2c(ifl)
          ifg = ifl + ifgBeg-1

          DO l = XCOORD,ZCOORD
            IF (grid%cofg(l,ic) < tripLoc(l)) THEN
              turb%bZofi(l,ZOF_LES_EDDYVIS,ifg) = 0._RFREAL
            ELSE
              turb%bZofi(l,ZOF_LES_EDDYVIS,ifg) = 1._RFREAL
            ENDIF
          ENDDO ! l
        ENDDO   ! ifl
      ENDDO     ! iPatch
#endif

    ENDIF ! triploc
  ENDIF   ! turbModel

! if applicable, initiate metric variables and utau of wlm

#ifdef RFLO
  DO iPatch=1,region%nPatches
    patch1 => region%levels(1)%patches(iPatch)
    patch  => region%levels(iLev)%patches(iPatch)
#endif
#ifdef RFLU
  DO iPatch=1,region%grid%nPatches
    patch  => region%patches(iPatch)
#endif

    doWlm = .false.
#ifdef RFLO
    IF (patch%bcType>=BC_NOSLIPWALL .AND. &
        patch%bcType<=BC_NOSLIPWALL+BC_RANGE) THEN ! my boundary type
      IF (patch%valBola%switches(WLM_INPUT_MODEL) /= WLM_MODEL_NOMODEL) THEN
        doWlm = .true.
      ENDIF
    ENDIF
#endif

    IF (doWlm .eqv. .true.) THEN

! --- get initial estimate of friction velocity utau
      CALL TURB_WlmInitia( region,patch )

#ifdef RFLO
! --- compute mapping coefficients from body fitted to cartesian and other metrics
      CALL TURB_FloWlmMetric( region,patch )

      IF (patch%valBola%distrib==BCDAT_DISTRIB) THEN
! ----- roughness distribution
!       CALL TURB_interpolate2Levels( )
        CALL ErrorStop( global,ERR_VAL_BCVAL,__LINE__,'No variable roughness yet'  )
      ELSE
! ----- distribution has constant value, so copied from a point on patch level 1
        patch%valBola%vals(:,WLM_VALS_ROUGH)= patch1%valBola%vals(1,WLM_VALS_ROUGH)
      ENDIF  ! distribution

      IF (patch%valBola%switches(WLM_INPUT_MODEL) == WLM_MODEL_EXTERN) THEN
! ----- wall stress distribution
!       CALL TURB_interpolate2Levels( )
      ENDIF
#endif

    ENDIF    ! doWlm
  ENDDO      ! iPatch

! finalize --------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_InitSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_InitSolution.F90,v $
! Revision 1.20  2009/08/26 12:28:52  mtcampbe
! Ported to Hera.   Fixed logical expression syntax errors.  Replaced all
! IF (logical_variable)  with IF (logical_variable .eqv. .true.) as
! consistent with the specification.  Also changed: IF( ASSOCIATED(expr) )
! to IF ( ASSOCIATED(expr) .eqv. .true. ).   Intel compilers produce code
! which silently fails for some mal-formed expressions, so these changes
! are a net which should ensure that they are evaluated as intended.
!
! Revision 1.19  2009/08/12 04:15:59  mtcampbe
! Major update, bugfix from Abe development, more propagation compatibility,
! some Rocstar IO changes, Ju's temporary clipping fix for turbulence. A bug
! fix for initialization IO.
!
! Revision 1.18  2009/06/29 17:15:32  juzhang
! initialization of RANS/DES added for restart from a laminar run
!
! Revision 1.17  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.16  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.15  2006/08/19 15:40:58  mparmar
! Removed bf2bg
!
! Revision 1.14  2006/01/30 23:06:15  wasistho
! removed comment after ifdef RFLU
!
! Revision 1.13  2006/01/17 17:51:40  wasistho
! applied tripping to all eddy viscosity models
!
! Revision 1.12  2006/01/12 23:46:55  wasistho
! POINTER to TARGET in rflu
!
! Revision 1.11  2006/01/12 09:48:49  wasistho
! enabled tripping fixed Smagorinsky
!
! Revision 1.10  2005/12/30 23:20:32  wasistho
! exclude rocflu from WLM treatment
!
! Revision 1.9  2005/12/29 19:48:01  wasistho
! modified rflu part
!
! Revision 1.8  2005/03/09 06:35:01  wasistho
! incorporated HDESSA
!
! Revision 1.7  2004/08/04 02:45:47  wasistho
! removed turb%avgCoI,J,K as it is defined as grid%c2fCoI,J,K
!
! Revision 1.6  2004/06/19 03:29:48  wasistho
! removed argument iReg in TURB_InitSolution
!
! Revision 1.5  2004/06/03 02:10:47  wasistho
! enabled non-uniform fix-Smagorinsky
!
! Revision 1.4  2004/03/20 03:28:29  wasistho
! prepared for RFLU
!
! Revision 1.3  2004/03/19 02:47:13  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.12  2004/02/26 21:26:29  wasistho
! initialize esg1Sum, esg4Sum
!
! Revision 1.11  2003/10/24 03:46:43  wasistho
! initiate mu_t to mu_l if RaNS is active
!
! Revision 1.10  2003/10/21 20:31:50  wasistho
! added dt relaxation in steady flow due to RANS source term
!
! Revision 1.9  2003/10/09 22:50:13  wasistho
! mv call to TURB_RansSAGetEddyVis from readSolution to initSolution
!
! Revision 1.8  2003/10/07 02:06:03  wasistho
! initial installation of RaNS-SA and DES
!
! Revision 1.7  2003/08/08 01:46:24  wasistho
! fixed turb. restart for steady flow
!
! Revision 1.6  2003/08/06 15:56:13  wasistho
! added vorticities computation
!
! Revision 1.5  2003/08/01 22:17:43  wasistho
! prepared rocturb for Genx
!
! Revision 1.4  2003/07/23 15:59:40  wasistho
! prepared more accurate rocturb restart
!
! Revision 1.3  2003/05/31 01:46:52  wasistho
! installed turb. wall layer model
!
! Revision 1.2  2002/10/16 07:48:19  wasistho
! Enable Fix Smagorinsky
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************







