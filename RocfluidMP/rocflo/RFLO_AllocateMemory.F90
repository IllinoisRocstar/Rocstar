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
! Purpose: allocate memory for all variables associated with the mixture
!          for all active regions on current processor.
!
! Description: none.
!
! Input: region = current region
!
! Output: region%mixt = mixture variables
!         region%grid = grid variables (also region%gridOld)
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_AllocateMemory.F90,v 1.24 2008/12/06 08:44:25 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_AllocateMemory( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModBndPatch, ONLY   : t_patch
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset, &
        RFLO_GetDimensDummy, RFLO_GetCellOffset, RFLO_GetDimensPhysNodes
  USE ModMixture, ONLY    : t_mixt
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), TARGET :: region

! ... loop variables
  INTEGER :: iLev, iPatch

! ... local variables
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: ibc, iec, ibn, ien, iCellOffset, ijCellOffset, &
             iNodeOffset, ijNodeOffset, h1, h2, errorFlag
  INTEGER :: n1, n2, ijBeg, ijEnd, iOff

  TYPE(t_level), POINTER  :: level
  TYPE(t_grid) , POINTER  :: grid, gridOld, gridOld2
  TYPE(t_mixt) , POINTER  :: mixt
  TYPE(t_global), POINTER :: global
  TYPE(t_patch) , POINTER :: patch

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_AllocateMemory',&
  'RFLO_AllocateMemory.F90' )

! loop over all grid levels

  DO iLev=1,region%nGridLevels

    level    => region%levels(iLev)
    grid     => region%levels(iLev)%grid
    gridOld  => region%levels(iLev)%gridOld
    gridOld2 => region%levels(iLev)%gridOld2
    mixt     => region%levels(iLev)%mixt

! - get cell and node dimensions

    CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( region,iLev,iCellOffset,ijCellOffset )
    ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCellOffset,ijCellOffset)
    iec = IndIJK(idcend,jdcend,kdcend,iCellOffset,ijCellOffset)

    CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetNodeOffset( region,iLev,iNodeOffset,ijNodeOffset )
    ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNodeOffset,ijNodeOffset)
    ien = IndIJK(idnend,jdnend,kdnend,iNodeOffset,ijNodeOffset)

    CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                  jpnbeg,jpnend,kpnbeg,kpnend )

! - grid coordinates

    ALLOCATE( grid%xyz(3,ibn:ien)  ,stat=errorFlag )
    ALLOCATE( grid%ijkDgen(ibn:ien),stat=errorFlag )
    IF (region%mixtInput%moveGrid) THEN
      ALLOCATE( gridOld%xyzOld(3,ibn:ien),stat=errorFlag )
      ALLOCATE( gridOld%xyz(3,ibn:ien)   ,stat=errorFlag )
      IF (iLev == 1) THEN     ! only on the finest grid
        ALLOCATE( grid%arcLen12(jpnbeg:jpnend,kpnbeg:kpnend),stat=errorFlag )
        ALLOCATE( grid%arcLen34(kpnbeg:kpnend,ipnbeg:ipnend),stat=errorFlag )
        ALLOCATE( grid%arcLen56(ipnbeg:ipnend,jpnbeg:jpnend),stat=errorFlag )
      ENDIF
      IF (global%moveGridScheme/=MOVEGRID_BLOCKS .AND. iLev==1) THEN
        ALLOCATE( grid%xyzOld(3,ibn:ien),stat=errorFlag )
      ENDIF
      IF ((global%moveGridScheme==MOVEGRID_FRAME .OR. &
           global%moveGridScheme==MOVEGRID_FOMS  .OR. &
           global%moveGridScheme==MOVEGRID_ELFRAME) .AND. iLev==1) THEN
        ALLOCATE( grid%nCorns(global%nRegions), stat=errorFlag )
      ENDIF
      IF (global%moveGridScheme==MOVEGRID_FOMS .AND. iLev==1) THEN
        ALLOCATE( grid%xyzOrth(3,ibn:ien), stat=errorFlag )
        ALLOCATE( grid%xyzTemp(3,ibn:ien), stat=errorFlag )
      ENDIF
      IF ((global%moveGridScheme==MOVEGRID_ELGLOBAL .OR. &
           global%moveGridScheme==MOVEGRID_ELFRAME) .AND. iLev==1) THEN
        ALLOCATE( grid%stu(      3,ibn:ien),stat=errorFlag )
        ALLOCATE( grid%stuOld(   3,ibn:ien),stat=errorFlag )
        ALLOCATE( grid%stui(     3,ibn:ien),stat=errorFlag )
        ALLOCATE( grid%stuj(     3,ibn:ien),stat=errorFlag )
        ALLOCATE( grid%stuk(     3,ibn:ien),stat=errorFlag )
        ALLOCATE( grid%pmat(   3,6,ibn:ien),stat=errorFlag )

        ALLOCATE( grid%aijk(       ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aimjk(      ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aipjk(      ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aijmk(      ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aijpk(      ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aijkm(      ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aijkp(      ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aimjmk(     ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aipjmk(     ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aimjpk(     ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aipjpk(     ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aimjkm(     ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aipjkm(     ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aimjkp(     ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aipjkp(     ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aijmkm(     ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aijpkm(     ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aijmkp(     ibn:ien),stat=errorFlag )
        ALLOCATE( grid%aijpkp(     ibn:ien),stat=errorFlag )

        DO iPatch=1,region%nPatches
          patch  => region%levels(iLev)%patches(iPatch)

          h1 = ABS(patch%l1end-patch%l1beg) + 2
          h2 = ABS(patch%l2end-patch%l2beg) + 2

! ------- elliptic PDE grid motion

          ALLOCATE( patch%st(     3,h1,h2),stat=errorFlag )
          ALLOCATE( patch%stOld(  3,h1,h2),stat=errorFlag )
          ALLOCATE( patch%sti(    3,h1,h2),stat=errorFlag )
          ALLOCATE( patch%stj(    3,h1,h2),stat=errorFlag )
          ALLOCATE( patch%stii(   3,h1,h2),stat=errorFlag )
          ALLOCATE( patch%stjj(   3,h1,h2),stat=errorFlag )
          ALLOCATE( patch%stij(   3,h1,h2),stat=errorFlag )
          ALLOCATE( patch%pfun( 2,3,h1,h2),stat=errorFlag )

          ALLOCATE( patch%aimjm( h1,h2),stat=errorFlag )
          ALLOCATE( patch%aijm(  h1,h2),stat=errorFlag )
          ALLOCATE( patch%aipjm( h1,h2),stat=errorFlag )
          ALLOCATE( patch%aimj(  h1,h2),stat=errorFlag )
          ALLOCATE( patch%aij(   h1,h2),stat=errorFlag )
          ALLOCATE( patch%aipj(  h1,h2),stat=errorFlag )
          ALLOCATE( patch%aimjp( h1,h2),stat=errorFlag )
          ALLOCATE( patch%aijp(  h1,h2),stat=errorFlag )
          ALLOCATE( patch%aipjp( h1,h2),stat=errorFlag )

        ENDDO      ! iPatch
      ENDIF        ! moveGridScheme
    ENDIF          ! moveGrid
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! - averaging coefficients, face vectors, volumes, grid speeds

    ALLOCATE( grid%c2fCoI(2,ibn:ien),stat=errorFlag )
    ALLOCATE( grid%c2fCoJ(2,ibn:ien),stat=errorFlag )
    ALLOCATE( grid%c2fCoK(2,ibn:ien),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    ALLOCATE( grid%c2eCoI(4,ibn:ien),stat=errorFlag )
    ALLOCATE( grid%c2eCoJ(4,ibn:ien),stat=errorFlag )
    ALLOCATE( grid%c2eCoK(4,ibn:ien),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    ALLOCATE( grid%si(4,ibn:ien),stat=errorFlag )
    ALLOCATE( grid%sj(4,ibn:ien),stat=errorFlag )
    ALLOCATE( grid%sk(4,ibn:ien),stat=errorFlag )
    ALLOCATE( grid%vol(ibc:iec) ,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    IF ((region%mixtInput%spaceOrder == DISCR_ORDER_4) .OR. &
        global%calcCellCtr) THEN
      ALLOCATE( grid%cofg(3,ibc:iec),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    ELSE
      NULLIFY( grid%cofg )
    ENDIF

    IF (global%calcFaceCtr) THEN
      ALLOCATE( grid%cfcI(3,ibn:ien),stat=errorFlag )
      ALLOCATE( grid%cfcJ(3,ibn:ien),stat=errorFlag )
      ALLOCATE( grid%cfcK(3,ibn:ien),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    ELSE
      NULLIFY( grid%cfcI, grid%cfcJ, grid%cfcK )
    ENDIF

    IF (region%mixtInput%moveGrid) THEN
      ALLOCATE( grid%siVel(ibn:ien)  ,stat=errorFlag )
      ALLOCATE( grid%sjVel(ibn:ien)  ,stat=errorFlag )
      ALLOCATE( grid%skVel(ibn:ien)  ,stat=errorFlag )
      ALLOCATE( gridOld%si(4,ibn:ien),stat=errorFlag )
      ALLOCATE( gridOld%sj(4,ibn:ien),stat=errorFlag )
      ALLOCATE( gridOld%sk(4,ibn:ien),stat=errorFlag )
      ALLOCATE( gridOld%vol(ibc:iec) ,stat=errorflag )
      IF (global%solverType==SOLV_IMPLICIT .AND. &
          global%tstepOrder==3 .AND. iLev == 1) THEN
        ALLOCATE( gridOld2%vol(ibc:iec),stat=errorflag )
      ENDIF
    ELSE
      ALLOCATE( grid%siVel(0:1),stat=errorFlag )
      ALLOCATE( grid%sjVel(0:1),stat=errorFlag )
      ALLOCATE( grid%skVel(0:1),stat=errorFlag )
    ENDIF
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    grid%siVel(:) = 0._RFREAL   ! zero out grid speeds
    grid%sjVel(:) = 0._RFREAL
    grid%skVel(:) = 0._RFREAL

! - time step & spectral radii

    ALLOCATE( level%dt(ibc:iec)   ,stat=errorFlag )
    ALLOCATE( mixt%srad(3,ibc:iec),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! - mixture variables

    ALLOCATE( mixt%cv   (CV_MIXT_NEQS,ibc:iec),stat=errorFlag )
    ALLOCATE( mixt%cvOld(CV_MIXT_NEQS,ibc:iec),stat=errorFlag )
    ALLOCATE( mixt%dv   (mixt%nDv    ,ibc:iec),stat=errorFlag )
    ALLOCATE( mixt%rhs  (CV_MIXT_NEQS,ibc:iec),stat=errorFlag )
    ALLOCATE( mixt%diss (CV_MIXT_NEQS,ibc:iec),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    IF (mixt%nTv > 0) THEN
      ALLOCATE( mixt%tv(mixt%nTv,ibc:iec),stat=errorFlag )
    ELSE
      NULLIFY( mixt%tv )
    ENDIF

    IF (mixt%nGrad > 0) THEN
      ALLOCATE( mixt%gradi(mixt%nGrad,ibn:ien),stat=errorFlag )
      ALLOCATE( mixt%gradj(mixt%nGrad,ibn:ien),stat=errorFlag )
      ALLOCATE( mixt%gradk(mixt%nGrad,ibn:ien),stat=errorFlag )
    ELSE
      NULLIFY( mixt%gradi )
      NULLIFY( mixt%gradj )
      NULLIFY( mixt%gradk )
    ENDIF

    IF (region%mixtInput%gasModel == GAS_MODEL_TCPERF) THEN
      ALLOCATE( mixt%gv(mixt%nGv,0:1),stat=errorFlag )
    ELSE
      ALLOCATE( mixt%gv(mixt%nGv,ibc:iec),stat=errorFlag )
    ENDIF

#ifdef STATS
    IF ((global%flowType == FLOW_UNSTEADY) .AND. &
        (global%doStat == ACTIVE) .AND. &
        (global%mixtNStat > 0)) THEN
      ALLOCATE( mixt%tav(global%mixtNStat,ibc:iec),stat=errorFlag )
    ELSE
      NULLIFY( mixt%tav )
    ENDIF
#endif

    IF (iLev>1 .AND. global%cycleType/=MGCYCLE_NO) THEN
      ALLOCATE( mixt%fterm(CV_MIXT_NEQS,ibc:iec),stat=errorFlag )
    ELSE
      NULLIFY( mixt%fterm )
    ENDIF

    IF (global%flowType == FLOW_UNSTEADY) THEN
      ALLOCATE( mixt%rhsSum(CV_MIXT_NEQS,ibc:iec),stat=errorFlag )
    ELSE
      NULLIFY( mixt%rhsSum )
    ENDIF

    IF (global%flowType==FLOW_UNSTEADY .AND. &
        global%solverType==SOLV_IMPLICIT .AND. &
        iLev == 1) THEN
      ALLOCATE( mixt%cvn  (CV_MIXT_NEQS,ibc:iec),stat=errorFlag )
      ALLOCATE( mixt%cvn1 (CV_MIXT_NEQS,ibc:iec),stat=errorFlag )
      ALLOCATE( mixt%cvn2 (CV_MIXT_NEQS,ibc:iec),stat=errorFlag )
      ALLOCATE( mixt%sDual(CV_MIXT_NEQS,ibc:iec),stat=errorFlag )
    ELSE
      NULLIFY( mixt%cvn   )
      NULLIFY( mixt%cvn1  )
      NULLIFY( mixt%cvn2  )
      NULLIFY( mixt%sDual )
    ENDIF

    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! - patch data allocation -----------------------------------------------------

    DO iPatch=1,region%nPatches
      patch  => region%levels(iLev)%patches(iPatch)

! --- aerodynamic coefficients

      IF (global%aeroCoeffs == ACTIVE) THEN
        n1    = ABS(patch%l1end-patch%l1beg)
        n2    = ABS(patch%l2end-patch%l2beg)
        iOff  = n1 + 1
        ijBeg = IndIJ( 0, 0,iOff)
        ijEnd = IndIJ(n1,n2,iOff)
      ELSE
        ijBeg = 0
        ijEnd = 1
      ENDIF
      ALLOCATE( patch%cp(ijBeg:ijEnd),   stat=errorFlag )
      ALLOCATE( patch%ch(ijBeg:ijEnd),   stat=errorFlag )
      ALLOCATE( patch%cf(3,ijBeg:ijEnd), stat=errorFlag )
      ALLOCATE( patch%forceCoeffs(3,FORCES_PRESS:FORCES_VISC), stat=errorFlag )
      ALLOCATE( patch%momentCoeffs(3,FORCES_PRESS:FORCES_VISC),stat=errorFlag )

      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! --- patch arclengths

      h1 = ABS(patch%l1end-patch%l1beg) + 2
      h2 = ABS(patch%l2end-patch%l2beg) + 2
      ALLOCATE( patch%arclen1( h2 ), stat=errorFlag )
      ALLOCATE( patch%arclen2( h1 ), stat=errorFlag )

      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    ENDDO ! iPatch

  ENDDO   ! iLev

! data allocation independent of grid level -----------------------------------

! coefficients of implicit residual smoothing

  IF ((global%flowType==FLOW_STEADY .OR. &
      (global%flowType==FLOW_UNSTEADY .AND. global%solverType==SOLV_IMPLICIT)) &
      .AND. region%mixtInput%smoocf>0.) THEN
    ALLOCATE( mixt%epsIrs(3,ibc:iec),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
  ELSE
    NULLIFY( mixt%epsIrs )
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_AllocateMemory

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_AllocateMemory.F90,v $
! Revision 1.24  2008/12/06 08:44:25  mtcampbe
! Updated license.
!
! Revision 1.23  2008/11/19 22:17:36  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.22  2006/03/25 23:25:27  wasistho
! adjusted size of all si,j,k
!
! Revision 1.21  2006/03/24 05:03:27  wasistho
! increased size of si, sj, sk
!
! Revision 1.20  2006/03/18 08:20:58  wasistho
! allocate patch%arclen1,2
!
! Revision 1.19  2006/03/16 22:02:08  wasistho
! bug fixed aero coeffs allocation
!
! Revision 1.18  2006/03/13 03:40:40  wasistho
! allocate aero coeffs
!
! Revision 1.17  2006/03/08 23:33:23  wasistho
! made gridOld%xyzOld available for all type gm
!
! Revision 1.16  2006/03/02 01:27:17  wasistho
! split movegrid_epde to elglobal and elframe
!
! Revision 1.15  2006/02/09 00:26:23  wasistho
! allocate nCorns for EPDE too
!
! Revision 1.14  2005/12/07 08:47:53  wasistho
! allocate arrays for surface mesh motion EPDE
!
! Revision 1.13  2005/12/03 09:33:00  wasistho
! allocate arrays for movegrid EPDE
!
! Revision 1.12  2005/11/28 20:04:14  wasistho
! allocate gridOld%xyzOld
!
! Revision 1.11  2005/11/11 07:18:18  wasistho
! allocate ijkDgen
!
! Revision 1.10  2005/10/31 21:09:36  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.9  2005/10/28 07:37:55  wasistho
! corrected nullify cfcI,J,K
!
! Revision 1.8  2005/10/27 05:12:01  wasistho
! allocate xyzOrth
!
! Revision 1.7  2005/10/20 06:49:47  wasistho
! allocate cfcI,J,K
!
! Revision 1.6  2005/06/10 19:32:40  wasistho
! move memory allocation of movegridFrame variables to modflo/RFLO_ModMoveGridFrame
!
! Revision 1.5  2005/06/02 03:21:32  wasistho
! shuffle MoveGridVms with MoveGridFrame
!
! Revision 1.4  2005/05/29 20:45:01  wasistho
! added regCornOrig
!
! Revision 1.3  2005/05/28 08:07:40  wasistho
! allocate data for moveGridFrame
!
! Revision 1.2  2005/05/21 05:08:29  wasistho
! changed == MOVEGRID_GLOBAL to /=MOVEGRID_BLOCKS
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.33  2004/08/03 00:33:03  wasistho
! added allocation of grid%c2eConI,J,K
!
! Revision 1.32  2004/08/02 19:32:24  wasistho
! changed grid%avgCo to grid%c2fCo
!
! Revision 1.31  2004/07/30 17:25:18  wasistho
! provide cell2face averaging coefficients
!
! Revision 1.30  2003/11/20 16:40:36  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.26  2003/10/03 20:18:07  wasistho
! added wall dist. criterion for cellcentroids allocation
!
! Revision 1.25  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.24  2003/08/11 21:51:18  jblazek
! Added basic global grid smoothing scheme.
!
! Revision 1.23  2003/07/03 21:48:45  jblazek
! Implemented dual-time stepping.
!
! Revision 1.22  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.21  2003/04/23 22:51:28  olawlor
! Added TARGET attribute to region variable, since
! pointers are cached into it.
!
! Revision 1.20  2003/03/14 22:05:11  jblazek
! Improved mesh motion algorithm - node movement exchaged between blocks.
!
! Revision 1.19  2002/11/02 01:56:04  wasistho
! Added TURB statistics
!
! Revision 1.18  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.17  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.16  2002/08/30 01:47:58  jblazek
! Added support for moving grids.
!
! Revision 1.15  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
! Revision 1.14  2002/06/14 21:38:45  wasistho
! Added time avg statistics
!
! Revision 1.13  2002/03/29 23:15:22  jblazek
! Corrected bug in MPI send.
!
! Revision 1.12  2002/03/22 19:26:23  jblazek
! Unused pointers are nullified.
!
! Revision 1.11  2002/03/18 23:11:32  jblazek
! Finished multiblock and MPI.
!
! Revision 1.10  2002/02/27 18:38:20  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
! Revision 1.9  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.8  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.7  2002/02/06 00:15:39  jblazek
! Improved injection BC. Added pointers to gradients.
!
! Revision 1.6  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.5  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.4  2002/01/10 00:02:07  jblazek
! Added calculation of mixture properties.
!
! Revision 1.3  2002/01/08 22:09:17  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.2  2001/12/22 00:09:39  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1  2001/12/11 21:59:29  jblazek
! memory allocation added.
!
!******************************************************************************







