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
! Purpose: read in user input related to inflow boundary condition.
!
! Description: present inflow bc is based on prescribed velocities and 
!              either temperature or pressure.
!
! Input: boundary condition file.
!
! Output: regions = BC data.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ReadBcInflowVelSection.F90,v 1.11 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReadBcInflowVelSection( regions,bcTitle )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadPatchSection, RFLO_ReadBcFromFile
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  INTEGER :: bcTitle

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  INTEGER, PARAMETER :: NVALS_MAX = 8

  CHARACTER(10)  :: keys(NVALS_MAX)
  CHARACTER(256) :: fname

  INTEGER :: nvals, brbeg, brend, prbeg, prend, distrib, switch
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: lbound, nijk

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL) :: vals(NVALS_MAX)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ReadBcInflowVelSection',&
  'RFLO_ReadBcInflowVelSection.F90' )

! specify keywords and search for them ----------------------------------------

  keys(1) = 'TYPE'
  keys(2) = 'VELX'
  keys(3) = 'VELY'
  keys(4) = 'VELZ'
  IF (bcTitle==BC_INFLOW_VELTEMP) THEN
    keys(5) = 'TEMP'
    keys(6) = 'PRESS'
  ELSEIF (bcTitle==BC_INFLOW_VELPRESS) THEN
    keys(5) = 'PRESS'
    keys(6) = 'TEMP'
  ELSE
    CALL ErrorStop( global,ERR_UNKNOWN_BC,&
    __LINE__ )
  ENDIF
  keys(7) = 'RECYCTURB'
  keys(8) = 'AMPLITUDE'

  nvals = NVALS_MAX

  CALL ReadPatchSection( global,IF_INPUT,nvals,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,fname,defined )


! module consistency check ----------------------------------------------------

  IF (vals(7) > 0.1_RFREAL) THEN
#ifdef STATS
    IF ((global%flowType == FLOW_UNSTEADY) .AND. &
        (global%doStat == ACTIVE) .AND. &
        (global%mixtNStat > 0)) THEN
    ELSE
      CALL ErrorStop( global,ERR_VAL_BCVAL,&
      __LINE__, &
                 'inflow bc parameter RECYCTURB > 0 inconsistent with STATS' )
    ENDIF
#else
    CALL ErrorStop( global,ERR_VAL_BCVAL,&
    __LINE__, &
                   'Activate STATS for inflow bc parameter RECYCTURB > 0' )
#endif
  ENDIF

! get switches & check if all necessary values defined ------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_INFLOW .AND. &
           patch%bcType<=BC_INFLOW+BC_RANGE)    .AND. &   ! my boundary type,
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        patch%bcType = bcTitle

        IF (patch%mixt%bcSet .eqv. .true.) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,&
          __LINE__,'Inflow boundary.' )

        patch%mixt%nSwitches = 3
        IF (patch%bcCoupled == BC_EXTERNAL) THEN   ! data from outside
          patch%mixt%distrib = BCDAT_DISTRIB    ! => always distribution
        ELSE
          patch%mixt%distrib = distrib
        ENDIF

        ALLOCATE( patch%mixt%switches(patch%mixt%nSwitches), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
        __LINE__ )

! ----- check if switch defined
        IF (defined(1).eqv..true.) THEN
          patch%mixt%switches(BCSWI_INFLOW_TYPE)   = BCOPT_SUBSONIC
          IF (vals(1) < 0.1) &
            patch%mixt%switches(BCSWI_INFLOW_TYPE) = BCOPT_SUPERSONIC
          IF (vals(1) > 1.9) &
            patch%mixt%switches(BCSWI_INFLOW_TYPE) = BCOPT_MIXED
        ELSE
          CALL ErrorStop( global,ERR_NO_BCSWITCH,&
          __LINE__,'(inflow type).' )
        ENDIF

        patch%mixt%switches(BCSWI_INFLOW_MODEL) = BCOPT_STEADY
        IF (defined(7).eqv..true.) THEN
          IF (vals(7) > 0.1) &
            patch%mixt%switches(BCSWI_INFLOW_MODEL) = BCOPT_UNSTEADY
        ENDIF

! ----- check if appropriate values specified
        IF (patch%mixt%switches(BCSWI_INFLOW_TYPE) == BCOPT_SUBSONIC) THEN
          IF (patch%mixt%distrib==BCDAT_CONSTANT .AND. &
              (.NOT. (defined(2).eqv..true.) .OR. &
               .NOT. (defined(3).eqv..true.) .OR. &
               .NOT. (defined(4).eqv..true.) .OR. &
               .NOT. (defined(5).eqv..true.))) CALL ErrorStop( global,ERR_BCVAL_MISSING,&
               __LINE__ )
        ENDIF
        IF (patch%mixt%switches(BCSWI_INFLOW_TYPE) == BCOPT_SUPERSONIC .OR. &
            patch%mixt%switches(BCSWI_INFLOW_TYPE) == BCOPT_MIXED) THEN
          IF (patch%mixt%distrib==BCDAT_CONSTANT .AND. &
              (.NOT. (defined(2).eqv..true.) .OR. &
               .NOT. (defined(3).eqv..true.) .OR. &
               .NOT. (defined(4).eqv..true.) .OR. &
               .NOT. (defined(5).eqv..true.) .OR. &
               .NOT. (defined(6).eqv..true.))) CALL ErrorStop( global,ERR_BCVAL_MISSING,&
               __LINE__ )
        ENDIF
        IF (defined(7).eqv..true.) THEN
          IF (patch%mixt%switches(BCSWI_INFLOW_MODEL) == BCOPT_UNSTEADY) THEN
            IF (defined(8).eqv..true.) THEN
              patch%mixt%amplitude = ABS( vals(8) )
            ELSE
              CALL ErrorStop( global,ERR_BCVAL_MISSING,&
              __LINE__, &
                           'AMPLITUDE should be specified in inflow with recycturb' )
            ENDIF
          ENDIF
        ENDIF

! ----- set flag to BC specified
        patch%mixt%bcSet = .true.

      ENDIF   ! my BC & processor, active
    ENDDO
  ENDDO

! copy values/distribution to variables ---------------------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_INFLOW .AND. &
           patch%bcType<=BC_INFLOW+BC_RANGE)    .AND. &   ! my boundary type,
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        switch = patch%mixt%switches(BCSWI_INFLOW_TYPE)
        IF (switch == BCOPT_SUBSONIC) THEN
! -----  BUGFIX - Didn't allocate enough for this BC (MTC)
          patch%mixt%nData = 5
        ELSE
          patch%mixt%nData = 5
        ENDIF

! ----- allocate memory for the values

        IF (patch%mixt%distrib == BCDAT_DISTRIB) THEN
          n1    = ABS(patch%l1end-patch%l1beg)
          n2    = ABS(patch%l2end-patch%l2beg)
          iOff  = n1 + 1
          ijBeg = IndIJ( 0, 0,iOff)
          ijEnd = IndIJ(n1,n2,iOff)
        ELSE
          ijBeg = 0
          ijEnd = 1
        ENDIF
        ALLOCATE( patch%mixt%vals(patch%mixt%nData,ijBeg:ijEnd), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
        __LINE__ )

! ----- distribution from file

        IF (patch%mixt%distrib==BCDAT_DISTRIB .AND. &
            patch%bcCoupled      /=BC_EXTERNAL  ) THEN
           WRITE(*,*) 'Reading more bc from file',fname
          CALL RFLO_ReadBcFromFile( global,fname,patch )

! ----- distribution from external source / constant value

        ELSE

          patch%mixt%vals(BCDAT_INFLOW_U,:) = vals(2)
          patch%mixt%vals(BCDAT_INFLOW_V,:) = vals(3)
          patch%mixt%vals(BCDAT_INFLOW_W,:) = vals(4)
          IF (bcTitle==BC_INFLOW_VELTEMP) THEN
            patch%mixt%vals(BCDAT_INFLOW_T,:) = vals(5)
            IF (defined(6).eqv..true.) THEN
               patch%mixt%vals(BCDAT_INFLOW_P,:) = vals(6)
            ENDIF
          ELSEIF (bcTitle==BC_INFLOW_VELPRESS) THEN
            patch%mixt%vals(BCDAT_INFLOW_P,:) = vals(5)
            IF (defined(6).eqv..true.) &
            patch%mixt%vals(BCDAT_INFLOW_T,:) = vals(6)
          ENDIF
        ENDIF  ! distribution?

      ENDIF    ! bcType, active region on my processor

    ENDDO      ! iPatch
  ENDDO        ! iReg

! search for max ijk-index for turbulence recycling ---------------------------

  IF (vals(7) > 0.1_RFREAL) THEN
    nijk = global%infloNijk
    DO iReg=brbeg,brend
      CALL RFLO_GetDimensPhysNodes( regions(iReg),1,ipnbeg,ipnend, &
                                    jpnbeg,jpnend,kpnbeg,kpnend )
      DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

        patch => regions(iReg)%levels(1)%patches(iPatch)

        IF ((patch%bcType>=BC_INFLOW .AND. &
             patch%bcType<=BC_INFLOW+BC_RANGE)    .AND. &   ! my boundary type,
            regions(iReg)%procid==global%myProcid .AND. &   ! region active and
            regions(iReg)%active==ACTIVE) THEN              ! on my processor
          lbound = patch%lbound
          IF (lbound==1 .OR. lbound==2) THEN
            nijk = ABS(ipnend-ipnbeg)+1
          ELSEIF (lbound==3 .OR. lbound==4) THEN 
            nijk = ABS(jpnend-jpnbeg)+1
          ELSEIF (lbound==5 .OR. lbound==6) THEN 
            nijk = ABS(kpnend-kpnbeg)+1
          ENDIF
        ENDIF  ! bcType
        global%infloNijk = MIN( global%infloNijk,nijk )
      ENDDO  ! iPatch
    ENDDO    ! iReg
    global%infloNijk = global%infloNijk/2+1
  ENDIF      ! recycTurb

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ReadBcInflowVelSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReadBcInflowVelSection.F90,v $
! Revision 1.11  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2008/10/23 18:20:53  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.8  2006/08/19 15:38:12  mparmar
! Renamed patch variables
!
! Revision 1.7  2006/04/13 21:37:43  wasistho
! take half of global%infloNijk
!
! Revision 1.6  2006/02/18 08:26:33  wasistho
! fixed reading temperature and pressure
!
! Revision 1.5  2005/11/19 03:25:36  wasistho
! bug fixed errorstop reading amplitude
!
! Revision 1.4  2005/11/18 07:18:42  wasistho
! added parameter AMPLITUDE
!
! Revision 1.3  2005/09/20 23:56:57  wasistho
! modified computation of global%infloNijk
!
! Revision 1.2  2005/09/20 23:03:22  wasistho
! added user input RECYCTURB
!
! Revision 1.1  2005/04/28 05:48:28  wasistho
! added velocity based inflow BC
!
!
!******************************************************************************







