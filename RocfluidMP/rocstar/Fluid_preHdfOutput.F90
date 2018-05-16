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
! Purpose: treatments before solution dump to Hdf files.
!
! Description: none.
!
! Input: globalGenx = global data structure (contains name of the window)
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: Fluid_preHdfOutput.F90,v 1.7 2008/12/06 08:44:00 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE Fluid_preHdfOutput( globalGenx )

  USE ModDataTypes
  USE ModRocstar, ONLY       : t_globalGenx
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModParameters

  IMPLICIT NONE
  INCLUDE 'comf90.h'

! ... parameters
  TYPE(t_globalGenx), POINTER  :: globalGenx
  TYPE(t_global), POINTER      :: global
  TYPE(t_region), POINTER      :: regions(:)

! ... loop variables
  INTEGER :: iReg, iStat

! ... local variables
  INTEGER :: iLev
  REAL(RFREAL) :: eps

!******************************************************************************

  global  => globalGenx%global

#ifdef RFLO
  regions => globalGenx%regions
#endif
#ifdef RFLU
  regions => globalGenx%levels(1)%regions
#endif

  CALL RegisterFunction( global,'Fluid_preHdfOutput',&
  'Fluid_preHdfOutput.F90' )

! set tav from accumulated values ---------------------------------------------
! to be moved to a wrapper

  eps = 100._RFREAL*EPSILON( 1._RFREAL )

#ifdef STATS

#ifdef RFLO
  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
      iLev = regions(iReg)%currLevel
#endif
#ifdef RFLU
  DO iReg=1,global%nRegionsLocal
#endif

      IF ((global%flowType==FLOW_UNSTEADY) .AND. (global%doStat==ACTIVE)) THEN
        IF (global%integrTime > eps) THEN
          IF (global%mixtNStat > 0) THEN
            DO iStat=1,global%mixtNStat

#ifdef RFLO
              regions(iReg)%levels(iLev)%mixt%tav(iStat,:) = &
              regions(iReg)%levels(iLev)%mixt%tav(iStat,:)/global%integrTime
#endif
#ifdef RFLU
              regions(iReg)%mixt%tav(iStat,:) = &
              regions(iReg)%mixt%tav(iStat,:)/global%integrTime
#endif

            ENDDO
          ENDIF  ! mixtNstat

#ifdef TURB
          IF ((global%turbActive .EQV. .true.) .AND. (global%turbNStat > 0)) THEN
            DO iStat=1,global%turbNStat

#ifdef RFLO
              regions(iReg)%levels(iLev)%turb%tav(iStat,:) = &
              regions(iReg)%levels(iLev)%turb%tav(iStat,:)/global%integrTime
#endif
#ifdef RFLU
              regions(iReg)%turb%tav(iStat,:) = &
              regions(iReg)%turb%tav(iStat,:)/global%integrTime
#endif

            ENDDO
          ENDIF  ! turbNstat
#endif

        ENDIF    ! integrTime
      ENDIF      ! unsteady and dostat

#ifdef RFLO
    ENDIF        ! region on this processor and active
  ENDDO          ! iReg
#endif
#ifdef RFLU
  ENDDO          ! iReg
#endif

#endif

  CALL DeregisterFunction( global )

END SUBROUTINE Fluid_preHdfOutput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: Fluid_preHdfOutput.F90,v $
! Revision 1.7  2008/12/06 08:44:00  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/01/07 10:19:46  wasistho
! added Rocflu treatment
!
! Revision 1.4  2005/12/08 19:56:44  wasistho
! modified purpose (cosmetics)
!
! Revision 1.3  2005/12/08 07:03:43  wasistho
! added ifdef TURB
!
! Revision 1.2  2005/12/08 02:04:57  wasistho
! added register deregister function
!
! Revision 1.1  2005/12/08 00:22:04  wasistho
! initial import
!
!
!******************************************************************************







