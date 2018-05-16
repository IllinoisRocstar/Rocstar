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
! Purpose: set dummycells to zero.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: regions%levels%mixt%rhs = residual in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ZeroDummyCellsMP.F90,v 1.5 2008/12/06 08:44:11 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ZeroDummyCellsMP( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_ZeroDummyCells
#endif
#ifdef RFLU
  USE ModInterfaces, ONLY: RFLU_ZeroVirtualCellVars
#endif
#ifdef TURB
  USE TURB_ModParameters
  USE ModInterfacesTurbulence, ONLY: TURB_RansZeroDummyCells
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), TARGET :: region

! ... local variables
#ifdef RFLO
  INTEGER :: iLev
#endif
  LOGICAL :: peulUsed,specUsed

  REAL(RFREAL), POINTER :: rhs(:,:), rhsPeul(:,:)

  TYPE(t_global), POINTER :: global
#ifdef RFLU
  TYPE(t_region), POINTER :: pRegion
#endif

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'ZeroDummyCellsMP',&
  'ZeroDummyCellsMP.F90' )

! get dimensions and pointers =================================================

  peulUsed = global%peulUsed
  specUsed = global%specUsed

#ifdef RFLO
  iLev = region%currLevel
  rhs  => region%levels(iLev)%mixt%rhs
#ifdef PEUL
  rhsPEul => region%levels(iLev)%peul%rhs
#endif
#endif

#ifdef RFLU
  pRegion => region
#endif

! zero out residuals in dummy cells -------------------------------------------

#ifdef RFLO
  CALL RFLO_ZeroDummyCells( region,rhs )
#ifdef PEUL
  IF (peulUsed) CALL RFLO_ZeroDummyCells( region,rhsPEul )
#endif
#ifdef TURB
  IF (region%mixtInput%flowModel == FLOW_NAVST .AND. &
      region%mixtInput%turbModel /= TURB_MODEL_NONE .AND. &
      region%turbInput%modelClass == MODEL_RANS) THEN
    CALL TURB_RansZeroDummyCells( region )
  ENDIF
#endif
#endif

#ifdef RFLU
  CALL RFLU_ZeroVirtualCellVars(pRegion,pRegion%mixt%rhs)
#ifdef SPEC
  IF ( specUsed .EQV. .TRUE. ) THEN
    CALL RFLU_ZeroVirtualCellVars(pRegion,pRegion%spec%rhs)
  END IF ! specUsed
#endif
#ifdef TURB
  IF ( pRegion%mixtInput%flowModel == FLOW_NAVST .AND. &
       pRegion%mixtInput%turbModel /= TURB_MODEL_NONE .AND. &
       pRegion%turbInput%modelClass == MODEL_RANS) THEN
    CALL RFLU_ZeroVirtualCellVars(pRegion,pRegion%turb%rhs)
  END IF ! pRegion%mixtInput%flowModel
#endif
#endif

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE ZeroDummyCellsMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ZeroDummyCellsMP.F90,v $
! Revision 1.5  2008/12/06 08:44:11  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:24  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/01/12 02:38:30  wasistho
! added modelClass condition in turb
!
! Revision 1.2  2005/05/16 20:40:32  haselbac
! Renamed RFLU_ZeroDummyCells, now use pRegion
!
! Revision 1.1  2004/12/01 16:52:32  haselbac
! Initial revision after changing case
!
! Revision 1.14  2004/07/28 15:29:18  jferry
! created global variable for spec use
!
! Revision 1.13  2004/03/19 02:40:15  wasistho
! renamed TURB_RFLO_RansZeroDummyCells to TURB_RansZeroDummyCells
!
! Revision 1.12  2004/03/11 03:33:06  wasistho
! changed rocturb nomenclature
!
! Revision 1.11  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.10  2004/03/02 21:50:44  jferry
! Corrected typo
!
! Revision 1.9  2003/11/25 21:01:46  haselbac
! Added support for rocspecies, changed TURB call
!
! Revision 1.8  2003/11/20 16:40:36  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/10/16 20:14:10  wasistho
! turbulence zero dummy cells computed through one routine
!
! Revision 1.4  2003/10/03 20:13:21  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.3  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.2  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.1  2003/03/28 19:48:44  fnajjar
! Initial import for RocfluidMP
!
!******************************************************************************







