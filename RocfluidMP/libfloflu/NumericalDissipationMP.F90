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
! Purpose: calculate numerical dissipation for RocfluidMP framework.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: regions%levels%mixt = new solution after one time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: NumericalDissipationMP.F90,v 1.3 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE NumericalDissipationMP( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE ModInterfaces, ONLY: NumericalDissipation

#ifdef PEUL
  USE ModInterfacesEulerian, ONLY : PEUL_CentralDissipation
#endif
#ifdef TURB
  USE ModInterfacesTurbulence, ONLY : TURB_RansNumericalDissipation
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... local variables
  INTEGER :: spaceDiscr,spaceOrder
  TYPE(t_global), POINTER :: global

#ifdef RFLU
  TYPE(t_region), POINTER :: pRegion
#endif

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'NumericalDissipationMP',&
  'NumericalDissipationMP.F90' )

! get parameters --------------------------------------------------------------

  spaceDiscr = region%mixtInput%spaceDiscr
  spaceOrder = region%mixtInput%spaceOrder

! compute numerical dissipation -----------------------------------------------

  CALL NumericalDissipation( region )

#ifdef PEUL
!  IF (global%peulUsed) CALL PEUL_NumericalDissipation( region )
  IF (global%peulUsed) CALL PEUL_CentralDissipation( region )
#endif

#ifdef TURB
  IF (region%mixtInput%flowModel == FLOW_NAVST .AND. &
      region%mixtInput%turbModel /= TURB_MODEL_NONE) &
    CALL TURB_RansNumericalDissipation( region )
#endif

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE NumericalDissipationMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: NumericalDissipationMP.F90,v $
! Revision 1.3  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:49:50  haselbac
! Initial revision after changing case
!
! Revision 1.14  2004/03/20 00:26:22  wasistho
! set turb_rflo_ransNumericalDiss to turb_ransNumerical..
!
! Revision 1.13  2004/03/11 03:32:55  wasistho
! changed rocturb nomenclature
!
! Revision 1.12  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.11  2004/01/29 22:52:43  haselbac
! Removed scalar dissipation routines
!
! Revision 1.10  2003/11/25 21:01:42  haselbac
! Added rocspecies support
!
! Revision 1.9  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/10/27 04:49:19  wasistho
! replace ransCentralDiss by ransNumericalDiss.
!
! Revision 1.5  2003/10/03 20:12:43  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.4  2003/10/01 23:52:09  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.3  2003/09/26 21:45:34  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.2  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.1  2003/03/28 19:46:30  fnajjar
! Initial import for RocfluidMP
!
!******************************************************************************







