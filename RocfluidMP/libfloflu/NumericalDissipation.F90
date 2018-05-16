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
! ******************************************************************************
!
! Purpose: Compute fluxes related to numerical dissipation.
!
! Description: 
!   None.
!
! Input: 
!   region     data of current region
!
! Output: 
!   None 
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: NumericalDissipation.F90,v 1.6 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2001-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE NumericalDissipation(region)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global  
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModParameters

#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_CentralDissipation, RFLO_RoeDissipFirst, &
                            RFLO_RoeDissipSecond
#endif

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), TARGET :: region

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: spaceDiscr, spaceOrder
  
  TYPE(t_global), POINTER :: global  
#ifdef RFLU
  TYPE(t_region), POINTER :: pRegion
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'NumericalDissipation',&
  'NumericalDissipation.F90' )

#ifdef RFLU
  pRegion => region
#endif

! ******************************************************************************
! Get dimensions and pointers
! ******************************************************************************

  spaceDiscr = region%mixtInput%spaceDiscr
  spaceOrder = region%mixtInput%spaceOrder

#ifdef RFLO
! ******************************************************************************
! 2nd-order central scheme 
! ******************************************************************************

  IF (spaceDiscr==DISCR_CEN_SCAL .AND. &
      (spaceOrder==DISCR_ORDER_1 .OR. spaceOrder==DISCR_ORDER_2)) THEN
    CALL RFLO_CentralDissipation( region )
  ENDIF
#endif

#ifdef RFLO
! ******************************************************************************
! Roe upwind scheme-
! ******************************************************************************

  IF (spaceDiscr == DISCR_UPW_ROE) THEN
    IF (spaceOrder == DISCR_ORDER_1) THEN
      CALL RFLO_RoeDissipFirst( region )
    ELSE IF (spaceOrder == DISCR_ORDER_2) THEN
      CALL RFLO_RoeDissipSecond( region )
    ENDIF
  ENDIF
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(region%global)

END SUBROUTINE NumericalDissipation

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: NumericalDissipation.F90,v $
! Revision 1.6  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/05/01 20:58:56  haselbac
! Not needed for RFLU
!
! Revision 1.3  2006/03/26 20:21:19  haselbac
! Rewrite bcos of of GL model
!
! Revision 1.2  2005/05/16 20:39:27  haselbac
! Changed args to pRegion, now USE flux moduless
!
! Revision 1.1  2004/12/01 16:49:45  haselbac
! Initial revision after changing case
!
! Revision 1.18  2004/01/29 22:52:42  haselbac
! Modified for unsteady flows
!
! Revision 1.17  2003/12/04 03:23:01  haselbac
! Added second-order schemes
!
! Revision 1.16  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.13  2003/05/29 17:28:42  jblazek
! Implemented Roe scheme.
!
! Revision 1.12  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.11  2003/04/29 22:52:35  jferry
! removed call to PEUL_CentralDissipation
!
! Revision 1.10  2003/03/11 16:04:19  jferry
! Enclosed USE statements for multi-physics routines within ifdefs
!
! Revision 1.9  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
! Revision 1.8  2002/09/09 15:19:31  haselbac
! mixtInput now under regions
!
! Revision 1.7  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.6  2002/07/29 17:13:08  jblazek
! Clean up after RFLU and TURB.
!
! Revision 1.5  2002/07/25 14:50:18  haselbac
! Modified logic because of optimal LES fluxes
!
! Revision 1.4  2002/05/04 16:36:16  haselbac
! Added RFLU statements
!
! Revision 1.3  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.2  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.1  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
!
! ******************************************************************************







