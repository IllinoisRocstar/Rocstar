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
! Purpose: delete data windows.
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
! $Id: Fluid_finalize.F90,v 1.3 2008/12/06 08:44:00 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE Fluid_finalize( globalGenx )

  USE ModRocstar, ONLY       : t_globalGenx
#ifdef RFLO  
  USE ModInterfaces, ONLY: RFLO_EndFlowSolver
#endif
#ifdef RFLU
  USE ModInterfaces, ONLY: RFLU_EndFlowSolver
#endif

  IMPLICIT NONE
  INCLUDE 'comf90.h'

! ... parameters
  TYPE(t_globalGenx), POINTER  :: globalGenx

!******************************************************************************

#ifdef RFLO
  CALL RFLO_EndFlowSolver( globalGenx%regions )
#endif
#ifdef RFLU
  CALL RFLU_EndFlowSolver(globalGenx%levels)
#endif

  CALL COM_delete_window( TRIM(globalGenx%global%winName)//'_surf' )
  CALL COM_delete_window( TRIM(globalGenx%global%winName)//'_vol' )

END SUBROUTINE Fluid_finalize

!******************************************************************************
!
! RCS Revision history:
!
! $Log: Fluid_finalize.F90,v $
! Revision 1.3  2008/12/06 08:44:00  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:23:41  haselbac
! Initial revision after changing case
!
! Revision 1.6  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.3  2002/10/07 14:09:53  haselbac
! Added RFLU/RFLO distinction
!
! Revision 1.2  2002/09/25 17:55:41  jblazek
! Added call to EndFlowSolver.
!
! Revision 1.1  2002/09/20 22:22:34  jblazek
! Finalized integration into GenX.
!
!******************************************************************************






