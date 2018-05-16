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
! Purpose: Define the derived data type encapsulating coloring.
!
! Description: None.
!
! Notes: None
!
! ******************************************************************************
!
! $Id $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE ModColoring

  USE ModDataTypes
  
  IMPLICIT NONE

  TYPE t_soc
    INTEGER :: nCellMembs
    INTEGER, DIMENSION(:), POINTER :: cellMembs
  END TYPE t_soc

END MODULE ModColoring

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModColoring.F90,v $
! Revision 1.3  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2005/08/17 20:04:15  hdewey2
! Initial revision
!
! ******************************************************************************






