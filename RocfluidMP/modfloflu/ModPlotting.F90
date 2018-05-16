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
! Purpose: Define data types for plotting variables.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ModPlotting.F90,v 1.3 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
! ******************************************************************************

MODULE ModPlotting

  USE ModDataTypes
  USE ModParameters

  IMPLICIT NONE

  TYPE t_plot
    CHARACTER(CHRLEN), DIMENSION(PV_XXXX_NVAR) :: pvNameLong,pvNameShort
    INTEGER :: nPv
    INTEGER :: pv2pvi(PV_XXXX_NVAR)
    INTEGER, DIMENSION(:), POINTER :: pvi2pv
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pv,pvVert
  END TYPE t_plot

END MODULE ModPlotting

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModPlotting.F90,v $
! Revision 1.3  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/03/19 21:39:13  haselbac
! Initial revision
!
! ******************************************************************************






