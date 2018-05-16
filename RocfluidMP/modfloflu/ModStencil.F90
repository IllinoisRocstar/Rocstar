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
! Purpose: Define the derived data type encapsulating stencils.
!
! Description: None.
!
! Notes: None
!
! ******************************************************************************
!
! $Id: ModStencil.F90,v 1.4 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

MODULE ModStencil

  USE ModDataTypes
  USE ModParameters
  
  IMPLICIT NONE

  TYPE t_stencil
    INTEGER :: nBFaceMembs,nCellMembs,nLayers
    INTEGER, DIMENSION(:), POINTER :: cellMembs
    INTEGER, DIMENSION(:,:), POINTER :: bFaceMembs,layerInfo
    REAL(RFREAL), DIMENSION(:), POINTER :: xyzMoms
  END TYPE t_stencil

  TYPE t_stencilInfo
    INTEGER :: nBFaceMembsMax,nCellMembsMax,nCellMembsMin,nLayersMax, & 
               orderNominal
  END TYPE t_stencilInfo

END MODULE ModStencil

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModStencil.F90,v $
! Revision 1.4  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/01/03 16:18:06  haselbac
! Changed t_stencil_info type, cosmetics
!
! Revision 1.1  2003/12/04 03:28:19  haselbac
! Initial revision
!
! ******************************************************************************






