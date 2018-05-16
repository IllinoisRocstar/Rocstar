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
! Purpose: define the overall data structure related to streams
!
! Description: none
!
! Notes: none
!
!******************************************************************************
!
! $Id: RVAV_ModDataStruct.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

MODULE RVAV_ModDataStruct

  USE ModDataTypes

  IMPLICIT NONE
  
! region related data for streams

  TYPE t_RVAV_region

    INTEGER :: iNodes, jNodes, kNodes, nVars
    REAL(RFREAL), DIMENSION(:,:), POINTER :: xyzS1,cvS1,dvS1
    REAL(RFREAL), DIMENSION(:,:), POINTER :: xyzS2,cvS2,dvS2
    REAL(RFREAL), DIMENSION(:,:,:,:), POINTER :: xyzAES2,cvAES2
   
  END TYPE t_RVAV_region

  
END MODULE RVAV_ModDataStruct

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_ModDataStruct.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1.1.1  2002/06/03 21:41:29  f-najjar
! Initial Import of RocVaV
!
!
!******************************************************************************






