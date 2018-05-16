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
! Purpose: Suite for vector and tensor/matrix operations
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModVectorTensor.F90,v 1.3 2008/12/06 08:44:17 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModVectorTensor

  USE ModGlobal, ONLY    : t_global 
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY      : t_grid
  USE ModBndPatch, ONLY  : t_patch
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLO_NormCrossProd

! private :
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModVectorTensor.F90,v $ $Revision: 1.3 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS

!******************************************************************************
!
! Purpose: perform cross product and normalize it
!
! Description: none.
!
! Input: s1, s2 = pair of vectors to be processed
!
! Output: s3 = s1 x s2
!
! Notes: none
!
!******************************************************************************

SUBROUTINE RFLO_NormCrossProd( s1,s2,s3 )

  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: s1(:), s2(:), s3(:)

! ... loop variables

! ... local variables
  REAL(RFREAL) :: smag

!******************************************************************************

  s3(XCOORD) = s1(YCOORD)*s2(ZCOORD)-s1(ZCOORD)*s2(YCOORD)
  s3(YCOORD) = s1(ZCOORD)*s2(XCOORD)-s1(XCOORD)*s2(ZCOORD)
  s3(ZCOORD) = s1(XCOORD)*s2(YCOORD)-s1(YCOORD)*s2(XCOORD)
  smag       = SQRT( s3(XCOORD)**2+s3(YCOORD)**2+s3(ZCOORD)**2 )
  s3(:)      = s3(:)/smag

END SUBROUTINE RFLO_NormCrossProd

! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModVectorTensor

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModVectorTensor.F90,v $
! Revision 1.3  2008/12/06 08:44:17  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:28  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2005/12/05 10:53:09  wasistho
! initial import
!
!
!
! ******************************************************************************






