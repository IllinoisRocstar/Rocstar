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
! Purpose: Invert general m-by-n matrix A using SVD
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   nRows       Number of rows of A
!   nCols       Number of columns of A
!   a           Matrix A, to be inverted
!
! Output:
!   aInv        Inverse of A   
!   sCount      Number of singular values below threshold 
!
! Notes:
!   1. Uses LAPACK routines to carry out SVD.
!
! ******************************************************************************
!
! $Id: RFLU_InvertMatrixSVD.F90,v 1.9 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_InvertMatrixSVD(global,nRows,nCols,a,aInv,sCount)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
 
  USE ModSortSearch
 
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: nRows,nCols
  INTEGER, INTENT(OUT), OPTIONAL :: sCount
  REAL(RFREAL) :: a(nRows,nCols),aInv(nCols,nRows)
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,i,info,j,workArrayRealSize
  INTEGER :: workArrayInt(8*nCols)  
  REAL(RFREAL) :: s(nCols)
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: workArrayReal
  REAL(RFREAL) :: sInv(nCols,nRows),u(nRows,nRows),vt(nCols,nCols), & 
                  v(nCols,nCols)

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InvertMatrixSVD.F90,v $ $Revision: 1.9 $'

  CALL RegisterFunction(global,'RFLU_InvertMatrixSVD',&
  'RFLU_InvertMatrixSVD.F90')

! ******************************************************************************
! Set work array size and allocate memory
! ******************************************************************************

  workArrayRealSize = 2*(4*nCols*nCols + nRows + 9*nCols)

  ALLOCATE(workArrayReal(workArrayRealSize),STAT=errorFlag)
  global%error = errorFlag    
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'workArrayReal')
  END IF ! global%error  

! ******************************************************************************
! Perform SVD
! ******************************************************************************
            
  CALL dgesdd('A',nRows,nCols,a,nRows,s,u,nRows,vt,nCols,workArrayReal, & 
              workArrayRealSize,workArrayInt,info)
  global%error = info
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_LAPACK_OUTPUT,__LINE__)
  END IF ! global%error

! ******************************************************************************
! Construct actual singular-value matrix
! ******************************************************************************

  sCount = 0
  
  DO i = 1,nCols ! Explicit loop to avoid Frost problem
    DO j = 1,nRows
      sInv(i,j) = 0.0_RFREAL
    END DO ! j
  END DO ! i

! Based on absolute size
!  DO i = 1,nCols  
!    IF ( s(i) > 100.0_RFREAL*EPSILON(1.0_RFREAL) ) THEN   
!      sInv(i,i) = 1.0_RFREAL/s(i)
!    ELSE 
!      sInv(i,i) = 0.0_RFREAL
!      sCount    = sCount + 1
!    END IF ! s
!  END DO ! i

! Based on ratio of successive singular values
  IF ( s(1) > 100.0_RFREAL*EPSILON(1.0_RFREAL) ) THEN 
    sInv(1,1) = 1.0_RFREAL/s(1)
  ELSE 
    sInv(1,1) = 0.0_RFREAL
    sCount = 1
  END IF ! s

  outerLoop: DO i = 2,nCols  
    IF ( s(i) > 0.01_RFREAL*s(i-1) ) THEN   
      sInv(i,i) = 1.0_RFREAL/s(i)
    ELSE 
      sInv(i,i) = 0.0_RFREAL
      sCount    = sCount + 1
      
      DO j = i+1,nCols
        sInv(j,j) = 0.0_RFREAL
        sCount = sCount + 1
      END DO ! j
      
      EXIT outerLoop
    END IF ! s
  END DO outerLoop
 
! ******************************************************************************
! Compute inverse
! ******************************************************************************

  aInv = MATMUL(TRANSPOSE(vt),MATMUL(sInv,TRANSPOSE(u))) ! Lapack

! ******************************************************************************
! End
! ******************************************************************************

  DEALLOCATE(workArrayReal,STAT=errorFlag)
  global%error = errorFlag    
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'workArrayReal')
  END IF ! global%error  

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InvertMatrixSVD

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InvertMatrixSVD.F90,v $
! Revision 1.9  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.6  2005/10/05 13:51:42  haselbac
! Cosmetics
!
! Revision 1.5  2003/12/04 03:23:56  haselbac
! Fixed bug, changed init, added new method of detecting singularity
!
! Revision 1.4  2003/07/22 01:56:50  haselbac
! Change to init of sInv and cosmetics
!
! Revision 1.3  2002/10/08 15:48:56  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.2  2002/09/09 14:15:01  haselbac
! global now under regions, bug fix: workArrayReal not deallocated
!
! Revision 1.1  2002/07/25 14:34:59  haselbac
! Initial revision
!
! ******************************************************************************







