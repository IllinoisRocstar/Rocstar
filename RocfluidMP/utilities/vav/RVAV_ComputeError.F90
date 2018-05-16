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
! Purpose: Computes Error Between Stream 1 and Stream 2 Data.
!
! Description: none.
!
! Input: Arrays evS1 and evS2 and their corresponding extents
!
! Output: L2 and L8 Norms
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RVAV_ComputeError.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RVAV_ComputeError( global,iCompare,iNodes,jNodes,kNodes )                              

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY : t_global
  USE ModParameters
  USE RVAV_ModGlobal
  USE RVAV_ModParameters
  IMPLICIT NONE
  
! ... parameter
  TYPE(t_global), POINTER :: global

  INTEGER :: iCompare, iNodes, jNodes, kNodes

! ... loop variables
  INTEGER :: i, j, k 

! ... local variables
  CHARACTER(CHRLEN) :: msg

  REAL(RFREAL) :: l2Norm, l8Norm
  REAL(RFREAL) :: l2NormS2, l8NormS2
  
!******************************************************************************

  CALL RegisterFunction( global, 'RVAV_ComputeError',&
  'RVAV_ComputeError.F90' )

! initialize variables

  l2Norm   =  0.0_RFREAL
  l8Norm   = -1.0E+30_RFREAL
  l2NormS2 =  0.0_RFREAL
  l8NormS2 = -1.0E+30_RFREAL
  
  DO k = 1, kNodes
    DO j = 1, jNodes
      DO i = 1, iNodes
        l2Norm = l2Norm + (globalRVAV%evS1(i,j,k) -globalRVAV%evS2(i,j,k))**2
        l8Norm = MAX(l8Norm,ABS(globalRVAV%evS1(i,j,k) - globalRVAV%evS2(i,j,k)) )
        l2NormS2 = l2NormS2 + globalRVAV%evS2(i,j,k)**2
        l8NormS2 = MAX(l8NormS2,ABS(globalRVAV%evS2(i,j,k))  )
      ENDDO ! i
    ENDDO ! j
  ENDDO ! k
  
  WRITE(*,1010) iCompare,SQRT(l2Norm)/SQRT(l2NormS2)*100.0_RFREAL, &
                         l8Norm/l8NormS2*100.0_RFREAL

  IF (global%verbLevel /= VERBOSE_NONE)  THEN
    WRITE(*,1000) iCompare,SQRT(l2Norm)/(REAL(iNodes*jNodes*kNodes, KIND=RFREAL)),l8Norm
  
    WRITE(200,*) 'ICompare = ',iCompare
  
    DO k = 1, kNodes
      DO j = 1, jNodes
        DO i = 1, iNodes
          WRITE(200,*) i,j,k,globalRVAV%evS1(i,j,k),globalRVAV%evS2(i,j,k), &
                             globalRVAV%evS1(i,j,k) -globalRVAV%evS2(i,j,k)
        ENDDO ! i
      ENDDO ! j
    ENDDO ! k
  END IF ! verbLevel 

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1000 FORMAT(' Comparison ',I5,', L2 Norm= ',1PE15.7,', L8 Norm =',1PE15.7)
1010 FORMAT(/,' RocVaV passed the test for Comparison ',I5,/,&
              '  Percentage of Normalized L2 Norm= ',1PE15.7,/,&
              '  Percentage of Normalized L8 Norm =',1PE15.7)
            
END SUBROUTINE RVAV_ComputeError

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_ComputeError.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:43:19  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.4  2002/09/10 00:01:45  f-najjar
! Variable global moved into regions()
!
! Revision 1.3  2002/06/19 14:40:15  f-najjar
! Included verbLevel calls for cleanup
!
! Revision 1.2  2002/06/19 13:59:56  f-najjar
! Computed Normalized L2 and L8 errors
!
! Revision 1.1.1.1  2002/06/03 21:41:29  f-najjar
! Initial Import of RocVaV
!
!******************************************************************************







