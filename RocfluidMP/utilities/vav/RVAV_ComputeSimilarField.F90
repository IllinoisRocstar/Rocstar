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
! Purpose: Computes the Similar Field for Stream 1 
!          for certain Analytical Solutions.
!
! Description: none.
!
! Input: Type of Solution, Array evS1 and their corresponding extents
!
! Output: Scale-Similar evS1
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RVAV_ComputeSimilarField.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RVAV_ComputeSimilarField( global, iNodes,jNodes,kNodes, &
                                     similarityType,variableIndex,ev )                              

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY : t_global
  USE ModParameters
  USE RVAV_ModGlobal
  USE RVAV_ModParameters
  IMPLICIT NONE
  
! ... parameter
  TYPE(t_global), POINTER :: global

  INTEGER :: iNodes,jNodes,kNodes
  INTEGER :: similarityType,variableIndex
  REAL(RFREAL), POINTER :: ev(:,:,:) ! this is an array for a single extracted
                                     ! variable with i,j,k dimensions

! ... loop variables
  INTEGER :: i, j, k 

! ... local variables
  CHARACTER(CHRLEN) :: msg

  REAL(RFREAL) :: uCenterLine,vInjection
  REAL(RFREAL) :: RhoRef, URef

!******************************************************************************

  CALL RegisterFunction( global, 'RVAV_ComputeSimilarField',&
  'RVAV_ComputeSimilarField.F90' )

  SELECT CASE (similarityType)
  
    CASE (RVAV_CULICK) 
    
      SELECT CASE (variableIndex)
         CASE (RVAV_U)
	     
             DO k = 1, kNodes
             DO i = 1, iNodes
               uCenterLine = globalRVAV%evS1(i,jNodes,k)
             DO j = 1, jNodes
               ev(i,j,k) = ev(i,j,k)/uCenterLine
             ENDDO ! j
             ENDDO ! i
             ENDDO ! k
         CASE (RVAV_V)
	     
             DO k = 1, kNodes
             DO i = 1, iNodes
               vInjection = globalRVAV%evS1(i,1,k)
             DO j = 1, jNodes
               ev(i,j,k) = ev(i,j,k)/vInjection
             ENDDO ! j
             ENDDO ! i
             ENDDO ! k 
      END SELECT ! variableIndex
      
    CASE (RVAV_BLASIUS)
    
      SELECT CASE (variableIndex)
         CASE (RVAV_RHO)
             DO k = 1, kNodes
             DO i = 1, iNodes
               RhoRef = globalRVAV%evS1(i,jNodes,k)
             DO j = 1, jNodes
               ev(i,j,k) = ev(i,j,k)/RhoRef
             ENDDO ! j
             ENDDO ! i
             ENDDO ! k 

         CASE (RVAV_U)
             DO k = 1, kNodes
             DO i = 1, iNodes
               URef = globalRVAV%evS1(i,jNodes,k)
             DO j = 1, jNodes
               ev(i,j,k) = ev(i,j,k)/URef
             ENDDO ! j
             ENDDO ! i
             ENDDO ! k

      END SELECT ! variableIndex
	
  END SELECT ! similarityType          

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RVAV_ComputeSimilarField

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_ComputeSimilarField.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:43:20  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.4  2002/10/21 18:48:55  f-najjar
! Bug fix for inconsistent calling sequence
!
! Revision 1.3  2002/09/10 00:01:45  f-najjar
! Variable global moved into regions()
!
! Revision 1.2  2002/06/24 15:49:33  f-najjar
! Included Similarity Scaling for Blasius
!
! Revision 1.1  2002/06/17 16:11:05  f-najjar
! Initial Import
!
!******************************************************************************







