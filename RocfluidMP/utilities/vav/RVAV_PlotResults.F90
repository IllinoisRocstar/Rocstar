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
! Purpose: Plots results of Stream 1 and Stream 2 Data.
!
! Description: none.
!
! Input: Arrays evS1 and evS2 and their corresponding extents
!
! Output: File 
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RVAV_PlotResults.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RVAV_PlotResults( global,iCompare,iNodes,jNodes,kNodes )

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY : t_global
  USE RVAV_ModGlobal
  USE RVAV_ModParameters
  IMPLICIT NONE

! ... parameter
  TYPE(t_global), POINTER :: global

  INTEGER :: iCompare, iNodes, jNodes, kNodes

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  CHARACTER(CHRLEN) ::  msg

  REAL(RFREAL) :: l2Norm, l8Norm

!******************************************************************************

  CALL RegisterFunction( global, 'RVAV_PlotResults',&
  'RVAV_PlotResults.F90' )

! initialize variables

  WRITE(200,'(A,I2)') 'ICompare = ',iCompare

  DO k = 1, kNodes
    DO j = 1, jNodes
      DO i = 1, iNodes
        WRITE(200,1000) i,j,k,globalRVAV%evS1(i,j,k),globalRVAV%evS2(i,j,k)
      ENDDO ! i
    ENDDO ! j
  ENDDO ! k

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1000 FORMAT(3(3X,3I5),2(3X,1PE15.7))

END SUBROUTINE RVAV_PlotResults

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_PlotResults.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:43:25  fnajjar
! Initial revision after changing case
!
! Revision 1.3  2003/05/15 02:57:08  jblazek
! Inlined index function.
!
! Revision 1.2  2002/09/10 00:01:45  f-najjar
! Variable global moved into regions()
!
! Revision 1.1.1.1  2002/06/03 21:41:29  f-najjar
! Initial Import of RocVaV
!
!******************************************************************************







