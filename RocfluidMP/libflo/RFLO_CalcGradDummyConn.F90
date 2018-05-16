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
! Purpose: determine gradients at the outerst dummy layers of connecting 
!          boundaries
!
! Description: copy defined gradients at the dummy points to the outerst dummy 
!              layers of connecting boundaries
!
! Input: region  = data of current region
!        lbound  = patch boundary index
!        idir,jdir,kdir = patch direction
!        indBeg,indEnd,jndBeg,jndEnd,kndBeg,kndEnd = indices of patch nodes
!        iBegV, iEndV = begin and end var index
!        iBegG, iEndG = begin and end gradient index
!  
! Output: gradi, gradj, gradk = gradients at dummy faces of the current patch
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_CalcGradDummyConn.F90,v 1.3 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CalcGradDummyConn( region,lbound, &
                                   idir  ,jdir  ,kdir  , &
                                   indBeg,indEnd,jndBeg,jndEnd,kndBeg,kndEnd, &
                                   iBegV ,iEndV ,iBegG ,iEndG , &
                                   gradi ,gradj ,gradk )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: lbound,idir,jdir,kdir
  INTEGER        :: indBeg, indEnd, jndBeg, jndEnd, kndBeg, kndEnd
  INTEGER        :: iBegV, iEndV, iBegG, iEndG
  REAL(RFREAL), POINTER :: gradi(:,:), gradj(:,:), gradk(:,:)

! ... loop variables
  INTEGER :: i, j, k, idum, jdum, kdum

! ... local variables
  INTEGER :: iLev, iNOff, ijNOff
  INTEGER :: ijkNI, ijkNJ, ijkNK, ijkND
  INTEGER :: iG(3), jG(3), kG(3)
  INTEGER :: iDumB(3), jDumB(3), kDumB(3), iDumE(3), jDumE(3), kDumE(3)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CalcGradDummyConn',&
       'RFLO_CalcGradDummyConn.F90' )

! get dimensions -------------------------------------------------------------

  iLev = region%currLevel
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  iDumB(:) = -idir*region%nDumCells
  jDumB(:) = -jdir*region%nDumCells
  kDumB(:) = -kdir*region%nDumCells
  iDumE(:) = iDumB(:)
  jDumE(:) = jDumB(:)
  kDumE(:) = kDumB(:)

  iG(:) = -idir*(region%nDumCells-1)
  jG(:) = -jdir*(region%nDumCells-1)
  kG(:) = -kdir*(region%nDumCells-1)
  IF (lbound == 2) THEN
    iDumB(2) = iDumB(2)+idir; iDumB(3) = iDumB(3)+idir
    iG(2) = iG(2)+idir; iG(3) = iG(3)+idir
  ELSE IF (lbound == 4) THEN
    jDumB(1) = jDumB(1)+jdir; jDumB(3) = jDumB(3)+jdir
    jG(1) = jG(1)+jdir; jG(3) = jG(3)+jdir
  ELSE IF (lbound == 6) THEN
    kDumB(1) = kDumB(1)+kdir; kDumB(2) = kDumB(2)+kdir
    kG(1) = kG(1)+kdir; kG(2) = kG(2)+kdir
  ENDIF

! 2D loop over patch nodes ----------------------------------------------------

  DO k=kndBeg,kndEnd
    DO j=jndBeg,jndEnd
      DO i=indBeg,indEnd

        ijkNI = IndIJK(i+iG(1),j+jG(1),k+kG(1),iNOff,ijNOff)  ! i reference
        ijkNJ = IndIJK(i+iG(2),j+jG(2),k+kG(2),iNOff,ijNOff)  ! j reference
        ijkNK = IndIJK(i+iG(3),j+jG(3),k+kG(3),iNOff,ijNOff)  ! k reference

! ----- 1D loop in direction normal to patch to define grads at dummy faces

        DO idum=iDumB(1),iDumE(1)
          DO jdum=jDumB(1),jDumE(1)
            DO kdum=kDumB(1),kDumE(1) 
              ijkND = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff)
              gradi(:,ijkND) = gradi(:,ijkNI)
            ENDDO
          ENDDO
        ENDDO

        DO idum=iDumB(2),iDumE(2)
          DO jdum=jDumB(2),jDumE(2)
            DO kdum=kDumB(2),kDumE(2) 
              ijkND  = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff)
              gradj(:,ijkND) = gradj(:,ijkNJ)
            ENDDO
          ENDDO
        ENDDO

        DO idum=iDumB(3),iDumE(3)
          DO jdum=jDumB(3),jDumE(3)
            DO kdum=kDumB(3),kDumE(3) 
              ijkND  = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff)
              gradk(:,ijkND) = gradk(:,ijkNK)
            ENDDO
          ENDDO
        ENDDO

      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_CalcGradDummyConn

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: RFLO_CalcGradDummyConn.F90,v $
! Revision 1.3  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.6  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.3  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.2  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/09/02 22:58:54  wasistho
! RFLO grad routines migrated from rocflo to libflo
!
! Revision 1.3  2002/07/22 22:59:11  jblazek
! Some more clean up.
!
! Revision 1.2  2002/07/19 23:43:21  wasistho
! made compliant with CODING RULE
!
!******************************************************************************







