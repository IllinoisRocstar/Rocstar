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
! Purpose: extrapolate gradients from the interior domain to the dummy points  
!          for symmetry boundary condition
!
! Description: dummy gradients mirrored from interior
!
! Input: region  = data of current region
!        lbound  = patch boundary index
!        idir,jdir,kdir = patch direction
!        inode,jnode,knode = patch nodes identifiers
!        ibeg,iend,jbeg,jend,kbeg,kend = patch cell indices
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
! $Id: RFLO_CalcGradDummySymm.F90,v 1.3 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CalcGradDummySymm( region,lbound, &
                                   idir  ,jdir  ,kdir  , &
                                   inode ,jnode ,knode , &
                                   ibeg  ,iend  ,jbeg  ,jend  ,kbeg  ,kend  , &
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
  INTEGER        :: lbound, idir, jdir, kdir, inode, jnode, knode
  INTEGER        :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER        :: indBeg, indEnd, jndBeg, jndEnd, kndBeg, kndEnd
  INTEGER        :: iBegV, iEndV, iBegG, iEndG
  REAL(RFREAL), POINTER :: gradi(:,:), gradj(:,:), gradk(:,:)

! ... loop variables
  INTEGER :: i, j, k, l, lx,ly,lz, idum, jdum, kdum

! ... local variables
  INTEGER :: iLev, iNOff, ijNOff
  INTEGER :: ijkNB, ijkNI, ijkNJ, ijkNK, ijkND, ijkNR, ijkNR1, ijkNR2, ijkNR3
  INTEGER :: im1jkNB, ijm1kNB, ijkm1NB
  INTEGER :: indxi(indBeg-1:indEnd), indxj(jndBeg-1:jndEnd)
  INTEGER :: indxk(kndBeg-1:kndEnd)
  INTEGER :: iG(3), jG(3), kG(3), iR(3), jR(3), kR(3)
  INTEGER :: iDumB(3), jDumB(3), kDumB(3), iDumE(3), jDumE(3), kDumE(3)
  INTEGER :: nvar

  REAL(RFREAL)  :: dabsi, dabsj, dabsk
  REAL(RFREAL)  :: bgradi(iBegG:iEndG)
  REAL(RFREAL)  :: bgradj(iBegG:iEndG)
  REAL(RFREAL)  :: bgradk(iBegG:iEndG)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CalcGradDummySymm',&
  'RFLO_CalcGradDummySymm.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  nvar = iEndV - iBegV + 1

  dabsi = REAL(ABS(idir))
  dabsj = REAL(ABS(jdir))
  dabsk = REAL(ABS(kdir))

  IF (lbound==1 .OR. lbound==3 .OR. lbound==5) THEN
    iDumB(:) = -idir*region%nDumCells
    jDumB(:) = -jdir*region%nDumCells
    kDumB(:) = -kdir*region%nDumCells
    iDumE(:) = -idir
    jDumE(:) = -jdir
    kDumE(:) = -kdir
  ELSE
    iDumB(:) = 0
    jDumB(:) = 0
    kDumB(:) = 0
    iDumE(:) = -idir*region%nDumCells
    jDumE(:) = -jdir*region%nDumCells
    kDumE(:) = -kdir*region%nDumCells
  ENDIF

  iG(:)=-abs(idir); jG(:)=-abs(jdir); kG(:)=-abs(kdir)
  iR(:)=0         ; jR(:)=0         ; kR(:)=0
  IF (lbound == 1) THEN
    iG(1)=0
    iR(2)=-1; iR(3)=-1
  ELSE IF (lbound == 3) THEN
    jG(2)=0
    jR(1)=-1; jR(3)=-1
  ELSE IF (lbound == 5) THEN
    kG(3)=0
    kR(1)=-1; kR(2)=-1
  ELSE IF (lbound == 2) THEN
    iG(1)=0
    iDumB(1) = iDumB(1)-idir
  ELSE IF (lbound == 4) THEN
    jG(2)=0
    jDumB(2) = jDumB(2)-jdir
  ELSE IF (lbound == 6) THEN
    kG(3)=0
    kDumB(3) = kDumB(3)-kdir
  ENDIF

! specify actual gradients at symmetry patch and denoted as reference gradients

  DO i=indBeg,indEnd
    indxi(i)=i
  Enddo
  DO j=jndBeg,jndEnd
    indxj(j)=j
  Enddo
  DO k=kndBeg,kndEnd
    indxk(k)=k
  Enddo
  indxi(indBeg-1)=indBeg
  indxj(jndBeg-1)=jndBeg
  indxk(kndBeg-1)=kndBeg

  DO k=kbeg+knode,kend+knode
    DO j=jbeg+jnode,jend+jnode
      DO i=ibeg+inode,iend+inode

        ijkNB   = IndIJK(i  ,j  ,k         ,iNOff,ijNOff)  ! ijk bnd nodes
        im1jkNB = IndIJK(indxi(i-1),j  ,k  ,iNOff,ijNOff)  ! im1jk bnd nodes
        ijm1kNB = IndIJK(i  ,indxj(j-1),k  ,iNOff,ijNOff)  ! ijm1k bnd nodes
        ijkm1NB = IndIJK(i  ,j  ,indxk(k-1),iNOff,ijNOff)  ! ijkm1 bnd nodes

        ijkNR1 = IndIJK(i+iR(1),j+jR(1),k+kR(1),iNOff,ijNOff)  ! reference nodes
        ijkNR2 = IndIJK(i+iR(2),j+jR(2),k+kR(2),iNOff,ijNOff)  ! reference nodes
        ijkNR3 = IndIJK(i+iR(3),j+jR(3),k+kR(3),iNOff,ijNOff)  ! reference nodes
        DO l=iBegG,iEndG
          gradi(l,ijkNR1)=dabsi*gradi(l,ijkNR1)+ &
                          dabsj*0.5_RFREAL*(gradj(l,ijkNB)+gradj(l,im1jkNB))+ &
                          dabsk*0.5_RFREAL*(gradk(l,ijkNB)+gradk(l,im1jkNB))

          gradj(l,ijkNR2)=dabsj*gradj(l,ijkNR2)+ &
                          dabsi*0.5_RFREAL*(gradi(l,ijkNB)+gradi(l,ijm1kNB))+ &
                          dabsk*0.5_RFREAL*(gradk(l,ijkNB)+gradk(l,ijm1kNB))

          gradk(l,ijkNR3)=dabsk*gradk(l,ijkNR3)+ &
                          dabsi*0.5_RFREAL*(gradi(l,ijkNB)+gradi(l,ijkm1NB))+ &
                          dabsj*0.5_RFREAL*(gradj(l,ijkNB)+gradj(l,ijkm1NB))
        ENDDO ! l
      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! fill in the last reference line of gradients ---------------------------------

  IF (lbound==1 .OR. lbound==2) THEN

    DO j=jbeg,jend
      ijkNB = IndIJK(indBeg,j,kend,iNOff,ijNOff)  ! boundary nodes
      ijkNR = IndIJK(indBeg+iR(3),j,kndEnd,iNOff,ijNOff)  ! reference nodes
      DO l=iBegG,iEndG
        gradk(l,ijkNR)=gradi(l,ijkNB)
      ENDDO
    ENDDO

    DO k=kbeg,kend
      ijkNB = IndIJK(indBeg,jend,k,iNOff,ijNOff)  ! boundary nodes
      ijkNR = IndIJK(indBeg+iR(2),jndEnd,k,iNOff,ijNOff)  ! reference nodes
      DO l=iBegG,iEndG
        gradj(l,ijkNR)=gradi(l,ijkNB)
      ENDDO
    ENDDO

  ELSEIF (lbound==3 .OR. lbound==4) THEN

    DO k=kbeg,kend
      ijkNB = IndIJK(iend,jndBeg,k,iNOff,ijNOff)  ! boundary nodes
      ijkNR = IndIJK(indEnd,jndBeg+jR(1),k,iNOff,ijNOff)  ! reference nodes
      DO l=iBegG,iEndG
        gradi(l,ijkNR)=gradj(l,ijkNB)
      ENDDO
    ENDDO

    DO i=ibeg,iend
      ijkNB = IndIJK(i,jndBeg,kend,iNOff,ijNOff)  ! boundary nodes
      ijkNR = IndIJK(i,jndBeg+jR(3),kndEnd,iNOff,ijNOff)  ! reference nodes
      DO l=iBegG,iEndG
        gradk(l,ijkNR)=gradj(l,ijkNB)
      ENDDO
    ENDDO

  ELSEIF (lbound==5 .OR. lbound==6) THEN

    DO j=jbeg,jend
      ijkNB = IndIJK(iend,j,kndBeg,iNOff,ijNOff)  ! boundary nodes
      ijkNR = IndIJK(indEnd,j,kndBeg+kR(1),iNOff,ijNOff)  ! reference nodes
      DO l=iBegG,iEndG
        gradi(l,ijkNR)=gradk(l,ijkNB)
      ENDDO
    END DO

    DO i=ibeg,iend
      ijkNB = IndIJK(i,jend,kndBeg,iNOff,ijNOff)  ! boundary nodes
      ijkNR = IndIJK(i,jndEnd,kndBeg+kR(2),iNOff,ijNOff)  ! reference nodes
      DO l=iBegG,iEndG
        gradj(l,ijkNR)=gradk(l,ijkNB)
      ENDDO
    END DO

  ENDIF

! specify gradients at dummy points by extrapolation (mirroring) ---------------
! 2D loop over patch nodes 

  DO k=kndBeg,kndEnd
    DO j=jndBeg,jndEnd
      DO i=indBeg,indEnd

! ----- 1D loop in direction normal to patch to define dummy face-gradi,j, and k 

        ijkNR1 = IndIJK(i+iR(1),j+jR(1),k+kR(1),iNOff,ijNOff)  ! reference nodes
        ijkNR2 = IndIJK(i+iR(2),j+jR(2),k+kR(2),iNOff,ijNOff)  ! reference nodes
        ijkNR3 = IndIJK(i+iR(3),j+jR(3),k+kR(3),iNOff,ijNOff)  ! reference nodes
        DO l=iBegG,iEndG
          bgradi(l)=gradi(l,ijkNR1)
          bgradj(l)=gradj(l,ijkNR2)
          bgradk(l)=gradk(l,ijkNR3)
        ENDDO

        DO idum=iDumB(1),iDumE(1)
          DO jdum=jDumB(1),jDumE(1)
            DO kdum=kDumB(1),kDumE(1) 
              ijkND = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff)  ! dummy points
              ijkNI = IndIJK(i+iG(1)-idum,j+jG(1)-jdum,k+kG(1)-kdum,iNOff,ijNOff)  ! interior
              DO l=iBegG,iEndG
                gradi(l,ijkND)=2.0_RFREAL*bgradi(l)-gradi(l,ijkNI)
              ENDDO
            ENDDO
          ENDDO  
        ENDDO    

        DO idum=iDumB(2),iDumE(2)
          DO jdum=jDumB(2),jDumE(2)
            DO kdum=kDumB(2),kDumE(2) 
              ijkND = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff) ! dummy points
              ijkNJ = IndIJK(i+iG(2)-idum,j+jG(2)-jdum,k+kG(2)-kdum,iNOff,ijNOff)  ! interior
              DO l=iBegG,iEndG
                gradj(l,ijkND)=2.0_RFREAL*bgradj(l)-gradj(l,ijkNJ)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        DO idum=iDumB(3),iDumE(3)
          DO jdum=jDumB(3),jDumE(3)
            DO kdum=kDumB(3),kDumE(3) 
              ijkND = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff) ! dummy points
              ijkNK = IndIJK(i+iG(3)-idum,j+jG(3)-jdum,k+kG(3)-kdum,iNOff,ijNOff)  ! interior
              DO l=iBegG,iEndG
                gradk(l,ijkND)=2.0_RFREAL*bgradk(l)-gradk(l,ijkNK)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_CalcGradDummySymm

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: RFLO_CalcGradDummySymm.F90,v $
! Revision 1.3  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.9  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.5  2002/10/14 23:47:39  wasistho
! Minor tunning to get speedup
!
! Revision 1.4  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.3  2002/10/03 02:54:49  wasistho
! modification through indxi,j,k
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
! Revision 1.2  2002/07/19 23:43:51  wasistho
! made compliant with CODING RULE
!
!******************************************************************************







