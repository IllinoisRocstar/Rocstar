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
! Purpose: Compute filter coefficients of face-face filtering for the case 
!          where the filtering direction is the same as the face direction
!          being treated.
!
! Description: The coefficients are obtained in the ijk direction.
!              FFLo if for low delta (1 or 2 grid-spacing) and FFHi is for 
!              high delta (2 or 4 grid-spacing).
!
! Input: global = global data only used by the register function
!        ijk    = ijk-face is being treated
!        i,j,kbeg,i,j,kend = cell index boundaries
!        minIdx, maxIdx    = begin and end index of temporary 1D array
!        segId             = face-width ID 
!        ..NOff            = cells and nodes offset
!        ds     = 1D temporary array
!        segm   = two components face-width
!
! Output: ffCofA = filter coefficients for the smaller filter-width
!         ffCofB = filter coefficients for the larger filter-width
!
! Notes: Mother routine = LesGenCoFF.
!
!******************************************************************************
!
! $Id: TURB_floLesGenCoFFUtil.F90,v 1.4 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloLesGenCoFFLo( global,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_global), POINTER :: global
  INTEGER               :: ijk,ibeg,iend,jbeg,jend,kbeg,kend,minIdx,maxIdx
  INTEGER               :: segId,iNOff,ijNOff
  REAL(RFREAL)          :: ds(minIdx:maxIdx)
  REAL(RFREAL), POINTER :: segm(:,:),ffCofA(:,:),ffCofB(:,:)

! ... loop variables
  INTEGER :: i, j, k, ijkN, ijkN1

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER           :: indi(ibeg:iend+1),indj(jbeg:jend+1),indk(kbeg:kend+1)
  REAL(RFREAL)      :: rWidth1, fac1, fac3

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floLesGenCoFFUtil.F90,v $'

  CALL RegisterFunction( global,'TURB_FloLesGenCoFFLo',&
  'TURB_floLesGenCoFFUtil.F90' )

! get indices ---------------------------------------------------------------

  DO i = ibeg,iend
    indi(i)=i
  ENDDO
  DO j = jbeg,jend
    indj(j)=j
  ENDDO
  DO k = kbeg,kend
    indk(k)=k
  ENDDO
  indi(iend+1)=iend
  indj(jend+1)=jend
  indk(kend+1)=kend

! start computations 

  IF (ijk==DIRI) THEN
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkN  = IndIJK(i     ,j     ,k        ,iNOff,ijNOff)
          ijkN1 = IndIJK(i     ,j     ,indk(k+1),iNOff,ijNOff)
          ds(i) = 0.5_RFREAL*(segm(segId,ijkN)+segm(segId,ijkN1))
        ENDDO
        DO i=1,2
          ds(ibeg-i) = ds(ibeg)
          ds(iend+i) = ds(iend)
        ENDDO
        DO i=ibeg,iend
          ijkN   = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          rWidth1= 1._RFREAL/(ds(i-1)+ds(i))
          fac1   = ds(i-1)*rWidth1
          fac3   = ds(i)*rWidth1
          ffCofA(1,ijkN) = 0.25_RFREAL*fac1
          ffCofA(2,ijkN) = 0.75_RFREAL
          ffCofA(3,ijkN) = 0.25_RFREAL*fac3
          ffCofB(1,ijkN) = 0.5_RFREAL*fac1
          ffCofB(2,ijkN) = 0.5_RFREAL
          ffCofB(3,ijkN) = 0.5_RFREAL*fac3
        ENDDO
      ENDDO 
    ENDDO     

  ELSEIF (ijk==DIRJ) THEN
    DO k=kbeg,kend
      DO i=ibeg,iend
        DO j=jbeg,jend
          ijkN  = IndIJK(i        ,j     ,k     ,iNOff,ijNOff)
          ijkN1 = IndIJK(indi(i+1),j     ,k     ,iNOff,ijNOff)
          ds(j) = 0.5_RFREAL*(segm(segId,ijkN)+segm(segId,ijkN1))
        ENDDO
        DO j=1,2
          ds(jbeg-j) = ds(jbeg)
          ds(jend+j) = ds(jend)
        ENDDO
        DO j=jbeg,jend
          ijkN   = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          rWidth1= 1._RFREAL/(ds(j-1)+ds(j))
          fac1   = ds(j-1)*rWidth1
          fac3   = ds(j)*rWidth1
          ffCofA(1,ijkN) = 0.25_RFREAL*fac1
          ffCofA(2,ijkN) = 0.75_RFREAL
          ffCofA(3,ijkN) = 0.25_RFREAL*fac3
          ffCofB(1,ijkN) = 0.5_RFREAL*fac1
          ffCofB(2,ijkN) = 0.5_RFREAL
          ffCofB(3,ijkN) = 0.5_RFREAL*fac3
        ENDDO
      ENDDO 
    ENDDO     

  ELSEIF (ijk==DIRK) THEN
    DO j=jbeg,jend
      DO i=ibeg,iend
        DO k=kbeg,kend
          ijkN  = IndIJK(i     ,j        ,k     ,iNOff,ijNOff)
          ijkN1 = IndIJK(i     ,indj(j+1),k     ,iNOff,ijNOff)
          ds(k) = 0.5_RFREAL*(segm(segId,ijkN)+segm(segId,ijkN1))
        ENDDO
        DO k=1,2
          ds(kbeg-k) = ds(kbeg)
          ds(kend+k) = ds(kend)
        ENDDO
        DO k=kbeg,kend
          ijkN   = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          rWidth1= 1._RFREAL/(ds(k-1)+ds(k))
          fac1   = ds(k-1)*rWidth1
          fac3   = ds(k)*rWidth1
          ffCofA(1,ijkN) = 0.25_RFREAL*fac1
          ffCofA(2,ijkN) = 0.75_RFREAL
          ffCofA(3,ijkN) = 0.25_RFREAL*fac3
          ffCofB(1,ijkN) = 0.5_RFREAL*fac1
          ffCofB(2,ijkN) = 0.5_RFREAL
          ffCofB(3,ijkN) = 0.5_RFREAL*fac3
        ENDDO
      ENDDO 
    ENDDO     
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloLesGenCoFFLo

!#############################################################################

SUBROUTINE TURB_FloLesGenCoFFHi( global,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
                       minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm,ffCofA,ffCofB )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_global), POINTER :: global
  INTEGER               :: ijk,ibeg,iend,jbeg,jend,kbeg,kend,minIdx,maxIdx
  INTEGER               :: segId,iNOff,ijNOff
  REAL(RFREAL)          :: ds(minIdx:maxIdx)
  REAL(RFREAL), POINTER :: segm(:,:),ffCofA(:,:),ffCofB(:,:)

! ... loop variables
  INTEGER :: i, j, k, ijkN, ijkN1

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER           :: indi(ibeg:iend+1),indj(jbeg:jend+1),indk(kbeg:kend+1)
  REAL(RFREAL)      :: rWidth1, fac1, fac2, fac3, fac4, fac5

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floLesGenCoFFUtil.F90,v $'

  CALL RegisterFunction( global,'TURB_FloLesGenCoFFHi',&
  'TURB_floLesGenCoFFUtil.F90' )

! get indices ---------------------------------------------------------------

  DO i = ibeg,iend
    indi(i)=i
  ENDDO
  DO j = jbeg,jend
    indj(j)=j
  ENDDO
  DO k = kbeg,kend
    indk(k)=k
  ENDDO
  indi(iend+1)=iend
  indj(jend+1)=jend
  indk(kend+1)=kend

! start computations

  IF (ijk==DIRI) THEN
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkN  = IndIJK(i     ,j     ,k        ,iNOff,ijNOff)
          ijkN1 = IndIJK(i     ,j     ,indk(k+1),iNOff,ijNOff)
          ds(i) = 0.5_RFREAL*(segm(segId,ijkN)+segm(segId,ijkN1))
        ENDDO
        DO i=1,2
          ds(ibeg-i) = ds(ibeg)
          ds(iend+i) = ds(iend)
        ENDDO
        DO i=ibeg,iend
          ijkN   = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          rWidth1= 1._RFREAL/(ds(i-1)+ds(i))
          fac1   = ds(i-1)*rWidth1
          fac3   = ds(i)*rWidth1
          ffCofA(1,ijkN) = 0.5_RFREAL*fac1
          ffCofA(2,ijkN) = 0.5_RFREAL
          ffCofA(3,ijkN) = 0.5_RFREAL*fac3

          rWidth1= 1._RFREAL/(ds(i-2)+ds(i-1)+ds(i)+ds(i+1))
          fac1   = ds(i-2)*rWidth1
          fac2   = (ds(i-2)+ds(i-1))*rWidth1
          fac3   = (ds(i-1)+ds(i))*rWidth1
          fac4   = (ds(i)+ds(i+1))*rWidth1
          fac5   = ds(i+1)*rWidth1
          ffCofB(1,ijkN) = 0.5_RFREAL*fac1
          ffCofB(2,ijkN) = 0.5_RFREAL*fac2
          ffCofB(3,ijkN) = 0.5_RFREAL*fac3
          ffCofB(4,ijkN) = 0.5_RFREAL*fac4
          ffCofB(5,ijkN) = 0.5_RFREAL*fac5
        ENDDO
      ENDDO 
    ENDDO     

  ELSEIF (ijk==DIRJ) THEN
    DO k=kbeg,kend
      DO i=ibeg,iend
        DO j=jbeg,jend
          ijkN  = IndIJK(i        ,j     ,k     ,iNOff,ijNOff)
          ijkN1 = IndIJK(indi(i+1),j     ,k     ,iNOff,ijNOff)
          ds(j) = 0.5_RFREAL*(segm(segId,ijkN)+segm(segId,ijkN1))
        ENDDO
        DO j=1,2
          ds(jbeg-j) = ds(jbeg)
          ds(jend+j) = ds(jend)
        ENDDO
        DO j=jbeg,jend
          ijkN   = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          rWidth1= 1._RFREAL/(ds(j-1)+ds(j))
          fac1   = ds(j-1)*rWidth1
          fac3   = ds(j)*rWidth1
          ffCofA(1,ijkN) = 0.5_RFREAL*fac1
          ffCofA(2,ijkN) = 0.5_RFREAL
          ffCofA(3,ijkN) = 0.5_RFREAL*fac3

          rWidth1= 1._RFREAL/(ds(j-2)+ds(j-1)+ds(j)+ds(j+1))
          fac1   = ds(j-2)*rWidth1
          fac2   = (ds(j-2)+ds(j-1))*rWidth1
          fac3   = (ds(j-1)+ds(j))*rWidth1
          fac4   = (ds(j)+ds(j+1))*rWidth1
          fac5   = ds(j+1)*rWidth1
          ffCofB(1,ijkN) = 0.5_RFREAL*fac1
          ffCofB(2,ijkN) = 0.5_RFREAL*fac2
          ffCofB(3,ijkN) = 0.5_RFREAL*fac3
          ffCofB(4,ijkN) = 0.5_RFREAL*fac4
          ffCofB(5,ijkN) = 0.5_RFREAL*fac5
        ENDDO
      ENDDO 
    ENDDO     

  ELSEIF (ijk==DIRK) THEN
    DO j=jbeg,jend
      DO i=ibeg,iend
        DO k=kbeg,kend
          ijkN  = IndIJK(i     ,j        ,k     ,iNOff,ijNOff)
          ijkN1 = IndIJK(i     ,indj(j+1),k     ,iNOff,ijNOff)
          ds(k) = 0.5_RFREAL*(segm(segId,ijkN)+segm(segId,ijkN1))
        ENDDO
        DO k=1,2
          ds(kbeg-k) = ds(kbeg)
          ds(kend+k) = ds(kend)
        ENDDO
        DO k=kbeg,kend
          ijkN   = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          rWidth1= 1._RFREAL/(ds(k-1)+ds(k))
          fac1   = ds(k-1)*rWidth1
          fac3   = ds(k)*rWidth1
          ffCofA(1,ijkN) = 0.5_RFREAL*fac1
          ffCofA(2,ijkN) = 0.5_RFREAL
          ffCofA(3,ijkN) = 0.5_RFREAL*fac3

          rWidth1= 1._RFREAL/(ds(k-2)+ds(k-1)+ds(k)+ds(k+1))
          fac1   = ds(k-2)*rWidth1
          fac2   = (ds(k-2)+ds(k-1))*rWidth1
          fac3   = (ds(k-1)+ds(k))*rWidth1
          fac4   = (ds(k)+ds(k+1))*rWidth1
          fac5   = ds(k+1)*rWidth1
          ffCofB(1,ijkN) = 0.5_RFREAL*fac1
          ffCofB(2,ijkN) = 0.5_RFREAL*fac2
          ffCofB(3,ijkN) = 0.5_RFREAL*fac3
          ffCofB(4,ijkN) = 0.5_RFREAL*fac4
          ffCofB(5,ijkN) = 0.5_RFREAL*fac5
        ENDDO
      ENDDO 
    ENDDO     
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloLesGenCoFFHi

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floLesGenCoFFUtil.F90,v $
! Revision 1.4  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.2  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!******************************************************************************








