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
! Purpose: Compute filter coefficients of cell-cell filtering.
!
! Description: The coefficients are obtained in the ijk direction.
!              CCLo if for low delta (1 or 2 grid-spacing) and CCHi is for 
!              high delta (2 or 4 grid-spacing).
!
! Input: global = global data only used by the register function
!        ijk    = ijk-face is being treated
!        i,j,kbeg,i,j,kend = cell index boundaries
!        minIdx, maxIdx    = begin and end index of temporary 1D array
!        segId             = face-width ID 
!        ..COff,..NOff     = cells and nodes offset
!        ds     = 1D temporary array
!        segm   = two components face-width
!
! Output: ccCofA = filter coefficients for the smaller filter-width
!         ccCofB = filter coefficients for the larger filter-width
!
! Notes: Mother routine = LesGenCoCC.
!
!******************************************************************************
!
! $Id: TURB_floLesGenCoCCUtil.F90,v 1.4 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloLesGenCoCCLo( global,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
                              minIdx,maxIdx,segId,iCOff,ijCOff,iNOff,ijNOff, &
                              ds,segm,ccCofA,ccCofB )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_global), POINTER :: global
  INTEGER               :: ijk,ibeg,iend,jbeg,jend,kbeg,kend,minIdx,maxIdx
  INTEGER               :: segId,iCOff,ijCOff,iNOff,ijNOff
  REAL(RFREAL)          :: ds(minIdx:maxIdx)
  REAL(RFREAL), POINTER :: segm(:,:),ccCofA(:,:),ccCofB(:,:)

! ... loop variables
  INTEGER :: i, j, k, ijkC, ijkN, ijkN1

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  REAL(RFREAL)      :: rWidth1, rWidth2, fac1, fac2, fac3

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floLesGenCoCCUtil.F90,v $'

  CALL RegisterFunction( global,'TURB_FloLesGenCoCCLo',&
  'TURB_floLesGenCoCCUtil.F90' )

! start computations --------------------------------------------------------

  IF (ijk==DIRI) THEN
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          ijkN1 = IndIJK(i     ,j     ,k+1   ,iNOff,ijNOff)
          ds(i) = 0.5_RFREAL*(segm(segId,ijkN)+segm(segId,ijkN1))
        ENDDO
        DO i=1,2
          ds(ibeg-i) = ds(ibeg)
          ds(iend+i) = ds(iend)
        ENDDO
        DO i=ibeg,iend
          ijkC   = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          rWidth1= 1._RFREAL/(ds(i-1)+ds(i))
          rWidth2= 1._RFREAL/(ds(i+1)+ds(i))
          fac1   = ds(i)*rWidth1
          fac2   = ds(i-1)*rWidth1+ds(i+1)*rWidth2+2._RFREAL
          fac3   = ds(i)*rWidth2
          ccCofA(1,ijkC) = 0.25_RFREAL*fac1
          ccCofA(2,ijkC) = 0.25_RFREAL*fac2
          ccCofA(3,ijkC) = 0.25_RFREAL*fac3
  
          rWidth1= 1._RFREAL/(ds(i-1)+2._RFREAL*ds(i)+ds(i+1))
          fac1   = (ds(i-1)+ds(i))*rWidth1
          fac3   = (ds(i+1)+ds(i))*rWidth1
          ccCofB(1,ijkC) = 0.5_RFREAL*fac1
          ccCofB(2,ijkC) = 0.5_RFREAL
          ccCofB(3,ijkC) = 0.5_RFREAL*fac3
        ENDDO
      ENDDO 
    ENDDO     

  ELSEIF (ijk==DIRJ) THEN
    DO k=kbeg,kend
      DO i=ibeg,iend
        DO j=jbeg,jend
          ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          ijkN1 = IndIJK(i+1   ,j     ,k     ,iNOff,ijNOff)
          ds(j) = 0.5_RFREAL*(segm(segId,ijkN)+segm(segId,ijkN1))
        ENDDO
        DO j=1,2
          ds(jbeg-j) = ds(jbeg)
          ds(jend+j) = ds(jend)
        ENDDO
        DO j=jbeg,jend
          ijkC   = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          rWidth1= 1._RFREAL/(ds(j-1)+ds(j))
          rWidth2= 1._RFREAL/(ds(j+1)+ds(j))
          fac1   = ds(j)*rWidth1
          fac2   = ds(j-1)*rWidth1+ds(j+1)*rWidth2+2._RFREAL
          fac3   = ds(j)*rWidth2
          ccCofA(1,ijkC) = 0.25_RFREAL*fac1
          ccCofA(2,ijkC) = 0.25_RFREAL*fac2
          ccCofA(3,ijkC) = 0.25_RFREAL*fac3
  
          rWidth1= 1._RFREAL/(ds(j-1)+2._RFREAL*ds(j)+ds(j+1))
          fac1   = (ds(j-1)+ds(j))*rWidth1
          fac3   = (ds(j+1)+ds(j))*rWidth1
          ccCofB(1,ijkC) = 0.5_RFREAL*fac1
          ccCofB(2,ijkC) = 0.5_RFREAL
          ccCofB(3,ijkC) = 0.5_RFREAL*fac3
        ENDDO
      ENDDO 
    ENDDO     

  ELSEIF (ijk==DIRK) THEN
    DO j=jbeg,jend
      DO i=ibeg,iend
        DO k=kbeg,kend
          ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          ijkN1 = IndIJK(i     ,j+1   ,k     ,iNOff,ijNOff)
          ds(k) = 0.5_RFREAL*(segm(segId,ijkN)+segm(segId,ijkN1))
        ENDDO
        DO k=1,2
          ds(kbeg-k) = ds(kbeg)
          ds(kend+k) = ds(kend)
        ENDDO
        DO k=kbeg,kend
          ijkC   = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          rWidth1= 1._RFREAL/(ds(k-1)+ds(k))
          rWidth2= 1._RFREAL/(ds(k+1)+ds(k))
          fac1   = ds(k)*rWidth1
          fac2   = ds(k-1)*rWidth1+ds(k+1)*rWidth2+2._RFREAL
          fac3   = ds(k)*rWidth2
          ccCofA(1,ijkC) = 0.25_RFREAL*fac1
          ccCofA(2,ijkC) = 0.25_RFREAL*fac2
          ccCofA(3,ijkC) = 0.25_RFREAL*fac3
  
          rWidth1= 1._RFREAL/(ds(k-1)+2._RFREAL*ds(k)+ds(k+1))
          fac1   = (ds(k-1)+ds(k))*rWidth1
          fac3   = (ds(k+1)+ds(k))*rWidth1
          ccCofB(1,ijkC) = 0.5_RFREAL*fac1
          ccCofB(2,ijkC) = 0.5_RFREAL
          ccCofB(3,ijkC) = 0.5_RFREAL*fac3
        ENDDO
      ENDDO 
    ENDDO     
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloLesGenCoCCLo

!#############################################################################

SUBROUTINE TURB_FloLesGenCoCCHi( global,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
                              minIdx,maxIdx,segId,iCOff,ijCOff,iNOff,ijNOff, &
                              ds,segm,ccCofA,ccCofB )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_global), POINTER :: global
  INTEGER               :: ijk,ibeg,iend,jbeg,jend,kbeg,kend,minIdx,maxIdx
  INTEGER               :: segId,iCOff,ijCOff,iNOff,ijNOff
  REAL(RFREAL)          :: ds(minIdx:maxIdx)
  REAL(RFREAL), POINTER :: segm(:,:),ccCofA(:,:),ccCofB(:,:)

! ... loop variables
  INTEGER :: i, j, k, ijkC, ijkN, ijkN1

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  REAL(RFREAL)      :: rWidth1, rWidth2, fac1, fac2, fac3, fac4, fac5

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floLesGenCoCCUtil.F90,v $'

  CALL RegisterFunction( global,'TURB_FloLesGenCoCCHi',&
  'TURB_floLesGenCoCCUtil.F90' )

! start computations --------------------------------------------------------

  IF (ijk==DIRI) THEN
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          ijkN1 = IndIJK(i     ,j     ,k+1   ,iNOff,ijNOff)
          ds(i) = 0.5_RFREAL*(segm(segId,ijkN)+segm(segId,ijkN1))
        ENDDO
        DO i=1,2
          ds(ibeg-i) = ds(ibeg)
          ds(iend+i) = ds(iend)
        ENDDO
        DO i=ibeg,iend
          ijkC   = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          rWidth1= 1._RFREAL/(ds(i-1)+2._RFREAL*ds(i)+ds(i+1))
          fac1   = (ds(i-1)+ds(i))*rWidth1
          fac3   = (ds(i+1)+ds(i))*rWidth1
          ccCofA(1,ijkC) = 0.5_RFREAL*fac1
          ccCofA(2,ijkC) = 0.5_RFREAL
          ccCofA(3,ijkC) = 0.5_RFREAL*fac3

          rWidth1= 1._RFREAL/ &
                   (ds(i-2)+ds(i+2)+2._RFREAL*(ds(i-1)+ds(i)+ds(i+1)))
          fac1   = (ds(i-2)+ds(i-1))*rWidth1  
          fac2   = (ds(i-2)+2._RFREAL*ds(i-1)+ds(i))*rWidth1  
          fac3   = (ds(i-1)+2._RFREAL*ds(i)+ds(i+1))*rWidth1  
          fac4   = (ds(i)+2._RFREAL*ds(i+1)+ds(i+2))*rWidth1  
          fac5   = (ds(i+1)+ds(i+2))*rWidth1  
          ccCofB(1,ijkC) = 0.5_RFREAL*fac1
          ccCofB(2,ijkC) = 0.5_RFREAL*fac2
          ccCofB(3,ijkC) = 0.5_RFREAL*fac3
          ccCofB(4,ijkC) = 0.5_RFREAL*fac4
          ccCofB(5,ijkC) = 0.5_RFREAL*fac5
        ENDDO
      ENDDO 
    ENDDO     

  ELSEIF (ijk==DIRJ) THEN
    DO k=kbeg,kend
      DO i=ibeg,iend
        DO j=jbeg,jend
          ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          ijkN1 = IndIJK(i+1   ,j     ,k     ,iNOff,ijNOff)
          ds(j) = 0.5_RFREAL*(segm(segId,ijkN)+segm(segId,ijkN1))
        ENDDO
        DO j=1,2
          ds(jbeg-j) = ds(jbeg)
          ds(jend+j) = ds(jend)
        ENDDO
        DO j=jbeg,jend
          ijkC   = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          rWidth1= 1._RFREAL/(ds(j-1)+2._RFREAL*ds(j)+ds(j+1))
          fac1   = (ds(j-1)+ds(j))*rWidth1
          fac3   = (ds(j+1)+ds(j))*rWidth1
          ccCofA(1,ijkC) = 0.5_RFREAL*fac1
          ccCofA(2,ijkC) = 0.5_RFREAL
          ccCofA(3,ijkC) = 0.5_RFREAL*fac3

          rWidth1= 1._RFREAL/ &
                   (ds(j-2)+ds(j+2)+2._RFREAL*(ds(j-1)+ds(j)+ds(j+1)))
          fac1   = (ds(j-2)+ds(j-1))*rWidth1  
          fac2   = (ds(j-2)+2._RFREAL*ds(j-1)+ds(j))*rWidth1  
          fac3   = (ds(j-1)+2._RFREAL*ds(j)+ds(j+1))*rWidth1  
          fac4   = (ds(j)+2._RFREAL*ds(j+1)+ds(j+2))*rWidth1  
          fac5   = (ds(j+1)+ds(j+2))*rWidth1  
          ccCofB(1,ijkC) = 0.5_RFREAL*fac1
          ccCofB(2,ijkC) = 0.5_RFREAL*fac2
          ccCofB(3,ijkC) = 0.5_RFREAL*fac3
          ccCofB(4,ijkC) = 0.5_RFREAL*fac4
          ccCofB(5,ijkC) = 0.5_RFREAL*fac5
        ENDDO
      ENDDO 
    ENDDO     

  ELSEIF (ijk==DIRK) THEN
    DO j=jbeg,jend
      DO i=ibeg,iend
        DO k=kbeg,kend
          ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          ijkN1 = IndIJK(i     ,j+1   ,k     ,iNOff,ijNOff)
          ds(k) = 0.5_RFREAL*(segm(segId,ijkN)+segm(segId,ijkN1))
        ENDDO
        DO k=1,2
          ds(kbeg-k) = ds(kbeg)
          ds(kend+k) = ds(kend)
        ENDDO
        DO k=kbeg,kend
          ijkC   = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          rWidth1= 1._RFREAL/(ds(k-1)+2._RFREAL*ds(k)+ds(k+1))
          fac1   = (ds(k-1)+ds(k))*rWidth1
          fac3   = (ds(k+1)+ds(k))*rWidth1
          ccCofA(1,ijkC) = 0.5_RFREAL*fac1
          ccCofA(2,ijkC) = 0.5_RFREAL
          ccCofA(3,ijkC) = 0.5_RFREAL*fac3

          rWidth1= 1._RFREAL/ &
                   (ds(k-2)+ds(k+2)+2._RFREAL*(ds(k-1)+ds(k)+ds(k+1)))
          fac1   = (ds(k-2)+ds(k-1))*rWidth1  
          fac2   = (ds(k-2)+2._RFREAL*ds(k-1)+ds(k))*rWidth1  
          fac3   = (ds(k-1)+2._RFREAL*ds(k)+ds(k+1))*rWidth1  
          fac4   = (ds(k)+2._RFREAL*ds(k+1)+ds(k+2))*rWidth1  
          fac5   = (ds(k+1)+ds(k+2))*rWidth1  
          ccCofB(1,ijkC) = 0.5_RFREAL*fac1
          ccCofB(2,ijkC) = 0.5_RFREAL*fac2
          ccCofB(3,ijkC) = 0.5_RFREAL*fac3
          ccCofB(4,ijkC) = 0.5_RFREAL*fac4
          ccCofB(5,ijkC) = 0.5_RFREAL*fac5
        ENDDO
      ENDDO 
    ENDDO     
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloLesGenCoCCHi

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floLesGenCoCCUtil.F90,v $
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








