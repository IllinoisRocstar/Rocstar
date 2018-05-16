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
! Purpose: Extrapolate coefficient from the interior domain or at the patch 
!          to the dummy points at physical and symmetry boundaries.
!
! Description: The averaging coefficients at dummy faces are obtained by
!              mirroring from the interior faces.
!
! Input: region   = current region data
!        lbound   = patch boundary ID
!        i,j,kdir = direction identifier
!        i,j,kndBeg = i,j,k begin index
!        i,j,kndEnd = i,j,k end index
!  
! Output: Averaging coefficients c2fCoI, c2fCoJ, c2fCoK at dummy faces.
!
! Notes: Mother routine = RFLO_C2fAvgCoeffsDummy.
!
!******************************************************************************
!
! $Id: RFLO_C2fAvgCoeffsDummyPhys.F90,v 1.3 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_C2fAvgCoeffsDummyPhys( region,lbound,idir,jdir,kdir, &
                                 indBeg,indEnd,jndBeg,jndEnd,kndBeg,kndEnd )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: lbound, idir, jdir, kdir
  INTEGER        :: indBeg, indEnd, jndBeg, jndEnd, kndBeg, kndEnd

! ... loop variables
  INTEGER :: i, j, k, idum, jdum, kdum

! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: iLev, iNOff, ijNOff
  INTEGER :: ijkNB, ijkNI, ijkNJ, ijkNK, ijkND
  INTEGER :: coIndxI(2), coIndxJ(2), coIndxK(2)
  INTEGER :: iG(3), jG(3), kG(3), iR(3), jR(3), kR(3)
  INTEGER :: iDumB(3), jDumB(3), kDumB(3), iDumE(3), jDumE(3), kDumE(3)

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RFLO_C2fAvgCoeffsDummyPhys',&
  'RFLO_C2fAvgCoeffsDummyPhys.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

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
  IF (lbound == 1) THEN
    iG(1)=0
  ELSE IF (lbound == 3) THEN
    jG(2)=0
  ELSE IF (lbound == 5) THEN
    kG(3)=0
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

  IF (lbound==1 .OR. lbound==2) THEN
    coIndxI(1) = 2
    coIndxI(2) = 1
    coIndxJ(1) = 1
    coIndxJ(2) = 2
    coIndxK(1) = 1
    coIndxK(2) = 2
  ELSEIF (lbound==3 .OR. lbound==4) THEN
    coIndxI(1) = 1
    coIndxI(2) = 2
    coIndxJ(1) = 2
    coIndxJ(2) = 1
    coIndxK(1) = 1
    coIndxK(2) = 2
  ELSEIF (lbound==5 .OR. lbound==6) THEN
    coIndxI(1) = 1
    coIndxI(2) = 2
    coIndxJ(1) = 1
    coIndxJ(2) = 2
    coIndxK(1) = 2
    coIndxK(2) = 1
  ENDIF

! specify averaging coefficients at dummy points by mirroring
! 2D loop over patch nodes 

  DO k=kndBeg,kndEnd
    DO j=jndBeg,jndEnd
      DO i=indBeg,indEnd

! ----- 1D loop in direction normal to patch to define dummy face values

        DO idum=iDumB(1),iDumE(1)
          DO jdum=jDumB(1),jDumE(1)
            DO kdum=kDumB(1),kDumE(1) 
              ijkND = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff)  ! dummy points
              ijkNI = IndIJK(i+iG(1)-idum,j+jG(1)-jdum,k+kG(1)-kdum,iNOff,ijNOff)
              region%levels(iLev)%grid%c2fCoI(coIndxI(2),ijkND) = &
              region%levels(iLev)%grid%c2fCoI(2,ijkNI)
              region%levels(iLev)%grid%c2fCoI(coIndxI(1),ijkND) = &
              1._RFREAL-region%levels(iLev)%grid%c2fCoI(coIndxI(2),ijkND)
            ENDDO
          ENDDO  
        ENDDO    

        DO idum=iDumB(2),iDumE(2)
          DO jdum=jDumB(2),jDumE(2)
            DO kdum=kDumB(2),kDumE(2) 
              ijkND = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff) ! dummy points
              ijkNJ = IndIJK(i+iG(2)-idum,j+jG(2)-jdum,k+kG(2)-kdum,iNOff,ijNOff)                      ! interior
              region%levels(iLev)%grid%c2fCoJ(coIndxJ(2),ijkND) = &
              region%levels(iLev)%grid%c2fCoJ(2,ijkNJ)
              region%levels(iLev)%grid%c2fCoJ(coIndxJ(1),ijkND) = &
              1._RFREAL-region%levels(iLev)%grid%c2fCoJ(coIndxJ(2),ijkND)
            ENDDO
          ENDDO
        ENDDO

        DO idum=iDumB(3),iDumE(3)
          DO jdum=jDumB(3),jDumE(3)
            DO kdum=kDumB(3),kDumE(3) 
              ijkND = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff) ! dummy points
              ijkNK = IndIJK(i+iG(3)-idum,j+jG(3)-jdum,k+kG(3)-kdum,iNOff,ijNOff)                      ! interior
              region%levels(iLev)%grid%c2fCoK(coIndxK(2),ijkND) = &
              region%levels(iLev)%grid%c2fCoK(2,ijkNK)
              region%levels(iLev)%grid%c2fCoK(coIndxK(1),ijkND) = &
              1._RFREAL-region%levels(iLev)%grid%c2fCoK(coIndxK(2),ijkND)
            ENDDO
          ENDDO
        ENDDO

      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_C2fAvgCoeffsDummyPhys

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: RFLO_C2fAvgCoeffsDummyPhys.F90,v $
! Revision 1.3  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.3  2004/08/03 00:52:13  wasistho
! changed avgCo to c2fCo in the description
!
! Revision 1.2  2004/08/02 19:32:48  wasistho
! changed grid%avgCo to grid%c2fCo
!
! Revision 1.1  2004/07/30 17:30:31  wasistho
! initial import routines starting with RFLO_c2fAvg...
!
!
!
!******************************************************************************







