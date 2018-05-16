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
! Purpose: Extrapolate face widths from the interior domain or at the patch 
!          to the dummy points at physical and symmetry boundaries.
!
! Description: The face widths at dummy faces are obtained by
!              mirroring from the interior faces.
!
! Input: region   = current region data
!        lbound   = patch boundary ID
!        i,j,kdir = direction identifier
!        i,j,kndBeg = i,j,k begin index
!        i,j,kndEnd = i,j,k end index
!  
! Output: Face widths workI, workJ, workK at dummy faces.
!
! Notes: Mother routine = TURB_FloFaceWidth.
!
!******************************************************************************
!
! $Id: TURB_floFaceWidthDummyPhys.F90,v 1.5 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloFaceWidthDummyPhys( region,lbound,idir,jdir,kdir, &
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
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: iLev, iNOff, ijNOff
  INTEGER :: ijkNB, ijkNI, ijkNJ, ijkNK, ijkND
  INTEGER :: iG(3), jG(3), kG(3), iR(3), jR(3), kR(3)
  INTEGER :: iDumB(3), jDumB(3), kDumB(3), iDumE(3), jDumE(3), kDumE(3)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floFaceWidthDummyPhys.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloFaceWidthDummyPhys',&
  'TURB_floFaceWidthDummyPhys.F90' )

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

! specify facewidths at dummy faces by mirroring
! 2D loop over patch nodes 

  DO k=kndBeg,kndEnd
    DO j=jndBeg,jndEnd
      DO i=indBeg,indEnd

! ----- 1D loop in direction normal to patch to define dummy face values

        DO idum=iDumB(1),iDumE(1)
          DO jdum=jDumB(1),jDumE(1)
            DO kdum=kDumB(1),kDumE(1) 
              ijkND = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff)  ! dummy points
              ijkNI = IndIJK(i+iG(1)-idum,j+jG(1)-jdum,k+kG(1)-kdum,iNOff,ijNOff) ! interior
              region%levels(iLev)%turb%workI(1,ijkND) = &
              region%levels(iLev)%turb%workI(1,ijkNI)
              region%levels(iLev)%turb%workI(2,ijkND) = &
              region%levels(iLev)%turb%workI(2,ijkNI)
            ENDDO
          ENDDO  
        ENDDO    

        DO idum=iDumB(2),iDumE(2)
          DO jdum=jDumB(2),jDumE(2)
            DO kdum=kDumB(2),kDumE(2) 
              ijkND = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff) ! dummy points
              ijkNJ = IndIJK(i+iG(2)-idum,j+jG(2)-jdum,k+kG(2)-kdum,iNOff,ijNOff) ! interior
              region%levels(iLev)%turb%workJ(1,ijkND) = &
              region%levels(iLev)%turb%workJ(1,ijkNJ)
              region%levels(iLev)%turb%workJ(2,ijkND) = &
              region%levels(iLev)%turb%workJ(2,ijkNJ)
            ENDDO
          ENDDO
        ENDDO

        DO idum=iDumB(3),iDumE(3)
          DO jdum=jDumB(3),jDumE(3)
            DO kdum=kDumB(3),kDumE(3) 
              ijkND = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff) ! dummy points
              ijkNK = IndIJK(i+iG(3)-idum,j+jG(3)-jdum,k+kG(3)-kdum,iNOff,ijNOff) ! interior
              region%levels(iLev)%turb%workK(1,ijkND) = &
              region%levels(iLev)%turb%workK(1,ijkNK)
              region%levels(iLev)%turb%workK(2,ijkND) = &
              region%levels(iLev)%turb%workK(2,ijkNK)
            ENDDO
          ENDDO
        ENDDO

      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloFaceWidthDummyPhys

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_floFaceWidthDummyPhys.F90,v $
! Revision 1.5  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/08/04 02:49:46  wasistho
! removed turb%avgCoI,J,K as it is defined as grid%c2fCoI,J,K
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
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







