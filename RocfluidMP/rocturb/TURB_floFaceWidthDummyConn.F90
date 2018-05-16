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
!          to the dummy points at connecting boundaries.
!
! Description: Dummy values of connecting bnd has been computed in prev. step 
!              up to second outerst layer. They are copied to the
!              outerst layer in this routine.
!
! Input: region   = current region data
!        lbound   = patch boundary ID
!        i,j,kdir = direction identifier
!        i,j,kndBeg = i,j,k begin index
!        i,j,kndEnd = i,j,k end index
!  
! Output: Face widths workI, workJ, workK at outerst layer of conn. bc.
!
! Notes: Mother routine = TURB_FloFaceWidth.
!        This routine is only called to complete the face widths at the 
!        right outerst faces (lbound = 2, 4 and 6).
!
!******************************************************************************
!
! $Id: TURB_floFaceWidthDummyConn.F90,v 1.5 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloFaceWidthDummyConn( region,lbound,idir,jdir,kdir, &
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
  INTEGER        :: lbound,idir,jdir,kdir
  INTEGER        :: indBeg, indEnd, jndBeg, jndEnd, kndBeg, kndEnd

! ... loop variables
  INTEGER :: i, j, k, idum, jdum, kdum

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global


  INTEGER :: iLev, iNOff, ijNOff
  INTEGER :: ijkNI, ijkNJ, ijkNK, ijkND
  INTEGER :: iG(3), jG(3), kG(3)
  INTEGER :: iDumB(3), jDumB(3), kDumB(3), iDumE(3), jDumE(3), kDumE(3)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floFaceWidthDummyConn.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloFaceWidthDummyConn',&
  'TURB_floFaceWidthDummyConn.F90' )

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
              ijkND  = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff)
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
              ijkND  = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff)
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

END SUBROUTINE TURB_FloFaceWidthDummyConn

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_floFaceWidthDummyConn.F90,v $
! Revision 1.5  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/08/04 02:49:32  wasistho
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
!
!******************************************************************************







