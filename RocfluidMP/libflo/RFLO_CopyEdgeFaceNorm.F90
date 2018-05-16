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
! Purpose: set values in edge and corner faces of a region
!
! Description: - every edge has 2 sides of dummy faces at which the values to 
!                be copied are known from previous treatments 
!              - face values at edges and corners are then obtaine from 
!                averaging of two known values at each side boundary
!              - face variables set in edges and corners are those whose faces
!                NORMAL to faces of the reference values 
!
! Input: region = region data.
!        iFBeg  = begin index of input face variables
!        iFEnd  = end index of input face variables
!        fvari,fvarj,fvark = input variables at face i, j and k
!
! Output: new values of face variables at edges and corners
!
! Notes: - corners are not treated separately but as elongation of edges;
!          hence only eight such extended edges are treated covering corners
!          at their ends
!        - this routine makes use of averaging routine AverageVecVar
!          in libfloflu
!
!******************************************************************************
!
! $Id: RFLO_CopyEdgeFaceNorm.F90,v 1.3 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CopyEdgeFaceNorm( region,iFBeg,iFEnd,fvari,fvarj,fvark )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : AverageVecVar, RFLO_GetDimensDummyNodes, &
                            RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset
  USE ModError
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: iFBeg, iFEnd
  REAL(RFREAL), POINTER :: fvari(:,:), fvarj(:,:), fvark(:,:)

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: iLev, iNOff, ijNOff, ijkD, ijkN1, ijkN2

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CopyEdgeFaceNorm',&
  'RFLO_CopyEdgeFaceNorm.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! edges 9, 10, 11, 12 ---------------------------------------------------------

  DO i=idnbeg,idnend
    DO j=jpnbeg-1,jdnbeg,-1         ! edges 9, 10
      DO k=kpnbeg-1,kdnbeg,-1
        ijkD  = IndIJK(i,j     ,k     ,iNOff,ijNOff)
        ijkN1 = IndIJK(i,jpnbeg,k     ,iNOff,ijNOff)
        ijkN2 = IndIJK(i,j     ,kpnbeg,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvarj)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvark)
      ENDDO
      DO k=kpnend,kdnend
        ijkD  = IndIJK(i,j     ,k      ,iNOff,ijNOff)
        ijkN1 = IndIJK(i,jpnbeg,k      ,iNOff,ijNOff)
        ijkN2 = IndIJK(i,j     ,kpnend-1,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvarj)
      ENDDO
      DO k=kpnend+1,kdnend
        ijkD  = IndIJK(i,j     ,k     ,iNOff,ijNOff)
        ijkN1 = IndIJK(i,jpnbeg,k     ,iNOff,ijNOff)
        ijkN2 = IndIJK(i,j     ,kpnend,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvark)
      ENDDO
    ENDDO

    DO j=jpnend,jdnend           ! edges 11, 12
      DO k=kpnbeg-1,kdnbeg,-1
        ijkD  = IndIJK(i,j     ,k     ,iNOff,ijNOff)
        ijkN1 = IndIJK(i,jpnend-1,k     ,iNOff,ijNOff)
        ijkN2 = IndIJK(i,j     ,kpnbeg,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvark)
      ENDDO
      DO k=kpnend+1,kdnend
        ijkD  = IndIJK(i,j     ,k     ,iNOff,ijNOff)
        ijkN1 = IndIJK(i,jpnend-1,k     ,iNOff,ijNOff)
        ijkN2 = IndIJK(i,j     ,kpnend,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvark)
      ENDDO
    ENDDO

    DO j=jpnend+1,jdnend           ! edges 11, 12
      DO k=kpnbeg-1,kdnbeg,-1
        ijkD  = IndIJK(i,j     ,k     ,iNOff,ijNOff)
        ijkN1 = IndIJK(i,jpnend,k     ,iNOff,ijNOff)
        ijkN2 = IndIJK(i,j     ,kpnbeg,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvarj)
      ENDDO
      DO k=kpnend,kdnend
        ijkD  = IndIJK(i,j     ,k     ,iNOff,ijNOff)
        ijkN1 = IndIJK(i,jpnend,k     ,iNOff,ijNOff)
        ijkN2 = IndIJK(i,j     ,kpnend-1,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvarj)
      ENDDO
    ENDDO
  ENDDO

! edges 2, 4, 6, 8 ------------------------------------------------------------

  DO j=jdnbeg,jdnend
    DO i=ipnbeg-1,idnbeg,-1        ! edges 2, 4
      DO k=kpnend,kdnend
        ijkD  = IndIJK(i  ,j,k  ,iNOff,ijNOff)
        ijkN1 = IndIJK(ipnbeg,j,k  ,iNOff,ijNOff)
        ijkN2 = IndIJK(i  ,j,kpnend-1,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvari)
      ENDDO
      DO k=kpnend+1,kdnend
        ijkD  = IndIJK(i  ,j,k  ,iNOff,ijNOff)
        ijkN1 = IndIJK(ipnbeg,j,k  ,iNOff,ijNOff)
        ijkN2 = IndIJK(i  ,j,kpnend,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvark)
      ENDDO
      DO k=kpnbeg-1,kdnbeg,-1
        ijkD  = IndIJK(i     ,j,k     ,iNOff,ijNOff)
        ijkN1 = IndIJK(ipnbeg,j,k     ,iNOff,ijNOff)
        ijkN2 = IndIJK(i     ,j,kpnbeg,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvari)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvark)
      ENDDO
    ENDDO

    DO i=ipnend,idnend           ! edges 6, 8
      DO k=kpnend+1,kdnend
        ijkD  = IndIJK(i     ,j,k     ,iNOff,ijNOff)
        ijkN1 = IndIJK(ipnend-1,j,k     ,iNOff,ijNOff)
        ijkN2 = IndIJK(i     ,j,kpnend,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvark)
      ENDDO
      DO k=kpnbeg-1,kdnbeg,-1
        ijkD  = IndIJK(i     ,j,k     ,iNOff,ijNOff)
        ijkN1 = IndIJK(ipnend-1,j,k     ,iNOff,ijNOff)
        ijkN2 = IndIJK(i     ,j,kpnbeg,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvark)
      ENDDO
    ENDDO

    DO i=ipnend+1,idnend           ! edges 6, 8
      DO k=kpnend,kdnend
        ijkD  = IndIJK(i     ,j,k     ,iNOff,ijNOff)
        ijkN1 = IndIJK(ipnend,j,k     ,iNOff,ijNOff)
        ijkN2 = IndIJK(i     ,j,kpnend-1,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvari)
      ENDDO
      DO k=kpnbeg-1,kdnbeg,-1
        ijkD  = IndIJK(i     ,j,k     ,iNOff,ijNOff)
        ijkN1 = IndIJK(ipnend,j,k     ,iNOff,ijNOff)
        ijkN2 = IndIJK(i     ,j,kpnbeg,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvari)
      ENDDO
    ENDDO
  ENDDO

! edges 1, 3, 5, 7 ------------------------------------------------------------

  DO k=kdnbeg,kdnend
    DO i=ipnbeg-1,idnbeg,-1        ! edges 1, 3
      DO j=jpnbeg-1,jdnbeg,-1
        ijkD  = IndIJK(i     ,j     ,k,iNOff,ijNOff)
        ijkN1 = IndIJK(ipnbeg,j     ,k,iNOff,ijNOff)
        ijkN2 = IndIJK(i     ,jpnbeg,k,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvari)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvarj)
      ENDDO
      DO j=jpnend,jdnend
        ijkD  = IndIJK(i     ,j     ,k,iNOff,ijNOff)
        ijkN1 = IndIJK(ipnbeg,j     ,k,iNOff,ijNOff)
        ijkN2 = IndIJK(i     ,jpnend-1,k,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvari)
      ENDDO
      DO j=jpnend+1,jdnend
        ijkD  = IndIJK(i     ,j     ,k,iNOff,ijNOff)
        ijkN1 = IndIJK(ipnbeg,j     ,k,iNOff,ijNOff)
        ijkN2 = IndIJK(i     ,jpnend,k,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvarj)
      ENDDO
    ENDDO

    DO i=ipnend,idnend           ! edges 5, 7
      DO j=jpnbeg-1,jdnbeg,-1
        ijkD  = IndIJK(i     ,j     ,k,iNOff,ijNOff)
        ijkN1 = IndIJK(ipnend-1,j     ,k,iNOff,ijNOff)
        ijkN2 = IndIJK(i     ,jpnbeg,k,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvarj)
      ENDDO
      DO j=jpnend+1,jdnend
        ijkD  = IndIJK(i     ,j       ,k,iNOff,ijNOff)
        ijkN1 = IndIJK(ipnend-1,j       ,k,iNOff,ijNOff)
        ijkN2 = IndIJK(i     ,jpnend,k,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvarj)
      ENDDO
    ENDDO
    DO i=ipnend+1,idnend           ! edges 5, 7
      DO j=jpnbeg-1,jdnbeg,-1
        ijkD  = IndIJK(i     ,j     ,k,iNOff,ijNOff)
        ijkN1 = IndIJK(ipnend,j     ,k,iNOff,ijNOff)
        ijkN2 = IndIJK(i     ,jpnbeg,k,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvari)
      ENDDO
      DO j=jpnend,jdnend
        ijkD  = IndIJK(i     ,j       ,k,iNOff,ijNOff)
        ijkN1 = IndIJK(ipnend,j       ,k,iNOff,ijNOff)
        ijkN2 = IndIJK(i     ,jpnend-1,k,iNOff,ijNOff)
        CALL AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvari)
      ENDDO
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_CopyEdgeFaceNorm

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CopyEdgeFaceNorm.F90,v $
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
! Revision 1.4  2002/07/30 02:43:09  wasistho
! AverageVecVar moved to libfloflu
!
! Revision 1.3  2002/07/22 22:59:11  jblazek
! Some more clean up.
!
! Revision 1.2  2002/07/19 23:45:00  wasistho
! made compliant with CODING RULE
!
! Revision 1.1  2002/05/21 02:05:34  wasistho
! add viscous terms
!
!******************************************************************************







