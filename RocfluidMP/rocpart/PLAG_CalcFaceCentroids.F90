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
! Purpose: calculate face centroids.
!
! Description: none.
!
! Input: region = current region.
!
! Output: plag%fc = face centroids.
!
! Notes: face centroids are defined for dummy cells.
!
!******************************************************************************
!
! $Id: PLAG_CalcFaceCentroids.F90,v 1.3 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CalcFaceCentroids( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, &
                            RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#define IndIJK(x,y,z,o1,o2) ((x) + ((y) -1)*o1 + ((z)-1)*o2)

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, iLev, j, k

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iCOff, ijCOff, ijkC, idnbeg, idnend, iNOff, ijNOff, ijkN, &
             jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER, DIMENSION(4) :: corner
  
  REAL(RFREAL)          :: xyzQuad(3,4)
  REAL(RFREAL), POINTER, DIMENSION(:,:)   :: xyz
  REAL(RFREAL), POINTER, DIMENSION(:,:,:) :: fc
  
  TYPE(t_global), POINTER :: global
  TYPE(t_plag),   POINTER :: pPlag
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CalcFaceCentroids.F90,v $ $Revision: 1.3 $'

  global => region%global
  
  CALL RegisterFunction( global,'PLAG_CalcFaceCentroids',&
  'PLAG_CalcFaceCentroids.F90' )
    
! Loop over all grid levels ---------------------------------------------------

  DO iLev=1,region%nGridLevels
  
! - Get dimensions ------------------------------------------------------------

    CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
    CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

! - Set pointers --------------------------------------------------------------

    pPlag => region%levels(iLev)%plag
    xyz   => region%levels(iLev)%grid%xyz
    fc    => pPlag%fc

! - i-direction ---------------------------------------------------------------

    DO k=kdnbeg,kdnend-1
      DO j=jdnbeg,jdnend-1
        DO i=idnbeg,idnend
          corner(1) = IndIJK(i,j  ,k  ,iNOff,ijNOff)
          corner(2) = IndIJK(i,j  ,k+1,iNOff,ijNOff)
          corner(3) = IndIJK(i,j+1,k+1,iNOff,ijNOff)
          corner(4) = IndIJK(i,j+1,k  ,iNOff,ijNOff)
          xyzQuad(1:3,1) = xyz(1:3,corner(1))
          xyzQuad(1:3,2) = xyz(1:3,corner(2))
          xyzQuad(1:3,3) = xyz(1:3,corner(3))
          xyzQuad(1:3,4) = xyz(1:3,corner(4))
          fc(1:3,ICOORD,corner(1)) = 0.25_RFREAL*     &
                                    (xyzQuad(1:3,1) + &
                                     xyzQuad(1:3,2) + &
                                     xyzQuad(1:3,3) + &
                                     xyzQuad(1:3,4) )    
        ENDDO
      ENDDO
    ENDDO

! - j-direction ---------------------------------------------------------------

    DO k=kdnbeg,kdnend-1
      DO j=jdnbeg,jdnend
        DO i=idnbeg,idnend-1
          corner(1) = IndIJK(i  ,j,k  ,iNOff,ijNOff)
          corner(2) = IndIJK(i+1,j,k  ,iNOff,ijNOff)
          corner(3) = IndIJK(i+1,j,k+1,iNOff,ijNOff)
          corner(4) = IndIJK(i  ,j,k+1,iNOff,ijNOff)
          xyzQuad(1:3,1) = xyz(1:3,corner(1))
          xyzQuad(1:3,2) = xyz(1:3,corner(2))
          xyzQuad(1:3,3) = xyz(1:3,corner(3))
          xyzQuad(1:3,4) = xyz(1:3,corner(4))
          fc(1:3,JCOORD,corner(1)) = 0.25_RFREAL*     &
                                    (xyzQuad(1:3,1) + &
                                     xyzQuad(1:3,2) + &
                                     xyzQuad(1:3,3) + &
                                     xyzQuad(1:3,4) ) 
        ENDDO
      ENDDO
    ENDDO

! - k-direction ---------------------------------------------------------------

    DO k=kdnbeg,kdnend
      DO j=jdnbeg,jdnend-1
        DO i=idnbeg,idnend-1
          corner(1) = IndIJK(i  ,j  ,k,iNOff,ijNOff)
          corner(2) = IndIJK(i  ,j+1,k,iNOff,ijNOff)
          corner(3) = IndIJK(i+1,j+1,k,iNOff,ijNOff)
          corner(4) = IndIJK(i+1,j  ,k,iNOff,ijNOff)
          xyzQuad(1:3,1) = xyz(1:3,corner(1))
          xyzQuad(1:3,2) = xyz(1:3,corner(2))
          xyzQuad(1:3,3) = xyz(1:3,corner(3))
          xyzQuad(1:3,4) = xyz(1:3,corner(4))
          fc(1:3,KCOORD,corner(1)) = 0.25_RFREAL*     &
                                    (xyzQuad(1:3,1) + &
                                     xyzQuad(1:3,2) + &
                                     xyzQuad(1:3,3) + &
                                     xyzQuad(1:3,4) ) 
        ENDDO
      ENDDO
    ENDDO

  ENDDO   ! iLev

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CalcFaceCentroids

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CalcFaceCentroids.F90,v $
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:02  fnajjar
! Initial revision after changing case
!
! Revision 1.2  2003/05/02 20:10:55  fnajjar
! Redefined IndIJK based on cpp for performance speedup
!
! Revision 1.1  2002/10/25 14:14:44  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************







