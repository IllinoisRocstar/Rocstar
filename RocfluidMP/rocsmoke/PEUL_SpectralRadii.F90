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
! Purpose: calculate spectral radii for smoke
!
! Description: none.
!
! Input: regions = data of all grid regions.
!
! Output: regions%levels%peul%srad = convective spectral radii for smoke
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_SpectralRadii.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_SpectralRadii( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset, &
        RFLO_GetNodeOffset, RFLO_copyVectorPatches, RFLO_copyVectorEdges, &
        RFLO_copyVectorCorners, RFLO_copyMatrixPatches, RFLO_copyMatrixEdges, &
        RFLO_copyMatrixCorners, RFLO_GetDimensPhys
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff, iNOff, ijNOff, ijkC, ijkN, ijkN1, indSvel

  REAL(RFREAL) :: rrho, u, v, w, sx, sy, sz, dS, sVel, vc
  REAL(RFREAL), POINTER :: cv(:,:)
  REAL(RFREAL), POINTER :: si(:,:), sj(:,:), sk(:,:), vol(:)
  REAL(RFREAL), POINTER :: siVel(:), sjVel(:), skVel(:)
  REAL(RFREAL), POINTER :: srad(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_SpectralRadii.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_SpectralRadii',&
  'PEUL_SpectralRadii.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend,jdcbeg,jdcend, &
                            kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  cv => region%levels(iLev)%mixt%cv

  si     => region%levels(iLev)%grid%si
  sj     => region%levels(iLev)%grid%sj
  sk     => region%levels(iLev)%grid%sk
  vol    => region%levels(iLev)%grid%vol
  siVel  => region%levels(iLev)%grid%siVel
  sjVel  => region%levels(iLev)%grid%sjVel
  skVel  => region%levels(iLev)%grid%skVel
  indSvel = region%levels(iLev)%grid%indSvel

  srad => region%levels(iLev)%peul%srad

! local time step, spectral radii ---------------------------------------------

  DO k=kdcbeg,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend
        ijkC = IndIJK(i,j,k,iCoff,ijCOff)
        ijkN = IndIJK(i,j,k,iNOff,ijNOff)
        rrho = 1._RFREAL/cv(CV_MIXT_DENS,ijkC)
        u    = cv(CV_MIXT_XMOM,ijkC)*rrho
        v    = cv(CV_MIXT_YMOM,ijkC)*rrho
        w    = cv(CV_MIXT_ZMOM,ijkC)*rrho

        ijkN1 = IndIJK(i+1,j,k,iNOff,ijNOff)
        sx    = 0.5_RFREAL*(si(XCOORD,ijkN)+si(XCOORD,ijkN1))
        sy    = 0.5_RFREAL*(si(YCOORD,ijkN)+si(YCOORD,ijkN1))
        sz    = 0.5_RFREAL*(si(ZCOORD,ijkN)+si(ZCOORD,ijkN1))
        dS    = SQRT(sx*sx+sy*sy+sz*sz)
        sVel  = 0.5_RFREAL*(siVel(ijkN*indSvel)+siVel(ijkN1*indSvel))
        vc    = sx*u + sy*v + sz*w - sVel*dS
        srad(ICOORD,ijkC) = ABS(vc)

        ijkN1 = IndIJK(i,j+1,k,iNOff,ijNOff)
        sx    = 0.5_RFREAL*(sj(XCOORD,ijkN)+sj(XCOORD,ijkN1))
        sy    = 0.5_RFREAL*(sj(YCOORD,ijkN)+sj(YCOORD,ijkN1))
        sz    = 0.5_RFREAL*(sj(ZCOORD,ijkN)+sj(ZCOORD,ijkN1))
        dS    = SQRT(sx*sx+sy*sy+sz*sz)
        sVel  = 0.5_RFREAL*(sjVel(ijkN*indSvel)+sjVel(ijkN1*indSvel))
        vc    = sx*u + sy*v + sz*w - sVel*dS
        srad(JCOORD,ijkC) = ABS(vc)

        ijkN1 = IndIJK(i,j,k+1,iNOff,ijNOff)
        sx    = 0.5_RFREAL*(sk(XCOORD,ijkN)+sk(XCOORD,ijkN1))
        sy    = 0.5_RFREAL*(sk(YCOORD,ijkN)+sk(YCOORD,ijkN1))
        sz    = 0.5_RFREAL*(sk(ZCOORD,ijkN)+sk(ZCOORD,ijkN1))
        dS    = SQRT(sx*sx+sy*sy+sz*sz)
        sVel  = 0.5_RFREAL*(skVel(ijkN*indSvel)+skVel(ijkN1*indSvel))
        vc    = sx*u + sy*v + sz*w - sVel*dS
        srad(KCOORD,ijkC) = ABS(vc)

      ENDDO  ! i
    ENDDO    ! j
  ENDDO      ! k

! treat dummy cells -----------------------------------------------------------

  CALL RFLO_copyMatrixPatches( iLev,region,srad )
  CALL RFLO_copyMatrixEdges( iLev,region,srad )
  CALL RFLO_copyMatrixCorners( iLev,region,srad )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_SpectralRadii

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_SpectralRadii.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:10:01  haselbac
! Initial revision after changing case
!
! Revision 1.2  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.1  2003/02/11 22:52:51  jferry
! Initial import of Rocsmoke
!
!
!******************************************************************************







