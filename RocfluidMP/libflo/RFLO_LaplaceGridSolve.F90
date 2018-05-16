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
! Purpose: conduct one Jacobi iteration to obtain new grid movements
!          in the flow domain (boundaries included).
!
! Description: none.
!
! Input: region = data of current region, old grid coordinates.
!
! Output: region%levels%grid%xyz = grid movements.
!
! Notes: on entry, xyz holds node coordinates from a previous smoothing
!        step. On exit however, xyz contains only the grid motion. 
!
!******************************************************************************
!
! $Id: RFLO_LaplaceGridSolve.F90,v 1.7 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_LaplaceGridSolve( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetDimensPhysNodes, &
                            RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend, iLev, iNOff, ijNOff
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, ibn, ien
  INTEGER :: ijk, ip1, ip2, im1, im2, jp1, jp2, jm1, jm2, kp1, kp2, km1, km2
  INTEGER :: moveBlock, method, in, ico, ic(8)

  LOGICAL, POINTER :: bndMoved(:)

  REAL(RFREAL) :: rxi, ryi, rzi, rxj, ryj, rzj, rxk, ryk, rzk, wi, wj, wk, &
                  si, sj, sk, sim, sjm, skm, sip, sjp, skp, d, p
  REAL(RFREAL) :: denom(8), dist(8,8)
  REAL(RFREAL), POINTER :: xyz(:,:), xyzOld(:,:), xyzOrig(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_LaplaceGridSolve',&
  'RFLO_LaplaceGridSolve.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = 1

  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend,jpnbeg,jpnend, &
                                kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

  xyz      => region%levels(iLev)%grid%xyz
  xyzOld   => region%levels(iLev)%grid%xyzOld
  xyzOrig  => region%levels(iLev)%gridOld%xyz
  bndMoved => region%levels(iLev)%grid%boundMoved
  p        =  region%global%moveGridWeight

! reset motion vectors --------------------------------------------------------

  DO ijk=ibn,ien
    xyzOld(XCOORD,ijk) = xyz(XCOORD,ijk) - xyzOrig(XCOORD,ijk)
    xyzOld(YCOORD,ijk) = xyz(YCOORD,ijk) - xyzOrig(YCOORD,ijk)
    xyzOld(ZCOORD,ijk) = xyz(ZCOORD,ijk) - xyzOrig(ZCOORD,ijk)
  ENDDO

! move block corners locally

  moveBlock = 0
  method    = 1

  IF (moveBlock==1) THEN
    ic(1) = IndIJK(ipnbeg ,jpnbeg ,kpnbeg ,iNOff,ijNOff)
    ic(2) = IndIJK(ipnbeg ,jpnbeg ,kpnend ,iNOff,ijNOff)
    ic(3) = IndIJK(ipnbeg ,jpnend ,kpnend ,iNOff,ijNOff)
    ic(4) = IndIJK(ipnbeg ,jpnend ,kpnbeg ,iNOff,ijNOff)
    ic(5) = IndIJK(ipnend ,jpnbeg ,kpnbeg ,iNOff,ijNOff)
    ic(6) = IndIJK(ipnend ,jpnbeg ,kpnend ,iNOff,ijNOff)
    ic(7) = IndIJK(ipnend ,jpnend ,kpnend ,iNOff,ijNOff)
    ic(8) = IndIJK(ipnend ,jpnend ,kpnbeg ,iNOff,ijNOff)

    DO ico = 1,8
      denom(ico) = 0._RFREAL
      DO in = 1,8
        IF (ico/=in) THEN
          dist(in,ico) = SQRT( (xyz(XCOORD,ic(in))-xyz(XCOORD,ic(ico)))**2 + &
                               (xyz(YCOORD,ic(in))-xyz(YCOORD,ic(ico)))**2 + &
                               (xyz(ZCOORD,ic(in))-xyz(ZCOORD,ic(ico)))**2 )
          dist(in,ico) = 1._RFREAL/dist(in,ico)
          denom(ico)   = denom(ico) + dist(in,ico)
        ENDIF
      ENDDO
    ENDDO

    IF (.NOT. bndMoved(1)) THEN
      IF ((.NOT. bndMoved(3)).AND.(.NOT. bndMoved(5))) THEN
        xyz(XCOORD,ic(1)) = 0._RFREAL
        xyz(YCOORD,ic(1)) = 0._RFREAL
        xyz(ZCOORD,ic(1)) = 0._RFREAL
        DO in=1,8
          IF (in/=1) THEN
            xyz(XCOORD,ic(1))= xyz(XCOORD,ic(1))+dist(in,1)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(1))= xyz(YCOORD,ic(1))+dist(in,1)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(1))= xyz(ZCOORD,ic(1))+dist(in,1)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(1)) = xyz(XCOORD,ic(1))/denom(1)
        xyzOld(YCOORD,ic(1)) = xyz(YCOORD,ic(1))/denom(1)
        xyzOld(ZCOORD,ic(1)) = xyz(ZCOORD,ic(1))/denom(1)
      ENDIF
      IF ((.NOT. bndMoved(3)).AND.(.NOT. bndMoved(6))) THEN
        xyz(XCOORD,ic(2)) = 0._RFREAL
        xyz(YCOORD,ic(2)) = 0._RFREAL
        xyz(ZCOORD,ic(2)) = 0._RFREAL
        DO in=1,8
          IF (in/=2) THEN
            xyz(XCOORD,ic(2))= xyz(XCOORD,ic(2))+dist(in,2)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(2))= xyz(YCOORD,ic(2))+dist(in,2)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(2))= xyz(ZCOORD,ic(2))+dist(in,2)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(2)) = xyz(XCOORD,ic(2))/denom(2)
        xyzOld(YCOORD,ic(2)) = xyz(YCOORD,ic(2))/denom(2)
        xyzOld(ZCOORD,ic(2)) = xyz(ZCOORD,ic(2))/denom(2)
      ENDIF
      IF ((.NOT. bndMoved(4)).AND.(.NOT. bndMoved(6))) THEN
        xyz(XCOORD,ic(3)) = 0._RFREAL
        xyz(YCOORD,ic(3)) = 0._RFREAL
        xyz(ZCOORD,ic(3)) = 0._RFREAL
        DO in=1,8
          IF (in/=3) THEN
            xyz(XCOORD,ic(3))= xyz(XCOORD,ic(3))+dist(in,3)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(3))= xyz(YCOORD,ic(3))+dist(in,3)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(3))= xyz(ZCOORD,ic(3))+dist(in,3)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(3)) = xyz(XCOORD,ic(3))/denom(3)
        xyzOld(YCOORD,ic(3)) = xyz(YCOORD,ic(3))/denom(3)
        xyzOld(ZCOORD,ic(3)) = xyz(ZCOORD,ic(3))/denom(3)
      ENDIF
      IF ((.NOT. bndMoved(4)).AND.(.NOT. bndMoved(5))) THEN
        xyz(XCOORD,ic(4)) = 0._RFREAL
        xyz(YCOORD,ic(4)) = 0._RFREAL
        xyz(ZCOORD,ic(4)) = 0._RFREAL
        DO in=1,8
          IF (in/=4) THEN
            xyz(XCOORD,ic(4))= xyz(XCOORD,ic(4))+dist(in,4)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(4))= xyz(YCOORD,ic(4))+dist(in,4)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(4))= xyz(ZCOORD,ic(4))+dist(in,4)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(4)) = xyz(XCOORD,ic(4))/denom(4)
        xyzOld(YCOORD,ic(4)) = xyz(YCOORD,ic(4))/denom(4)
        xyzOld(ZCOORD,ic(4)) = xyz(ZCOORD,ic(4))/denom(4)
      ENDIF
    ENDIF

    IF (.NOT. bndMoved(2)) THEN
      IF ((.NOT. bndMoved(3)).AND.(.NOT. bndMoved(5))) THEN
        xyz(XCOORD,ic(5)) = 0._RFREAL
        xyz(YCOORD,ic(5)) = 0._RFREAL
        xyz(ZCOORD,ic(5)) = 0._RFREAL
        DO in=1,8
          IF (in/=5) THEN
            xyz(XCOORD,ic(5))= xyz(XCOORD,ic(5))+dist(in,5)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(5))= xyz(YCOORD,ic(5))+dist(in,5)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(5))= xyz(ZCOORD,ic(5))+dist(in,5)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(5)) = xyz(XCOORD,ic(5))/denom(5)
        xyzOld(YCOORD,ic(5)) = xyz(YCOORD,ic(5))/denom(5)
        xyzOld(ZCOORD,ic(5)) = xyz(ZCOORD,ic(5))/denom(5)
      ENDIF
      IF ((.NOT. bndMoved(3)).AND.(.NOT. bndMoved(6))) THEN
        xyz(XCOORD,ic(6)) = 0._RFREAL
        xyz(YCOORD,ic(6)) = 0._RFREAL
        xyz(ZCOORD,ic(6)) = 0._RFREAL
        DO in=1,8
          IF (in/=6) THEN
            xyz(XCOORD,ic(6))= xyz(XCOORD,ic(6))+dist(in,6)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(6))= xyz(YCOORD,ic(6))+dist(in,6)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(6))= xyz(ZCOORD,ic(6))+dist(in,6)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(6)) = xyz(XCOORD,ic(6))/denom(6)
        xyzOld(YCOORD,ic(6)) = xyz(YCOORD,ic(6))/denom(6)
        xyzOld(ZCOORD,ic(6)) = xyz(ZCOORD,ic(6))/denom(6)
      ENDIF
      IF ((.NOT. bndMoved(4)).AND.(.NOT. bndMoved(6))) THEN
        xyz(XCOORD,ic(7)) = 0._RFREAL
        xyz(YCOORD,ic(7)) = 0._RFREAL
        xyz(ZCOORD,ic(7)) = 0._RFREAL
        DO in=1,8
          IF (in/=7) THEN
            xyz(XCOORD,ic(7))= xyz(XCOORD,ic(7))+dist(in,7)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(7))= xyz(YCOORD,ic(7))+dist(in,7)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(7))= xyz(ZCOORD,ic(7))+dist(in,7)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(7)) = xyz(XCOORD,ic(7))/denom(7)
        xyzOld(YCOORD,ic(7)) = xyz(YCOORD,ic(7))/denom(7)
        xyzOld(ZCOORD,ic(7)) = xyz(ZCOORD,ic(7))/denom(7)
      ENDIF
      IF ((.NOT. bndMoved(4)).AND.(.NOT. bndMoved(5))) THEN
        xyz(XCOORD,ic(8)) = 0._RFREAL
        xyz(YCOORD,ic(8)) = 0._RFREAL
        xyz(ZCOORD,ic(8)) = 0._RFREAL
        DO in=1,8
          IF (in/=8) THEN
            xyz(XCOORD,ic(8))= xyz(XCOORD,ic(8))+dist(in,8)*xyzOld(XCOORD,ic(in))
            xyz(YCOORD,ic(8))= xyz(YCOORD,ic(8))+dist(in,8)*xyzOld(YCOORD,ic(in))
            xyz(ZCOORD,ic(8))= xyz(ZCOORD,ic(8))+dist(in,8)*xyzOld(ZCOORD,ic(in))
          ENDIF
        ENDDO
        xyzOld(XCOORD,ic(8)) = xyz(XCOORD,ic(8))/denom(8)
        xyzOld(YCOORD,ic(8)) = xyz(YCOORD,ic(8))/denom(8)
        xyzOld(ZCOORD,ic(8)) = xyz(ZCOORD,ic(8))/denom(8)
      ENDIF
    ENDIF   ! bndMoved
  ENDIF     ! moveBlock corners

! compute new coordinates -----------------------------------------------------

  IF (method==1) THEN

! - New weighted Laplacian Smoothing of Bono Wasistho:

    DO k=kpnbeg,kpnend
     DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        ijk  = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ip1  = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
        ip2  = IndIJK(i+2,j  ,k  ,iNOff,ijNOff)
        im1  = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        im2  = IndIJK(i-2,j  ,k  ,iNOff,ijNOff)
        jp1  = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
        jp2  = IndIJK(i  ,j+2,k  ,iNOff,ijNOff)
        jm1  = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        jm2  = IndIJK(i  ,j-2,k  ,iNOff,ijNOff)
        kp1  = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
        kp2  = IndIJK(i  ,j  ,k+2,iNOff,ijNOff)
        km1  = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)
        km2  = IndIJK(i  ,j  ,k-2,iNOff,ijNOff)

        sim = SQRT( (xyzOrig(XCOORD,im1)+xyzOld(XCOORD,im1)-xyzOrig(XCOORD,ijk)- &
                     xyzOld(XCOORD,ijk))**2 + &
                    (xyzOrig(YCOORD,im1)+xyzOld(YCOORD,im1)-xyzOrig(YCOORD,ijk)- &
                     xyzOld(YCOORD,ijk))**2 + &
                    (xyzOrig(ZCOORD,im1)+xyzOld(ZCOORD,im1)-xyzOrig(ZCOORD,ijk)- &
                     xyzOld(ZCOORD,ijk))**2 )**p       
        sip = SQRT( (xyzOrig(XCOORD,ip1)+xyzOld(XCOORD,ip1)-xyzOrig(XCOORD,ijk)- &
                     xyzOld(XCOORD,ijk))**2 + &
                    (xyzOrig(YCOORD,ip1)+xyzOld(YCOORD,ip1)-xyzOrig(YCOORD,ijk)- &
                     xyzOld(YCOORD,ijk))**2 + &
                    (xyzOrig(ZCOORD,ip1)+xyzOld(ZCOORD,ip1)-xyzOrig(ZCOORD,ijk)- &
                     xyzOld(ZCOORD,ijk))**2 )**p

        sjm = SQRT( (xyzOrig(XCOORD,jm1)+xyzOld(XCOORD,jm1)-xyzOrig(XCOORD,ijk)- &
                     xyzOld(XCOORD,ijk))**2 + &
                    (xyzOrig(YCOORD,jm1)+xyzOld(YCOORD,jm1)-xyzOrig(YCOORD,ijk)- &
                     xyzOld(YCOORD,ijk))**2 + &
                    (xyzOrig(ZCOORD,jm1)+xyzOld(ZCOORD,jm1)-xyzOrig(ZCOORD,ijk)- &
                     xyzOld(ZCOORD,ijk))**2 )**p       
        sjp = SQRT( (xyzOrig(XCOORD,jp1)+xyzOld(XCOORD,jp1)-xyzOrig(XCOORD,ijk)- &
                     xyzOld(XCOORD,ijk))**2 + &
                    (xyzOrig(YCOORD,jp1)+xyzOld(YCOORD,jp1)-xyzOrig(YCOORD,ijk)- &
                     xyzOld(YCOORD,ijk))**2 + &
                    (xyzOrig(ZCOORD,jp1)+xyzOld(ZCOORD,jp1)-xyzOrig(ZCOORD,ijk)- &
                     xyzOld(ZCOORD,ijk))**2 )**p

        skm = SQRT( (xyzOrig(XCOORD,km1)+xyzOld(XCOORD,km1)-xyzOrig(XCOORD,ijk)- &
                     xyzOld(XCOORD,ijk))**2 + &
                    (xyzOrig(YCOORD,km1)+xyzOld(YCOORD,km1)-xyzOrig(YCOORD,ijk)- &
                     xyzOld(YCOORD,ijk))**2 + &
                    (xyzOrig(ZCOORD,km1)+xyzOld(ZCOORD,km1)-xyzOrig(ZCOORD,ijk)- &
                     xyzOld(ZCOORD,ijk))**2 )**p       
        skp = SQRT( (xyzOrig(XCOORD,kp1)+xyzOld(XCOORD,kp1)-xyzOrig(XCOORD,ijk)- &
                     xyzOld(XCOORD,ijk))**2 + &
                    (xyzOrig(YCOORD,kp1)+xyzOld(YCOORD,kp1)-xyzOrig(YCOORD,ijk)- &
                     xyzOld(YCOORD,ijk))**2 + &
                    (xyzOrig(ZCOORD,kp1)+xyzOld(ZCOORD,kp1)-xyzOrig(ZCOORD,ijk)- &
                     xyzOld(ZCOORD,ijk))**2 )**p

        sim = 1._RFREAL/sim       
        sip = 1._RFREAL/sip
        sjm = 1._RFREAL/sjm       
        sjp = 1._RFREAL/sjp
        skm = 1._RFREAL/skm       
        skp = 1._RFREAL/skp

        rxi = sim*xyzOld(XCOORD,im1) + sip*xyzOld(XCOORD,ip1)
        ryi = sim*xyzOld(YCOORD,im1) + sip*xyzOld(YCOORD,ip1)
        rzi = sim*xyzOld(ZCOORD,im1) + sip*xyzOld(ZCOORD,ip1)

        rxj = sjm*xyzOld(XCOORD,jm1) + sjp*xyzOld(XCOORD,jp1)
        ryj = sjm*xyzOld(YCOORD,jm1) + sjp*xyzOld(YCOORD,jp1)
        rzj = sjm*xyzOld(ZCOORD,jm1) + sjp*xyzOld(ZCOORD,jp1)

        rxk = skm*xyzOld(XCOORD,km1) + skp*xyzOld(XCOORD,kp1)
        ryk = skm*xyzOld(YCOORD,km1) + skp*xyzOld(YCOORD,kp1)
        rzk = skm*xyzOld(ZCOORD,km1) + skp*xyzOld(ZCOORD,kp1)

        d = 1._RFREAL/(sim+sip+sjm+sjp+skm+skp)

        xyz(XCOORD,ijk) = (rxi+rxj+rxk)*d
        xyz(YCOORD,ijk) = (ryi+ryj+ryk)*d
        xyz(ZCOORD,ijk) = (rzi+rzj+rzk)*d

      ENDDO   ! i
     ENDDO    ! j
    ENDDO     ! k

  ELSEIF (method==2) THEN

! - Old method by Jiri Blazek

    DO k=kpnbeg,kpnend
     DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        ijk  = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        ip1  = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
        ip2  = IndIJK(i+2,j  ,k  ,iNOff,ijNOff)
        im1  = IndIJK(i-1,j  ,k  ,iNOff,ijNOff)
        im2  = IndIJK(i-2,j  ,k  ,iNOff,ijNOff)
        jp1  = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
        jp2  = IndIJK(i  ,j+2,k  ,iNOff,ijNOff)
        jm1  = IndIJK(i  ,j-1,k  ,iNOff,ijNOff)
        jm2  = IndIJK(i  ,j-2,k  ,iNOff,ijNOff)
        kp1  = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
        kp2  = IndIJK(i  ,j  ,k+2,iNOff,ijNOff)
        km1  = IndIJK(i  ,j  ,k-1,iNOff,ijNOff)
        km2  = IndIJK(i  ,j  ,k-2,iNOff,ijNOff)

        rxi = xyzOld(XCOORD,im1) + xyzOld(XCOORD,ip1)
        ryi = xyzOld(YCOORD,im1) + xyzOld(YCOORD,ip1)
        rzi = xyzOld(ZCOORD,im1) + xyzOld(ZCOORD,ip1)

        rxj = xyzOld(XCOORD,jm1) + xyzOld(XCOORD,jp1)
        ryj = xyzOld(YCOORD,jm1) + xyzOld(YCOORD,jp1)
        rzj = xyzOld(ZCOORD,jm1) + xyzOld(ZCOORD,jp1)

        rxk = xyzOld(XCOORD,km1) + xyzOld(XCOORD,kp1)
        ryk = xyzOld(YCOORD,km1) + xyzOld(YCOORD,kp1)
        rzk = xyzOld(ZCOORD,km1) + xyzOld(ZCOORD,kp1)

        si = (xyzOld(XCOORD,ip1)+xyzOrig(XCOORD,ip1)- &
              xyzOld(XCOORD,ijk)-xyzOrig(XCOORD,ijk))**2 + &
             (xyzOld(YCOORD,ip1)+xyzOrig(YCOORD,ip1)- &
              xyzOld(YCOORD,ijk)-xyzOrig(YCOORD,ijk))**2 + &
             (xyzOld(ZCOORD,ip1)+xyzOrig(ZCOORD,ip1)- &
              xyzOld(ZCOORD,ijk)-xyzOrig(ZCOORD,ijk))**2 + &
             (xyzOld(XCOORD,im1)+xyzOrig(XCOORD,im1)- &
              xyzOld(XCOORD,ijk)-xyzOrig(XCOORD,ijk))**2 + &
             (xyzOld(YCOORD,im1)+xyzOrig(YCOORD,im1)- &
              xyzOld(YCOORD,ijk)-xyzOrig(YCOORD,ijk))**2 + &
             (xyzOld(ZCOORD,im1)+xyzOrig(ZCOORD,im1)- &
              xyzOld(ZCOORD,ijk)-xyzOrig(ZCOORD,ijk))**2 + &
             (xyzOld(XCOORD,ip2)+xyzOrig(XCOORD,ip2)- &
              xyzOld(XCOORD,ip1)-xyzOrig(XCOORD,ip1))**2 + &
             (xyzOld(YCOORD,ip2)+xyzOrig(YCOORD,ip2)- &
              xyzOld(YCOORD,ip1)-xyzOrig(YCOORD,ip1))**2 + &
             (xyzOld(ZCOORD,ip2)+xyzOrig(ZCOORD,ip2)- &
              xyzOld(ZCOORD,ip1)-xyzOrig(ZCOORD,ip1))**2 + &
             (xyzOld(XCOORD,im2)+xyzOrig(XCOORD,im2)- &
              xyzOld(XCOORD,im1)-xyzOrig(XCOORD,im1))**2 + &
             (xyzOld(YCOORD,im2)+xyzOrig(YCOORD,im2)- &
              xyzOld(YCOORD,im1)-xyzOrig(YCOORD,im1))**2 + &
             (xyzOld(ZCOORD,im2)+xyzOrig(ZCOORD,im2)- &
              xyzOld(ZCOORD,im1)-xyzOrig(ZCOORD,im1))**2

        sj = (xyzOld(XCOORD,jp1)+xyzOrig(XCOORD,jp1)- &
              xyzOld(XCOORD,ijk)-xyzOrig(XCOORD,ijk))**2 + &
             (xyzOld(YCOORD,jp1)+xyzOrig(YCOORD,jp1)- &
              xyzOld(YCOORD,ijk)-xyzOrig(YCOORD,ijk))**2 + &
             (xyzOld(ZCOORD,jp1)+xyzOrig(ZCOORD,jp1)- &
              xyzOld(ZCOORD,ijk)-xyzOrig(ZCOORD,ijk))**2 + &
             (xyzOld(XCOORD,jm1)+xyzOrig(XCOORD,jm1)- &
              xyzOld(XCOORD,ijk)-xyzOrig(XCOORD,ijk))**2 + &
             (xyzOld(YCOORD,jm1)+xyzOrig(YCOORD,jm1)- &
              xyzOld(YCOORD,ijk)-xyzOrig(YCOORD,ijk))**2 + &
             (xyzOld(ZCOORD,jm1)+xyzOrig(ZCOORD,jm1)- &
              xyzOld(ZCOORD,ijk)-xyzOrig(ZCOORD,ijk))**2 + &
             (xyzOld(XCOORD,jp2)+xyzOrig(XCOORD,jp2)- &
              xyzOld(XCOORD,jp1)-xyzOrig(XCOORD,jp1))**2 + &
             (xyzOld(YCOORD,jp2)+xyzOrig(YCOORD,jp2)- &
              xyzOld(YCOORD,jp1)-xyzOrig(YCOORD,jp1))**2 + &
             (xyzOld(ZCOORD,jp2)+xyzOrig(ZCOORD,jp2)- &
              xyzOld(ZCOORD,jp1)-xyzOrig(ZCOORD,jp1))**2 + &
             (xyzOld(XCOORD,jm2)+xyzOrig(XCOORD,jm2)- &
              xyzOld(XCOORD,jm1)-xyzOrig(XCOORD,jm1))**2 + &
             (xyzOld(YCOORD,jm2)+xyzOrig(YCOORD,jm2)- &
              xyzOld(YCOORD,jm1)-xyzOrig(YCOORD,jm1))**2 + &
             (xyzOld(ZCOORD,jm2)+xyzOrig(ZCOORD,jm2)- &
              xyzOld(ZCOORD,jm1)-xyzOrig(ZCOORD,jm1))**2

        sk = (xyzOld(XCOORD,kp1)+xyzOrig(XCOORD,kp1)- &
              xyzOld(XCOORD,ijk)-xyzOrig(XCOORD,ijk))**2 + &
             (xyzOld(YCOORD,kp1)+xyzOrig(YCOORD,kp1)- &
              xyzOld(YCOORD,ijk)-xyzOrig(YCOORD,ijk))**2 + &
             (xyzOld(ZCOORD,kp1)+xyzOrig(ZCOORD,kp1)- &
              xyzOld(ZCOORD,ijk)-xyzOrig(ZCOORD,ijk))**2 + &
             (xyzOld(XCOORD,km1)+xyzOrig(XCOORD,km1)- &
              xyzOld(XCOORD,ijk)-xyzOrig(XCOORD,ijk))**2 + &
             (xyzOld(YCOORD,km1)+xyzOrig(YCOORD,km1)- &
              xyzOld(YCOORD,ijk)-xyzOrig(YCOORD,ijk))**2 + &
             (xyzOld(ZCOORD,km1)+xyzOrig(ZCOORD,km1)- &
              xyzOld(ZCOORD,ijk)-xyzOrig(ZCOORD,ijk))**2 + &
             (xyzOld(XCOORD,kp2)+xyzOrig(XCOORD,kp2)- &
              xyzOld(XCOORD,kp1)-xyzOrig(XCOORD,kp1))**2 + &
             (xyzOld(YCOORD,kp2)+xyzOrig(YCOORD,kp2)- &
              xyzOld(YCOORD,kp1)-xyzOrig(YCOORD,kp1))**2 + &
             (xyzOld(ZCOORD,kp2)+xyzOrig(ZCOORD,kp2)- &
              xyzOld(ZCOORD,kp1)-xyzOrig(ZCOORD,kp1))**2 + &
             (xyzOld(XCOORD,km2)+xyzOrig(XCOORD,km2)- &
              xyzOld(XCOORD,km1)-xyzOrig(XCOORD,km1))**2 + &
             (xyzOld(YCOORD,km2)+xyzOrig(YCOORD,km2)- &
              xyzOld(YCOORD,km1)-xyzOrig(YCOORD,km1))**2 + &
             (xyzOld(ZCOORD,km2)+xyzOrig(ZCOORD,km2)- &
              xyzOld(ZCOORD,km1)-xyzOrig(ZCOORD,km1))**2

        wi  = 1._RFREAL/(si)**p
        wj  = 1._RFREAL/(sj)**p
        wk  = 1._RFREAL/(sk)**p
        d   = 2._RFREAL*(wi+wj+wk)

        xyz(XCOORD,ijk) = (wi*rxi+wj*rxj+wk*rxk)/d
        xyz(YCOORD,ijk) = (wi*ryi+wj*ryj+wk*ryk)/d
        xyz(ZCOORD,ijk) = (wi*rzi+wj*rzj+wk*rzk)/d
      ENDDO   ! i
     ENDDO    ! j
    ENDDO     ! k
  ENDIF       ! method

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_LaplaceGridSolve

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_LaplaceGridSolve.F90,v $
! Revision 1.7  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2005/06/03 01:56:48  wasistho
! gather specifying moveBlock and method on top
!
! Revision 1.4  2005/06/02 04:55:25  wasistho
! provide moveBlock and method as internal options
!
! Revision 1.3  2005/05/28 08:54:02  wasistho
! use im1,ip1,etc i.o. im2,ip2 for weights
!
! Revision 1.2  2005/05/28 05:40:45  wasistho
! move block corners and new weighted Laplacian smoothing
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.5  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.2  2003/08/25 21:51:23  jblazek
! Full version of global grid motion scheme.
!
! Revision 1.1  2003/08/11 21:51:17  jblazek
! Added basic global grid smoothing scheme.
!
!******************************************************************************







