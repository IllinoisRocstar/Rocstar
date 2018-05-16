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
! Purpose: conduct linear interpolation to remesh blocks with inverted cells
!
! Description: none.
!
! Input: region = data of current region
!
! Output: region%levels%grid%xyz = new mesh in regions with inverted cells
!
! Notes: none
!
!******************************************************************************
!
! $Id: RFLO_GridRemesh.F90,v 1.6 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_GridRemesh( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset, &
                            RFLO_GetCellOffset, RFLO_GenerateCoarseGrids, &
                            RFLO_CopyGeometryDummy, RFLO_ExtrapolateGeometry 
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, iter

! ... local variables
  INTEGER :: iLev, iNOff, ijNOff, iCOff, ijCOff
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: ib, ie, ir, jb, je, jr, kb, ke, kr
  INTEGER :: ijk, ijkC, ip1, im1, jp1, jm1, kp1, km1
  INTEGER :: ia1, is1, ja1, js1, ka1, ks1

  REAL(RFREAL) :: rxi, ryi, rzi, rxj, ryj, rzj, rxk, ryk, rzk, rd, p
  REAL(RFREAL) :: rmi, rpi, rmj, rpj, rmk, rpk, volMax
  REAL(RFREAL) :: qmi, qpi, qmj, qpj, qmk, qpk
  REAL(RFREAL) :: rati, ratj, ratk, dsi, dsj, dsk, dli, dlj, dlk, wi, wj, wk
  REAL(RFREAL), POINTER :: xyz(:,:), xyzOrig(:,:), xyzOld(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_GridRemesh',&
  'RFLO_GridRemesh.F90' )

! write to stdout -------------------------------------------------------------

  IF (region%global%verbLevel /= VERBOSE_NONE) THEN
    WRITE(STDOUT,1000) SOLVER_NAME,region%iRegionGlobal
  ENDIF    ! verbLevel

! get dimensions and pointers -------------------------------------------------

  iLev = 1

  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend,jpnbeg,jpnend, &
                                kpnbeg,kpnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  xyz     => region%levels(iLev)%grid%xyz
  xyzOrig => region%levels(iLev)%gridOld%xyz
  xyzOld  => region%levels(iLev)%grid%xyzOld
  p       =  region%global%moveGridWeight

  rd     =  1._RFREAL/3._RFREAL
  volMax = -1.e+30_RFREAL

! determine remeshing direction -----------------------------------------------

  dsi = 0._RFREAL
  DO i=ipnbeg,ipnend-1
    ijk  = IndIJK(i  ,jpnbeg ,kpnbeg ,iNOff,ijNOff)
    ip1  = IndIJK(i+1,jpnbeg ,kpnbeg ,iNOff,ijNOff)

    dsi= dsi+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                   (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                   (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

    ijk  = IndIJK(i  ,jpnbeg ,kpnend ,iNOff,ijNOff)
    ip1  = IndIJK(i+1,jpnbeg ,kpnend ,iNOff,ijNOff)

    dsi= dsi+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                   (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                   (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

    ijk  = IndIJK(i  ,jpnend ,kpnend ,iNOff,ijNOff)
    ip1  = IndIJK(i+1,jpnend ,kpnend ,iNOff,ijNOff)

    dsi= dsi+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                   (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                   (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

    ijk  = IndIJK(i  ,jpnend ,kpnbeg ,iNOff,ijNOff)
    ip1  = IndIJK(i+1,jpnend ,kpnbeg ,iNOff,ijNOff)

    dsi= dsi+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                   (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                   (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )
  ENDDO

  dli = 0._RFREAL

  ijk  = IndIJK(ipnbeg ,jpnbeg ,kpnbeg ,iNOff,ijNOff)
  ip1  = IndIJK(ipnend ,jpnbeg ,kpnbeg ,iNOff,ijNOff)
  dli= dli+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                 (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                 (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

  ijk  = IndIJK(ipnbeg ,jpnbeg ,kpnend ,iNOff,ijNOff)
  ip1  = IndIJK(ipnend ,jpnbeg ,kpnend ,iNOff,ijNOff)
  dli= dli+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                 (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                 (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

  ijk  = IndIJK(ipnbeg ,jpnend ,kpnend ,iNOff,ijNOff)
  ip1  = IndIJK(ipnend ,jpnend ,kpnend ,iNOff,ijNOff)
  dli= dli+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                 (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                 (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

  ijk  = IndIJK(ipnbeg ,jpnend ,kpnbeg ,iNOff,ijNOff)
  ip1  = IndIJK(ipnend ,jpnend ,kpnbeg ,iNOff,ijNOff)
  dli= dli+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                 (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                 (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

! --------------------

  dsj = 0._RFREAL
  DO j=jpnbeg,jpnend-1
    ijk  = IndIJK(ipnbeg ,j    ,kpnbeg ,iNOff,ijNOff)
    ip1  = IndIJK(ipnbeg ,j+1  ,kpnbeg ,iNOff,ijNOff)

    dsj= dsj+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                   (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                   (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

    ijk  = IndIJK(ipnbeg ,j    ,kpnend ,iNOff,ijNOff)
    ip1  = IndIJK(ipnbeg ,j+1  ,kpnend ,iNOff,ijNOff)

    dsj= dsj+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                   (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                   (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

    ijk  = IndIJK(ipnend ,j    ,kpnend ,iNOff,ijNOff)
    ip1  = IndIJK(ipnend ,j+1  ,kpnend ,iNOff,ijNOff)

    dsj= dsj+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                   (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                   (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

    ijk  = IndIJK(ipnend ,j    ,kpnbeg ,iNOff,ijNOff)
    ip1  = IndIJK(ipnend ,j+1  ,kpnbeg ,iNOff,ijNOff)

    dsj= dsj+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                   (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                   (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )
  ENDDO

  dlj = 0._RFREAL

  ijk  = IndIJK(ipnbeg ,jpnbeg ,kpnbeg ,iNOff,ijNOff)
  ip1  = IndIJK(ipnbeg ,jpnend ,kpnbeg ,iNOff,ijNOff)
  dlj= dlj+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                 (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                 (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

  ijk  = IndIJK(ipnbeg ,jpnbeg ,kpnend ,iNOff,ijNOff)
  ip1  = IndIJK(ipnbeg ,jpnend ,kpnend ,iNOff,ijNOff)
  dlj= dlj+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                 (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                 (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

  ijk  = IndIJK(ipnend ,jpnbeg ,kpnend ,iNOff,ijNOff)
  ip1  = IndIJK(ipnend ,jpnend ,kpnend ,iNOff,ijNOff)
  dlj= dlj+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                 (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                 (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

  ijk  = IndIJK(ipnend ,jpnbeg ,kpnbeg ,iNOff,ijNOff)
  ip1  = IndIJK(ipnend ,jpnend ,kpnbeg ,iNOff,ijNOff)
  dlj= dlj+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                 (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                 (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

! -----------------

  dsk = 0._RFREAL
  DO k=kpnbeg,kpnend-1
    ijk  = IndIJK(ipnbeg ,jpnbeg ,k   ,iNOff,ijNOff)
    ip1  = IndIJK(ipnbeg ,jpnbeg ,k+1 ,iNOff,ijNOff)

    dsk= dsk+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                   (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                   (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

    ijk  = IndIJK(ipnend ,jpnbeg ,k   ,iNOff,ijNOff)
    ip1  = IndIJK(ipnend ,jpnbeg ,k+1 ,iNOff,ijNOff)

    dsk= dsk+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                   (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                   (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

    ijk  = IndIJK(ipnend ,jpnend ,k   ,iNOff,ijNOff)
    ip1  = IndIJK(ipnend ,jpnend ,k+1 ,iNOff,ijNOff)

    dsk= dsk+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                   (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                   (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

    ijk  = IndIJK(ipnbeg ,jpnend ,k   ,iNOff,ijNOff)
    ip1  = IndIJK(ipnbeg ,jpnend ,k+1 ,iNOff,ijNOff)

    dsk= dsk+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                   (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                   (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )
  ENDDO

  dlk = 0._RFREAL

  ijk  = IndIJK(ipnbeg ,jpnbeg ,kpnbeg ,iNOff,ijNOff)
  ip1  = IndIJK(ipnbeg ,jpnbeg ,kpnend ,iNOff,ijNOff)
  dlk= dlk+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                 (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                 (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

  ijk  = IndIJK(ipnend ,jpnbeg ,kpnbeg ,iNOff,ijNOff)
  ip1  = IndIJK(ipnend ,jpnbeg ,kpnend ,iNOff,ijNOff)
  dlk= dlk+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                 (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                 (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

  ijk  = IndIJK(ipnend ,jpnend ,kpnbeg ,iNOff,ijNOff)
  ip1  = IndIJK(ipnend ,jpnend ,kpnend ,iNOff,ijNOff)
  dlk= dlk+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                 (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                 (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

  ijk  = IndIJK(ipnend ,jpnend ,kpnbeg ,iNOff,ijNOff)
  ip1  = IndIJK(ipnend ,jpnend ,kpnend ,iNOff,ijNOff)
  dlk= dlk+ SQRT((xyzOrig(XCOORD,ip1)- xyzOrig(XCOORD,ijk))**2 + &
                 (xyzOrig(YCOORD,ip1)- xyzOrig(YCOORD,ijk))**2 + &
                 (xyzOrig(ZCOORD,ip1)- xyzOrig(ZCOORD,ijk))**2 )

! linear interolate for new mesh ---------------------------------------------

  rati = dli/dsi
  ratj = dlj/dsj
  ratk = dlk/dsk
  IF (rati > ratj .AND. rati > ratk .AND. ratj > ratk) THEN
    wi = 1._RFREAL
    wj = 0._RFREAL
    wk = 0._RFREAL
  ELSEIF (rati > ratj .AND. rati > ratk .AND. ratj < ratk) THEN
    wi = 1._RFREAL
    wj = 0._RFREAL
    wk = 0._RFREAL
  ELSEIF (ratj > rati .AND. ratj > ratk .AND. rati > ratk) THEN
    wi = 0._RFREAL
    wj = 1._RFREAL
    wk = 0._RFREAL
  ELSEIF (ratj > rati .AND. ratj > ratk .AND. rati < ratk) THEN
    wi = 0._RFREAL
    wj = 1._RFREAL
    wk = 0._RFREAL
  ELSEIF (ratk > rati .AND. ratk > ratj .AND. rati > ratj) THEN
    wi = 0._RFREAL
    wj = 0._RFREAL
    wk = 1._RFREAL
  ELSEIF (ratk > rati .AND. ratk > ratj .AND. rati < ratj) THEN
    wi = 0._RFREAL
    wj = 0._RFREAL
    wk = 1._RFREAL
  ENDIF

  DO k=kpnbeg+1,kpnend-1
    DO j=jpnbeg+1,jpnend-1
      DO i=ipnbeg+1,ipnend-1
        ijk  = IndIJK(i     ,j  ,k  ,iNOff,ijNOff)
        ip1  = IndIJK(ipnend,j  ,k  ,iNOff,ijNOff)
        im1  = IndIJK(ipnbeg,j  ,k  ,iNOff,ijNOff)
        jp1  = IndIJK(i  ,jpnend,k  ,iNOff,ijNOff)
        jm1  = IndIJK(i  ,jpnbeg,k  ,iNOff,ijNOff)
        kp1  = IndIJK(i  ,j  ,kpnend,iNOff,ijNOff)
        km1  = IndIJK(i  ,j  ,kpnbeg,iNOff,ijNOff)

        ia1  = IndIJK(i+1   ,j  ,k  ,iNOff,ijNOff)
        is1  = IndIJK(i-1   ,j  ,k  ,iNOff,ijNOff)
        ja1  = IndIJK(i  ,j+1   ,k  ,iNOff,ijNOff)
        js1  = IndIJK(i  ,j-1   ,k  ,iNOff,ijNOff)
        ka1  = IndIJK(i  ,j  ,k+1   ,iNOff,ijNOff)
        ks1  = IndIJK(i  ,j  ,k-1   ,iNOff,ijNOff)
 
        qpi= (xyzOrig(XCOORD,ia1)- xyzOrig(XCOORD,ijk))**2 + &
             (xyzOrig(YCOORD,ia1)- xyzOrig(YCOORD,ijk))**2 + &
             (xyzOrig(ZCOORD,ia1)- xyzOrig(ZCOORD,ijk))**2
        qmi= (xyzOrig(XCOORD,is1)- xyzOrig(XCOORD,ijk))**2 + &
             (xyzOrig(YCOORD,is1)- xyzOrig(YCOORD,ijk))**2 + &
             (xyzOrig(ZCOORD,is1)- xyzOrig(ZCOORD,ijk))**2
 
        qpj= (xyzOrig(XCOORD,ja1)- xyzOrig(XCOORD,ijk))**2 + &
             (xyzOrig(YCOORD,ja1)- xyzOrig(YCOORD,ijk))**2 + &
             (xyzOrig(ZCOORD,ja1)- xyzOrig(ZCOORD,ijk))**2
        qmj= (xyzOrig(XCOORD,js1)- xyzOrig(XCOORD,ijk))**2 + &
             (xyzOrig(YCOORD,js1)- xyzOrig(YCOORD,ijk))**2 + &
             (xyzOrig(ZCOORD,js1)- xyzOrig(ZCOORD,ijk))**2
 
        qpk= (xyzOrig(XCOORD,ka1)- xyzOrig(XCOORD,ijk))**2 + &
             (xyzOrig(YCOORD,ka1)- xyzOrig(YCOORD,ijk))**2 + &
             (xyzOrig(ZCOORD,ka1)- xyzOrig(ZCOORD,ijk))**2
        qmk= (xyzOrig(XCOORD,ks1)- xyzOrig(XCOORD,ijk))**2 + &
             (xyzOrig(YCOORD,ks1)- xyzOrig(YCOORD,ijk))**2 + &
             (xyzOrig(ZCOORD,ks1)- xyzOrig(ZCOORD,ijk))**2
       
        qpi = SQRT( qpi )
        qmi = SQRT( qmi )
        qpj = SQRT( qpj )
        qmj = SQRT( qmj )
        qpk = SQRT( qpk )
        qmk = SQRT( qmk )

        rpi = REAL(ipnend-i)
        rmi = REAL(i-ipnbeg)
        rpj = REAL(jpnend-j)
        rmj = REAL(j-jpnbeg)
        rpk = REAL(kpnend-k)
        rmk = REAL(k-kpnbeg)

        rxi = (rpi*xyzOrig(XCOORD,im1) + rmi*xyzOrig(XCOORD,ip1))/(rpi+rmi)+&
              (qpj*xyzOrig(XCOORD,js1) + qmj*xyzOrig(XCOORD,ja1))/(qpj+qmj)+&
              (qpk*xyzOrig(XCOORD,ks1) + qmk*xyzOrig(XCOORD,ka1))/(qpk+qmk)
        ryi = (rpi*xyzOrig(YCOORD,im1) + rmi*xyzOrig(YCOORD,ip1))/(rpi+rmi)+&
              (qpj*xyzOrig(YCOORD,js1) + qmj*xyzOrig(YCOORD,ja1))/(qpj+qmj)+&
              (qpk*xyzOrig(YCOORD,ks1) + qmk*xyzOrig(YCOORD,ka1))/(qpk+qmk)
        rzi = (rpi*xyzOrig(ZCOORD,im1) + rmi*xyzOrig(ZCOORD,ip1))/(rpi+rmi)+&
              (qpj*xyzOrig(ZCOORD,js1) + qmj*xyzOrig(ZCOORD,ja1))/(qpj+qmj)+&
              (qpk*xyzOrig(ZCOORD,ks1) + qmk*xyzOrig(ZCOORD,ka1))/(qpk+qmk)

        rxj = (rpj*xyzOrig(XCOORD,jm1) + rmj*xyzOrig(XCOORD,jp1))/(rpj+rmj)+&
              (qpi*xyzOrig(XCOORD,is1) + qmi*xyzOrig(XCOORD,ia1))/(qpi+qmi)+&
              (qpk*xyzOrig(XCOORD,ks1) + qmk*xyzOrig(XCOORD,ka1))/(qpk+qmk)
        ryj = (rpj*xyzOrig(YCOORD,jm1) + rmj*xyzOrig(YCOORD,jp1))/(rpj+rmj)+&
              (qpi*xyzOrig(YCOORD,is1) + qmi*xyzOrig(YCOORD,ia1))/(qpi+qmi)+&
              (qpk*xyzOrig(YCOORD,ks1) + qmk*xyzOrig(YCOORD,ka1))/(qpk+qmk)
        rzj = (rpj*xyzOrig(ZCOORD,jm1) + rmj*xyzOrig(ZCOORD,jp1))/(rpj+rmj)+&
              (qpi*xyzOrig(ZCOORD,is1) + qmi*xyzOrig(ZCOORD,ia1))/(qpi+qmi)+&
              (qpk*xyzOrig(ZCOORD,ks1) + qmk*xyzOrig(ZCOORD,ka1))/(qpk+qmk)

        rxk = (rpk*xyzOrig(XCOORD,km1) + rmk*xyzOrig(XCOORD,kp1))/(rpk+rmk)+&
              (qpi*xyzOrig(XCOORD,is1) + qmi*xyzOrig(XCOORD,ia1))/(qpi+qmi)+&
              (qpj*xyzOrig(XCOORD,js1) + qmj*xyzOrig(XCOORD,ja1))/(qpj+qmj)
        ryk = (rpk*xyzOrig(YCOORD,km1) + rmk*xyzOrig(YCOORD,kp1))/(rpk+rmk)+&
              (qpi*xyzOrig(YCOORD,is1) + qmi*xyzOrig(YCOORD,ia1))/(qpi+qmi)+&
              (qpj*xyzOrig(YCOORD,js1) + qmj*xyzOrig(YCOORD,ja1))/(qpj+qmj)
        rzk = (rpk*xyzOrig(ZCOORD,km1) + rmk*xyzOrig(ZCOORD,kp1))/(rpk+rmk)+&
              (qpi*xyzOrig(ZCOORD,is1) + qmi*xyzOrig(ZCOORD,ia1))/(qpi+qmi)+&
              (qpj*xyzOrig(ZCOORD,js1) + qmj*xyzOrig(ZCOORD,ja1))/(qpj+qmj)

        xyz(XCOORD,ijk) = rd*(wi*rxi+wj*rxj+wk*rxk)/(wi+wj+wk)
        xyz(YCOORD,ijk) = rd*(wi*ryi+wj*ryj+wk*ryk)/(wi+wj+wk)
        xyz(ZCOORD,ijk) = rd*(wi*rzi+wj*rzj+wk*rzk)/(wi+wj+wk)

      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! update coarse grids and dummy cells

  CALL RFLO_GenerateCoarseGrids( region )   ! coarsen finest grid
  CALL RFLO_CopyGeometryDummy( region )     ! copy to dummy nodes
  CALL RFLO_ExtrapolateGeometry( region )   ! extrapolate

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

1000 FORMAT(A,1X,'Remesh region: ',I6)

END SUBROUTINE RFLO_GridRemesh

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_GridRemesh.F90,v $
! Revision 1.6  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/06/01 07:12:20  wasistho
! added 4 adjacent neighbours in averaging process
!
! Revision 1.3  2005/05/28 05:42:50  wasistho
! cosmetics
!
! Revision 1.2  2005/05/27 08:39:36  wasistho
! write to stdout when remeshing
!
! Revision 1.1  2005/05/27 01:53:41  wasistho
! added rflo_gridremesh
!
!
!******************************************************************************







