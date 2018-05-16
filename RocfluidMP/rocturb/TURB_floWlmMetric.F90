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
! Purpose: Compute metric coefficients needed in wlm computations.
!
! Description: Metric quantities stored in data structure are,
!              - transformation matrix coeffs. from body fitted to Cartesian
!              - grid spacing in body fitted coord. in wall parallel 
!                directions at first (or user selected) points from the wall
!              - wall normal distance to first points from the wall.
!
! Input: region  = data of current region.
!        patch   = current patch.
!
! Output: metric quantities described above.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_floWlmMetric.F90,v 1.5 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloWlmMetric( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices,RFLO_GetPatchDirection, &
                            RFLO_GetNodeOffset, RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k, l

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: ilev, lbound, nRef, isub, jsub, ksub
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: inode, jnode, knode, idir, jdir, kdir, ir, jr, kr
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ix,jx,kx, iz,jz,kz, ixz,jxz,kxz
  INTEGER :: ijkN1B, ijkN1, ijkN1x, ijkN1z, ijkVal
  INTEGER :: ijkN2B, ijkN2, ijkN2x, ijkN2z
  INTEGER :: ijkN1Bx, ijkN1Bz, ijkN1Bxz, ijkN2Bx, ijkN2Bz, ijkN2Bxz
  INTEGER :: ijkN0B, ijkN0Bx, ijkN0Bz, ijkN0Bxz

  REAL(RFREAL)          :: sxi(3), set(3), szt(3), fCenb(3), fCene(3), xyz(3,4)
  REAL(RFREAL)          :: modSxi, modSet, modSzt, sgn
  REAL(RFREAL), POINTER :: sfxi(:,:), sfet(:,:), sfzt(:,:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floWlmMetric.F90,v $ $Revision: 1.5 $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloWlmMetric',&
  'TURB_floWlmMetric.F90' )

! get dimensions, parameters and pointers -------------------------------------

  ilev    =  region%currLevel
  lbound  =  patch%lbound
  nRef    =  patch%valBola%switches(WLM_INPUT_REFPOINT) - 1

  CALL RFLO_GetPatchIndices( region,patch,ilev,ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetNodeOffset( region,ilev,iNOff,ijNOff )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )

  sgn  = +1._RFREAL
  inode = 0
  jnode = 0
  knode = 0
  IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
    sgn   = -1._RFREAL
    inode = -idir
    jnode = -jdir
    knode = -kdir
  ENDIF
  
  IF (lbound==1 .OR. lbound==2) THEN
    ir = idir*nRef
    jr = 0
    kr = 0
    ix = 0
    jx = 0
    kx = 1
    iz = 0
    jz = 1
    kz = 0
    sfxi => region%levels(ilev)%grid%sk
    sfet => region%levels(ilev)%grid%si
    sfzt => region%levels(ilev)%grid%sj
  ELSEIF (lbound==3 .OR. lbound==4) THEN
    ir = 0
    jr = jdir*nRef
    kr = 0
    ix = 1
    jx = 0
    kx = 0
    iz = 0
    jz = 0
    kz = 1
    sfxi => region%levels(ilev)%grid%si
    sfet => region%levels(ilev)%grid%sj
    sfzt => region%levels(ilev)%grid%sk
  ELSEIF (lbound==5 .OR. lbound==6) THEN
    ir = 0
    jr = 0
    kr = kdir*nRef
    ix = 0
    jx = 1
    kx = 0
    iz = 1
    jz = 0
    kz = 0
    sfxi => region%levels(ilev)%grid%sj
    sfet => region%levels(ilev)%grid%sk
    sfzt => region%levels(ilev)%grid%si
  ENDIF
  ixz = ix+iz
  jxz = jx+jz
  kxz = kx+kz

! Get mapping coefficients : define xi,eta,zeta normal vectors from i,j,k 
! face vectors at cell center. The fact that i,j,k face vectors points to 
! negative i,j,k directions are taken into account. Next, get finite 
! difference spacings

  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,iend

! ----- get mapping coefficients

        isub = i+ir
        jsub = j+jr
        ksub = k+kr

        ijkN1B = IndIJK(isub+inode,jsub+jnode,ksub+knode,iNOff,ijNOff) !bnd nodes
        ijkN1  = IndIJK(isub      ,jsub      ,ksub      ,iNOff,ijNOff)  
        ijkN1x = IndIJK(isub+ix   ,jsub+jx   ,ksub+kx   ,iNOff,ijNOff)  
        ijkN1z = IndIJK(isub+iz   ,jsub+jz   ,ksub+kz   ,iNOff,ijNOff)  

        ijkN2B = &
        IndIJK(isub+inode+idir,jsub+jnode+jdir,ksub+knode+kdir,iNOff,ijNOff) 

        IF (lbound==1 .OR. lbound==2) THEN
          ijkVal = IndIJ(j-jbeg ,k-kbeg ,jend-jbeg+1)
        ELSEIF (lbound==3 .OR. lbound==4) THEN
          ijkVal = IndIJ(k-kbeg ,i-ibeg ,kend-kbeg+1)
        ELSEIF (lbound==5 .OR. lbound==6) THEN
          ijkVal = IndIJ(i-ibeg ,j-jbeg ,iend-ibeg+1)
        ENDIF

        DO l=XCOORD,ZCOORD
          sxi(l) = -0.5_RFREAL*(sfxi(l,ijkN1)+sfxi(l,ijkN1x))
          set(l) = -sgn*0.5_RFREAL*(sfet(l,ijkN1B)+sfet(l,ijkN2B))
          szt(l) = -0.5_RFREAL*(sfzt(l,ijkN1)+sfzt(l,ijkN1z))
        ENDDO

        modSxi = SQRT( sxi(1)*sxi(1)+sxi(2)*sxi(2)+sxi(3)*sxi(3) )
        modSet = SQRT( set(1)*set(1)+set(2)*set(2)+set(3)*set(3) )
        modSzt = SQRT( szt(1)*szt(1)+szt(2)*szt(2)+szt(3)*szt(3) )

        patch%valBola%vals(ijkVal,WLM_VALS_XIX) = sxi(1)/modSxi
        patch%valBola%vals(ijkVal,WLM_VALS_ETX) = set(1)/modSet
        patch%valBola%vals(ijkVal,WLM_VALS_ZTX) = szt(1)/modSzt

        patch%valBola%vals(ijkVal,WLM_VALS_XIY) = sxi(2)/modSxi
        patch%valBola%vals(ijkVal,WLM_VALS_ETY) = set(2)/modSet
        patch%valBola%vals(ijkVal,WLM_VALS_ZTY) = szt(2)/modSzt

        patch%valBola%vals(ijkVal,WLM_VALS_XIZ) = sxi(3)/modSxi
        patch%valBola%vals(ijkVal,WLM_VALS_ETZ) = set(3)/modSet
        patch%valBola%vals(ijkVal,WLM_VALS_ZTZ) = szt(3)/modSzt

! ----- get finite difference spacings dXi, dZeta, respectively

        isub = i+ir+inode
        jsub = j+jr+jnode
        ksub = k+kr+knode

        ijkN1Bx = IndIJK(isub+ix  ,jsub+jx  ,ksub+kx  ,iNOff,ijNOff) 
        ijkN1Bz = IndIJK(isub+iz  ,jsub+jz  ,ksub+kz  ,iNOff,ijNOff) 
        ijkN1Bxz= IndIJK(isub+ixz ,jsub+jxz ,ksub+kxz ,iNOff,ijNOff)

        ijkN2Bx = &
        IndIJK(isub+idir+ix  ,jsub+jdir+jx  ,ksub+kdir+kx  ,iNOff,ijNOff)  
        ijkN2Bz = &
        IndIJK(isub+idir+iz  ,jsub+jdir+jz  ,ksub+kdir+kz  ,iNOff,ijNOff)  
        ijkN2Bxz= &
        IndIJK(isub+idir+ixz ,jsub+jdir+jxz ,ksub+kdir+kxz ,iNOff,ijNOff)  

        DO l=XCOORD,ZCOORD
          xyz(l,1) = region%levels(ilev)%grid%xyz(l,ijkN1B)
          xyz(l,2) = region%levels(ilev)%grid%xyz(l,ijkN2B)
          xyz(l,3) = region%levels(ilev)%grid%xyz(l,ijkN1Bz)
          xyz(l,4) = region%levels(ilev)%grid%xyz(l,ijkN2Bz)
          fCenb(l) = 0.25_RFREAL*(xyz(l,1) + xyz(l,2) + xyz(l,3) + xyz(l,4))

          xyz(l,1) = region%levels(ilev)%grid%xyz(l,ijkN1Bx)
          xyz(l,2) = region%levels(ilev)%grid%xyz(l,ijkN2Bx)
          xyz(l,3) = region%levels(ilev)%grid%xyz(l,ijkN1Bxz)
          xyz(l,4) = region%levels(ilev)%grid%xyz(l,ijkN2Bxz)
          fCene(l) = 0.25_RFREAL*(xyz(l,1) + xyz(l,2) + xyz(l,3) + xyz(l,4))
        ENDDO
        patch%valBola%vals(ijkVal,WLM_VALS_DXI) = SQRT( (fCene(1)-fCenb(1))**2+ &
                                                        (fCene(2)-fCenb(2))**2+ &
                                                        (fCene(3)-fCenb(3))**2 )

        DO l=XCOORD,ZCOORD
          xyz(l,1) = region%levels(ilev)%grid%xyz(l,ijkN1B)
          xyz(l,2) = region%levels(ilev)%grid%xyz(l,ijkN2B)
          xyz(l,3) = region%levels(ilev)%grid%xyz(l,ijkN1Bx)
          xyz(l,4) = region%levels(ilev)%grid%xyz(l,ijkN2Bx)
          fCenb(l) = 0.25_RFREAL*(xyz(l,1) + xyz(l,2) + xyz(l,3) + xyz(l,4))

          xyz(l,1) = region%levels(ilev)%grid%xyz(l,ijkN1Bz)
          xyz(l,2) = region%levels(ilev)%grid%xyz(l,ijkN2Bz)
          xyz(l,3) = region%levels(ilev)%grid%xyz(l,ijkN1Bxz)
          xyz(l,4) = region%levels(ilev)%grid%xyz(l,ijkN2Bxz)
          fCene(l) = 0.25_RFREAL*(xyz(l,1) + xyz(l,2) + xyz(l,3) + xyz(l,4))
        ENDDO
        patch%valBola%vals(ijkVal,WLM_VALS_DZT) = SQRT( (fCene(1)-fCenb(1))**2+ &
                                                        (fCene(2)-fCenb(2))**2+ &
                                                        (fCene(3)-fCenb(3))**2 )

! ----- get wall distance wDist

        isub = i+inode
        jsub = j+jnode
        ksub = k+knode

        ijkN0B   = IndIJK(isub     ,jsub     ,ksub     ,iNOff,ijNOff)
        ijkN0Bx  = IndIJK(isub+ix  ,jsub+jx  ,ksub+kx  ,iNOff,ijNOff)
        ijkN0Bz  = IndIJK(isub+iz  ,jsub+jz  ,ksub+kz  ,iNOff,ijNOff)
        ijkN0Bxz = IndIJK(isub+ixz ,jsub+jxz ,ksub+kxz ,iNOff,ijNOff)

        DO l=XCOORD,ZCOORD
          xyz(l,1) = region%levels(ilev)%grid%xyz(l,ijkN1B)
          xyz(l,2) = region%levels(ilev)%grid%xyz(l,ijkN1Bx)
          xyz(l,3) = region%levels(ilev)%grid%xyz(l,ijkN1Bz)
          xyz(l,4) = region%levels(ilev)%grid%xyz(l,ijkN1Bxz)
          fCenb(l) = 0.25_RFREAL*(xyz(l,1) + xyz(l,2) + xyz(l,3) + xyz(l,4))
          xyz(l,1) = region%levels(ilev)%grid%xyz(l,ijkN2B)
          xyz(l,2) = region%levels(ilev)%grid%xyz(l,ijkN2Bx)
          xyz(l,3) = region%levels(ilev)%grid%xyz(l,ijkN2Bz)
          xyz(l,4) = region%levels(ilev)%grid%xyz(l,ijkN2Bxz)
          fCene(l) = 0.25_RFREAL*(xyz(l,1) + xyz(l,2) + xyz(l,3) + xyz(l,4)) 
          fCene(l) = 0.50_RFREAL*(fCene(l)+fCenb(l))
          xyz(l,1) = region%levels(ilev)%grid%xyz(l,ijkN0B)
          xyz(l,2) = region%levels(ilev)%grid%xyz(l,ijkN0Bx)
          xyz(l,3) = region%levels(ilev)%grid%xyz(l,ijkN0Bz)
          xyz(l,4) = region%levels(ilev)%grid%xyz(l,ijkN0Bxz)
          fCenb(l) = 0.25_RFREAL*(xyz(l,1) + xyz(l,2) + xyz(l,3) + xyz(l,4))
          set(l)   = -sgn*0.5_RFREAL*(sfet(l,ijkN1B)+sfet(l,ijkN2B))
        ENDDO
        modSet = SQRT( set(1)*set(1)+set(2)*set(2)+set(3)*set(3) )
        patch%valBola%vals(ijkVal,WLM_VALS_WDIST) = &
                      ABS( (fCene(1)-fCenb(1))*set(1)+ &
                           (fCene(2)-fCenb(2))*set(2)+ &
                           (fCene(3)-fCenb(3))*set(3) )/modSet

      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloWlmMetric

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_floWlmMetric.F90,v $
! Revision 1.5  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/03/19 02:53:16  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:46  wasistho
! changed turb nomenclature
!
! Revision 1.3  2004/02/14 03:42:33  wasistho
! added new WLM parameter: reference point
!
!
!******************************************************************************







