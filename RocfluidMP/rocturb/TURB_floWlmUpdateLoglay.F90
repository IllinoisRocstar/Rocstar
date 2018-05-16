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
! Purpose: Model wall stresses based on logarithmic assumption at constant p.
!
! Description: Wall stresses are derived from modeling terms in BL equations.
!              Only viscous term is modeled based on log layer assumption at 
!              zero pressure gradient. BL convective term, pressure gradient
!              and time derivative term are not considered. Surface roughnes
!              is however taken into account. The method is inspired by paper
!              of Hoffman and Benocci, "Approximate wall BC for LES", 5th
!              Advances in Turbulence, Siena, Italy 1994. Surface roughness
!              were not considered in their model.
!
! Input: region  = data of current region.
!        patch   = current patch.
!
! Output: total and wall parallel components of wall stresses in body fitted 
!         coordinate.
!
! Notes: the pressure gradient and unsteady terms are modeled if BNDLAY 
!        (boundary layer) model is selected and added to the model viscous 
!        term computed here.
!
!******************************************************************************
!
! $Id: TURB_floWlmUpdateLoglay.F90,v 1.5 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloWlmUpdateLoglay( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: ilev, lbound, idir, jdir, kdir, iRef, jRef, kRef, nRef
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iCOff, ijCOff
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag
  INTEGER :: ijkC, ijkCip, ijkCim,ijkCjp, ijkCjm,ijkCkp, ijkCkm, ijkVal

  REAL(RFREAL)              :: dens, wdist, rVonk, consB, ksize, abVel, utau
  REAL(RFREAL)              :: yPls, ksPls, ksPlsMax, func, dfunc, tauWall
  REAL(RFREAL)              :: minUt, maxUt, relaxCoef 
  REAL(RFREAL), POINTER     :: cv(:,:), dv(:,:), tv(:,:), vals(:,:)
  REAL(RFREAL), ALLOCATABLE :: mul(:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floWlmUpdateLoglay.F90,v $ $Revision: 1.5 $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloWlmUpdateLoglay',&
  'TURB_floWlmUpdateLoglay.F90' )

! get dimensions and parameters -----------------------------------------------

  ilev   =  region%currLevel
  lbound =  patch%lbound
  nRef   =  patch%valBola%switches(WLM_INPUT_REFPOINT) - 1

  CALL RFLO_GetPatchIndices( region,patch,ilev,ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )

  relaxCoef = 0.2_RFREAL   ! coef. for utau iterat. in log-law (usr inp. later)

! get pointers and allocate temporary workspace -------------------------------

  cv    => region%levels(ilev)%mixt%cv
  dv    => region%levels(ilev)%mixt%dv
  tv    => region%levels(ilev)%mixt%tv
  vals  => patch%valBola%vals 

  n1    = ABS(patch%l1end-patch%l1beg)
  n2    = ABS(patch%l2end-patch%l2beg)
  iOff  = n1 + 1
  ijBeg = IndIJ( 0, 0,iOff)
  ijEnd = IndIJ(n1,n2,iOff)

  ALLOCATE( mul(ijBeg:ijEnd) )

! begin and end indices for i, j and k directions

! IF (lbound==1 .OR. lbound==2) THEN
!   indxb = jbeg
!   indxe = jend
!   jndxb = kbeg
!   jndxe = kend
! ELSEIF (lbound==3 .OR. lbound==4) THEN
!   indxb = kbeg
!   indxe = kend
!   jndxb = ibeg
!   jndxe = iend
! ELSE
!   indxb = ibeg
!   indxe = iend
!   jndxb = jbeg
!   jndxe = jend
! ENDIF

! Note that indx correspond to zeta-coordinate and jndx to xi.
! Get p gradient in xi and zeta direction and velocities depending on lbound.

  IF (lbound==1 .OR. lbound==2) THEN

! - j-dir == indx-dir == zeta-dir;  k-dir == jndx-dir == xi-dir
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend

          iRef  = i+idir*nRef
          ijkC  = IndIJK(iRef   ,j      ,k      ,iCOff,ijCOff)
          ijkCjp= IndIJK(iRef   ,j+1    ,k      ,iCOff,ijCOff)
          ijkCjm= IndIJK(iRef   ,j-1    ,k      ,iCOff,ijCOff)
          ijkCkp= IndIJK(iRef   ,j      ,k+1    ,iCOff,ijCOff)
          ijkCkm= IndIJK(iRef   ,j      ,k-1    ,iCOff,ijCOff)
          ijkVal= IndIJ(j-jbeg ,k-kbeg ,jend-jbeg+1)

          vals(ijkVal,WLM_VALS_DPDXI) = &
                          (dv(DV_MIXT_PRES,ijkCkp)-dv(DV_MIXT_PRES,ijkCkm))* &
                          0.5_RFREAL/vals(ijkVal,WLM_VALS_DXI)
          vals(ijkVal,WLM_VALS_DPDZT) = & 
                          (dv(DV_MIXT_PRES,ijkCjp)-dv(DV_MIXT_PRES,ijkCjm))* &
                          0.5_RFREAL/vals(ijkVal,WLM_VALS_DZT)
          vals(ijkVal,WLM_VALS_XIV) = &
                          dv(DV_MIXT_UVEL,ijkC)*vals(ijkVal,WLM_VALS_XIX) + &
                          dv(DV_MIXT_VVEL,ijkC)*vals(ijkVal,WLM_VALS_XIY) + &
                          dv(DV_MIXT_WVEL,ijkC)*vals(ijkVal,WLM_VALS_XIZ) 
          vals(ijkVal,WLM_VALS_ETV) = &
                          dv(DV_MIXT_UVEL,ijkC)*vals(ijkVal,WLM_VALS_ETX) + &
                          dv(DV_MIXT_VVEL,ijkC)*vals(ijkVal,WLM_VALS_ETY) + &
                          dv(DV_MIXT_WVEL,ijkC)*vals(ijkVal,WLM_VALS_ETZ) 
          vals(ijkVal,WLM_VALS_ZTV) = &
                          dv(DV_MIXT_UVEL,ijkC)*vals(ijkVal,WLM_VALS_ZTX) + &
                          dv(DV_MIXT_VVEL,ijkC)*vals(ijkVal,WLM_VALS_ZTY) + &
                          dv(DV_MIXT_WVEL,ijkC)*vals(ijkVal,WLM_VALS_ZTZ) 
          vals(ijkVal,WLM_VALS_DENS) = cv(CV_MIXT_DENS,ijkC)
          mul(ijkVal)   = tv(TV_MIXT_MUEL,ijkC)
        ENDDO
      ENDDO
    ENDDO

  ELSEIF (lbound==3 .OR. lbound==4) THEN

! - k-dir == indx-dir == zeta-dir;  i-dir == jndx-dir == xi-dir
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend

          jRef  = j+jdir*nRef
          ijkC  = IndIJK(i      ,jRef   ,k      ,iCOff,ijCOff)
          ijkCip= IndIJK(i+1    ,jRef   ,k      ,iCOff,ijCOff)
          ijkCim= IndIJK(i-1    ,jRef   ,k      ,iCOff,ijCOff)
          ijkCkp= IndIJK(i      ,jRef   ,k+1    ,iCOff,ijCOff)
          ijkCkm= IndIJK(i      ,jRef   ,k-1    ,iCOff,ijCOff)
          ijkVal= IndIJ(k-kbeg ,i-ibeg ,kend-kbeg+1)

          vals(ijkVal,WLM_VALS_DPDXI) = &
                          (dv(DV_MIXT_PRES,ijkCip)-dv(DV_MIXT_PRES,ijkCim))* &
                          0.5_RFREAL/vals(ijkVal,WLM_VALS_DXI)
          vals(ijkVal,WLM_VALS_DPDZT) = & 
                          (dv(DV_MIXT_PRES,ijkCkp)-dv(DV_MIXT_PRES,ijkCkm))* &
                          0.5_RFREAL/vals(ijkVal,WLM_VALS_DZT)
          vals(ijkVal,WLM_VALS_XIV) = &
                          dv(DV_MIXT_UVEL,ijkC)*vals(ijkVal,WLM_VALS_XIX) + &
                          dv(DV_MIXT_VVEL,ijkC)*vals(ijkVal,WLM_VALS_XIY) + &
                          dv(DV_MIXT_WVEL,ijkC)*vals(ijkVal,WLM_VALS_XIZ) 
          vals(ijkVal,WLM_VALS_ETV) = &
                          dv(DV_MIXT_UVEL,ijkC)*vals(ijkVal,WLM_VALS_ETX) + &
                          dv(DV_MIXT_VVEL,ijkC)*vals(ijkVal,WLM_VALS_ETY) + &
                          dv(DV_MIXT_WVEL,ijkC)*vals(ijkVal,WLM_VALS_ETZ) 
          vals(ijkVal,WLM_VALS_ZTV) = &
                          dv(DV_MIXT_UVEL,ijkC)*vals(ijkVal,WLM_VALS_ZTX) + &
                          dv(DV_MIXT_VVEL,ijkC)*vals(ijkVal,WLM_VALS_ZTY) + &
                          dv(DV_MIXT_WVEL,ijkC)*vals(ijkVal,WLM_VALS_ZTZ) 
          vals(ijkVal,WLM_VALS_DENS) = cv(CV_MIXT_DENS,ijkC)
          mul(ijkVal)   = tv(TV_MIXT_MUEL,ijkC)
        ENDDO
      ENDDO
    ENDDO

  ELSE

! - i-dir == indx-dir == zeta-dir;  j-dir == jndx-dir == xi-dir
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend

          kRef  = k+kdir*nRef
          ijkC  = IndIJK(i      ,j      ,kRef   ,iCOff,ijCOff)
          ijkCip= IndIJK(i+1    ,j      ,kRef   ,iCOff,ijCOff)
          ijkCim= IndIJK(i-1    ,j      ,kRef   ,iCOff,ijCOff)
          ijkCjp= IndIJK(i      ,j+1    ,kRef   ,iCOff,ijCOff)
          ijkCjm= IndIJK(i      ,j-1    ,kRef   ,iCOff,ijCOff)
          ijkVal= IndIJ(i-ibeg ,j-jbeg ,iend-ibeg+1)

          vals(ijkVal,WLM_VALS_DPDXI) = &
                          (dv(DV_MIXT_PRES,ijkCjp)-dv(DV_MIXT_PRES,ijkCjm))* &
                          0.5_RFREAL/vals(ijkVal,WLM_VALS_DXI)
          vals(ijkVal,WLM_VALS_DPDZT) = & 
                          (dv(DV_MIXT_PRES,ijkCip)-dv(DV_MIXT_PRES,ijkCim))* &
                          0.5_RFREAL/vals(ijkVal,WLM_VALS_DZT)
          vals(ijkVal,WLM_VALS_XIV) = &
                          dv(DV_MIXT_UVEL,ijkC)*vals(ijkVal,WLM_VALS_XIX) + &
                          dv(DV_MIXT_VVEL,ijkC)*vals(ijkVal,WLM_VALS_XIY) + &
                          dv(DV_MIXT_WVEL,ijkC)*vals(ijkVal,WLM_VALS_XIZ) 
          vals(ijkVal,WLM_VALS_ETV) = &
                          dv(DV_MIXT_UVEL,ijkC)*vals(ijkVal,WLM_VALS_ETX) + &
                          dv(DV_MIXT_VVEL,ijkC)*vals(ijkVal,WLM_VALS_ETY) + &
                          dv(DV_MIXT_WVEL,ijkC)*vals(ijkVal,WLM_VALS_ETZ) 
          vals(ijkVal,WLM_VALS_ZTV) = &
                          dv(DV_MIXT_UVEL,ijkC)*vals(ijkVal,WLM_VALS_ZTX) + &
                          dv(DV_MIXT_VVEL,ijkC)*vals(ijkVal,WLM_VALS_ZTY) + &
                          dv(DV_MIXT_WVEL,ijkC)*vals(ijkVal,WLM_VALS_ZTZ) 
          vals(ijkVal,WLM_VALS_DENS) = cv(CV_MIXT_DENS,ijkC)
          mul(ijkVal)   = tv(TV_MIXT_MUEL,ijkC)
        ENDDO
      ENDDO
    ENDDO

  ENDIF

! compute tau-wall from total contributions

  rVonk    = 2.5_RFREAL
  consB    = 8.5_RFREAL
  ksPlsMax = 70._RFREAL
  vals(:,WLM_VALS_TAUUX:WLM_VALS_TAUWZ) = 0._RFREAL

  DO ijkVal=ijBeg,ijEnd

    wdist = vals(ijkVal,WLM_VALS_WDIST)
    ksize = vals(ijkVal,WLM_VALS_ROUGH)
    dens  = vals(ijkVal,WLM_VALS_DENS)
    abVel = SQRT( vals(ijkVal,WLM_VALS_XIV)**2+vals(ijkVal,WLM_VALS_ZTV)**2 )

    IF (ksize > REAL_SMALL) THEN

! --- surface roughness may be large enough to be taken into account

      utau  = rVonk*LOG( wdist/ksize ) + consB
      utau  = abVel/utau   

! --- utau at LHS above is large for high ksize, but turns to negative when
!     ksize becomes even larger. Theoretically there is maximum value of utau, 
!     i.e. when ksPls reaches value 70. This contrain is applied here.

      maxUt = ksPlsMax*mul(ijkVal)/(ksize*dens) 
      utau  = MIN( utau, maxUt )          ! for large ks but utau still > 0.
      IF (utau < 0._RFREAL) utau = maxUt  ! for large ks that utau < 0.

      ksPls = utau*ksize*dens/mul(ijkVal)
      IF (ksPls >= 5._RFREAL) THEN
!wlmCheckprobe-----------------------------------------
!        yPls  = ksPls*wdist/ksize
!        write(*,*)1,ijkVal,ksPls,yPls,utau,utau/abVel
!------------------------------------------------------
        GOTO 30
      ELSE
        GOTO 20
      ENDIF
    ENDIF

20  CONTINUE

! - surface hidraulically smooth, find utau using Newton-Raphson iteration

!    minUt = 1.E-8_RFREAL*abVel   ! constrain for low utau (near separation)
!    utau  = MAX( minUt,vals(ijkVal,WLM_VALS_UTAU) ) ! initial estimate

!    func = 1._RFREAL
!    DO WHILE (ABS( func ) > STOP_VALUE )
!      yPls  = utau*wdist*dens/mul(ijkVal)
!      func  = rVonk*LOG( yPls ) + 5.5_RFREAL - abVel/utau
!      dfunc = (rVonk*yPls + abVel)/(utau*utau)
!      utau  = MAX( minUt , (utau  - relaxCoef*func/dfunc) )
!wlmCheckprobe---------------------
!      write(*,*) ijkVal,func,yPls
!----------------------------------
!    ENDDO 

! - or using rough estimate (replaced by Dean`s correlation later)
    
    utau = 0.05_RFREAL*abVel
    yPls = utau*wdist*dens/mul(ijkVal)
!wlmCheckprobe---------------------
!    write(*,*) ijkVal,utau,wdist,dens,mul(ijkVal),yPls
!----------------------------------

! - check roughness height in wall unit

    ksPls = utau*ksize*dens/mul(ijkVal)
!wlmCheckprobe-------------------------------------
!    write(*,*)2,ijkVal,ksPls,yPls,utau,utau/abVel
!--------------------------------------------------
30  CONTINUE

! - wall stress in xi and zeta direction and store utau for next stage

    tauWall = dens*utau*utau
    vals(ijkVal,WLM_VALS_TAUUY) = vals(ijkVal,WLM_VALS_XIV)/abVel*tauWall
    vals(ijkVal,WLM_VALS_TAUWY) = vals(ijkVal,WLM_VALS_ZTV)/abVel*tauWall

    vals(ijkVal,WLM_VALS_TAUVX) = vals(ijkVal,WLM_VALS_TAUUY)     
    vals(ijkVal,WLM_VALS_TAUVZ) = vals(ijkVal,WLM_VALS_TAUWY)      

    vals(ijkVal,WLM_VALS_UTAU)  = utau

! - store total wall stress in heat-flux array for heat transfer comp. later
    vals(ijkVal,WLM_VALS_HFLUX) = tauWall
!checkprobe--------------------------------------------------------------------
!    write(*,*) region%procId,patch%lbound,ijkVal,vals(ijkVal,WLM_VALS_TAUUY), &
!               vals(ijkVal,WLM_VALS_TAUWY),utau/abVel
!------------------------------------------------------------------------------   

  ENDDO       ! ijkVal

! deallocate temporary arrays

  DEALLOCATE( mul )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloWlmUpdateLoglay

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_floWlmUpdateLoglay.F90,v $
! Revision 1.5  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/03/19 02:53:24  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:46  wasistho
! changed turb nomenclature
!
! Revision 1.4  2004/02/14 03:43:30  wasistho
! added new WLM parameter: reference point
!
!
!******************************************************************************







