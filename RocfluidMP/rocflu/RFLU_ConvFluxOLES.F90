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
! Purpose: Compute convective fluxes using optimal LES approach.
!
! Description: None.
!
! Input: region = data of current region.
!
! Output: region%mixt%rhs = convective fluxes added to the residual.
!
! Notes: 
!   1. No boundary fluxes because present implementation works only for 
!      isotropic turbulence
!   2. No grid speeds implemented.
!
!******************************************************************************
!
! $Id: RFLU_ConvFluxOLES.F90,v 1.7 2008/12/06 08:44:29 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ConvFluxOLES(region)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid  
  USE ModError
  USE ModParameters
  
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i,ifc,ifcp,iPatch,j,k

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: c1,c2,ic1l,ic1g,ic2l,ic2g,nCells,netMassFluxCntr
  INTEGER, DIMENSION(:), POINTER :: f2fp
  INTEGER, DIMENSION(:,:), POINTER :: fs,f2c
  REAL(RFREAL) :: el,er,irl,irr,netMassFlux,nm,nx,ny,nz,ql,qr,rl,rm,rr,vcont, & 
                  ul,ur,vl,vr,wl,wr,pl,pr
  REAL(RFREAL) :: fc(5),fl(3),fq(3),v1(3),v2(3)
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,dv,rhs,fn
  REAL(RFREAL), DIMENSION(:,:,:,:), POINTER :: wtl
  REAL(RFREAL), DIMENSION(:,:,:,:,:,:), POINTER :: wtq
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ConvFluxOLES.F90,v $ $Revision: 1.7 $'

  global => region%global

  CALL RegisterFunction(global,'RFLU_ConvFluxOLES',&
  'RFLU_ConvFluxOLES.F90')

! get dimensions and pointers -------------------------------------------------

  cv  => region%mixt%cv
  dv  => region%mixt%dv
  rhs => region%mixt%rhs
  
  f2c  => region%grid%f2c
  fn   => region%grid%fn

  fs   => region%grid%fsOLES 
  f2fp => region%grid%f2fpOLES
  wtl  => region%grid%wtLinOLES
  wtq  => region%grid%wtQuadOLES

  netMassFluxCntr = 0
  netMassFlux     = 0.0_RFREAL

! stationary grid -------------------------------------------------------------
! flux (except through boundary)

  nCells = SIZE(fs,1) ! Get size of optimal LES stencil

  DO ifc = 1,region%grid%nFaces
    ifcp = f2fp(ifc) ! Get prototype face
  
    c1 = f2c(1,ifc)
    c2 = f2c(2,ifc)
                        
    nx = fn(XCOORD,ifc)
    ny = fn(YCOORD,ifc)
    nz = fn(ZCOORD,ifc)
    nm = fn(XYZMAG,ifc)
  
    rl  = cv(CV_MIXT_DENS,c1)
    irl = 1.0_RFREAL/rl
    
    ul  = cv(CV_MIXT_XMOM,c1)*irl
    vl  = cv(CV_MIXT_YMOM,c1)*irl
    wl  = cv(CV_MIXT_ZMOM,c1)*irl
    el  = cv(CV_MIXT_ENER,c1)*irl
    pl  = dv(DV_MIXT_PRES,c1)
    
    ql  = ul*nx + vl*ny + wl*nz
    
    rr  = cv(CV_MIXT_DENS,c2)
    irr = 1.0_RFREAL/rr
    
    rm  = 0.5_RFREAL*(rl + rr)
    
    ur  = cv(CV_MIXT_XMOM,c2)*irr
    vr  = cv(CV_MIXT_YMOM,c2)*irr
    wr  = cv(CV_MIXT_ZMOM,c2)*irr
    er  = cv(CV_MIXT_ENER,c2)*irr
    pr  = dv(DV_MIXT_PRES,c2)
    
    qr  = ur*nx + vr*ny + wr*nz    
   
! - Fluxes: mass and energy
      
    fc(1) = 0.5_RFREAL*(ql* rl          + qr* rr         )*nm
    fc(5) = 0.5_RFREAL*(ql*(rl*el + pl) + qr*(rr*er + pr))*nm    
    
! - Fluxes: momenta    

    fl(:) = 0.0_RFREAL
    fq(:) = 0.0_RFREAL

    DO i = 1,3  
      DO ic1l = 1,nCells           
        ic1g = fs(ic1l,ifc)

        v1(1:3) = cv(CV_MIXT_XMOM:CV_MIXT_ZMOM,ic1g)/cv(CV_MIXT_DENS,ic1g)

        DO j = 1,3
          fl(i) = fl(i) + wtl(i,j,ic1l,ifcp)*v1(j)          
        END DO ! j
      
        DO ic2l = 1,nCells
          ic2g = fs(ic2l,ifc)
          v2(1:3) = cv(CV_MIXT_XMOM:CV_MIXT_ZMOM,ic2g)/cv(CV_MIXT_DENS,ic2g)
        
          DO j = 1,3
            DO k = 1,3              
              fq(i) = fq(i) + wtq(i,j,k,ic1l,ic2l,ifcp)*v1(j)*v2(k)
            END DO ! k
          END DO ! j
         
        END DO ! ic2l
      END DO ! ic1l
         
      fl(i) = rm*fl(i)*nm
      fq(i) = rm*fq(i)*nm*SIGN(1.0_RFREAL,fn(ifcp,ifc))     
    
      fc(i+1) = fl(i) + fq(i)
    END DO ! i     

! - Add pressure fluxes (simple central difference)

    fc(2) = fc(2) + 0.5_RFREAL*(pl + pr)*nx*nm
    fc(3) = fc(3) + 0.5_RFREAL*(pl + pr)*ny*nm    
    fc(4) = fc(4) + 0.5_RFREAL*(pl + pr)*nz*nm

! - Add/subtract from residual
    
    rhs(CV_MIXT_DENS,c1) = rhs(CV_MIXT_DENS,c1) + fc(1)
    rhs(CV_MIXT_XMOM,c1) = rhs(CV_MIXT_XMOM,c1) + fc(2)
    rhs(CV_MIXT_YMOM,c1) = rhs(CV_MIXT_YMOM,c1) + fc(3)
    rhs(CV_MIXT_ZMOM,c1) = rhs(CV_MIXT_ZMOM,c1) + fc(4)
    rhs(CV_MIXT_ENER,c1) = rhs(CV_MIXT_ENER,c1) + fc(5)

    rhs(CV_MIXT_DENS,c2) = rhs(CV_MIXT_DENS,c2) - fc(1)
    rhs(CV_MIXT_XMOM,c2) = rhs(CV_MIXT_XMOM,c2) - fc(2)
    rhs(CV_MIXT_YMOM,c2) = rhs(CV_MIXT_YMOM,c2) - fc(3)
    rhs(CV_MIXT_ZMOM,c2) = rhs(CV_MIXT_ZMOM,c2) - fc(4)
    rhs(CV_MIXT_ENER,c2) = rhs(CV_MIXT_ENER,c2) - fc(5) 
    
! - Periodic boundary fluxes

    IF ( MAX(c1,c2) > region%grid%nCells ) THEN
      netMassFluxCntr = netMassFluxCntr + 1 
      netMassFlux = netMassFlux + fc(1)
    END IF ! MAX        
  END DO ! ifc

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ConvFluxOLES

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ConvFluxOLES.F90,v $
! Revision 1.7  2008/12/06 08:44:29  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2003/06/04 22:52:30  haselbac
! Removed unnecessary interface
!
! Revision 1.4  2003/04/10 14:38:36  haselbac
! Changed interface statement
!
! Revision 1.3  2003/01/28 14:33:45  haselbac
! Use parameters in fn
!
! Revision 1.2  2002/09/09 15:46:05  haselbac
! complete basic coding, global now under regions, several bug fixes
!
! Revision 1.1  2002/07/25 14:18:28  haselbac
! Initial revision
!
!******************************************************************************







