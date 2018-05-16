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
! ******************************************************************************
!
! Purpose: Compute central convective fluxes for mixture using second-order 
!   accurate approximation through a patch by using an average of variables.
!
! Description: None.
!
! Input:
!   pRegion      Pointer to region
!   pPatch       Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_CentralSecondPatch_GL.F90,v 1.5 2008/12/06 08:44:29 mtcampbe Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CentralSecondPatch_GL(pRegion,pPatch)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
  USE ModMPI

  USE RFLU_ModGridSpeedUtils, ONLY: RFLU_DescaleGridSpeed
  USE RFLU_ModRindStates, ONLY: RFLU_SetRindStateFarfieldPerf, &
                                RFLU_SetRindStateInjectPerf, & 
                                RFLU_SetRindStateSlipWallPerf

  USE ModInterfaces, ONLY: BcondInflowPerf_GL, &
                           BcondOutflowPerf_GL, &
                           MixtGasLiq_C, &
                           MixtGasLiq_P, &
                           MixtLiq_C2_Bp, &
                           MixtPerf_C2_GRT, &
                           MixtPerf_Cv_CpR, &
                           MixtLiq_D_DoBpPPoBtTTo, &
                           MixtPerf_D_PRT, &                           
                           MixtPerf_G_CpR, &
                           MixtPerf_Eo_DGPUVW, &
                           MixtPerf_P_DEoGVm2, &
                           MixtPerf_P_DRT, &                           
                           MixtPerf_R_M, &
                           MixtPerf_R_CpG, &
                           RFLU_DecidePrint, &
                           RFLU_PrintLocInfo

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: corrFlag,decidePrintFlag
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER, PARAMETER :: MAX_INOUTFLOW_LOCS = 10
  INTEGER :: c1,bcOptFixed,bcOptType,distrib,ifl,indCp,indGs,indMf,indSd, &
             indMol,nLocs
  INTEGER :: loc(MAX_INOUTFLOW_LOCS,MIN_VAL:MAX_VAL)
  REAL(RFREAL) :: Bgl2,Bll2,Bp,Bt,Bvl2,Cgl2,Cll2,cml,cp,cvg,cvl,Cvl2, &
                  cvm,cvv,dx,dy,dz,fs,fsu,iCpRef,irl,mf,mm,nm,nx,ny,nz, &
                  pl,po,pr,pRef,press, ql,qr,rl,rel,rer,Rg,rgl,rgpgl, &
                  rgpgr,rhogr,rholr,rhovr,rll,rlpll,ro,rr,rRef,rul,rur, &
                  Rv,rvl,rvpvl,rvpvr,rvr,rwl,rwr, rYgl,rYll,rYvl,temp,tl, & 
                  to,tr,ul,ur,vfgr,vflr,vfvr,vl,Vl2,vm2,vr,vRef,wl,wr,xc, &
                  yc,zc,Ygl,Ygr,Yll,Ylr,Yvl,Yvr,vfgl,vfll,vfvl
  REAL(RFREAL) :: flx(5)
  REAL(RFREAL), DIMENSION(:), POINTER :: mfMixt
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,dv,gv,pCvSpec,rhs,vals,sd
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad,pGradSpec 
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_CentralSecondPatch_GL.F90,v $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CentralSecondPatch_GL',&
  'RFLU_CentralSecondPatch_GL.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  indCp  = pRegion%mixtInput%indCp
  indMf  = pRegion%mixtInput%indMfMixt
  indSd  = pRegion%mixtInput%indSd
  indMol = pRegion%mixtInput%indMol

  indGs = pRegion%grid%indGs

  distrib = pPatch%mixt%distrib

  cv  =>  pRegion%mixt%cv
  dv  =>  pRegion%mixt%dv
  gv  =>  pRegion%mixt%gv
  grad => pRegion%mixt%gradCell
  rhs =>  pRegion%mixt%rhs
  sd  =>  pRegion%mixt%sd

  pCvSpec => pRegion%spec%cv
  pGradSpec => pRegion%spec%gradCell

  mfMixt => pPatch%mfMixt

  nLocs = 0

  decidePrintFlag = RFLU_DecidePrint(global)

  pRef = global%refPressure
  rRef = global%refDensity
  vRef = global%refVelocity

  iCpRef = 2.0_RFREAL/(rRef*vRef*vRef)

! ******************************************************************************
! Define constants
! ******************************************************************************

  ro  = global%refDensityLiq
  po  = global%refPressLiq 
  to  = global%refTempLiq
  Bp  = global%refBetaPLiq
  Bt  = global%refBetaTLiq
  cvl = global%refCvLiq

  Rg  = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
  cvg = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht,Rg)

  Rv  = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
  cvv = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht,Rv)

! ******************************************************************************
! Select boundary type
! ******************************************************************************

  SELECT CASE ( pPatch%bcType )

! ==============================================================================
!   Inflow based on velocity and temperature
! ==============================================================================

    CASE ( BC_INFLOW_VELTEMP ) 
      vals       => pPatch%mixt%vals
      bcOptType  =  pPatch%mixt%switches(BCSWI_INFLOW_TYPE)

      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)
        nm = pPatch%fn(XYZMAG,ifl)

        xc = pPatch%fc(XCOORD,ifl)
        yc = pPatch%fc(YCOORD,ifl)
        zc = pPatch%fc(ZCOORD,ifl)

        fs = pPatch%gs(indGs*ifl)

        rl  = cv(CV_MIXT_DENS,c1)
        irl = 1.0_RFREAL/rl

        ul  = cv(CV_MIXT_XMOM,c1)*irl
        vl  = cv(CV_MIXT_YMOM,c1)*irl
        wl  = cv(CV_MIXT_ZMOM,c1)*irl
        tl  = dv(DV_MIXT_TEMP,c1)

        Ygl = pCvSpec(1,c1)*irl
        Yvl = pCvSpec(2,c1)*irl

        dx  = xc - pRegion%grid%cofg(XCOORD,c1)
        dy  = yc - pRegion%grid%cofg(YCOORD,c1)
        dz  = zc - pRegion%grid%cofg(ZCOORD,c1)

        rl = rl + grad(XCOORD,GRC_MIXT_DENS,c1)*dx &
                + grad(YCOORD,GRC_MIXT_DENS,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_DENS,c1)*dz
        ul = ul + grad(XCOORD,GRC_MIXT_XVEL,c1)*dx &
                + grad(YCOORD,GRC_MIXT_XVEL,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_XVEL,c1)*dz
        vl = vl + grad(XCOORD,GRC_MIXT_YVEL,c1)*dx &
                + grad(YCOORD,GRC_MIXT_YVEL,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_YVEL,c1)*dz
        wl = wl + grad(XCOORD,GRC_MIXT_ZVEL,c1)*dx &
                + grad(YCOORD,GRC_MIXT_ZVEL,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_ZVEL,c1)*dz
        tl = tl + grad(XCOORD,GRC_MIXT_TEMP,c1)*dx &
                + grad(YCOORD,GRC_MIXT_TEMP,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_TEMP,c1)*dz

        Ygl = Ygl + pGradSpec(XCOORD,1,c1)*dx &
                  + pGradSpec(YCOORD,1,c1)*dy &
                  + pGradSpec(ZCOORD,1,c1)*dz
        Yvl = Yvl + pGradSpec(XCOORD,2,c1)*dx &
                  + pGradSpec(YCOORD,2,c1)*dy &
                  + pGradSpec(ZCOORD,2,c1)*dz
        
        IF ( Ygl > 1.0_RFREAL) THEN
          Ygl = 1.0_RFREAL
          Yvl = 0.0_RFREAL
        ELSE IF ( Ygl < 0.0_RFREAL ) THEN
          Ygl = 0.0_RFREAL
        ELSE IF ( Yvl > 1.0_RFREAL ) THEN
          Yvl = 1.0_RFREAL
          Ygl = 0.0_RFREAL
        ELSE IF ( Yvl < 0.0_RFREAL ) THEN
          Yvl = 0.0_RFREAL
        END IF ! Ygl

        Yll = 1.0_RFREAL - Ygl - Yvl

        IF ( Yll > 1.0_RFREAL) THEN
          Yll = 1.0_RFREAL
        ELSE IF ( Yll < 0.0_RFREAL) THEN
          Yll = 0.0_RFREAL
        END IF ! Yll
 
        Vl2 = ul*ul + vl*vl + wl*wl

        Cll2 = MixtLiq_C2_Bp(Bp)
        Cvl2 = MixtPerf_C2_GRT(1.0_RFREAL,Rv,tl)
        Cgl2 = MixtPerf_C2_GRT(1.0_RFREAL,Rg,tl)

        rYgl = rl*Ygl
        rYvl = rl*Yvl
        rYll = rl*Yll 

        pl = MixtGasLiq_P(rYll,rYvl,rYgl,Cll2,Cvl2,Cgl2,rl, &
                          ro,po,to,Bp,Bt,tl)

        rll = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,pl,po,tl,to)
        rvl = MixtPerf_D_PRT(pl,Rv,tl)
        rgl = MixtPerf_D_PRT(pl,Rg,tl)

        vfgl = rYgl/rgl
        vfvl = rYvl/rvl
        vfll = rYll/rll    

        cvm = (rYll*cvl + rYvl*cvv + rYgl*cvg)/rl

        Bll2 = -Bt/Bp
        Bvl2 = rvl*Rv
        Bgl2 = rgl*Rg

        cml = MixtGasLiq_C(cvm,rl,pl,rll,rvl,rgl,vfll,vfvl,vfgl, &
                           Cll2,Cvl2,Cgl2,Bll2,Bvl2,Bgl2)
        rul = rl*ul
        rvl = rl*vl
        rwl = rl*wl
        rel = rl*cvm*tl + 0.5*rl*Vl2

        rgpgl = rYgl 
        rvpvl = rYvl
        rlpll = rYll  

        ur    = vals(BCDAT_INFLOW_U,distrib*ifl)
        vr    = vals(BCDAT_INFLOW_V,distrib*ifl)
        wr    = vals(BCDAT_INFLOW_W,distrib*ifl)
        temp  = vals(BCDAT_INFLOW_T,distrib*ifl)

        Ygr   = pPatch%spec%vals(1,distrib*ifl)  
        Yvr   = pPatch%spec%vals(2,distrib*ifl)
        Ylr   = 1.0_RFREAL - Ygr - Yvr

        IF ( bcOptType /= BCOPT_SUBSONIC ) THEN
          press  = vals(BCDAT_INFLOW_P, distrib*ifl)
          rholr  = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,press,po,temp,to)
          rhogr  = MixtPerf_D_PRT(press,Rg,temp)
          rhovr  = MixtPerf_D_PRT(press,Rv,temp)
          rr     = 1.0_RFREAL/(Ygr/rhogr + Yvr/rhovr + Ylr/rholr) 
          vfgr   = (rr*Ygr)/rhogr
          vfvr   = (rr*Yvr)/rhovr
          vflr   = (rr*Ylr)/rholr 
        ELSE
          press  = pl+(SQRT(ur*ur+vr*vr+wr*wr) &
                 - SQRT(ul*ul+vl*vl+wl*wl))*rl*cml  
          rholr  = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,press,po,temp,to)
          rhogr  = MixtPerf_D_PRT(press,Rg,temp)
          rhovr  = MixtPerf_D_PRT(press,Rv,temp)
          rr     = 1.0_RFREAL/(Ygr/rhogr + Yvr/rhovr + Ylr/rholr)
          vfgr   = (rr*Ygr)/rhogr
          vfvr   = (rr*Yvr)/rhovr
          vflr   = (rr*Ylr)/rholr
        END IF ! bcOptType

        CALL BcondInflowPerf_GL(bcOptType,ro,po,to,Bp,Bt,cvl,cvv,cvg,Rg,Rv, &
                                ur,vr,wr,vfgr,vfvr,vflr,temp,press,nx,ny, &  
                                nz,rl,rul,rvl,rwl,rel,rgpgl,rvpvl,pl,rr, &  
                                rur,rvr,rwr,rer,rgpgr,rvpvr,pr)

        ql = (rul*nx + rvl*ny + rwl*nz)/rl - fs
        qr = (rur*nx + rvr*ny + rwr*nz)/rr - fs

        flx(1) = 0.5_RFREAL*(ql* rl                + qr* rr               )*nm
        flx(2) = 0.5_RFREAL*(ql* rul       + pl*nx + qr* rur       + pr*nx)*nm
        flx(3) = 0.5_RFREAL*(ql* rvl       + pl*ny + qr* rvr       + pr*ny)*nm
        flx(4) = 0.5_RFREAL*(ql* rwl       + pl*nz + qr* rwr       + pr*nz)*nm
        flx(5) = 0.5_RFREAL*(ql*(rel + pl) + pl*fs + qr*(rer + pr) + pr*fs)*nm

        mfMixt(indMf*ifl) = flx(1)
        pPatch%cp(ifl) = iCpRef*(0.5_RFREAL*(pl + pr) - pRef)

        rhs(CV_MIXT_DENS,c1) = rhs(CV_MIXT_DENS,c1) + flx(1)
        rhs(CV_MIXT_XMOM,c1) = rhs(CV_MIXT_XMOM,c1) + flx(2)
        rhs(CV_MIXT_YMOM,c1) = rhs(CV_MIXT_YMOM,c1) + flx(3)
        rhs(CV_MIXT_ZMOM,c1) = rhs(CV_MIXT_ZMOM,c1) + flx(4)
        rhs(CV_MIXT_ENER,c1) = rhs(CV_MIXT_ENER,c1) + flx(5)

        sd(SD_XMOM,c1*indSd) = sd(SD_XMOM,c1*indSd) &
                             + 0.5_RFREAL*(rul/rl + rur/rr)*flx(1)
        sd(SD_YMOM,c1*indSd) = sd(SD_YMOM,c1*indSd) &
                             + 0.5_RFREAL*(rvl/rl + rvr/rr)*flx(1)
        sd(SD_ZMOM,c1*indSd) = sd(SD_ZMOM,c1*indSd) &
                             + 0.5_RFREAL*(rwl/rl + rwr/rr)*flx(1)

        IF ( pRegion%irkStep == 1 ) THEN
          global%massIn  = global%massIn  - MIN(flx(1),0.0_RFREAL)
          global%massOut = global%massOut + MAX(flx(1),0.0_RFREAL)
        END IF ! pRegion

        IF ( (global%checkLevel == CHECK_HIGH) .AND. &
             (global%verbLevel >= VERBOSE_HIGH) .AND. &
             (global%myProcid == MASTERPROC) .AND. &
             (decidePrintFlag .EQV. .TRUE.) ) THEN
          IF ( flx(1) > 0.0_RFREAL ) THEN
            nLocs = nLocs + 1

            IF ( nLocs == 1 ) THEN
              global%warnCounter = global%warnCounter + 1

              WRITE(STDOUT,'(A,1X,A,1X,I9)') SOLVER_NAME, &
                    '*** WARNING *** Outflow detected at inflow boundary!'
              WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                               pRegion%iRegionGlobal
              IF ( global%flowType == FLOW_UNSTEADY ) THEN
                WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME, &
                                                    'Current time:', &
                                                    global%currentTime
              ELSE
                WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME, &
                                                 'Current iteration number:', &
                                                 global%currentIter
              END IF ! global%flowType

              WRITE(STDOUT,'(A,3X,A,1X,I1)') SOLVER_NAME, &
                                             'Runge-Kutta stage:', &
                                             pRegion%irkStep
            END IF ! nLocs

            IF ( nLocs <= MAX_INOUTFLOW_LOCS ) THEN
              loc(nLocs,MIN_VAL:MAX_VAL) = c1
            END IF ! nLocs
          END IF ! flx(1)
        END IF ! global%checkLevel
      END DO ! ifl

! ==============================================================================
!   Outflow
! ==============================================================================

    CASE ( BC_OUTFLOW )
      bcOptType = pPatch%mixt%switches(BCSWI_OUTFLOW_TYPE)

      IF ( bcOptType /= BCOPT_SUPERSONIC ) THEN
        vals => pPatch%mixt%vals
      END IF ! bcOptType

      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)
        nm = pPatch%fn(XYZMAG,ifl)

        xc = pPatch%fc(XCOORD,ifl)
        yc = pPatch%fc(YCOORD,ifl)
        zc = pPatch%fc(ZCOORD,ifl)

        fs = pPatch%gs(indGs*ifl)

        rl  = cv(CV_MIXT_DENS,c1)
        irl = 1.0_RFREAL/rl

        ul  = cv(CV_MIXT_XMOM,c1)*irl
        vl  = cv(CV_MIXT_YMOM,c1)*irl
        wl  = cv(CV_MIXT_ZMOM,c1)*irl
        tl  = dv(DV_MIXT_TEMP,c1)

        Ygl = pCvSpec(1,c1)*irl
        Yvl = pCvSpec(2,c1)*irl

        dx  = xc - pRegion%grid%cofg(XCOORD,c1)
        dy  = yc - pRegion%grid%cofg(YCOORD,c1)
        dz  = zc - pRegion%grid%cofg(ZCOORD,c1)

        rl = rl + grad(XCOORD,GRC_MIXT_DENS,c1)*dx &
                + grad(YCOORD,GRC_MIXT_DENS,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_DENS,c1)*dz
        ul = ul + grad(XCOORD,GRC_MIXT_XVEL,c1)*dx &
                + grad(YCOORD,GRC_MIXT_XVEL,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_XVEL,c1)*dz
        vl = vl + grad(XCOORD,GRC_MIXT_YVEL,c1)*dx &
                + grad(YCOORD,GRC_MIXT_YVEL,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_YVEL,c1)*dz
        wl = wl + grad(XCOORD,GRC_MIXT_ZVEL,c1)*dx &
                + grad(YCOORD,GRC_MIXT_ZVEL,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_ZVEL,c1)*dz
        tl = tl + grad(XCOORD,GRC_MIXT_TEMP,c1)*dx &
                + grad(YCOORD,GRC_MIXT_TEMP,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_TEMP,c1)*dz

        Ygl = Ygl + pGradSpec(XCOORD,1,c1)*dx &
                  + pGradSpec(YCOORD,1,c1)*dy &
                  + pGradSpec(ZCOORD,1,c1)*dz
        Yvl = Yvl + pGradSpec(XCOORD,2,c1)*dx &
                  + pGradSpec(YCOORD,2,c1)*dy &
                  + pGradSpec(ZCOORD,2,c1)*dz

        IF ( Ygl > 1.0_RFREAL ) THEN
          Ygl = 1.0_RFREAL
          Yvl = 0.0_RFREAL
        ELSE IF ( Ygl < 0.0_RFREAL ) THEN
          Ygl = 0.0_RFREAL
        ELSE IF ( Yvl > 1.0_RFREAL ) THEN
          Yvl = 1.0_RFREAL
          Ygl = 0.0_RFREAL
        ELSE IF ( Yvl < 0.0_RFREAL ) THEN
          Yvl = 0.0_RFREAL
        END IF ! Ygl
 
        Yll = 1.0_RFREAL - Ygl - Yvl

        IF ( Yll > 1.0_RFREAL) THEN
          Yll = 1.0_RFREAL
        ELSE IF ( Yll < 0.0_RFREAL) THEN
          Yll = 0.0_RFREAL
        END IF ! Yll

        Vl2 = ul*ul + vl*vl + wl*wl

        Cll2 = MixtLiq_C2_Bp(Bp)
        Cvl2 = MixtPerf_C2_GRT(1.0_RFREAL,Rv,tl)
        Cgl2 = MixtPerf_C2_GRT(1.0_RFREAL,Rg,tl)

        rYgl = rl*Ygl
        rYvl = rl*Yvl
        rYll = rl*Yll

        pl = MixtGasLiq_P(rYll,rYvl,rYgl,Cll2,Cvl2,Cgl2,rl,ro,po,to,Bp,Bt,tl)

        cvm = (rYll*cvl + rYvl*cvv + rYgl*cvg)/rl

        rul = rl*ul
        rvl = rl*vl
        rwl = rl*wl
        rel = rl*cvm*tl + 0.5*rl*Vl2       
        rgpgl = rYgl
        rvpvl = rYvl
        rlpll = rYll


        IF ( bcOptType /= BCOPT_SUPERSONIC ) THEN
          pr = vals(BCDAT_OUTFLOW_PRESS,distrib*ifl)
        ELSE
          pr = pl
        END IF ! bcOptType

        CALL BcondOutflowPerf_GL(bcOptType,ro,po,to,Bp,Bt,cvl,cvv, &
                                 cvg,Rg,Rv,pr,nx,ny,nz,rl,rul,rvl, &
                                 rwl,rel,rgpgl,rvpvl,pl,rr,rur,rvr, &
                                 rwr,rer,rgpgr,rvpvr)

        ql = (rul*nx + rvl*ny + rwl*nz)/rl - fs
        qr = (rur*nx + rvr*ny + rwr*nz)/rr - fs

        flx(1) = 0.5_RFREAL*(ql* rl                + qr* rr               )*nm
        flx(2) = 0.5_RFREAL*(ql* rul       + pl*nx + qr* rur       + pr*nx)*nm
        flx(3) = 0.5_RFREAL*(ql* rvl       + pl*ny + qr* rvr       + pr*ny)*nm
        flx(4) = 0.5_RFREAL*(ql* rwl       + pl*nz + qr* rwr       + pr*nz)*nm 
        flx(5) = 0.5_RFREAL*(ql*(rel + pl) + pl*fs + qr*(rer + pr) + pr*fs)*nm

        mfMixt(indMf*ifl) = flx(1)
        pPatch%cp(ifl) = iCpRef*(0.5_RFREAL*(pl + pr) - pRef)

        rhs(CV_MIXT_DENS,c1) = rhs(CV_MIXT_DENS,c1) + flx(1)
        rhs(CV_MIXT_XMOM,c1) = rhs(CV_MIXT_XMOM,c1) + flx(2)
        rhs(CV_MIXT_YMOM,c1) = rhs(CV_MIXT_YMOM,c1) + flx(3)
        rhs(CV_MIXT_ZMOM,c1) = rhs(CV_MIXT_ZMOM,c1) + flx(4)
        rhs(CV_MIXT_ENER,c1) = rhs(CV_MIXT_ENER,c1) + flx(5)

        sd(SD_XMOM,c1*indSd) = sd(SD_XMOM,c1*indSd) &
                             + 0.5_RFREAL*(rul/rl + rur/rr)*flx(1)
        sd(SD_YMOM,c1*indSd) = sd(SD_YMOM,c1*indSd) &
                             + 0.5_RFREAL*(rvl/rl + rvr/rr)*flx(1)
        sd(SD_ZMOM,c1*indSd) = sd(SD_ZMOM,c1*indSd) &
                             + 0.5_RFREAL*(rwl/rl + rwr/rr)*flx(1)

        IF ( pRegion%irkStep == 1 ) THEN
          global%massIn  = global%massIn  - MIN(flx(1),0.0_RFREAL)
          global%massOut = global%massOut + MAX(flx(1),0.0_RFREAL)
        END IF ! pRegion

        IF ( (global%checkLevel == CHECK_HIGH) .AND. &
             (global%verbLevel >= VERBOSE_HIGH) .AND. &
             (global%myProcid == MASTERPROC) .AND. &
             (decidePrintFlag .EQV. .TRUE.) ) THEN
          IF ( flx(1) < 0.0_RFREAL ) THEN
            nLocs = nLocs + 1

            IF ( nLocs == 1 ) THEN
              global%warnCounter = global%warnCounter + 1

              WRITE(STDOUT,'(A,1X,A,1X,I9)') SOLVER_NAME, &
                    '*** WARNING *** Inflow detected at outflow boundary!'
              WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                               pRegion%iRegionGlobal
              IF ( global%flowType == FLOW_UNSTEADY ) THEN
                WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME, &
                                                    'Current time:', &
                                                    global%currentTime
              ELSE
                WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME, &
                                                 'Current iteration number:', &
                                                 global%currentIter
              END IF ! global%flowType

              WRITE(STDOUT,'(A,3X,A,1X,I1)') SOLVER_NAME, &
                                             'Runge-Kutta stage:', &
                                             pRegion%irkStep
            END IF ! nLocs

            IF ( nLocs <= MAX_INOUTFLOW_LOCS ) THEN
              loc(nLocs,MIN_VAL:MAX_VAL) = c1
            END IF ! nLocs
          END IF ! flx(1)
        END IF ! global%checkLevel
      END DO ! ifl

! ------------------------------------------------------------------------------
!     Write info on inflow at outflow boundary
! ------------------------------------------------------------------------------

      IF ( (global%checkLevel == CHECK_HIGH) .AND. &
           (global%verbLevel >= VERBOSE_HIGH) .AND. &
           (global%myProcid == MASTERPROC) .AND. &
           (decidePrintFlag .EQV. .TRUE.) .AND. &
           (nLocs > 0) ) THEN
        IF ( nLocs > MAX_INOUTFLOW_LOCS ) THEN
           WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A,1X,I9,1X,A)') SOLVER_NAME, &
                 'Only wrote the first',MAX_INOUTFLOW_LOCS,'of',nLocs, &
                 'outflow faces with inflow.'
          CALL RFLU_PrintLocInfo(pRegion,loc,MAX_INOUTFLOW_LOCS, &
                                 LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
        ELSE
          CALL RFLU_PrintLocInfo(pRegion,loc(1:nLocs,MIN_VAL:MAX_VAL),nLocs, &
                                 LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
        END IF ! nLocs
      END IF ! global%verbLevel

! ==============================================================================
!   Slip wall (weak imposition)
! ==============================================================================

    CASE ( BC_SLIPWALL )
      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)
        nm = pPatch%fn(XYZMAG,ifl)

        xc = pPatch%fc(XCOORD,ifl)
        yc = pPatch%fc(YCOORD,ifl)
        zc = pPatch%fc(ZCOORD,ifl)

        fs  = pPatch%gs(indGs*ifl)
        fsu = RFLU_DescaleGridSpeed(pRegion,fs)

        rl  = cv(CV_MIXT_DENS,c1)
        irl = 1.0_RFREAL/rl

        ul  = cv(CV_MIXT_XMOM,c1)*irl
        vl  = cv(CV_MIXT_YMOM,c1)*irl
        wl  = cv(CV_MIXT_ZMOM,c1)*irl
        tl  = dv(DV_MIXT_TEMP,c1)

        Ygl = pCvSpec(1,c1)*irl
        Yvl = pCvSpec(2,c1)*irl

        dx  = xc - pRegion%grid%cofg(XCOORD,c1)
        dy  = yc - pRegion%grid%cofg(YCOORD,c1)
        dz  = zc - pRegion%grid%cofg(ZCOORD,c1)

        rl = rl + grad(XCOORD,GRC_MIXT_DENS,c1)*dx &
                + grad(YCOORD,GRC_MIXT_DENS,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_DENS,c1)*dz
        ul = ul + grad(XCOORD,GRC_MIXT_XVEL,c1)*dx &
                + grad(YCOORD,GRC_MIXT_XVEL,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_XVEL,c1)*dz
        vl = vl + grad(XCOORD,GRC_MIXT_YVEL,c1)*dx &
                + grad(YCOORD,GRC_MIXT_YVEL,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_YVEL,c1)*dz
        wl = wl + grad(XCOORD,GRC_MIXT_ZVEL,c1)*dx &
                + grad(YCOORD,GRC_MIXT_ZVEL,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_ZVEL,c1)*dz
        tl = tl + grad(XCOORD,GRC_MIXT_TEMP,c1)*dx &
                + grad(YCOORD,GRC_MIXT_TEMP,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_TEMP,c1)*dz

        Ygl = Ygl + pGradSpec(XCOORD,1,c1)*dx &
                  + pGradSpec(YCOORD,1,c1)*dy &
                  + pGradSpec(ZCOORD,1,c1)*dz
        Yvl = Yvl + pGradSpec(XCOORD,2,c1)*dx &
                  + pGradSpec(YCOORD,2,c1)*dy &
                  + pGradSpec(ZCOORD,2,c1)*dz
  
        IF ( Ygl > 1.0_RFREAL ) THEN
          Ygl = 1.0_RFREAL
          Yvl = 0.0_RFREAL
        ELSE IF ( Ygl < 0.0_RFREAL ) THEN
          Ygl = 0.0_RFREAL
        ELSE IF ( Yvl > 1.0_RFREAL ) THEN
          Yvl = 1.0_RFREAL
          Ygl = 0.0_RFREAL
        ELSE IF ( Yvl < 0.0_RFREAL ) THEN
          Yvl = 0.0_RFREAL
        END IF ! Ygl
 
        Yll = 1.0_RFREAL - Ygl - Yvl

        IF ( Yll > 1.0_RFREAL) THEN
          Yll = 1.0_RFREAL
        ELSE IF ( Yll < 0.0_RFREAL) THEN
          Yll = 0.0_RFREAL
        END IF ! Yll

        Cll2 = MixtLiq_C2_Bp(Bp)
        Cvl2 = MixtPerf_C2_GRT(1.0_RFREAL,Rv,tl)
        Cgl2 = MixtPerf_C2_GRT(1.0_RFREAL,Rg,tl)

        rYgl = rl*Ygl
        rYvl = rl*Yvl
        rYll = rl*Yll 

        pl = MixtGasLiq_P(rYll,rYvl,rYgl,Cll2,Cvl2,Cgl2,rl,ro,po,to,Bp,Bt,tl)

        cvm = (rYll*cvl + rYvl*cvv + rYgl*cvg)/rl

        rul = rl*ul
        rvl = rl*vl
        rwl = rl*wl

        cp = gv(GV_MIXT_CP ,indCp *c1)
        mm = gv(GV_MIXT_MOL,indMol*c1)
   
        CALL RFLU_SetRindStateSlipWallPerf(cp,mm,nx,ny,nz,rl,rul, &
                                           rvl,rwl,fsu,pl)

        flx(2) = pl*nx*nm 
        flx(3) = pl*ny*nm
        flx(4) = pl*nz*nm
        flx(5) = pl*fs*nm

        mfMixt(indMf*ifl) = 0.0_RFREAL
        pPatch%cp(ifl) = iCpRef*(pl - pRef)

        rhs(CV_MIXT_XMOM,c1) = rhs(CV_MIXT_XMOM,c1) + flx(2)
        rhs(CV_MIXT_YMOM,c1) = rhs(CV_MIXT_YMOM,c1) + flx(3)
        rhs(CV_MIXT_ZMOM,c1) = rhs(CV_MIXT_ZMOM,c1) + flx(4)
        rhs(CV_MIXT_ENER,c1) = rhs(CV_MIXT_ENER,c1) + flx(5)
      END DO ! ifl

! ==============================================================================
!   No-slip wall
! ==============================================================================

    CASE ( BC_NOSLIPWALL )
      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)
        nm = pPatch%fn(XYZMAG,ifl)

        xc = pPatch%fc(XCOORD,ifl)
        yc = pPatch%fc(YCOORD,ifl)
        zc = pPatch%fc(ZCOORD,ifl)

        fs = pPatch%gs(indGs*ifl)

        rl  = cv(CV_MIXT_DENS,c1)
        tl  = dv(DV_MIXT_TEMP,c1)
        Ygl = pCvSpec(1,c1)*irl
        Yvl = pCvSpec(2,c1)*irl

        dx  = xc - pRegion%grid%cofg(XCOORD,c1)
        dy  = yc - pRegion%grid%cofg(YCOORD,c1)
        dz  = zc - pRegion%grid%cofg(ZCOORD,c1)

        rl = rl + grad(XCOORD,GRC_MIXT_DENS,c1)*dx &
                + grad(YCOORD,GRC_MIXT_DENS,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_DENS,c1)*dz
        tl = tl + grad(XCOORD,GRC_MIXT_TEMP,c1)*dx &
                + grad(YCOORD,GRC_MIXT_TEMP,c1)*dy &
                + grad(ZCOORD,GRC_MIXT_TEMP,c1)*dz

        Ygl = Ygl + pGradSpec(XCOORD,1,c1)*dx &
                  + pGradSpec(YCOORD,1,c1)*dy &
                  + pGradSpec(ZCOORD,1,c1)*dz
        Yvl = Yvl + pGradSpec(XCOORD,2,c1)*dx &
                  + pGradSpec(YCOORD,2,c1)*dy &
                  + pGradSpec(ZCOORD,2,c1)*dz
      
        IF ( Ygl > 1.0_RFREAL ) THEN
          Ygl = 1.0_RFREAL
          Yvl = 0.0_RFREAL
        ELSE IF ( Ygl < 0.0_RFREAL ) THEN
          Ygl = 0.0_RFREAL
        ELSE IF ( Yvl > 1.0_RFREAL ) THEN
          Yvl = 1.0_RFREAL
          Ygl = 0.0_RFREAL
        ELSE IF ( Yvl < 0.0_RFREAL ) THEN
          Yvl = 0.0_RFREAL
        END IF ! Ygl
 
        Yll = 1.0_RFREAL - Ygl - Yvl

        IF ( Yll > 1.0_RFREAL) THEN
          Yll = 1.0_RFREAL
        ELSE IF ( Yll < 0.0_RFREAL) THEN
          Yll = 0.0_RFREAL
        END IF ! Yll

        Cll2 = MixtLiq_C2_Bp(Bp)
        Cvl2 = MixtPerf_C2_GRT(1.0_RFREAL,Rv,tl)
        Cgl2 = MixtPerf_C2_GRT(1.0_RFREAL,Rg,tl)

        rYgl = rl*Ygl
        rYvl = rl*Yvl
        rYll = rl*Yll 

        pl = MixtGasLiq_P(rYll,rYvl,rYgl,Cll2,Cvl2,Cgl2,rl,ro,po,to,Bp,Bt,tl)

        flx(2) = pl*nx*nm 
        flx(3) = pl*ny*nm
        flx(4) = pl*nz*nm
        flx(5) = pl*fs*nm

        mfMixt(indMf*ifl) = 0.0_RFREAL
        pPatch%cp(ifl) = iCpRef*(pl - pRef)

        rhs(CV_MIXT_XMOM,c1) = rhs(CV_MIXT_XMOM,c1) + flx(2)
        rhs(CV_MIXT_YMOM,c1) = rhs(CV_MIXT_YMOM,c1) + flx(3)
        rhs(CV_MIXT_ZMOM,c1) = rhs(CV_MIXT_ZMOM,c1) + flx(4)
        rhs(CV_MIXT_ENER,c1) = rhs(CV_MIXT_ENER,c1) + flx(5)
      END DO ! ifl

! ==============================================================================
!   Boundaries for which fluxes must not or need not be computed
! ==============================================================================

    CASE ( BC_PERIODIC, & 
           BC_SYMMETRY, & 
           BC_VIRTUAL )

! ==============================================================================
!   Default
! ==============================================================================

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pPatch%bcType

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CentralSecondPatch_GL

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CentralSecondPatch_GL.F90,v $
! Revision 1.5  2008/12/06 08:44:29  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/08/19 15:39:56  mparmar
! Renamed patch variables
!
! Revision 1.2  2006/04/13 20:02:27  haselbac
! Added periodic and virtual boundary CASE statement
!
! Revision 1.1  2006/03/26 20:21:03  haselbac
! Initial revision
!
! ******************************************************************************







