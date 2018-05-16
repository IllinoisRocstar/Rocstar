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
! Purpose: Calculate max. allowable local/global time step in the case
!          of viscous flow.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_TimeStepViscous.F90,v 1.16 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_TimeStepViscous(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModError
  USE ModParameters

  USE RFLU_ModGridSpeedUtils, ONLY: RFLU_DescaleGridSpeed
  USE RFLU_ModRindStates, ONLY: RFLU_SetRindStateInjectPerf  
  
  USE ModInterfaces, ONLY: MixtPerf_C_GRT, &
                           MixtPerf_G_CpR, & 
                           MixtPerf_R_M
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: c1,c2,distrib,gasModel,icg,ifg,ifl,indCp,indGs,indMol,iPatch, & 
             turbModel
  INTEGER :: dtMinLoc(1)
  INTEGER, DIMENSION(:), POINTER :: bf2c   
  INTEGER, DIMENSION(:,:), POINTER :: f2c 
  REAL(RFREAL) :: as,a1,a2,cp1,cp2,cp,dn,dtMin,dummyRealOut,dx,dy,dz,fmue,fs, & 
                  fsu,gam,Hl,ir1,ir2,minj,mm,mm1,mm2,nm,nx,ny,nz,pl,pr,prLam, &
                  prTurb,rgas,rl,rr,r1,r2,srad,tinj,uinj,ul,u1,u2,vinj,vl,v1, & 
                  v2,winj,wl,w1,w2,vc                 
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pTv,pVals
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_TimeStepViscous.F90,v $ $Revision: 1.16 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TimeStepViscous',&
  'RFLU_TimeStepViscous.F90')

! ******************************************************************************
! Get dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indGs = pGrid%indGs

  pCv => pRegion%mixt%cv
  pDv => pRegion%mixt%dv
  pGv => pRegion%mixt%gv
  pTv => pRegion%mixt%tv  
  
  prLam  = pRegion%mixtInput%prLam
  prTurb = pRegion%mixtInput%prTurb

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol
  
  turbModel = pRegion%mixtInput%turbModel
  gasModel  = pRegion%mixtInput%gasModel

! ******************************************************************************
! Initialize time step
! ******************************************************************************

  DO icg = 1,pGrid%nCellsTot ! Explicit loop to avoid Frost problem
    pRegion%dt(icg) = 0.0_RFREAL 
  END DO ! icg

! ******************************************************************************
! Local time step
! ******************************************************************************

! ==============================================================================
! Interior faces
! ==============================================================================

  DO ifg = 1,pGrid%nFaces
    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg)
    nm = pGrid%fn(XYZMAG,ifg)                 

    dx = pGrid%cofg(XCOORD,c2) - pGrid%cofg(XCOORD,c1)
    dy = pGrid%cofg(YCOORD,c2) - pGrid%cofg(YCOORD,c1)
    dz = pGrid%cofg(ZCOORD,c2) - pGrid%cofg(ZCOORD,c1)            
    dn = SQRT(dx*dx + dy*dy + dz*dz)

    fs = pGrid%gs(indGs*ifg)

    r1 = pCv(CV_MIXT_DENS,c1)
    r2 = pCv(CV_MIXT_DENS,c2)        

    ir1 = 1.0_RFREAL/r1
    ir2 = 1.0_RFREAL/r2

    u1 = pCv(CV_MIXT_XMOM,c1)*ir1
    v1 = pCv(CV_MIXT_YMOM,c1)*ir1
    w1 = pCv(CV_MIXT_ZMOM,c1)*ir1
    a1 = pDv(DV_MIXT_SOUN,c1)

    u2 = pCv(CV_MIXT_XMOM,c2)*ir2
    v2 = pCv(CV_MIXT_YMOM,c2)*ir2
    w2 = pCv(CV_MIXT_ZMOM,c2)*ir2
    a2 = pDv(DV_MIXT_SOUN,c2)

    vc = 0.5_RFREAL*((u1+u2)*nx + (v1+v2)*ny + (w1+w2)*nz)
    as = 0.5_RFREAL*(a1+a2)

#ifdef TURB
    IF ( turbModel == ACTIVE ) THEN
      fmue = 0.5_RFREAL*(pTv(TV_MIXT_MUEL,c1) + pTv(TV_MIXT_MUEL,c2))/prLam  & 
           + 0.5*RFREAL*(pTv(TV_MIXT_MUET,c1) + pTv(TV_MIXT_MUET,c2))/prTurb
    ELSE 
      fmue = 0.5_RFREAL*(pTv(TV_MIXT_MUEL,c1) + pTv(TV_MIXT_MUEL,c2))/prLam
    END IF ! turbModel
#else
    fmue = 0.5_RFREAL*(pTv(TV_MIXT_MUEL,c1) + pTv(TV_MIXT_MUEL,c2))/prLam
#endif      

    srad = (ABS(vc - fs) + as)*nm & 
         + 4.0_RFREAL*fmue/(0.5_RFREAL*(r1 + r2))*nm/dn

    pRegion%dt(c1) = pRegion%dt(c1) + srad
    pRegion%dt(c2) = pRegion%dt(c2) + srad
  END DO ! ifg

! ==============================================================================
! Boundary faces
! ==============================================================================

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)

    SELECT CASE ( pPatch%bcType ) 
    
! ------------------------------------------------------------------------------
!     Injection boundary
! ------------------------------------------------------------------------------   
    
      CASE ( BC_INJECTION:BC_INJECTION+BC_RANGE ) 
        pVals => pPatch%mixt%vals     
         
        distrib = pPatch%mixt%distrib
        
        DO ifl = 1,pPatch%nBFaces
          c1 = pPatch%bf2c(ifl)

          nx = pPatch%fn(XCOORD,ifl)
          ny = pPatch%fn(YCOORD,ifl)
          nz = pPatch%fn(ZCOORD,ifl)
          nm = pPatch%fn(XYZMAG,ifl)

          dx = pGrid%cofg(XCOORD,c1) - pPatch%fc(XCOORD,ifl)
          dy = pGrid%cofg(YCOORD,c1) - pPatch%fc(YCOORD,ifl)
          dz = pGrid%cofg(ZCOORD,c1) - pPatch%fc(ZCOORD,ifl)
          dn = SQRT(dx*dx + dy*dy + dz*dz)

          fs  = pPatch%gs(indGs*ifl)
          fsu = RFLU_DescaleGridSpeed(pRegion,fs)

          rl = pCv(CV_MIXT_DENS,c1)
          pl = pDv(DV_MIXT_PRES,c1)

          minj = pVals(BCDAT_INJECT_MFRATE,distrib*ifl)

          IF ( minj > 0.0_RFREAL ) THEN ! Surface burning
            tinj = pVals(BCDAT_INJECT_TEMP,distrib*ifl)

            IF ( gasModel == GAS_MODEL_TCPERF ) THEN
              cp = pGv(GV_MIXT_CP ,indCp *c1)        
              mm = pGv(GV_MIXT_MOL,indMol*c1)

              CALL RFLU_SetRindStateInjectPerf(cp,mm,nx,ny,nz,minj,tinj,pl, &
                                               fsu,rl,ul,vl,wl,Hl)
  
              r1   = rl
              vc   = ul*nx + vl*ny + wl*nz  
              rgas = MixtPerf_R_M(mm)
              gam  = MixtPerf_G_CpR(cp,rgas)              
              a1   = MixtPerf_C_GRT(gam,rgas,tinj)
            ELSE ! Defensive programming
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)              
            END IF ! gasModel
          ELSE ! Surface NOT burning
            r1 = rl
            vc = 0.0_RFREAL
            a1 = pDv(DV_MIXT_SOUN,c1)                        
          END IF ! minj
                  
#ifdef TURB
          IF ( turbModel == ACTIVE ) THEN
            fmue = pTv(TV_MIXT_MUEL,c1)/prLam + pTv(TV_MIXT_MUET,c1)/prTurb
          ELSE 
            fmue = pTv(TV_MIXT_MUEL,c1)/prLam
          END IF ! turbModel
#else
          fmue = pTv(TV_MIXT_MUEL,c1)/prLam
#endif  
                  
          pRegion%dt(c1) = pRegion%dt(c1) & 
                         + (ABS(vc - fsu) + a1)*nm + 4.0_RFREAL*fmue/r1*nm/dn
        END DO ! ifl

! ------------------------------------------------------------------------------
!     All other boundaries
! ------------------------------------------------------------------------------   
 
      CASE DEFAULT        
        DO ifl = 1,pPatch%nBFaces 
          c1 = pPatch%bf2c(ifl)

          nx = pPatch%fn(XCOORD,ifl)
          ny = pPatch%fn(YCOORD,ifl)
          nz = pPatch%fn(ZCOORD,ifl)
          nm = pPatch%fn(XYZMAG,ifl)

          dx = pGrid%cofg(XCOORD,c1) - pPatch%fc(XCOORD,ifl)
          dy = pGrid%cofg(YCOORD,c1) - pPatch%fc(YCOORD,ifl)
          dz = pGrid%cofg(ZCOORD,c1) - pPatch%fc(ZCOORD,ifl)
          dn = SQRT(dx*dx + dy*dy + dz*dz)

          fs  = pPatch%gs(indGs*ifl)
          fsu = RFLU_DescaleGridSpeed(pRegion,fs)           

          r1 = pCv(CV_MIXT_DENS,c1)        

          ir1 = 1.0_RFREAL/r1

          u1 = pCv(CV_MIXT_XMOM,c1)*ir1
          v1 = pCv(CV_MIXT_YMOM,c1)*ir1
          w1 = pCv(CV_MIXT_ZMOM,c1)*ir1
          a1 = pDv(DV_MIXT_SOUN,c1)

          vc = u1*nx + v1*ny + w1*nz

#ifdef TURB
          IF ( turbModel == ACTIVE ) THEN
            fmue = pTv(TV_MIXT_MUEL,c1)/prLam + pTv(TV_MIXT_MUET,c1)/prTurb
          ELSE 
            fmue = pTv(TV_MIXT_MUEL,c1)/prLam
          END IF ! turbModel
#else
          fmue = pTv(TV_MIXT_MUEL,c1)/prLam
#endif  

          pRegion%dt(c1) = pRegion%dt(c1) & 
                         + (ABS(vc - fsu) + a1)*nm + 4.0_RFREAL*fmue/r1*nm/dn
        END DO ! ifl                         
    END SELECT ! pPatch%bcType 
  END DO ! iPatch

! ******************************************************************************
! Finalize computation of local time step. NOTE set time step in dummy cells 
! to crazy value, that way code will blow up immediately should it be used by
! accident. 
! ******************************************************************************

  DO icg = 1,pGrid%nCells  
    pRegion%dt(icg) = pGrid%vol(icg)/pRegion%dt(icg)
  END DO ! icg

  DO icg = pGrid%nCells+1,pGrid%nCellsTot
    pRegion%dt(icg) = REAL(ABS(CRAZY_VALUE_INT),KIND=RFREAL)
  END DO ! icg

! ******************************************************************************
! For unsteady flows, determine minimum step
! ******************************************************************************

  IF ( global%flowType == FLOW_UNSTEADY ) THEN 
    dtMin    = MINVAL(pRegion%dt(1:pGrid%nCells))    
    dtMinLoc = MINLOC(pRegion%dt(1:pGrid%nCells))
    
    pRegion%dtMin    = dtMin
    pRegion%dtMinLoc = dtMinLoc(1)    
  END IF ! global%flowType

! ******************************************************************************
! Finalize
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_TimeStepViscous

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_TimeStepViscous.F90,v $
! Revision 1.16  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.15  2008/11/19 22:17:43  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.14  2006/08/19 15:40:02  mparmar
! Renamed patch variables
!
! Revision 1.13  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.12  2005/10/31 21:09:37  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.11  2005/05/28 18:03:56  haselbac
! Removed unnecessary setting of rgas and gamma
!
! Revision 1.10  2004/10/19 19:29:34  haselbac
! Now receive pRegion, partial rewrite, cosmetics
!
! Revision 1.9  2004/01/11 02:17:34  jiao
! Eliminated some redundant trailing spaces that made some lines too long.
! This changed was necessary to compile with NAG F90 compiler.
!
! Revision 1.8  2003/12/04 03:30:08  haselbac
! Changed timestep computation for viscous flows
!
! Revision 1.7  2003/10/15 02:43:57  haselbac
! Added code to determine location of minimum dt
!
! Revision 1.6  2003/08/28 00:58:43  haselbac
! Fixed bug in TURB section for timestep
!
! Revision 1.5  2003/04/18 20:00:09  haselbac
! Clean-up, added explicit initialization (prevent Frost-problem)
!
! Revision 1.4  2003/01/28 14:50:33  haselbac
! Use parameters in fn
!
! Revision 1.3  2002/11/27 20:27:45  haselbac
! Proper contribution from injecting boundaries
!
! Revision 1.2  2002/11/08 21:35:36  haselbac
! Fixed bug in use of grid speed
!
! Revision 1.1  2002/09/09 16:28:02  haselbac
! Initial revision
!
! ******************************************************************************







