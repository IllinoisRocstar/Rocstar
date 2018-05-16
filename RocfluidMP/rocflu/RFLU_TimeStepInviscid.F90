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
!          of inviscid flow.
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_TimeStepInviscid.F90,v 1.23 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_TimeStepInviscid(pRegion)

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

#ifdef PETSC
  USE RFLU_ModPETScNewtonKrylov, ONLY: RFLU_PETSC_GetDtScale
#endif
  
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
  INTEGER :: c1,c2,distrib,gasModel,icg,ifg,ifl,indCp,indGs,indMol,iPatch
  INTEGER :: dtMinLoc(1) 
  REAL(RFREAL) :: as,a1,a2,cp,dtMin,dtScale,dummyRealout,fs,fsu,gam,Hl,ir1, &
                  ir2,minj,mm,nm,nx,ny,nz,pl,pr,rgas,rl,rr,rur,rvr,rwr,srad, &
                  tinj,uinj,ul,u1,u2,vc,vinj,vl,vm2,v1,v2,winj,wl,w1,w2
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pVals
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_TimeStepInviscid.F90,v $ $Revision: 1.23 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_TimeStepInviscid',&
  'RFLU_TimeStepInviscid.F90')

! ******************************************************************************
! Get dimensions and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  indCp    = pRegion%mixtInput%indCp
  indMol   = pRegion%mixtInput%indMol
  gasModel = pRegion%mixtInput%gasModel

  pCv => pRegion%mixt%cv
  pDv => pRegion%mixt%dv
  pGv => pRegion%mixt%gv

  indGs = pGrid%indGs
  
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

  DO ifg = 1,pRegion%grid%nFaces

    c1 = pGrid%f2c(1,ifg)
    c2 = pGrid%f2c(2,ifg)

    nx = pGrid%fn(XCOORD,ifg)
    ny = pGrid%fn(YCOORD,ifg)
    nz = pGrid%fn(ZCOORD,ifg)
    nm = pGrid%fn(XYZMAG,ifg)                  

    fs  = pGrid%gs(indGs*ifg)
    fsu = RFLU_DescaleGridSpeed(pRegion,fs)    

    ir1 = 1.0_RFREAL/pCv(CV_MIXT_DENS,c1)
    ir2 = 1.0_RFREAL/pCv(CV_MIXT_DENS,c2)

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

    srad = (ABS(vc - fsu) + as)*nm

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

          fs  = pPatch%gs(indGs*ifl)
          fsu = RFLU_DescaleGridSpeed(pRegion,fs)

          pl = pDv(DV_MIXT_PRES,c1)

          minj = pVals(BCDAT_INJECT_MFRATE,distrib*ifl)

          IF ( minj > 0.0_RFREAL ) THEN ! Surface burning
            tinj = pVals(BCDAT_INJECT_TEMP,distrib*ifl)          
          
            IF ( gasModel == GAS_MODEL_TCPERF ) THEN
              cp = pGv(GV_MIXT_CP ,indCp *c1)        
              mm = pGv(GV_MIXT_MOL,indMol*c1)

              CALL RFLU_SetRindStateInjectPerf(cp,mm,nx,ny,nz,minj,tinj,pl, &
                                               fsu,rl,ul,vl,wl,Hl)
  
              vc   = ul*nx + vl*ny + wl*nz  
              rgas = MixtPerf_R_M(mm)
              gam  = MixtPerf_G_CpR(cp,rgas)              
              a1   = MixtPerf_C_GRT(gam,rgas,tinj)
            ELSE ! Defensive programming
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)              
            END IF ! gasModel
          ELSE ! Surface NOT burning
            vc = fsu
            a1 = pDv(DV_MIXT_SOUN,c1)
          END IF ! minj
                  
          pRegion%dt(c1) = pRegion%dt(c1) + (ABS(vc - fsu) + a1)*nm
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

          fs  = pPatch%gs(indGs*ifl)
          fsu = RFLU_DescaleGridSpeed(pRegion,fs)           

          ir1 = 1.0_RFREAL/pCv(CV_MIXT_DENS,c1)

          u1 = pCv(CV_MIXT_XMOM,c1)*ir1
          v1 = pCv(CV_MIXT_YMOM,c1)*ir1
          w1 = pCv(CV_MIXT_ZMOM,c1)*ir1
          a1 = pDv(DV_MIXT_SOUN,c1)

          vc = u1*nx + v1*ny + w1*nz

          pRegion%dt(c1) = pRegion%dt(c1) + (ABS(vc - fsu) + a1)*nm
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
! For Newton-Krylov solvers, ramp the timestep
! ******************************************************************************
 
  IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN
    dtScale = 1.0
#ifdef PETSC
    CALL RFLU_PETSC_GetDtScale(global,dtScale)
#endif
    DO icg = 1,pGrid%nCells
      pRegion%dt(icg) = dtScale*pRegion%dt(icg)
    ENDDO ! icg
  END IF ! global%solverType

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

END SUBROUTINE RFLU_TimeStepInviscid

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_TimeStepInviscid.F90,v $
! Revision 1.23  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.22  2008/11/19 22:17:43  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.21  2006/08/19 15:40:00  mparmar
! Renamed patch variables
!
! Revision 1.20  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.19  2006/02/08 21:34:46  hdewey2
! Set dtScale so not undefined when running implicit solver wo PETSC
!
! Revision 1.18  2005/11/14 17:00:58  haselbac
! Removed DEBUG statement
!
! Revision 1.17  2005/10/31 21:09:37  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.16  2005/09/23 14:56:01  haselbac
! Bug fix: Missing ifdef PETSC
!
! Revision 1.15  2005/09/22 17:12:33  hdewey2
! Changed implicit timestep ramping to call RFLU_PETSC_GetDtScale to get the 
! timestep scaling.
!
! Revision 1.14  2005/08/17 20:28:24  hdewey2
! Added timestep ramping for the implicit Newton-Krylov solver
!
! Revision 1.13  2004/10/19 19:29:26  haselbac
! Now receive pRegion, partial rewrite, cosmetics
!
! Revision 1.12  2004/01/11 02:17:34  jiao
! Eliminated some redundant trailing spaces that made some lines too long.
! This changed was necessary to compile with NAG F90 compiler.
!
! Revision 1.11  2003/10/15 02:43:57  haselbac
! Added code to determine location of minimum dt
!
! Revision 1.10  2003/09/04 20:19:37  haselbac
! Removed temporary fix for Rocburn problem in GENx
!
! Revision 1.9  2003/08/22 20:31:56  haselbac
! Added temporary fix for ignition problems
!
! Revision 1.8  2003/04/18 20:00:09  haselbac
! Clean-up, added explicit initialization (prevent Frost-problem)
!
! Revision 1.7  2003/03/15 18:58:15  haselbac
! Clean-up, cosmetics
!
! Revision 1.6  2003/01/28 14:50:06  haselbac
! Use parameters in fn
!
! Revision 1.5  2002/11/27 20:27:45  haselbac
! Proper contribution from injecting boundaries
!
! Revision 1.4  2002/11/08 21:34:49  haselbac
! Fixed bug in use of grid speed
!
! Revision 1.3  2002/09/09 15:51:56  haselbac
! global and mixtInput under regions, add grid speeds
!
! Revision 1.2  2002/07/25 18:29:48  haselbac
! Bug fix, bug only shows for cases with dummy cells
!
! Revision 1.1  2002/05/04 17:02:00  haselbac
! Initial revision
!
! ******************************************************************************







