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
! Purpose: Initialize boundary-condition data using hard-coded values.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine assumes a perfect gas.
!
! ******************************************************************************
!
! $Id: RFLU_InitBcDataHardCode.F90,v 1.7 2008/12/06 08:44:56 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_InitBcDataHardCode(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModParameters
    
  USE RFLU_ModExactFlow, ONLY: RFLU_ComputeExactFlowProudman, & 
                               RFLU_ComputeExactFlowRingleb, & 
                               RFLU_ComputeExactFlowSsVortex      
  USE RFLU_ModFlowHardCode, ONLY: RFLU_GetParamsHardCodeProudman, & 
                                  RFLU_GetParamsHardCodeRingleb, & 
                                  RFLU_GetParamsHardCodeSsVortex    
    
  USE ModInterfaces, ONLY: MixtPerf_C_DGP, &
                           MixtPerf_P_DRT, &
                           MixtPerf_R_CpG, &
                           MixtPerf_T_DPR 
         
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

  CHARACTER(CHRLEN) :: errorString,RCSIdentString
  INTEGER :: ifl,iPatch
  REAL(RFREAL) :: betah,cpGas,d,dInc,gGas,height,mach,Mi,mInj,p,pTot,rGas, &
                  ri,t,tTot,u,v,vInj,vMag,w,x,y
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InitBcDataHardCode.F90,v $ $Revision: 1.7 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_InitBcDataHardCode', &
                        'RFLU_InitBcDataHardCode.F90')
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Initializing boundary-condition data from hard code...'
                             
    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,A)') SOLVER_NAME,'Case: ',TRIM(global%casename)
    END IF ! global%verbLevel                         
  END IF ! global%verbLevel

! ==============================================================================
! Set pointers
! ==============================================================================

  pGrid => pRegion%grid

! ==============================================================================
! Set constants
! ==============================================================================

  cpGas = global%refCp
  gGas  = global%refGamma  
  rGas  = MixtPerf_R_CpG(cpGas,gGas)
  
! ******************************************************************************
! Initialize data fields based on user input
! ******************************************************************************

  SELECT CASE ( global%casename )
  
! ==============================================================================
!   Proudman-Culick flow. NOTE this problem is two-dimensional and assumed to 
!   lie in the x-y plane, and that the injection boundary is located at 
!   y = -height.
! ==============================================================================  
  
    CASE ( "onera_c0", "onera_c0_2d_100x50" )
      CALL RFLU_GetParamsHardCodeProudman(dInc,mInj,vInj,pTot)
    
      height = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVertTot))
      
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)
        
        IF ( pPatch%bcType == BC_INJECTION ) THEN
          IF ( pPatch%mixt%distrib == BCDAT_DISTRIB ) THEN  
            DO ifl = 1,pPatch%nBFaces
              x = pPatch%fc(XCOORD,ifl)
              y = pPatch%fc(YCOORD,ifl)

              CALL RFLU_ComputeExactFlowProudman(global,x,y,height,dInc, &
                                                 vInj,pTot,d,u,v,w,p)

              t = MixtPerf_T_DPR(d,p,rGas)

              pPatch%mixt%vals(BCDAT_INJECT_MFRATE,ifl) = d*v
              pPatch%mixt%vals(BCDAT_INJECT_TEMP  ,ifl) = t                                  
            END DO ! ifl
          ELSE 
            WRITE(errorString,'(A,1X,I3)') 'Patch:',iPatch
            CALL ErrorStop(global,ERR_DISTRIB_INVALID,__LINE__, &
                           TRIM(errorString))
          END IF ! pPatch%mixt%distrib
        ELSE IF ( pPatch%bcType == BC_OUTFLOW ) THEN
          IF ( pPatch%mixt%distrib == BCDAT_DISTRIB ) THEN  
            DO ifl = 1,pPatch%nBFaces
              x = pPatch%fc(XCOORD,ifl)
              y = pPatch%fc(YCOORD,ifl)

              CALL RFLU_ComputeExactFlowProudman(global,x,y,height,dInc, &
                                                 vInj,pTot,d,u,v,w,p)

              pPatch%mixt%vals(BCDAT_OUTFLOW_PRESS,ifl) = p
            END DO ! ifl     
          ELSE 
            WRITE(errorString,'(A,1X,I3)') 'Patch:',iPatch
            CALL ErrorStop(global,ERR_DISTRIB_INVALID,__LINE__, &
                           TRIM(errorString))          
          END IF ! pPatch%mixt%distrib              
        END IF ! pPatch%bcType      
      END DO ! iPatch  
  
! ==============================================================================
!   Ringleb flow. NOTE this problem is two-dimensional and assumed to lie in 
!   the x-y plane and that the exact solution is restricted to gamma = 1.4.
! ==============================================================================  
  
    CASE ( "ringleb" )
      CALL RFLU_GetParamsHardCodeRingleb(pTot,tTot)
          
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)
        
        IF ( pPatch%bcType == BC_INFLOW_TOTANG ) THEN
          IF ( pPatch%mixt%distrib == BCDAT_DISTRIB ) THEN  
            DO ifl = 1,pPatch%nBFaces
              x = pPatch%fc(XCOORD,ifl)
              y = pPatch%fc(YCOORD,ifl)

              CALL RFLU_ComputeExactFlowRingleb(x,y,rGas,pTot,tTot,d,u,v,w,p)

              betah = ATAN2(v,u) ! w assumed to be zero

              pPatch%mixt%vals(BCDAT_INFLOW_PTOT ,ifl) = pTot
              pPatch%mixt%vals(BCDAT_INFLOW_TTOT ,ifl) = tTot 
              pPatch%mixt%vals(BCDAT_INFLOW_BETAH,ifl) = betah
              pPatch%mixt%vals(BCDAT_INFLOW_BETAV,ifl) = 0.0_RFREAL            
            END DO ! ifl
          ELSE 
            WRITE(errorString,'(A,1X,I3)') 'Patch:',iPatch
            CALL ErrorStop(global,ERR_DISTRIB_INVALID,__LINE__, &
                           TRIM(errorString))
          END IF ! pPatch%mixt%distrib
        ELSE IF ( pPatch%bcType == BC_OUTFLOW ) THEN
          IF ( pPatch%mixt%distrib == BCDAT_DISTRIB ) THEN  
            DO ifl = 1,pPatch%nBFaces
              x = pPatch%fc(XCOORD,ifl)
              y = pPatch%fc(YCOORD,ifl)

              CALL RFLU_ComputeExactFlowRingleb(x,y,rGas,pTot,tTot,d,u,v,w,p)

              pPatch%mixt%vals(BCDAT_OUTFLOW_PRESS,ifl) = p
            END DO ! ifl     
          ELSE 
            WRITE(errorString,'(A,1X,I3)') 'Patch:',iPatch
            CALL ErrorStop(global,ERR_DISTRIB_INVALID,__LINE__, &
                           TRIM(errorString))          
          END IF ! pPatch%mixt%distrib              
        END IF ! pPatch%bcType      
      END DO ! iPatch
      
! ==============================================================================
!   Supersonic vortex flow. NOTE this problem is two-dimensional and assumed 
!   to lie in the x-y plane. 
! ==============================================================================  
  
    CASE ( "ssvorth20x5l1"   ,"ssvortp20x5l1",    & 
           "ssvorth20x5l3"   ,"ssvortp20x5l3",    &
           "ssvorth40x10l1"  ,"ssvortp40x10l1",   & 
           "ssvorth40x10l3"  ,"ssvortp40x10l3",   & 
           "ssvorth80x20l1"  ,"ssvortp80x20l1",   & 
           "ssvorth80x20l3"  ,"ssvortp80x20l3",   & 
           "ssvorth160x40l1" ,"ssvortp160x40l1",  &
           "ssvorth160x40l3" ,"ssvortp160x40l3",  &
           "ssvorth320x80l1" ,"ssvortp320x80l1",  &
           "ssvorth320x80l3" ,"ssvortp320x80l3",  &
           "ssvorth640x160l1","ssvortp640x160l1", &
           "ssvorth640x160l3","ssvortp640x160l3" )
      CALL RFLU_GetParamsHardCodeSsVortex(ri,Mi,pTot,tTot)
           
      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)
        
        IF ( pPatch%bcType == BC_INFLOW_TOTANG ) THEN
          IF ( pPatch%mixt%distrib == BCDAT_DISTRIB ) THEN 
            DO ifl = 1,pPatch%nBFaces
              x = pPatch%fc(XCOORD,ifl)
              y = pPatch%fc(YCOORD,ifl)

              CALL RFLU_ComputeExactFlowSsVortex(x,y,gGas,rGas,ri,Mi,pTot, & 
                                                 tTot,d,u,v,w,p)

              vMag  = SQRT(u*u + v*v) ! w assumed to be zero
              betah = ASIN(v/vMag) 
              mach  = vMag/MixtPerf_C_DGP(d,gGas,p)

              pPatch%mixt%vals(BCDAT_INFLOW_PTOT ,ifl) = pTot
              pPatch%mixt%vals(BCDAT_INFLOW_TTOT ,ifl) = tTot 
              pPatch%mixt%vals(BCDAT_INFLOW_BETAH,ifl) = betah
              pPatch%mixt%vals(BCDAT_INFLOW_BETAV,ifl) = 0.0_RFREAL
              pPatch%mixt%vals(BCDAT_INFLOW_MACH ,ifl) = mach                       
            END DO ! ifl  
          ELSE 
            WRITE(errorString,'(A,1X,I3)') 'Patch:',iPatch
            CALL ErrorStop(global,ERR_DISTRIB_INVALID,__LINE__, &
                           TRIM(errorString))            
          END IF ! pPatch%mixt%distrib        
        END IF ! pPatch%bcType      
      END DO ! iPatch     
        
! ==============================================================================
!   Default - must be due to input error
! ==============================================================================  
            
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
  END SELECT ! global%casename

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Initializing boundary-condition data from hard code done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InitBcDataHardCode

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InitBcDataHardCode.F90,v $
! Revision 1.7  2008/12/06 08:44:56  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:18:06  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/08/19 15:41:09  mparmar
! Renamed patch variables
!
! Revision 1.4  2006/03/08 23:39:39  haselbac
! Added 1 layer ssvort cases
!
! Revision 1.3  2005/10/09 15:40:19  haselbac
! Added 2d C0 case
!
! Revision 1.2  2005/04/27 02:14:52  haselbac
! Adapted to changes in inflow treatment
!
! Revision 1.1  2005/04/15 15:08:13  haselbac
! Initial revision
!
! Revision 1.5  2004/10/19 19:30:53  haselbac
! Removed setting of relative velocity on injection boundaries
!
! Revision 1.4  2004/07/06 15:15:53  haselbac
! Adapted to changes in libflu and modflu, cosmetics
!
! Revision 1.3  2004/02/23 23:05:10  haselbac
! Added Proudman solution for ONERA C0 case
!
! Revision 1.2  2004/02/13 03:05:54  haselbac
! Added more casenames
!
! Revision 1.1  2004/01/29 22:58:49  haselbac
! Initial revision
!
! ******************************************************************************







