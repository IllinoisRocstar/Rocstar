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
! Purpose: Initialize flow field in a region from scratch.
!
! Description: Use user input to define initial flow field.
!
! Input: 
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_InitFlowScratch.F90,v 1.5 2008/12/06 08:44:56 mtcampbe Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_InitFlowScratch(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  
  USE ModInterfaces, ONLY: MixtGasLiq_Eo_CvmTVm2, &
                           MixtLiq_D_DoBpPPoBtTTo, &
                           MixtPerf_Cv_CpR, &
                           MixtPerf_D_PRT, &
                           MixtPerf_Eo_DGPUVW, &
                           MixtPerf_G_CpR, & 
                           MixtPerf_R_M, &                             
                           RFLU_GetCvLoc
  
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
  INTEGER :: cvMixtPres,cvMixtXVel,cvMixtYVel,cvMixtZVel,icg,indCp,indMol
  REAL(RFREAL) :: cp,Cvm,cvg,cvv,Eo,g,gc,mw,r,rhog,rhol,rhov,Rg,Rv,Vm2,Yg,Yl,Yv
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pGv
  TYPE(t_global), POINTER :: global
  TYPE(t_mixt_input), POINTER :: pMixtInput

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InitFlowScratch.F90,v $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_InitFlowScratch', &
                        'RFLU_InitFlowScratch.F90')
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Initializing flow field from scratch...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pCv        => pRegion%mixt%cv
  pGv        => pRegion%mixt%gv  
  pMixtInput => pRegion%mixtInput
  
  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol  
  
! ******************************************************************************
! Initialize flow field based on user input and fluid model
! ******************************************************************************

  SELECT CASE ( pMixtInput%fluidModel ) 
 
! ==============================================================================  
!   Incompressible fluid model  
! ==============================================================================  

    CASE ( FLUID_MODEL_INCOMP )     
      pRegion%mixt%cvState = CV_MIXT_STATE_PRIM

      cvMixtXVel = RFLU_GetCvLoc(global,pMixtInput%fluidModel,CV_MIXT_XVEL)
      cvMixtYVel = RFLU_GetCvLoc(global,pMixtInput%fluidModel,CV_MIXT_YVEL)
      cvMixtZVel = RFLU_GetCvLoc(global,pMixtInput%fluidModel,CV_MIXT_ZVEL)
      cvMixtPres = RFLU_GetCvLoc(global,pMixtInput%fluidModel,CV_MIXT_PRES)                

      DO icg = 1,pRegion%grid%nCellsTot
        pCv(cvMixtXVel,icg) = pRegion%mixtInput%iniVelX
        pCv(cvMixtYVel,icg) = pRegion%mixtInput%iniVelY
        pCv(cvMixtZVel,icg) = pRegion%mixtInput%iniVelZ
        pCv(cvMixtPres,icg) = pRegion%mixtInput%iniPress 
      END DO ! icg   
      
! ==============================================================================  
!   Compressible fluid model  
! ==============================================================================  
    
    CASE ( FLUID_MODEL_COMP )
      pRegion%mixt%cvState = CV_MIXT_STATE_CONS
          
      SELECT CASE ( pRegion%mixtInput%gasModel ) 
      
! ------------------------------------------------------------------------------
!       Gas models without liquid phase
! ------------------------------------------------------------------------------      
      
        CASE ( GAS_MODEL_TCPERF, & 
               GAS_MODEL_TPERF, & 
               GAS_MODEL_MIXT_TCPERF, & 
               GAS_MODEL_MIXT_TPERF, & 
               GAS_MODEL_MIXT_PSEUDO )     
          DO icg = 1,pRegion%grid%nCellsTot
            mw = pGv(GV_MIXT_MOL,indMol*icg)
            cp = pGv(GV_MIXT_CP ,indCp *icg)

            gc = MixtPerf_R_M(mw)
            g  = MixtPerf_G_CpR(cp,gc)

            r = pRegion%mixtInput%iniDens

            pCv(CV_MIXT_DENS,icg) = r
            pCv(CV_MIXT_XMOM,icg) = r*pRegion%mixtInput%iniVelX
            pCv(CV_MIXT_YMOM,icg) = r*pRegion%mixtInput%iniVelY
            pCv(CV_MIXT_ZMOM,icg) = r*pRegion%mixtInput%iniVelZ

            Eo = MixtPerf_Eo_DGPUVW(pCv(CV_MIXT_DENS,icg),g,    & 
                                    pRegion%mixtInput%iniPress, & 
                                    pRegion%mixtInput%iniVelX,  & 
                                    pRegion%mixtInput%iniVelY,  & 
                                    pRegion%mixtInput%iniVelZ)

            pCv(CV_MIXT_ENER,icg) = pCv(CV_MIXT_DENS,icg)*Eo
          END DO ! icg    

! ------------------------------------------------------------------------------
!       Mixture of gas, liquid, and vapor
! ------------------------------------------------------------------------------      

        CASE ( GAS_MODEL_MIXT_GASLIQ ) 
          Rg  = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
          cvg = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht,Rg)

          Rv  = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
          cvv = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht,Rv)

          rhol = MixtLiq_D_DoBpPPoBtTTo(global%refDensityLiq, &
                                        global%refBetaPLiq, &
                                        global%refBetaTLiq, &
                                        pRegion%mixtInput%iniPress, &
                                        global%refPressLiq, &
                                        pRegion%mixtInput%iniTemp, &
                                        global%refTempLiq)

          DO icg = 1,pRegion%grid%nCellsTot
            Yg = pRegion%specInput%specType(1)%initVal
            Yv = pRegion%specInput%specType(2)%initVal 
            Yl = 1.0_RFREAL - Yg - Yv        

            rhog = MixtPerf_D_PRT(pRegion%mixtInput%iniPress,Rg, &
                                  pRegion%mixtInput%iniTemp)
            rhov = MixtPerf_D_PRT(pRegion%mixtInput%iniPress,Rv, &
                                  pRegion%mixtInput%iniTemp) 

            r = 1.0_RFREAL/(Yg/rhog + Yv/rhov + Yl/rhol)

            pCv(CV_MIXT_DENS,icg) = r
            pCv(CV_MIXT_XMOM,icg) = r*pRegion%mixtInput%iniVelX
            pCv(CV_MIXT_YMOM,icg) = r*pRegion%mixtInput%iniVelY
            pCv(CV_MIXT_ZMOM,icg) = r*pRegion%mixtInput%iniVelZ        

            Cvm = Yg*cvg + Yv*cvv + Yl*global%refCvLiq 

            Vm2 = pRegion%mixtInput%iniVelX*pRegion%mixtInput%iniVelX &
                + pRegion%mixtInput%iniVelY*pRegion%mixtInput%iniVelY &
                + pRegion%mixtInput%iniVelZ*pRegion%mixtInput%iniVelZ

            pCv(CV_MIXT_ENER,icg) = pCv(CV_MIXT_DENS,icg)* &
             MixtGasLiq_Eo_CvmTVm2(Cvm,pRegion%mixtInput%iniTemp,Vm2) 
          END DO ! icg

! ------------------------------------------------------------------------------
!       Default
! ------------------------------------------------------------------------------      

        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%gasModel
      
! ==============================================================================  
!   Default
! ==============================================================================  

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%fluidModel
  
! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Initializing flow field from scratch done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InitFlowScratch

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InitFlowScratch.F90,v $
! Revision 1.5  2008/12/06 08:44:56  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:18:06  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/03/26 20:22:26  haselbac
! Added support for GL model
!
! Revision 1.2  2005/11/10 02:43:23  haselbac
! Added support for variable properties
!
! Revision 1.1  2005/04/15 15:08:18  haselbac
! Initial revision
!
! ******************************************************************************







