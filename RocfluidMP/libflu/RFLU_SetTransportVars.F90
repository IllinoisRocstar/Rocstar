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
! Purpose: Set transport variables.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!  icgBeg       Beginning cell index
!  icgEnd       Ending cell index
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_SetTransportVars.F90,v 1.6 2008/12/06 08:44:13 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetTransportVars(pRegion,icgBeg,icgEnd)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
                                 
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(IN) :: icgBeg,icgEnd
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: icg,indCp
  REAL(RFREAL) :: absTemp,iPrLam,iPrLamLiq,kLiq,kGas,kVap,muLiq,muGas,muVap,&
                  refTemp,refVisc,refViscLiq,suthCoef,s1,s2,s3,s4,term,y,z
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pTv
#ifdef SPEC
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvSpec
#endif
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetTransportVars.F90,v $ $Revision: 1.6 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SetTransportVars',&
  'RFLU_SetTransportVars.F90')

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************
 
  pCv => pRegion%mixt%cv
  pDv => pRegion%mixt%dv 
  pGv => pRegion%mixt%gv 
  pTv => pRegion%mixt%tv

#ifdef SPEC
  pCvSpec => pRegion%spec%cv
#endif

  indCp = pRegion%mixtInput%indCp

  iPrLam   = 1.0_RFREAL/pRegion%mixtInput%prLam
  refVisc  = pRegion%mixtInput%refVisc
  refTemp  = pRegion%mixtInput%refTemp
  suthCoef = pRegion%mixtInput%suthCoef

! ******************************************************************************  
! Compute dependent variables
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%fluidModel )
  
! ==============================================================================  
!   Incompressible fluid model
! ==============================================================================  
   
    CASE ( FLUID_MODEL_INCOMP ) 

! ==============================================================================  
!   Compressible fluid model
! ==============================================================================  

    CASE ( FLUID_MODEL_COMP )             
      SELECT CASE ( pRegion%mixtInput%gasModel ) 
      
! ------------------------------------------------------------------------------
!       Thermally and calorically perfect or thermally perfect gas
! ------------------------------------------------------------------------------      
            
        CASE ( GAS_MODEL_TCPERF, & 
               GAS_MODEL_TPERF )        
          SELECT CASE ( pRegion%mixtInput%viscModel ) 
      
! -------- Sutherland model ----------------------------------------------------
      
            CASE ( VISC_SUTHR )
              DO icg = icgBeg,icgEnd
                term = SQRT(pDv(DV_MIXT_TEMP,icg)/suthCoef) & 
                      *(1.0_RFREAL + refTemp/suthCoef)/ & 
                       (1.0_RFREAL + refTemp/pDv(DV_MIXT_TEMP,icg))

                pTv(TV_MIXT_MUEL,icg) = term*refVisc
                pTv(TV_MIXT_TCOL,icg) = pGv(GV_MIXT_CP  ,icg*indCp)* & 
                                        pTv(TV_MIXT_MUEL,icg      )*iPrLam
              END DO ! icg

! --------- Constant viscosity model -------------------------------------------
                
            CASE ( VISC_FIXED ) 
              DO icg = icgBeg,icgEnd
                pTv(TV_MIXT_MUEL,icg) = refVisc
                pTv(TV_MIXT_TCOL,icg) = pGv(GV_MIXT_CP  ,icg*indCp) & 
                                       *pTv(TV_MIXT_MUEL,icg      )*iPrLam
              END DO ! icg           
        
! --------- Antibes model ------------------------------------------------------
              
            CASE ( VISC_ANTIB ) 
              s1 = 110.0_RFREAL
              s2 = 120.0_RFREAL
              s3 = 230.0_RFREAL
              s4 = 1.0_RFREAL/refTemp*SQRT(s2/refTemp)*(refTemp+s1)/s3

              IF ( refTemp <= s2 ) THEN
                DO icg = icgBeg,icgEnd
                  absTemp = ABS(pDv(DV_MIXT_TEMP,icg))

                  IF ( absTemp < s2 ) THEN
                    term = absTemp/refTemp
                  ELSE
                    term = SQRT(absTemp/s2)*(absTemp/refTemp)*(s3/(absTemp+s1))
                  END IF ! absTemp

                  pTv(TV_MIXT_MUEL,icg) = refVisc*term
                  pTv(TV_MIXT_TCOL,icg) = pGv(GV_MIXT_CP  ,icg*indCp) & 
                                         *pTv(TV_MIXT_MUEL,icg      )*iPrLam
                END DO ! icg
              ELSE 
                DO icg = icgBeg,icgEnd
                  absTemp = ABS(pDv(DV_MIXT_TEMP,icg))

                  IF ( absTemp < s2 ) THEN
                    term = s4
                  ELSE
                    term = SQRT(absTemp/refTemp)*(absTemp/refTemp) & 
                          *(refTemp+s1)/(absTemp+s1)
                  END IF ! absTemp

                  pTv(TV_MIXT_MUEL,icg) = term*refVisc
                  pTv(TV_MIXT_TCOL,icg) = pGv(GV_MIXT_CP  ,icg*indCp) & 
                                         *pTv(TV_MIXT_MUEL,icg      )*iPrLam
                END DO ! icg      
              END IF ! refTemp         

! --------- Default ------------------------------------------------------------

            CASE DEFAULT 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pRegion%mixtInput%viscModel
      
! ------------------------------------------------------------------------------
!       Gas-liquid mixture fluid model
! ------------------------------------------------------------------------------

        CASE ( GAS_MODEL_MIXT_GASLIQ )
#ifdef SPEC        
          DO icg = icgBeg,icgEnd

! --------- Use Sutherland formula for gas -------------------------------------

            term = SQRT(pDv(DV_MIXT_TEMP,icg)/suthCoef) & 
                   *(1.0_RFREAL + refTemp/suthCoef)/ & 
                    (1.0_RFREAL + refTemp/pDv(DV_MIXT_TEMP,icg))

            muGas = term*refVisc
            kGas  = pRegion%specInput%specType(1)%pMaterial%spht*muGas*iPrLam
   
! --------- Use Sutherland formula for vapor -----------------------------------

            term = SQRT(pDv(DV_MIXT_TEMP,icg)/suthCoef) & 
                   *(1.0_RFREAL + refTemp/suthCoef)/ & 
                    (1.0_RFREAL + refTemp/pDv(DV_MIXT_TEMP,icg)) 

            muVap = term*refVisc   
            kVap  = pRegion%specInput%specType(2)%pMaterial%spht*muVap*iPrLam

! --------- Compute liquid viscosity and conductivity --------------------------

            refViscLiq = 1.0E-02_RFREAL ! NOTE hard-coded for now
            z = refTemp/pDv(DV_MIXT_TEMP,icg)
            y = -1.704_RFREAL -5.306_RFREAL*z + 7.003_RFREAL*z*z

            muLiq = refViscLiq*EXP(y) 
            kliq  = 0.585_RFREAL

! --------- Compute mixture viscosity and conductivity --------------------------

            pTv(TV_MIXT_MUEL,icg) = & 
              ((pCv(CV_MIXT_DENS,icg)- pCvSpec(1,icg) &
              - pCvSpec(2,icg))*muLiq + pCvSpec(1,icg)*muGas & 
              + pCvSpec(2,icg)*muVap)/pCv(CV_MIXT_DENS,icg)
            pTv(TV_MIXT_TCOL,icg) = & 
              ((pCv(CV_MIXT_DENS,icg)- pCvSpec(1,icg) &
              - pCvSpec(2,icg))*kLiq + pCvSpec(1,icg)*kGas & 
              + pCvSpec(2,icg)*kVap)/pCv(CV_MIXT_DENS,icg)
          END DO ! icg
#else
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
#endif          
            
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

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetTransportVars

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetTransportVars.F90,v $
! Revision 1.6  2008/12/06 08:44:13  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.3  2006/03/26 20:21:45  haselbac
! Added support for GL model
!
! Revision 1.2  2005/04/15 15:06:22  haselbac
! Added range arguments, removed pGrid declaration
!
! Revision 1.1  2004/11/14 20:02:49  haselbac
! Initial revision
!
! ******************************************************************************







