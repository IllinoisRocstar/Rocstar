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
! Purpose: Calculate cell gradients for higher-order schemes in RocfluidMP
!   framework.
!
! Description: None.
!
! Input:
!   region         Data of current region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: CellGradientsMP.F90,v 1.12 2008/12/06 08:44:08 mtcampbe Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE CellGradientsMP(region)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region

#ifdef RFLU
  USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, &
                               RFLU_ConvertCvPrim2Cons, &
                               RFLU_ScalarConvertCvCons2Prim, &
                               RFLU_ScalarConvertCvPrim2Cons
  USE RFLU_ModDifferentiationCells
  USE RFLU_ModLimiters, ONLY: RFLU_CreateLimiter, &
                              RFLU_ComputeLimiterBarthJesp, &
                              RFLU_ComputeLimiterVenkat, &
                              RFLU_LimitGradCells, &
                              RFLU_LimitGradCellsSimple, &
                              RFLU_DestroyLimiter
  USE RFLU_ModWENO, ONLY: RFLU_WENOGradCellsWrapper, & 
                          RFLU_WENOGradCellsXYZWrapper                              
#endif

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), TARGET :: region

! ==============================================================================
! Locals
! ==============================================================================

  TYPE(t_global), POINTER :: global

#ifdef RFLU
  INTEGER :: varInfoMixt(CV_MIXT_DENS:CV_MIXT_PRES)
  TYPE(t_region), POINTER :: pRegion
#ifdef SPEC
  INTEGER :: errorFlag,iSpec
  INTEGER, DIMENSION(:), ALLOCATABLE :: varInfoSpec
#endif
#endif

! ******************************************************************************
! Start
! ******************************************************************************

  global => region%global

  CALL RegisterFunction(global,'CellGradientsMP',&
  'CellGradientsMP.F90')

#ifdef RFLU
! ******************************************************************************
! Compute cell gradients
! ******************************************************************************

  pRegion => region

  IF ( pRegion%mixtInput%spaceOrder > 1 ) THEN

! ==============================================================================
!   Mixture
! ==============================================================================

! ------------------------------------------------------------------------------ 
!   Convert to primitive state vector
! ------------------------------------------------------------------------------

    SELECT CASE ( pRegion%mixtInput%gasModel )     
      CASE ( GAS_MODEL_TCPERF,      &
             GAS_MODEL_TPERF,       & 
             GAS_MODEL_MIXT_TCPERF, & 
             GAS_MODEL_MIXT_TPERF,  & 
             GAS_MODEL_MIXT_PSEUDO ) 
        CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)

        varInfoMixt(CV_MIXT_DENS) = V_MIXT_DENS
        varInfoMixt(CV_MIXT_XVEL) = V_MIXT_XVEL
        varInfoMixt(CV_MIXT_YVEL) = V_MIXT_YVEL       
        varInfoMixt(CV_MIXT_ZVEL) = V_MIXT_ZVEL
        varInfoMixt(CV_MIXT_PRES) = V_MIXT_PRES           
      CASE ( GAS_MODEL_MIXT_GASLIQ )
        CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWT)
        
        varInfoMixt(CV_MIXT_DENS) = V_MIXT_DENS
        varInfoMixt(CV_MIXT_XVEL) = V_MIXT_XVEL
        varInfoMixt(CV_MIXT_YVEL) = V_MIXT_YVEL       
        varInfoMixt(CV_MIXT_ZVEL) = V_MIXT_ZVEL
        varInfoMixt(CV_MIXT_TEMP) = V_MIXT_TEMP               
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%gasModel

! ------------------------------------------------------------------------------ 
!   Compute gradients
! ------------------------------------------------------------------------------
                                                               
    CALL RFLU_ComputeGradCellsWrapper(pRegion,CV_MIXT_DENS,CV_MIXT_PRES, &
                                      GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                      varInfoMixt,pRegion%mixt%cv, &
                                      pRegion%mixt%gradCell)
                                      
! ------------------------------------------------------------------------------ 
!   Modify gradients
! ------------------------------------------------------------------------------

    SELECT CASE ( pRegion%mixtInput%reconst ) 
      CASE ( RECONST_NONE ) 
      CASE ( RECONST_WENO_SIMPLE ) 
        CALL RFLU_WENOGradCellsWrapper(pRegion,GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                       pRegion%mixt%gradCell)
        CALL RFLU_LimitGradCellsSimple(pRegion,CV_MIXT_DENS,CV_MIXT_PRES, &
                                      GRC_MIXT_DENS,GRC_MIXT_PRES, & 
                                      pRegion%mixt%cv,pRegion%mixt%cvInfo, &
                                      pRegion%mixt%gradCell)
      CASE ( RECONST_WENO_XYZ ) 
        CALL RFLU_WENOGradCellsXYZWrapper(pRegion,GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                          pRegion%mixt%gradCell)
        CALL RFLU_LimitGradCellsSimple(pRegion,CV_MIXT_DENS,CV_MIXT_PRES, &
                                      GRC_MIXT_DENS,GRC_MIXT_PRES, & 
                                      pRegion%mixt%cv,pRegion%mixt%cvInfo, &
                                      pRegion%mixt%gradCell)
      CASE ( RECONST_LIM_BARTHJESP )  
        CALL RFLU_CreateLimiter(pRegion,GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                pRegion%mixt%lim)
        CALL RFLU_ComputeLimiterBarthJesp(pRegion,CV_MIXT_DENS,CV_MIXT_PRES, &
                                          GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                          pRegion%mixt%cv, &
                                          pRegion%mixt%gradCell, &
                                          pRegion%mixt%lim)
        CALL RFLU_LimitGradCells(pRegion,GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                 pRegion%mixt%gradCell,pRegion%mixt%lim)
        CALL RFLU_DestroyLimiter(pRegion,pRegion%mixt%lim)
      CASE ( RECONST_LIM_VENKAT )  
        CALL RFLU_CreateLimiter(pRegion,GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                pRegion%mixt%lim)
        CALL RFLU_ComputeLimiterVenkat(pRegion,CV_MIXT_DENS,CV_MIXT_PRES, &
                                       GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                       pRegion%mixt%cv,pRegion%mixt%gradCell, &
                                       pRegion%mixt%lim)
        CALL RFLU_LimitGradCells(pRegion,GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                 pRegion%mixt%gradCell,pRegion%mixt%lim)
        CALL RFLU_DestroyLimiter(pRegion,pRegion%mixt%lim)
      CASE DEFAULT 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%reconst

! ------------------------------------------------------------------------------ 
!   Convert back to primitive state vector
! ------------------------------------------------------------------------------

    CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)

#ifdef SPEC
! ==============================================================================
!   Species
! ==============================================================================

    IF ( global%specUsed .EQV. .TRUE. ) THEN
    
! ------------------------------------------------------------------------------ 
!     Convert to primitive state vector
! ------------------------------------------------------------------------------
    
      CALL RFLU_ScalarConvertCvCons2Prim(pRegion,pRegion%spec%cv, &
                                         pRegion%spec%cvState)

      ALLOCATE(varInfoSpec(pRegion%specInput%nSpecies),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'varInfoSpec')
      END IF ! global%error

      DO iSpec = 1,pRegion%specInput%nSpecies
        varInfoSpec(iSpec) = V_SPEC_VAR1 + iSpec - 1
      END DO ! iSpec

! ------------------------------------------------------------------------------ 
!     Compute gradients
! ------------------------------------------------------------------------------
                                                               
      CALL RFLU_ComputeGradCellsWrapper(pRegion,1,pRegion%specInput%nSpecies, &
                                        1,pRegion%specInput%nSpecies, &
                                        varInfoSpec,pRegion%spec%cv, &
                                        pRegion%spec%gradCell)

      DEALLOCATE(varInfoSpec,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'varInfoSpec')
      END IF ! global%error

! ------------------------------------------------------------------------------ 
!     Modify gradients
! ------------------------------------------------------------------------------

      SELECT CASE ( pRegion%mixtInput%reconst ) 
        CASE ( RECONST_NONE )
        CASE ( RECONST_WENO_SIMPLE ) 
          CALL RFLU_WENOGradCellsWrapper(pRegion,1,pRegion%specInput%nSpecies, &
                                         pRegion%spec%gradCell)
          CALL RFLU_LimitGradCellsSimple(pRegion,1, & 
                                         pRegion%specInput%nSpecies,1, & 
                                         pRegion%specInput%nSpecies, & 
                                         pRegion%spec%cv,pRegion%spec%cvInfo, & 
                                         pRegion%spec%gradCell)
        CASE ( RECONST_WENO_XYZ ) 
          CALL RFLU_WENOGradCellsXYZWrapper(pRegion,1, & 
                                            pRegion%specInput%nSpecies, &
                                            pRegion%spec%gradCell)
          CALL RFLU_LimitGradCellsSimple(pRegion,1, & 
                                         pRegion%specInput%nSpecies,1, & 
                                         pRegion%specInput%nSpecies, & 
                                         pRegion%spec%cv,pRegion%spec%cvInfo, & 
                                         pRegion%spec%gradCell)
        CASE ( RECONST_LIM_BARTHJESP )  
          CALL RFLU_CreateLimiter(pRegion,1,pRegion%specInput%nSpecies, &
                                  pRegion%spec%lim)
          CALL RFLU_ComputeLimiterBarthJesp(pRegion,1, & 
                                            pRegion%specInput%nSpecies,1, & 
                                            pRegion%specInput%nSpecies, & 
                                            pRegion%spec%cv, & 
                                            pRegion%spec%gradCell, & 
                                            pRegion%spec%lim)
          CALL RFLU_LimitGradCells(pRegion,1,pRegion%specInput%nSpecies, &
                                   pRegion%spec%gradCell,pRegion%spec%lim)
          CALL RFLU_DestroyLimiter(pRegion,pRegion%spec%lim)
        CASE ( RECONST_LIM_VENKAT )  
          CALL RFLU_CreateLimiter(pRegion,1,pRegion%specInput%nSpecies, &
                                  pRegion%spec%lim)
          CALL RFLU_ComputeLimiterVenkat(pRegion,1, & 
                                         pRegion%specInput%nSpecies,1, & 
                                         pRegion%specInput%nSpecies, & 
                                         pRegion%spec%cv, & 
                                         pRegion%spec%gradCell, & 
                                         pRegion%spec%lim)
          CALL RFLU_LimitGradCells(pRegion,1,pRegion%specInput%nSpecies, &
                                   pRegion%spec%gradCell,pRegion%spec%lim)
          CALL RFLU_DestroyLimiter(pRegion,pRegion%spec%lim)
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! pRegion%mixtInput%reconst

      CALL RFLU_ScalarConvertCvPrim2Cons(pRegion,pRegion%spec%cv, &
                                         pRegion%spec%cvState)
    END IF ! global%specUsed
#endif
  END IF ! pRegion%mixtInput%spaceOrder
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE CellGradientsMP

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: CellGradientsMP.F90,v $
! Revision 1.12  2008/12/06 08:44:08  mtcampbe
! Updated license.
!
! Revision 1.11  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.10  2006/04/27 15:10:04  haselbac
! Added VENKAT limiter option
!
! Revision 1.9  2006/04/15 16:52:30  haselbac
! Added RECONST_NONE option in CASE statements
!
! Revision 1.8  2006/04/07 14:37:20  haselbac
! Changes bcos of new WENO module and renamed routines
!
! Revision 1.7  2006/03/26 20:21:11  haselbac
! Introduced GL model, affects cv conversion
!
! Revision 1.6  2006/01/06 22:03:36  haselbac
! Added call to wrapper
!
! Revision 1.5  2005/10/27 18:52:58  haselbac
! Adapted to new cell grad routines
!
! Revision 1.4  2005/10/05 13:45:07  haselbac
! Revamped gradient computation
!
! Revision 1.3  2005/07/11 20:24:36  haselbac
! Added calls to WENO scheme with separate weighting in xyz directions
!
! Revision 1.2  2005/07/11 19:21:15  mparmar
! Added limiter reconstruction option
!
! Revision 1.1  2004/12/01 16:48:20  haselbac
! Initial revision after changing case
!
! Revision 1.3  2004/07/28 15:29:18  jferry
! created global variable for spec use
!
! Revision 1.2  2004/01/29 22:52:37  haselbac
! Added second-order for species
!
! Revision 1.1  2003/12/04 03:22:54  haselbac
! Initial revision
!
! ******************************************************************************







