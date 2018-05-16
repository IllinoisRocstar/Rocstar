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
! Purpose: Suite of routines to convert state vector.
!
! Description: None.
!
! Notes: 
!
! ******************************************************************************
!
! $Id: RFLU_ModBoundConvertCv.F90,v 1.3 2008/12/06 08:44:20 mtcampbe Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModBoundConvertCv

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
    
  IMPLICIT NONE
   
  PRIVATE
  PUBLIC :: RFLU_BXV_ConvertCvCons2Prim, & 
            RFLU_BXV_ConvertCvPrim2Cons
   
  SAVE    
     
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS




! ******************************************************************************
!
! Purpose: Convert conserved state vector to primitive variables.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to data of current region
!   pPatch              Pointer to data of current patch
!   cvStateFuture       Future state of conserved variables
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_BXV_ConvertCvCons2Prim(pRegion,pPatch,cvStateFuture)

  USE ModInterfaces, ONLY: MixtPerf_R_M, &
                           MixtPerf_T_DPR

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: cvStateFuture
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: icg,ifl,indMol
  REAL(RFREAL) :: gc,ir,mw,p,r
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_ConvertCvCons2Prim',&
  'RFLU_ModBoundConvertCv.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pCv   => pPatch%mixt%cv
  pDv   => pPatch%mixt%dv

  pGv   => pRegion%mixt%gv ! pGv taken from interior domain

  indMol = pRegion%mixtInput%indMol

! ******************************************************************************
! Actual conversion
! ******************************************************************************

  SELECT CASE ( pPatch%mixt%cvState )

! ==============================================================================
!   Convert from conservative to primitive form
! ==============================================================================

    CASE ( CV_MIXT_STATE_CONS )
      SELECT CASE ( cvStateFuture )

! ------------------------------------------------------------------------------
!       Convert to duvwp form
! ------------------------------------------------------------------------------

        CASE ( CV_MIXT_STATE_DUVWP )
          pPatch%mixt%cvState = CV_MIXT_STATE_DUVWP

          DO ifl = 1,pPatch%nBFaces
            ir = 1.0_RFREAL/pCv(CV_MIXT_DENS,ifl)

            pCv(CV_MIXT_XVEL,ifl) = ir*pCv(CV_MIXT_XMOM,ifl)
            pCv(CV_MIXT_YVEL,ifl) = ir*pCv(CV_MIXT_YMOM,ifl)
            pCv(CV_MIXT_ZVEL,ifl) = ir*pCv(CV_MIXT_ZMOM,ifl)
            pCv(CV_MIXT_PRES,ifl) = pDv(DV_MIXT_PRES,ifl)
          END DO ! ifl

! ------------------------------------------------------------------------------
!       Convert to duvwt form
! ------------------------------------------------------------------------------

        CASE (CV_MIXT_STATE_DUVWT)
          pPatch%mixt%cvState = CV_MIXT_STATE_DUVWT

          SELECT CASE ( pRegion%mixtInput%fluidModel )

! --------- Compressible fluid model -----------------------------------------

            CASE ( FLUID_MODEL_COMP )
              SELECT CASE ( pRegion%mixtInput%gasModel )

! ------------- TC perfect gas or mixture thereof, pseudo-gas

                CASE ( GAS_MODEL_TCPERF, &
                       GAS_MODEL_MIXT_TCPERF, &
                       GAS_MODEL_MIXT_PSEUDO )
                  DO ifl = 1,pPatch%nBFaces
                    icg = pPatch%bf2c(ifl)

                    r  = pCv(CV_MIXT_DENS,ifl)
                    p  = pDv(DV_MIXT_PRES,ifl)
                    ir = 1.0_RFREAL/r

                    pCv(CV_MIXT_XVEL,ifl) = ir*pCv(CV_MIXT_XMOM,ifl)
                    pCv(CV_MIXT_YVEL,ifl) = ir*pCv(CV_MIXT_YMOM,ifl)
                    pCv(CV_MIXT_ZVEL,ifl) = ir*pCv(CV_MIXT_ZMOM,ifl)

                    mw = pGv(GV_MIXT_MOL,indMol*icg)
                    gc = MixtPerf_R_M(mw)

                    pCv(CV_MIXT_TEMP,ifl) = MixtPerf_T_DPR(r,p,gc)
                  END DO ! ifl

! ------------- Gas-liquid mixture

                CASE ( GAS_MODEL_MIXT_GASLIQ )
                  DO ifl = 1,pPatch%nBFaces
                    icg = pPatch%bf2c(ifl)

                    ir = 1.0_RFREAL/pCv(CV_MIXT_DENS,ifl)

                    pCv(CV_MIXT_XVEL,ifl) = ir*pCv(CV_MIXT_XMOM,ifl)
                    pCv(CV_MIXT_YVEL,ifl) = ir*pCv(CV_MIXT_YMOM,ifl)
                    pCv(CV_MIXT_ZVEL,ifl) = ir*pCv(CV_MIXT_ZMOM,ifl)

                    pCv(CV_MIXT_TEMP,ifl) = pDv(DV_MIXT_TEMP,ifl)
                  END DO ! ifl
                CASE DEFAULT
                  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
              END SELECT ! pRegion%mixtInput%gasModel

! --------- Default ----------------------------------------------------------

            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pRegion%mixtInput%fluidModel

! ------------------------------------------------------------------------------
!       Default
! ------------------------------------------------------------------------------

        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! cvStateFuture

! ==============================================================================
!   Error - invalid input
! ==============================================================================

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pPatch%mixt%cvStateFuture

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_ConvertCvCons2Prim






! ******************************************************************************
!
! Purpose: Convert primitive state vector to consverved variables.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to data of current region
!   pPatch              Pointer to data of current patch
!   cvStateFuture       Future state of conserved variables
!
! Output: None.
!
! Notes:
!   1. Strictly speaking, cvStateFuture is not needed (there is only one
!      state for conserved variables), but kept for consistency with
!      RFLU_ConvertCvCons2Prim.
!
! ******************************************************************************

SUBROUTINE RFLU_BXV_ConvertCvPrim2Cons(pRegion,pPatch,cvStateFuture)

  USE ModInterfaces, ONLY: MixtPerf_Cv_CpR, &
                           MixtGasLiq_Eo_CvmTVm2, &
                           MixtPerf_Eo_DGPUVW, &
                           MixtPerf_G_CpR, &
                           MixtPerf_R_M

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: cvStateFuture
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: icg,ifl,indCp,indMol
  REAL(RFREAL) :: cp,Cvm,cvg,cvl,cvv,g,gc,mw,p,r,Rg,Rv,rYg,rYl,rYv,t,u,v, &
                  Vm2,w
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BXV_ConvertCvPrim2Cons',&
  'RFLU_ModBoundConvertCv.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pCv   => pPatch%mixt%cv
  pDv   => pPatch%mixt%dv

  pGv   => pRegion%mixt%gv ! pGv is taken from interior cells

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

! ******************************************************************************
! Actual conversion
! ******************************************************************************

  IF ( pPatch%mixt%cvState == CV_MIXT_STATE_DUVWP .OR. &
       pPatch%mixt%cvState == CV_MIXT_STATE_DUVWT ) THEN

! ==============================================================================
!   Convert from primitive to conservative form
! ==============================================================================

    SELECT CASE ( cvStateFuture )
      CASE ( CV_MIXT_STATE_CONS )
        pPatch%mixt%cvState = CV_MIXT_STATE_CONS

        SELECT CASE ( pRegion%mixtInput%gasModel )

! ------------------------------------------------------------------------------
!         Thermally and calorically perfect gas or pseudo-gas
! ------------------------------------------------------------------------------

          CASE ( GAS_MODEL_TCPERF, &
                 GAS_MODEL_MIXT_TCPERF, &
                 GAS_MODEL_MIXT_PSEUDO )
            DO ifl = 1,pPatch%nBFaces
              icg = pPatch%bf2c(ifl)

              r = pCv(CV_MIXT_DENS,ifl)
              u = pCv(CV_MIXT_XVEL,ifl)
              v = pCv(CV_MIXT_YVEL,ifl)
              w = pCv(CV_MIXT_ZVEL,ifl)
              p = pDv(DV_MIXT_PRES,ifl)

              pCv(CV_MIXT_XMOM,ifl) = r*u
              pCv(CV_MIXT_YMOM,ifl) = r*v
              pCv(CV_MIXT_ZMOM,ifl) = r*w

              cp = pGv(GV_MIXT_CP,indCp*icg)
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)

              pCv(CV_MIXT_ENER,ifl) = r*MixtPerf_Eo_DGPUVW(r,g,p,u,v,w)
            END DO ! icg

! ------------------------------------------------------------------------------
!         Mixture of gas, liquid, and vapor
! ------------------------------------------------------------------------------

          CASE ( GAS_MODEL_MIXT_GASLIQ )
            CALL ErrorStop(global,ERR_GASMODEL_INVALID,__LINE__, &
                           'Can only be used with species module.')

! ------------------------------------------------------------------------------
!         Other or invalid gas models
! ------------------------------------------------------------------------------

          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pRegion%mixtInput

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! cvStateFuture

! ==============================================================================
! Error - invalid input
! ==============================================================================

  ELSE
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END IF ! pPatch%mixt%cvState

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BXV_ConvertCvPrim2Cons






  
! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLU_ModBoundConvertCv


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModBoundConvertCv.F90,v $
! Revision 1.3  2008/12/06 08:44:20  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:31  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2006/08/19 15:37:43  mparmar
! Initial revision
!
! ******************************************************************************








