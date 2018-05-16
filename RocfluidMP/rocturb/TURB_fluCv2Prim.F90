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
! Purpose: Convert conserved state vector to primitive variables.
!
! Description: None.
!
! Input:  region        = Pointer to data of current region
!         cvStateFuture = Future state of conserved variables
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
  
SUBROUTINE TURB_FluCv2Prim( pRegion,cvStateFuture )

  USE ModDataTypes
  USE ModGlobal, ONLY    : t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY      : t_grid
  USE ModInterfaces, ONLY: MixtPerf_R_M, MixtPerf_T_DPR
  USE ModParameters
  USE ModError
  IMPLICIT NONE


! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region), TARGET :: pRegion
  INTEGER, INTENT(IN) :: cvStateFuture

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: ic,indMol
  REAL(RFREAL) :: ir,mol,p,r,rgas
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
  TYPE(t_grid),   POINTER :: pGrid
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: TURB_fluCv2Prim.F90,v $ $Revision: 1.4 $'

  global => pRegion%global

  CALL RegisterFunction(global,'TURB_FluCv2Prim',&
  'TURB_fluCv2Prim.F90')

! *****************************************************************************
! Set pointers and variables
! *****************************************************************************

  pGrid => pRegion%grid
  pCv   => pRegion%mixt%cv
  pDv   => pRegion%mixt%dv
  pGv   => pRegion%mixt%gv

  indMol = pRegion%mixtInput%indMol

! *****************************************************************************
! Actual conversion
! *****************************************************************************

  SELECT CASE (pRegion%mixt%cvState)

! =============================================================================
!   Convert from conservative to primitive form 
! =============================================================================

    CASE (CV_MIXT_STATE_CONS) 
      SELECT CASE (cvStateFuture)

! -----------------------------------------------------------------------------      
!       Convert to duvwp form        
! -----------------------------------------------------------------------------

        CASE (CV_MIXT_STATE_DUVWP)
          pRegion%mixt%cvState = CV_MIXT_STATE_DUVWP          
          
          DO ic = 1,pGrid%nCellsTot
            ir = 1.0_RFREAL/pCv(CV_MIXT_DENS,ic)

            pCv(CV_MIXT_XVEL,ic) = ir*pCv(CV_MIXT_XMOM,ic)
            pCv(CV_MIXT_YVEL,ic) = ir*pCv(CV_MIXT_YMOM,ic)
            pCv(CV_MIXT_ZVEL,ic) = ir*pCv(CV_MIXT_ZMOM,ic)

            pCv(CV_MIXT_PRES,ic) = pDv(DV_MIXT_PRES,ic)
          END DO ! ic

! -----------------------------------------------------------------------------          
!       Convert to duvwt form             
! -----------------------------------------------------------------------------

        CASE (CV_MIXT_STATE_DUVWT) 
          pRegion%mixt%cvState = CV_MIXT_STATE_DUVWT                
                            
          SELECT CASE (pRegion%mixtInput%gasModel)

! --------- Perfect gas          

            CASE (GAS_MODEL_TCPERF)
              DO ic = 1,pGrid%nCellsTot                
                r  = pCv(CV_MIXT_DENS,ic)
                p  = pDv(DV_MIXT_PRES,ic)                
                ir = 1.0_RFREAL/r

                pCv(CV_MIXT_XVEL,ic) = ir*pCv(CV_MIXT_XMOM,ic)
                pCv(CV_MIXT_YVEL,ic) = ir*pCv(CV_MIXT_YMOM,ic)
                pCv(CV_MIXT_ZVEL,ic) = ir*pCv(CV_MIXT_ZMOM,ic)

                mol  = pGv(GV_MIXT_MOL,indMol*ic)
                rgas = MixtPerf_R_M(mol)

                pCv(CV_MIXT_TEMP,ic) = MixtPerf_T_DPR(r,p,rgas)
              END DO ! ic
                                
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pRegion%mixtInput
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)   
      END SELECT ! cvStateFuture

! =============================================================================
!   Error - invalid input
! =============================================================================

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixt%cvStates

! *****************************************************************************
!   End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE TURB_FluCv2Prim


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: TURB_fluCv2Prim.F90,v $
!   Revision 1.4  2008/12/06 08:44:44  mtcampbe
!   Updated license.
!
!   Revision 1.3  2008/11/19 22:17:56  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.2  2005/10/31 21:09:37  haselbac
!   Changed specModel and SPEC_MODEL_NONE
!
!   Revision 1.1  2004/03/29 21:10:01  wasistho
!   add flu routines
!
!
! ******************************************************************************







