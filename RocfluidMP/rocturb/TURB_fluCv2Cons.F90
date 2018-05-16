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
! Purpose: Convert primitive state vector to consverved variables.
!
! Description: None.
!
! Input: region        = Pointer to data of current region
!        cvStateFuture = Future state of conserved variables
!
! Output: None.
!
! Notes: 
!   1. Strictly speaking, cvStateFuture is not needed (there is only one
!      state for conserved variables), but kept for consistency with 
!      TURB_FluCv2Prim.
!
!******************************************************************************
  
SUBROUTINE TURB_FluCv2Cons(pRegion,cvStateFuture)

  USE ModDataTypes
  USE ModGlobal, ONLY    : t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY      : t_grid
  USE ModInterfaces, ONLY: MixtPerf_Eo_DGPUVW, MixtPerf_G_CpR, MixtPerf_R_M
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
  INTEGER :: ic,indCp,indMol
  REAL(RFREAL) :: cp,g,mol,p,r,rgas,u,v,w
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv
  TYPE(t_grid),   POINTER :: pGrid
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: TURB_fluCv2Cons.F90,v $ $Revision: 1.4 $'

  global => pRegion%global

  CALL RegisterFunction(global,'TURB_FluCv2Cons',&
  'TURB_fluCv2Cons.F90')

! *****************************************************************************
! Set pointers and variables
! *****************************************************************************

  pGrid => pRegion%grid
  pCv   => pRegion%mixt%cv
  pDv   => pRegion%mixt%dv
  pGv   => pRegion%mixt%gv

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol

! *****************************************************************************
! Actual conversion
! *****************************************************************************
     
  IF ( pRegion%mixt%cvState == CV_MIXT_STATE_DUVWP .OR. & 
       pRegion%mixt%cvState == CV_MIXT_STATE_DUVWT ) THEN 
    
! =============================================================================
!   Convert from primitive to conservative form 
! =============================================================================

    SELECT CASE (cvStateFuture)
      CASE (CV_MIXT_STATE_CONS)
        pRegion%mixt%cvState = CV_MIXT_STATE_CONS          

        SELECT CASE (pRegion%mixtInput%gasModel)

! ------- Perfect gas

          CASE (GAS_MODEL_TCPERF)
            DO ic = 1,pGrid%nCellsTot
              r = pCv(CV_MIXT_DENS,ic)
              u = pCv(CV_MIXT_XVEL,ic)
              v = pCv(CV_MIXT_YVEL,ic)
              w = pCv(CV_MIXT_ZVEL,ic)
              p = pDv(DV_MIXT_PRES,ic)

              pCv(CV_MIXT_XMOM,ic) = r*u
              pCv(CV_MIXT_YMOM,ic) = r*v
              pCv(CV_MIXT_ZMOM,ic) = r*w

              cp   = pGv(GV_MIXT_CP,indCp*ic)
              mol  = pGv(GV_MIXT_MOL,indMol*ic)                  
              rgas = MixtPerf_R_M(mol)
              g    = MixtPerf_G_CpR(cp,rgas)

              pCv(CV_MIXT_ENER,ic) = r*MixtPerf_Eo_DGPUVW(r,g,p,u,v,w)
            END DO ! ic 

! ------- Other or invalid gas models                

          CASE DEFAULT 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! pRegion%mixtInput

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! cvStateFuture

! =============================================================================
! Error - invalid input
! =============================================================================

  ELSE
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END IF ! pRegion%mixt%cvState

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE TURB_FluCv2Cons


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: TURB_fluCv2Cons.F90,v $
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







