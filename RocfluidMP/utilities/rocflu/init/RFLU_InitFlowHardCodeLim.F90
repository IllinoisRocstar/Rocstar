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
! Purpose: Initialize part of flow field in a region using hard-coded values.
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
! $Id: RFLU_InitFlowHardCodeLim.F90,v 1.4 2008/12/06 08:44:56 mtcampbe Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_InitFlowHardCodeLim(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModParameters

  USE ModInterfaces, ONLY: MixtPerf_Eo_DGPUVW, &
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
  INTEGER :: icg,indCp,indMol
  REAL(RFREAL) :: cp,d,mw,p,g,gc,u,v,w,x
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pGv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt_input), POINTER :: pMixtInput

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InitFlowHardCodeLim.F90,v $ $Revision: 1.4 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_InitFlowHardCodeLim', &
                        'RFLU_InitFlowHardCodeLim.F90')
 
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Initializing flow field from limited hard code...'
                             
    IF ( global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,A)') SOLVER_NAME,'Case: ',TRIM(global%casename)
    END IF ! global%verbLevel                         
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid      => pRegion%grid
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
      
! TEMPORARY, to be replaced by proper code
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
! END TEMPORARY

! ==============================================================================
!   Compressible fluid model
! ==============================================================================  
    
    CASE ( FLUID_MODEL_COMP )
      pRegion%mixt%cvState = CV_MIXT_STATE_CONS
    
      SELECT CASE ( global%casename )

! ------------------------------------------------------------------------------
!       Cylinder diffracting shock
! ------------------------------------------------------------------------------

        CASE ( "cylds" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN 
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3 
              v = 0.0_RFREAL 
              w = 0.0_RFREAL 
              p = pMixtInput%prepRealVal4 
            
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)
        
              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)            
            
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
            END IF ! x
          END DO ! icg           
  
! ------------------------------------------------------------------------------
!       Skews diffracting shock
! ------------------------------------------------------------------------------

! ----- Shock Mach number of 3.0 -----------------------------------------------  

        CASE ( "skews_ms2p0","skews_ms3p0","skews_ms4p0" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN 
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3 
              v = 0.0_RFREAL 
              w = 0.0_RFREAL 
              p = pMixtInput%prepRealVal4 

              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)
        
              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)  

              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
            END IF ! x
          END DO ! icg           
  
! ------------------------------------------------------------------------------
!       Sphere diffracting shock
! ------------------------------------------------------------------------------

        CASE ( "sphds" )
          DO icg = 1,pGrid%nCellsTot
            x = pGrid%cofg(XCOORD,icg)

            IF ( x < pMixtInput%prepRealVal1 ) THEN 
              d = pMixtInput%prepRealVal2
              u = pMixtInput%prepRealVal3 
              v = 0.0_RFREAL 
              w = 0.0_RFREAL 
              p = pMixtInput%prepRealVal4 
            
              mw = pGv(GV_MIXT_MOL,indMol*icg)
              cp = pGv(GV_MIXT_CP ,indCp *icg)
        
              gc = MixtPerf_R_M(mw)
              g  = MixtPerf_G_CpR(cp,gc)            
            
              pCv(CV_MIXT_DENS,icg) = d
              pCv(CV_MIXT_XMOM,icg) = d*u
              pCv(CV_MIXT_YMOM,icg) = d*v
              pCv(CV_MIXT_ZMOM,icg) = d*w
              pCv(CV_MIXT_ENER,icg) = d*MixtPerf_Eo_DGPUVW(d,g,p,u,v,w)               
            END IF ! x
          END DO ! icg           
  
! ------------------------------------------------------------------------------
!       Default - must be due to input error
! ------------------------------------------------------------------------------

        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
      END SELECT ! global%casename
      
! ==============================================================================
!   Default
! ==============================================================================  
    
    CASE DEFAULT 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__) 
  END SELECT ! pMixtInput%fluidModel

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
          'Initializing flow field from limited hard code done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InitFlowHardCodeLim

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InitFlowHardCodeLim.F90,v $
! Revision 1.4  2008/12/06 08:44:56  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:06  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/03/26 20:22:23  haselbac
! Removed error trap for GL model
!
! Revision 1.1  2005/11/10 02:54:20  haselbac
! Renamed to shorten name
!
! Revision 1.6  2005/09/13 21:38:47  haselbac
! Removed hardcoded gamma value
!
! Revision 1.5  2005/07/05 19:49:04  mparmar
! Added init for diffracting shock over sphere
!
! Revision 1.4  2005/06/16 20:57:11  haselbac
! Now use gGas for Skews cases instead of hardcoded 1.4
!
! Revision 1.3  2005/06/16 20:55:54  haselbac
! Added new Skews cases
!
! Revision 1.2  2005/04/22 15:22:34  haselbac
! Changed message printed to screen
!
! Revision 1.1  2005/04/15 16:20:22  haselbac
! Initial revision
!
! ******************************************************************************







