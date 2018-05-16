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
! Purpose: Check for posivity of variables
!
! Description: None.
!
! Input: 
!   pRegion     Region data
!
! Output: None.
!
! Notes: 
!   1. Compute and check pressure here because, it being computed later 
!      through a call to MixtureProperties, a negative value would not be
!      detected and could cause negative pressure and temperature.
!
! ******************************************************************************
!
! $Id: RFLU_CheckPositivity_GL.F90,v 1.4 2008/12/06 08:44:11 mtcampbe Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CheckPositivity_GL(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModParameters
  USE ModMPI
  
  USE ModInterfaces, ONLY: MixtGasLiq_P, &
                           MixtLiq_C2_Bp, &
                           MixtPerf_C2_GRT, &
                           MixtPerf_Cv_CpR, &
                           MixtPerf_R_M, &
                           MixtPerf_T_CvEoVm2, &
                           RFLU_PrintLocInfo
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER, PARAMETER :: MAX_NEGATIVE_LOCS = 10
  INTEGER :: icg,indCp,indMol,nLocs
  INTEGER :: loc(MAX_NEGATIVE_LOCS,MIN_VAL:MAX_VAL)
  REAL(RFREAL) :: Bp,Bt,Cg2,cpGas,cpVap,Cl2,Cv2,cvg,cvl,Cvm,cvv,Eo,p,po, &
                  rGas,rho,rhoYg,rhoYl,rhoYv,ro,rrho,rVap,t,to,u,v,Vm2,w
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pCvSpec,pDv,pGv
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_CheckPositivity_GL.F90,v $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CheckPositivity_GL',&
  'RFLU_CheckPositivity_GL.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

  pCv => pRegion%mixt%cv
  pDv => pRegion%mixt%dv
  pGv => pRegion%mixt%gv

#ifdef SPEC
  pCvSpec => pRegion%spec%cv
#endif

  nLocs = 0

  indCp  = pRegion%mixtInput%indCp
  indMol = pRegion%mixtInput%indMol 

  ro  = global%refDensityLiq
  po  = global%refPressLiq
  to  = global%refTempLiq
  Bp  = global%refBetaPLiq
  Bt  = global%refBetaTLiq
  cvl = global%refCvLiq

  rGas = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
  cvg  = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht,rGas)

  rVap = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
  cvv  = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht,rVap)

! ******************************************************************************
! Loop over cells and check for positivity
! ******************************************************************************

  DO icg = 1,pGrid%nCells
    rho   = pCv(CV_MIXT_DENS,icg)
    rrho  = 1.0_RFREAL/rho
    rhoYg = pCvSpec(1,icg)
    rhoYv = pCvSpec(2,icg)
    rhoYl = rho - rhoYg - rhoYv

    u  = rrho*pCv(CV_MIXT_XMOM,icg)
    v  = rrho*pCv(CV_MIXT_YMOM,icg)
    w  = rrho*pCv(CV_MIXT_ZMOM,icg)        
    Eo = rrho*pCv(CV_MIXT_ENER,icg)

    Vm2 = u*u + v*v + w*w
    Cvm = (rhoYl*cvl + rhoYv*cvv + rhoYg*cvg)/rho
    t   = MixtPerf_T_CvEoVm2(Cvm,Eo,Vm2)

    Cl2 = MixtLiq_C2_Bp(Bp) 
    Cv2 = MixtPerf_C2_GRT(1.0_RFREAL,rGas,t)
    Cg2 = MixtPerf_C2_GRT(1.0_RFREAL,rVap,t)
 
    p = MixtGasLiq_P(rhoYl,rhoYv,rhoYg,Cl2,Cv2,Cg2,rho,ro,po,to,Bp,Bt,t)

    IF ( (rho <= 0.0_RFREAL) .OR. (p <= 0.0_RFREAL) ) THEN
        nLocs = nLocs + 1   

      IF ( nLocs == 1 ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
              'Negative positive-definite variables detected!'
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Module: Mixture.'        
              
        IF ( global%flowType == FLOW_UNSTEADY ) THEN 
          WRITE(STDOUT,'(A,3X,A,1X,1PE12.5)') SOLVER_NAME,'Current time:', &
                                              global%currentTime              
        ELSE 
          WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME, &
                                         'Current iteration number:', &
                                         global%currentIter           
        END IF ! global%flowType                 
                                            
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal 
        WRITE(STDOUT,'(A,6X,A,6(1X,A))') SOLVER_NAME,'#', &
                                         '   Density   ', &
                                         '  x-velocity ', &
                                         '  y-velocity ', &
                                         '  z-velocity ', &
                                         '   Pressure  ', &
                                         ' Temperature '       
      END IF ! nLocs

      IF ( nLocs <= MAX_NEGATIVE_LOCS ) THEN 
        WRITE(STDOUT,'(A,4X,I3,6(1X,E13.6))') SOLVER_NAME,nLocs, & 
                                              rho,u,v,w,p,t
        loc(nLocs,MIN_VAL:MAX_VAL) = icg                                   
      END IF ! nLocs
    END IF ! pCv    
  END DO ! icg

! ******************************************************************************
! Write out message and call error handling routine
! ******************************************************************************

  IF ( nLocs > 0 ) THEN 
    IF ( nLocs > MAX_NEGATIVE_LOCS ) THEN 
       WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A,1X,I9,1X,A)') SOLVER_NAME, & 
             'Only wrote the first',MAX_NEGATIVE_LOCS,'of',nLocs, & 
             'cells with negative positive-definite variables.'    
      CALL RFLU_PrintLocInfo(pRegion,loc,MAX_NEGATIVE_LOCS, & 
                             LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
    ELSE 
      CALL RFLU_PrintLocInfo(pRegion,loc(1:nLocs,MIN_VAL:MAX_VAL),nLocs, & 
                             LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
    END IF ! nLocs
    
    CALL ErrorStop(global,ERR_NEGATIVE_POSDEF,__LINE__)   
  END IF ! nLocs

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CheckPositivity_GL

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CheckPositivity_GL.F90,v $
! Revision 1.4  2008/12/06 08:44:11  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.1  2006/03/26 20:20:56  haselbac
! Initial revision
!
! ******************************************************************************







