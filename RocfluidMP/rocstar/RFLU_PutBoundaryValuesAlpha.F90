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
! Purpose: Put outgoing data into buffers for GENX.
!
! Description: None.
!
! Input: 
!   region     	Dimensions of patches, types of BCs, flow variables
!
! Output: None.
!
! Notes: 
!   1. The normal depends also on the fractional time alpha, but Rocflu does 
!      not move the mesh inside the Runge-Kutta stages, so copying the normal
!      vector here is not necessary.
!
! ******************************************************************************
!
! $Id: RFLU_PutBoundaryValuesAlpha.F90,v 1.17 2008/12/06 08:44:01 mtcampbe Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PutBoundaryValuesAlpha(region)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModParameters
  USE ModMPI  

  USE RFLU_ModGridSpeedUtils, ONLY: RFLU_DescaleGridSpeed
  USE RFLU_ModRindStates, ONLY: RFLU_SetRindStateInjectPerf
  
  IMPLICIT NONE

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region), TARGET :: region
  LOGICAL :: initialize

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: c1,distrib,gasModel,ifl,indCp,indGs,indMol,iPatch
  REAL(RFREAL) :: cp,dummyReal,fs,fsu,minj,mm,nx,ny,nz,pl,rl,tinj
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pMixtVals
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = & 
   '$RCSfile: RFLU_PutBoundaryValuesAlpha.F90,v $ $Revision: 1.17 $'

  global => region%global

  CALL RegisterFunction(global,'RFLU_PutBoundaryValuesAlpha',&
  'RFLU_PutBoundaryValuesAlpha.F90')

! *****************************************************************************
! Set pointers
! *****************************************************************************

  pRegion => region
  pGrid   => pRegion%grid
  pCv     => pRegion%mixt%cv
  pDv     => pRegion%mixt%dv
  pGv     => pRegion%mixt%gv

  indGs    = region%grid%indGs
  indCp    = region%mixtInput%indCp
  indMol   = region%mixtInput%indMol
  gasModel = region%mixtInput%gasModel

! *****************************************************************************
! Loop over ALL interacting boundaries 
! *****************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A,1X,I3)') SOLVER_NAME,'Putting alpha-dependent '// & 
                                   'values into GENX buffers...'
  END IF ! global%myProcid

  DO iPatch=1,pGrid%nPatches
    pPatch => region%patches(iPatch)

    distrib = pPatch%mixt%distrib

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%checkLevel == CHECK_HIGH .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      IF ( pPatch%bcCoupled == BC_NOT_COUPLED ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A)') SOLVER_NAME,'Patch:',iPatch, & 
	                                    '(not interacting)'
      ELSE IF ( pPatch%bcCoupled == BC_NOT_BURNING ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A)') SOLVER_NAME,'Patch:',iPatch, & 
	                                    '(interacting)'      
      ELSE IF ( pPatch%bcCoupled == BC_BURNING ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A)') SOLVER_NAME,'Patch:',iPatch, & 
	                                    '(interacting and burning)'          
      END IF ! pPatch%bcCoupled                 
    END IF ! global%myProcid 

! =============================================================================
!   Outgoing data (initialized in RFLU_InitGENXInterface)
! =============================================================================

! -----------------------------------------------------------------------------
!   Interacting non-burning boundary
! -----------------------------------------------------------------------------

    IF ( pPatch%bcCoupled == BC_NOT_BURNING ) THEN
      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        pPatch%rhofAlp(ifl) = pCv(CV_MIXT_DENS,c1)  
        
        pPatch%nfAlp(XCOORD,ifl) = pPatch%fn(XCOORD,ifl)
        pPatch%nfAlp(YCOORD,ifl) = pPatch%fn(YCOORD,ifl)
        pPatch%nfAlp(ZCOORD,ifl) = pPatch%fn(ZCOORD,ifl)          
      END DO ! ifl 

      DO ifl = pPatch%nBFaces+1,pPatch%nBFacesTot
        pPatch%rhofAlp(ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        
        pPatch%nfAlp(XCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%nfAlp(YCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%nfAlp(ZCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      END DO ! ifl 

      IF ( global%myProcid == MASTERPROC .AND. & 
           global%checkLevel == CHECK_HIGH .AND. & 
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        IF ( pPatch%nBFaces > 0 .AND. global%verbLevel >= VERBOSE_HIGH) THEN
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Minimum/maximum values:'     
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'rhofAlp: ', & 
                MINVAL(pPatch%rhofAlp(1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%rhofAlp(1:pPatch%nBFaces))
        END IF ! pPatch%nBFaces
      END IF ! global%myProcid
    END IF ! pPatch%bcCoupled

! -----------------------------------------------------------------------------
!   Interacting burning boundary - NOTE call to BcondInjectionPerf only done to
!   set rhofAlp, other quantities are not used or set, hence set rhoVrel to
!   crazy values and use dummyRealOut
! -----------------------------------------------------------------------------

    IF ( pPatch%bcCoupled == BC_BURNING ) THEN
      pMixtVals => pPatch%mixt%vals

      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        nx = pPatch%fn(XCOORD,ifl)
        ny = pPatch%fn(YCOORD,ifl)
        nz = pPatch%fn(ZCOORD,ifl)        
        
        fs  = pPatch%gs(indGs*ifl)
        fsu = RFLU_DescaleGridSpeed(pRegion,fs)
		
        IF ( gasModel == GAS_MODEL_TCPERF ) THEN
          cp = pGv(GV_MIXT_CP ,indCp *c1)
          mm = pGv(GV_MIXT_MOL,indMol*c1)

          minj = pMixtVals(BCDAT_INJECT_MFRATE,distrib*ifl)
	
          IF ( minj > 0.0_RFREAL ) THEN ! Surface burning
            tinj = pMixtVals(BCDAT_INJECT_TEMP,distrib*ifl)

            pl = pDv(DV_MIXT_PRES,c1)

            CALL RFLU_SetRindStateInjectPerf(cp,mm,nx,ny,nz,minj,tinj,pl, &
                                             fsu,rl,dummyReal,dummyReal, &
		  			     dummyReal,dummyReal)					     				     
          ELSE ! Surface NOT burning
            rl = pCv(CV_MIXT_DENS,c1)					                     
          END IF ! minj
        ELSE 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! gasModel

        pPatch%rhofAlp(ifl) = rl
        
        pPatch%nfAlp(XCOORD,ifl) = nx
        pPatch%nfAlp(YCOORD,ifl) = ny
        pPatch%nfAlp(ZCOORD,ifl) = nz         
      END DO ! ifl 

      DO ifl = pPatch%nBFaces+1,pPatch%nBFacesTot
        pPatch%rhofAlp(ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        
        pPatch%nfAlp(XCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%nfAlp(YCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%nfAlp(ZCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      END DO ! ifl

      IF ( global%myProcid == MASTERPROC .AND. & 
           global%checkLevel == CHECK_HIGH .AND. & 
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        IF ( pPatch%nBFaces > 0 .AND. global%verbLevel >= VERBOSE_HIGH) THEN
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Minimum/maximum values:'     
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'rhofAlp: ', & 
                MINVAL(pPatch%rhofAlp(1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%rhofAlp(1:pPatch%nBFaces))				
        END IF ! pPatch%nBFaces
      END IF ! global%myProcid     
    END IF ! pPatch%bcCoupled    
  END DO ! iPatch

! *****************************************************************************
! End
! *****************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A,1X,I3)') SOLVER_NAME,'Putting alpha-dependent '// & 
                                   'values into GENX buffers done.'
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PutBoundaryValuesAlpha

! *****************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PutBoundaryValuesAlpha.F90,v $
! Revision 1.17  2008/12/06 08:44:01  mtcampbe
! Updated license.
!
! Revision 1.16  2008/11/19 22:17:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.15  2006/08/19 15:37:56  mparmar
! Renamed patch variables
!
! Revision 1.14  2005/10/31 21:09:33  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.13  2004/10/19 19:35:59  haselbac
! Changed setting of quantities on burning boundaries, cosmetics
!
! Revision 1.12  2004/08/13 22:35:03  haselbac
! Bug fix: pMixtVals only set for burning patches
!
! Revision 1.11  2003/09/04 20:19:36  haselbac
! Removed temporary fix for Rocburn problem in GENx
!
! Revision 1.10  2003/08/23 21:41:31  haselbac
! Fixed comment
!
! Revision 1.9  2003/08/22 20:31:55  haselbac
! Added temporary fix for ignition problems
!
! Revision 1.8  2003/08/20 02:09:58  haselbac
! Changed verbosity conditions to reduce solver output in GENx runs
!
! Revision 1.7  2003/05/13 23:44:30  haselbac
! Added init of dummy faces
!
! Revision 1.6  2003/04/12 21:34:12  haselbac
! Made parallel: Added MASTERPROC checks
!
! Revision 1.5  2003/03/05 20:37:16  jiao
! ACH: Changed call to BcondInjectionPerf, see comment
!
! Revision 1.4  2003/02/24 23:42:40  haselbac
! Added setting of nfAlp at request of Jim for PC iterations
!
! Revision 1.3  2003/01/28 18:27:30  haselbac
! Cosmetics only
!
! Revision 1.2  2002/11/27 20:23:30  haselbac
! Fixed bug in setting of variables on burning boundaries
!
! Revision 1.1  2002/10/12 19:06:13  haselbac
! Initial revision
!
! *****************************************************************************







