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
!   1. bcFlag is initialized in RFLU_InitGENXInterface, because Jim needs that
!      value when COM_init_dataitem is called.
!
! ******************************************************************************
!
! $Id: RFLU_PutBoundaryValues.F90,v 1.21 2008/12/06 08:44:01 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PutBoundaryValues(region)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModParameters
  USE ModMPI  

  USE RFLU_ModGridSpeedUtils, ONLY: RFLU_DescaleGridSpeed
  USE RFLU_ModRindStates, ONLY: RFLU_SetRindStateInjectPerf, & 
                                RFLU_SetRindStateSlipWallPerf
  
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

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: c1,distrib,gasModel,ifl,indCp,indGs,indMol,iPatch,iv,ivg
  REAL(RFREAL) :: chRef,cp,dummyReal,fs,fsu,minj,mm,nx,ny,nz,pl,rl,rRef,rul, &
                  rvl,rwl,tinj,tl,vRef
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pDv,pGv,pMixtVals
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PutBoundaryValues.F90,v $ $Revision: 1.21 $'

  global => region%global

  CALL RegisterFunction(global,'RFLU_PutBoundaryValues',&
  'RFLU_PutBoundaryValues.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pRegion => region
  pGrid   => pRegion%grid
  pCv     => pRegion%mixt%cv
  pDv     => pRegion%mixt%dv
  pGv     => pRegion%mixt%gv

  indGs    = region%grid%indGs
  indCp    = region%mixtInput%indCp
  indMol   = region%mixtInput%indMol
  gasModel = region%mixtInput%gasModel

  rRef = global%refDensity
  vRef = global%refVelocity
  
  chRef = 0.5_RFREAL*(rRef*vRef*vRef*vRef)

! ******************************************************************************
! Loop over ALL boundaries 
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A,1X,I3)') SOLVER_NAME, & 
                                   'Putting values into GENX buffers...'
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

! ==============================================================================
!   Surface grid
! ==============================================================================

    DO iv = 1,pPatch%nBVert
      ivg = pPatch%bv(iv)

      pPatch%xyz(XCOORD,iv) = pGrid%xyz(XCOORD,ivg)
      pPatch%xyz(YCOORD,iv) = pGrid%xyz(YCOORD,ivg)
      pPatch%xyz(ZCOORD,iv) = pGrid%xyz(ZCOORD,ivg)        
    END DO ! iv

    DO iv = pPatch%nBVert+1,pPatch%nBVertTot
      pPatch%xyz(XCOORD,iv) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pPatch%xyz(YCOORD,iv) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      pPatch%xyz(ZCOORD,iv) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    END DO ! iv    

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%checkLevel == CHECK_HIGH .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      IF ( pPatch%nBVert > 0 ) THEN 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Minimum/maximum values:'    
        WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'xyz.x:  ', & 
              MINVAL(pPatch%xyz(XCOORD,1:pPatch%nBVert)), & 
              MAXVAL(pPatch%xyz(XCOORD,1:pPatch%nBVert))
        WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'xyz.y:  ', & 
              MINVAL(pPatch%xyz(YCOORD,1:pPatch%nBVert)), & 
              MAXVAL(pPatch%xyz(YCOORD,1:pPatch%nBVert))
        WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'xyz.z:  ', & 
              MINVAL(pPatch%xyz(ZCOORD,1:pPatch%nBVert)), & 
              MAXVAL(pPatch%xyz(ZCOORD,1:pPatch%nBVert))
      END IF ! pPatch%nBVert
    END IF ! global%myProcid

! ==============================================================================
!   Outgoing data 
! ==============================================================================
    
! ------------------------------------------------------------------------------
!   Interacting non-burning boundary
! ------------------------------------------------------------------------------

    IF ( pPatch%bcCoupled == BC_NOT_BURNING ) THEN
      DO ifl = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifl)

        pPatch%rhofAlp(ifl)      = pCv(CV_MIXT_DENS,c1)
        pPatch%pf(ifl)           = pDv(DV_MIXT_PRES,c1)  

        pPatch%qc(ifl)           = chRef*pPatch%ch(ifl)
        pPatch%qr(ifl)           = 0.0_RFREAL ! TEMPORARY
        
        pPatch%nfAlp(XCOORD,ifl) = pPatch%fn(XCOORD,ifl)
        pPatch%nfAlp(YCOORD,ifl) = pPatch%fn(YCOORD,ifl)
        pPatch%nfAlp(ZCOORD,ifl) = pPatch%fn(ZCOORD,ifl)
    
        pPatch%tracf(XCOORD,ifl) = pPatch%pf(ifl)*pPatch%nfAlp(XCOORD,ifl)   
        pPatch%tracf(YCOORD,ifl) = pPatch%pf(ifl)*pPatch%nfAlp(YCOORD,ifl)
        pPatch%tracf(ZCOORD,ifl) = pPatch%pf(ifl)*pPatch%nfAlp(ZCOORD,ifl)     
      END DO ! ifl     
      
      DO ifl = pPatch%nBFaces+1,pPatch%nBFacesTot
        pPatch%rhofAlp(ifl)      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%pf(ifl)           = REAL(CRAZY_VALUE_INT,KIND=RFREAL)

        pPatch%qc(ifl)           = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%qr(ifl)           = REAL(CRAZY_VALUE_INT,KIND=RFREAL)

        pPatch%nfAlp(XCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%nfAlp(YCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%nfAlp(ZCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    
        pPatch%tracf(XCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%tracf(YCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%tracf(ZCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      END DO ! ifl

      IF ( global%myProcid == MASTERPROC .AND. & 
           global%checkLevel == CHECK_HIGH .AND. & 
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        IF ( pPatch%nBFaces > 0 ) THEN
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'rhofAlp:', & 
                MINVAL(pPatch%rhofAlp(1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%rhofAlp(1:pPatch%nBFaces))       
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'pf:     ', & 
                MINVAL(pPatch%pf(1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%pf(1:pPatch%nBFaces)) 
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'qc:     ', & 
                MINVAL(pPatch%qc(1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%qc(1:pPatch%nBFaces))
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'qr:     ', & 
                MINVAL(pPatch%qr(1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%qr(1:pPatch%nBFaces))                 
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'nfAlp.x:', & 
                MINVAL(pPatch%nfAlp(XCOORD,1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%nfAlp(XCOORD,1:pPatch%nBFaces)) 
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'nfAlp.y:', & 
                MINVAL(pPatch%nfAlp(YCOORD,1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%nfAlp(YCOORD,1:pPatch%nBFaces))  
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'nfAlp.z:', & 
                MINVAL(pPatch%nfAlp(ZCOORD,1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%nfAlp(ZCOORD,1:pPatch%nBFaces))  
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'tracf.x:', & 
                MINVAL(pPatch%tracf(XCOORD,1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%tracf(XCOORD,1:pPatch%nBFaces)) 
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'tracf.y:', & 
                MINVAL(pPatch%tracf(YCOORD,1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%tracf(YCOORD,1:pPatch%nBFaces))  
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'tracf.z:', & 
                MINVAL(pPatch%tracf(ZCOORD,1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%tracf(ZCOORD,1:pPatch%nBFaces)) 
        END IF ! pPatch%nBFaces
      END IF ! global%myProcid      
    END IF ! pPatch%bcCoupled

! ------------------------------------------------------------------------------
!   Interacting burning boundary
! ------------------------------------------------------------------------------

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
					     
            tl = tinj					     
          ELSE ! Surface NOT burning
            rl  = pCv(CV_MIXT_DENS,c1)	
            rul = pCv(CV_MIXT_XMOM,c1)
            rvl = pCv(CV_MIXT_YMOM,c1)
            rwl = pCv(CV_MIXT_ZMOM,c1)	  
            pl  = pDv(DV_MIXT_PRES,c1)
		
            CALL RFLU_SetRindStateSlipWallPerf(cp,mm,nx,ny,nz,rl,rul,rvl,rwl, &
                                               fsu,pl)
					     
            tl = pDv(DV_MIXT_TEMP,c1)					                     
          END IF ! minj
        ELSE 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! gasModel
    
        pPatch%rhofAlp(ifl)      = rl
        pPatch%pf(ifl)           = pl
        pPatch%tempf(ifl)        = tl  

        pPatch%qc(ifl)           = chRef*pPatch%ch(ifl)
        pPatch%qr(ifl)           = 0.0_RFREAL ! TEMPORARY
        
        pPatch%nfAlp(XCOORD,ifl) = nx
        pPatch%nfAlp(YCOORD,ifl) = ny
        pPatch%nfAlp(ZCOORD,ifl) = nz
    
        pPatch%tracf(XCOORD,ifl) = pPatch%pf(ifl)*pPatch%nfAlp(XCOORD,ifl)   
        pPatch%tracf(YCOORD,ifl) = pPatch%pf(ifl)*pPatch%nfAlp(YCOORD,ifl)
        pPatch%tracf(ZCOORD,ifl) = pPatch%pf(ifl)*pPatch%nfAlp(ZCOORD,ifl)    
      END DO ! ifl

      DO ifl = pPatch%nBFaces+1,pPatch%nBFacesTot
        pPatch%rhofAlp(ifl)      = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%pf(ifl)           = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%tempf(ifl)        = REAL(CRAZY_VALUE_INT,KIND=RFREAL)

        pPatch%qc(ifl)           = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%qr(ifl)           = REAL(CRAZY_VALUE_INT,KIND=RFREAL)

        pPatch%nfAlp(XCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%nfAlp(YCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%nfAlp(ZCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    
        pPatch%tracf(XCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%tracf(YCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
        pPatch%tracf(ZCOORD,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
      END DO ! ifl

      IF ( global%myProcid == MASTERPROC .AND. & 
           global%checkLevel == CHECK_HIGH .AND. & 
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        IF ( pPatch%nBFaces > 0 ) THEN
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'rhofAlp:', & 
                MINVAL(pPatch%rhofAlp(1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%rhofAlp(1:pPatch%nBFaces))     
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'pf:     ', & 
                MINVAL(pPatch%pf(1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%pf(1:pPatch%nBFaces))  
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'tempf:  ', & 
                MINVAL(pPatch%tempf(1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%tempf(1:pPatch%nBFaces))
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'qc:     ', & 
                MINVAL(pPatch%qc(1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%qc(1:pPatch%nBFaces))
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'qr:     ', & 
                MINVAL(pPatch%qr(1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%qr(1:pPatch%nBFaces))                 
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'nfAlp.x:', & 
                MINVAL(pPatch%nfAlp(XCOORD,1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%nfAlp(XCOORD,1:pPatch%nBFaces)) 
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'nfAlp.y:', & 
                MINVAL(pPatch%nfAlp(YCOORD,1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%nfAlp(YCOORD,1:pPatch%nBFaces))  
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'nfAlp.z:', & 
                MINVAL(pPatch%nfAlp(ZCOORD,1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%nfAlp(ZCOORD,1:pPatch%nBFaces))  
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'tracf.x:', & 
                MINVAL(pPatch%tracf(XCOORD,1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%tracf(XCOORD,1:pPatch%nBFaces)) 
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'tracf.y:', & 
                MINVAL(pPatch%tracf(YCOORD,1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%tracf(YCOORD,1:pPatch%nBFaces))  
          WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'tracf.z:', & 
                MINVAL(pPatch%tracf(ZCOORD,1:pPatch%nBFaces)), & 
                MAXVAL(pPatch%tracf(ZCOORD,1:pPatch%nBFaces))				
        END IF ! pPatch%nBFaces
      END IF ! global%myProcid
    END IF ! pPatch%bcType
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A,1X,I3)') SOLVER_NAME, &
                                   'Putting values into GENX buffers done.'
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PutBoundaryValues

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PutBoundaryValues.F90,v $
! Revision 1.21  2008/12/06 08:44:01  mtcampbe
! Updated license.
!
! Revision 1.20  2008/11/19 22:17:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.19  2006/08/19 15:37:54  mparmar
! Renamed patch variables
!
! Revision 1.18  2005/10/31 21:09:33  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.17  2005/10/14 14:02:20  haselbac
! Added setting of heat fluxes
!
! Revision 1.16  2005/09/27 14:00:10  haselbac
! Added setting of qc
!
! Revision 1.15  2004/10/19 19:23:09  haselbac
! Removed printing of bFlag extrema and use of relative velocity, cosmetics
!
! Revision 1.14  2004/08/13 22:35:03  haselbac
! Bug fix: pMixtVals only set for burning patches
!
! Revision 1.13  2003/09/04 20:19:36  haselbac
! Removed temporary fix for Rocburn problem in GENx
!
! Revision 1.12  2003/08/22 20:31:55  haselbac
! Added temporary fix for ignition problems
!
! Revision 1.11  2003/08/20 02:09:58  haselbac
! Changed verbosity conditions to reduce solver output in GENx runs
!
! Revision 1.10  2003/05/13 23:44:05  haselbac
! Added init of dummy faces and vertices
!
! Revision 1.9  2003/04/24 15:36:57  haselbac
! Deleted setting of bFlag and logical input argument
!
! Revision 1.8  2003/04/12 21:34:12  haselbac
! Made parallel: Added MASTERPROC checks
!
! Revision 1.7  2003/01/28 18:27:53  haselbac
! Cosmetics only
!
! Revision 1.6  2002/11/27 20:23:30  haselbac
! Fixed bug in setting of variables on burning boundaries
!
! Revision 1.5  2002/10/19 16:11:50  haselbac
! rhof set of interacting patches, cosmetic changes to output
!
! Revision 1.4  2002/10/12 19:08:33  haselbac
! Cosmetic changes only
!
! Revision 1.3  2002/10/12 14:48:51  haselbac
! Corrected bcFlag bug, some cosmetic changes
!
! Revision 1.2  2002/10/07 14:09:32  haselbac
! Changed ModBndpatch to ModBndPatch
!
! Revision 1.1  2002/10/05 18:28:18  haselbac
! Initial revision
!
! ******************************************************************************







