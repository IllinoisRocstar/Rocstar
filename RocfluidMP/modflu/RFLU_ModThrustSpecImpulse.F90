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
! Purpose: Collection of routines to compute thrust and specific impulse.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModThrustSpecImpulse.F90,v 1.4 2008/12/06 08:44:24 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModThrustSpecImpulse

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$$'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_TSI_ComputeGlobalThrustSI, &  
            RFLU_TSI_PrintGlobalVals, & 
            RFLU_TSI_WriteGlobalVals
            
! ==============================================================================
! Private functions
! ==============================================================================
  
  PRIVATE :: RFLU_TSI_CloseGlobalVals, &
             RFLU_TSI_OpenGlobalVals

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  





! ******************************************************************************
!
!   Purpose: Close thrust and specific impulse files.
!
!   Description: None.
!
!   Input:
!     global            Pointer to global data
!
!   Notes:
!
! ******************************************************************************

  SUBROUTINE RFLU_TSI_CloseGlobalVals(global,iPatch,iFile1,iFile2)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

! *****************************************************************************
!   Declarations and definitions
! *****************************************************************************

! =============================================================================
!   Arguments
! =============================================================================

    INTEGER, INTENT(IN) :: iPatch,iFile1,iFile2
    TYPE(t_global), POINTER :: global
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName1,iFileName2
    INTEGER :: errorFlag 
   
! *****************************************************************************
!   Start
! *****************************************************************************
 
    CALL RegisterFunction(global,'RFLU_TSI_CloseGlobalVals',&
  'RFLU_ModThrustSpecImpulse.F90')

! ==============================================================================
!   Opening thrust and specific impulse file
! ==============================================================================

    IF ( global%verbLevel > VERBOSE_LOW ) THEN       
      WRITE(STDOUT,'(A,1X,A,I4.4)') SOLVER_NAME,'Closing thrust and '// &
                               'specific impulse files for iPatch =',iPatch 
    END IF ! global%verbLevel
      
! ------------------------------------------------------------------------------
!   Build file name
! ------------------------------------------------------------------------------

    WRITE(iFileName1,'(A,I4.4)') TRIM(global%outDir)// & 
                                 TRIM(global%casename)//'.thr_',iPatch

    WRITE(iFileName2,'(A,I4.4)') TRIM(global%outDir)// & 
                                 TRIM(global%casename)//'.isp_',iPatch

! ------------------------------------------------------------------------------
!   Open file
! ------------------------------------------------------------------------------

    CLOSE(iFile1,IOSTAT=errorFlag)
    IF (global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '// & 
                     TRIM(iFileName1))
    END IF ! global%error
    CLOSE(iFile2,IOSTAT=errorFlag)
    IF (global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '// & 
                     TRIM(iFileName2))
    END IF ! global%error

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_LOW ) THEN       
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing thrust and '// & 
                               'specific impulse files done.'             
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_TSI_CloseGlobalVals










! *******************************************************************************
!
! Purpose: Compute thrust and specific impulse coefficients.
!
! Description: None.
!
! Input:
!   regions             Region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_TSI_ComputeGlobalThrustSI(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,iPatch,iPatchGlobal
    REAL(RFREAL) :: forceCoeff2Thrust,thrust2Isp,massCoeffIn,massCoeffOut, &
                    massIn,massOut
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pForceCoeff,pForceCoeffVac
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pSpecImpulse,pSpecImpulseVac, &
                                             pThrust,pThrustVac
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_TSI_ComputeGlobalThrustSI',&
  'RFLU_ModThrustSpecImpulse.F90')

    pForceCoeff     => pRegion%forceCoeffsGlobal
    pForceCoeffVac  => pRegion%forceVacCoeffsGlobal
    pSpecImpulse    => pRegion%specImpulseGlobal
    pSpecImpulseVac => pRegion%specImpulseVacGlobal
    pThrust         => pRegion%thrustGlobal
    pThrustVac      => pRegion%thrustVacGlobal
        
! ******************************************************************************
!   Compute thrust
! ******************************************************************************
        
! ==============================================================================
!   De-normalize mass coefficients
! ==============================================================================

    massCoeffIn  = 0.0_RFREAL
    massCoeffOut = 0.0_RFREAL
    massIn       = 0.0_RFREAL
    massOut      = 0.0_RFREAL

    DO iPatch = 1,pRegion%global%nPatches
      massCoeffIn  = massCoeffIn  + pRegion%massCoeffsGlobal(MASS_IN ,iPatch)
      massCoeffOut = massCoeffOut + pRegion%massCoeffsGlobal(MASS_OUT,iPatch)
    END DO ! iPatch

    massIn  = massCoeffIn*global%refDensity*global%refVelocity*global%forceRefArea
    massOut = massCoeffOut*global%refDensity*global%refVelocity*global%forceRefArea

! ==============================================================================
!   De-normalize force coefficients to compute thrusts
! ==============================================================================

    forceCoeff2Thrust  = 0.5_RFREAL*global%refDensity*global%refVelocity &
                                   *global%refVelocity*global%forceRefArea

    DO iPatch = 1,pRegion%global%nPatches
      pThrust(XCOORD,iPatch) = forceCoeff2Thrust &
                             *(pForceCoeff(XCOORD,COMP_MOM ,iPatch) &
                             + pForceCoeff(XCOORD,COMP_PRES,iPatch) &
                             + pForceCoeff(XCOORD,COMP_VISC,iPatch))
      pThrust(YCOORD,iPatch) = forceCoeff2Thrust &
                             *(pForceCoeff(YCOORD,COMP_MOM ,iPatch) &
                             + pForceCoeff(YCOORD,COMP_PRES,iPatch) &
                             + pForceCoeff(YCOORD,COMP_VISC,iPatch))
      pThrust(ZCOORD,iPatch) = forceCoeff2Thrust &
                             *(pForceCoeff(ZCOORD,COMP_MOM ,iPatch) &
                             + pForceCoeff(ZCOORD,COMP_PRES,iPatch) &
                             + pForceCoeff(ZCOORD,COMP_VISC,iPatch))
      
      pThrustVac(XCOORD,iPatch) = forceCoeff2Thrust &
                                *(pForceCoeffVac(XCOORD,COMP_MOM ,iPatch) &
                                + pForceCoeffVac(XCOORD,COMP_PRES,iPatch) &
                                + pForceCoeffVac(XCOORD,COMP_VISC,iPatch))
      pThrustVac(YCOORD,iPatch) = forceCoeff2Thrust &
                                *(pForceCoeffVac(YCOORD,COMP_MOM ,iPatch) &
                                + pForceCoeffVac(YCOORD,COMP_PRES,iPatch) &
                                + pForceCoeffVac(YCOORD,COMP_VISC,iPatch))
      pThrustVac(ZCOORD,iPatch) = forceCoeff2Thrust &
                                *(pForceCoeffVac(ZCOORD,COMP_MOM ,iPatch) &
                                + pForceCoeffVac(ZCOORD,COMP_PRES,iPatch) &
                                + pForceCoeffVac(ZCOORD,COMP_VISC,iPatch))
    END DO ! iPatch

! ==============================================================================
!   Normalize thrusts to compute specific impulse
! ==============================================================================

    thrust2Isp  = 1.0_RFREAL/(massIn*global%gravity)

    DO iPatch = 1,pRegion%global%nPatches
      pSpecImpulse(XCOORD,iPatch) = thrust2Isp*pThrust(XCOORD,iPatch)
      pSpecImpulse(YCOORD,iPatch) = thrust2Isp*pThrust(YCOORD,iPatch)
      pSpecImpulse(ZCOORD,iPatch) = thrust2Isp*pThrust(ZCOORD,iPatch)
      
      pSpecImpulseVac(XCOORD,iPatch) = thrust2Isp*pThrustVac(XCOORD,iPatch)
      pSpecImpulseVac(YCOORD,iPatch) = thrust2Isp*pThrustVac(YCOORD,iPatch)
      pSpecImpulseVac(ZCOORD,iPatch) = thrust2Isp*pThrustVac(ZCOORD,iPatch)
    END DO ! iPatch

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_TSI_ComputeGlobalThrustSI






  


! ******************************************************************************
!
!   Purpose: Open thrust and specific impulse files.
!
!   Description: None.
!
!   Input:
!     global            Pointer to global data
!
!   Notes:
!
! ******************************************************************************

  SUBROUTINE RFLU_TSI_OpenGlobalVals(global,iPatch,iFile1,iFile2)

    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    USE ModParameters
    USE ModError
    USE ModMPI

! *****************************************************************************
!   Declarations and definitions
! *****************************************************************************

! =============================================================================
!   Arguments
! =============================================================================

    INTEGER, INTENT(IN) :: iPatch
    TYPE(t_global), POINTER :: global
    
    INTEGER, INTENT(OUT) :: iFile1,iFile2

! ==============================================================================
!   Locals
! ==============================================================================
    
    LOGICAL :: fileExists1,fileExists2
    CHARACTER(CHRLEN) :: iFileName1,iFileName2
    INTEGER :: errorFlag 
   
! *****************************************************************************
!   Start
! *****************************************************************************
 
    CALL RegisterFunction(global,'RFLU_TSI_OpenGlobalVals',&
  'RFLU_ModThrustSpecImpulse.F90')

! ==============================================================================
!   Opening thrust and specific impulse file
! ==============================================================================

    IF ( global%verbLevel > VERBOSE_LOW ) THEN       
      WRITE(STDOUT,'(A,1X,A,I4.4)') SOLVER_NAME,'Opening thrust and '// &
                               'specific impulse files for iPatch =',iPatch 
    END IF ! global%verbLevel
      
! ------------------------------------------------------------------------------
!   Build file name
! ------------------------------------------------------------------------------

    iFile1 = IF_THRUST
    iFile2 = IF_ISP

    WRITE(iFileName1,'(A,I4.4)') TRIM(global%outDir)// & 
                                 TRIM(global%casename)//'.thr_',iPatch

    WRITE(iFileName2,'(A,I4.4)') TRIM(global%outDir)// & 
                                 TRIM(global%casename)//'.isp_',iPatch

! ------------------------------------------------------------------------------
!   Open file
! ------------------------------------------------------------------------------

    IF ( global%restartFromScratch .EQV. .FALSE. ) THEN
      INQUIRE(FILE=iFileName1,EXIST=fileExists1)

      IF ( fileExists1 .EQV. .TRUE. ) THEN 
        OPEN(iFile1,FILE=iFileName1,FORM='FORMATTED',STATUS='OLD', &
             POSITION='APPEND',IOSTAT=errorFlag)
      ELSE 
        OPEN(iFile1,FILE=iFileName1,FORM='FORMATTED',STATUS='NEW', &
             IOSTAT=errorFlag)
      END IF ! fileExists1

      INQUIRE(FILE=iFileName2,EXIST=fileExists2)

      IF ( fileExists2 .EQV. .TRUE. ) THEN 
        OPEN(iFile2,FILE=iFileName2,FORM='FORMATTED',STATUS='OLD', &
             POSITION='APPEND',IOSTAT=errorFlag)
      ELSE 
        OPEN(iFile1,FILE=iFileName2,FORM='FORMATTED',STATUS='NEW', &
             IOSTAT=errorFlag)
      END IF ! fileExists2
    ELSE
      IF ( global%thrustWriteCntr == 1 ) THEN 
        OPEN(iFile1,FILE=iFileName1,FORM='FORMATTED',STATUS='UNKNOWN', &
             IOSTAT=errorFlag)
        OPEN(iFile2,FILE=iFileName2,FORM='FORMATTED',STATUS='UNKNOWN', &
             IOSTAT=errorFlag)
      ELSE 
        OPEN(iFile1,FILE=iFileName1,FORM='FORMATTED',STATUS='OLD', &
             POSITION='APPEND',IOSTAT=errorFlag)          
        OPEN(iFile2,FILE=iFileName2,FORM='FORMATTED',STATUS='OLD', &
             POSITION='APPEND',IOSTAT=errorFlag)          
      END IF ! global%thrustWriteCntr
    END IF ! global      

    global%error = errorFlag

    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '// & 
                     TRIM(iFileName1)//' and '//TRIM(iFileName2))
    END IF ! global%error

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_LOW ) THEN       
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening thrust and '// & 
                               'specific impulse files done.'             
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_TSI_OpenGlobalVals











! *******************************************************************************
!
! Purpose: Print thrust and specific impulse.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine should only be called for the master process. It will not 
!      do anything if called by any other process (safeguard).
!
! ******************************************************************************

  SUBROUTINE RFLU_TSI_PrintGlobalVals(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: iPatch 
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_TSI_PrintGlobalVals',&
  'RFLU_ModThrustSpecImpulse.F90')
                
! ******************************************************************************
!   Print global force coefficients
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing thrust and '// &
                               'specific impulse...' 
      
      DO iPatch = 1,pRegion%global%nPatches
        IF ( pRegion%thrustFlagsGlobal(iPatch) .EQV. .TRUE. ) THEN
          WRITE(STDOUT,'(A,5X,A,1X,I2)') SOLVER_NAME,'Patch:',iPatch
          WRITE(STDOUT,'(A,7X,A,7X,A,6X,A,8X,A,8X,A)') SOLVER_NAME, & 
                'Component','Thrust','Thrust vacuum','Isp','Isp vacuum'
          WRITE(STDOUT,'(A,7X,A,4(1X,E13.6))') SOLVER_NAME,'x-direction:', & 
                pRegion%thrustGlobal(XCOORD,iPatch), & 
                pRegion%thrustVacGlobal(XCOORD,iPatch), &
                pRegion%specImpulseGlobal(XCOORD,iPatch), & 
                pRegion%specImpulseVacGlobal(XCOORD,iPatch)
          WRITE(STDOUT,'(A,7X,A,4(1X,E13.6))') SOLVER_NAME,'y-direction:', & 
                pRegion%thrustGlobal(YCOORD,iPatch), & 
                pRegion%thrustVacGlobal(YCOORD,iPatch), &
                pRegion%specImpulseGlobal(YCOORD,iPatch), & 
                pRegion%specImpulseVacGlobal(YCOORD,iPatch)
          WRITE(STDOUT,'(A,7X,A,4(1X,E13.6))') SOLVER_NAME,'z-direction:', & 
                pRegion%thrustGlobal(ZCOORD,iPatch), & 
                pRegion%thrustVacGlobal(ZCOORD,iPatch), &
                pRegion%specImpulseGlobal(ZCOORD,iPatch), & 
                pRegion%specImpulseVacGlobal(ZCOORD,iPatch)
        END IF ! pRegion%thrustFlagsGlobal
      END DO ! iPatch     
 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing thrust and '// & 
                               'specific impulse done.'             
    END IF ! global%myProcid

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_TSI_PrintGlobalVals








! *******************************************************************************
!
! Purpose: Write thrust and specific impulse to file.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine should only be called for the master process. It will not 
!      do anything if called by any other process (safeguard).
!   2. At present, each patch is written to separate file. This may change in
!      future with definition of superpatches which are treated together. 
!
! ******************************************************************************

  SUBROUTINE RFLU_TSI_WriteGlobalVals(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    LOGICAL :: fileExists1,fileExists2
    CHARACTER(CHRLEN) :: iFileName1,iFileName2
    INTEGER :: errorFlag,iFile1,iFile2,iPatch 
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'RFLU_TSI_WriteGlobalVals',&
  'RFLU_ModThrustSpecImpulse.F90')

! ******************************************************************************
!   Increment counter
! ******************************************************************************

    global%thrustWriteCntr = global%thrustWriteCntr + 1
                
! ******************************************************************************
!   Write global force coefficients
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC ) THEN
      IF ( global%verbLevel > VERBOSE_LOW ) THEN   
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing thrust and '// &
                                 'specific impulse...' 
      END IF ! global%verbLevel

! ==============================================================================
!   Loop over patches - DONT LOOP OVER PATCH, just write for a given patch
! ==============================================================================

      DO iPatch = 1,global%nPatches
        IF ( pRegion%thrustFlagsGlobal(iPatch) .EQV. .TRUE. ) THEN
          IF ( global%verbLevel > VERBOSE_LOW ) THEN   
            WRITE(STDOUT,'(A,3X,A,1X,I2)') SOLVER_NAME,'Patch:',iPatch 
          END IF ! global%verbLevel

! ------------------------------------------------------------------------------
!   Open file
! ------------------------------------------------------------------------------

          CALL RFLU_TSI_OpenGlobalVals(global,iPatch,iFile1,iFile2)

! ------------------------------------------------------------------------------
!   Write to file
! ------------------------------------------------------------------------------

          IF ( global%flowType == FLOW_STEADY ) THEN                               
            WRITE(iFile1,'(I6,6(1X,E13.6))') global%currentIter, & 
                  pRegion%thrustGlobal(XCOORD,iPatch), & 
                  pRegion%thrustVacGlobal(XCOORD,iPatch), &
              
                  pRegion%thrustGlobal(YCOORD,iPatch), & 
                  pRegion%thrustVacGlobal(YCOORD,iPatch), &
  
                  pRegion%thrustGlobal(ZCOORD,iPatch), & 
                  pRegion%thrustVacGlobal(ZCOORD,iPatch)
          ELSE 
            WRITE(iFile1,'(1PE12.5,6(1X,E13.6))') global%currentTime, & 
                  pRegion%thrustGlobal(XCOORD,iPatch), & 
                  pRegion%thrustVacGlobal(XCOORD,iPatch), &
                
                  pRegion%thrustGlobal(YCOORD,iPatch), & 
                  pRegion%thrustVacGlobal(YCOORD,iPatch), &
  
                  pRegion%thrustGlobal(ZCOORD,iPatch), & 
                  pRegion%thrustVacGlobal(ZCOORD,iPatch)
          END IF ! global%flowType                                     
  
          IF ( global%flowType == FLOW_STEADY ) THEN                               
            WRITE(iFile2,'(I6,6(1X,E13.6))') global%currentIter, & 
                  pRegion%specImpulseGlobal(XCOORD,iPatch), & 
                  pRegion%specImpulseVacGlobal(XCOORD,iPatch), &
                
                  pRegion%specImpulseGlobal(YCOORD,iPatch), & 
                  pRegion%specImpulseVacGlobal(YCOORD,iPatch), &
  
                  pRegion%specImpulseGlobal(ZCOORD,iPatch), & 
                  pRegion%specImpulseVacGlobal(ZCOORD,iPatch)
          ELSE 
            WRITE(iFile2,'(1PE12.5,6(1X,E13.6))') global%currentTime, & 
                  pRegion%specImpulseGlobal(XCOORD,iPatch), & 
                  pRegion%specImpulseVacGlobal(XCOORD,iPatch), &
                
                  pRegion%specImpulseGlobal(YCOORD,iPatch), & 
                  pRegion%specImpulseVacGlobal(YCOORD,iPatch), &
  
                  pRegion%specImpulseGlobal(ZCOORD,iPatch), & 
                  pRegion%specImpulseVacGlobal(ZCOORD,iPatch)
          END IF ! global%flowType                                     

! ------------------------------------------------------------------------------
!   Close file
! ------------------------------------------------------------------------------

          CALL RFLU_TSI_CloseGlobalVals(global,iPatch,iFile1,iFile2)
        END IF ! pRegion%thrustFlagsGlobal
      END DO! iPatch        
            
! ******************************************************************************
!   End  
! ******************************************************************************

      IF ( global%verbLevel > VERBOSE_LOW ) THEN       
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing thrust and '// & 
                                 'specific impulse done.'             
      END IF ! global%verbLevel
    END IF ! global%myProcid

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_TSI_WriteGlobalVals






END MODULE RFLU_ModThrustSpecImpulse

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModThrustSpecImpulse.F90,v $
! Revision 1.4  2008/12/06 08:44:24  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:35  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/12/19 19:13:42  rfiedler
! Increased verbosity level required for writing to screen dump to > VERBOSE_LOW.
!
! Revision 1.1  2006/10/20 21:16:32  mparmar
! Initial import
!
!
! ******************************************************************************











