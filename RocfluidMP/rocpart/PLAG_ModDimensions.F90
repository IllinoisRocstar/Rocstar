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
! Purpose: Suite of routines to compute, read, and write dimension files.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: PLAG_ModDimensions.F90,v 1.7 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_ModDimensions

  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global 
  USE ModBndPatch, ONLY: t_patch
  USE ModBorder, ONLY: t_border
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY: t_plag,t_plag_input
  USE ModMPI
  
  USE PLAG_ModParameters
    
  IMPLICIT NONE
     
  PRIVATE
  PUBLIC :: PLAG_CalcNPclsGlobal, &
            PLAG_PrintDimensions, &
            PLAG_PrintNPclsGlobal, & 
            PLAG_RFLU_ReadDimensions, & 
            PLAG_SetDimensions, & 
            PLAG_SetMaxDimensions, & 
            PLAG_RFLU_WriteDimensions

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Private
! ==============================================================================
   
  CHARACTER(CHRLEN), PRIVATE :: RCSIdentString = &
    '$RCSfile: PLAG_ModDimensions.F90,v $ $Revision: 1.7 $' 
    
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS





! ******************************************************************************
!
! Purpose: Calculate the total number of particles for all regions.
!
! Description: None.
!
! Input:
!   regions     Pointer to all regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_CalcNPclsGlobal(regions)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: regions(:)
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iReg,iRegGlobal,nPclsGlobal,nPclsLocal

    TYPE(t_global), POINTER :: global    
    TYPE(t_region), POINTER :: pRegion
    TYPE(t_plag),   POINTER :: pPlag 

! ******************************************************************************
!   Start
! ******************************************************************************

    global => regions(0)%global

    CALL RegisterFunction(global,'PLAG_CalcNPclsGlobal',&
  'PLAG_ModDimensions.F90')

! ******************************************************************************
!   Initialize variables
! ******************************************************************************

    nPclsLocal = 0

    DO iReg = 0,global%nRegionsLocal
      pRegion => regions(iReg)

      pRegion%plag%nPclsGlobal = 0
    END DO ! iReg

! ******************************************************************************
!   Compute sum of particles for all regions on the same processor. NOTE cannot
!   include region 0 in this sum because it would pick up the CRAZY_VALUE_INT
!   value assigned below.
! ******************************************************************************

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)

      nPclsLocal = nPclsLocal + pRegion%plag%nPcls
    END DO ! iReg
    
! ******************************************************************************
!   Determine global sum of particles over all processors
! ******************************************************************************

! ==============================================================================
!   Perform reduction operation. NOTE need to include region index 0 
!    to make sure that this works properly for serial runs.
! ==============================================================================

    CALL MPI_Allreduce(nPclsLocal,nPclsGlobal,1,MPI_INTEGER,MPI_SUM, &
                       global%mpiComm,errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
    END IF ! global%errorFlag
    
! ******************************************************************************
!   Store the global sum of particles. NOTE copy into serial nPcls for master
!   process, but for others set to CRAZY_VALUE_INT.
! ******************************************************************************

    DO iReg = 0,global%nRegionsLocal
      pRegion => regions(iReg)

      pRegion%plag%nPclsGlobal = nPclsGlobal 
    END DO ! iReg

    IF ( global%myProcid == MASTERPROC ) THEN 
      regions(0)%plag%nPcls = regions(0)%plag%nPclsGlobal
    ELSE 
      regions(0)%plag%nPcls = CRAZY_VALUE_INT
    END IF ! global%myProcid

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_CalcNPclsGlobal








! ******************************************************************************
!
! Purpose: Print the number of particles.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_PrintDimensions(pRegion)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion
       
! ==============================================================================
!   Local variables
! ==============================================================================

    TYPE(t_global), POINTER :: global    

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_PrintDimensions',&
  'PLAG_ModDimensions.F90')

! ******************************************************************************
!   Print information
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing number of particles...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
      IF ( global%flowType == FLOW_UNSTEADY ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%flowType

      WRITE(STDOUT,'(A,3X,A,9X,I8)') SOLVER_NAME,'Number of particles:', &
                                     pRegion%plag%nPcls
      WRITE(STDOUT,'(A,3X,A,1X,I8)') SOLVER_NAME,'Maximum number of '// &
                                     'particles:',pRegion%plag%nPclsMax

      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing number of particles done.'
    END IF ! global%verbLevel

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_PrintDimensions








! ******************************************************************************
!
! Purpose: Print the total number of particles in all regions.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_PrintNPclsGlobal(pRegion)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion
       
! ==============================================================================
!   Local variables
! ==============================================================================

    TYPE(t_global), POINTER :: global    

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_PrintNPclsGlobal',&
  'PLAG_ModDimensions.F90')

! ******************************************************************************
!   Print information
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing total number of '// & 
                               'particles...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
      IF ( global%flowType == FLOW_UNSTEADY ) THEN
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', &
                                            global%currentTime
      END IF ! global%flowType

      WRITE(STDOUT,'(A,3X,A,1X,I8)') SOLVER_NAME,'Number of particles:', &
                                     pRegion%plag%nPclsGlobal

      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing total number of '// &
                               'particles done.'
    END IF ! global%verbLevel

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_PrintNPclsGlobal







! ******************************************************************************
!
! Purpose: Read dimensions for Lagrangian particles.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLU_ReadDimensions(pRegion)

    USE RFLU_ModGrid

    USE ModBuildFileNames, ONLY: BuildFileNameBasic, &
                                 BuildFileNameUnsteady  

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion  
   
! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: errorString,iFileName,RCSIdentString,sectionString, & 
                         timeString1,timeString2
    INTEGER :: errorFlag,dummy,iPatch,iFile,loopCounter,nCont,nPclsMax
    REAL(RFREAL) :: currentTime
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_RFLU_ReadDimensions',&
  'PLAG_ModDimensions.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading particle dimensions...'
    END IF ! global%verbLevel

! ==============================================================================
!   Build file name
! ==============================================================================

    iFile = IF_DIMS

    IF ( global%flowType == FLOW_UNSTEADY ) THEN 
#ifndef GENX
      currentTime = global%currentTime
#else
      IF ( global%timeStamp > 0.0_RFREAL ) THEN
        global%warnCounter = global%warnCounter + 1    

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_NONE ) THEN        
          WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'*** WARNING ***', & 
                                        'Hard-code - read file from time zero.'
        END IF ! global%myProcid
      END IF ! global%timeStamp

      currentTime = 0.0_RFREAL ! Hard-code for GENX restart
#endif        

! TEMPORARY
      IF ( pRegion%iRegionGlobal == 0 ) THEN 
        currentTime = 0.0_RFREAL
      END IF ! pRegion%iRegionGlobal
! END TEMPORARY

      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.pdim', & 
                                 pRegion%iRegionGlobal,currentTime, & 
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN        
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal          
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                            global%currentTime   
      END IF ! global%verbLevel
    ELSE 
      CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.pdim', & 
                              pRegion%iRegionGlobal,iFileName) 

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal
      END IF ! global%verbLevel
    END IF ! global

! ==============================================================================
!   Open file
! ==============================================================================

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD", & 
         IOSTAT=errorFlag)   
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error
    
! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel 

    READ(iFile,'(A)') sectionString  
    IF ( TRIM(sectionString) /= '# ROCPART dimensions file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM  
    
! ==============================================================================
!   Rest of file
! ==============================================================================

    pGrid => pRegion%grid  

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile,'(A)') sectionString

      SELECT CASE ( TRIM(sectionString) ) 

! ------------------------------------------------------------------------------
!       Actual number of particles
! ------------------------------------------------------------------------------ 
    
        CASE ( '# Actual number of particles' )     
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Actual number of particles...'
          END IF ! global%verbLevel          

          READ(iFile,'(I8)') pRegion%plag%nPcls   

! ------------------------------------------------------------------------------
!       Maximum number of particles
! ------------------------------------------------------------------------------ 
        
        CASE ( '# Maximum number of particles' )     
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Maximum number of particles...'
          END IF ! global%verbLevel          

          READ(iFile,'(I8)') pRegion%plag%nPclsMax     

! ------------------------------------------------------------------------------
!       Number of constituents
! ------------------------------------------------------------------------------ 
        
        CASE ( '# Number of constituents' )     
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Number of constituents...'
          END IF ! global%verbLevel          

          READ(iFile,'(I8)') nCont

          IF ( nCont /= pRegion%plagInput%nCont ) THEN 
            WRITE(errorString,'(A,1X,I2,1X,A,1X,I2)') 'Specified:',nCont, & 
                  'but expected:',pRegion%plagInput%nCont        
            CALL ErrorStop(global,ERR_PLAG_NCONT_INVALID,__LINE__)          
          END IF ! nCont              
        
! ------------------------------------------------------------------------------
!       Next particle identifier
! ------------------------------------------------------------------------------ 
        
        CASE ( '# Next particle identifier' )     
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Next particle identifier...'
          END IF ! global%verbLevel          

          READ(iFile,'(I8)') pRegion%plag%nextIdNumber         
        
! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------ 
      
        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel           
      
          EXIT
      
! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------ 
      
        CASE DEFAULT
          IF ( global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
          END IF ! verbosityLevel           

          CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)      

      END SELECT ! TRIM

! ------------------------------------------------------------------------------
!     Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------  
  
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter

    END DO ! <empty> 

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%myProcid == MASTERPROC .AND. &
         global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error 
   
! ******************************************************************************
!   End
! ******************************************************************************
    
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading particle dimensions done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  
  
  END SUBROUTINE PLAG_RFLU_ReadDimensions








! ******************************************************************************
!
! Purpose: Set dimensions for Lagrangian particle datastructure 
!          at initialization stage.
!
! Description: None.
!
! Input:
!   pRegion        Pointer to region
!   nPcls          Number of particles 
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_SetDimensions(pRegion,nPcls)

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: nPcls
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_SetDimensions',&
  'PLAG_ModDimensions.F90')

! ******************************************************************************
!   Set number of particles
! ******************************************************************************

    pRegion%plag%nPcls = nPcls

! ******************************************************************************
!   End
! ******************************************************************************
    
    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_SetDimensions








! ******************************************************************************
!
! Purpose: Set maximum dimensions for Lagrangian particle datastructure.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_SetMaxDimensions(pRegion)
  
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

    INTEGER :: nPcls,nPclsMax
    REAL(RFREAL), PARAMETER :: PLAG_REALLOC_FACT = 1.20_RFREAL   
    TYPE(t_global), POINTER :: global 

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_SetMaxDimensions',&
  'PLAG_ModDimensions.F90')

! ******************************************************************************
!   Set variables 
! ******************************************************************************

    nPcls    = pRegion%plag%nPcls
    nPclsMax = pRegion%plag%nPclsMax

! ******************************************************************************
!   Set maximum dimension of particle datastructure. First ensure that 
!   a minimum value is preserved for null particle field.
! ******************************************************************************
  
    nPclsMax = NINT(PLAG_REALLOC_FACT*REAL(nPcls,KIND=RFREAL))

! ******************************************************************************
!   Update maximum value of particle datastructure
! ******************************************************************************

    pRegion%plag%nPclsMax = MAX(nPclsMax,NPCLS_TOT_MIN)

! ******************************************************************************
!   End
! ******************************************************************************
    
    CALL DeregisterFunction(global)  
  
  END SUBROUTINE PLAG_SetMaxDimensions






! ******************************************************************************
!
! Purpose: Write dimensions for Lagrangian particles.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLU_WriteDimensions(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic, &
                                 BuildFileNameUnsteady  

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion  
   
! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,RCSIdentString,sectionString, & 
                         timeString1,timeString2
    INTEGER :: errorFlag,dummy,iFile
    REAL(RFREAL) :: currentTime
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_RFLU_WriteDimensions',&
  'PLAG_ModDimensions.F90')
    
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing particle dimensions...'
    END IF ! global%verbLevel

! ==============================================================================
!   Build file name
! ==============================================================================

    iFile = IF_DIMS

    IF ( global%flowType == FLOW_UNSTEADY ) THEN 
#ifndef GENX
      currentTime = global%currentTime
#else
      IF ( global%timeStamp > 0.0_RFREAL ) THEN
        global%warnCounter = global%warnCounter + 1    

        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_NONE ) THEN        
          WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'*** WARNING ***', & 
                                        'Hard-code - read file from time zero.'
        END IF ! global%myProcid
      END IF ! global%timeStamp

      currentTime = 0.0_RFREAL ! Hard-code for GENX restart
#endif        
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.pdim', & 
                                 pRegion%iRegionGlobal,currentTime, & 
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN        
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal          
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                            global%currentTime   
      END IF ! global%verbLevel
    ELSE 
      CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.pdim', & 
                              pRegion%iRegionGlobal,iFileName) 

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal
      END IF ! global%verbLevel
    END IF ! global

! ==============================================================================
!   Open file
! ==============================================================================

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", & 
         IOSTAT=errorFlag)   
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel 

    sectionString = '# ROCPART dimensions file'  
    WRITE(iFile,'(A)') TRIM(sectionString) 

! ==============================================================================
!   Write dimensions
! ==============================================================================

    sectionString = '# Actual number of particles'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(I8)') pRegion%plag%nPcls   

    sectionString = '# Maximum number of particles'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(I8)') pRegion%plag%nPclsMax  

    sectionString = '# Number of constituents'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(I8)') pRegion%plagInput%nCont  

    sectionString = '# Next particle identifier'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(I8)') pRegion%plag%nextIdNumber      
    
! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile,'(A)') TRIM(sectionString) 
 
! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%myProcid == MASTERPROC .AND. &
         global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error
   
! ******************************************************************************
!   End
! ******************************************************************************
  
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing dimensions done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
  
  END SUBROUTINE PLAG_RFLU_WriteDimensions








! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE PLAG_ModDimensions


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModDimensions.F90,v $
! Revision 1.7  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2007/03/31 23:55:43  haselbac
! Removed Tot from names, bug fix in CalcNPclsGlobal, removed tabs
!
! Revision 1.4  2007/03/27 00:55:33  haselbac
! Substantial changes and additions
!
! Revision 1.3  2007/03/21 21:33:10  fnajjar
! Activated IO of nPclsTotGlobal when doPrint is on
!
! Revision 1.2  2007/03/20 22:03:27  fnajjar
! Included PLAG_CalcnPclsTotGlobal routine
!
! Revision 1.1  2007/03/20 17:39:12  fnajjar
! Initial import
!
! ******************************************************************************













