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
! Purpose: Collection of routines related to probes.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModProbes.F90,v 1.6 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModProbes

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_CloseProbeFiles, &
            RFLU_DecideWriteProbes, & 
            RFLU_FindProbeCells, & 
            RFLU_OpenProbeFiles, & 
            RFLU_PrintProbeInfo
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModProbes.F90,v $ $Revision: 1.6 $' 
              
! *****************************************************************************
! Routines
! *****************************************************************************

  CONTAINS
  
  
  
  
  


! ******************************************************************************
!
! Purpose: Close probe files. 
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

SUBROUTINE RFLU_CloseProbeFiles(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE (t_region), POINTER :: pRegion

! ******************************************************************************
! Locals
! ******************************************************************************

  INTEGER :: errorFlag,iProbe
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CloseProbeFiles',&
  'RFLU_ModProbes.F90')

! ******************************************************************************
! Loop over probes
! ******************************************************************************

  DO iProbe = 1,global%nProbes

! ==============================================================================
!   Check if probe region on current processor
! ==============================================================================

    IF ( global%probePos(iProbe,PROBE_REGION) == pRegion%iRegionGlobal ) THEN
      CLOSE(IF_PROBE+iProbe-1,IOSTAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__)
      END IF ! global%error
    END IF ! global
  END DO ! iProbe

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CloseProbeFiles
  
  
  
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Determine whether to write probe data to files.
!
! Description: None.
!
! Input:
!   global                      Pointer to global data
!
! Output: 
!   RFLU_DecideWrite = .TRUE.   If should write probe data to files
!   RFLU_DecideWrite = .FALSE.  If should not write probe data to files
!
! Notes: None.
!
! ******************************************************************************

LOGICAL FUNCTION RFLU_DecideWriteProbes(global)
      
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: logical1,logical2,logical3

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_DecideWriteProbes',&
  'RFLU_ModProbes.F90')

! ******************************************************************************
! Initialize
! ******************************************************************************

  RFLU_DecideWriteProbes = .FALSE.
  
! ******************************************************************************
! Determine whether should print to screen
! ******************************************************************************

! ==============================================================================
! Unsteady flow
! ==============================================================================

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    logical1 = ABS(global%timeSinceProbe-global%probeSaveTime) & 
             < 0.1_RFREAL*global%dtMin  
    logical2 = (global%timeSinceProbe > global%probeSaveTime)    
    logical3 = (global%iterSinceRestart == 1)
  
    IF ( logical1 .OR. logical2 .OR. logical3 ) THEN    
      RFLU_DecideWriteProbes = .TRUE.
    END IF ! logical1

! ==============================================================================
! Steady flow
! ==============================================================================

  ELSE    
    RFLU_DecideWriteProbes = (MOD(global%currentIter,global%probeSaveIter) == 0)

    IF ( global%currentIter == 1 ) THEN 
      RFLU_DecideWriteProbes = .TRUE.
    END IF ! global%currentIter
  END IF ! global%flowType

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_DecideWriteProbes 
  
  
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Determine cell centroids which are closest to probe coordinates.
!
! Description: None.
!
! Input: 
!  pRegion      Pointer to region
!
! Output: None.
!
! Notes: 
!   1. It is not enough to simply get closest cell from Octree, because 
!      proximity to cell centroid does not guarantee that the specified 
!      location is actually inside the cell.
!   2. If the probe location is very close to a vertex, a very large number
!      of cells may have to be searched for tetrahedral grids...
!
! ******************************************************************************

SUBROUTINE RFLU_FindProbeCells(pRegion)

  USE RFLU_ModGeometryTools, ONLY: RFLU_TestInBoundBox
  USE RFLU_ModInCellTest, ONLY: RFLU_ICT_TestInCell
  USE RFLU_ModOctree
   
  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ===============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  INTEGER :: ccSize,errorFlag,icg,iProbe
  INTEGER, DIMENSION(:), ALLOCATABLE :: cc
  REAL(RFREAL) :: delFrac,xDel,xp,xMax,xMin,yDel,yp,yMax,yMin,zDel,zp,zMax,zMin
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_FindProbeCells',&
  'RFLU_ModProbes.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Finding cells containing probes...'    
  END IF ! global%verbLevel

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pGrid => pRegion%grid

  delFrac = 0.01_RFREAL
  
  ccSize = MIN(100,pGrid%nCells) ! Must be larger than unity 

  ALLOCATE(cc(ccSize),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cc')
  END IF ! global%error

! ******************************************************************************  
! Build bounding box
! ******************************************************************************

  xMin = MINVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  xMax = MAXVAL(pGrid%xyz(XCOORD,1:pGrid%nVert))
  yMin = MINVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
  yMax = MAXVAL(pGrid%xyz(YCOORD,1:pGrid%nVert))
  zMin = MINVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))
  zMax = MAXVAL(pGrid%xyz(ZCOORD,1:pGrid%nVert))

  xDel = xMax - xMin 
  yDel = yMax - yMin
  zDel = zMax - zMin

  xMin = xMin - delFrac*xDel
  xMax = xMax + delFrac*xDel 
  yMin = yMin - delFrac*yDel
  yMax = yMax + delFrac*yDel 
  zMin = zMin - delFrac*zDel
  zMax = zMax + delFrac*zDel  

! ******************************************************************************  
! Build Octree with cell centroids: NOTE only for actual cells
! ******************************************************************************

  CALL RFLU_CreateOctree(global,pGrid%nCells)
  
  CALL RFLU_BuildOctree(pGrid%cofg(XCOORD,1:pGrid%nCells), & 
                        pGrid%cofg(YCOORD,1:pGrid%nCells), & 
                        pGrid%cofg(ZCOORD,1:pGrid%nCells), & 
                        xMin,xMax,yMin,yMax,zMin,zMax)

! ******************************************************************************
! Loop over probes
! ******************************************************************************

  DO iProbe = 1,global%nProbes

! ==============================================================================  
!   Test probe location against bounding box of partition
! ==============================================================================  
  
    xp = global%probeXyz(iProbe,1)
    yp = global%probeXyz(iProbe,2)
    zp = global%probeXyz(iProbe,3)

    IF ( RFLU_TestInBoundBox(global,xp,yp,zp,xMin,xMax,yMin,yMax, & 
                             zMin,zMax) .EQV. .TRUE. ) THEN  
                               
! ------------------------------------------------------------------------------  
!     Query octree to get closest cells  
! ------------------------------------------------------------------------------ 
    
      CALL RFLU_QueryOctree(xp,yp,zp,ccSize,cc)
        
! ------------------------------------------------------------------------------
!     Test cells obtained from Octree if they contain specified location
! ------------------------------------------------------------------------------        
           
      cellLoop: DO icg = 1,ccSize                 
        IF ( RFLU_ICT_TestInCell(pRegion,xp,yp,zp,cc(icg)) .EQV. .TRUE. ) THEN
          global%probePos(iProbe,PROBE_REGION) = pRegion%iRegionGlobal
          global%probePos(iProbe,PROBE_CELL)   = cc(icg)

          EXIT cellLoop
        END IF ! RFLU_ICT_TestInCell        
      END DO cellLoop
    END IF ! RFLU_TestInBoundNox 
  END DO ! iProbe

! ******************************************************************************
! Destroy Octree and deallocate temporary memory
! ******************************************************************************

  CALL RFLU_DestroyOctree(global)

  DEALLOCATE(cc,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cc')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN          
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Finding cells containing probes done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_FindProbeCells  
  
  
  
  
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Open probe files and position files at right line when restarting. 
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

SUBROUTINE RFLU_OpenProbeFiles(pRegion)

  USE ModBuildFileNames, ONLY: BuildFileNamePlain
  USE ModTools, ONLY: FloatLess

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE (t_region), POINTER :: pRegion

! ******************************************************************************
! Locals
! ******************************************************************************

  LOGICAL :: fileAppend,fileExists
  CHARACTER(CHRLEN+9) :: fname
  INTEGER :: errorFlag,iProbe,probeIter
  REAL(RFREAL) :: probeTime
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_OpenProbeFiles',&
  'RFLU_ModProbes.F90')

! ******************************************************************************
! Initialize
! ******************************************************************************
  
  probeTime = HUGE(1.0_RFREAL)

! ******************************************************************************
! Loop over probes
! ******************************************************************************

  DO iProbe = 1,global%nProbes

! ==============================================================================
!   Check if probe region on current processor
! ==============================================================================

    IF ( global%probePos(iProbe,PROBE_REGION) == pRegion%iRegionGlobal ) THEN

! ------------------------------------------------------------------------------
!     Generate file name
! ------------------------------------------------------------------------------

      WRITE(fname,'(A,I4.4)') TRIM(global%outDir)// & 
                              TRIM(global%casename)//'.prb_',iProbe

! ------------------------------------------------------------------------------
!     Open file with appropriate options
! ------------------------------------------------------------------------------

      IF ( (global%flowType    == FLOW_UNSTEADY .AND. &
            global%currentTime > 0.0_RFREAL) .OR.  &
           (global%flowType    == FLOW_STEADY .AND. & 
            global%currentIter > 1) ) THEN
        INQUIRE(FILE=fname,EXIST=fileExists)

        IF ( fileExists .EQV. .TRUE. ) THEN
          fileAppend = .TRUE.

          OPEN(IF_PROBE+iProbe-1,FILE=fname,FORM='FORMATTED',STATUS='OLD', &
               POSITION='APPEND',IOSTAT=errorFlag)
        ELSE
          fileAppend = .FALSE.
        
          OPEN(IF_PROBE+iProbe-1,FILE=fname,FORM='FORMATTED',STATUS='UNKNOWN', &
               IOSTAT=errorFlag)
        END IF ! fileExists
      ELSE
        fileAppend = .FALSE.

        OPEN(IF_PROBE+iProbe-1,FILE=fname,FORM='FORMATTED',STATUS='UNKNOWN', &
             IOSTAT=errorFlag)
      END IF ! global

      global%error = errorFlag
      IF (global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname))
      END IF ! global%error
 
! ------------------------------------------------------------------------------
!     Back up to right time when appending to existing file. NOTE if encounter
!     errors when backspacing or reading file (for example because it may be 
!     empty), just quit loop instead of treating error properly.
! ------------------------------------------------------------------------------

      IF ( fileAppend .EQV. .TRUE. ) THEN
      
! ----- Unsteady flow ----------------------------------------------------------
      
        IF ( global%flowType == FLOW_UNSTEADY ) THEN 
          emptyLoopUnsteady: DO 
            BACKSPACE(IF_PROBE+iProbe-1,IOSTAT=errorFlag)
            IF ( errorFlag /= ERR_NONE ) THEN 
              EXIT emptyLoopUnsteady
            END IF ! errorFlag
          
            READ(IF_PROBE+iProbe-1,*,IOSTAT=errorFlag) probeTime
            IF ( errorFlag /= ERR_NONE ) THEN 
              EXIT emptyLoopUnsteady
            END IF ! errorFlag

            IF ( FloatLess(probeTime,global%currentTime) .EQV. .TRUE. ) THEN           
              EXIT emptyLoopUnsteady
            ELSE 
              BACKSPACE(IF_PROBE+iProbe-1,IOSTAT=errorFlag)
              IF ( errorFlag /= ERR_NONE ) THEN 
                EXIT emptyLoopUnsteady
              END IF ! errorFlag              
            END IF ! probeTime
          END DO emptyLoopUnsteady

! ----- Steady flow ------------------------------------------------------------

        ELSE 
          emptyLoopSteady: DO 
            BACKSPACE(IF_PROBE+iProbe-1,IOSTAT=errorFlag)
            IF ( errorFlag /= ERR_NONE ) THEN 
              EXIT emptyLoopSteady
            END IF ! errorFlag
          
            READ(IF_PROBE+iProbe-1,*,IOSTAT=errorFlag) probeIter
            IF ( errorFlag /= ERR_NONE ) THEN 
              EXIT emptyLoopSteady
            END IF ! errorFlag

            IF ( probeIter < global%currentIter ) THEN           
              EXIT emptyLoopSteady
            ELSE 
              BACKSPACE(IF_PROBE+iProbe-1,IOSTAT=errorFlag)
              IF ( errorFlag /= ERR_NONE ) THEN 
                EXIT emptyLoopSteady
              END IF ! errorFlag
            END IF ! probeTime
          END DO emptyLoopSteady         
        END IF ! global%flowType 
      END IF ! fileAppend
    END IF ! global
  END DO ! iProbe

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_OpenProbeFiles
  
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Print information on probe cells.
!
! Description: None.
!
! Input: 
!  global               Pointer to global data
!
! Output: None.
!
! Notes: 
!   1. It is important to note that the use of the max reduction requires the 
!      global%probePos array to be initialized with negative integers (see 
!      readProbeSection.F90). Hence if the entries of a given probe are still
!      equal to the initial values, know that probe was not located in any of
!      the regions.
!
! ******************************************************************************

SUBROUTINE RFLU_PrintProbeInfo(global)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================
! Arguments 
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: errorString
  INTEGER :: errorFlag,iProbe
  INTEGER, DIMENSION(:), ALLOCATABLE :: globalVals,localVals

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_PrintProbeInfo',&
  'RFLU_ModProbes.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing probe information...'
  END IF ! global%verbLevel

! ******************************************************************************
! Gather information
! ******************************************************************************

! ==============================================================================
! Allocate temporary memory
! ==============================================================================

  ALLOCATE(localVals(global%nProbes),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'localVals')
  END IF ! global%error

  ALLOCATE(globalVals(global%nProbes),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'globalVals')
  END IF ! global%error

! ==============================================================================
! Peform reduction operations
! ==============================================================================
  
  DO iProbe = 1,global%nProbes
    localVals(iProbe) = global%probePos(iProbe,PROBE_REGION)
  END DO ! iProbe
   
  CALL MPI_Allreduce(localVals,globalVals,global%nProbes,MPI_INTEGER,MPI_MAX, &
                     global%mpiComm,errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
  END IF ! global%errorFlag
   
  DO iProbe = 1,global%nProbes
    global%probePos(iProbe,PROBE_REGION) = globalVals(iProbe)
  END DO ! iProbe   

  DO iProbe = 1,global%nProbes
    localVals(iProbe) = global%probePos(iProbe,PROBE_CELL)
  END DO ! iProbe
   
  CALL MPI_Allreduce(localVals,globalVals,global%nProbes,MPI_INTEGER,MPI_MAX, &
                     global%mpiComm,errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
  END IF ! global%errorFlag
  
  DO iProbe = 1,global%nProbes
    global%probePos(iProbe,PROBE_CELL) = globalVals(iProbe)
  END DO ! iProbe 

! ==============================================================================
! Deallocate temporary memory
! ==============================================================================

  DEALLOCATE(localVals,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'localVals')
  END IF ! global%error

  DEALLOCATE(globalVals,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'globalVals')
  END IF ! global%error

! ******************************************************************************
! If probe position values are still equal to initial values, could not find 
! a cell containing probe position
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC ) THEN 
    DO iProbe = 1,global%nProbes    
      IF ( (global%probePos(iProbe,PROBE_REGION) == CRAZY_VALUE_INT) .OR. & 
           (global%probePos(iProbe,PROBE_CELL  ) == CRAZY_VALUE_INT) ) THEN
        WRITE(errorString,'(A,1X,I3)') 'Probe:',iProbe
        CALL ErrorStop(global,ERR_PROBE_LOCATION,__LINE__,TRIM(errorString))
      END IF ! global%probePos
    END DO ! iProbe
  END IF ! global%myProcid

! ******************************************************************************
! Print information
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,5X,A,3X,A,8X,A)') SOLVER_NAME,'#','Region','Cell'
       
    DO iProbe = 1,global%nProbes
      WRITE(STDOUT,'(A,3X,I3,3X,I6,3X,I9)') SOLVER_NAME,iProbe, &
        global%probePos(iProbe,PROBE_REGION),global%probePos(iProbe,PROBE_CELL)         
    END DO ! iProbe
  END IF ! global%myProcid  

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing probe information done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PrintProbeInfo  
  
  
  
  
  
  
! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModProbes


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModProbes.F90,v $
! Revision 1.6  2008/12/06 08:44:23  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:34  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/12/15 13:25:08  haselbac
! Fixed bug in format statement, found by ifort
!
! Revision 1.3  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.2  2005/12/24 21:30:45  haselbac
! Adapted to changes in ICT
!
! Revision 1.1  2005/04/29 12:41:04  haselbac
! Initial revision
!
! ******************************************************************************
  











