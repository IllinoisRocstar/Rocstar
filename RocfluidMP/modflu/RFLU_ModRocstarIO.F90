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
! Purpose: Collection of routines related to GENX I/O.
!
! Description: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModGENXIO.F90,v 1.15 2008/12/06 08:44:21 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModRocstarIO

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid 
  USE ModMPI 

  USE ModBuildFileNames, ONLY: BuildRegionIdString

  USE RFLU_ModRocstarUtils

  USE RFLU_ModInterfacesLibrary, ONLY: RFLU_TestIsFirstRegion

  IMPLICIT NONE

  INCLUDE 'comf90.h'

  PRIVATE
  PUBLIC :: RFLU_GENX_BuildRocinPaneStrings, & 
            RFLU_GENX_CloseRocinCtrlFiles, &
            RFLU_GENX_DecideReadFile, & 
            RFLU_GENX_DecideWriteFile, &                         
            RFLU_GENX_GetDataFlow, &
            RFLU_GENX_GetDataGSpeedsSurf, &
            RFLU_GENX_GetDataGSpeedsVol, & 
            RFLU_GENX_GetDataInterf, &
            RFLU_GENX_GetDimensions, &
            RFLU_GENX_GetDimensionsDerived, &  
            RFLU_GENX_GetGlobalData, &  
            RFLU_GENX_GetGrid, &
            RFLU_GENX_OpenRocinCtrlFiles, &
            RFLU_GENX_PutDataFlow, &
            RFLU_GENX_PutDataGSpeedsSurf, &
            RFLU_GENX_PutDataGSpeedsVol, &
            RFLU_GENX_PutDataInterf, &                           
            RFLU_GENX_PutGrid, & 
            RFLU_GENX_PutGridSurf, & 
            RFLU_GENX_PutGridVol, & 
            RFLU_GENX_ReadWindow, & 
            RFLU_GENX_WriteRocinCtrlFiles

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================  
! Private
! ==============================================================================  

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModRocstarIO.F90,v $ $Revision: 1.15 $'       
  
! ==============================================================================  
! Public
! ==============================================================================  
  
  INTEGER, PARAMETER, PUBLIC :: GENX_WINDOW_TYPE_SURF = 1, & 
                                GENX_WINDOW_TYPE_VOL  = 2 

! ******************************************************************************
! Contained routines
! ******************************************************************************

  CONTAINS 

 
  
  
  


! ******************************************************************************
!
! Purpose: Build pane strings for Rocin control file.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!
! Output: 
!   paneStringVol       String with pane ids of volume data
!   paneStringSurf      String with pane ids of surface data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_BuildRocinPaneStrings(pRegion,paneStringVol, & 
                                             paneStringSurf)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    CHARACTER(*), INTENT(INOUT) :: paneStringSurf,paneStringVol
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  

    CHARACTER(CHRLEN) :: tempString
    INTEGER :: iPatch,paneId    
    TYPE(t_global), POINTER :: global    
        
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_GENX_BuildRocinPaneStrings', & 
                          'RFLU_ModRocstarIO.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building Rocin pane strings...'
    END IF ! global%verbLevel
        
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid

! ******************************************************************************
!   Write file
! ******************************************************************************

    CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,0,paneId)

    WRITE(tempString,*) paneId    
    WRITE(paneStringVol,'(A)') TRIM(paneStringVol)//' '// & 
                               ADJUSTL(TRIM(tempString))

    DO iPatch = 1,pRegion%grid%nPatches
      CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,iPatch,paneId)      

      WRITE(tempString,*) paneId    
      WRITE(paneStringSurf,'(A)') TRIM(paneStringSurf)//' '// & 
                                  ADJUSTL(TRIM(tempString))    
    END DO ! iPatch
  
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building Rocin pane strings done.'
    END IF ! global%verbLevel
      
    CALL DeregisterFunction(global)  
  
  END SUBROUTINE RFLU_GENX_BuildRocinPaneStrings  
  
  
  
  
  
  

! ******************************************************************************
!
! Purpose: Close control files for Rocin.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_CloseRocinCtrlFiles(global)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_global), POINTER :: global

! ==============================================================================  
!   Locals 
! ==============================================================================  

    CHARACTER(CHRLEN) :: iFileName
    CHARACTER(GENX_TIME_STRING_LEN) :: timeString
    INTEGER :: errorFlag    
        
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_GENX_CloseRocinCtrlFiles', & 
                          'RFLU_ModRocstarIO.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing Rocin control files...'
    END IF ! global%verbLevel
     
! ******************************************************************************
!   Close files
! ******************************************************************************

! ==============================================================================  
!   Volume file 
! ==============================================================================  
 
    CALL RFLU_GENX_BuildTimeString(global%currentTime,timeString)
    iFileName = './Rocflu/Rocin/fluid_in_'//TRIM(timeString)//'.txt'

    CLOSE(IF_CTRL_VOL,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName) 
    END IF ! global%error

! ==============================================================================  
!   Surface file 
! ==============================================================================  

    CALL RFLU_GENX_BuildTimeString(global%currentTime,timeString)
    iFileName = './Rocflu/Rocin/ifluid_in_'//TRIM(timeString)//'.txt'

    CLOSE(IF_CTRL_SURF,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName) 
    END IF ! global%error
 
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing Rocin control files done.'
    END IF ! global%verbLevel
      
    CALL DeregisterFunction(global)  
  
  END SUBROUTINE RFLU_GENX_CloseRocinCtrlFiles
  
  
  



! ******************************************************************************
!
! Purpose: Decide whether GENX files should be read.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  LOGICAL FUNCTION RFLU_GENX_DecideReadFile(global)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_GENX_DecideReadFile',&
  'RFLU_ModRocstarIO.F90')
    
! ******************************************************************************
!   Decide whether need to read files
! ******************************************************************************
    
    IF ( global%moduleType == MODULE_TYPE_PART ) THEN
      RFLU_GENX_DecideReadFile = .TRUE.            
    ELSE IF ( global%moduleType == MODULE_TYPE_SOLVER ) THEN 
      RFLU_GENX_DecideReadFile = .TRUE.
    ELSE IF ( global%moduleType == MODULE_TYPE_POSTPROC ) THEN 
      RFLU_GENX_DecideReadFile = .TRUE. 
    ELSE IF ( global%moduleType == MODULE_TYPE_INIT ) THEN 
      RFLU_GENX_DecideReadFile = .TRUE.    
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! global%moduleType

! ******************************************************************************
!   End
! ******************************************************************************
  
    CALL DeregisterFunction(global)  
  
  END FUNCTION RFLU_GENX_DecideReadFile






! ******************************************************************************
!
! Purpose: Decide whether files should be written.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  LOGICAL FUNCTION RFLU_GENX_DecideWriteFile(global)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_GENX_DecideWriteFile',&
  'RFLU_ModRocstarIO.F90')
    
! ******************************************************************************
!   Decide whether need to write files
! ******************************************************************************
    
    IF ( global%moduleType == MODULE_TYPE_PART ) THEN
      RFLU_GENX_DecideWriteFile = .TRUE.
    ELSE IF ( global%moduleType == MODULE_TYPE_SOLVER ) THEN 
      RFLU_GENX_DecideWriteFile = .FALSE.
    ELSE IF ( global%moduleType == MODULE_TYPE_POSTPROC ) THEN 
      RFLU_GENX_DecideWriteFile = .TRUE.  
    ELSE IF ( global%moduleType == MODULE_TYPE_INIT ) THEN 
      RFLU_GENX_DecideWriteFile = .TRUE.    
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! global%moduleType

! ******************************************************************************
!   End
! ******************************************************************************
  
    CALL DeregisterFunction(global)  
  
  END FUNCTION RFLU_GENX_DecideWriteFile







! ******************************************************************************
!
! Purpose: Get flow data through Roccom.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetDataFlow(pRegion)
  
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
    
    CHARACTER(CHRLEN) :: winName,winNameIn
    INTEGER :: handleIn,handleObtain,handleOut
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global

    handleObtain = global%handleObtain
    
    pRegion%mixt%cvState = CV_MIXT_STATE_CONS

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting flow data...'
    END IF ! global%verbLevel 
        
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid

! ******************************************************************************
!   Get flow data
! ******************************************************************************

    winNameIn = global%volWinNameInput
    winName   = global%volWinName  
    
! ==============================================================================
!   Conserved variables
! ==============================================================================    
    
    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.rhof')
    handleOut = COM_get_dataitem_handle(TRIM(winName)//'.rhof')
    CALL COM_call_function(handleObtain,2,handleIn,handleOut)

    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.rhovf')
    handleOut = COM_get_dataitem_handle(TRIM(winName)//'.rhovf')
    CALL COM_call_function(handleObtain,2,handleIn,handleOut)

    handleIn  = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.rhoEf')    
    handleOut = COM_get_dataitem_handle(TRIM(winName)//'.rhoEf')
    CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    
! ==============================================================================
!   Dependent variables NOT obtained from Roccom because will be set later
! ==============================================================================    

! ******************************************************************************
!   End
! ******************************************************************************
 
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting flow data done.'
    END IF ! global%verbLevel   
  
  END SUBROUTINE RFLU_GENX_GetDataFlow 








! ******************************************************************************
!
! Purpose: Get surface grid speed data through Roccom.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetDataGSpeedsSurf(pRegion)
  
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
    
    CHARACTER(CHRLEN) :: winName,winNameIn
    INTEGER :: handleIn,handleObtain,handleOut
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global

    handleObtain = global%handleObtain
    
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting surface grid speeds...'
    END IF ! global%verbLevel      
    
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid

! ******************************************************************************
!   Get grid speeds
! ******************************************************************************
        
    winNameIn = global%surfWinNameInput
    winName   = global%surfWinName  
    
    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.gs')
    IF ( handleIn > 0 ) THEN    
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.gs')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! handleIn

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting surface grid speeds done.'
    END IF ! global%verbLevel  
  
  END SUBROUTINE RFLU_GENX_GetDataGSpeedsSurf








! ******************************************************************************
!
! Purpose: Get volume grid speed data through Roccom.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetDataGSpeedsVol(pRegion)
  
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
    
    CHARACTER(CHRLEN) :: winName,winNameIn
    INTEGER :: handleIn,handleObtain,handleOut
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global

    handleObtain = global%handleObtain
    
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting volume grid speeds...'
    END IF ! global%verbLevel      
    
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid

! ******************************************************************************
!   Get grid speeds
! ******************************************************************************
    
    winNameIn = global%volWinNameInput
    winName   = global%volWinName  
    
    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.gs')    
    IF ( handleIn > 0 ) THEN
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.gs')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! handleIn
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting volume grid speeds done.'
    END IF ! global%verbLevel  
  
  END SUBROUTINE RFLU_GENX_GetDataGSpeedsVol





! ******************************************************************************
!
! Purpose: Get interface data through Roccom.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetDataInterf(pRegion)
  
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
    
    CHARACTER(CHRLEN) :: winName,winNameIn
    INTEGER :: handleIn,handleObtain,handleOut
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global

    handleObtain = global%handleObtain

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting interface data...'
    END IF ! global%verbLevel  
      
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid

! ******************************************************************************
!   Get interface data
! ******************************************************************************

    winNameIn = global%surfWinNameInput
    winName   = global%surfWinName  
    
! ==============================================================================
!   Input data
! ==============================================================================    
    
    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.du_alp')
    IF ( handleIn > 0 ) THEN 
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.du_alp')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! handleIn

    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.rhofvf_alp')
    IF ( handleIn > 0 ) THEN 
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.rhofvf_alp')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! handleIn

    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.mdot_alp')
    IF ( handleIn > 0 ) THEN 
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.mdot_alp')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! handleIn

    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.Tflm_alp')
    IF ( handleIn > 0 ) THEN 
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.Tflm_alp')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! handleIn
    
    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.Tb_alp')
    IF ( handleIn > 0 ) THEN 
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.Tb_alp')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! handleIn
    
! ==============================================================================
!   Output data
! ==============================================================================    
        
    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.nfAlp')
    IF ( handleIn > 0 ) THEN 
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.nfAlp')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! handleIn
    
    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.rhofAlp')
    IF ( handleIn > 0 ) THEN 
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.rhofAlp')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! handleIn
    
    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.pf')
    IF ( handleIn > 0 ) THEN
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.pf')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! handleIn
    
    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.tf')
    IF ( handleIn > 0 ) THEN
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.tf')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! handleIn
    
    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.qc')
    IF ( handleIn > 0 ) THEN
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.qc')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! handleIn

    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.qr')
    IF ( handleIn > 0 ) THEN
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.qr')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! handleIn

    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.tempf')
    IF ( handleIn > 0 ) THEN
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.tempf')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! handleIn
    
    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.bFlag')
    IF ( handleIn > 0 ) THEN
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.bFlag')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! handleIn    
    
! ******************************************************************************
!   End
! ******************************************************************************
  
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting interface data done.'
    END IF ! global%verbLevel  
  
  END SUBROUTINE RFLU_GENX_GetDataInterf







! ******************************************************************************
!
! Purpose: Get dimensions through Roccom.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetDimensions(pRegion)
  
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
    
    CHARACTER(CHRLEN) :: connTable,winNameIn
    CHARACTER, DIMENSION(:), POINTER :: connTables
    INTEGER :: handle,iBeg,iConnTable,iEnd,iPane,iPatch,nChars,nConnTables, & 
               nActual,nDummy,nPanes,nVirtual,paneId
    INTEGER, DIMENSION(:), POINTER :: paneList        
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global

    pGrid => pRegion%grid

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting dimensions...'
    END IF ! global%verbLevel 
        
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid

! ******************************************************************************
!   Get dimensions
! ******************************************************************************

! ==============================================================================
!   Volume
! ==============================================================================    
    
    winNameIn = global%volWinNameInput  

    CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,0,paneId)    
           
! ------------------------------------------------------------------------------
!   Vertices
! ------------------------------------------------------------------------------    
    
    CALL COM_get_size(TRIM(winNameIn)//'.nc',paneId,pGrid%nVert,nVirtual)
    
    pGrid%nVertTot = pGrid%nVert + nVirtual    

! ------------------------------------------------------------------------------
!   Connectivity tables
! ------------------------------------------------------------------------------    
    
    pGrid%nTets = 0 
    pGrid%nHexs = 0 
    pGrid%nPris = 0 
    pGrid%nPyrs = 0 
    
    pGrid%nTetsTot = 0
    pGrid%nHexsTot = 0
    pGrid%nPrisTot = 0
    pGrid%nPyrsTot = 0 
    
    CALL COM_get_connectivities(TRIM(winNameIn),paneId,nConnTables,connTables)
    
    iBeg = 1
    
    DO iConnTable = 1,nConnTables 
      iEnd = iBeg
      nChars = 0
      
      DO 
        IF ( connTables(iEnd) /= ' ' ) THEN 
          nChars = nChars + 1     
          connTable(nChars:nChars) = connTables(iEnd)
          iEnd = iEnd + 1
        ELSE 
          EXIT 
        END IF ! connTables
        
        IF ( iEnd > UBOUND(connTables,1) ) THEN 
          EXIT 
        END IF ! iEnd
      END DO ! <emptyLoop>

      SELECT CASE ( connTable(1:nChars) )       
        CASE ( ':T4:real' ) 
          CALL COM_get_size(TRIM(winNameIn)//'.:T4:real',paneId,pGrid%nTets, & 
                            nDummy)
        CASE ( ':T4:virtual' ) 
          CALL COM_get_size(TRIM(winNameIn)//'.:T4:virtual',paneId, &
                            pGrid%nTetsTot,nDummy)
        CASE ( ':H8:real' ) 
          CALL COM_get_size(TRIM(winNameIn)//'.:H8:real',paneId,pGrid%nHexs, &
                            nDummy)     
        CASE ( ':H8:virtual' ) 
          CALL COM_get_size(TRIM(winNameIn)//'.:H8:virtual',paneId, &
                            pGrid%nHexsTot,nDummy)      
        CASE ( ':W6:real' ) 
          CALL COM_get_size(TRIM(winNameIn)//'.:W6:real',paneId,pGrid%nPris, &
                            nDummy)     
        CASE ( ':W6:virtual' )  
          CALL COM_get_size(TRIM(winNameIn)//'.:W6:virtual',paneId, &
                            pGrid%nPrisTot,nDummy)      
        CASE ( ':P5:real' ) 
          CALL COM_get_size(TRIM(winNameIn)//'.:P5:real',paneId,pGrid%nPyrs, &
                            nDummy)     
        CASE ( ':P5:virtual' )
          CALL COM_get_size(TRIM(winNameIn)//'.:P5:virtual',paneId, &
                            pGrid%nPyrsTot,nDummy)      
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! connTable            
    END DO ! iConnTable     
                  
    pGrid%nTetsTot = pGrid%nTetsTot + pGrid%nTets
    pGrid%nHexsTot = pGrid%nHexsTot + pGrid%nHexs
    pGrid%nPrisTot = pGrid%nPrisTot + pGrid%nPris
    pGrid%nPyrsTot = pGrid%nPyrsTot + pGrid%nPyrs                                        
                                
    CALL COM_free_buffer(connTables)
                  
! ==============================================================================
!   Surface
! ==============================================================================    

    winNameIn = global%surfWinNameInput  

! -----------------------------------------------------------------------------
!   Get number of patches. NOTE get panes of surface window only.
! -----------------------------------------------------------------------------

    CALL COM_get_panes(TRIM(winNameIn),nPanes,paneList)
    
! DEBUG
!    WRITE(*,*) 'nPanes:',nPanes
!    
!    pGrid%nPatches = nPanes
!    
!    DO iPane = 1,nPanes
!      WRITE(*,*) iPane,paneList(iPane)
!    END DO ! iPane
! END DEBUG    
    
    CALL COM_free_buffer(paneList)    
    
! -----------------------------------------------------------------------------
!   Loop over patches
! -----------------------------------------------------------------------------
    
! DEBUG
!    WRITE(*,*) 'nPatches:',pGrid%nPatches
! END DEBUG    
    
    DO iPatch = 1,pGrid%nPatches
      pGrid%patchDimens(PATCH_DIMENS_NBTRIS    ,iPatch) = 0 
      pGrid%patchDimens(PATCH_DIMENS_NBTRISTOT ,iPatch) = 0       
      pGrid%patchDimens(PATCH_DIMENS_NBQUADS   ,iPatch) = 0 
      pGrid%patchDimens(PATCH_DIMENS_NBQUADSTOT,iPatch) = 0       
        
      CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,iPatch,paneId)
      
      CALL COM_get_connectivities(TRIM(winNameIn),paneId,nConnTables, & 
                                  connTables)
        
! --- Connectivity tables -----------------------------------------------------
    
      iBeg = 1
    
      DO iConnTable = 1,nConnTables 
        iEnd = iBeg
        nChars = 0

        DO 
          IF ( connTables(iEnd) /= ' ' ) THEN 
            nChars = nChars + 1   
            connTable(nChars:nChars) = connTables(iEnd)
            iEnd = iEnd + 1
          ELSE 
            EXIT 
          END IF ! connTables

          IF ( iEnd > UBOUND(connTables,1) ) THEN 
            EXIT 
          END IF ! iEnd
        END DO ! <emptyLoop>
        
! DEBUG
!      WRITE(*,*) connTable(1:nChars)
! END DEBUG     
        
        SELECT CASE ( connTable(1:nChars) )       
          CASE ( ':t3:real' ) 
            CALL COM_get_size(TRIM(winNameIn)//'.:t3:real',paneId,nActual, & 
                              nDummy)
            pGrid%patchDimens(PATCH_DIMENS_NBTRIS,iPatch) = nActual
          CASE ( ':t3:virtual' ) 
            CALL COM_get_size(TRIM(winNameIn)//'.:t3:virtual',paneId,nVirtual, & 
                              nDummy)
            pGrid%patchDimens(PATCH_DIMENS_NBTRISTOT,iPatch) = nVirtual     
          CASE ( ':q4:real' ) 
            CALL COM_get_size(TRIM(winNameIn)//'.:q4:real',paneId,nActual, &  
                              nDummy)
            pGrid%patchDimens(PATCH_DIMENS_NBQUADS,iPatch) = nActual    
          CASE ( ':q4:virtual' ) 
            CALL COM_get_size(TRIM(winNameIn)//'.:q4:virtual',paneId,nVirtual, & 
                              nDummy)
            pGrid%patchDimens(PATCH_DIMENS_NBQUADSTOT,iPatch) = nVirtual
          CASE DEFAULT 
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! connTable  
      END DO ! iConnTable      
                    
      pGrid%patchDimens(PATCH_DIMENS_NBTRISTOT,iPatch) & 
        = pGrid%patchDimens(PATCH_DIMENS_NBTRIS   ,iPatch) & 
        + pGrid%patchDimens(PATCH_DIMENS_NBTRISTOT,iPatch) 
                                            
      pGrid%patchDimens(PATCH_DIMENS_NBQUADSTOT,iPatch) & 
        = pGrid%patchDimens(PATCH_DIMENS_NBQUADS   ,iPatch) & 
        + pGrid%patchDimens(PATCH_DIMENS_NBQUADSTOT,iPatch)
        
      CALL COM_free_buffer(connTables)
      
! DEBUG      
!      WRITE(*,*) pGrid%patchDimens(PATCH_DIMENS_NBTRIS    ,iPatch), &  
!                 pGrid%patchDimens(PATCH_DIMENS_NBTRISTOT ,iPatch), &      
!                 pGrid%patchDimens(PATCH_DIMENS_NBQUADS   ,iPatch), &  
!                 pGrid%patchDimens(PATCH_DIMENS_NBQUADSTOT,iPatch)    
! END DEBUG              
    END DO ! iPatch
                        
! ******************************************************************************
!   End
! ******************************************************************************
 
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting dimensions done.'
    END IF ! global%verbLevel   
  
  END SUBROUTINE RFLU_GENX_GetDimensions








! ******************************************************************************
!
! Purpose: Get derived dimensions through Roccom.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Need to get derived dimensions such as number of vertices on a given 
!      patch before registration. 
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetDimensionsDerived(pRegion)
  
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
    
    CHARACTER(CHRLEN) :: winNameIn
    INTEGER :: iPatch,nVertTot,nVertVirtual,paneId
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global

    pGrid => pRegion%grid

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting derived dimensions...'
    END IF ! global%verbLevel 
        
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid

! ******************************************************************************
!   Get dimensions
! ******************************************************************************

    winNameIn = global%surfWinNameInput  
    
! ==============================================================================
!   Loop over patches
! ==============================================================================    
    
    DO iPatch = 1,pGrid%nPatches    
      pPatch => pRegion%patches(iPatch)
      
      CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,iPatch,paneId)

      CALL COM_get_size(TRIM(winNameIn)//'.nc',paneId,nVertTot,nVertVirtual)
 
      pPatch%nBVert    = nVertTot - nVertVirtual
      pPatch%nBVertTot = nVertTot
    END DO ! iPatch
        
! ******************************************************************************
!   End
! ******************************************************************************
 
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting derived dimensions done.'
    END IF ! global%verbLevel   
  
  END SUBROUTINE RFLU_GENX_GetDimensionsDerived



  






! ******************************************************************************
!
! Purpose: Get grid through Roccom.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetGrid(pRegion)
  
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
    
    CHARACTER(CHRLEN) :: winName,winNameIn 
    INTEGER :: handleIn,handleObtain,handleOut
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
 
    pGrid => pRegion%grid

    handleObtain = global%handleObtain

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting grid...'
    END IF ! global%verbLevel
        
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid

! ******************************************************************************
!   Get volume grid data
! ******************************************************************************

    winNameIn = global%volWinNameInput
    winName   = global%volWinName  

! ==============================================================================
!   Coordinates and connectivity
! ==============================================================================

    handleIn  = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.pmesh')
    handleOut = COM_get_dataitem_handle(TRIM(winName)//'.pmesh')
    CALL COM_call_function(handleObtain,2,handleIn,handleOut)
  
! ******************************************************************************
!   Get surface grid data
! ******************************************************************************

    winNameIn = global%surfWinNameInput
    winName   = global%surfWinName  

! ==============================================================================
!   Coordinates and predefined connectivity
! ==============================================================================

    handleIn  = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.pmesh')
    handleOut = COM_get_dataitem_handle(TRIM(winName)//'.pmesh')
    CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    
! ==============================================================================
!   Other quantities. NOTE IF statement for cnstr_type to retain backward-
!   compatibility.
! ==============================================================================

    IF ( pGrid%nPatches > 0 ) THEN 
      handleIn  = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.bcflag')
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.bcflag')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)

      handleIn  = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.patchNo')
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.patchNo')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    
! Use MOVEDIR from the *.bc file to get cnstr_type; don't use HDF the values
!
!      handleIn  = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.cnstr_type')
!      IF ( handleIn > 0 ) THEN 
!        handleOut = COM_get_dataitem_handle(TRIM(winName)//'.cnstr_type')
!        CALL COM_call_function(handleObtain,2,handleIn,handleOut)    
!      END IF ! handleIn

      handleIn  = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.t3g:real')
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.t3g:real')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)

      handleIn  = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.t3g:virtual')
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.t3g:virtual')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)

      handleIn  = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.q4g:real')
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.q4g:real')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)

      handleIn  = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.q4g:virtual')
      handleOut = COM_get_dataitem_handle(TRIM(winName)//'.q4g:virtual')
      CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    END IF ! pGrid%nPatches

! ******************************************************************************
!   End
! ******************************************************************************
  
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting grid done.'
    END IF ! global%verbLevel  
  
  END SUBROUTINE RFLU_GENX_GetGrid  







! ******************************************************************************
!
! Purpose: Open control files for Rocin.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_OpenRocinCtrlFiles(global)
    
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_global), POINTER :: global

! ==============================================================================  
!   Locals 
! ==============================================================================  

    CHARACTER(CHRLEN) :: iFileName
    CHARACTER(GENX_TIME_STRING_LEN) :: timeString
    INTEGER :: errorFlag    
        
! ******************************************************************************
!   Start
! ******************************************************************************
    CALL RegisterFunction(global,'RFLU_GENX_OpenRocinCtrlFiles', & 
                          'RFLU_ModRocstarIO.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening Rocin control files...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Open files
! ******************************************************************************

! ==============================================================================  
!   Volume file 
! ==============================================================================  
        
    CALL RFLU_GENX_BuildTimeString(global%currentTime,timeString)
    iFileName = TRIM(global%outDirHDF)//'fluid_in_'//TRIM(timeString)//'.txt'

    OPEN(IF_CTRL_VOL,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", & 
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName) 
    END IF ! global%error    

! ==============================================================================  
!   Surface file 
! ==============================================================================  
     
    CALL RFLU_GENX_BuildTimeString(global%currentTime,timeString)
    iFileName = TRIM(global%outDirHDF)//'ifluid_in_'//TRIM(timeString)//'.txt'

    OPEN(IF_CTRL_SURF,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", & 
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName) 
    END IF ! global%error    
     
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Opening Rocin control files done.'
    END IF ! global%verbLevel
      
    CALL DeregisterFunction(global)  
  
  END SUBROUTINE RFLU_GENX_OpenRocinCtrlFiles







! ******************************************************************************
!
! Purpose: 
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_PutDataFlow(pRegion)
  
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
    
    CHARACTER(CHRLEN) :: fileNameMixt,fileNameGrid,fileStemMixt,fileStemGrid, &
                         matName,regIdString,winName
    CHARACTER(GENX_TIME_STRING_LEN) :: timeString
    INTEGER :: communicator,handleAttr,handleSetOption,handleAddAttr, & 
               handlePutAttr,nArgs,paneId,sz,ng
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

! ******************************************************************************
!   Set variables
! ******************************************************************************

    communicator = global%communicator

    winName = TRIM(global%volWinName)
    matName = 'fluid_vol'
    
    nArgs = 7

! ******************************************************************************
!   Build time and file strings
! ******************************************************************************

    CALL BuildRegionIdString(global,pRegion%iRegionGlobal,regIdString)
    CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,0,paneId)
    CALL RFLU_GENX_BuildTimeString(global%currentTime,timeString)

    CALL RFLU_GENX_GetFileStemGridVol(fileStemGrid)
    fileNameGrid = TRIM(fileStemGrid)//'_'//TRIM(timeString)//'_'// & 
                   TRIM(regIdString)                      

    CALL RFLU_GENX_GetFileStemMixtVol(fileStemMixt)
    fileNameMixt = TRIM(global%outDirHDF)//TRIM(fileStemMixt)//'_'// & 
                   TRIM(timeString)//'_'//TRIM(regIdString)         
        
! ******************************************************************************
!   Get function handles and set options
! ******************************************************************************

    handlePutAttr   = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.put_dataitem')
    handleAddAttr   = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.add_dataitem')    
    handleSetOption = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.set_option')

    CALL COM_call_function(handleSetOption,2,'rankwidth','0')
    ! MS
    !CALL COM_call_function(handleSetOption,2,'format','HDF')
    ! MS End

! ******************************************************************************
!   Write data through Rocout commands 
! ******************************************************************************
    
! ==============================================================================  
!   Conserved variables
! ==============================================================================  

    handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.rhof')
    CALL COM_get_size(TRIM(winName)//'.rhof', paneId, sz, ng)
    CALL COM_call_function(handlePutAttr,nArgs,TRIM(fileNameMixt), & 
                           handleAttr,TRIM(matName),timeString, & 
                           TRIM(fileNameGrid),communicator,paneId)
    
    handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.rhovf')
    CALL COM_call_function(handleAddAttr,nArgs,TRIM(fileNameMixt), & 
                           handleAttr,TRIM(matName),timeString, & 
                           TRIM(fileNameGrid),communicator,paneId)

    handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.rhoEf')
    CALL COM_call_function(handleAddAttr,nArgs,TRIM(fileNameMixt), & 
                           handleAttr,TRIM(matName),timeString, & 
                           TRIM(fileNameGrid),communicator,paneId)

! ==============================================================================  
!   Dependent variables
! ==============================================================================  

    handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.pf')
    CALL COM_call_function(handleAddAttr,nArgs,TRIM(fileNameMixt), & 
                           handleAttr,TRIM(matName),timeString, & 
                           TRIM(fileNameGrid),communicator,paneId)

    handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.Tf')
    CALL COM_call_function(handleAddAttr,nArgs,TRIM(fileNameMixt), & 
                           handleAttr,TRIM(matName),timeString, & 
                           TRIM(fileNameGrid),communicator,paneId)

    handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.af')
    CALL COM_call_function(handleAddAttr,nArgs,TRIM(fileNameMixt), & 
                           handleAttr,TRIM(matName),timeString, & 
                           TRIM(fileNameGrid),communicator,paneId)

! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_GENX_PutDataFlow









! ******************************************************************************
!
! Purpose: Write surface grid speeds through Rocout.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_PutDataGSpeedsSurf(pRegion)
  
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
    
    CHARACTER(CHRLEN) :: fileNameGrid,fileNameGSp,fileStemGrid,fileStemGSp, & 
                         matName,regIdString,winName
    CHARACTER(GENX_TIME_STRING_LEN) :: timeString
    INTEGER :: communicator,handleAttr,handleSetOption,handleAddAttr, & 
               handlePutAttr,handleWriteAttr,iPatch,nArgs,paneId
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global
    
    pGrid => pRegion%grid

    CALL RegisterFunction(global,'RFLU_GENX_PutDataGSpeedsSurf',&
  'RFLU_ModRocstarIO.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing surface grid speeds through Rocout...'
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid

! ******************************************************************************
!   Set variables
! ******************************************************************************

    communicator = global%communicator

    winName = TRIM(global%surfWinName)
    matName = 'fluid_surf'

    nArgs = 7

! ******************************************************************************
!   Build strings and file name 
! ******************************************************************************

    CALL BuildRegionIdString(global,pRegion%iRegionGlobal,regIdString)
    CALL RFLU_GENX_BuildTimeString(global%currentTime,timeString)

    CALL RFLU_GENX_GetFileStemGridSurf(fileStemGrid)    
    fileNameGrid = TRIM(fileStemGrid)//'_'//TRIM(timeString)//'_'// &
                   TRIM(regIdString)  

    CALL RFLU_GENX_GetFileStemGSpSurf(fileStemGSp)
    fileNameGSp = TRIM(global%outDirHDF)//TRIM(fileStemGSp)//'_'// & 
                  TRIM(timeString)//'_'//TRIM(regIdString)

! ******************************************************************************
!   Get function handles and set options
! ******************************************************************************

    handlePutAttr   = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.put_dataitem')
    handleAddAttr   = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.add_dataitem')    
    handleSetOption = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.set_option')

    CALL COM_call_function(handleSetOption,2,'rankwidth','0')
    ! MS
    !CALL COM_call_function(handleSetOption,2,'format','HDF')
    ! MS End

! ******************************************************************************
!   Write grid file through Rocout commands
! ******************************************************************************

! ==============================================================================  
!   Loop over patches
! ==============================================================================  

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,iPatch,paneId)

      IF ( iPatch == 1 ) THEN 
        handleWriteAttr = handlePutAttr
      ELSE 
        handleWriteAttr = handleAddAttr
      END IF ! iPatch

      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.gs')
      CALL COM_call_function(handleWriteAttr,nArgs,TRIM(fileNameGSp), &
                             handleAttr,TRIM(matName),timeString, &
                             TRIM(fileNameGrid),communicator,paneId)
    END DO ! iPatch
 
! ******************************************************************************
!   End
! ******************************************************************************
  
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
            'Writing surface grid speeds through Rocout done.'
    END IF ! global%myProcid  
  
    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_GENX_PutDataGSpeedsSurf









! ******************************************************************************
!
! Purpose: 
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_PutDataGSpeedsVol(pRegion)
  
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
    
    CHARACTER(CHRLEN) :: fileNameGrid,fileNameGSp,fileStemGrid,fileStemGSp, &
                         matName,regIdString,winName
    CHARACTER(GENX_TIME_STRING_LEN) :: timeString
    INTEGER :: communicator,handleAttr,handlePutAttr,handleSetOption,nArgs, & 
               paneId
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_GENX_PutDataGSpeedsVol',&
  'RFLU_ModRocstarIO.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing volume grid speeds through Rocout...'
    END IF ! global%myProcid
    
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid

! ******************************************************************************
!   Set variables
! ******************************************************************************

    communicator = global%communicator

    winName = TRIM(global%volWinName)
    matName = 'fluid_vol'

    nArgs = 7

! ******************************************************************************
!   Build time and file strings
! ******************************************************************************

    CALL BuildRegionIdString(global,pRegion%iRegionGlobal,regIdString)
    CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,0,paneId)
    CALL RFLU_GENX_BuildTimeString(global%currentTime,timeString)

    CALL RFLU_GENX_GetFileStemGridVol(fileStemGrid)
    fileNameGrid = TRIM(fileStemGrid)//'_'//TRIM(timeString)//'_'// & 
                   TRIM(regIdString)                      

    CALL RFLU_GENX_GetFileStemGSpVol(fileStemGSp)
    fileNameGSp = TRIM(global%outDirHDF)//TRIM(fileStemGSp)//'_'// & 
                  TRIM(timeString)//'_'//TRIM(regIdString)         
        
! ******************************************************************************
!   Get function handles and set options
! ******************************************************************************

    handlePutAttr   = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.put_dataitem')
    handleSetOption = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.set_option')

    CALL COM_call_function(handleSetOption,2,'rankwidth','0')
    ! MS
    !CALL COM_call_function(handleSetOption,2,'format','HDF')
    ! MS End

! ******************************************************************************
!   Write data through Rocout commands 
! ******************************************************************************
    
    handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.gs')
    CALL COM_call_function(handlePutAttr,nArgs,TRIM(fileNameGsp), & 
                           handleAttr,TRIM(matName),timeString, & 
                           TRIM(fileNameGrid),communicator,paneId)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
            'Writing volume grid speeds through Rocout done.'
    END IF ! global%myProcid  
    
    CALL DeregisterFunction(global)    
    
  END SUBROUTINE RFLU_GENX_PutDataGSpeedsVol








! ******************************************************************************
!
! Purpose: Write interface data through Roccom.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_PutDataInterf(pRegion)
  
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
    
    CHARACTER(CHRLEN) :: fileNameGrid,fileNameMixt,fileStemGrid,fileStemMixt, & 
                         matName,regIdString,winName
    CHARACTER(GENX_TIME_STRING_LEN) :: timeString
    INTEGER :: communicator,handleAddAttr,handleAttr,handlePutAttr, & 
               handleSetOption,handleWriteAttr,iPatch,nArgs,paneId
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing interface data through Roccom...'
    END IF ! global%verbLevel  

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid

! ******************************************************************************
!   Set variables
! ******************************************************************************

    communicator = global%communicator

    winName = TRIM(global%surfWinName)
    matName = 'fluid_surf'

    nArgs = 7

! ******************************************************************************
!   Build strings and file name 
! ******************************************************************************

    CALL BuildRegionIdString(global,pRegion%iRegionGlobal,regIdString)
    CALL RFLU_GENX_BuildTimeString(global%currentTime,timeString)

    CALL RFLU_GENX_GetFileStemGridSurf(fileStemGrid)    
    fileNameGrid = TRIM(fileStemGrid)//'_'//TRIM(timeString)//'_'// & 
                   TRIM(regIdString)  

    CALL RFLU_GENX_GetFileStemMixtSurf(fileStemMixt)
    fileNameMixt = TRIM(global%outDirHDF)//TRIM(fileStemMixt)//'_'// & 
                   TRIM(timeString)//'_'//TRIM(regIdString)

! ******************************************************************************
!   Get function handles and set options
! ******************************************************************************

    handlePutAttr   = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.put_dataitem')
    handleAddAttr   = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.add_dataitem')    
    handleSetOption = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.set_option')

    CALL COM_call_function(handleSetOption,2,'rankwidth','0')
    ! MS
    !CALL COM_call_function(handleSetOption,2,'format','HDF')
    ! MS End

! ******************************************************************************
!   Write data through Rocout commands
! ******************************************************************************

! ==============================================================================  
!   Loop over patches
! ==============================================================================  

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,iPatch,paneId)

      IF ( iPatch == 1 ) THEN 
        handleWriteAttr = handlePutAttr
      ELSE 
        handleWriteAttr = handleAddAttr
      END IF ! iPatch

! ------------------------------------------------------------------------------
!     Input data
! ------------------------------------------------------------------------------

      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.du_alp')
      CALL COM_call_function(handleWriteAttr,nArgs,TRIM(fileNameMixt), &
                             handleAttr,TRIM(matName),timeString, &
                             TRIM(fileNameGrid),communicator,paneId)
                             
      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.rhofvf_alp')
      CALL COM_call_function(handleWriteAttr,nArgs,TRIM(fileNameMixt), &
                             handleAttr,TRIM(matName),timeString, &
                             TRIM(fileNameGrid),communicator,paneId)
                             
      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.mdot_alp')
      CALL COM_call_function(handleWriteAttr,nArgs,TRIM(fileNameMixt), &
                             handleAttr,TRIM(matName),timeString, &
                             TRIM(fileNameGrid),communicator,paneId)
                             
      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.Tflm_alp')
      CALL COM_call_function(handleWriteAttr,nArgs,TRIM(fileNameMixt), &
                             handleAttr,TRIM(matName),timeString, &
                             TRIM(fileNameGrid),communicator,paneId)

! ------------------------------------------------------------------------------
!     Output data
! ------------------------------------------------------------------------------

      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.nfAlp')
      CALL COM_call_function(handleWriteAttr,nArgs,TRIM(fileNameMixt), &
                             handleAttr,TRIM(matName),timeString, &
                             TRIM(fileNameGrid),communicator,paneId)
                             
      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.rhofAlp')
      CALL COM_call_function(handleWriteAttr,nArgs,TRIM(fileNameMixt), &
                             handleAttr,TRIM(matName),timeString, &
                             TRIM(fileNameGrid),communicator,paneId)
                             
      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.pf')
      CALL COM_call_function(handleWriteAttr,nArgs,TRIM(fileNameMixt), &
                             handleAttr,TRIM(matName),timeString, &
                             TRIM(fileNameGrid),communicator,paneId)                             

      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.tf')
      CALL COM_call_function(handleWriteAttr,nArgs,TRIM(fileNameMixt), &
                             handleAttr,TRIM(matName),timeString, &
                             TRIM(fileNameGrid),communicator,paneId)                                                                                       

      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.qc')
      CALL COM_call_function(handleWriteAttr,nArgs,TRIM(fileNameMixt), &
                             handleAttr,TRIM(matName),timeString, &
                             TRIM(fileNameGrid),communicator,paneId)
                             
      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.qr')
      CALL COM_call_function(handleWriteAttr,nArgs,TRIM(fileNameMixt), &
                             handleAttr,TRIM(matName),timeString, &
                             TRIM(fileNameGrid),communicator,paneId)
                             
      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.tempf')
      CALL COM_call_function(handleWriteAttr,nArgs,TRIM(fileNameMixt), &
                             handleAttr,TRIM(matName),timeString, &
                             TRIM(fileNameGrid),communicator,paneId)
                             
      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.bFlag')
      CALL COM_call_function(handleWriteAttr,nArgs,TRIM(fileNameMixt), &
                             handleAttr,TRIM(matName),timeString, &
                             TRIM(fileNameGrid),communicator,paneId)                             
    END DO ! iPatch      

! ******************************************************************************
!   End
! ******************************************************************************
  
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing interface data through Roccom done.'
    END IF ! global%verbLevel  
  
  END SUBROUTINE RFLU_GENX_PutDataInterf








! ******************************************************************************
!
! Purpose: Write volume and surface grids through Rocout.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_PutGrid(pRegion)
  
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
    
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
        
    pGrid => pRegion%grid

! ******************************************************************************
!   Write volume and surface grid
! ******************************************************************************

    CALL RFLU_GENX_PutGridVol(pRegion)

    IF ( pGrid%nPatches > 0 ) THEN 
      CALL RFLU_GENX_PutGridSurf(pRegion)
    END IF ! pGrid%nPatches

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_PutGrid
  
  
  
  
  
  



! ******************************************************************************
!
! Purpose: Write surface grid through Rocout.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_PutGridSurf(pRegion)
  
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
    
    CHARACTER(CHRLEN) :: fileName,fileStem,matName,regIdString,winName
    CHARACTER(GENX_TIME_STRING_LEN) :: timeString
    INTEGER :: communicator,handleAttr,handleSetOption,handleAddAttr, & 
               handlePutAttr,handleWriteAttr,iPatch,nArgs,paneId, handleWrtAttr
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global
    
    pGrid => pRegion%grid

    CALL RegisterFunction(global,'RFLU_GENX_PutGridSurf',&
  'RFLU_ModRocstarIO.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing surface grid through Rocout...'
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid

! ******************************************************************************
!   Set variables
! ******************************************************************************

    communicator = global%communicator

    winName = TRIM(global%surfWinName)
    matName = 'fluid_surf'

    nArgs = 7

! ******************************************************************************
!   Build strings and file name 
! ******************************************************************************

    CALL BuildRegionIdString(global,pRegion%iRegionGlobal,regIdString)
    CALL RFLU_GENX_BuildTimeString(global%currentTime,timeString)
    CALL RFLU_GENX_GetFileStemGridSurf(fileStem)
    
    fileName = TRIM(global%outDirHDF)//TRIM(fileStem)//'_'// & 
               TRIM(timeString)//'_'//TRIM(regIdString)  

! ******************************************************************************
!   Get function handles and set options
! ******************************************************************************

    handlePutAttr   = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.put_dataitem')
    handleAddAttr   = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.add_dataitem')    
    handleWrtAttr   = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.write_dataitem')    
    handleSetOption = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.set_option')

    CALL COM_call_function(handleSetOption,2,'rankwidth','0')
    ! MS
    !CALL COM_call_function(handleSetOption,2,'format','HDF')
    ! MS End

! ******************************************************************************
!   Write grid file through Rocout commands
! ******************************************************************************

! ==============================================================================  
!   Loop over patches
! ==============================================================================  

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,iPatch,paneId)

! ------------------------------------------------------------------------------
!     Coordinates, predefined connectivity lists, and pconn
! ------------------------------------------------------------------------------

      IF ( iPatch == 1 ) THEN 
        handleWriteAttr = handlePutAttr
      ELSE 
        handleWriteAttr = handleAddAttr
      END IF ! iPatch
      !handleWriteAttr = handleWrtAttr

      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.pmesh')
      CALL COM_call_function(handleWriteAttr,nArgs,TRIM(fileName), &
                             handleAttr,TRIM(matName),timeString, &
                             '',communicator,paneId) 

! ------------------------------------------------------------------------------
!     Non-predefined connectivity lists. NOTE always need to be called, even if 
!     do not actually exist.
! ------------------------------------------------------------------------------

      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.t3g:real')
      CALL COM_call_function(handleAddAttr,nArgs,TRIM(fileName), &
                             handleAttr,TRIM(matName),timeString, &
                             '',communicator,paneId)
        
      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.t3g:virtual')
      CALL COM_call_function(handleAddAttr,nArgs,TRIM(fileName), &
                             handleAttr,TRIM(matName),timeString, &
                             '',communicator,paneId)

      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.q4g:real')
      CALL COM_call_function(handleAddAttr,nArgs,TRIM(fileName), &
                             handleAttr,TRIM(matName),timeString, &
                             '',communicator,paneId)

      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.q4g:virtual')
      CALL COM_call_function(handleAddAttr,nArgs,TRIM(fileName), &
                             handleAttr,TRIM(matName),timeString, &
                             '',communicator,paneId)

! ------------------------------------------------------------------------------
!     Patch quantities
! ------------------------------------------------------------------------------

      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.bcflag')
      CALL COM_call_function(handleAddAttr,nArgs,TRIM(fileName), &
                             handleAttr,TRIM(matName),timeString, &
                             '',communicator,paneId)

      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.patchNo')
      CALL COM_call_function(handleAddAttr,nArgs,TRIM(fileName), &
                             handleAttr,TRIM(matName),timeString, &
                             '',communicator,paneId)
                             
      handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.cnstr_type')
      CALL COM_call_function(handleAddAttr,nArgs,TRIM(fileName), &
                             handleAttr,TRIM(matName),timeString, &
                             '',communicator,paneId)                             
    END DO ! iPatch
 
! ******************************************************************************
!   End
! ******************************************************************************
  
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing surface grid through Rocout done.'
    END IF ! global%myProcid  
  
    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_GENX_PutGridSurf
  
  
  
  
  
  
  
  

! ******************************************************************************
!
! Purpose: Write volume grid through Rocout.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_PutGridVol(pRegion)
    
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
    
    CHARACTER(CHRLEN) :: fileName,fileStem,matName,regIdString,winName
    CHARACTER(GENX_TIME_STRING_LEN) :: timeString
    INTEGER :: communicator,handleAttr,handleSetOption,handleAddAttr, &
               handlePutAttr,nArgs,paneId
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start
! ******************************************************************************
    
    global => pRegion%global
    
    pGrid => pRegion%grid

    CALL RegisterFunction(global,'RFLU_GENX_PutGridVol',&
  'RFLU_ModRocstarIO.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing volume grid through Rocout...'
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid

! ******************************************************************************
!   Set variables
! ******************************************************************************

    communicator = global%communicator
    
    winName = TRIM(global%volWinName)
    matName = 'fluid_vol'
        
    nArgs = 7        
        
! ******************************************************************************
!   Build strings and file name 
! ******************************************************************************

    CALL BuildRegionIdString(global,pRegion%iRegionGlobal,regIdString)
    CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,0,paneId)
    CALL RFLU_GENX_BuildTimeString(global%currentTime,timeString)

    CALL RFLU_GENX_GetFileStemGridVol(fileStem)
    fileName = TRIM(global%outDirHDF)//TRIM(fileStem)//'_'// & 
               TRIM(timeString)//'_'//TRIM(regIdString)

! ******************************************************************************
!   Get function handles and set options
! ******************************************************************************

    handlePutAttr   = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.put_dataitem')
    handleAddAttr   = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.add_dataitem')    
    handleSetOption = COM_get_function_handle(TRIM(global%winNameOut)// & 
                                              '.set_option')

    CALL COM_call_function(handleSetOption,2,'rankwidth','0')
    ! MS
    !CALL COM_call_function(handleSetOption,2,'format','HDF')
    ! MS End

! ******************************************************************************
!   Write grid file through Rocout commands
! ******************************************************************************

! ==============================================================================  
!   Coordinates and connectivity
! ==============================================================================  
    handleAttr = COM_get_dataitem_handle(TRIM(winName)//'.pmesh')
    ! MS
    CALL COM_call_function(handlePutAttr,nArgs,TRIM(fileName),handleAttr, & 
                           TRIM(matName),timeString,'', &
                           communicator,paneId)
    ! MS End
    ! Original
    !CALL COM_call_function(handlePutAttr,nArgs,TRIM(fileName),handleAttr, & 
    !                       TRIM(matName),timeString,TRIM(fileName), &
    !                       communicator,paneId)
    ! Original End
                           
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing volume grid through Rocout done.'
    END IF ! global%myProcid

    CALL DeregisterFunction(global)
      
  END SUBROUTINE RFLU_GENX_PutGridVol  








! ******************************************************************************
!
! Purpose: Read window.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   windowType          Window type
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_ReadWindow(pRegion,windowType)
    
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    INTEGER, INTENT(IN) :: windowType
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: fileName,fileStem,matName,regIdString,winName, &
                         winNameIn
    CHARACTER(GENX_TIME_STRING_LEN) :: timeString
    INTEGER :: communicator,handleReadWin,paneId
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start
! ******************************************************************************
    
    global => pRegion%global
    
    pGrid => pRegion%grid

    CALL RegisterFunction(global,'RFLU_GENX_ReadWindow',&
  'RFLU_ModRocstarIO.F90')

! ******************************************************************************
!   Set variables
! ******************************************************************************

    communicator = global%communicator
                
! ******************************************************************************
!   Build strings and file name 
! ******************************************************************************

    CALL BuildRegionIdString(global,pRegion%iRegionGlobal,regIdString)
    CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,0,paneId)
    CALL RFLU_GENX_BuildTimeString(global%currentTime,timeString)

    IF ( windowType == GENX_WINDOW_TYPE_SURF ) THEN 
      winNameIn = global%surfWinNameInput
      winName   = global%surfWinName

      CALL RFLU_GENX_GetFileStemGridSurf(fileStem)
    ELSE IF ( windowType == GENX_WINDOW_TYPE_VOL ) THEN 
      winNameIn = global%volWinNameInput
      winName   = global%volWinName

      CALL RFLU_GENX_GetFileStemGridVol(fileStem)   
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! windowType

    fileName = TRIM(global%outDirHDF)//TRIM(fileStem)//'_'// & 
               TRIM(timeString)//'_'//TRIM(regIdString)//'*'

! ******************************************************************************
!   Get function handles and set options
! ******************************************************************************

    handleReadWin = COM_get_function_handle(TRIM(global%winNameIn)// & 
                                            '.read_window')

! ******************************************************************************
!   Read window
! ******************************************************************************

    CALL COM_call_function(handleReadWin,2,TRIM(fileName), & 
                           TRIM(winNameIn))
                           
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
      
  END SUBROUTINE RFLU_GENX_ReadWindow  









! ******************************************************************************
!
! Purpose: Write control file for Rocin.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   paneStringVol       String with pane ids for volume data
!   paneStringSurf      String with pane ids for surface data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_WriteRocinCtrlFiles(global,iProc,paneStringVol, & 
                                           paneStringSurf)
      
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    CHARACTER(*) :: paneStringSurf,paneStringVol
    INTEGER, INTENT(IN) :: iProc
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    CHARACTER(5) :: iRegString
    CHARACTER(CHRLEN) :: fileName,fileNames,fileStem
    CHARACTER(GENX_TIME_STRING_LEN) :: timeString       
    INTEGER :: iReg
    REAL(RFREAL) :: currentTime       
        
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_GENX_WriteRocinCtrlFiles', & 
                          'RFLU_ModRocstarIO.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing Rocin control files...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Get time string
! ******************************************************************************

    currentTime = 0.0_RFREAL ! Hardcoded, required by Rocin for correct starts
    
    CALL RFLU_GENX_BuildTimeString(currentTime,timeString)
          
! ******************************************************************************
!   Write file for volume data
! ******************************************************************************

! ==============================================================================
!   Processor id
! ==============================================================================

    WRITE(IF_CTRL_VOL,'(A,1X,I6)') '@Proc:',iProc-1       

! ==============================================================================
!   File names
! ==============================================================================

    fileNames = ''
    
    CALL RFLU_GENX_GetFileStemVol(fileStem)

    DO iReg = 1,global%nRegionsLocal 
      WRITE(iRegString,'(I5.5)') global%regMap(iReg)

      fileName  = TRIM(fileStem)//'*'//timeString//'_'//iRegString//'*'
      fileNames = TRIM(fileNames)//' '//TRIM(fileName)
    END DO ! iReg

    WRITE(IF_CTRL_VOL,'(A,1X,A)') '@Files:',TRIM(fileNames)

! ==============================================================================
!   Pane list
! ==============================================================================

    WRITE(IF_CTRL_VOL,'(A,1X,A)') '@Panes:',TRIM(paneStringVol)
 
! ******************************************************************************
!   Write file for surface data
! ******************************************************************************

! ==============================================================================
!   Processor id
! ==============================================================================

    WRITE(IF_CTRL_SURF,'(A,1X,I6)') '@Proc:',iProc-1       

! ==============================================================================
!   File names
! ==============================================================================

    fileNames = ''
    
    CALL RFLU_GENX_GetFileStemSurf(fileStem)

    DO iReg = 1,global%nRegionsLocal 
      WRITE(iRegString,'(I5.5)') global%regMap(iReg)

      fileName  = TRIM(fileStem)//'*'//timeString//'_'//iRegString//'*'
      fileNames = TRIM(fileNames)//' '//TRIM(fileName)
    END DO ! iReg

    WRITE(IF_CTRL_SURF,'(A,1X,A)') '@Files:',TRIM(fileNames)

! ==============================================================================
!   Pane list
! ==============================================================================

    WRITE(IF_CTRL_SURF,'(A,1X,A)') '@Panes:',TRIM(paneStringSurf)
 
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing Rocin control files done.'
    END IF ! global%verbLevel
      
    CALL DeregisterFunction(global)  
  
  END SUBROUTINE RFLU_GENX_WriteRocinCtrlFiles







! ******************************************************************************
!
! Purpose: Get Global data through Roccom.
!
! Description: None.
!
! Input:
!   pGlobal		Pointer to Global structure
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_GetGlobalData(Global)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_global), POINTER :: Global

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winName,winNameIn
    Double Precision :: ZoomLocal
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting Global data...'
    END IF ! global%verbLevel 
        

    winNameIn = global%surfWinNameInput
    winName   = global%surfWinName 

    CALL COM_set_size( TRIM(winName)//'.zoomFact',0,1)
    CALL COM_set_array( TRIM(winName)//'.zoomFact',0, global%Zoomfactor)

    RETURN
  END SUBROUTINE RFLU_GENX_GetGlobalData
    

END MODULE RFLU_ModRocstarIO

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModRocstarIO.F90,v $
! Revision 1.15  2008/12/06 08:44:21  mtcampbe
! Updated license.
!
! Revision 1.14  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.13  2007/04/14 14:12:37  mtcampbe
! Mods for TZ
!
! Revision 1.12  2006/11/07 18:13:12  mtcampbe
! Commented the line to obtain cntstr_type from Roccom
!
! Revision 1.11  2006/08/08 17:23:35  rfiedler
! Use MOVEDIR from *.bc to get cnstr_type, not the HDF values.
!
! Revision 1.10  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.9  2005/10/14 14:16:08  haselbac
! Added getting of Tb_alp
!
! Revision 1.8  2005/06/17 21:31:23  haselbac
! Bug fix in RFLU_GENX_GetGrid: Added IF statement when getting patch data
!
! Revision 1.7  2005/06/10 18:04:37  haselbac
! Bug fix for backward-compatibility in GENx with cnstr_type
!
! Revision 1.6  2005/06/09 20:20:14  haselbac
! Added cnstr_type to IO routines
!
! Revision 1.5  2005/04/15 15:06:56  haselbac
! Added routine to read window, cosmetics
!
! Revision 1.4  2004/12/10 15:33:46  haselbac
! Added reading/writing of patchNo dataitem
!
! Revision 1.3  2004/11/03 17:02:48  haselbac
! Removal of vertex and cell flag IO
!
! Revision 1.2  2004/10/22 14:01:29  haselbac
! Removed trailing dot to match new hdf file names
!
! Revision 1.1  2004/10/19 19:27:26  haselbac
! Initial revision
!
! ******************************************************************************

















