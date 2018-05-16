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
! Purpose: Suite for restart info routines
!
! Description: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModRestartInfo.F90,v 1.5 2008/12/06 08:44:17 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModRestartInfo

  USE ModGlobal, ONLY    : t_global 
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY      : t_grid
  USE ModBndPatch, ONLY  : t_patch
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLO_ReadRestartInfo, &
            RFLO_WriteRestartInfo

! private : RFLO_OpenRestartInfo
!           RFLO_CloseRestartInfo
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModRestartInfo.F90,v $ $Revision: 1.5 $'
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS

! ******************************************************************************
!
! Purpose: Open restart file.
!
! Description: None.
!
! Input: 
!   global              Pointer to global data
!   filePosition        Position at which file to be opened (start or end)
!
! Output: 
!   fileExists          Logical indicating whether file exists
!
! Notes: 
!   1. The filePosition parameter is needed because the restart info file
!      is opened with two goals. The first is to open it with the goal of
!      getting the last output iteration or time (ReadRestartInfo). 
!      The second is to be able to write to it by appending additional 
!      lines (WriteRestartInfo).
!   2. The fileExists parameter is needed because if the file exists
!      on restarting the code, it will need to be read to determine the 
!      last output iteration or time.
!
! ******************************************************************************

SUBROUTINE RFLO_OpenRestartInfo(global,filePosition,fileExists)
  
  USE ModBuildFileNames, ONLY: BuildFileNamePlain

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: filePosition
  TYPE(t_global), POINTER :: global

  LOGICAL, INTENT(OUT) :: fileExists

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: iFileName
  INTEGER :: errorFlag

! *****************************************************************************
! Start
! *****************************************************************************

  CALL RegisterFunction(global,'RFLO_OpenRestartInfo',&
       'RFLO_ModRestartInfo.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Opening restart info file...'
  END IF ! global%verbLevel

! =============================================================================
! Open file
! =============================================================================

  CALL BuildFileNamePlain(global,FILEDEST_OUTDIR,'.rin',iFileName)

  INQUIRE(FILE=iFileName,EXIST=fileExists)
    
  IF ( fileExists .EQV. .TRUE. ) THEN
    IF ( filePosition == FILE_POSITION_START ) THEN 
      OPEN(IF_RESTINFO,FILE=iFileName,FORM='FORMATTED',STATUS='OLD', &
           IOSTAT=errorFlag) 
    ELSE IF ( filePosition == FILE_POSITION_END ) THEN 
      OPEN(IF_RESTINFO,FILE=iFileName,FORM='FORMATTED',STATUS='OLD', &
           POSITION='APPEND',IOSTAT=errorFlag)     
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,&
           __LINE__)
    END IF ! filePosition
  ELSE
    OPEN(IF_RESTINFO,FILE=iFileName,FORM='FORMATTED',STATUS='NEW', &
         IOSTAT=errorFlag)    
  END IF ! file

  global%error = errorFlag
  IF ( global%error /= 0 ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,&
         __LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! *****************************************************************************
! End
! *****************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Opening restart info file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLO_OpenRestartInfo


! ******************************************************************************
!
! Purpose: Close restart info file.
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

SUBROUTINE RFLO_CloseRestartInfo(global)

  USE ModBuildFileNames, ONLY: BuildFileNamePlain  

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_global), POINTER :: global

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: iFileName
  INTEGER :: errorFlag
  
! *****************************************************************************
! Start
! *****************************************************************************

  CALL RegisterFunction(global,'RFLO_CloseRestartInfo',&
       'RFLO_ModRestartInfo.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing restart info file...'
  END IF ! global%verbLevel

! =============================================================================
! Close file
! =============================================================================

  CALL BuildFileNamePlain(global,FILEDEST_OUTDIR,'.rin',iFileName)

  CLOSE(IF_RESTINFO,IOSTAT=errorFlag)                        
  global%error = errorFlag
  IF ( global%error /= 0 ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,&
         __LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! *****************************************************************************
! End
! *****************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing restart info file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLO_CloseRestartInfo


! ******************************************************************************
!
! Purpose: Read restart info file to get last output iteration or time.
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

SUBROUTINE RFLO_ReadRestartInfo(global)

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

  LOGICAL :: fileExists
  INTEGER :: dummyInteger,errorFlag
  REAL(RFREAL) :: dummyRFReal

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLO_ReadRestartInfo',&
       'RFLO_ModRestartInfo.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Reading restart info file...'
  END IF ! global%verbLevel

! ==============================================================================
! Open file
! ==============================================================================

  CALL RFLO_OpenRestartInfo(global,FILE_POSITION_START,fileExists)

! ==============================================================================
! Set or read restart iteration or time
! ==============================================================================

  IF ( global%flowType == FLOW_STEADY ) THEN ! steady flow
    global%currentIter = 0

    IF ( fileExists .EQV. .TRUE. ) THEN
      DO 
        READ(IF_RESTINFO,*,IOSTAT=errorFlag) dummyInteger
        
        IF ( errorFlag /= ERR_NONE ) THEN 
          EXIT
        ELSE 
          global%currentIter = dummyInteger
        END IF ! errorFlag
      END DO ! <empty>
    END IF ! fileExists
  ELSE ! unsteady flow
    global%currentTime = 0.0_RFREAL
  
    IF ( fileExists .EQV. .TRUE. ) THEN
      DO 
        READ(IF_RESTINFO,*,IOSTAT=errorFlag) dummyRFReal
        
        IF ( errorFlag /= ERR_NONE ) THEN 
          EXIT
        ELSE 
          global%currentTime = dummyRFReal
          global%timeStamp   = global%currentTime
        END IF ! errorFlag
      END DO ! <empty>
    END IF ! fileExists
  END IF ! global%flowType

! ==============================================================================
! Write info
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    IF ( global%flowType == FLOW_STEADY ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME, & 
                                       'Restart iteration:',global%currentIter
    ELSE 
      WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME, & 
                                          'Restart time:',global%currentTime      
    END IF ! global%flowType
  END IF ! global%myProcid

! ==============================================================================
! Close file
! ==============================================================================

  CALL RFLO_CloseRestartInfo(global)        

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Reading restart info file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLO_ReadRestartInfo


!******************************************************************************
!
! Purpose: Write iteration or time to restart info file.
!
! Description: None.
!
! Input: 
!   global              Pointer to global data
!
! Output: None.
!
! Notes: 
!   1. Only master process writes to file. 
!   2. Every time write to file, open and close it to make sure always have
!      latest data in file.
!
!******************************************************************************

SUBROUTINE RFLO_WriteRestartInfo(global)

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_global), POINTER :: global

! =============================================================================
! Locals
! =============================================================================

  LOGICAL :: dummyLogical

! *****************************************************************************
! Start
! *****************************************************************************

  CALL RegisterFunction(global,'RFLO_WriteRestartInfo',&
       'RFLO_ModRestartInfo.F90')

! =============================================================================
! Write restart info to file
! =============================================================================

  IF ( global%myProcid == MASTERPROC ) THEN 
    CALL RFLO_OpenRestartInfo(global,FILE_POSITION_END,dummyLogical)

    IF ( global%flowType == FLOW_STEADY ) THEN 
      WRITE(IF_RESTINFO,*) global%currentIter
    ELSE 
      WRITE(IF_RESTINFO,*) global%currentTime
    END IF ! global%flowType

    CALL RFLO_CloseRestartInfo(global)
  END IF ! global%myProcid

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLO_WriteRestartInfo


! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModRestartInfo

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModRestartInfo.F90,v $
! Revision 1.5  2008/12/06 08:44:17  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:28  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/07 15:19:17  haselbac
! Removed tabs
!
! Revision 1.2  2006/02/04 03:59:11  wasistho
! global%timeStamp copied from currentTime in ReadRestartInfo
!
! Revision 1.1  2006/02/01 20:03:51  wasistho
! initial import
!
!
!
! ******************************************************************************










