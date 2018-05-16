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
!
! $Id: RFLU_ReadRestartInfo.F90,v 1.6 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ReadRestartInfo(global)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModMPI
  
  USE ModInterfaces, ONLY: RFLU_CloseRestartInfo, & 
                           RFLU_OpenRestartInfo

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
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: dummyInteger,errorFlag
  REAL(RFREAL) :: dummyRFReal

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ReadRestartInfo.F90,v $ $Revision: 1.6 $'

  CALL RegisterFunction(global,'RFLU_ReadRestartInfo',&
  'RFLU_ReadRestartInfo.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Reading restart info file...'
  END IF ! global%verbLevel

! ==============================================================================
! Open file
! ==============================================================================

  CALL RFLU_OpenRestartInfo(global,FILE_POSITION_START,fileExists)

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

  CALL RFLU_CloseRestartInfo(global)        

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Reading restart info file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ReadRestartInfo


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ReadRestartInfo.F90,v $
! Revision 1.6  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.3  2004/06/25 20:04:51  haselbac
! Removed some code and put into RFLU_SetRestartTimeFlag
!
! Revision 1.2  2004/06/16 20:00:33  haselbac
! Added variables to make file handling easier, cosmetics
!
! Revision 1.1  2003/06/20 22:32:30  haselbac
! Initial revision
!
! ******************************************************************************







