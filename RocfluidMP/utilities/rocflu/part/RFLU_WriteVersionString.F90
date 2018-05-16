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
!******************************************************************************
!
! Purpose: Write file with version string.
!
! Description: None.
!
! Input: 
!   global              Pointer to global type
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_WriteVersionString.F90,v 1.4 2008/12/06 08:45:03 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_WriteVersionString(global)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModMPI
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName,moduleString,RCSIdentString,versionString
  INTEGER :: errorFlag,versionWidth

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_WriteVersionString.F90,v $ $Revision: 1.4 $'

  CALL RegisterFunction(global,'RFLU_WriteVersionString',&
  'RFLU_WriteVersionString.F90')

! ******************************************************************************
! Build strings
! ******************************************************************************

  CALL BuildVersionString(versionString)

  moduleString = 'rflupart'
  iFileName    = TRIM(global%outDir)//'rflupart.vrs'  
  
! ******************************************************************************
! Open file, write string, close file
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC ) THEN    
    OPEN(IF_VRS,FILE=iFileName,FORM='FORMATTED',STATUS='UNKNOWN', &
         POSITION='APPEND',IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__)
    END IF ! global%error  
      
    WRITE(IF_VRS,'(A)') TRIM(moduleString)//" "//versionString
    
    CLOSE(IF_VRS,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__)
    END IF ! global%error                                    
  END IF ! global%myProcid

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_WriteVersionString

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_WriteVersionString.F90,v $
! Revision 1.4  2008/12/06 08:45:03  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/04/21 01:39:57  haselbac
! Changed file name, should have been changed before
!
! Revision 1.1  2005/04/15 15:09:20  haselbac
! Initial revision
!
! Revision 1.1  2003/03/31 16:06:31  haselbac
! Initial revision
!
! ******************************************************************************







