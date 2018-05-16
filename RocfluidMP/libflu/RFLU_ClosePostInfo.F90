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
! Purpose: Close file with information for rflupost.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ClosePostInfo.F90,v 1.5 2008/12/06 08:44:11 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ClosePostInfo(global)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY : t_global
  USE ModMPI
  USE ModParameters

  USE ModBuildFileNames, ONLY: BuildFileNamePlain  
  
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

  CHARACTER(CHRLEN) :: iFileName,RCSIdentString
  INTEGER :: errorFlag

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ClosePostInfo.F90,v $ $Revision: 1.5 $'

  CALL RegisterFunction(global,'RFLU_ClosePostInfo',&
  'RFLU_ClosePostInfo.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing post-processor info file...'
  END IF ! global%verbLevel

! ******************************************************************************
! Close file
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC ) THEN
    CALL BuildFileNamePlain(global,FILEDEST_OUTDIR,'.inf',iFileName)

    CLOSE(IF_POSTINFO,IOSTAT=errorFlag)                        
    global%error = errorFlag
    IF ( global%error /= 0 ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error
  END IF ! global%myProcid

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                            'Closing post-processor info file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ClosePostInfo

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ClosePostInfo.F90,v $
! Revision 1.5  2008/12/06 08:44:11  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/06/16 20:00:15  haselbac
! Added use of ModBuildFileNames, cosmetics
!
! Revision 1.2  2003/06/20 22:33:33  haselbac
! Corrected comment
!
! Revision 1.1  2003/06/04 22:24:29  haselbac
! Initial revision
!
! ******************************************************************************







