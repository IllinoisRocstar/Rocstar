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
! Purpose: Open file for total mass check.
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
! $Id: RFLU_OpenTotalMassFile.F90,v 1.5 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_OpenTotalMassFile(global)

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
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: fname,RCSIdentString
  INTEGER :: errorFlag

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_OpenTotalMassFile.F90,v $ $Revision: 1.5 $'

  CALL RegisterFunction(global,'RFLU_OpenTotalMassFile',&
  'RFLU_OpenTotalMassFile.F90')

! ==============================================================================
! Open file
! ==============================================================================

  IF ( global%myProcid == MASTERPROC ) THEN
    CALL BuildFileNamePlain(global,FILEDEST_OUTDIR,'.mass',fname)
  
! ==============================================================================
!   Append to existing file (restart) or create new file
! ==============================================================================

    IF ( ( global%flowType == FLOW_UNSTEADY .AND. & 
           global%currentTime > 0.0_RFREAL ) .OR. &
         ( global%flowType == FLOW_STEADY   .AND. & 
           global%currentIter>1 ) ) THEN
      OPEN(IF_MASS,FILE=fname,FORM='FORMATTED',STATUS='OLD', &
                   POSITION='APPEND',IOSTAT=errorFlag)                     
    ELSE
      OPEN(IF_MASS,FILE=fname,FORM='FORMATTED',STATUS='UNKNOWN', &
                   IOSTAT=errorFlag)
    END IF ! global
   
    global%error = errorFlag
    IF ( global%error /= 0 ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname))
    END IF ! global%error
  END IF ! global%myProcid

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_OpenTotalMassFile

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_OpenTotalMassFile.F90,v $
! Revision 1.5  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:43  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/06/16 20:01:14  haselbac
! Added use of ModBuildFileNames, cosmetics
!
! Revision 1.2  2003/01/28 14:45:46  haselbac
! Use common building of file name
!
! Revision 1.1  2002/11/08 21:55:20  haselbac
! Initial revision
!
! ******************************************************************************







