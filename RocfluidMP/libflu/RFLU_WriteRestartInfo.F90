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
!
! $Id: RFLU_WriteRestartInfo.F90,v 1.4 2008/12/06 08:44:13 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_WriteRestartInfo(global)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_CloseRestartInfo,RFLU_OpenRestartInfo

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
  CHARACTER(CHRLEN) :: RCSIdentString

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_WriteRestartInfo.F90,v $ $Revision: 1.4 $'

  CALL RegisterFunction(global,'RFLU_WriteRestartInfo',&
  'RFLU_WriteRestartInfo.F90')

! =============================================================================
! Write restart info to file
! =============================================================================

  IF ( global%myProcid == MASTERPROC ) THEN 
    CALL RFLU_OpenRestartInfo(global,FILE_POSITION_END,dummyLogical)

    IF ( global%flowType == FLOW_STEADY ) THEN 
      WRITE(IF_RESTINFO,*) global%currentIter
    ELSE 
      WRITE(IF_RESTINFO,*) global%currentTime
    END IF ! global%flowType

    CALL RFLU_CloseRestartInfo(global)
  END IF ! global%myProcid

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_WriteRestartInfo


!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_WriteRestartInfo.F90,v $
! Revision 1.4  2008/12/06 08:44:13  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.1  2003/06/20 22:32:30  haselbac
! Initial revision
!
!******************************************************************************







