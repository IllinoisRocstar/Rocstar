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
! Purpose: Print dots on a line to indicate progress of a length task.
!
! Description: None.
!
! Input: 
!   nowValue    Current value
!   endValye    End value
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_PrintProgressDots.F90,v 1.6 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_PrintProgressDots(nowValue,endValue)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: global 
  USE ModParameters
  
  IMPLICIT NONE

! ... parameters
  INTEGER, INTENT(IN) :: endValue,nowValue

! ... loop variables


! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  REAL(RFREAL) :: pDone

!******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PrintProgressDots.F90,v $ $Revision: 1.6 $'

  CALL RegisterFunction('RFLU_PrintProgressDots',&
  'RFLU_PrintProgressDots.F90')

! start -----------------------------------------------------------------------

  pDone = 100.0_RFREAL*nowValue/REAL(endValue,KIND=RFREAL)
    
  IF ( INT(pDone) >= 10*global%progressCounter ) THEN 
    WRITE(STDOUT,'(A)',ADVANCE="NO") '.'
    global%progressCounter = global%progressCounter + 1
  END IF ! NINT
  
  IF ( nowValue == endValue ) THEN 
    global%progressCounter = 1
    WRITE(STDOUT,'(A)') ' '
  END IF ! nowValue

! end -------------------------------------------------------------------------

  CALL DeregisterFunction()

END SUBROUTINE 

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PrintProgressDots.F90,v $
! Revision 1.6  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2003/05/16 02:27:44  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.3  2003/03/15 17:06:38  haselbac
! Added KIND qualifyer
!
! Revision 1.2  2002/10/27 18:54:05  haselbac
! Removed tabs
!
! Revision 1.1  2002/07/25 14:34:59  haselbac
! Initial revision
!
!******************************************************************************







