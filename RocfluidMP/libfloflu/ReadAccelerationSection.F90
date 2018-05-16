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
! Purpose: Read in user input related to acceleration terms.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReadAccelerationSection.F90,v 1.6 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReadAccelerationSection(global)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  
  USE ModInterfaces, ONLY: ReadSection  
  
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

  LOGICAL :: defined(5)
  CHARACTER(10) :: keys(5)
  REAL(RFREAL) :: vals(5)

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction( global,'ReadAccelerationSection',&
  'ReadAccelerationSection.F90' )

! ******************************************************************************
! Specify keywords and search for them
! ******************************************************************************

  keys(1) = 'TYPE'
  keys(2) = 'ACCELX'
  keys(3) = 'ACCELY'
  keys(4) = 'ACCELZ'
  keys(5) = 'GRAVITY'
  
  CALL ReadSection(global,IF_INPUT,5,keys,vals,defined)

! ******************************************************************************
! Set values
! ******************************************************************************

  IF ( defined(1) .EQV. .TRUE. ) THEN
    IF ( NINT(vals(1)) == 1 ) THEN 
      global%accelOn = .TRUE.
    ELSE 
      global%accelOn = .FALSE.
    END IF ! NINT(vals)
  ELSE 
    global%accelOn = .FALSE.
  END IF ! defined

  IF ( defined(2) .EQV. .TRUE. ) THEN
    global%accelX = vals(2)
  END IF ! defined
  
  IF ( defined(3) .EQV. .TRUE. ) THEN
    global%accelY = vals(3)
  END IF ! defined
  
  IF ( defined(4) .EQV. .TRUE. ) THEN
    global%accelZ = vals(4)
  END IF ! defined

  IF ( defined(5) .EQV. .TRUE. ) THEN
    global%gravity = vals(5)
  END IF ! defined

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE ReadAccelerationSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadAccelerationSection.F90,v $
! Revision 1.6  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/10/20 21:20:19  mparmar
! Added reading of gravity
!
! Revision 1.3  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.2  2005/11/10 01:55:42  haselbac
! Added Rocflu support, clean-up
!
! Revision 1.1  2004/12/01 16:50:09  haselbac
! Initial revision after changing case
!
! Revision 1.4  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.1  2003/08/28 20:05:38  jblazek
! Added acceleration terms.
!
! ******************************************************************************







