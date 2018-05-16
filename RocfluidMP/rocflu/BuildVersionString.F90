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
! Purpose: Build version string for printing in header.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   versionString       String containing version number and date
!
! Notes: 
!   1. The strings are NOT to be edited by anyone except the main code 
!      developer. 
!   2. Marks Rocbuild program will edit the build string to insert a 
!      build number or identifier
!
! ******************************************************************************
!
! $Id: BuildVersionString.F90,v 1.2562 2008/12/06 08:44:29 mtcampbe Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE BuildVersionString(versionString)

  USE ModDataTypes
  USE ModError
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
!*******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(CHRLEN), INTENT(OUT) :: versionString

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(LEN=2) :: major,minor,patch
  CHARACTER(LEN=3) :: build
  CHARACTER(CHRLEN) :: date,RCSIdentString

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: BuildVersionString.F90,v $ $Revision: 1.2562 $'

! ==============================================================================
! Set strings: DO NOT EDIT UNLESS YOU ARE MAIN DEVELOPER
! ==============================================================================

  major = '12'
  minor = '16'
  patch = '0'
  build = '1059' ! to be edited by Rocbuild

  date  = '04/04/07'

! ==============================================================================
! Write into version string
! ==============================================================================

  WRITE(versionString,'(A)') TRIM(major)//'.'//TRIM(minor)//'.'//TRIM(patch)
  WRITE(versionString,'(A)') 'Version: '//TRIM(versionString)//'-'//TRIM(build)
  WRITE(versionString,'(A)') TRIM(versionString)//', Date: '//TRIM(date)

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE BuildVersionString






