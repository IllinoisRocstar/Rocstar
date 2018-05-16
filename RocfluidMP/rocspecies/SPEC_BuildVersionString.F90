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
!   versionString       String containing version number and date.
!
! Notes: 
!   1. The strings are NOT to be edited by anyone except the developer of 
!      this physical module. 
!
! ******************************************************************************
!
! $Id: SPEC_BuildVersionString.F90,v 1.23 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_BuildVersionString( versionString )

  USE ModDataTypes
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(*) :: versionString

! ==============================================================================
! Locals 
! ==============================================================================

  CHARACTER(LEN=2)  :: major, minor, patch
  CHARACTER(CHRLEN) :: date,RCSIdentString

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSFile$ $Revision: 1.23 $'

! ==============================================================================
! Set strings: DO NOT EDIT UNLESS YOU ARE ROCSPECIES DEVELOPER
! ==============================================================================

  major = '2'
  minor = '7'
  patch = '1'

  date  = '04/05/07'

! ==============================================================================
! Write into version string
! ==============================================================================

  WRITE(versionString,'(A)') TRIM(major)//'.'//TRIM(minor)//'.'//TRIM(patch)
  WRITE(versionString,'(A)') 'Version: '//TRIM(versionString)
  WRITE(versionString,'(A)') TRIM(versionString)//', Date: '//TRIM(date)

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE SPEC_BuildVersionString






