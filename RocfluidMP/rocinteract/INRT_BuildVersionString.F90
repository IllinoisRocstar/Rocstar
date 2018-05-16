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
! Purpose: Build version string for printing in header.
!
! Description: none.
!
! Input: none.
!
! Output: 
!   versionString = string containing version number and date.
!
! Notes: 
!   1. The strings are NOT to be edited by anyone except the developer of 
!      this physical module. 
!
!******************************************************************************
!
! $Id: INRT_BuildVersionString.F90,v 1.6 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_BuildVersionString( versionString )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  CHARACTER(*) :: versionString

! ... local variables
  CHARACTER(LEN=2)  :: major, minor, patch
  CHARACTER(CHRLEN) :: date

!******************************************************************************
! set strings: DO NOT EDIT UNLESS YOU ARE ROCINTERACT DEVELOPER

  major = '2'
  minor = '0'
  patch = '1'

  date  = '03/08/07'

! write into string

  WRITE(versionString,'(A)') TRIM(major)//'.'//TRIM(minor)//'.'//TRIM(patch)
  WRITE(versionString,'(A)') 'Version: '//TRIM(versionString)
  WRITE(versionString,'(A)') TRIM(versionString)//', Date: '//TRIM(date)

END SUBROUTINE INRT_BuildVersionString

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_BuildVersionString.F90,v $
! Revision 1.6  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2007/03/08 14:59:11  fnajjar
! Updated patch version and date
!
! Revision 1.3  2007/03/07 22:16:18  fnajjar
! Updated major version and date
!
! Revision 1.2  2004/12/01 22:02:03  fnajjar
! Updated minor version and date
!
! Revision 1.1  2004/12/01 00:04:32  wasistho
! added BuildVersionString
!
!
!******************************************************************************






