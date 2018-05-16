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
!   1. The strings are NOT to be edited by anyone except the main code 
!      developer. 
!   2. Marks Rocbuild program will edit the build string to insert a 
!      build number or identifier
!
!******************************************************************************
!
! $Id: POST_BuildVersionString.F90,v 1.3 2008/12/06 08:44:49 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE BuildVersionString( versionString )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  CHARACTER(*) :: versionString

! ... local variables
  CHARACTER(LEN=2)  :: major, minor, patch
  CHARACTER(LEN=3)  :: build
  CHARACTER(CHRLEN) :: date

!******************************************************************************
! set strings: DO NOT EDIT UNLESS YOU ARE MAIN DEVELOPER

  major = '2'
  minor = '3'
  patch = '0'
  build = '0' ! to be edited by Rocbuild

  date  = '09/04/03'

! write into string

  WRITE(versionString,'(A)') TRIM(major)//'.'//TRIM(minor)//'.'//TRIM(patch)
  WRITE(versionString,'(A)') 'Version: '//TRIM(versionString)//'-'//TRIM(build)
  WRITE(versionString,'(A)') TRIM(versionString)//', Date: '//TRIM(date)

END SUBROUTINE BuildVersionString

!******************************************************************************
!
! RCS Revision history:
!
! $Log: POST_BuildVersionString.F90,v $
! Revision 1.3  2008/12/06 08:44:49  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:59  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/03 02:30:43  wasistho
! added prefix
!
! Revision 1.10  2003/09/19 22:38:11  jblazek
! Added turbulence input.
!
! Revision 1.9  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.8  2002/10/19 00:40:31  jblazek
! Added utility (rflosurf) to write out surface grids for GenX.
!
! Revision 1.7  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.6  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.5  2002/08/07 21:07:55  jblazek
! Executables are left in the top-level directory instead of $(HOME)/bin.
!
! Revision 1.4  2002/07/20 00:42:05  jblazek
! Added ASCII Tecplot format.
!
! Revision 1.3  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.2  2002/06/22 01:13:38  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.1  2002/06/14 17:16:41  jblazek
! Added version string.
!
!******************************************************************************






