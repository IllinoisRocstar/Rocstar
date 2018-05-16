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
!******************************************************************************
!
! $Id: BuildVersionString.F90,v 1.9 2008/12/06 08:44:54 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE BuildVersionString(versionString)

  USE ModDataTypes
  USE ModError
  
  IMPLICIT NONE

! ... parameters
  CHARACTER(CHRLEN), INTENT(OUT) :: versionString

! ... loop variables


! ... local variables
  CHARACTER(LEN=2) :: major,minor,patch
  CHARACTER(LEN=3) :: build
  CHARACTER(CHRLEN) :: date,RCSIdentString

!******************************************************************************

  RCSIdentString = '$RCSfile: BuildVersionString.F90,v $ $Revision: 1.9 $'

! Set strings: DO NOT EDIT UNLESS YOU ARE MAIN DEVELOPER ----------------------

  major = '0'
  minor = '1'
  patch = '6'
  build = '0' ! to be edited by Rocbuild

  date  = '05/15/03'

! Write into string -----------------------------------------------------------

  WRITE(versionString,'(A)') TRIM(major)//'.'//TRIM(minor)//'.'//TRIM(patch)
  WRITE(versionString,'(A)') 'Version: '//TRIM(versionString)//'-'//TRIM(build)  
  WRITE(versionString,'(A)') TRIM(versionString)//', Date: '//TRIM(date)

! Done ------------------------------------------------------------------------

END SUBROUTINE BuildVersionString

!******************************************************************************
!
! RCS Revision history:
!
! $Log: BuildVersionString.F90,v $
! Revision 1.9  2008/12/06 08:44:54  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:18:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2003/05/16 02:30:49  haselbac
! Updated patch level and date
!
! Revision 1.6  2003/04/17 00:17:47  haselbac
! Updated patch level and date
!
! Revision 1.5  2003/04/12 22:21:12  haselbac
! Updated patch level and date
!
! Revision 1.4  2003/04/10 14:51:15  haselbac
! Updated patch level and date
!
! Revision 1.3  2003/04/07 14:29:05  haselbac
! Updated patch level and date
!
! Revision 1.2  2003/04/02 17:35:02  haselbac
! Updated patch level and date
!
! Revision 1.1  2003/04/01 17:02:46  haselbac
! Initial revision
!
! Revision 1.3  2003/03/31 16:24:38  haselbac
! Updated patch level and date
!
! Revision 1.2  2003/03/20 20:45:16  haselbac
! Updated patch level and date
!
! Revision 1.1.1.1  2003/03/19 17:22:22  haselbac
! Initial revision
!
! Revision 1.35  2003/03/15 19:20:16  haselbac
! Updated major patch level and date
!
! Revision 1.34  2003/02/25 21:49:19  haselbac
! Changed patch level and date
!
! Revision 1.33  2003/02/20 19:50:47  haselbac
! Updated patch level and date
!
! Revision 1.32  2003/02/07 23:11:21  haselbac
! Updated patch level and date
!
! Revision 1.31  2003/02/06 19:33:23  haselbac
! Updated patch level and date
!
! Revision 1.30  2003/02/02 21:14:21  haselbac
! Updated patch level and date
!
! Revision 1.29  2003/02/01 00:30:35  haselbac
! Incremented patch level
!
! Revision 1.28  2003/01/31 13:54:54  haselbac
! Updated patch level and date
!
! Revision 1.27  2003/01/30 19:08:21  haselbac
! Updated patch level and date
!
! Revision 1.26  2003/01/28 21:00:31  haselbac
! Updated version number and date
!
! Revision 1.25  2003/01/09 15:32:37  haselbac
! Updated patch level and date
!
! Revision 1.24  2003/01/08 21:09:35  haselbac
! Changed patch level and date
!
! Revision 1.23  2002/12/03 23:02:10  haselbac
! Updated patch level and date
!
! Revision 1.22  2002/11/27 22:37:41  haselbac
! Updated patch level
!
! Revision 1.21  2002/11/27 20:46:46  haselbac
! Updated patch level and date
!
! Revision 1.20  2002/11/08 21:39:03  haselbac
! Updated patch level and date
!
! Revision 1.19  2002/10/27 19:26:38  haselbac
! Updated minor version number, patch level and date
!
! Revision 1.18  2002/10/21 14:34:58  haselbac
! Updated patch level and date
!
! Revision 1.17  2002/10/19 16:51:26  haselbac
! Updated patch level and date
!
! Revision 1.16  2002/10/17 20:11:34  haselbac
! Changed patch level and date
!
! Revision 1.15  2002/10/17 14:42:21  haselbac
! Updated patch level and date
!
! Revision 1.14  2002/10/12 15:30:56  haselbac
! Updated patch level and date
!
! Revision 1.13  2002/10/09 20:48:45  haselbac
! Updated patch level and date
!
! Revision 1.12  2002/10/08 21:52:11  haselbac
! Incremented patch level and date
!
! Revision 1.11  2002/10/07 14:25:10  haselbac
! Updated patch level and date
!
! Revision 1.10  2002/10/05 19:54:18  haselbac
! Incremented patch level and date
!
! Revision 1.9  2002/09/10 20:33:53  haselbac
! Incremented patch level and date
!
! Revision 1.8  2002/09/09 16:35:17  haselbac
! Changed patch number and date
!
! Revision 1.7  2002/07/25 16:22:37  haselbac
! Changed patch number and date
!
! Revision 1.6  2002/06/27 16:03:52  haselbac
! Changed patch number and date
!
! Revision 1.5  2002/06/17 21:52:24  haselbac
! Changed patch level and date
!
! Revision 1.4  2002/06/14 20:44:36  haselbac
! Changed patch level and date
!
! Revision 1.3  2002/06/10 21:42:33  haselbac
! Changed date
!
! Revision 1.2  2002/06/10 21:41:52  haselbac
! Incremented patch level
!
! Revision 1.1  2002/06/10 21:36:31  haselbac
! Initial revision
!
!******************************************************************************






