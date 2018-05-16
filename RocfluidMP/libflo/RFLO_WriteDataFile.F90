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
! Purpose: write real or integer data to a file (grid or solution).
!
! Description: file contains the following subroutines:
!
! - RFLO_WriteDataFileReal = real variable
! - RFLO_WriteDataFileInt  = integer variable
!
! Input: fileId  = I/O channel
!        form    = type of file (ASCII, binary, HDF)
!        nDim1/2 = 1st and 2nd dimension of var.
!        var     = real variable to read
!        ivar    = integer variable to read.
!
! Output: to file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_WriteDataFile.F90,v 1.3 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_WriteDataFileInt( global,fileId,form,nDim1,nDim2,ivar )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  INTEGER :: fileId, form, nDim1, nDim2

  INTEGER :: ivar(:,:)

  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: i1, i2

!******************************************************************************

  CALL RegisterFunction( global,'RFLO_WriteDataFileInt',&
  'RFLO_WriteDataFile.F90' )

! write

  IF (form == FORMAT_ASCII) THEN
    WRITE(fileId,*,err=10) ((ivar(i1,i2), i2=1,nDim2), i1=1,nDim1)
  ELSE IF (form == FORMAT_BINARY) THEN
    WRITE(fileId,err=10)   ((ivar(i1,i2), i2=1,nDim2), i1=1,nDim1)
  ELSE
    CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,&
    __LINE__ )
  ENDIF

! finalize, handle errors

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_WRITE,&
  __LINE__ )

999  CONTINUE
END SUBROUTINE RFLO_WriteDataFileInt

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_WriteDataFileReal( global,fileId,form,nDim1,nDim2,var )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  INTEGER :: fileId, form, nDim1, nDim2

  REAL(RFREAL) :: var(:,:)

  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER :: i1, i2

!******************************************************************************

  CALL RegisterFunction( global,'RFLO_WriteDataFileReal',&
  'RFLO_WriteDataFile.F90' )

! write

  IF (form == FORMAT_ASCII) THEN
    WRITE(fileId,*,err=10) ((var(i1,i2), i2=1,nDim2), i1=1,nDim1)
  ELSE IF (form == FORMAT_BINARY) THEN
    WRITE(fileId,err=10)   ((var(i1,i2), i2=1,nDim2), i1=1,nDim1)
  ELSE
    CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,&
    __LINE__ )
  ENDIF

! finalize, handle errors

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_WRITE,&
  __LINE__ )

999  CONTINUE
END SUBROUTINE RFLO_WriteDataFileReal

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_WriteDataFile.F90,v $
! Revision 1.3  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:21  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:17  wasistho
! lower to upper case
!
! Revision 1.4  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.3  2002/10/23 18:43:24  jblazek
! Changed temporary pointer arrays into allocatable arrays
! in grid and solution I/O routines.
!
! Revision 1.2  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/06/07 16:40:36  jblazek
! Grid & solution for all regions in one file.
!
!******************************************************************************








