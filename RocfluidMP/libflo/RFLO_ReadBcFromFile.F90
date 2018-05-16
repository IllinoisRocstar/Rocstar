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
! Purpose: read boundary condition data from a file.
!
! Description: none.
!
! Input: global = global variables (needed for error function)
!        fname  = file name
!        patch  = BC patch for which the data is to be read in.
!
! Output: patch%mixt%vals = BC data for the mixture.
!
! Notes: currently only the mixture BC data are read in.
!
!******************************************************************************
!
! $Id: RFLO_ReadBcFromFile.F90,v 1.6 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReadBcFromFile( global,fname,patch )

  USE ModDataTypes
  USE ModGlobal, ONLY   : t_global
  USE ModBndPatch, ONLY : t_patch
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  CHARACTER(*) :: fname

  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch

! ... loop variables
  INTEGER :: iReg, iPatch, n, i, j, ij

! ... local variables
  INTEGER :: n1, n2, iOff, nf1, nf2, errorFlag

!******************************************************************************

  CALL RegisterFunction( global,'RFLO_ReadBcFromFile',&
  'RFLO_ReadBcFromFile.F90' )

! dimensions

  n1   = ABS(patch%l1end-patch%l1beg)
  n2   = ABS(patch%l2end-patch%l2beg)
  iOff = n1 + 1

! read the file file

  WRITE(fname,'(A)') TRIM(global%inDir)//TRIM(fname)
  OPEN(IF_DISTR,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,&
    __LINE__,'File: '//TRIM(fname) )

  READ(IF_DISTR,*,err=10,end=10) nf1,nf2
  IF (nf1/=n1+1 .OR. nf2/=n2+1) &
    CALL ErrorStop( global,ERR_PATCH_DIMENS,&
    __LINE__ )

  DO n=1,patch%mixt%nData
    DO j=0,n2
      DO i=0,n1
        ij = IndIJ(i,j,iOff)
        READ(IF_DISTR,*,err=10,end=10) patch%mixt%vals(n,ij)
      ENDDO
    ENDDO
  ENDDO

  CLOSE(IF_DISTR,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_CLOSE,&
    __LINE__,'File: '//TRIM(fname) )

GOTO 999

! error handling

10  CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,&
  __LINE__,'File: '//TRIM(fname) )

999 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ReadBcFromFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReadBcFromFile.F90,v $
! Revision 1.6  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/08/19 15:38:09  mparmar
! Renamed patch variables
!
! Revision 1.3  2005/05/04 19:03:26  wasistho
! included dir-name in fname
!
! Revision 1.2  2005/05/03 03:20:52  wasistho
! added GOTO 999
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.2  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/19 00:40:30  jblazek
! Added utility (rflosurf) to write out surface grids for GenX.
!
! Revision 1.3  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.2  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.1  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
!******************************************************************************







