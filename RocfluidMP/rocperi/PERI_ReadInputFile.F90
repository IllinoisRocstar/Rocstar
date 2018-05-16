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
! Purpose: call section readers pertinent to PERI (done on all processors).
!
! Description: none.
!
! Input: regions = user input file of all regions.
!
! Output: regions = PERI input parameters.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PERI_ReadInputFile.F90,v 1.4 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_ReadInputFile( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY      : t_region
  USE ModGlobal, ONLY          : t_global
  USE PERI_ModInterfaces, ONLY : PERI_ReadPeriSection
  USE ModError
  USE ModParameters
  
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  CHARACTER(2*CHRLEN+4) :: fname
  CHARACTER(256)        :: line
  INTEGER               :: errorFlag

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_ReadInputFile.F90,v $'

  global => regions(1)%global
  CALL RegisterFunction( global,'PERI_ReadInputFile',&
  'PERI_ReadInputFile.F90' )

! open file

  fname = TRIM(global%inDir)//TRIM(global%casename)//'.inp'
  OPEN(IF_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! read file looking for keywords

  DO
    READ(IF_INPUT,'(A256)',err=10,end=86) line

    SELECT CASE(TRIM(line))
    CASE ('# PERIODICFLOW')
      CALL PERI_ReadPeriSection( regions )
    END SELECT
  ENDDO

86  CONTINUE

! close file ------------------------------------------------------------------

  CLOSE(IF_INPUT,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )

  GOTO 999

! error handling

10  CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) )

999 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE PERI_ReadInputFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_ReadInputFile.F90,v $
! Revision 1.4  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/06/11 21:49:40  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/06/08 23:56:56  wasistho
! changed nomenclature
!
! Revision 1.1.1.1  2003/03/29 03:36:30  wasistho
! install ROCPERI
!
!
!******************************************************************************







