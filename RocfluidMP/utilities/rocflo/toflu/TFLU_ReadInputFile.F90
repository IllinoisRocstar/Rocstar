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
! Purpose: read in user input from .inp file.
!
! Description: none.
!
! Input: regions
!
! Output: all regions: - input and output grid format
!                      - grid level
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TFLU_ReadInputFile.F90,v 1.4 2008/12/06 08:44:53 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadInputFile( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE TFLU_ModInterfaces, ONLY : ReadFormatsSection, ReadMultigridSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  CHARACTER(2*CHRLEN+4) :: fname
  CHARACTER(256)        :: line

  INTEGER :: errorFlag

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'ReadInputFile',&
  'TFLU_ReadInputFile.F90' )

! open file

  fname = TRIM(global%inDir)//TRIM(global%casename)//'.inp'
  OPEN(IF_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! read file looking for keywords

  DO
    READ(IF_INPUT,'(A256)',err=10,end=999) line
    SELECT CASE(TRIM(line))
    CASE ('# FORMATS')
      CALL ReadFormatsSection( global )

    CASE ('# MULTIGRID')
      CALL ReadMultigridSection( global )
    END SELECT
  ENDDO

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

END SUBROUTINE ReadInputFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TFLU_ReadInputFile.F90,v $
! Revision 1.4  2008/12/06 08:44:53  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/03 03:44:20  wasistho
! rflo_modinterfacestoflu to tflu_modinterfaces
!
! Revision 1.1  2004/12/03 02:59:30  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:58:20  wasistho
! lower to upper case
!
! Revision 1.1.1.1  2004/08/17 01:41:39  wasistho
! initial checkin
!
!
!******************************************************************************







