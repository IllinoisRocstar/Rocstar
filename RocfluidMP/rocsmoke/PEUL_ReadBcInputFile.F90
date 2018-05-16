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
! Purpose: read in user input related to PEUL boundary conditions
!          (done on all processors).
!
! Description: none.
!
! Input: boundary condition file.
!
! Output: regions = BC data for all regions.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_ReadBcInputFile.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_ReadBcInputFile( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters

  USE PEUL_ModInterfaces, ONLY : PEUL_ReadBcInflowSection, &
                                 PEUL_ReadBcFarfSection, &
                                 PEUL_ReadBcInjectSection
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  CHARACTER(CHRLEN)     :: RCSIdentString
  CHARACTER(2*CHRLEN+9) :: fname
  CHARACTER(256)        :: line

  INTEGER :: errorFlag

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_ReadBcInputFile.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PEUL_ReadBcInputFile',&
  'PEUL_ReadBcInputFile.F90' )

! open file

  WRITE(fname,'(A)') TRIM(global%inDir)//TRIM(global%casename)//'.bc'
  OPEN(IF_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! read file looking for keywords

  DO
    READ(IF_INPUT,'(A256)',err=10,end=86) line

    SELECT CASE(TRIM(line))

    CASE ('# BC_PEUL_INFLOW')
      CALL PEUL_ReadBcInflowSection( regions )

    CASE ('# BC_PEUL_FARF')
      CALL PEUL_ReadBcFarfSection( regions )

    CASE ('# BC_PEUL_INJECT')
      CALL PEUL_ReadBcInjectSection( regions )

    END SELECT
  ENDDO

86   CONTINUE

! close file ------------------------------------------------------------------

  CLOSE(IF_INPUT,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )

! finalization & error handling -----------------------------------------------

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) )

999  CONTINUE

END SUBROUTINE PEUL_ReadBcInputFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_ReadBcInputFile.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:43  haselbac
! Initial revision after changing case
!
! Revision 1.3  2003/02/26 23:38:30  jferry
! eliminated end=999 convention to ensure that files get closed
!
! Revision 1.2  2003/02/12 23:34:48  jferry
! Replaced [io]stat=global%error with local errorFlag
!
! Revision 1.1  2003/02/11 22:52:51  jferry
! Initial import of Rocsmoke
!
!
!******************************************************************************







