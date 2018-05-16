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
! Purpose: read in user input related to boundary conditions pertinent to 
!          rocturb (done on all processors).
!
! Description: none.
!
! Input: boundary condition file.
!
! Output: regions = TURB-BC data for all regions.
!
! Notes: none
!
!******************************************************************************
!
! $Id: TURB_ReadBcInputFile.F90,v 1.6 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_ReadBcInputFile( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE TURB_ModInterfaces, ONLY : TURB_CoWlmReadBcSection
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  CHARACTER(CHRLEN)     :: RCSIdentString
  CHARACTER(2*CHRLEN+9) :: fname
  CHARACTER(256)        :: line
  TYPE(t_global), POINTER :: global

  INTEGER :: errorFlag

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_ReadBcInputFile.F90,v $'

  global => regions(1)%global

  CALL RegisterFunction( global,'TURB_ReadBcInputFile',&
  'TURB_ReadBcInputFile.F90' )

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Entering TURB_ReadBcInputFile...'
  END IF ! global%verbLevel

! open file

  WRITE(fname,'(A)') TRIM(global%inDir)//TRIM(global%casename)//'.bc'
  OPEN(IF_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! read file looking for keywords

  DO
    READ(IF_INPUT,'(A256)',err=10,end=99) line
    SELECT CASE(TRIM(line))

    CASE ('# BC_WALLMODEL')
      CALL TURB_CoWlmReadBcSection( regions )

    END SELECT
  ENDDO

99  CONTINUE
  CLOSE(IF_INPUT,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )

! finalization & error handling -----------------------------------------------

  GOTO 999

10  CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) )

999 CONTINUE

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Leaving TURB_ReadBcInputFile.'
  END IF ! global%verbLevel

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_ReadBcInputFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_ReadBcInputFile.F90,v $
! Revision 1.6  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/02/04 04:59:42  wasistho
! added enter and leave statements
!
! Revision 1.3  2004/03/13 03:15:43  wasistho
! get rid of flo/flu identifier in TURB_Co.. routines
!
! Revision 1.2  2004/03/08 23:31:05  wasistho
! changed turb nomenclature
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.1  2003/05/31 01:48:23  wasistho
! installed turb. wall layer model
!
!
!******************************************************************************







