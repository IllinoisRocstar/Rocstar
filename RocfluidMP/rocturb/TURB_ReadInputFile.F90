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
! Purpose: Read in user input within TURB section (done on all processors).
!
! Description: none.
!
! Input: regions = user input file of all regions.
!
! Output: regions = turbulence input parameters.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_ReadInputFile.F90,v 1.7 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_ReadInputFile( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY      : t_region
  USE ModGlobal, ONLY          : t_global
  USE TURB_ModInterfaces, ONLY : TURB_ReadTurbSection
  USE ModTurbulence
  USE ModError
  USE ModMPI
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

  RCSIdentString = '$RCSfile: TURB_ReadInputFile.F90,v $'

  global => regions(1)%global
  CALL RegisterFunction( global,'TURB_ReadInputFile',&
  'TURB_ReadInputFile.F90' )

! open file -------------------------------------------------------------------

     fname = TRIM(global%inDir)//TRIM(global%casename)//'.inp'
     OPEN(IF_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag)
     global%error = errorFlag
     IF (global%error /= 0) THEN
        CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )
     ENDIF


! read file looking for keywords

  DO
    READ(IF_INPUT,'(A256)',err=10,end=86) line
        
    SELECT CASE(TRIM(line))
    CASE ('# TURBULENCE')
      CALL TURB_ReadTurbSection( regions )
    END SELECT
  ENDDO


86   CONTINUE

! close file ------------------------------------------------------------------

  CLOSE(IF_INPUT,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )

! finalization & error handling -----------------------------------------------

  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) )

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_ReadInputFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_ReadInputFile.F90,v $
! Revision 1.7  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2008/10/23 18:20:57  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.4  2006/04/07 15:19:25  haselbac
! Removed tabs
!
! Revision 1.3  2006/02/04 04:59:19  wasistho
! added enter and leave statements
!
! Revision 1.2  2004/03/19 02:48:02  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.4  2003/02/26 23:38:30  jferry
! eliminated end=999 convention to ensure that files get closed
!
! Revision 1.3  2002/10/16 02:23:35  wasistho
! Changed global%error flag
!
! Revision 1.2  2002/10/16 00:06:04  jblazek
! Added input directory to the file name.
!
! Revision 1.1  2002/10/14 23:55:30  wasistho
! Install Rocturb
!
!******************************************************************************







