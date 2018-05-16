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
! Purpose: read in user input from initialization file.
!
! Description: none.
!
! Input: regions
!
! Output: regions
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PREP_ReadInputFile.F90,v 1.4 2008/12/06 08:44:50 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadInputFile( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE PREP_ModInterfaces, ONLY : ReadFormatsSection, ReadReferenceSection, &
           ReadInitflowSection, ReadMultigridSection, ReadTimestepSection
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
  'PREP_ReadInputFile.F90' )

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

    CASE ('# REFERENCE')
      CALL ReadReferenceSection( global )

    CASE ('# INITFLOW')
      CALL ReadInitflowSection( regions )

    CASE ('# MULTIGRID')
      CALL ReadMultigridSection( global )

    CASE ('# TIMESTEP')
      CALL ReadTimestepSection( global )
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
! $Log: PREP_ReadInputFile.F90,v $
! Revision 1.4  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/03 03:30:19  wasistho
! rflo_modinterfacesprep to prep_modinterfaces
!
! Revision 1.1  2004/12/03 02:20:08  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:40:49  wasistho
! lower to upper case
!
! Revision 1.11  2004/07/23 04:32:31  wasistho
! Genx: readin from Rocin, standalone: read .inp file i.o. command line input
!
! Revision 1.10  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.9  2003/04/10 03:44:38  jblazek
! Merged .ini file into .inp file.
!
! Revision 1.8  2003/03/20 22:27:56  haselbac
! Renamed ModInterfaces
!
! Revision 1.7  2003/03/20 19:44:22  haselbac
! Corrected mistake in phased check-in
!
! Revision 1.6  2003/03/20 19:35:43  haselbac
! Modified RegFun call to avoid probs with long 'PREP_ReadInputFile.F90' names
!
! Revision 1.5  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.4  2002/09/20 22:22:37  jblazek
! Finalized integration into GenX.
!
! Revision 1.3  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/02/21 23:25:07  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/02 15:57:08  jblazek
! Added flow initialization.
!
!******************************************************************************








