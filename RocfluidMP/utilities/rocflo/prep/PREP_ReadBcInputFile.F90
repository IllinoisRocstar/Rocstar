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
! Purpose: read in user input related to boundary conditions
!          (done on all processors).
!
! Description: none.
!
! Input: boundary condition file.
!
! Output: regions = BC data for all regions.
!
! Notes: 
!
!******************************************************************************
!
! $Id: PREP_ReadBcInputFile.F90,v 1.5 2008/12/06 08:44:50 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadBcInputFile( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  CHARACTER(2*CHRLEN+9) :: fname
  CHARACTER(256)        :: line

  INTEGER :: errorFlag

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'ReadBcInputFile',&
  'PREP_ReadBcInputFile.F90' )

! open file

  WRITE(fname,'(A)') TRIM(global%inDir)//TRIM(global%casename)//'.bc'
  OPEN(IF_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! allocate variables

  ALLOCATE( global%prepBcDefined(BC_CODE_MAX) ,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! read file looking for keywords

  global%prepBcDefined(:) = .FALSE.
  
  DO
    READ(IF_INPUT,'(A256)',err=10,end=86) line
    SELECT CASE(TRIM(line))

    CASE ('# BC_SLIPW')
      global%prepBcDefined(BC_SLIPWALL) = .TRUE.

    CASE ('# BC_NOSLIP')
      global%prepBcDefined(BC_NOSLIPWALL) = .TRUE.

! Warning: keep this temporarily for backward compatibility
    CASE ('# BC_INFLOW')
      global%prepBcDefined(BC_INFLOW) = .TRUE.
! END Warning
    
    CASE ('# BC_INFLOW_TOTANG')
      global%prepBcDefined(BC_INFLOW_TOTANG) = .TRUE.
    
    CASE ('# BC_INFLOW_VELTEMP')
      global%prepBcDefined(BC_INFLOW_VELTEMP) = .TRUE.
    
    CASE ('# BC_INFLOW_VELPRESS')
      global%prepBcDefined(BC_INFLOW_VELPRESS) = .TRUE.

    CASE ('# BC_OUTFLOW')
      global%prepBcDefined(BC_OUTFLOW) = .TRUE.
    
    CASE ('# BC_FARF')
      global%prepBcDefined(BC_FARFIELD) = .TRUE.
    
    CASE ('# BC_INJECT')
      global%prepBcDefined(BC_INJECTION) = .TRUE.
    
    CASE ('# BC_INJECT_APN')
      global%prepBcDefined(BC_INJECTION) = .TRUE.

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

END SUBROUTINE ReadBcInputFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PREP_ReadBcInputFile.F90,v $
! Revision 1.5  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/01/29 09:34:50  wasistho
! added injection_apn
!
! Revision 1.2  2005/04/29 03:31:10  wasistho
! added distribution bc file generator
!
! Revision 1.1  2004/12/03 02:20:08  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:40:49  wasistho
! lower to upper case
!
! Revision 1.1  2004/07/27 20:29:47  wasistho
! added readBcInputFile and checkBcValidity
!
!
!******************************************************************************







