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
! ******************************************************************************
!
! Purpose: Read in user input related to time-dependent boundary conditions
!   (done on all processors).
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Reads TBCs for all modules
!
! ******************************************************************************
!
! $Id: RFLU_ReadTbcInputFile.F90,v 1.5 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ReadTbcInputFile(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  USE ModMPI

  USE ModInterfaces, ONLY : RFLU_ReadTbcSection

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN+9) :: fname
  CHARACTER(256)      :: line

  INTEGER :: errorFlag,loopCounter

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction( global,'RFLU_ReadTbcInputFile',&
  'RFLU_ReadTbcInputFile.F90' )

! ******************************************************************************
! Open file
! ******************************************************************************

  WRITE(fname,'(A)') TRIM(global%inDir)//TRIM(global%casename)//'.bc'
  OPEN(IF_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! read file looking for keywords

  loopCounter = 0

  DO
    READ(IF_INPUT,'(A256)',ERR=10) line
    
    SELECT CASE( TRIM(line) )

    CASE ('# TBC_PIECEWISE')
      CALL RFLU_ReadTbcSection( pRegion, TBC_PIECEWISE )

    CASE ('# TBC_SINUSOIDAL')
      CALL RFLU_ReadTbcSection( pRegion, TBC_SINUSOIDAL )

    CASE ('# TBC_STOCHASTIC')
      CALL RFLU_ReadTbcSection( pRegion, TBC_STOCHASTIC )

    CASE ('# TBC_WHITENOISE')
      CALL RFLU_ReadTbcSection( pRegion, TBC_WHITENOISE )

    CASE ('# END')
      EXIT
    END SELECT ! TRIM(line)

    loopCounter = loopCounter + 1 ! Prevent infinite loop
    IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
      CALL ErrorStop(global,ERR_INFINITE_LOOP ,__LINE__)
    END IF ! loopCounter
  ENDDO ! <empty>

! close file

  CLOSE(IF_INPUT,IOSTAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname))
  END IF ! global%error

! finalization & error handling -----------------------------------------------

  GOTO 999

! error handling

10  CONTINUE
  CALL ErrorStop(global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname))

999 CONTINUE

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading Rocflu time-dependent '// & 
                             'boundary condition file done.'
  END IF ! global%verbLevel

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ReadTbcInputFile

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ReadTbcInputFile.F90,v $
! Revision 1.5  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/04/27 02:08:02  haselbac
! Cosmetics only
!
! Revision 1.2  2003/06/10 22:54:42  jferry
! Added Piecewise TBC
!
! Revision 1.1  2003/06/04 20:05:54  jferry
! re-worked implementation of TBCs in unstructured code
!
! ******************************************************************************







