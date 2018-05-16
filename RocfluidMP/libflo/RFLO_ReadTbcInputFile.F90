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
! Purpose: read in user input related to time-dependent boundary conditions
!          (done on all processors).
!
! Description: none.
!
! Input: boundary condition file.
!
! Output: regions = TBC data for all regions.
!
! Notes: reads TBCs for all modules
!
!******************************************************************************
!
! $Id: RFLO_ReadTbcInputFile.F90,v 1.3 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReadTbcInputFile( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_ReadTbcSection
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

  CALL RegisterFunction( global,'RFLO_ReadTbcInputFile',&
  'RFLO_ReadTbcInputFile.F90' )

! open file

  WRITE(fname,'(A)') TRIM(global%inDir)//TRIM(global%casename)//'.bc'
  OPEN(IF_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,&
    __LINE__,'File: '//TRIM(fname) )

! read file looking for keywords

  DO
    READ(IF_INPUT,'(A256)',err=10,end=86) line
    SELECT CASE(TRIM(line))

    CASE ('# TBC_PIECEWISE')
      CALL RFLO_ReadTbcSection( regions, TBC_PIECEWISE )

    CASE ('# TBC_SINUSOIDAL')
      CALL RFLO_ReadTbcSection( regions, TBC_SINUSOIDAL )

    CASE ('# TBC_STOCHASTIC')
      CALL RFLO_ReadTbcSection( regions, TBC_STOCHASTIC )

    CASE ('# TBC_WHITENOISE')
      CALL RFLO_ReadTbcSection( regions, TBC_WHITENOISE )

    END SELECT
  ENDDO

86   CONTINUE

! close file ------------------------------------------------------------------

  CLOSE(IF_INPUT,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_CLOSE,&
    __LINE__,'File: '//TRIM(fname) )

! finalization & error handling -----------------------------------------------

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
 CALL ErrorStop( global,ERR_FILE_READ,&
 __LINE__,'File: '//TRIM(fname) )

999  CONTINUE

END SUBROUTINE RFLO_ReadTbcInputFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReadTbcInputFile.F90,v $
! Revision 1.3  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:21  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:17  wasistho
! lower to upper case
!
! Revision 1.7  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.4  2003/06/10 22:54:42  jferry
! Added Piecewise TBC
!
! Revision 1.3  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.2  2003/02/26 23:38:30  jferry
! eliminated end=999 convention to ensure that files get closed
!
! Revision 1.1  2003/02/11 22:30:21  jferry
! Re-worked BC and TBC input routines to add multi-physics capability
!
!
!******************************************************************************







