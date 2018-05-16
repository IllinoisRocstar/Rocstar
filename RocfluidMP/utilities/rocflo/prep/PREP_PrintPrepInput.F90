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
! Purpose: print on screen prep input for checking purposes.
!
! Description: none.
!
! Input: regions = user input.
!
! Output: to standard output.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PREP_PrintPrepInput.F90,v 1.3 2008/12/06 08:44:50 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PrintPrepInput( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'PrintPrepInput',&
  'PREP_PrintPrepInput.F90' )

! start -----------------------------------------------------------------------

  WRITE(STDOUT,1000) SOLVER_NAME

  WRITE(STDOUT,1005) SOLVER_NAME//' Actual rocflo-prep input selected:'
  WRITE(STDOUT,1030) SOLVER_NAME//' '
  WRITE(STDOUT,1015) SOLVER_NAME//'   Grid level            ',global%startLevel

  IF (global%gridFormat == FORMAT_ASCII) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   Input grid format      =  ASCII'
  ELSEIF (global%gridFormat == FORMAT_BINARY) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   Input grid format      =  binary'
  ELSEIF (global%gridFormat == FORMAT_HDF) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   Input grid format      =  HDF'
  ENDIF

  IF (global%solutFormat == FORMAT_ASCII) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   Output solution format =  ASCII'
  ELSEIF (global%solutFormat == FORMAT_BINARY) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   Output solution format =  binary'
  ELSEIF (global%solutFormat == FORMAT_HDF) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   Output solution format =  HDF'
  ENDIF

  IF (global%flowType == FLOW_STEADY) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   Flow type              =  steady'
  ELSE
    WRITE(STDOUT,1030) SOLVER_NAME//'   Flow type              =  unsteady'
  ENDIF

  WRITE(STDOUT,1020) SOLVER_NAME//'   Reference gamma       ',global%refGamma

! repeat for all regions ------------------------------------------------------

  IF (global%verbLevel < VERBOSE_HIGH) GOTO 888

  DO iReg=1,global%nRegions
    WRITE(STDOUT,1025) SOLVER_NAME,iReg

! - dimensions

    WRITE(STDOUT,1005) SOLVER_NAME//'     Dimensions:'
    WRITE(STDOUT,1015) SOLVER_NAME//'       icells      ',regions(iReg)%levels(1)%grid%ipc
    WRITE(STDOUT,1015) SOLVER_NAME//'       jcells      ',regions(iReg)%levels(1)%grid%jpc
    WRITE(STDOUT,1015) SOLVER_NAME//'       kcells      ',regions(iReg)%levels(1)%grid%kpc
    WRITE(STDOUT,1015) SOLVER_NAME//'       dummy cells ',regions(iReg)%nDumCells

! - initial solution

    WRITE(STDOUT,1005) SOLVER_NAME//'     Initial solution:'
    WRITE(STDOUT,1020) SOLVER_NAME//'       density     ',regions(iReg)%mixtInput%iniDens
    WRITE(STDOUT,1020) SOLVER_NAME//'       x-velocity  ',regions(iReg)%mixtInput%iniVelX
    WRITE(STDOUT,1020) SOLVER_NAME//'       y-velocity  ',regions(iReg)%mixtInput%iniVelY
    WRITE(STDOUT,1020) SOLVER_NAME//'       z-velocity  ',regions(iReg)%mixtInput%iniVelZ
    WRITE(STDOUT,1020) SOLVER_NAME//'       pressure    ',regions(iReg)%mixtInput%iniPress

  ENDDO   ! iReg

! finish ----------------------------------------------------------------------

888  CONTINUE

  WRITE(STDOUT,1035) SOLVER_NAME

  CALL DeregisterFunction( global )

1000 FORMAT(/,A,1X,80('-'))
1005 FORMAT(/,A)
1010 FORMAT(A,' = ',I2)
1015 FORMAT(A,' = ',I8)
1020 FORMAT(A,' = ',E12.5)
1025 FORMAT(/,A,'   Region ',I6,':')
1030 FORMAT(A)
1035 FORMAT(/,A,1X,80('-'),/)

END SUBROUTINE PrintPrepInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PREP_PrintPrepInput.F90,v $
! Revision 1.3  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/03 02:20:08  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:40:49  wasistho
! lower to upper case
!
! Revision 1.2  2004/07/28 01:47:08  wasistho
! cosmetics and cleanup
!
! Revision 1.1  2004/07/27 03:39:05  wasistho
! initial import genxWriteRocinout and printPrepInput
!
!
!******************************************************************************







