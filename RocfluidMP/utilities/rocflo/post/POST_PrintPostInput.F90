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
! Purpose: print on screen post input for checking purposes.
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
! $Id: POST_PrintPostInput.F90,v 1.3 2008/12/06 08:44:49 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PrintPostInput( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables

! ... local variables
  TYPE(t_global), POINTER :: global
  INTEGER :: iReg, iLev, l, h

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'PrintPostInput',&
  'POST_PrintPostInput.F90' )

! obtain parameters -----------------------------------------------------------

  iReg = region%iRegionGlobal
  iLev = region%currLevel

! start screen output ---------------------------------------------------------

  IF (iReg == 1) THEN

    WRITE(STDOUT,1000) SOLVER_NAME
    WRITE(STDOUT,1005) SOLVER_NAME//' Actual rocflo-post input selected:'
    WRITE(STDOUT,1030) SOLVER_NAME//' '
    WRITE(STDOUT,1015) SOLVER_NAME//'   Grid level            ',global%startLevel

    IF (global%flowType == FLOW_STEADY) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//'   Flow type              =  steady'
    ELSE
      WRITE(STDOUT,1030) SOLVER_NAME//'   Flow type              =  unsteady'
    ENDIF

    IF (global%postPlotType == PLOT_GRID_ONLY) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//'   Plot type              =  grid only'
    ELSE
      WRITE(STDOUT,1030) SOLVER_NAME//'   Plot type              =  grid+solution'
    ENDIF

    IF (global%flowType == FLOW_STEADY) THEN
      WRITE(STDOUT,1015) SOLVER_NAME//'   At iteration          ',global%postIter
    ELSE
      WRITE(STDOUT,1020) SOLVER_NAME//'   At time               ',global%postTime
    ENDIF

    IF (global%postOutFmt == PLOT_FMT_GENERIC) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//'   Output plot format     =  generic'
    ELSEIF (global%postOutFmt == PLOT_FMT_TECPLOT) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//'   Output plot format     =  Tecplot binary'
    ELSEIF (global%postOutFmt == PLOT_FMT_TECASCII) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//'   Output plot format     =  Tecplot ASCII'
    ENDIF

    WRITE(STDOUT,1030) SOLVER_NAME// &
      '   Gas properties         =  see Ref., Visc. and Material Section in .inp'

    IF ((global%postStatsFlag .EQV. .FALSE.) .AND. &
        (global%postTurbFlag  .EQV. .FALSE.) .AND. &
        (global%postPlagFlag  .EQV. .FALSE.) .AND. &
        (global%postRadiFlag  .EQV. .FALSE.) .AND. &
        (global%postSpecFlag  .EQV. .FALSE.)) THEN
      WRITE(STDOUT,1005) SOLVER_NAME//'   Additional modules postprocessed: none'
    ELSE
      WRITE(STDOUT,1005) SOLVER_NAME//'   Additional modules postprocessed:'
      IF (global%postStatsFlag) &
        WRITE(STDOUT,1030) SOLVER_NAME//'     Statistics solution'
      IF (global%postTurbFlag) &
        WRITE(STDOUT,1030) SOLVER_NAME//'     Turbulence solution'
      IF (global%postPlagFlag) &
        WRITE(STDOUT,1030) SOLVER_NAME//'     Particles solution'
      IF (global%postRadiFlag) &
        WRITE(STDOUT,1030) SOLVER_NAME//'     Radiation solution'
      IF (global%postSpecFlag) &
        WRITE(STDOUT,1030) SOLVER_NAME//'     Species solution'
    ENDIF ! module flag
  ENDIF   ! iReg=1

! grid and solution ----------------------------------------------------------

  IF (global%verbLevel < VERBOSE_HIGH) GOTO 888

  WRITE(STDOUT,1025) SOLVER_NAME,iReg

! dimensions

  WRITE(STDOUT,1005) SOLVER_NAME//'     Dimensions:'
  WRITE(STDOUT,1015) SOLVER_NAME//'       icells      ',region%levels(iLev)%grid%ipc
  WRITE(STDOUT,1015) SOLVER_NAME//'       jcells      ',region%levels(iLev)%grid%jpc
  WRITE(STDOUT,1015) SOLVER_NAME//'       kcells      ',region%levels(iLev)%grid%kpc
  WRITE(STDOUT,1015) SOLVER_NAME//'       dummy cells ',region%nDumCells

! minmax solution

  IF (global%postPlotType == PLOT_GRID_FLOW) THEN

    l = LBOUND( region%levels(iLev)%mixt%cv,2 )
    h = UBOUND( region%levels(iLev)%mixt%cv,2 )

    WRITE(STDOUT,1005) SOLVER_NAME//'     Minmax solution, dummy layers included:'
    WRITE(STDOUT,1022) SOLVER_NAME//'       density     ', &
                 MINVAL( region%levels(iLev)%mixt%cv(CV_MIXT_DENS,l:h) ), &
                 MAXVAL( region%levels(iLev)%mixt%cv(CV_MIXT_DENS,l:h) )
    WRITE(STDOUT,1022) SOLVER_NAME//'       x-velocity  ', &
                 MINVAL( region%levels(iLev)%mixt%dv(DV_MIXT_UVEL,l:h) ), &
                 MAXVAL( region%levels(iLev)%mixt%dv(DV_MIXT_UVEL,l:h) )
    WRITE(STDOUT,1022) SOLVER_NAME//'       y-velocity  ', &
                 MINVAL( region%levels(iLev)%mixt%dv(DV_MIXT_VVEL,l:h) ), &
                 MAXVAL( region%levels(iLev)%mixt%dv(DV_MIXT_VVEL,l:h) )
    WRITE(STDOUT,1022) SOLVER_NAME//'       z-velocity  ', &
                 MINVAL( region%levels(iLev)%mixt%dv(DV_MIXT_WVEL,l:h) ), &
                 MAXVAL( region%levels(iLev)%mixt%dv(DV_MIXT_WVEL,l:h) )
    WRITE(STDOUT,1022) SOLVER_NAME//'       pressure    ', &
                 MINVAL( region%levels(iLev)%mixt%dv(DV_MIXT_PRES,l:h) ), &
                 MAXVAL( region%levels(iLev)%mixt%dv(DV_MIXT_PRES,l:h) )
  ENDIF ! pltType

! finish ----------------------------------------------------------------------

888  CONTINUE

  IF (global%verblevel >= VERBOSE_HIGH .OR. iReg == 1) &
    WRITE(STDOUT,1035) SOLVER_NAME

  CALL DeregisterFunction( global )

1000 FORMAT(/,A,1X,80('-'))
1005 FORMAT(/,A)
1010 FORMAT(A,' = ',I2)
1015 FORMAT(A,' = ',I8)
1020 FORMAT(A,' = ',E12.5)
1022 FORMAT(A,' = ',2E12.5)
1025 FORMAT(/,A,'   Region ',I6,':')
1030 FORMAT(A)
1035 FORMAT(/,A,1X,80('-'))

END SUBROUTINE PrintPostInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: POST_PrintPostInput.F90,v $
! Revision 1.3  2008/12/06 08:44:49  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:59  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/03 02:03:16  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:32:01  wasistho
! lower to upper case
!
! Revision 1.2  2004/08/27 04:01:29  wasistho
! improved (cosmetics) screen output
!
! Revision 1.1  2004/07/28 01:51:28  wasistho
! initial import printPostInput
!
!
!******************************************************************************







