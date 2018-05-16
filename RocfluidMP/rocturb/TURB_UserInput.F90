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
! Purpose: Read TURB user input and store it in TURB data structure.
!
! Description: TURB simulation parameters are first set to default values.
!              They are then overruled by user input parameters. Based
!              on this, derived TURB variables are defined. Finally,
!              necessary checkings are performed against inconsistencies.
!
! Input: regions = user input in all regions
!
! Output: regions = TURB user input, parameters and derived variables
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_UserInput.F90,v 1.9 2009/08/26 12:28:53 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_UserInput( regions ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY      : t_region
  USE ModGlobal, ONLY          : t_global
  USE TURB_ModInterfaces, ONLY : TURB_CheckParamInput, &
                              TURB_DerivedInputValues, TURB_InitInputValues, &
                              TURB_ReadBcInputFile,    TURB_ReadInputFile
  USE ModTurbulence
  USE ModError
  USE ModMPI
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  INTEGER :: sndInteg, rcvInteg
  CHARACTER(2*CHRLEN+4) :: fname
  INTEGER :: errorFlag
!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_UserInput.F90,v $'

  global => regions(1)%global
  CALL RegisterFunction( global,'TURB_UserInput',&
  'TURB_UserInput.F90' )

! initialize parameters --------------------------------------------------

  CALL TURB_InitInputValues( regions )

! read user input



  CALL TURB_ReadInputFile( regions )
  
! global test for active turbulence and assign global values

#ifdef RFLO
  DO iReg = 1,global%nRegions
#endif
#ifdef RFLU
  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
#endif
    IF ((regions(iReg)%mixtInput%flowModel == FLOW_NAVST) .AND. &
        (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE)) THEN
      global%turbActive = .TRUE.
 
      IF ((regions(iReg)%mixtInput%turbModel == TURB_MODEL_SA) .OR. &
          (regions(iReg)%mixtInput%turbModel == TURB_MODEL_DESSA) .OR. &
          (regions(iReg)%mixtInput%turbModel == TURB_MODEL_HDESSA)) THEN
        global%turbCalcWDist = .TRUE.
      ENDIF
      global%calcCellCtr = global%turbCalcWDist
    ENDIF
  ENDDO  ! iReg

#ifdef RFLU
#ifdef MPI
! broadcast .true. global values 

!  sndInteg = 0
!  IF (global%turbActive .eqv. .true.) THEN
!    sndInteg = 1 
!    CALL MPI_ALLREDUCE( sndInteg,rcvInteg,1,MPI_INTEGER, MPI_SUM, &
!                        global%mpiComm, global%mpierr )
!  ENDIF
!  IF (rcvInteg > 0) global%turbActive = .TRUE.

!  sndInteg = 0
!  IF (global%turbCalcWDist .eqv. .true.) THEN
!    sndInteg = 1 
!    CALL MPI_ALLREDUCE( sndInteg,rcvInteg,1,MPI_INTEGER, MPI_SUM, &
!                        global%mpiComm, global%mpierr )
!  ENDIF
!  IF (rcvInteg > 0) global%turbCalcWDist = .TRUE.
!  global%calcCellCtr = global%turbCalcWDist
#endif
#endif

  IF (global%turbActive .eqv. .true.) THEN

! - read boundary conditions

    CALL TURB_ReadBcInputFile( regions )

! - set model & turbulence parameters from user input

    CALL TURB_DerivedInputValues( regions )

! - check user input and parameters set in the code

    CALL TURB_CheckParamInput( regions )

  ENDIF  ! turbActive

! finalize --------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_UserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_UserInput.F90,v $
! Revision 1.9  2009/08/26 12:28:53  mtcampbe
! Ported to Hera.   Fixed logical expression syntax errors.  Replaced all
! IF (logical_variable)  with IF (logical_variable .eqv. .true.) as
! consistent with the specification.  Also changed: IF( ASSOCIATED(expr) )
! to IF ( ASSOCIATED(expr) .eqv. .true. ).   Intel compilers produce code
! which silently fails for some mal-formed expressions, so these changes
! are a net which should ensure that they are evaluated as intended.
!
! Revision 1.8  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2008/10/23 18:20:57  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.5  2005/12/29 19:45:37  wasistho
! removed CHARM stuff
!
! Revision 1.4  2005/04/15 15:07:36  haselbac
! Removed Charm/FEM stuff
!
! Revision 1.3  2005/03/09 06:36:01  wasistho
! incorporated HDESSA
!
! Revision 1.2  2004/03/19 02:48:20  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.4  2003/10/07 02:07:48  wasistho
! initial installation of RaNS-SA and DES
!
! Revision 1.3  2003/08/01 22:17:20  wasistho
! prepared rocturb for Genx
!
! Revision 1.2  2003/05/31 01:48:23  wasistho
! installed turb. wall layer model
!
! Revision 1.1  2002/10/14 23:55:30  wasistho
! Install Rocturb
!
!
!******************************************************************************







