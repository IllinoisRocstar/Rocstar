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
! Purpose: update solution/quantities pertinent to turbulence if desired
!
! Description: Depending user input parameters, different quantities are 
!              updated, such as vorticity components.
!
! Input: region = data of current region
!
! Output: specified TURB variables updated
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_SolutionUpdate.F90,v 1.7 2009/08/26 12:28:53 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_SolutionUpdate( region,istage,ibc,iec ) ! PUBLIC

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModTurbulence, ONLY : t_turb
  USE ModInterfaces,      ONLY : RkUpdateGeneric
  USE TURB_ModInterfaces, ONLY : TURB_CalcVortic, TURB_RansEmsUpdate
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
#ifdef RFLO
  TYPE(t_region) :: region
#endif
#ifdef RFLU
  TYPE(t_region), TARGET :: region
#endif
  INTEGER :: istage, ibc, iec

! ... local variables
  TYPE(t_global), POINTER :: global
  TYPE(t_turb),   POINTER :: turb
  INTEGER :: iLev, nCv
  LOGICAL :: compVort

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'TURB_SolutionUpdate',&
  'TURB_SolutionUpdate.F90' )

! get pointers and parameters

  nCv  =  region%turbInput%nCv 
#ifdef RFLO
  iLev =  region%currLevel
  turb => region%levels(iLev)%turb
#endif
#ifdef RFLU
  turb => region%turb
#endif

! update RaNS or DES solution -------------------------------------------------

  IF (region%turbInput%modelClass == MODEL_RANS) THEN
    IF (global%flowType == FLOW_UNSTEADY .AND. &
        global%solverType == SOLV_EXPLICIT) THEN
      CALL RkUpdateGeneric( region,VAR_TYPE_CELL,istage,ibc,iec,1,nCv, &
                            turb%cv,turb%cvOld,turb%rhs,turb%rhsSum )
    ELSE
      CALL TURB_RansEmsUpdate( region )
    ENDIF
  ENDIF

! update vorticities as post processing for LES -------------------------------

  IF (region%turbInput%modelClass == MODEL_LES) THEN

    compVort = .FALSE.

    IF (region%turbInput%calcVort == CALCVORT_FDT .AND. &
        region%irkStep == global%nrkSteps) THEN
      compVort = .TRUE.
    ENDIF
    IF (region%turbInput%calcVort == CALCVORT_SDT .AND. &
        region%irkStep == global%nrkSteps .AND. &
        (global%currentTime+global%dtMin) >= global%dTimeSystem) THEN
      compVort = .TRUE.
    ENDIF

    IF (compVort .eqv. .true.) CALL TURB_CalcVortic( region )
  ENDIF

! finalize -------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_SolutionUpdate

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_SolutionUpdate.F90,v $
! Revision 1.7  2009/08/26 12:28:53  mtcampbe
! Ported to Hera.   Fixed logical expression syntax errors.  Replaced all
! IF (logical_variable)  with IF (logical_variable .eqv. .true.) as
! consistent with the specification.  Also changed: IF( ASSOCIATED(expr) )
! to IF ( ASSOCIATED(expr) .eqv. .true. ).   Intel compilers produce code
! which silently fails for some mal-formed expressions, so these changes
! are a net which should ensure that they are evaluated as intended.
!
! Revision 1.6  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2004/12/04 06:08:42  wasistho
! removed unsteady restriction for vort.calculation
!
! Revision 1.3  2004/11/30 23:57:10  wasistho
! bugfix, split rflo/rflu for turb pointer
!
! Revision 1.2  2004/11/17 23:44:08  wasistho
! used generic RK-update for rocturb
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.5  2004/01/21 03:52:32  wasistho
! modify due to newly unsteady explicit and implicit (dual t-s) option
!
! Revision 1.4  2003/10/16 20:17:22  wasistho
! installed RaNS in steady state flow (Exp.Mult.Stg)
!
! Revision 1.3  2003/10/07 02:07:20  wasistho
! initial installation of RaNS-SA and DES
!
! Revision 1.2  2003/08/14 01:47:35  wasistho
! added vorticity call per system timestep
!
! Revision 1.1  2003/08/06 15:58:41  wasistho
! added vorticities computation
!
!
!
!******************************************************************************







