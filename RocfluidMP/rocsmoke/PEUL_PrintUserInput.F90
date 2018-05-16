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
! Purpose: write out user input for Eulerian particles for checking purposes.
!
! Description: none.
!
! Input: region = user input.
!
! Output: to standard output.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_PrintUserInput.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_PrintUserInput( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartEul,    ONLY : t_peul_ptype
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(IN) :: region

! ... loop variables
  INTEGER :: ipt

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_global),     POINTER :: global
  TYPE(t_peul_ptype), POINTER :: ptype

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_PrintUserInput.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_PrintUserInput',&
  'PEUL_PrintUserInput.F90' )

! begin -----------------------------------------------------------------------

  WRITE(STDOUT,*)
  WRITE(STDOUT,1010) SOLVER_NAME//'       Number of Smoke particle types', &
    region%peulInput%nPtypes
  WRITE(STDOUT,1020) SOLVER_NAME//'       Smoothing coefficient', &
    region%peulInput%smoocf

  DO ipt = 1,region%peulInput%nPtypes

    WRITE(STDOUT,*)

    IF (region%peulInput%nPtypes > 1) &
      WRITE(STDOUT,1010) SOLVER_NAME//'       *** Particle Type',ipt

    ptype => region%peulInput%ptypes(ipt)

    SELECT CASE (ptype%negReport)
    CASE (PEUL_NEG_REPORT_NONE)
      WRITE(STDOUT,1030) SOLVER_NAME//'         report negative values =  NO'
    CASE (PEUL_NEG_REPORT_USED)
      WRITE(STDOUT,1030) SOLVER_NAME//'         report negative values =  YES'
    CASE DEFAULT
      CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )
    END SELECT ! ptype%negReport

    SELECT CASE (ptype%clipModel)
    CASE (PEUL_CLIP_MODEL_NONE)
      WRITE(STDOUT,1030) SOLVER_NAME//'         clip   negative values =  NO'
    CASE (PEUL_CLIP_MODEL_USED)
      WRITE(STDOUT,1030) SOLVER_NAME//'         clip   negative values =  YES'
    CASE DEFAULT
      CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )
    END SELECT ! ptype%clipModel

    SELECT CASE (ptype%methodV)
    CASE (PEUL_METHV_FLUIDVEL)
      WRITE(STDOUT,1030) SOLVER_NAME//'         smoke method =  '// &
        'Fluid Velocity'
    CASE (PEUL_METHV_EQEUL)
      WRITE(STDOUT,1030) SOLVER_NAME//'         smoke method =  '// &
        'Equilibrium Eulerian Velocity'
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! ptype%methodV

    WRITE(STDOUT,1030) SOLVER_NAME//'         material =  '// &
      TRIM(ptype%material%name)
    WRITE(STDOUT,1021) SOLVER_NAME//'         diameter',ptype%diam
    WRITE(STDOUT,1021) SOLVER_NAME//'         puff fac',ptype%puff
    WRITE(STDOUT,1021) SOLVER_NAME//'         eff dens',ptype%denseff
    WRITE(STDOUT,1021) SOLVER_NAME//'         tauVcoef',ptype%tauVcoef

! - Sc not actually used
!    IF (region%mixtInput%flowModel /= FLOW_EULER) &
!      WRITE(STDOUT,1020) SOLVER_NAME//'         Sc      ',ptype%Sc

    IF (region%mixtInput%spaceDiscr == DISCR_CEN_SCAL) THEN
! - k2 not actually used
!      WRITE(STDOUT,1020) SOLVER_NAME//'         k2      ',ptype%vis2
      WRITE(STDOUT,1020) SOLVER_NAME//'         1/k4    ',1._RFREAL/ptype%vis4
    ENDIF

  ENDDO ! ipt

  WRITE(STDOUT,*)

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1000 FORMAT(/,A,1X,80('-'))
1005 FORMAT(/,A)
1010 FORMAT(A,' =',I3)
1015 FORMAT(A,' =',I9)
1020 FORMAT(A,' =',ES15.5)
1021 FORMAT(A,' =',EN15.5)
1025 FORMAT(/,A,' Region ',I6,':')
1030 FORMAT(A)
1035 FORMAT(/,A,1X,80('-'),/)

END SUBROUTINE PEUL_PrintUserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_PrintUserInput.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:38  haselbac
! Initial revision after changing case
!
! Revision 1.6  2004/05/03 15:09:42  jferry
! added equilibrium Eulerian capability for smoke
!
! Revision 1.5  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.4  2003/03/24 23:30:53  jferry
! overhauled rocinteract to allow interaction design to use user input
!
! Revision 1.3  2003/03/11 16:04:57  jferry
! Created data type for material properties
!
! Revision 1.2  2003/02/12 19:03:13  jferry
! corrected a comment
!
! Revision 1.1  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







