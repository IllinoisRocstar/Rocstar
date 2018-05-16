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
! Purpose: set values derived from user input for Eulerian particles.
!
! Description: none.
!
! Input: regions = input parameters for all regions.
!
! Output: regions = numerical parameters and no. of equations for each region.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_DerivedInputValues.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_DerivedInputValues( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartEul,    ONLY : t_peul_input, t_peul_ptype, t_peul
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iLev, iPtype

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_global),     POINTER :: global
  TYPE(t_peul_input), POINTER :: input
  TYPE(t_peul_ptype), POINTER :: ptype
  TYPE(t_peul),       POINTER :: peul

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_DerivedInputValues.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PEUL_DerivedInputValues',&
  'PEUL_DerivedInputValues.F90' )

! global values ---------------------------------------------------------------

! (none currently)

! region related values -------------------------------------------------------

  DO iReg=LBOUND(regions,1),UBOUND(regions,1)

    input => regions(iReg)%peulInput

! - number of variables

    DO iLev=1,regions(iReg)%nGridLevels

      peul => regions(iReg)%levels(iLev)%peul

      peul%nCv1 = CV_PEUL_NEQS
      peul%nDv1 = 0
      peul%nTv1 = 0

      peul%nCv  = input%nPtypes * peul%nCv1
      peul%nDv  = input%nPtypes * peul%nDv1
      peul%nTv  = input%nPtypes * peul%nTv1

    ENDDO ! iLev

! - Check peulInput quantities (none currently)

!-  Fill auxillary peulInput variables (none currently)

    DO iPtype = 1,input%nPtypes

      ptype => input%ptypes(iPtype)

! --- Check peulInput%ptypes quantities

      IF (.NOT. ASSOCIATED(ptype%material)) &
        CALL ErrorStop( global,ERR_VAL_UNDEFINED,__LINE__ )

      IF (ptype%material%molw <= 0._RFREAL) &
        CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )

      IF (ptype%material%molw > 10._RFREAL) THEN
        WRITE(STDOUT,*)
        WRITE(STDOUT,*) '### WARNING: Molecular weight suspiciously large: ', &
          ptype%material%molw
        WRITE(STDOUT,*) '###          E.g., in SI units carbon is 0.012'
        WRITE(STDOUT,*)
      ENDIF

      IF (ptype%material%dens <= 0._RFREAL) &
        CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )

      IF (ptype%material%spht <= 0._RFREAL) &
        CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )

      IF (ptype%diam <= 0._RFREAL) &
        CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )

      IF (ptype%diam > 0.099_RFREAL) THEN
        WRITE(STDOUT,*)
        WRITE(STDOUT,*) '### WARNING: Particle diameter suspiciously large: ', &
          ptype%diam
        WRITE(STDOUT,*) '###          Input units are meters'
        WRITE(STDOUT,*)
      ENDIF

      IF (ptype%puff < 1._RFREAL) THEN
        IF (ptype%puff < 0.999_RFREAL) THEN
          CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )
        ELSE
          ptype%puff = 1._RFREAL
        ENDIF
      ENDIF

! --- Sc not actually used
!      IF (ptype%Sc <= 0._RFREAL) &
!        CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )

! --- k2 not actually used
!      IF (ptype%vis2 < 0._RFREAL) &
!        CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )

      IF (ptype%vis4 <= 0._RFREAL) &
        CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )

      ptype%vis4 = 1._RFREAL/ptype%vis4

      SELECT CASE (ptype%negReport)
      CASE (PEUL_NEG_REPORT_NONE)
      CASE (PEUL_NEG_REPORT_USED)
      CASE DEFAULT
        CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )
      END SELECT ! ptype%negReport

      SELECT CASE (ptype%clipModel)
      CASE (PEUL_CLIP_MODEL_NONE)
      CASE (PEUL_CLIP_MODEL_USED)
      CASE DEFAULT
        CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )
      END SELECT ! ptype%clipModel

      SELECT CASE (ptype%methodV)

      CASE (PEUL_METHV_FLUIDVEL)

      CASE (PEUL_METHV_EQEUL)

! ----- require Navier-Stokes for Eq. Eul. because it uses the velocity
! ----- gradient information computed in ViscousFluxes

        IF (regions(iReg)%mixtInput%flowModel /= FLOW_NAVST) THEN
          WRITE(STDOUT,*) '### PEUL: Eq. Eul. method requires Navier-Stokes'
          CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )
        END IF ! flowModel

        regions(iReg)%mixtInput%computeTv = .TRUE.
        DO iLev=1,regions(iReg)%nGridLevels
          IF (regions(iReg)%levels(iLev)%mixt%nTv < 2) THEN
            regions(iReg)%levels(iLev)%mixt%nTv = 2
          END IF ! nTv
        END DO ! iLev
        IF (regions(iReg)%mixtInput%viscModel < VISC_SUTHR .OR.  &
            regions(iReg)%mixtInput%viscModel > VISC_ANTIB) THEN
          CALL ErrorStop(global,ERR_UNKNOWN_VISCMODEL,__LINE__)
        END IF ! viscModel

      CASE DEFAULT
        CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )

      END SELECT ! ptype%methodV

! --- Fill auxillary peulInput%ptypes variables

      ptype%denseff = ptype%material%dens / ptype%puff
      ptype%voleff  = ptype%diam**3 * global%pi / 6._RFREAL
      ptype%volmat  = ptype%voleff / ptype%puff

! --- Note: include pressure correction: (1 - beta) factor
      ptype%tauVcoef = (ptype%denseff - 1._RFREAL)*ptype%diam**2 / 18._RFREAL

      ptype%maxConc = 0.7_RFREAL*ptype%denseff

    ENDDO ! iPtype
  ENDDO  ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_DerivedInputValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_DerivedInputValues.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:32  haselbac
! Initial revision after changing case
!
! Revision 1.6  2004/08/02 16:38:40  jferry
! Added pressure correction to tau coefficient
!
! Revision 1.5  2004/07/30 22:47:36  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.4  2004/05/03 15:09:42  jferry
! added equilibrium Eulerian capability for smoke
!
! Revision 1.3  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.2  2003/03/11 16:04:57  jferry
! Created data type for material properties
!
! Revision 1.1  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







