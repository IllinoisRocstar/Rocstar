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
! Purpose: Set values derived from user input.
!
! Description: Global and user input parameters are used to derive derived 
!              variables which are mostly components of periInput data type 
!              (region%periInput) as rocperi data structure does not have
!              its own global data type. Derivations for CPR are based on 1D 
!              analysis can be found in AIAA-86-1447 (Traineau, Hervat and 
!              Kuentsmann). Derivations for channel flow follow Dean`s and
!              theoretical correlations.
!
! Input: regions = Input parameters for all regions.
!
! Output: regions = Derived variables stored as periInput data.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PERI_DerivedInputValues.F90,v 1.5 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_DerivedInputValues( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPeriodic, ONLY   : t_peri_input
  USE ModBndPatch, ONLY   : t_patch
  USE ModError
  USE ModParameters
  USE PERI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER     :: global
  TYPE(t_peri_input), POINTER :: input
  TYPE(t_patch), POINTER      :: patches(:)

! ... cpr local variables
  REAL(RFREAL) :: rgas, gamma, gogampls, cprEps, headPres, headTemp, headSv 
  REAL(RFREAL) :: minj, denom, chokLen, delta, mfRat, phi, meanPgrad

! ... channel local variables
  REAL(RFREAL) :: ucUbulk, skinFriction, tauWall
  LOGICAL      :: pgTauWall

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_DerivedInputValues.F90,v $'

  global => regions(1)%global
  CALL RegisterFunction( global,'PERI_DerivedInputValues',&
  'PERI_DerivedInputValues.F90' )

  pgTauWall = .false.

#ifdef RFLO
  DO iReg=1,global%nRegions

    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

      patches => regions(iReg)%levels(1)%patches(:)
#endif
#ifdef RFLU
  DO iReg = LBOUND(regions,1),UBOUND(regions,1)

      patches => regions(iReg)%patches(:)
#endif

      input => regions(iReg)%periInput
      IF (input%flowKind==PERI_FLOW_CPR) THEN

! ----- get cpr head end pressure and temperature

        input%headPres = global%refPressure
        gamma          = global%refGamma
        rgas           = global%refCp*(1._RFREAL - 1._RFREAL/gamma)
        input%headTemp = global%refPressure/global%refDensity/rgas

! ----- derive initial mean pressure gradient from 1D analysis

        gogampls = gamma / (gamma + 1._RFREAL)
        cprEps   = input%cprEpsilon
        headPres = input%headPres
        headTemp = input%headTemp
        headSv   = SQRT( gamma*rgas*headTemp )
        IF (patches(3)%bcType >= BC_INJECTION .AND. &
            patches(3)%bcType <= BC_INJECTION+ &
            BC_RANGE) THEN
          minj = patches(3)%mixt%vals(BCDAT_INJECT_MFRATE,0)
        ELSE
          CALL ErrorStop( regions(iReg)%global,ERR_PERI_CPRBC,__LINE__, &
                          'cpr injection bc should be on patch 3 and 4' )
        ENDIF
        denom    = headSv*minj*SQRT( 2._RFREAL*(gamma+1._RFREAL) )
        choklen  = headPres*gamma/denom
        delta    = global%refLength
        mfRat    = 1._RFREAL/cprEps/chokLen
        phi      = SQRT(1._RFREAL - mfRat*mfRat)
        meanPgrad= -gogampls*mfRat/phi*headPres/(chokLen*delta)  ! dP/dx

! ----- cpr mean pressure gradient (in 'slow' coordinate)
        meanPgrad= meanPgrad/cprEps                              ! dP/dx_s

! ----- assign number of cpr variables and cpr mean values to permanent data str.

        input%nVar      = CPR_NVAR
        input%minjRate  = minj
        input%meanPgrad = meanPgrad
        input%bulkmFlux = 2._RFREAL*delta*minj/cprEps

      ELSEIF (input%flowKind==PERI_FLOW_CHANNEL) THEN

! ----- meaning global reference parameters
!       global%refDensity = wall density
!       global%refRenum   = bulk Reynolds number
!       global%refVelocity= estimate bulk velocity
!       global%refLength  = half channel height

        IF (global%refRenum > CNL_CRITREYN) THEN  !turbulent
! ------- uc/ubulk, from Dean`s correlation:
          ucUbulk = 1.28_RFREAL*(2._RFREAL*global%refRenum)**(-0.0116_RFREAL)
          skinFriction = 0.073_RFREAL*(2._RFREAL*global%refRenum)**(-0.25_RFREAL)
        ELSE  ! laminar
! ------- uc/ubulk, from laminar correlation:
          ucUbulk = 1.5_RFREAL
          skinFriction = 6._RFREAL/(global%refRenum*ucUbulk*ucUbulk)
        ENDIF

! ----- initial estimate physical parameters:
        input%cnlCvel = ucUbulk*global%refVelocity
        input%cnlUtau = global%refVelocity*(0.5_RFREAL*skinFriction)**0.5_RFREAL
        tauWall       = global%refDensity*input%cnlUtau**2._RFREAL
        input%cnlRetau= input%cnlUtau/global%refVelocity*global%refRenum

! ----- initial pressure gradient and reference mass flux:
        input%meanPgrad = -tauWall/global%refLength
        input%bulkmFlux = 2._RFREAL*global%refLength*global%refVelocity* &
                          global%refDensity

! ----- pressure gradient computation method:
        IF (input%pgradType == CNL_PGRAD_TAUWALL) pgTauWall = .TRUE.

      ELSEIF (input%flowKind==PERI_FLOW_BOLA) THEN

      ENDIF ! flowKind

#ifdef RFLO
    ENDIF  ! region active and my processor
#endif
  ENDDO   ! iReg

! global pressure gradient computation method:

#ifdef RFLO
  DO iReg=1,global%nRegions

    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
#endif
#ifdef RFLU
  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
#endif

      input => regions(iReg)%periInput
      IF (input%flowKind==PERI_FLOW_CHANNEL) THEN
        IF (pgTauWall) input%pgradType = CNL_PGRAD_TAUWALL
      ENDIF ! flowKind
#ifdef RFLO
    ENDIF   ! region active and my processor
#endif
  ENDDO     ! iReg

! Channel halfwidth delta and bulk mass flux are overruled in PERI_InitSolution
! with the actual values.

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PERI_DerivedInputValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_DerivedInputValues.F90,v $
! Revision 1.5  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/08/19 15:40:09  mparmar
! Renamed patch variables
!
! Revision 1.2  2004/06/11 21:49:40  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/06/08 23:56:56  wasistho
! changed nomenclature
!
! Revision 1.5  2003/10/14 21:26:18  wasistho
! replace Uc by Ub in uTau = f(Cf) formulation
!
! Revision 1.4  2003/09/18 01:57:13  wasistho
! added ijksplit and pgradType in PERI_PgradUpdate
!
! Revision 1.3  2003/04/03 00:33:01  wasistho
! correct initial mflux channel flow
!
! Revision 1.2  2003/04/02 01:48:48  wasistho
! minimize CPR user input
!
! Revision 1.1.1.1  2003/03/29 03:36:30  wasistho
! install ROCPERI
!
!
!
!******************************************************************************







