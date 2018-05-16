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
! Purpose: invoke breakup model case for Lagrangian particles.
!
! Description: none.
!
! Input: region  = current region.
!
! Output: region%levels(iLev)%plag%cv
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_CalcBreakup.F90,v 1.3 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CalcBreakup( region, iReg )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartLag,    ONLY : t_plag, t_plag_input
  USE ModError
  USE PLAG_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

  INTEGER, INTENT(IN) :: iReg

! ... loop variables
  INTEGER :: iCont, iPcls

! ... local variables
  CHARACTER(CHRLEN)  :: RCSIdentString

  INTEGER :: breakupModel, breakupWebSwi, nCont, nPcls
#ifdef RFLO
  INTEGER :: iLev
#endif
  INTEGER, POINTER, DIMENSION(:) :: pCvPlagMass, pDvPlagVolu

  REAL(RFREAL) :: breakupFac, breakupFacR, densG, diamL, diamLSplit, &
                  oneThird, pi, presG, relVelMagL, surfTensL,        &
                  surfTensSum, voluSum, voluSumR, weberL, weberCrit

  REAL(RFREAL),          DIMENSION(3)   :: relVel
  REAL(RFREAL), POINTER, DIMENSION(:)   :: pSurfTens
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv,pCv, pDv

  TYPE(t_plag),   POINTER :: pPlag
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CalcBreakup.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PLAG_CalcBreakup',&
  'PLAG_CalcBreakup.F90' )

! begin =======================================================================

! Check if there are any particles

#ifdef RFLO
  iLev  = region%currLevel
  nPcls = 0
  IF (global%plagUsed) nPcls = region%levels(iLev)%plag%nPcls
#endif
#ifdef RFLU
  nPcls = 0
  IF (global%plagUsed) nPcls = region%plag%nPcls
#endif
  IF (nPcls < 1) GO TO 999

! Set pointers ----------------------------------------------------------------

#ifdef RFLO
  pPlag => region%levels(iLev)%plag
#endif
#ifdef RFLU
  pPlag => region%plag
#endif

  pArv  => pPlag%arv
  pCv   => pPlag%cv
  pDv   => pPlag%dv

  pCvPlagMass => pPlag%cvPlagMass
  pDvPlagVolu => pPlag%dvPlagVolu
  pSurfTens   => region%plagInput%surftens

! Get dimensions --------------------------------------------------------------

  pi = global%pi
  oneThird = 1.0_RFREAL/3.0_RFREAL

  nCont  = region%plagInput%nCont

  breakupModel  = region%plagInput%breakupModel
  breakupWebSwi = region%plagInput%breakupWebSwi
  breakupFac    = region%plagInput%breakupFac
  breakupFacR   = 1.0_RFREAL/breakupFac

! Set appropriate coefficients pertinent to breakup Model ---------------------

  SELECT CASE (breakupModel)

    CASE (PLAG_BREAKUP_MODEL1)
      weberCrit = 10.0_RFREAL

    CASE DEFAULT
      CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

  END SELECT ! breakupModel

! Loop over all the particles -------------------------------------------------

  DO iPcls = 1,nPcls

! - Extract gas properties ----------------------------------------------------

    densG = pDv(DV_PLAG_DENSMIXT,iPcls)

! - Compute Particle surface tension ------------------------------------------

    voluSum  = SUM( pDv(pDvPlagVolu(:),iPcls) )
    voluSumR = 1.0_RFREAL/voluSum

    surfTensSum = SUM( pDv(pDvPlagVolu(:),iPcls) * pSurfTens(:) )

    surfTensL =  surfTensSum * voluSumR

    diamL = pDv(DV_PLAG_DIAM,iPcls)

! - Compute relative velocities and its magnitude ----------------------------

    relVel(1) = pDv(DV_PLAG_UVELMIXT,iPcls)-pDv(DV_PLAG_UVEL,iPcls)
    relVel(2) = pDv(DV_PLAG_VVELMIXT,iPcls)-pDv(DV_PLAG_VVEL,iPcls)
    relVel(3) = pDv(DV_PLAG_WVELMIXT,iPcls)-pDv(DV_PLAG_WVEL,iPcls)

    relVelMagL  = relVel(1)*relVel(1) &
                + relVel(2)*relVel(2) &
                + relVel(3)*relVel(3)

! - Compute weber number ------------------------------------------------------

    weberL = densG * diamL * relVelMagL / surfTensL

! - Check if critical weber number is met -------------------------------------

    IF ( weberL >= weberCrit ) THEN

#ifdef PLAG_DEBUG
      WRITE(*,'(A,3X,I3,3X,I4,3X,1PE12.5)') &
      'PLAG_CalcBreakup-Critical We reached: iReg, iPcls, We = ',&
       iReg, iPcls, weberL
#endif

! -- Redefine breakupFactor based on critical Weber, if needed ----------------

      IF ( breakupWebSwi == PLAG_BREAKUP_WEBSWI1 ) THEN
        breakupFac  = ( densG * diamL * relVelMagL /( surfTensL *weberCrit ) ) **3
        breakupFacR = 1.0_RFREAL/breakupFac

#ifdef PLAG_DEBUG
      WRITE(*,'(A,3X,I3,3X,I4,3X,1PE12.5)') &
      'PLAG_CalcBreakup-Breakup Switch Active: iReg, iPcls, breakupFac = ',&
       iReg, iPcls, breakupFac
#endif

      END IF ! breakupSwi

! -- Update cv and arv values --------------------------------------------------

      DO iCont = 1, nCont
       pCv(pCvPlagMass(iCont),iPcls) = pCv(pCvPlagMass(iCont),iPcls)*breakupFacR
      END DO ! iCont

      pCv(CV_PLAG_XMOM,iPcls) = pCv(CV_PLAG_XMOM,iPcls) * breakupFacR
      pCv(CV_PLAG_YMOM,iPcls) = pCv(CV_PLAG_YMOM,iPcls) * breakupFacR
      pCv(CV_PLAG_ZMOM,iPcls) = pCv(CV_PLAG_ZMOM,iPcls) * breakupFacR
      pCv(CV_PLAG_ENER,iPcls) = pCv(CV_PLAG_ENER,iPcls) * breakupFacR
      pCv(CV_PLAG_ENERVAPOR,iPcls) = pCv(CV_PLAG_ENERVAPOR,iPcls) * breakupFacR

      pArv(ARV_PLAG_SPLOAD,iPcls) = pArv(ARV_PLAG_SPLOAD,iPcls)* breakupFac
    END IF ! weberL

  END DO ! iPcls

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CalcBreakup

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CalcBreakup.F90,v $
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:00  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2004/03/25 21:16:43  jferry
! fixed Vapor Energy bug
!
! Revision 1.4  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.3  2004/02/06 21:18:20  fnajjar
! Initial Integration of Rocpart with Rocflu
!
! Revision 1.2  2003/09/15 20:26:54  fnajjar
! Corrected breakupFac and removed cubeRootFac
!
! Revision 1.1  2003/09/13 20:15:13  fnajjar
! Initialimport of breakup model
!
!******************************************************************************







