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
! Purpose: compute interaction source for drag forces on Lagrangian particles.
!
! Description: none.
!
! Input: region  = current region.
!
! Output: region%levels(iLev)%plag%inrtSources
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_CalcDrag.F90,v 1.5 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_CalcDrag( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartLag,    ONLY : t_plag
  USE ModError
  USE ModParameters
  USE INRT_ModParameters

#ifdef PLAG
  USE PLAG_ModParameters
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! ... loop variables
  INTEGER :: iPcls

! ... local variables
  CHARACTER(CHRLEN)  :: RCSIdentString

  INTEGER :: dragModel, iCont, nCont, nPcls
#ifdef RFLO
  INTEGER :: iLev
#endif
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass

  REAL(RFREAL) :: diamL, massL, mixtVolR, pi, psiL, &
                  relVelMagL, reyL, tauLR
  REAL(RFREAL),          DIMENSION(3)   :: relVel, accelL
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv, pDv, pTv

  TYPE(t_plag)  , POINTER :: pPlag
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_CalcDrag.F90,v $ $Revision: 1.5 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_CalcDrag',&
  'INRT_CalcDrag.F90' )

#ifdef PLAG
! begin -----------------------------------------------------------------------

! Check if there are any particles

  nPcls = 0

#ifdef RFLO
  iLev  = region%currLevel
  IF (global%plagUsed) nPcls = region%levels(iLev)%plag%nPcls
#endif
#ifdef RFLU
  IF (global%plagUsed) nPcls = region%plag%nPcls
#endif

  IF (nPcls < 1) GO TO 999

! Get dimensions --------------------------------------------------------------

  pi = global%pi

  nCont     = region%plagInput%nCont
  dragModel = region%inrtInput%inrts(INRT_TYPE_DRAG)%switches( &
    INRT_SWI_DRAG_MODEL)

! Set pointers ----------------------------------------------------------------

#ifdef RFLO
  pPlag     => region%levels(iLev)%plag
#endif
#ifdef RFLU
  pPlag     => region%plag
#endif

  pCv       => pPlag%cv
  pDv       => pPlag%dv
  pTv       => pPlag%tv

  pCvPlagMass => pPlag%cvPlagMass

! Loop over all the particles -------------------------------------------------

  DO iPcls = 1,nPcls

    diamL = pDv(DV_PLAG_DIAM,iPcls)

    massL = SUM( pCv(pCvPlagMass(:),iPcls) )

    tauLR = 3.0_RFREAL*pi*pTv(TV_PLAG_MUELMIXT,iPcls)* &
            diamL/massL

    relVel(1)   = pDv(DV_PLAG_UVELMIXT,iPcls)-pDv(DV_PLAG_UVEL,iPcls)
    relVel(2)   = pDv(DV_PLAG_VVELMIXT,iPcls)-pDv(DV_PLAG_VVEL,iPcls)
    relVel(3)   = pDv(DV_PLAG_WVELMIXT,iPcls)-pDv(DV_PLAG_WVEL,iPcls)

    SELECT CASE (dragModel)

      CASE (INRT_DRAG_MODEL_STOKES)

        psiL = 1.0_RFREAL

      CASE (INRT_DRAG_MODEL_SN)

        relVelMagL  = SQRT( relVel(1)*relVel(1)+ &
                            relVel(2)*relVel(2)+ &
                            relVel(3)*relVel(3)  )

        reyL = diamL * relVelMagL * pDv(DV_PLAG_DENSMIXT,iPcls) / &
                                    pTv(TV_PLAG_MUELMIXT,iPcls)

        psiL = 1.0_RFREAL + 0.15_RFREAL* (reyL**0.687_RFREAL)

      CASE (INRT_DRAG_MODEL_SMRFLD)

        relVelMagL  = SQRT( relVel(1)*relVel(1)+ &
                            relVel(2)*relVel(2)+ &
                            relVel(3)*relVel(3)  )

        reyL = diamL * relVelMagL * pDv(DV_PLAG_DENSMIXT,iPcls) / &
                                    pTv(TV_PLAG_MUELMIXT,iPcls)

        psiL = 112.0_RFREAL/24.0_RFREAL *reyL**(+0.02_RFREAL)

      CASE DEFAULT
        CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

    END SELECT ! dragModel

    accelL(1) = psiL*relVel(1)*tauLR
    accelL(2) = psiL*relVel(2)*tauLR
    accelL(3) = psiL*relVel(3)*tauLR

! - fill the effect of the particles (L) on the gas (G).
! - the source terms will be added to the particle Rhs.

    pPlag%inrtSources(INRT_DRAG_L_XMOM_G,iPcls) = -massL*accelL(1)
    pPlag%inrtSources(INRT_DRAG_L_YMOM_G,iPcls) = -massL*accelL(2)
    pPlag%inrtSources(INRT_DRAG_L_ZMOM_G,iPcls) = -massL*accelL(3)

  END DO ! iPcls

! finalize --------------------------------------------------------------------

999  CONTINUE
#endif
  CALL DeregisterFunction( global )

END SUBROUTINE INRT_CalcDrag

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_CalcDrag.F90,v $
! Revision 1.5  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/03/08 14:59:36  fnajjar
! Fixed bug in Sommerfeld drag law for psiL
!
! Revision 1.2  2007/03/07 22:16:47  fnajjar
! Added Sommerfeld drag law
!
! Revision 1.1  2004/12/01 21:56:14  fnajjar
! Initial revision after changing case
!
! Revision 1.6  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.5  2004/01/31 03:59:22  haselbac
! Initial integration for Rocflu and Rocpart
!
! Revision 1.4  2003/04/03 21:10:17  jferry
! implemented additional safety checks for rocinteract
!
! Revision 1.3  2003/04/02 22:32:03  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.2  2003/03/11 16:09:39  jferry
! Added comments
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************







