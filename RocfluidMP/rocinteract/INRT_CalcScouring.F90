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
! Purpose: compute interaction sources on Lagrangian particles
!          for the scouring process.
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
! $Id: INRT_CalcScouring.F90,v 1.4 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_CalcScouring( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
#ifdef RFLO
#ifdef PEUL
  USE ModPartEul,    ONLY : t_peul
#endif
#endif
#ifdef RFLU
  USE ModSpecies,    ONLY : t_spec
#endif
  USE ModPartLag,    ONLY : t_plag
  USE ModMixture,    ONLY : t_mixt
  USE ModInteract
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
  INTEGER :: iEdge, iCont, iPcls

! ... local variables
  CHARACTER(CHRLEN)  :: RCSIdentString

  INTEGER :: iCell, indPeul0, iPeul, nCont, nEdges, nPcls
#ifdef RFLO
  INTEGER :: iLev
#endif
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass
  INTEGER, POINTER, DIMENSION(:,:) :: pAivL

  REAL(RFREAL) :: captureArea, coeffScour, diamL, mDotDepo, &
                  oneFourth, peulConc, pi, relVelMagL
  REAL(RFREAL),          DIMENSION(3)   :: velL, velS
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCvS, pDvL
#ifdef RFLO
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pDvG
#endif
#ifdef RFLU
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCvG
#endif

  TYPE(t_inrt_input),    POINTER :: pInputInrt
  TYPE(t_inrt_interact), POINTER :: pInrtScour
  TYPE(t_plag),          POINTER :: pPlag
#ifdef RFLO
#ifdef PEUL
  TYPE(t_peul),          POINTER :: pPeul
#endif
#endif
#ifdef RFLU
  TYPE(t_spec),          POINTER :: pPeul
#endif
  TYPE(t_mixt),          POINTER :: pMixt
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_CalcScouring.F90,v $ $Revision: 1.4 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_CalcScouring',&
  'INRT_CalcScouring.F90' )

#ifdef PLAG
! begin =======================================================================

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

#ifdef RFLU
! Check that have primitive state vector --------------------------------------

  IF ( region%mixt%cvState /= CV_MIXT_STATE_DUVWP ) THEN
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! region%mixt%cvState
#endif

! Set pointers ----------------------------------------------------------------

#ifdef RFLO
  pPlag  => region%levels(iLev)%plag
#ifdef PEUL
  pPeul  => region%levels(iLev)%peul
#endif
  pMixt  => region%levels(iLev)%mixt
#endif
#ifdef RFLU
  pPlag  => region%plag
  pPeul  => region%spec
  pMixt  => region%mixt
#endif

  pInputInrt => region%inrtInput

  pAivL  => pPlag%aiv
#ifdef RFLO
#ifdef PEUL
  pCvS   => pPeul%cv
#endif
#endif
#ifdef RFLU
  pCvS   => pPeul%cv
#endif
  pDvL   => pPlag%dv

#ifdef RFLO
  pDvG   => pMixt%dv
#endif
#ifdef RFLU
  pCvG   => pMixt%cv
#endif

  pInrtScour  => pInputInrt%inrts(INRT_TYPE_SCOURING)

! Get dimensions --------------------------------------------------------------

  oneFourth = 1.0_RFREAL/4.0_RFREAL
  pi = global%pi

  nCont  = region%plagInput%nCont
  nEdges = pInrtScour%nEdges

  indPeul0 = pInputInrt%indPeul0

! Loop over all the particles -------------------------------------------------

  DO iPcls = 1,nPcls

    diamL  = pDvL(DV_PLAG_DIAM,iPcls)

    velL(1)   = pDvL(DV_PLAG_UVEL,iPcls)
    velL(2)   = pDvL(DV_PLAG_VVEL,iPcls)
    velL(3)   = pDvL(DV_PLAG_WVEL,iPcls)

    iCell = pAivL(AIV_PLAG_ICELLS,iPcls)

! - the computation of relVelMagL will have to be moved inside the iEdge
! - loop when a smoke velocity not equal to fluid velocity is implemented.

#ifdef RFLO
    velS(1)   = pDvG(DV_MIXT_UVEL,iCell)
    velS(2)   = pDvG(DV_MIXT_VVEL,iCell)
    velS(3)   = pDvG(DV_MIXT_WVEL,iCell)
#endif
#ifdef RFLU
    velS(1)   = pCvG(CV_MIXT_XVEL,iCell)
    velS(2)   = pCvG(CV_MIXT_YVEL,iCell)
    velS(3)   = pCvG(CV_MIXT_ZVEL,iCell)
#endif

    relVelMagL  = SQRT( ( velL(1)-velS(1) )**2  &
                      + ( velL(2)-velS(2) )**2  &
                      + ( velL(3)-velS(3) )**2  )

    DO iEdge = 1,nEdges

      iPeul = pInrtScour%edges(iEdge)%iNode(1) - indPeul0
      peulConc = pCvS(iPeul,iCell)

!      Use definitions like these when smoke velocity comes to exist
!      velS(1)   = pDvS(DV_PEUL_UVEL,iCell)
!      velS(2)   = pDvS(DV_PEUL_VVEL,iCell)
!      velS(3)   = pDvS(DV_PEUL_WVEL,iCell)
!
!    relVelMagL  = SQRT( ( velL(1)-velS(1) )**2  &
!                      + ( velL(2)-velS(2) )**2  &
!                      + ( velL(3)-velS(3) )**2  )

      coeffScour = pInrtScour%data(INRT_DAT_SCOURING_COEF0 + iEdge)

      captureArea = oneFourth * pi * diamL**2 * coeffScour

      mDotDepo = peulConc * captureArea * relVelMagL

      pPlag%inrtSources(iEdge,iPcls) = mDotDepo

    ENDDO ! iEdge
  ENDDO ! iPcls

! finalize --------------------------------------------------------------------

999  CONTINUE
#endif
  CALL DeregisterFunction( global )

END SUBROUTINE INRT_CalcScouring

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_CalcScouring.F90,v $
! Revision 1.4  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/02/15 20:18:11  wasistho
! put peul within ifdef
!
! Revision 1.1  2004/12/01 21:56:16  fnajjar
! Initial revision after changing case
!
! Revision 1.7  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.6  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.5  2004/02/02 22:52:21  haselbac
! Fixed bug: Wrong parameter subscript
!
! Revision 1.4  2004/01/31 03:59:22  haselbac
! Initial integration for Rocflu and Rocpart
!
! Revision 1.3  2003/04/03 22:52:56  fnajjar
! Included correct pointer for InputInrt
!
! Revision 1.2  2003/04/03 21:10:18  jferry
! implemented additional safety checks for rocinteract
!
! Revision 1.1  2003/04/03 16:19:28  fnajjar
! Initial Import of routines for burning and scouring
!
!******************************************************************************







