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
!          for the burning case.
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
! $Id: INRT_CalcBurning.F90,v 1.4 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_CalcBurning( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartLag,    ONLY : t_plag
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
  INTEGER :: iCont, iPcls, iPeulOutEdge

! ... local variables
  CHARACTER(CHRLEN)  :: RCSIdentString

  INTEGER :: burnModel, iContIn, iContOut, nCont, nEdges, &
             nPeulOutEdges, nPeulOxEdges, iPeulOx, nPcls, iCell
#ifdef RFLO
  INTEGER :: iLev
#endif
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv

  LOGICAL :: oxUsed, sendToVapor

  REAL(RFREAL) :: coeffHeatLatent, coeffLatentHeat, coeffPeul, densAl,  &
                  densAl2O3, densRatio, diamL, diffRel, dOxH2RdOxnoH2,  &
                  expDiamL, expHermsen, expPresG, expTempG, expXiEffG,  &
                  hCond, hEvap, hReac, hSolid, heatLatent, mDotBurn,    &
                  mFracL, molwAl, molwAl2O3, molwRatio, presG, tempG,   &
                  volFracL, xiCO2, xiEffG, xiO2, xiH2, xiH2O, massMixt, &
                  massOx, sendTemp

  REAL(RFREAL), POINTER, DIMENSION(:)   :: pDens, pMolw
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv, pDv, pCvMixt, pCvPeul

  TYPE(t_inrt_input),    POINTER :: pInputInrt
  TYPE(t_inrt_interact), POINTER :: pInrtBurn
  TYPE(t_plag),          POINTER :: pPlag
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_CalcBurning.F90,v $ $Revision: 1.4 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_CalcBurning',&
  'INRT_CalcBurning.F90' )

#ifdef PLAG
! begin =======================================================================

! Set grid-dependent pointers -------------------------------------------------

#ifdef RFLO
  iLev        =  region%currLevel
  pCvMixt     => region%levels(iLev)%mixt%cv
#ifdef PEUL
  pCvPeul     => region%levels(iLev)%peul%cv
#endif
  pPlag       => region%levels(iLev)%plag
#endif

#ifdef RFLU
  pCvMixt     => region%mixt%cv
  pCvPeul     => region%spec%cv
  pPlag       => region%plag
#endif

! Check if there are any particles

  nPcls = 0
  IF (global%plagUsed) nPcls = pPlag%nPcls

  IF (nPcls < 1) GO TO 999

! Set pointers ----------------------------------------------------------------

  pAiv        => pPlag%aiv
  pCv         => pPlag%cv
  pDv         => pPlag%dv

  pCvPlagMass => pPlag%cvPlagMass

  pDens       => region%plagInput%dens
  pMolw       => region%plagInput%molw

  pInputInrt  => region%inrtInput
  pInrtBurn   => pInputInrt%inrts(INRT_TYPE_BURNING)

! Get dimensions --------------------------------------------------------------

! Set parameters for oxidizer use

  nPeulOxEdges = 0
  oxUsed = (pInrtBurn%switches(INRT_SWI_BURNING_OXUSED) /= 0)

  IF (oxUsed) THEN

    nPeulOxEdges = 1
! - Extract index of oxidizer smoke type
    iPeulOx = pInrtBurn%edges(INRT_BURNING_S_MASS_X0 + &
      nPeulOxEdges)%iNode(1) - pInputInrt%indPeul0

  END IF ! oxUsed

! Set other parameters

  nCont = region%plagInput%nCont

  nEdges        = pInrtBurn%nEdges
  nPeulOutEdges = nEdges - nPeulOxEdges - INRT_BURNING_NEDGES0

  burnModel       = pInrtBurn%switches(INRT_SWI_BURNING_MODEL)
  coeffLatentHeat = pInrtBurn%data(INRT_DAT_BURNING_HEAT_COEF)

  mFracL = pInrtBurn%data(INRT_DAT_BURNING_MFRC_PLAG)

! Extract index of AL (In) and AL2O3 (Out). "In" material represents ----------
! the material that will burn and "out" is the burning product

  iContIn  = pInrtBurn%edges(INRT_BURNING_L_MASS_X)%iNode(1) &
           - pInputInrt%indPlag0
  iContOut = pInrtBurn%edges(INRT_BURNING_X_MASS_L + nPeulOxEdges)%iNode(2) &
           - pInputInrt%indPlag0

  densAl    = pDens(iContIn)
  densAl2O3 = pDens(iContOut)
  densRatio = densAL/densAl2O3

  molwAl    = pMolw(iContIn)
  molwAl2O3 = pMolw(iContOut)
  molwRatio = molwAl2O3/(2.0_RFREAL*molwAl)

! Set appropriate coefficients pertinent to burn Model ------------------------

  SELECT CASE (burnModel)

    CASE (INRT_BURNING_MODEL_BECKSTEAD)

! --- Exponents from Hermsens Model -------------------------------------------

      expHermsen = 1.9_RFREAL
      expDiamL   = 3.0_RFREAL-expHermsen

! --- Exponents for Burning Rate Variables ------------------------------------

      expTempG  = 1.57_RFREAL
      expPresG  = 0.20_RFREAL
      expXiEffG = 0.39_RFREAL

! --- Compute Measure of Availability of Oxidizing Species --------------------

      xiCO2 = 0.20_RFREAL
      xiO2  = 0.02_RFREAL
      xiH2O = 0.20_RFREAL
      xiH2  = 0.20_RFREAL

! --- Will override this value of xiEffG if an oxidizer species is used

      xiEffG = ( xiO2 + 0.58_RFREAL * xiH2O + 0.22_RFREAL * xiCO2 )

! --- Compute Relative Diffusivity Coefficient, D_rel -------------------------

! --- Define Ratio of Diffusivity of Oxidants with and without H2 -------------

      dOxH2RdOxnoH2 = 3.7_RFREAL

      diffRel = 1.0_RFREAL + xiH2 * ( dOxH2RdOxnoH2 - 1.0_RFREAL )

! --- Latent Heat for Energy Source -------------------------------------------

      hEvap  = 10896.0_RFREAL *1000.0_RFREAL
      hReac  =  9543.0_RFREAL *1000.0_RFREAL
      hCond  = 29767.0_RFREAL *1000.0_RFREAL
      hSolid =     0.0_RFREAL

      heatLatent = -hEvap + hReac + hCond + hSolid

! --- Scale Latent heat by appropriate coefficient ----------------------------

      heatLatent = heatLatent*coeffLatentHeat

    CASE DEFAULT
      CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

  END SELECT ! burnModel

  sendToVapor = (pInrtBurn%switches(INRT_SWI_BURNING_VAPOR_METH) /= &
                 INRT_BURNING_VAPOR_METH_NONE)

  sendTemp    = pInrtBurn%data(INRT_DAT_BURNING_VAPOR_TEMP)

  IF (sendTemp < 500._RFREAL .OR. sendTemp > 10000._RFREAL) THEN
    CALL ErrorStop( global,ERR_INRT_BADVAL,__LINE__ )
  END IF ! sendTemp

! Loop over all the particles -------------------------------------------------

  DO iPcls = 1,nPcls

    diamL = pDv(DV_PLAG_DIAM,iPcls)
    tempG = pDv(DV_PLAG_TEMPMIXT,iPcls)
    presG = pDV(DV_PLAG_PRESMIXT,iPcls)

    SELECT CASE (burnModel)

      CASE (INRT_BURNING_MODEL_BECKSTEAD)

! ----- Compute Particle Volume Fraction --------------------------------------

        volFracL =  pCv(pCvPlagMass(iContIn),iPcls)/             &
                  ( pCv(pCvPlagMass(iContIn),iPcls)              &
                  + pCv(pCvPlagMass(iContOut),iPcls) * densRatio )

! ----- If oxidizer is used, compute xiEffG -----------------------------------

        IF (oxUsed) THEN

! ------- Compute gas and oxidizer mass in this cell

          iCell = pAiv(AIV_PLAG_ICELLS,iPcls)

          massMixt = pCvMixt(CV_MIXT_DENS,iCell)
#ifdef RFLO
#ifdef PEUL
          massOx   = pCvPeul(iPeulOx,iCell)
#endif
#endif
#ifdef RFLU
          massOx   = pCvPeul(iPeulOx,iCell)
#endif

          IF (massOx <= 0._RFREAL) THEN

            xiEffG = 0._RFREAL

          ELSE

! --------- Notes: this assumes that the molecular weights of the gas
! --------- and oxidizer are the same.  Technically, xiEffG should be
! --------- massOx / massMixt, but that would be too dangerous.

            xiEffG = massOx / (massMixt + massOx)

          END IF ! massOx

        END IF ! oxUsed

! ----- Compute Mass Burning Rate of Al ---------------------------------------

        mDotBurn = 2.885E-13_RFREAL * densAL * ( tempG**expTempG ) &
                     * ( xiEffG**expXiEffG ) * ( presG**expPresG ) &
                     * (  diamL**expDiamL  ) * diffRel * volFracL

      CASE DEFAULT
        CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

    END SELECT ! burnModel

! - Fill the interaction source terms -----------------------------------------

    pPlag%inrtSources(INRT_BURNING_G_MASS_X,iPcls) = (molwRatio - 1.0_RFREAL) &
                                                   * mDotBurn

    pPlag%inrtSources(INRT_BURNING_L_MASS_X,iPcls) = mDotBurn

    IF (oxUsed) THEN
      pPlag%inrtSources(INRT_BURNING_S_MASS_X0 + nPeulOxEdges,iPcls) = &
        (molwRatio - 1.0_RFREAL) * mDotBurn
    END IF ! oxUsed

    IF (sendToVapor .AND. tempG > sendTemp) THEN
      pPlag%inrtSources(INRT_BURNING_X_ENER_G  + nPeulOxEdges,iPcls) = &
        0.0_RFREAL
      pPlag%inrtSources(INRT_BURNING_X_ENER_LV + nPeulOxEdges,iPcls) = &
        heatLatent*mDotBurn
    ELSE
      pPlag%inrtSources(INRT_BURNING_X_ENER_G  + nPeulOxEdges,iPcls) = &
        heatLatent*mDotBurn
      pPlag%inrtSources(INRT_BURNING_X_ENER_LV + nPeulOxEdges,iPcls) = &
        0.0_RFREAL
    END IF ! sendToVapor

    pPlag%inrtSources(INRT_BURNING_X_MASS_G + nPeulOxEdges,iPcls) = 0.0_RFREAL

    pPlag%inrtSources(INRT_BURNING_X_MASS_L + nPeulOxEdges,iPcls) = &
      mFracL*molwRatio*mDotBurn

    DO iPeulOutEdge = 1,nPeulOutEdges
      coeffPeul = (1.0_RFREAL-mFracL)*molwRatio * &
                  pInrtBurn%data(INRT_DAT_BURNING_MFRC_PEUL0 + iPeulOutEdge)

      pPlag%inrtSources(INRT_BURNING_X_MASS_S0 + nPeulOxEdges + iPeulOutEdge, &
        iPcls) = coeffPeul*mDotBurn
    END DO ! iPeulOutEdge

  END DO ! iPcls

! finalize --------------------------------------------------------------------

999  CONTINUE
#endif
  CALL DeregisterFunction( global )

END SUBROUTINE INRT_CalcBurning

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_CalcBurning.F90,v $
! Revision 1.4  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/02/15 20:18:17  wasistho
! put peul within ifdef
!
! Revision 1.1  2004/12/01 21:56:13  fnajjar
! Initial revision after changing case
!
! Revision 1.9  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.8  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.7  2004/03/02 21:47:29  jferry
! Added After Update interactions
!
! Revision 1.6  2004/01/31 03:59:22  haselbac
! Initial integration for Rocflu and Rocpart
!
! Revision 1.5  2003/09/19 20:35:26  jferry
! Implemented oxidizer species for burning interaction
!
! Revision 1.4  2003/04/04 16:27:39  jferry
! fixed inconsistency in use of burn rate information
!
! Revision 1.3  2003/04/03 22:52:04  fnajjar
! Bug fix for INRT data
!
! Revision 1.2  2003/04/03 21:10:17  jferry
! implemented additional safety checks for rocinteract
!
! Revision 1.1  2003/04/03 16:19:28  fnajjar
! Initial Import of routines for burning and scouring
!
!******************************************************************************







