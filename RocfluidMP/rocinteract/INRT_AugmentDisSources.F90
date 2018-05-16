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
! Purpose: converts the transfers of primary quantities along Edges for each
!          particle into augmentations of RHS terms for all quantities, for
!          an interaction involving Lagrangian particles
!
! Description: none.
!
! Input: iInrt = index of interaction
!
! Output: augments region%levels(iLev)%...%rhs structures
!
! Notes:
!
!   The RHS structures use opposite sign as the input source structure
!
!   For efficiency, this routine requires Nodes to be stored in this order:
!     Mixture, Lagrangian particle, Eulerian particle, Internal
!
!   The energy corresponding to a gas mass is actually taken to be the
!   enthalpy because energy added in reactions are measured as enthalpies
!
!******************************************************************************
!
! $Id: INRT_AugmentDisSources.F90,v 1.5 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_AugmentDisSources( region,iInrt )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInteract,   ONLY : t_inrt_input,t_inrt_interact
  USE ModMixture,    ONLY : t_mixt
  USE ModPartLag,    ONLY : t_plag
#ifdef RFLO
#ifdef PEUL
  USE ModPartEul,    ONLY : t_peul
#endif
#endif
#ifdef RFLU
  USE ModSpecies,    ONLY : t_spec
#endif
  USE ModError
  USE ModParameters
  USE INRT_ModParameters
#ifdef PLAG
  USE PLAG_ModParameters
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region
  INTEGER,        INTENT(IN)            :: iInrt

! ... loop variables
  INTEGER :: iPcls,iPlag,iPeul,iEdge

! ... local variables
  INTEGER, PARAMETER :: MAX_NODES = 101

  CHARACTER(CHRLEN)  :: RCSIdentString

  LOGICAL :: computeAux

  INTEGER :: nPcls,nPlag,nPeul,nIntl,nInputEdges,nNodes,nEdges
  INTEGER :: indMixt,indPlag0,indPeul0,indIntl
  INTEGER :: indPlagVapor,indPlagn,indPeuln
  INTEGER :: indCp,gasModel,ic,iNod,off1beg,ic1beg
#ifdef RFLO
  INTEGER :: iLev
#endif
  INTEGER :: tEdge,iNode(2),token(2)
  INTEGER, POINTER :: pCvPlagMass(:), aiv(:,:)

  REAL(RFREAL) :: temp1,temp2,spht,massdot,hcapdot,enerdot
  REAL(RFREAL) :: kinedot1,kinedot2,thrmdot1,thrmdot2,enerdot1,enerdot2
  REAL(RFREAL) :: intlMass,intlEner,intlTemp,intlHcap,contFac
  REAL(RFREAL) :: sphtPlag(MAX_NODES),sphtPeul(MAX_NODES)
  REAL(RFREAL), DIMENSION(3) :: velo1,velo2,intlVelo,intlMome
  REAL(RFREAL), DIMENSION(3) :: momedot,momedot1,momedot2
  REAL(RFREAL)               :: src(MAX_NODES,5)
  REAL(RFREAL), DIMENSION(:),   POINTER :: p1begMixtTemp,pPlagTemp
  REAL(RFREAL), DIMENSION(:,:), POINTER :: p1begMixtVelo,pPlagVelo
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pMixtRhs,pPlagRhs,pPeulRhs
  REAL(RFREAL), DIMENSION(:,:), POINTER :: gv,arv,primary

  TYPE(t_inrt_input),    POINTER :: input
  TYPE(t_inrt_interact), POINTER :: inrt
  TYPE(t_mixt),          POINTER :: pMixt
#ifdef RFLO
#ifdef PLAG
  TYPE(t_plag),          POINTER :: pPlag
#endif
#ifdef PEUL
  TYPE(t_peul),          POINTER :: pPeul
#endif
#endif
#ifdef RFLU
  TYPE(t_plag),          POINTER :: pPlag
  TYPE(t_spec),          POINTER :: pPeul
#endif
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_AugmentDisSources.F90,v $ $Revision: 1.5 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_AugmentDisSources',&
  'INRT_AugmentDisSources.F90' )

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

#ifdef RFLU
! Check that have primitive state vector --------------------------------------

  IF ( region%mixt%cvState /= CV_MIXT_STATE_DUVWP ) THEN
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! region%mixt%cvState
#endif

! initialize interaction constants and pointers

  input => region%inrtInput
  inrt  => input%inrts(iInrt)

  nPlag = input%nPlag
  nPeul = input%nPeul

  nIntl       = inrt%nIntl
  nInputEdges = inrt%nInputEdges

  indMixt  = input%indMixt
  indPlag0 = input%indPlag0
  indPeul0 = input%indPeul0
  indIntl  = input%indIntl

  indPlagVapor = input%indPlagVapor

  indPlagn = indPlag0 + nPlag
  indPeuln = indPeul0 + nPeul

  nNodes = input%nNodes
  IF (nNodes > MAX_NODES .OR. nPlag > MAX_NODES .OR. nPeul > MAX_NODES) &
    CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

  nEdges = inrt%nEdges

  computeAux = input%computeAux

! initialize data constants and pointers

#ifdef RFLO
  pMixt => region%levels(iLev)%mixt
  pPlag => region%levels(iLev)%plag
#ifdef PEUL
  pPeul => region%levels(iLev)%peul
#endif
  indCp =  pMixt%indCp
#endif
#ifdef RFLU
  pMixt => region%mixt
  pPlag => region%plag
  pPeul => region%spec

  indCp =  region%mixtInput%indCp
#endif

  pCvPlagMass => pPlag%cvPlagMass
  aiv => pPlag%aiv
  arv => pPlag%arv

! Constructing pointers from sections of arrays can be dangerous
! Extensive checks are employed to ensure correct behavior

! The string "1beg" is used to identify arrays beginning at 1 instead of ibc

#ifdef RFLO
  p1begMixtVelo => pMixt%dv(DV_MIXT_UVEL:DV_MIXT_WVEL,:)
  off1beg = LBOUND(p1begMixtVelo,2) - LBOUND(pMixt%dv,2)
#endif
#ifdef RFLU
  p1begMixtVelo => pMixt%cv(CV_MIXT_XVEL:CV_MIXT_ZVEL,:)
  off1beg = LBOUND(p1begMixtVelo,2) - LBOUND(pMixt%cv,2)
#endif

  IF (UBOUND(p1begMixtVelo,1) /= 3) THEN
    CALL ErrorStop( global,ERR_INRT_INDEXRANGE,__LINE__ )
  ENDIF

  p1begMixtTemp => pMixt%dv(DV_MIXT_TEMP,:)

  IF (off1beg /= LBOUND(p1begMixtTemp,1) - LBOUND(pMixt%dv,2)) THEN
    CALL ErrorStop( global,ERR_INRT_INDEXRANGE,__LINE__ )
  ENDIF

! The string "1beg" is not needed here:  Plag arrays already begin at 1

  pPlagVelo => pPlag%dv(DV_PLAG_UVEL:DV_PLAG_WVEL,:)

  IF (UBOUND(pPlagVelo,1) /= 3) THEN
    CALL ErrorStop( global,ERR_INRT_INDEXRANGE,__LINE__ )
  ENDIF

  IF (LBOUND(pPlagVelo,2) /= LBOUND(pPlag%dv,2)) THEN
    CALL ErrorStop( global,ERR_INRT_INDEXRANGE,__LINE__ )
  ENDIF

  pPlagTemp => pPlag%dv(DV_PLAG_TEMP,:)

  IF (LBOUND(pPlagTemp,1) /= LBOUND(pPlag%dv,2)) THEN
    CALL ErrorStop( global,ERR_INRT_INDEXRANGE,__LINE__ )
  ENDIF

  pMixtRhs => pMixt%rhs
  pPlagRhs => pPlag%rhs
#ifdef RFLO
#ifdef PEUL
  pPeulRhs => pPeul%rhs
#endif
#endif
#ifdef RFLU
  pPeulRhs => pPeul%rhs
#endif
  gasModel =  region%mixtInput%gasModel
  gv       => pMixt%gv

  IF (nPlag > 0) sphtPlag(1:nPlag) = region%plagInput%spht(1:nPlag)

  DO iPeul = 1,nPeul
#ifdef RFLO
#ifdef PEUL
    sphtPeul(iPeul) = region%peulInput%ptypes(iPeul)%material%spht
#endif
#endif
! use spht for now:  should be spht at constant pressure (or volume?)
#ifdef RFLU
    sphtPeul(iPeul) = region%specInput%specType(iPeul)%pMaterial%spht
#endif
  END DO ! iPeul

  primary => pPlag%inrtSources

  DO iPcls = 1,nPcls

    ic = aiv(AIV_PLAG_ICELLS,iPcls)
    ic1beg = ic + off1beg

! ----------------------
! - Compute source terms
! ----------------------

    src(:nNodes,:) = 0._RFREAL
    intlHcap       = 0._RFREAL

    DO iEdge = 1, nEdges

      tEdge = inrt%edges(iEdge)%tEdge
      iNode = inrt%edges(iEdge)%iNode
      token = inrt%edges(iEdge)%token

      SELECT CASE (tEdge)

! --- Mass Edge

      CASE (INRT_EDGE_MASS)

        massdot = primary(iEdge,iPcls)
        IF (massdot == 0._RFREAL) CYCLE ! do not waste time with null transfer

! ----- BEGIN INLINE: makeSpecificHeat(iNode,spht) ------------------------BEGI

! ----- Nodes at both ends should have the same specific heat, so use either.
! ----- Restriction: cannot use an internal Node.

        IF (iNode(1) == indIntl) THEN
          iNod = iNode(2)
        ELSE
          iNod = iNode(1)
        END IF ! iNode(1)

        IF (iNod == indMixt) THEN

          SELECT CASE (gasModel)

          CASE (GAS_MODEL_TCPERF)
            spht = gv(GV_MIXT_CP,indCp*ic) ! specific heat at constant pressure

          CASE DEFAULT
            CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

          END SELECT ! gasModel

        ELSE IF (iNod <= indPlagn) THEN
          spht = sphtPlag(iNod-indPlag0)

        ELSE IF (iNod == indIntl .OR. iNod == indPlagVapor) THEN
          CALL ErrorStop( global,ERR_INRT_INDEXRANGE,__LINE__ )

        ELSE
          spht = sphtPeul(iNod-indPeul0)

        END IF ! iNod

! ----- END   INLINE: makeSpecificHeat(iNode,spht) ------------------------ENDI

        hcapdot = spht * massdot

! ----- BEGIN INLINE: makeVelocity(iNode(1),velo1) ------------------------BEGI
! ----- BEGIN INLINE: makeTemperature(iNode(1),temp1) ---------------------BEGI

        iNod = iNode(1)

        IF (iNod == indMixt) THEN
          velo1 = p1begMixtVelo(:,ic1beg)
          temp1 = p1begMixtTemp(ic1beg)
        ELSE IF (iNod <= indPlagn) THEN
          velo1 = pPlagVelo(:,iPcls)
          temp1 = pPlagTemp(iPcls)
        ELSE IF (iNod == indPlagVapor) THEN
          CALL ErrorStop( global,ERR_INRT_INDEXRANGE,__LINE__ )
        ELSE IF (iNod <= indPeuln) THEN
          velo1 = p1begMixtVelo(:,ic1beg) ! Sets Smoke velocity = Fluid velocity
          temp1 = p1begMixtTemp(ic1beg)   ! Sets Smoke temp     = Fluid temp
        ELSE
          velo1 = intlVelo
          temp1 = intlTemp
        END IF ! iNod

! ----- END   INLINE: makeVelocity(iNode(1),velo1) ------------------------ENDI
! ----- END   INLINE: makeTemperature(iNode(1),temp1) ---------------------ENDI

! ----- BEGIN INLINE: makeVelocity(iNode(2),velo2) ------------------------BEGI
! ----- BEGIN INLINE: makeTemperature(iNode(2),temp2) ---------------------BEGI

        iNod = iNode(2)

        IF (iNod == indMixt) THEN
          velo2 = p1begMixtVelo(:,ic1beg)
          temp2 = p1begMixtTemp(ic1beg)
        ELSE IF (iNod <= indPlagn) THEN
          velo2 = pPlagVelo(:,iPcls)
          temp2 = pPlagTemp(iPcls)
        ELSE IF (iNod == indPlagVapor) THEN
          CALL ErrorStop( global,ERR_INRT_INDEXRANGE,__LINE__ )
        ELSE IF (iNod <= indPeuln) THEN
          velo2 = p1begMixtVelo(:,ic1beg) ! Sets Smoke velocity = Fluid velocity
          temp2 = p1begMixtTemp(ic1beg)   ! Sets Smoke temp     = Fluid temp
        ELSE
          velo2 = 0._RFREAL ! intlVelo not defined yet
          temp2 = 0._RFREAL ! intlTemp not defined yet
! ------- If downwind Node is Internal, then augment its heat capacity
          intlHcap = intlHcap + hcapdot
        END IF ! iNod

! ----- END   INLINE: makeVelocity(iNode(2),velo2) ------------------------ENDI
! ----- END   INLINE: makeTemperature(iNode(2),temp2) ---------------------BEGI

        IF (computeAux) THEN

          momedot1 = massdot * velo1
          momedot2 = massdot * velo2

          kinedot1 = 0.5_RFREAL*DOT_PRODUCT(velo1,momedot1)
          kinedot2 = 0.5_RFREAL*DOT_PRODUCT(velo2,momedot2)

          thrmdot1 = hcapdot * temp1
          thrmdot2 = hcapdot * temp2

        ELSE

          momedot1 = 0._RFREAL
          momedot2 = 0._RFREAL

          kinedot1 = 0._RFREAL
          kinedot2 = 0._RFREAL

          thrmdot1 = 0._RFREAL
          thrmdot2 = 0._RFREAL

        ENDIF

! ----- Compute sources for upwind Node

        SELECT CASE(token(1))

        CASE(INRT_PERM_PMASS)
          src(iNode(1),1)   = src(iNode(1),1)   - massdot
          src(iNode(1),2:4) = src(iNode(1),2:4) - momedot1
          src(iNode(1),5  ) = src(iNode(1),5  ) - kinedot1 - thrmdot1

        CASE(INRT_PERM_PMOME)
          src(iNode(1),1)   = src(iNode(1),1)   - massdot
          src(iNode(1),2:4) = src(iNode(1),2:4) - momedot2
          src(iNode(1),5  ) = src(iNode(1),5  ) - kinedot2 - thrmdot1

        CASE(INRT_PERM_PALL)
          src(iNode(1),1)   = src(iNode(1),1)   - massdot
          src(iNode(1),2:4) = src(iNode(1),2:4) - momedot2
          src(iNode(1),5  ) = src(iNode(1),5  ) - kinedot2 - thrmdot2

        END SELECT ! token(1)

! ----- Compute sources for downwind Node

        SELECT CASE(token(2))

        CASE(INRT_PERM_PMASS)
          src(iNode(2),1)   = src(iNode(2),1)   + massdot
          src(iNode(2),2:4) = src(iNode(2),2:4) + momedot2
          src(iNode(2),5  ) = src(iNode(2),5  ) + kinedot2 + thrmdot2

        CASE(INRT_PERM_PMOME)
          src(iNode(2),1)   = src(iNode(2),1)   + massdot
          src(iNode(2),2:4) = src(iNode(2),2:4) + momedot1
          src(iNode(2),5  ) = src(iNode(2),5  ) + kinedot1 + thrmdot2

        CASE(INRT_PERM_PALL)
          src(iNode(2),1)   = src(iNode(2),1)   + massdot
          src(iNode(2),2:4) = src(iNode(2),2:4) + momedot1
          src(iNode(2),5  ) = src(iNode(2),5  ) + kinedot1 + thrmdot1

        END SELECT ! token(2)

! --- Momentum Edge

      CASE (INRT_EDGE_MOME)

        momedot = primary(iEdge:iEdge+2,iPcls)

! ----- BEGIN INLINE: makeVelocity(iNode(1),velo1) ------------------------BEGI

        iNod = iNode(1)

        IF (iNod == indMixt) THEN
          velo1 = p1begMixtVelo(:,ic1beg)
        ELSE IF (iNod <= indPlagn) THEN
          velo1 = pPlagVelo(:,iPcls)
        ELSE IF (iNod == indPlagVapor) THEN
          CALL ErrorStop( global,ERR_INRT_INDEXRANGE,__LINE__ )
        ELSE IF (iNod <= indPeuln) THEN
          velo1 = p1begMixtVelo(:,ic1beg) ! Sets Smoke velocity = Fluid velocity
        ELSE
          velo1 = intlVelo
        END IF ! iNod

! ----- END   INLINE: makeVelocity(iNode(1),velo1) ------------------------ENDI

! ----- BEGIN INLINE: makeVelocity(iNode(2),velo2) ------------------------BEGI

        iNod = iNode(2)

        IF (iNod == indMixt) THEN
          velo2 = p1begMixtVelo(:,ic1beg)
        ELSE IF (iNod <= indPlagn) THEN
          velo2 = pPlagVelo(:,iPcls)
        ELSE IF (iNod == indPlagVapor) THEN
          CALL ErrorStop( global,ERR_INRT_INDEXRANGE,__LINE__ )
        ELSE IF (iNod <= indPeuln) THEN
          velo2 = p1begMixtVelo(:,ic1beg) ! Sets Smoke velocity = Fluid velocity
        ELSE
          velo2 = 0._RFREAL ! intlVelo not defined yet
        END IF ! iNod

! ----- END   INLINE: makeVelocity(iNode(2),velo2) ------------------------ENDI

        IF (computeAux) THEN

          enerdot1 = DOT_PRODUCT(momedot,velo1)
          enerdot2 = DOT_PRODUCT(momedot,velo2)

        ELSE

          enerdot1 = 0._RFREAL
          enerdot2 = 0._RFREAL

        ENDIF

! ----- Compute sources for upwind Node

        SELECT CASE(token(1))

        CASE(INRT_PERM_PMOME)
          src(iNode(1),2:4) = src(iNode(1),2:4) - momedot
          src(iNode(1),5  ) = src(iNode(1),5  ) - enerdot1

        CASE(INRT_PERM_PALL)
          src(iNode(1),2:4) = src(iNode(1),2:4) - momedot
          src(iNode(1),5  ) = src(iNode(1),5  ) - enerdot2

        END SELECT ! token(1)

! ----- Compute sources for downwind Node

        SELECT CASE(token(2))

        CASE(INRT_PERM_PMOME)
          src(iNode(2),2:4) = src(iNode(2),2:4) + momedot
          src(iNode(2),5)   = src(iNode(2),5)   + enerdot2

        CASE(INRT_PERM_PALL)
          src(iNode(2),2:4) = src(iNode(2),2:4) + momedot
          src(iNode(2),5)   = src(iNode(2),5)   + enerdot1

        END SELECT ! token(2)

! --- Dummy momentum Edge

      CASE (INRT_EDGE_MOME_DUM)
        CYCLE

! --- Energy Edge

      CASE (INRT_EDGE_ENER)

        enerdot = primary(iEdge,iPcls)

        IF (token(1) == INRT_PERM_PALL) &
          src(iNode(1),5) = src(iNode(1),5) - enerdot

        IF (token(2) == INRT_PERM_PALL) &
          src(iNode(2),5) = src(iNode(2),5) + enerdot

! --- Ghost mass Edge

      CASE (INRT_EDGE_MASS_GHO)

! ----- Downwind node must have Block Token.  Check value of upwind node

        IF (token(1) >= INRT_PERM_PMASS) THEN

          massdot = primary(iEdge,iPcls)
          src(iNode(1),1) = src(iNode(1),1) - massdot

        END IF ! token(1)

      CASE DEFAULT
        CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

      END SELECT ! tEdge

! --- Compute velocity and temperature at internal Node if necessary

      IF (nIntl > 0) THEN
        IF (iEdge == nInputEdges) THEN

          intlMass = src(indIntl,1)
          intlMome = src(indIntl,2:4)
          intlEner = src(indIntl,5)

          IF (intlMass > 0._RFREAL) THEN
            intlVelo = intlMome / intlMass
          ELSE
            intlVelo = 0._RFREAL
          END IF ! intlMass

          IF (intlHcap > 0._RFREAL) THEN
            intlTemp = (intlEner - 0.5_RFREAL*DOT_PRODUCT(intlMome, &
                        intlVelo)) / intlHcap
          ELSE
            intlTemp = 0._RFREAL
          END IF ! intlHcap

        END IF ! iEdge
      END IF ! nIntl

    END DO ! iEdge

! ------------------
! - Augment RHS data
! ------------------

! - Continuum conversion factor: converts from single particle value to
! - superparticle value

! - Note: does not incorporate volume because rhs values are not per volume
! - for either particles or continuua, in contrast to cv values, which are
! - not per volume for particles, but are per volume for continuua.

    contFac = arv(ARV_PLAG_SPLOAD,iPcls)

! - Augment Gas Sources

    pMixtRhs(  CV_MIXT_DENS             ,ic) = &
      pMixtRhs(CV_MIXT_DENS             ,ic) - contFac*src(indMixt,1  )

    pMixtRhs(  CV_MIXT_XMOM:CV_MIXT_ZMOM,ic) = &
      pMixtRhs(CV_MIXT_XMOM:CV_MIXT_ZMOM,ic) - contFac*src(indMixt,2:4)

    pMixtRhs(  CV_MIXT_ENER             ,ic) = &
      pMixtRhs(CV_MIXT_ENER             ,ic) - contFac*src(indMixt,5  )

! - Augment Lagrangian Particle Sources

    DO iPlag = 1, nPlag

      pPlagRhs(  pCvPlagMass(iPlag)       ,iPcls) = &
        pPlagRhs(pCvPlagMass(iPlag)       ,iPcls) - src(indPlag0+iPlag,1)

      pPlagRhs(  CV_PLAG_XMOM:CV_PLAG_ZMOM,iPcls) = &
        pPlagRhs(CV_PLAG_XMOM:CV_PLAG_ZMOM,iPcls) - src(indPlag0+iPlag,2:4)

      pPlagRhs(  CV_PLAG_ENER             ,iPcls) = &
        pPlagRhs(CV_PLAG_ENER             ,iPcls) - src(indPlag0+iPlag,5)

    END DO ! iPlag

! - Augment Lagrangian Particle Vapor Energy

    pPlagRhs(  CV_PLAG_ENERVAPOR,iPcls) = &
      pPlagRhs(CV_PLAG_ENERVAPOR,iPcls) - src(indPlagVapor,5)

! - Augment Smoke Sources
#ifdef RFLO
#ifdef PEUL
    DO iPeul = 1, nPeul
      pPeulRhs(iPeul,ic) = pPeulRhs(iPeul,ic) - contFac * src(indPeul0+iPeul,1)
    END DO ! iPeul
#endif
#endif
#ifdef RFLU
    DO iPeul = 1, nPeul
      pPeulRhs(iPeul,ic) = pPeulRhs(iPeul,ic) - contFac * src(indPeul0+iPeul,1)
    END DO ! iPeul
#endif
  END DO ! iPcls

! finalize --------------------------------------------------------------------

999  CONTINUE
#endif
  CALL DeregisterFunction( global )

END SUBROUTINE INRT_AugmentDisSources

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_AugmentDisSources.F90,v $
! Revision 1.5  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/02/15 20:18:30  wasistho
! put peul within ifdef
!
! Revision 1.2  2005/10/31 21:09:37  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/12/01 21:56:10  fnajjar
! Initial revision after changing case
!
! Revision 1.15  2004/07/28 15:42:12  jferry
! deleted defunct constructs: useDetangle, useSmokeDrag, useSmokeHeatTransfer
!
! Revision 1.14  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.13  2004/03/25 21:14:53  jferry
! fixed pointer offset bug
!
! Revision 1.12  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.11  2004/03/02 21:47:29  jferry
! Added After Update interactions
!
! Revision 1.10  2004/01/31 03:59:22  haselbac
! Initial integration for Rocflu and Rocpart
!
! Revision 1.9  2003/09/19 20:35:26  jferry
! Implemented oxidizer species for burning interaction
!
! Revision 1.8  2003/05/08 17:17:14  jferry
! changed energy associated with mass to enthalpy
!
! Revision 1.7  2003/05/07 15:13:10  jferry
! Rearranged for efficiency
!
! Revision 1.6  2003/04/09 15:02:39  jferry
! removed erroneous volume normalization for continuum rhs
!
! Revision 1.5  2003/04/03 21:10:17  jferry
! implemented additional safety checks for rocinteract
!
! Revision 1.4  2003/04/02 22:32:03  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.3  2003/03/24 23:30:52  jferry
! overhauled rocinteract to allow interaction design to use user input
!
! Revision 1.2  2003/03/11 16:05:54  jferry
! Created data type for material properties
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************







