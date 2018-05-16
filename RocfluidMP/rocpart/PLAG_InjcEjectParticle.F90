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
! Purpose: inject particle based on pool volume
!          for the multiphase injection algorithm.
!
! Description: none.
!
! Input: region    = current region
!        iReg      = current region number.
!
! Output: regions(iReg)%levels%plag = injected plag values for cv, aiv, arv
!                                     of current region.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_InjcEjectParticle.F90,v 1.8 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_InjcEjectParticle( region, iReg )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag, t_plag_input, t_tile_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModRandom, ONLY     : Rand1Uniform
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset
  USE PLAG_ModInterfaces, ONLY : PLAG_InjcMakeParticle, PLAG_InjcSetInjection
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  USE INRT_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

  INTEGER :: iReg

! ... loop variables
  INTEGER :: i, iCont, iPatch, j, k

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: bcType, burnStat, ibeg, iCOff, idir, iend, ijCOff, ijNOff, &
             ijkC, ijkN,iLev, injcDiamDist, inode, iNOff, iPcls, iTile, &
             jbeg, jdir, jend, jnode,kbeg, kdir, kend, knode,           &
             lbound, n1, n2, nextIdNumber, nCont, nPatches, nOff

  INTEGER :: iFilePlag
  INTEGER :: iterCellLocate
  INTEGER :: ejecModel
  INTEGER :: nPcls,nPclsMax
  INTEGER, PARAMETER :: ITER_CELL_LOCATE_MAX = 20

  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass, pCvTileMass
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv, pAivOld

  LOGICAL :: injectQ
  LOGICAL :: cellLocate

  REAL(RFREAL) :: area, cellHeight, injcBeta, pi, plagVolRatio, &
                  plagMomeNrm, poolVolSum, sgn, sxn,syn,szn,   &
                  tInjcCoeff, tInjcSum, volMeanPart
  REAL(RFREAL) :: meanSuperParticleVolume, spload
  REAL(RFREAL) :: injcBetaFac, injcBetaFacInv
  REAL(RFREAL) :: poolVolOldSum, poolVolFinal, remainingVolume
  REAL(RFREAL) :: currentSuperParticleVolume, poolVolume
  REAL(RFREAL) :: poolExcess, pExcessS, possibleExcess
  REAL(RFREAL) :: countdown, countdownNext, deltaVolume, randUnif
  REAL(RFREAL) :: currentParticleVolume
  REAL(RFREAL) :: poolVolCurr

  REAL(RFREAL),          DIMENSION(3)   :: poolxyz, sNormal

  REAL(RFREAL), POINTER, DIMENSION(:)   :: pVol
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: sFace, pXyz

  REAL(RFREAL), POINTER, DIMENSION(:,:)   :: pArv, pCvPlag, pCvTile, &
                                             pCvOldTile, pDvTile
  REAL(RFREAL), POINTER, DIMENSION(:,:,:) :: pFc

  TYPE(t_patch),     POINTER :: pPatch
  TYPE(t_tile_plag), POINTER :: pTilePlag
  TYPE(t_plag),      POINTER :: pPlag
  TYPE(t_global),    POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InjcEjectParticle.F90,v $ $Revision: 1.8 $'

  global => region%global

  CALL RegisterFunction( global, 'PLAG_InjcEjectParticle',&
  'PLAG_InjcEjectParticle.F90' )

! Get dimensions --------------------------------------------------------------

  iLev      = region%currLevel
  nPatches  = region%nPatches

  nCont     = region%plagInput%nCont
  pi        = global%pi

  injcDiamDist  = region%plagInput%injcDiamDist
  injcBeta      = region%plagInput%injcBeta

  volMeanPart   = pi/6.0_RFREAL * region%plagInput%injcDiamMean**3

  ejecModel     = region%plagInput%ejecModel
  spLoad        = region%plagInput%spLoad

  IF ( ejecModel == PLAG_EJEC_CRE ) THEN
    meanSuperParticleVolume = volMeanPart *spLoad

    injcBetaFac = 2.0_RFREAL *injcBeta *meanSuperParticleVolume**2
    injcBetaFacInv = 1.0_RFREAL/injcBetaFac
  ENDIF ! ejecModel

! Set pointers ----------------------------------------------------------------

  pXyz       => region%levels(iLev)%grid%xyz

  pPlag       => region%levels(iLev)%plag
  pCvPlag     => pPlag%cv
  pCvPlagMass => pPlag%cvPlagMass
  pAiv        => pPlag%aiv
  pAivOld     => pPlag%aivOld
  pArv        => pPlag%arv

  pVol        => region%levels(iLev)%grid%vol
  pFc         => region%levels(iLev)%plag%fc

! ******************************************************************************
! Trap error if nPcls exceeds maximum datastructure dimension, nPclsMax
! ******************************************************************************

  nPcls    = pPlag%nPcls
  nPclsMax = region%plagInput%nPclsMax
  
  IF ( nPcls >= nPclsMax ) THEN
    WRITE(STDOUT,*) ' PLAG_InjcEjectParticle: Datastructure Dimension Exceeded ',&
                      nPcls,nPclsMax
    CALL ErrorStop( global,ERR_PLAG_MEMOVERFLOW,__LINE__ )

#ifdef MPI
    CALL MPI_Barrier( global%mpiComm,global%mpierr )
    IF ( global%mpierr /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

    CALL MPI_Finalize(global%mpierr)
    IF ( global%mpierr /= ERR_NONE ) &
      CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

  END IF ! nPcls

! Loop over patches -----------------------------------------------------------

  DO iPatch=1,nPatches

    pPatch  => region%levels(iLev)%patches(iPatch)

    bcType = pPatch%bcType

! - Select injection boundary condition ---------------------------------------

    IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) THEN

! - Get dimensions ------------------------------------------------------------

      lbound = pPatch%lbound

      CALL RFLO_GetPatchIndices( region,pPatch,iLev, &
                                 ibeg,iend,jbeg,jend,kbeg,kend )
      CALL RFLO_GetPatchDirection( pPatch,idir,jdir,kdir )
      CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
      CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

      nOff        = ABS(pPatch%l1end-pPatch%l1beg) + 1

! - Select right face vector and make it point inwards ------------------------

      sgn   = -1._RFREAL
      inode = 0
      jnode = 0
      knode = 0
      IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
        sgn   = +1._RFREAL
        inode = -idir
        jnode = -jdir
        knode = -kdir
      ENDIF ! lbound

! - Get the appropriate face vector -------------------------------------------

      IF (lbound==1 .OR. lbound==2) sFace => region%levels(iLev)%plag%si
      IF (lbound==3 .OR. lbound==4) sFace => region%levels(iLev)%plag%sj
      IF (lbound==5 .OR. lbound==6) sFace => region%levels(iLev)%plag%sk

! - Loop over patch -----------------------------------------------------------

      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkC  = IndIJK(i,j,k,iCOff,ijCOff)
            ijkN  = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)

            area  = SQRT(sFace(XCOORD,ijkN)*sFace(XCOORD,ijkN)+ &
                         sFace(YCOORD,ijkN)*sFace(YCOORD,ijkN)+ &
                         sFace(ZCOORD,ijkN)*sFace(ZCOORD,ijkN))

            sxn   = sgn*sFace(XCOORD,ijkN)/area
            syn   = sgn*sFace(YCOORD,ijkN)/area
            szn   = sgn*sFace(ZCOORD,ijkN)/area

            IF      (lbound==1 .OR. lbound==2) THEN
              n1 = j - jbeg + 1
              n2 = k - kbeg + 1
            ELSE IF (lbound==3 .OR. lbound==4) THEN
              n1 = k - kbeg + 1
              n2 = i - ibeg + 1
            ELSE IF (lbound==5 .OR. lbound==6) THEN
              n1 = i - ibeg + 1
              n2 = j - jbeg + 1
            ENDIF ! lbound

! -- Select tile pointers -----------------------------------------------------

            iTile  = IndIJ(n1,n2,nOff)
            pTilePlag => pPatch%tilePlag

            pCvTile     => pTilePlag%cv
            pCvTileMass => pTilePlag%cvTileMass
            pDvTile     => pTilePlag%dv

            pCvOldTile  => pTilePlag%cvOld

!------------------------------------------------------------------------------
!           Invoke ejection model 
!------------------------------------------------------------------------------

            SELECT CASE(ejecModel)

!------------------------------------------------------------------------------
!             Ejection Model 1
!------------------------------------------------------------------------------

              CASE(PLAG_EJEC_MODEL1)
                poolVolSum = SUM ( pCvTile( pCvTileMass(:),iTile ) / &
                                   region%plagInput%dens(:) )

                pDvTile(DV_TILE_POOLVOLD,iTile) = poolVolSum

                tInjcCoeff = injcBeta * volMeanPart

                tInjcSum   = 0.0_RFREAL

                CALL PLAG_InjcSetInjection( region, pTilePlag, iTile,         &
                                            tInjcCoeff, tInjcSum, poolVolSum, &
                                            injectQ, plagVolRatio )

! -- If injectQ set to True, PLAG datastructure for new injected particle -----

                DO WHILE ( injectQ )
                  CALL PLAG_RFLO_EjectParticle( region,pPlag,pTilePlag,iTile,&
                                                lbound,iNOff,ijNOff,i,j,k,   &
                                                sxn,syn,szn,area,plagVolRatio )

! --- Update tile diameter and superparticle loading factor -------------------

                  CALL PLAG_InjcMakeParticle( region, injcDiamDist,          &
                                              pDvTile(DV_TILE_DIAM,  iTile), &
                                              pDvTile(DV_TILE_SPLOAD,iTile)  )

! --- Check if another particle could be injected -----------------------------

                  CALL PLAG_InjcSetInjection( region, pTilePlag, iTile,         &
                                              tInjcCoeff, tInjcSum, poolVolSum, &
                                              injectQ, plagVolRatio  )

                END DO ! injectQ

!------------------------------------------------------------------------------
!             Ejection Model based on Conservative Random Ejection (CRE)
!------------------------------------------------------------------------------

              CASE(PLAG_EJEC_CRE)

! ------------- poolVolSum is the pool volume after the current timestep is finished
                poolVolSum    = SUM ( pCvTile( pCvTileMass(:),iTile ) / &
                                      region%plagInput%dens(:) )
! ------------- poolVolOldSum is the pool volume before the current timestep began
                poolVolOldSum = SUM ( pCvOldTile( pCvTileMass(:),iTile ) / &
                                      region%plagInput%dens(:) )

! ------------- poolVolume is built up incrementally with volume injection,
! -------------   and is decremented when (super)particles are ejected
                poolVolume   = poolVolOldSum

! ------------- poolVolFinal starts with all volume injection having occurred,
! -------------   and is decremented when (super)particles are ejected
                poolVolFinal = poolVolSum

! ------------- poolVolCurr starts with all volume injection having occurred,
! -------------   and is decremented when (super)particles are ejected
                poolVolCurr = poolVolSum

! ------------- remainingVolume is the volume to be added to the pool over the
! -------------   current time step: it decreases to 0 as all volume is added
                remainingVolume = poolVolSum -poolVolOldSum

                creEjectLoop: DO
                  currentSuperParticleVolume = pi/6.0_RFREAL *spLoad &
                                             * pDvTile(DV_TILE_DIAM,iTile)**3.0_RFREAL

                  currentParticleVolume = pi/6.0_RFREAL &
                                        * pDvTile(DV_TILE_DIAM,iTile)**3.0_RFREAL

! --------------- poolExcess is the volume that would remain in the pool after ejection
!                  poolExcess = poolVolume -currentSuperParticleVolume

! TEMPORARY
                  poolExcess = poolVolCurr -currentSuperParticleVolume
! END TEMPORARY

! --------------- possibleExcess is the volume that would remain in the pool at the end of
! ---------------   the timestep after the current particle is (and no others are) ejected
                  possibleExcess = poolExcess + remainingVolume

! --------------- pExcessS is an auxiliary variable that occurs in a few formulas
                  IF ( poolExcess > 0.0_RFREAL ) THEN
                    pExcessS = poolExcess**2
                  ELSE
                    pExcessS = 0.0_RFREAL
                  ENDIF

! --------------- countdown is the internal ejection clock for the current particle:  it
! ---------------   starts to tick down when pool volume > current (super)particle volume
                  countdown = pDvTile(DV_TILE_COUNTDOWN,iTile)

! --------------- countdownNext is what the clock could tick down to within remainingVolume
                  IF ( possibleExcess > 0.0_RFREAL ) THEN
                    countdownNext = countdown + injcBetaFacInv &
                                  * (pExcessS -possibleExcess**2.0_RFREAL)
                  ELSE
                    countdownNext = countdown
                  ENDIF

!------------------------------------------------------------------------------
!                 Pool volume too small to eject particle
!------------------------------------------------------------------------------

! --------------- if the clock cannot tick down to zero within remainingVolume, no more
! ---------------   particles will be ejected in this timestep

                  IF ( countdownNext> 0.0_RFREAL ) THEN
                    poolVolume = poolVolume + remainingVolume ! not needed except as check
                    pDvTile(DV_TILE_COUNTDOWN,iTile) = countdownNext
                    EXIT creEjectLoop

!------------------------------------------------------------------------------
!                 Pool volume big enough to eject particle
!------------------------------------------------------------------------------

! --------------- if the clock can tick down below zero, then we need to find when it
! ---------------   reaches zero, and then eject the particle
                  ELSE

! ----------------- deltaVolume is the volume that must be added to the pool (out of
! -----------------   remainingVolume) in order to bring countdown exactly to zero
                    deltaVolume = SQRT(injcBetaFac*countdown +pExcessS) -poolExcess

! ----------------- therefore deltaVolume is added to the pool, the countdown hits zero,
! -----------------   and the particle is ejected
                    poolVolume  = poolVolume +deltaVolume -currentSuperParticleVolume

! ----------------- the deltaVolume added to the pool is taken from remainingVolume
                    remainingVolume = remainingVolume - deltaVolume

! ----------------- plagVolRatio must be in terms of poolVolFinal:  the pool cv variables
! -----------------   are expressed in terms of values at the end of the time step, not
! -----------------   in terms of the intermediate values like poolVolume

                    plagVolRatio = currentParticleVolume / poolVolFinal

! ----------------- poolVolFinal, like the cv variables, is affected by particle ejection
                    poolVolFinal = poolVolFinal - currentSuperParticleVolume

! ----------------- eject particle

                    CALL PLAG_RFLO_EjectParticle( region,pPlag,pTilePlag,iTile, &
                                                  lbound,iNOff,ijNOff,i,j,k,    &
                                                  sxn,syn,szn,area,plagVolRatio )

! ----------------- make new particle diameter and superparticle loading factor
                    CALL PLAG_InjcMakeParticle( region, injcDiamDist,          &
                                                pDvTile(DV_TILE_DIAM,  iTile), &
                                                pDvTile(DV_TILE_SPLOAD,iTile)  )

! TEMPORARY
! ----------------- Compute current pool volume after being decreased 
! ----------------- in PLAG_EjectParticle     
                    poolVolCurr  = SUM ( pCvTile( pCvTileMass(:),iTile ) / &
                                         region%plagInput%dens(:) )

! END TEMPORARY

! ----------------- set countdown for new particle
                    randUnif = Rand1Uniform(region%randData)
                    IF ( randUnif <= 0.0_RFREAL) THEN
! ------------------- treats randUnif as 1.9e-22
                      pDvTile(DV_TILE_COUNTDOWN,iTile) = 50.0_RFREAL
                    ELSE
                      pDvTile(DV_TILE_COUNTDOWN,iTile) = -LOG(randUnif)
                    END IF ! randUnif

                  END IF ! countdownNext
                END DO creEjectLoop

! ------------- Temporary check: poolVolume and poolVolFinal should now agree

!                print *, 'CRE check: iTile ', iTile,&
!                  (poolVolume - poolVolFinal) / poolVolFinal

            END SELECT ! ejecModel

          END DO ! i
        END DO ! j
      END DO ! k

      iFilePlag = 200+iReg

!      WRITE(STDOUT,'(A,I4,I8)') '     iReg: nPcls = ',iReg, pPlag%nPcls
!      WRITE(iFilePlag,'(1PE12.5,10(2X,I8))') global%currentTime+global%dtMin,iReg, &
!      pPlag%nPcls,&
!      pPatch%tilePlag%nPclsInjc(1),pPatch%tilePlag%nPclsInjc(10),  &
!!      pPatch%tilePlag%nPclsInjc(20),pPatch%tilePlag%nPclsInjc(50), &
!      SUM(pTilePlag%nPclsInjc(:))

! -- Reset injected particle size ---------------------------------------------

      pPatch%tilePlag%nPclsInjc(:) = 0

    ENDIF  ! bcType

  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

CONTAINS

!******************************************************************************

   SUBROUTINE PLAG_RFLO_EjectParticle(region,pPlag,pTilePlag,iTile,lbound, &
                                      iNOff,ijNOff,i,j,k,sxn,syn,szn,area, &
                                      plagVolRatio)

!... parameters
     INTEGER,      INTENT(IN)  :: lbound,iNOff,ijNOff,iTile,i,j,k

     REAL(RFREAL), INTENT(IN)  :: sxn,syn,szn,area,plagVolRatio
     
     TYPE(t_region) :: region
     TYPE(t_tile_plag), POINTER :: pTilePlag
     TYPE(t_plag),      POINTER :: pPlag     

!... local variables

     INTEGER :: burnStat,iCont,iterCellLocate,nCont
     INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass, pCvTileMass
     INTEGER, POINTER, DIMENSION(:,:) :: pAiv, pAivOld

     REAL(RFREAL) :: cellHeight 
     REAL(RFREAL), DIMENSION(3) :: poolxyz,sNormal
     REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv, pCvPlag, pCvTile, &
                                              pCvOldTile, pDvTile

!------------------------------------------------------------------------------
!    Get dimensions and set pointers
!------------------------------------------------------------------------------
 
     nCont = region%plagInput%nCont
     
     pCvPlag     => pPlag%cv
     pCvPlagMass => pPlag%cvPlagMass
     pAiv        => pPlag%aiv
     pAivOld     => pPlag%aivOld
     pArv        => pPlag%arv
     
     pCvTile     => pTilePlag%cv
     pCvTileMass => pTilePlag%cvTileMass
     pDvTile     => pTilePlag%dv

     pCvOldTile  => pTilePlag%cvOld

!------------------------------------------------------------------------------
!    Initial burning status of particles 
!------------------------------------------------------------------------------

! Currently all particles start off burning if the burning interaction is used

     IF (region%inrtInput%inrts(INRT_TYPE_BURNING)%used) THEN
       burnStat = INRT_BURNSTAT_ON
     ELSE
       burnStat = INRT_BURNSTAT_OFF
     ENDIF

!------------------------------------------------------------------------------
!    Increment PLAG datastructure for new injected particle 
!------------------------------------------------------------------------------

      pPlag%nPcls = pPlag%nPcls + 1
      iPcls       = pPlag%nPcls

      pTilePlag%nPclsInjc(iTile)   = pTilePlag%nPclsInjc(iTile) + 1

!------------------------------------------------------------------------------
!     Keep track of global count in particle ejected in region 
!------------------------------------------------------------------------------

      pPlag%nextIdNumber = pPlag%nextIdNumber + 1
      nextIdNumber       = pPlag%nextIdNumber

!------------------------------------------------------------------------------
!     Compute PLAG cv datastructure for new injected particle 
!------------------------------------------------------------------------------

       DO iCont = 1, nCont
         pCvPlag(pCvPlagMass(iCont),iPcls) = pCvTile(pCvTileMass(iCont),iTile) * &
                                                    plagVolRatio
       END DO ! iCont

       pCvPlag(CV_PLAG_ENER,iPcls) = pCvTile(CV_TILE_ENER,  iTile) *plagVolRatio

       plagMomeNrm                 = pCvTile(CV_TILE_MOMNRM,iTile) *plagVolRatio
       pCvPlag(CV_PLAG_XMOM,iPcls) = plagMomeNrm*sxn
       pCvPlag(CV_PLAG_YMOM,iPcls) = plagMomeNrm*syn
       pCvPlag(CV_PLAG_ZMOM,iPcls) = plagMomeNrm*szn

       pCvPlag(CV_PLAG_ENERVAPOR,iPcls) = 0._RFREAL

!------------------------------------------------------------------------------
!      Determine particle position 
!      include a test to insure particle is in Tile Cell
!------------------------------------------------------------------------------

       sNormal(XCOORD) = sxn
       sNormal(YCOORD) = syn
       sNormal(ZCOORD) = szn
       cellHeight      = pVol(ijkC)/area

       iterCellLocate = 1

199    CONTINUE
         CALL PLAG_InjcSetPositions( lbound,iNOff,ijNOff,i,j,k, &
                                     sNormal,cellHeight,poolxyz )

         CALL PLAG_InjcTestCell( region,iNOff,ijNOff,i,j,k, &
                                 poolxyz,cellLocate         )

         IF ( cellLocate .EQV. .FALSE. ) THEN
            iterCellLocate = iterCellLocate+1
            IF ( iterCellLocate > ITER_CELL_LOCATE_MAX ) THEN
              WRITE(*,*) ' PLAG_InjcEjectParticle: Unable to create a particle after ',&
                           ITER_CELL_LOCATE_MAX, ' iterations'

              CALL ErrorStop( global,ERR_PLAG_CELLINDEX,__LINE__ )
             END IF ! iterCellLocate
             GOTO 199
          ENDIF ! cellLocate

          pCvPlag(CV_PLAG_XPOS,iPcls) = poolxyz(XCOORD)
          pCvPlag(CV_PLAG_YPOS,iPcls) = poolxyz(YCOORD)
          pCvPlag(CV_PLAG_ZPOS,iPcls) = poolxyz(ZCOORD)

!------------------------------------------------------------------------------
!         Set aiv and arv
!------------------------------------------------------------------------------

          pAiv(AIV_PLAG_PIDINI,iPcls) = nextIdNumber
          pAiv(AIV_PLAG_REGINI,iPcls) = iReg
          pAiv(AIV_PLAG_REGCRT,iPcls) = iReg
          pAiv(AIV_PLAG_ICELLS,iPcls) = ijkC
          pAiv(AIV_PLAG_INDEXI,iPcls) = i
          pAiv(AIV_PLAG_INDEXJ,iPcls) = j
          pAiv(AIV_PLAG_INDEXK,iPcls) = k
          pAiv(AIV_PLAG_BURNSTAT,iPcls) = burnStat
          pAiv(AIV_PLAG_STATUS,iPcls) = PLAG_STATUS_KEEP

          pArv(ARV_PLAG_SPLOAD,iPcls) = pDvTile(DV_TILE_SPLOAD,iTile)

!------------------------------------------------------------------------------
!         Load at injection aivOld array
!------------------------------------------------------------------------------

          pAivOld(:,iPcls) =    pAiv(:,iPcls)

!------------------------------------------------------------------------------
!         Set time of particle injection 
!            timeinjectPart = global%currentTime+tInjcSum*global%dtMin
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!         Decrease pool variables
!------------------------------------------------------------------------------

          DO iCont = 1, nCont
            pCvTile( pCvTileMass(iCont),iTile) = pCvTile(pCvTileMass(iCont),iTile) -  &
                                                 pDvTile(DV_TILE_SPLOAD,    iTile) *  &
                                                 pCvPlag(pCvPlagMass(iCont),iPcls)
          END DO ! iCont

          pCvTile(CV_TILE_MOMNRM,iTile) = pCvTile(CV_TILE_MOMNRM,iTile) - &
                                          pDvTile(DV_TILE_SPLOAD,iTile) * &
                                          plagMomeNrm

          pCvTile(CV_TILE_ENER  ,iTile) = pCvTile(CV_TILE_ENER  ,iTile) - &
                                          pDvTile(DV_TILE_SPLOAD,iTile) * &
                                          pCvPlag(CV_PLAG_ENER,iPcls)


   END SUBROUTINE PLAG_RFLO_EjectParticle
!******************************************************************************

!******************************************************************************
   SUBROUTINE PLAG_InjcSetPositions( lbound,iNOff, ijNOff,i,j,k,    &
                                     sNormal,cellHeight,posPlag )

!... parameters
     INTEGER,      INTENT(IN)  :: lbound,iNOff,ijNOff,i,j,k

     REAL(RFREAL), INTENT(IN)                :: cellHeight
     REAL(RFREAL), DIMENSION(3), INTENT(IN)  :: sNormal
     REAL(RFREAL), DIMENSION(3), INTENT(OUT) :: posPlag

!... local variables
     INTEGER,  DIMENSION(4)  :: corner

     LOGICAL                 :: useTriangle1

     REAL(RFREAL), PARAMETER :: HEIGHT_FRACTION = 1.0E-3_RFREAL
     REAL(RFREAL) :: areaTriangle1, areaTriangle2, xRand, yRand, zRand

     REAL(RFREAL),  DIMENSION(ZCOORD)   :: faceCentroid,v1,v2
     REAL(RFREAL),  DIMENSION(ZCOORD,4) :: pNode

!******************************************************************************

! Select the coordinates of the four points on the injecting face -------------
! Use unshifted indices to get correct extents
!  Note: pNode(1) = pNode00, pNode(2) = pNode10
!        pNode(3) = pNode01, pNode(4) = pNode11

     SELECT CASE (lbound)
       CASE(1,2)
         corner(1) = IndIJK(i+inode  ,j+jnode  ,k+knode  ,iNOff,ijNOff)
         corner(2) = IndIJK(i+inode  ,j+jnode+1,k+knode  ,iNOff,ijNOff)
         corner(3) = IndIJK(i+inode  ,j+jnode  ,k+knode+1,iNOff,ijNOff)
         corner(4) = IndIJK(i+inode  ,j+jnode+1,k+knode+1,iNOff,ijNOff)

         pNode(XCOORD:ZCOORD,1) = pXyz(XCOORD:ZCOORD,corner(1))
         pNode(XCOORD:ZCOORD,2) = pXyz(XCOORD:ZCOORD,corner(2))
         pNode(XCOORD:ZCOORD,3) = pXyz(XCOORD:ZCOORD,corner(3))
         pNode(XCOORD:ZCOORD,4) = pXyz(XCOORD:ZCOORD,corner(4))

         faceCentroid(XCOORD:ZCOORD) = (/pFc(XCOORD,ICOORD,corner(1)), &
                                         pFc(YCOORD,ICOORD,corner(1)), &
                                         pFc(ZCOORD,ICOORD,corner(1))/)

       CASE(3,4)
         corner(1) = IndIJK(i+inode  ,j+jnode  ,k+knode  ,iNOff,ijNOff)
         corner(2) = IndIJK(i+inode+1,j+jnode  ,k+knode  ,iNOff,ijNOff)
         corner(3) = IndIJK(i+inode  ,j+jnode  ,k+knode+1,iNOff,ijNOff)
         corner(4) = IndIJK(i+inode+1,j+jnode  ,k+knode+1,iNOff,ijNOff)

         pNode(XCOORD:ZCOORD,1) = pXyz(XCOORD:ZCOORD,corner(1))
         pNode(XCOORD:ZCOORD,2) = pXyz(XCOORD:ZCOORD,corner(2))
         pNode(XCOORD:ZCOORD,3) = pXyz(XCOORD:ZCOORD,corner(3))
         pNode(XCOORD:ZCOORD,4) = pXyz(XCOORD:ZCOORD,corner(4))

         faceCentroid(XCOORD:ZCOORD) = (/pFc(XCOORD,JCOORD,corner(1)), &
                                         pFc(YCOORD,JCOORD,corner(1)), &
                                         pFc(ZCOORD,JCOORD,corner(1))/)

       CASE(5,6)
         corner(1) = IndIJK(i+inode  ,j+jnode  ,k+knode  ,iNOff,ijNOff)
         corner(2) = IndIJK(i+inode+1,j+jnode  ,k+knode  ,iNOff,ijNOff)
         corner(3) = IndIJK(i+inode  ,j+jnode+1,k+knode  ,iNOff,ijNOff)
         corner(4) = IndIJK(i+inode+1,j+jnode+1,k+knode  ,iNOff,ijNOff)

         pNode(XCOORD:ZCOORD,1) = pXyz(XCOORD:ZCOORD,corner(1))
         pNode(XCOORD:ZCOORD,2) = pXyz(XCOORD:ZCOORD,corner(2))
         pNode(XCOORD:ZCOORD,3) = pXyz(XCOORD:ZCOORD,corner(3))
         pNode(XCOORD:ZCOORD,4) = pXyz(XCOORD:ZCOORD,corner(4))

         faceCentroid(XCOORD:ZCOORD) = (/pFc(XCOORD,KCOORD,corner(1)), &
                                         pFc(YCOORD,KCOORD,corner(1)), &
                                         pFc(ZCOORD,KCOORD,corner(1))/)


     END SELECT ! lbound

! Partition the face into two triangles and compute the areas -----------------
!  of their projections normal to sNormal -------------------------------------

     v1(1:3) = pNode(XCOORD:ZCOORD,2)-pNode(XCOORD:ZCOORD,1)
     v2(1:3) = pNode(XCOORD:ZCOORD,3)-pNode(XCOORD:ZCOORD,1)

     areaTriangle1 = 0.5_RFREAL*ABS( sNormal(XCOORD)*(v1(2)*v2(3)-v1(3)*v2(2)) &
                                   + sNormal(YCOORD)*(v1(3)*v2(1)-v1(1)*v2(3)) &
                                   + sNormal(ZCOORD)*(v1(1)*v2(2)-v1(2)*v2(1)) )

     v1(XCOORD:ZCOORD) = pNode(XCOORD:ZCOORD,2)-pNode(XCOORD:ZCOORD,4)
     v2(XCOORD:ZCOORD) = pNode(XCOORD:ZCOORD,3)-pNode(XCOORD:ZCOORD,4)

     areaTriangle2 = 0.5_RFREAL*ABS( sNormal(XCOORD)*(v1(2)*v2(3)-v1(3)*v2(2)) &
                                   + sNormal(YCOORD)*(v1(3)*v2(1)-v1(1)*v2(3)) &
                                   + sNormal(ZCOORD)*(v1(1)*v2(2)-v1(2)*v2(1)) )

! Select which triangle to select a point in ----------------------------------

     xRand = Rand1Uniform(region%randData)
     useTriangle1 = ( (areaTriangle1+areaTriangle2)*xRand < areaTriangle1 )

! Select a random point within the appropriate triangle -----------------------

     xRand = Rand1Uniform(region%randData)
     yRand = Rand1Uniform(region%randData)

! - reflect back into triangle ------------------------------------------------

     IF( xRand+yRand  > 1.0_RFREAL ) THEN
       xRand = 1.0_RFREAL-xRand
       yRand = 1.0_RFREAL-yRand
     ENDIF ! xRand

     zRand = 1.0_RFREAL -(xRand+yRand)

     IF ( useTriangle1 ) THEN
       posPlag(XCOORD:ZCOORD) = xRand*pNode(XCOORD:ZCOORD,1) &
                              + yRand*pNode(XCOORD:ZCOORD,2) &
                              + zRand*pNode(XCOORD:ZCOORD,3)
     ELSE
       posPlag(XCOORD:ZCOORD) = xRand*pNode(XCOORD:ZCOORD,4) &
                              + yRand*pNode(XCOORD:ZCOORD,2) &
                              + zRand*pNode(XCOORD:ZCOORD,3)
     ENDIF ! useTriangle1

! Adjust posPlag to be on tile surface ----------------------------------------

     posPlag(XCOORD:ZCOORD) = posPlag(XCOORD:ZCOORD)                   &
                            - DOT_PRODUCT(     posPlag(XCOORD:ZCOORD)- &
                                          faceCentroid(XCOORD:ZCOORD), &
                                               sNormal(XCOORD:ZCOORD)) &
                            * sNormal(XCOORD:ZCOORD)

! Add tiny offset to put posPlag inside cell

     posPlag(XCOORD:ZCOORD) = posPlag(XCOORD:ZCOORD) &
                            + HEIGHT_FRACTION*cellHeight*sNormal(XCOORD:ZCOORD)

   END SUBROUTINE PLAG_InjcSetPositions
!******************************************************************************

   SUBROUTINE PLAG_InjcTestCell( region,iNOff,ijNOff,i,j,k,posPlag,cellLocate )

!... parameters
     TYPE(t_region) :: region

     INTEGER,      INTENT(IN)  :: iNOff,ijNOff,i,j,k

     LOGICAL,      INTENT(OUT) :: cellLocate

     REAL(RFREAL), DIMENSION(3), INTENT(IN) :: posPlag

! ... loop variables
     INTEGER :: mbound

!... local variables
     INTEGER :: ijkNR,ijkNRI,ijkNRJ,ijkNRK,nbound
     INTEGER, DIMENSION(6) :: inCellFlag

     REAL(RFREAL)               :: dpFace, rsgn
     REAL(RFREAL), DIMENSION(3) :: diffPos, faceCentroid, pSFace
     REAL(RFREAL), POINTER, DIMENSION(:,:) :: pSNormal

!******************************************************************************

! Get dimensions --------------------------------------------------------------

     nbound = 6
     inCellFlag(:) = 0
     cellLocate = .FALSE.

     ijkNR   = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
     ijkNRI  = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
     ijkNRJ  = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
     ijkNRK  = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)

! Loop over all cell faces ----------------------------------------------------

     DO mbound = 1, nBound

       inCellFlag(mbound) = 0

       SELECT CASE (mbound)

! - i-face check --------------------------------------------------------------

         CASE(1)
           pSNormal    => region%levels(iLev)%plag%si

           rsgn        = -1._RFREAL
           pSFace(1:3) = (/rsgn*pSNormal(XCOORD,ijkNR), &
                           rsgn*pSNormal(YCOORD,ijkNR), &
                           rsgn*pSNormal(ZCOORD,ijkNR)/)
           faceCentroid(1:3) = (/pFc(XCOORD,ICOORD,ijkNR), &
                                 pFc(YCOORD,ICOORD,ijkNR), &
                                 pFc(ZCOORD,ICOORD,ijkNR)/)
         CASE(2)
           pSNormal    => region%levels(iLev)%plag%si

           rsgn        = +1._RFREAL
           pSFace(1:3) = (/rsgn*pSNormal(XCOORD,ijkNRI), &
                           rsgn*pSNormal(YCOORD,ijkNRI), &
                           rsgn*pSNormal(ZCOORD,ijkNRI)/)
           faceCentroid(1:3) = (/pFc(XCOORD,ICOORD,ijkNRI), &
                                 pFc(YCOORD,ICOORD,ijkNRI), &
                                 pFc(ZCOORD,ICOORD,ijkNRI)/)

! - j-face check --------------------------------------------------------------

         CASE(3)
           pSNormal    => region%levels(iLev)%plag%sj

           rsgn        = -1._RFREAL
           pSFace(1:3) = (/rsgn*pSNormal(XCOORD,ijkNR), &
                           rsgn*pSNormal(YCOORD,ijkNR), &
                           rsgn*pSNormal(ZCOORD,ijkNR)/)
           faceCentroid(1:3) = (/pFc(XCOORD,JCOORD,ijkNR), &
                                 pFc(YCOORD,JCOORD,ijkNR), &
                                 pFc(ZCOORD,JCOORD,ijkNR)/)

         CASE(4)
           pSNormal    => region%levels(iLev)%plag%sj

           rsgn        = +1._RFREAL
           pSFace(1:3) = (/rsgn*pSNormal(XCOORD,ijkNRJ), &
                           rsgn*pSNormal(YCOORD,ijkNRJ), &
                           rsgn*pSNormal(ZCOORD,ijkNRJ)/)
           faceCentroid(1:3) = (/pFc(XCOORD,JCOORD,ijkNRJ), &
                                 pFc(YCOORD,JCOORD,ijkNRJ), &
                                 pFc(ZCOORD,JCOORD,ijkNRJ)/)

! - k-face check -------------------------------------------------------------

         CASE(5)
           pSNormal    => region%levels(iLev)%plag%sk

           rsgn        = -1._RFREAL
           pSFace(1:3) = (/rsgn*pSNormal(XCOORD,ijkNR), &
                           rsgn*pSNormal(YCOORD,ijkNR), &
                           rsgn*pSNormal(ZCOORD,ijkNR)/)
           faceCentroid(1:3) = (/pFc(XCOORD,KCOORD,ijkNR), &
                                 pFc(YCOORD,KCOORD,ijkNR), &
                                 pFc(ZCOORD,KCOORD,ijkNR)/)

         CASE(6)
           pSNormal    => region%levels(iLev)%plag%sk

           rsgn        = +1._RFREAL
           pSFace(1:3) = (/rsgn*pSNormal(XCOORD,ijkNRK), &
                           rsgn*pSNormal(YCOORD,ijkNRK), &
                           rsgn*pSNormal(ZCOORD,ijkNRK)/)
           faceCentroid(1:3) = (/pFc(XCOORD,KCOORD,ijkNRK), &
                                 pFc(YCOORD,KCOORD,ijkNRK), &
                                 pFc(ZCOORD,KCOORD,ijkNRK)/)

       END SELECT ! mbound

! - Compute position vector difference ----------------------------------------
!   and perform dot product with face vectors ---------------------------------
!   need to check (r_p - r_fc). n_fc > 0 --------------------------------------

       diffPos(1:3) = posPlag(1:3)-faceCentroid(1:3)
       dpFace = DOT_PRODUCT( pSFace,diffPos )

       IF ( dpFace >= 0.0_RFREAL ) inCellFlag(mbound) = 1

     ENDDO ! mbound

! Perform test for incell location --------------------------------------------

     IF ( SUM( inCellFlag(1:6) ) == 6 ) cellLocate = .TRUE.

   END SUBROUTINE PLAG_InjcTestCell
!******************************************************************************


END SUBROUTINE PLAG_InjcEjectParticle

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InjcEjectParticle.F90,v $
! Revision 1.8  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2007/03/06 23:13:13  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.5  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.4  2005/02/03 23:24:04  fnajjar
! Added error trap if nPcls exceeds nPclsTot
!
! Revision 1.3  2005/01/21 14:39:06  fnajjar
! Bug fix to invoke region instead of pRegion
!
! Revision 1.2  2005/01/20 15:36:19  fnajjar
! Bug fix to properly bypass negative pool volume
!
! Revision 1.1  2004/12/01 20:57:40  fnajjar
! Initial revision after changing case
!
! Revision 1.18  2004/07/23 22:43:17  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.17  2004/06/30 15:43:57  fnajjar
! Included kernel for CRE model
!
! Revision 1.16  2004/06/16 23:06:33  fnajjar
! Renamed variabled for CRE kernel
!
! Revision 1.15  2004/04/09 23:12:11  fnajjar
! Added status to injected particle datastructure
!
! Revision 1.14  2004/03/25 21:16:06  jferry
! made initial BurnStatus depend on whether burning interaction is used
!
! Revision 1.13  2004/03/02 21:47:30  jferry
! Added After Update interactions
!
! Revision 1.12  2003/11/21 22:35:51  fnajjar
! Update Random Number Generator
!
! Revision 1.11  2003/11/03 21:21:51  fnajjar
! Changed definition of face vectors pointing to PLAG datastructure
!
! Revision 1.10  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.9  2003/05/06 23:43:58  fnajjar
! Commented off I/O and moved it to PLAG_appendDataFromBuffers
!
! Revision 1.8  2003/04/18 23:13:34  fnajjar
! Bug fix for incorrect definition of normal vector to pSFace
!
! Revision 1.7  2003/04/18 19:23:18  fnajjar
! Redefined dpFace to be a true dot product using FORTRAN90 intrinisic
!
! Revision 1.6  2003/04/16 22:56:38  fnajjar
! Included improved positioning kernel for complex geometry
!
! Revision 1.3  2003/02/04 19:06:46  f-najjar
! Commented write statement for out-of-range tiles
!
! Revision 1.2  2003/01/10 19:22:26  f-najjar
! Included iReg in calling sequence
!
! Revision 1.1  2002/10/25 14:16:31  f-najjar
! Initial Import of Rocpart
!
!******************************************************************************







