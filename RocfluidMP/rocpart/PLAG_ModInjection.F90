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
! ******************************************************************************
!
! Purpose: Suite of routines for stochastic injection algorithm.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: PLAG_ModInjection.F90,v 1.11 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_ModInjection
  
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModDataTypes
  USE ModGlobal,     ONLY: t_global
  USE ModGrid,       ONLY: t_grid
  USE ModPartLag,    ONLY: t_plag, t_plag_input, t_tile_plag
  USE ModBndPatch,   ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModRandom,     ONLY: Rand1Uniform

  USE PLAG_ModParameters
  USE INRT_ModParameters

  USE PLAG_ModInterfaces, ONLY : PLAG_InjcMakeParticle, &
                                 PLAG_InjcSetInjection

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: PLAG_EjectParticle, &
            PLAG_InvokeEjecModel1, &
            PLAG_InvokeConsRandEjec
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: PLAG_ModInjection.F90,v $ $Revision: 1.11 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS



!******************************************************************************
!
! Purpose: Fills PLAG datastructure for ejected particle.
!
! Description: none.
!
! Input:
!   pGlobal      Pointer to global data
!   pPlag        Pointer to particle data
!   pTilePlag    Pointer to tile particle data
!   iTile        Tile index
!   icg          Cell index
!   burnStat     Particle burning status
!   nVert        Number of vertices on tile element surface
!   normDirFlag  Normal direction flag
!   xyzVertex    Coordinates of element vertices
!   fn           Face normal coordinates
!   fc           Face centroids coordinates
!   cellHeight   Cell height
!   plagVolRatio Volume ratio particles get from tiles
!
! Output:
!   pPlag      Plag values for cv, aiv, arv of current region.
!   pTilePlag  Tile values for cv, dv of current region.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE PLAG_EjectParticle( pRegion,pPlag,pTilePlag,      &
                               iTile,icg,burnStat,           &
                               nVert,normDirFlag,xyzVertex,  &
                               fn,fc,cellHeight,plagVolRatio )

#ifdef RFLU
  USE RFLU_ModInCellTest, ONLY: RFLU_ICT_TestInCell
#endif

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: burnStat,icg,iTile,normDirFlag,nVert

  REAL(RFREAL),                      INTENT(IN) :: cellHeight,plagVolRatio
  REAL(RFREAL), DIMENSION(ZCOORD),   INTENT(IN) :: fc,fn
  REAL(RFREAL), DIMENSION(ZCOORD,4), INTENT(IN) :: xyzVertex     
   
  TYPE(t_plag),      POINTER :: pPlag
  TYPE(t_tile_plag), POINTER :: pTilePlag
  TYPE(t_region), POINTER  :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER, PARAMETER :: ITER_CELL_LOCATE_MAX = 20
  INTEGER :: iCont,iCvPlagMass,iCvTileMass,iPcl,iReg,iterCellLocate, &
             nCont,nextIdNumber,nPcls,nPclsMax
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass, pCvTileMass
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv

  LOGICAL :: cellLocate
   
  REAL(RFREAL) :: plagMomeNrm,sxn,syn,szn,xLoc,yLoc,zLoc
  REAL(RFREAL), DIMENSION(3) :: poolxyz
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv,pCvPlag,pCvTile,pDvTile
  
  TYPE(t_global),    POINTER :: global 
! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ModInjection.F90,v $ $Revision: 1.11 $'

  global => pRegion%global

  CALL RegisterFunction( global, 'PLAG_EjectParticle',&
  'PLAG_ModInjection.F90' )

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pCvPlag     => pPlag%cv
  pCvPlagMass => pPlag%cvPlagMass
  pAiv        => pPlag%aiv
  pArv        => pPlag%arv

  pCvTile     => pTilePlag%cv
  pCvTileMass => pTilePlag%cvTileMass
  pDvTile     => pTilePlag%dv

  iReg     = pRegion%iRegionGlobal
  nCont    = pRegion%plagInput%nCont

! ******************************************************************************
! Load normals to point inwards
!   Note: RFLU defines normals pointing outwards.
! ******************************************************************************

  sxn = fn(XCOORD) *normDirFlag
  syn = fn(YCOORD) *normDirFlag
  szn = fn(ZCOORD) *normDirFlag

! ******************************************************************************
! Increment PLAG datastructure for new injected particle
! ******************************************************************************

  pPlag%nPcls = pPlag%nPcls + 1
  iPcl        = pPlag%nPcls

  pTilePlag%nPclsInjc(iTile) = pTilePlag%nPclsInjc(iTile) + 1

! ******************************************************************************
! Trap error if nPcls exceeds maximum datastructure dimension, nPclsMax
! ******************************************************************************

  nPcls    = pPlag%nPcls
  nPclsMax = pRegion%plag%nPclsMax 

  IF ( nPcls >= nPclsMax ) THEN
    WRITE(STDOUT,*) ' PLAG_EjectParticle: Datastructure Dimension Exceeded ',&
                      nPcls,nPclsMax
    CALL ErrorStop( global,ERR_PLAG_MEMOVERFLOW,__LINE__ )
  END IF ! nPcls

! ******************************************************************************
! Keep track of global count in particle ejected in region
! ******************************************************************************

  pPlag%nextIdNumber = pPlag%nextIdNumber + 1
  nextIdNumber       = pPlag%nextIdNumber

! ******************************************************************************
! Compute PLAG cv datastructure for new injected particle
! ******************************************************************************

  DO iCont = 1, nCont
    iCvPlagMass = pCvPlagMass(iCont)
    iCvTileMass = pCvTileMass(iCont)
    pCvPlag(iCvPlagMass,iPcl) = pCvTile(iCvTileMass,iTile)    *plagVolRatio
  END DO ! iCont

  pCvPlag(CV_PLAG_ENER,iPcl)  = pCvTile(CV_TILE_ENER,  iTile) *plagVolRatio

  plagMomeNrm                 = pCvTile(CV_TILE_MOMNRM,iTile) *plagVolRatio
  pCvPlag(CV_PLAG_XMOM,iPcl)  = plagMomeNrm*sxn
  pCvPlag(CV_PLAG_YMOM,iPcl)  = plagMomeNrm*syn
  pCvPlag(CV_PLAG_ZMOM,iPcl)  = plagMomeNrm*szn

  pCvPlag(CV_PLAG_ENERVAPOR,iPcl) = 0._RFREAL

! ******************************************************************************
! Determine particle position
! include a test to insure particle is in tile cell
!  Note: RFLU defines normals pointing outwards.
! ******************************************************************************

  iterCellLocate = 1

199 CONTINUE

  CALL PLAG_InjcSetPoolPos( pRegion,nVert,normDirFlag,xyzVertex, &
                            fc(XCOORD:ZCOORD),fn(XCOORD:ZCOORD), &
                            cellHeight,poolxyz )

  xLoc = poolxyz(XCOORD)
  yLoc = poolxyz(YCOORD)
  zLoc = poolxyz(ZCOORD)

  cellLocate = RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,icg)

! TEMPORARY
!   PRINT*,' nVert       = ', nVert
!   PRINT*,' icg         = ' ,icg
!   PRINT*,' iTile       = ' ,iTile
!   PRINT*,' normDirFlag = ', normDirFlag
!   PRINT*,' fc          = ', fc(XCOORD:ZCOORD)
!   PRINT*,' fn          = ', fn(XCOORD:ZCOORD)
!   PRINT*,' cellHeight  = ', cellHeight
!   PRINT*,' poolxyz     = ', poolxyz
!   PRINT*,' cellLocate  = ', cellLocate
! END TEMPORARY

  IF (cellLocate .EQV. .FALSE.) THEN
    iterCellLocate = iterCellLocate+1
    IF ( iterCellLocate > ITER_CELL_LOCATE_MAX ) THEN
       WRITE(STDOUT,*) ' PLAG_EjectParticle: Unable to create a particle after ',&
                         ITER_CELL_LOCATE_MAX, ' iterations'

       CALL ErrorStop( global,ERR_PLAG_CELLINDEX,__LINE__ )
    END IF ! iterCellLocate
    GOTO 199
  ENDIF ! cellLocate

  pCvPlag(CV_PLAG_XPOS,iPcl) = poolxyz(XCOORD)
  pCvPlag(CV_PLAG_YPOS,iPcl) = poolxyz(YCOORD)
  pCvPlag(CV_PLAG_ZPOS,iPcl) = poolxyz(ZCOORD)

! ******************************************************************************
! Set aiv and arv
! ******************************************************************************

  pAiv(AIV_PLAG_PIDINI,iPcl) = nextIdNumber
  pAiv(AIV_PLAG_REGINI,iPcl) = iReg
  pAiv(AIV_PLAG_ICELLS,iPcl) = icg
  pAiv(AIV_PLAG_STATUS,iPcl) = PLAG_STATUS_KEEP
  pAiv(AIV_PLAG_BURNSTAT,iPcl) = burnStat

  pArv(ARV_PLAG_SPLOAD,iPcl) = pDvTile(DV_TILE_SPLOAD,iTile)

! ******************************************************************************
! Decrease pool variables
! ******************************************************************************

  DO iCont = 1, nCont
    iCvPlagMass = pCvPlagMass(iCont)
    iCvTileMass = pCvTileMass(iCont)
    pCvTile(iCvTileMass,iTile) = pCvTile(iCvTileMass,   iTile) -  &
                                 pDvTile(DV_TILE_SPLOAD,iTile) *  &
                                 pCvPlag(iCvPlagMass,iPcl)
  END DO ! iCont

  pCvTile(CV_TILE_MOMNRM,iTile) = pCvTile(CV_TILE_MOMNRM,iTile) - &
                                  pDvTile(DV_TILE_SPLOAD,iTile) * &
                                  plagMomeNrm

  pCvTile(CV_TILE_ENER  ,iTile) = pCvTile(CV_TILE_ENER  ,iTile) - &
                                  pDvTile(DV_TILE_SPLOAD,iTile) * &
                                  pCvPlag(CV_PLAG_ENER,iPcl)

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_EjectParticle




! ******************************************************************************
!
! Purpose: Determine coordinates for injected particle on patch surface
!   based on random number generator.
!
! Description: None.
!
! Input:
!   pRegion       Pointer to region data
!   nVert         Number of vertices
!   normDirFlag   Flag setting normal direction
!   xyzVertex     Coordinates of patch
!   faceCentroid  Coordinates of face centroid
!   sNormal       Vector of face normal
!   cellHeight    Height of cell to which patch is attached
!
! Output:
!   posPlag       Coordinates of injected particle
!
! Notes:
!   1. Normals need to be set pointing inward
!   2. normDirFlag is set to -1 for RFLU and +1 for RFLO
!   3. RFLU always set normals pointing outwards of the computational domain
!   4. RFLU numbers vertices in the counterclockwise direction as 1-2-3-4
!
! ******************************************************************************

  SUBROUTINE PLAG_InjcSetPoolPos( pRegion,nVert,normDirFlag,xyzVertex,    &
                                  faceCentroid,sNormal,cellHeight,posPlag )
    
    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================
    
    INTEGER,      INTENT(IN) :: nVert,normDirFlag
       
    REAL(RFREAL), INTENT(IN) :: cellHeight
    REAL(RFREAL), DIMENSION(ZCOORD),   INTENT(IN)  :: faceCentroid,sNormal
    REAL(RFREAL), DIMENSION(ZCOORD,4), INTENT(IN)  :: xyzVertex
    REAL(RFREAL), DIMENSION(ZCOORD),   INTENT(OUT) :: posPlag
     
    TYPE(t_region), POINTER  :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    LOGICAL :: useTriangle1
     
    REAL(RFREAL), PARAMETER :: HEIGHT_FRACTION = 1.0E-3_RFREAL
    REAL(RFREAL) :: areaTriangle1,areaTriangle2,xRand,xTriRand,yRand,zRand
    REAL(RFREAL), DIMENSION(ZCOORD) :: sNormalInward, v1,v2  

    TYPE(t_global), POINTER :: global
      
! ******************************************************************************
!   Start
! ******************************************************************************

    RCSIdentString = '$RCSfile: PLAG_ModInjection.F90,v $ $Revision: 1.11 $'

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_InjcSetPoolPos',&
  'PLAG_ModInjection.F90' )

! ****************************************************************************** 
!   Set normals to point inwards 
! ****************************************************************************** 

    sNormalInward(XCOORD:ZCOORD) = normDirFlag  *sNormal(XCOORD:ZCOORD) 

! ****************************************************************************** 
!   Select a random point within the appropriate triangle 
! ****************************************************************************** 

    xTriRand = Rand1Uniform(pRegion%randData)   
    xRand    = Rand1Uniform(pRegion%randData) 
    yRand    = Rand1Uniform(pRegion%randData) 

! ******************************************************************************     
!   Reflect back into triangle
! ******************************************************************************

    IF ( xRand+yRand  > 1.0_RFREAL ) THEN    
      xRand = 1.0_RFREAL-xRand
      yRand = 1.0_RFREAL-yRand
    ENDIF ! xRand
     
    zRand = 1.0_RFREAL -(xRand+yRand)

! ******************************************************************************  
!   Setup triangles depending on element type
! ******************************************************************************  

    SELECT CASE(nVert)     

! ==============================================================================
!     For triangular elements, always select triangle1
! ==============================================================================

      CASE (3)
        useTriangle1 = .TRUE.        

! ==============================================================================
!     Partition the face into two triangles and compute the areas 
!     of their projections normal to sNormal
! ==============================================================================

      CASE (4)
        v1(XCOORD:ZCOORD) = xyzVertex(XCOORD:ZCOORD,2) -xyzVertex(XCOORD:ZCOORD,1)
        v2(XCOORD:ZCOORD) = xyzVertex(XCOORD:ZCOORD,3) -xyzVertex(XCOORD:ZCOORD,1)
    
        areaTriangle1 = 0.5_RFREAL*ABS( sNormal(XCOORD)*(v1(2)*v2(3)-v1(3)*v2(2)) &
                                      + sNormal(YCOORD)*(v1(3)*v2(1)-v1(1)*v2(3)) &
                                      + sNormal(ZCOORD)*(v1(1)*v2(2)-v1(2)*v2(1)) )  

        v1(XCOORD:ZCOORD) = xyzVertex(XCOORD:ZCOORD,3) -xyzVertex(XCOORD:ZCOORD,1)
        v2(XCOORD:ZCOORD) = xyzVertex(XCOORD:ZCOORD,4) -xyzVertex(XCOORD:ZCOORD,1)

        areaTriangle2 = 0.5_RFREAL*ABS( sNormal(XCOORD)*(v1(2)*v2(3)-v1(3)*v2(2)) &
                                      + sNormal(YCOORD)*(v1(3)*v2(1)-v1(1)*v2(3)) &
                                      + sNormal(ZCOORD)*(v1(1)*v2(2)-v1(2)*v2(1)) )

! ------------------------------------------------------------------------------
!     Select which triangle to select a point in 
! ------------------------------------------------------------------------------

        useTriangle1 = ( (areaTriangle1+areaTriangle2)*xTriRand < areaTriangle1 ) 
    END SELECT ! nVert 

! ******************************************************************************  
!   Set particle position on triangular tile
! ****************************************************************************** 
     
    IF ( useTriangle1 ) THEN
      posPlag(XCOORD:ZCOORD) = xRand*xyzVertex(XCOORD:ZCOORD,1) &
                             + yRand*xyzVertex(XCOORD:ZCOORD,2) &
                             + zRand*xyzVertex(XCOORD:ZCOORD,3) 
    ELSE
      posPlag(XCOORD:ZCOORD) = xRand*xyzVertex(XCOORD:ZCOORD,1) &
                             + yRand*xyzVertex(XCOORD:ZCOORD,3) &
                             + zRand*xyzVertex(XCOORD:ZCOORD,4) 
    ENDIF ! useTriangle1

! ****************************************************************************** 
!   Adjust posPlag to be on tile surface      
! ******************************************************************************     

    posPlag(XCOORD:ZCOORD) = posPlag(XCOORD:ZCOORD)                   &
                           - DOT_PRODUCT(     posPlag(XCOORD:ZCOORD)- &
                                         faceCentroid(XCOORD:ZCOORD), &
                                        sNormalInward(XCOORD:ZCOORD)) &
                           * sNormalInward(XCOORD:ZCOORD) 

! ******************************************************************************     
!   Add tiny offset to put posPlag inside cell
! ******************************************************************************      

    posPlag(XCOORD:ZCOORD) = posPlag(XCOORD:ZCOORD)      &
                           + HEIGHT_FRACTION*cellHeight  &
                           *sNormalInward(XCOORD:ZCOORD) 

! *****************************************************************************
!   End
! *****************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_InjcSetPoolPos



!******************************************************************************
!
! Purpose: Invoke ejection model1 to generate particles on tile surfaces.
!
! Description: none.
!
! Input:
!   pRegion      Pointer to region data
!   pPlag        Pointer to particle data
!   pTilePlag    Pointer to tile particle data
!   iTile        Tile identifier
!   icg          Cell index
!   burnStat     Particle burning status
!   nVert        Number of vertices on tile element surface
!   normDirFlag  Normal direction flag
!   xyzVertex    Coordinates of element vertices
!   volMeanPart  Mean particle volume
!   fn           Face normal coordinates
!   fc           Face centroids coordinates
!   cellHeight   Cell height
!
! Output:
!   pPlag      Plag values for cv, aiv, arv of current region.
!   pTilePlag  Tile values for cv, dv of current region.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE PLAG_InvokeEjecModel1( pRegion,pPlag,pTilePlag,      &
                                  iTile,icg,burnStat,nVert,     &
                                  normDirFlag,xyzVertex,        &
                                  volMeanPart,fn,fc,cellHeight  )

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: burnStat,icg,iTile,normDirFlag,nVert

  REAL(RFREAL),                      INTENT(IN) :: cellHeight,volMeanPart
  REAL(RFREAL), DIMENSION(ZCOORD),   INTENT(IN) :: fc,fn
  REAL(RFREAL), DIMENSION(ZCOORD,4), INTENT(IN) :: xyzVertex     

  TYPE(t_region),    POINTER :: pRegion  
  TYPE(t_plag),      POINTER :: pPlag
  TYPE(t_tile_plag), POINTER :: pTilePlag

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: injcDiamDist
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass, pCvTileMass

  LOGICAL :: injectQ

  REAL(RFREAL) :: injcBeta,pi,plagVolRatio,poolVolSum,tInjcCoeff,tInjcSum
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCvTile,pDvTile

  TYPE(t_global),    POINTER :: global 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ModInjection.F90,v $ $Revision: 1.11 $'

  global => pRegion%global

  CALL RegisterFunction( global, 'PLAG_InvokeEjecModel1',&
  'PLAG_ModInjection.F90' )

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pCvTile     => pTilePlag%cv
  pCvTileMass => pTilePlag%cvTileMass
  pDvTile     => pTilePlag%dv

  injcDiamDist = pRegion%plagInput%injcDiamDist
  injcBeta     = pRegion%plagInput%injcBeta

! ******************************************************************************
! Begin injection algorithm
! ******************************************************************************
  
  poolVolSum = SUM ( pCvTile( pCvTileMass(:),iTile ) / &
                     pRegion%plagInput%dens(:) )

  pDvTile(DV_TILE_POOLVOLD,iTile) = poolVolSum

  tInjcCoeff = injcBeta *volMeanPart

  tInjcSum   = 0.0_RFREAL

  CALL PLAG_InjcSetInjection( pRegion,pTilePlag,iTile,        &
                              tInjcCoeff,tInjcSum,poolVolSum, &
                              injectQ,plagVolRatio            )

! ******************************************************************************
! If injectQ set to True, PLAG datastructure for new injected particle
! ******************************************************************************

  DO WHILE ( injectQ )

! TEMPORARY
!   PRINT*,'Ejection Model =', pRegion%plagInput%ejecModel
!   PRINT*,' IFC = ' ,iTile
!   PRINT*,'injectQ     = ',injectQ
!   PRINT*,' plagVolRatio = ',plagVolRatio
!   PRINT*,' tInjcSum   = ',tInjcSum
!   PRINT*,' poolVolSum         = ',poolVolSum
!   PRINT*,' pDvTile(DV_TILE_DIAM,iTile)        = ',pDvTile(DV_TILE_DIAM,iTile)
!   PRINT*,' pDvTile(DV_TILE_COUNTDOWN,iTile)   = ',pDvTile(DV_TILE_COUNTDOWN,iTile)
! END TEMPORARY

! ==============================================================================
!  Eject particle
! ==============================================================================

    CALL PLAG_EjectParticle( pRegion,pPlag,pTilePlag,       &
                             iTile,icg,burnStat,            &
                             nVert,normDirFlag,xyzVertex,   &
                             fn,fc,cellHeight,plagVolRatio  )

! ==============================================================================
!   Update tile diameter and superparticle loading factor
! ==============================================================================

    CALL PLAG_InjcMakeParticle( pRegion, injcDiamDist,         &
                                pDvTile(DV_TILE_DIAM,  iTile), &
                                pDvTile(DV_TILE_SPLOAD,iTile)  )

! ==============================================================================
!   Check if another particle could be injected
! ==============================================================================

    CALL PLAG_InjcSetInjection( pRegion,pTilePlag,iTile,        &
                                tInjcCoeff,tInjcSum,poolVolSum, &
                                injectQ,plagVolRatio            )

  END DO ! injectQ

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_InvokeEjecModel1



!******************************************************************************
!
! Purpose: Invoke conservative random ejection model
!          to generate particles on tile surfaces.
!
! Description: none.
!
! Input:
!   pRegion      Pointer to region data
!   pPlag        Pointer to particle data
!   pTilePlag    Pointer to tile particle data
!   iTile        Tile identifier
!   icg          Cell index
!   burnStat     Particle burning status
!   nVert        Number of vertices on tile element surface
!   normDirFlag  Normal direction flag
!   xyzVertex    Coordinates of element vertices
!   volMeanPart  Mean particle volume
!   fn           Face normal coordinates
!   fc           Face centroids coordinates
!   cellHeight   Cell height
!
! Output:
!   pPlag      Plag values for cv, aiv, arv of current region.
!   pTilePlag  Tile values for cv, dv of current region.
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE PLAG_InvokeConsRandEjec( pRegion,pPlag,pTilePlag,      &
                                    iTile,icg,burnStat,nVert,     &
                                    normDirFlag,xyzVertex,        &
                                    volMeanPart,fn,fc,cellHeight  )

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: burnStat,icg,iTile,normDirFlag,nVert

  REAL(RFREAL),                      INTENT(IN) :: cellHeight,volMeanPart
  REAL(RFREAL), DIMENSION(ZCOORD),   INTENT(IN) :: fc,fn
  REAL(RFREAL), DIMENSION(ZCOORD,4), INTENT(IN) :: xyzVertex     

  TYPE(t_region),    POINTER :: pRegion  
  TYPE(t_plag),      POINTER :: pPlag
  TYPE(t_tile_plag), POINTER :: pTilePlag

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: ejecModel,injcDiamDist
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass, pCvTileMass

  REAL(RFREAL) :: injcBeta,pi,plagVolRatio,poolVolSum,tInjcCoeff,tInjcSum
  REAL(RFREAL) :: meanSuperParticleVolume, spload
  REAL(RFREAL) :: injcBetaFac, injcBetaFacInv
  REAL(RFREAL) :: poolVolOldSum, poolVolFinal, remainingVolume
  REAL(RFREAL) :: currentSuperParticleVolume, poolVolume
  REAL(RFREAL) :: poolExcess, pExcessS, possibleExcess
  REAL(RFREAL) :: countdown, countdownNext, deltaVolume, randUnif
  REAL(RFREAL) :: currentParticleVolume
  REAL(RFREAL) :: poolVolCurr
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCvTile,pCvOldTile,pDvTile

  TYPE(t_global),    POINTER :: global 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ModInjection.F90,v $ $Revision: 1.11 $'

  global => pRegion%global

  CALL RegisterFunction( global, 'PLAG_InvokeConsRandEjec',&
  'PLAG_ModInjection.F90' )

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pCvTile     => pTilePlag%cv
  pCvOldTile  => pTilePlag%cvOld
  pCvTileMass => pTilePlag%cvTileMass
  pDvTile     => pTilePlag%dv

 
  injcDiamDist = pRegion%plagInput%injcDiamDist
  injcBeta     = pRegion%plagInput%injcBeta
  
  pi        = global%pi
  spload       = pRegion%plagInput%spload
  meanSuperParticleVolume = volMeanPart *spLoad

  injcBetaFac = 2.0_RFREAL *injcBeta *meanSuperParticleVolume**2
  injcBetaFacInv = 1.0_RFREAL/injcBetaFac

! ******************************************************************************
! Begin injection algorithm
! ******************************************************************************

! ==============================================================================
! poolVolSum is the pool volume after the current timestep is finished
! ==============================================================================

  poolVolSum    = SUM ( pCvTile( pCvTileMass(:),iTile ) / &
                        pRegion%plagInput%dens(:) )

! ==============================================================================
! poolVolOldSum is the pool volume before the current timestep began
! ==============================================================================

  poolVolOldSum = SUM ( pCvOldTile( pCvTileMass(:),iTile ) / &
                        pRegion%plagInput%dens(:) )

! ==============================================================================
! poolVolume is built up incrementally with volume injection,
!  and is decremented when (super)particles are ejected
! ==============================================================================

  poolVolume   = poolVolOldSum

! ==============================================================================
! poolVolFinal starts with all volume injection having occurred,
!  and is decremented when (super)particles are ejected
! ==============================================================================

  poolVolFinal = poolVolSum

! ==============================================================================
! poolVolCurr starts with all volume injection having occurred,
!  and is decremented when (super)particles are ejected
! ==============================================================================

  poolVolCurr = poolVolSum

! ==============================================================================
! remainingVolume is the volume to be added to the pool over the
!  current time step: it decreases to 0 as all volume is added
! ==============================================================================

  remainingVolume = poolVolSum -poolVolOldSum

! TEMPORARY
!   PRINT*,'Ejection Model =',  pRegion%plagInput%ejecModel
!   PRINT*,' IFC = ' ,iTile
!   PRINT*,' plagVolRatio = ',plagVolRatio
!   PRINT*,' poolVolSum         = ',poolVolSum
!   PRINT*,' remainingVolume    = ',remainingVolume
!   PRINT*,' pDvTile(DV_TILE_DIAM,iTile)        = ',pDvTile(DV_TILE_DIAM,iTile)
!   PRINT*,' pDvTile(DV_TILE_COUNTDOWN,iTile)   = ',pDvTile(DV_TILE_COUNTDOWN,iTile)
! END TEMPORARY

! ==============================================================================
! Start the ejection loop for CRE model
! ==============================================================================

  creEjectLoop: DO
    currentSuperParticleVolume = pi/6.0_RFREAL *spLoad &
                               * pDvTile(DV_TILE_DIAM,iTile)**3.0_RFREAL

    currentParticleVolume = pi/6.0_RFREAL &
                          * pDvTile(DV_TILE_DIAM,iTile)**3.0_RFREAL

! ------------------------------------------------------------------------------
!   poolExcess is the volume that would remain in the pool after ejection
! ------------------------------------------------------------------------------

!    poolExcess = poolVolume -currentSuperParticleVolume

! TEMPORARY
    poolExcess = poolVolCurr -currentSuperParticleVolume
! END TEMPORARY

! ------------------------------------------------------------------------------
!   possibleExcess is the volume that would remain in the pool at the end of
!     the timestep after the current particle is (and no others are) ejected
! ------------------------------------------------------------------------------

    possibleExcess = poolExcess + remainingVolume

! ------------------------------------------------------------------------------
!   pExcessS is an auxiliary variable that occurs in a few formulas
! ------------------------------------------------------------------------------

    IF ( poolExcess > 0.0_RFREAL ) THEN
      pExcessS = poolExcess**2
    ELSE
      pExcessS = 0.0_RFREAL
    ENDIF ! poolExcess

! ------------------------------------------------------------------------------
!   countdown is the internal ejection clock for the current particle:  it
!    starts to tick down when pool volume > current (super)particle volume
! ------------------------------------------------------------------------------
 
    countdown = pDvTile(DV_TILE_COUNTDOWN,iTile)

! ------------------------------------------------------------------------------
!   countdownNext is what the clock could tick down to within remainingVolume
! ------------------------------------------------------------------------------

    IF ( possibleExcess > 0.0_RFREAL ) THEN
      countdownNext = countdown + injcBetaFacInv &
                    * (pExcessS -possibleExcess**2.0_RFREAL)
    ELSE
      countdownNext = countdown
    ENDIF ! possibleExcess
 
! ------------------------------------------------------------------------------
!   Pool volume too small to eject particle
! ------------------------------------------------------------------------------

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     if the clock cannot tick down to zero within remainingVolume, no more
!       particles will be ejected in this timestep
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IF ( countdownNext> 0.0_RFREAL ) THEN
      poolVolume = poolVolume + remainingVolume ! not needed except as check
      pDvTile(DV_TILE_COUNTDOWN,iTile) = countdownNext
      EXIT creEjectLoop

!-------------------------------------------------------------------------------
!   Pool volume big enough to eject particle
!-------------------------------------------------------------------------------

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     if the clock can tick down below zero, then we need to find when it
!       reaches zero, and then eject the particle
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ELSE

! --- deltaVolume is the volume that must be added to the pool (out of
!     remainingVolume) in order to bring countdown exactly to zero

      deltaVolume = SQRT(injcBetaFac*countdown +pExcessS) -poolExcess

! --- therefore deltaVolume is added to the pool, the countdown hits zero,
!      and the particle is ejected
      
      poolVolume  = poolVolume +deltaVolume -currentSuperParticleVolume

! --- the deltaVolume added to the pool is taken from remainingVolume
      
      remainingVolume = remainingVolume -deltaVolume

! --- plagVolRatio must be in terms of poolVolFinal:  the pool cv variables
!     are expressed in terms of values at the end of the time step, not
!     in terms of the intermediate values like poolVolume

      plagVolRatio = currentParticleVolume / poolVolFinal

! --- poolVolFinal, like the cv variables, is affected by particle ejection
      
      poolVolFinal = poolVolFinal -currentSuperParticleVolume

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Eject particle
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CALL PLAG_EjectParticle( pRegion,pPlag,pTilePlag,       &
                               iTile,icg,burnStat,            &
                               nVert,normDirFlag,xyzVertex,   &
                               fn,fc,cellHeight,plagVolRatio  )

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Update tile diameter and superparticle loading factor
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CALL PLAG_InjcMakeParticle( pRegion, injcDiamDist,         &
                                  pDvTile(DV_TILE_DIAM,  iTile), &
                                  pDvTile(DV_TILE_SPLOAD,iTile)  )

! TEMPORARY
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Compute current pool volume after being decreased in PLAG_EjectParticle
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
       poolVolCurr  = SUM ( pCvTile( pCvTileMass(:),iTile ) / &
                            pRegion%plagInput%dens(:) )

! END TEMPORARY

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Set countdown for new particle
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      randUnif = Rand1Uniform(pRegion%randData)
      
      IF ( randUnif <= 0.0_RFREAL) THEN
! ----- treats randUnif as 1.9e-22
        pDvTile(DV_TILE_COUNTDOWN,iTile) = 50.0_RFREAL
      ELSE
        pDvTile(DV_TILE_COUNTDOWN,iTile) = -LOG(randUnif)
      END IF ! randUnif

    END IF ! countdownNext
  END DO creEjectLoop

!------------------------------------------------------------------------------
! Temporary check: poolVolume and poolVolFinal should now agree
!------------------------------------------------------------------------------

  IF ( ABS(poolVolume - poolVolFinal) / poolVolFinal > 1.0E-10_RFREAL) &
  PRINT*,'CRE check: iTile ', iTile,(poolVolume - poolVolFinal) / poolVolFinal  

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_InvokeConsRandEjec

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE PLAG_ModInjection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModInjection.F90,v $
! Revision 1.11  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2007/03/06 23:15:32  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.8  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.7  2005/12/30 16:26:25  fnajjar
! Bug fix for vertex and triangle defs in posPlag
!
! Revision 1.6  2005/12/24 21:35:49  haselbac
! Adapted to changes in ICT, cosmetics
!
! Revision 1.5  2005/01/20 15:36:19  fnajjar
! Bug fix to properly bypass negative pool volume
!
! Revision 1.4  2004/07/28 18:56:54  fnajjar
! Included error trap for memory overflow
!
! Revision 1.3  2004/07/12 15:52:38  fnajjar
! Included interface and routine for conservative random ejection model
!
! Revision 1.2  2004/07/12 14:43:17  fnajjar
! Added pertinent routines to prepare for CRE ejection model
!
! Revision 1.1  2004/03/08 22:35:37  fnajjar
! Initial import of Injection Module
!
! ******************************************************************************










