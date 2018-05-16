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
! Purpose: main injection driver specific to RFLU
!          for the multiphase injection algorithm.
!
! Description: none.
!
! Input:
!   pRegion    Pointer to region data
!
! Output:
!   pPlag      Plag values for cv, aiv, arv of current region.
!   pTilePlag  Tile values for cv, dv of current region.
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: PLAG_RFLU_InjectionDriver.F90,v 1.14 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_InjectionDriver(pRegion)

  USE ModDataTypes
  USE ModPartLag,    ONLY : t_plag, t_plag_input, t_tile_plag
  USE ModBndPatch,   ONLY : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModGrid,       ONLY : t_grid
  USE ModError
  USE ModParameters

  USE PLAG_ModParameters
  USE INRT_ModParameters
  USE PLAG_ModInjection

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: bcType,burnStat,ejecModel,icg,ifc,injcDiamDist,iPatch,&
             iTile,iv1,iv2,iv3,iv4,nextIdNumber,normDirFlag,nPatches,nVert
  INTEGER :: iFilePlagTile
  INTEGER, POINTER, DIMENSION(:)   :: bf2c,pCvTileMass
  INTEGER, POINTER, DIMENSION(:,:) :: bf2v

  REAL(RFREAL) :: cellHeight,injcBeta,pi,volMeanPart
  REAL(RFREAL) :: injcBetaFac,injcBetaFacInv,meanSuperParticleVolume,spload

  REAL(RFREAL),          DIMENSION(3,4) :: xyzVertex

  REAL(RFREAL), POINTER, DIMENSION(:)   :: pVol
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: fc,fn,pXyz

  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCvTile,pDvTile

  TYPE(t_grid),      POINTER :: pGrid
  TYPE(t_global),    POINTER :: global
  TYPE(t_patch),     POINTER :: pPatch
  TYPE(t_plag),      POINTER :: pPlag
  TYPE(t_tile_plag), POINTER :: pTilePlag

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_InjectionDriver.F90,v $ $Revision: 1.14 $'

  global => pRegion%global

  CALL RegisterFunction( global, 'PLAG_RFLU_InjectionDriver',&
  'PLAG_RFLU_InjectionDriver.F90' )

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

  nPatches = pGrid%nPatches

  pi       = global%pi
  normDirFlag = -1

  injcDiamDist  = pRegion%plagInput%injcDiamDist
  injcBeta      = pRegion%plagInput%injcBeta
  volMeanPart   = pi/6.0_RFREAL * pRegion%plagInput%injcDiamMean**3
  ejecModel     = pRegion%plagInput%ejecModel
  spLoad        = pRegion%plagInput%spLoad

  IF ( ejecModel == PLAG_EJEC_CRE ) THEN
    meanSuperParticleVolume = volMeanPart *spLoad

    injcBetaFac    = 2.0_RFREAL *injcBeta *meanSuperParticleVolume**2
    injcBetaFacInv = 1.0_RFREAL/injcBetaFac
  ENDIF ! ejecModel

  pXyz        => pGrid%xyz

  pVol        => pGrid%vol

! ******************************************************************************
! Initial burning status of particles
! ******************************************************************************

! Currently all particles start off burning if the burning interaction is used

  IF (pRegion%inrtInput%inrts(INRT_TYPE_BURNING)%used) THEN
    burnStat = INRT_BURNSTAT_ON
  ELSE
    burnStat = INRT_BURNSTAT_OFF
  ENDIF

! ******************************************************************************
! Loop over patches
! ******************************************************************************

  DO iPatch=1,nPatches

    pPatch => pRegion%patches(iPatch)

    bcType =  pPatch%bcType

    bf2c   => pPatch%bf2c
    bf2v   => pPatch%bf2v

    fc     => pPatch%fc
    fn     => pPatch%fn

!===============================================================================
!   Select injection boundary condition
!===============================================================================

    IF ( (bcType >= BC_INJECTION .AND. bcType <= BC_INJECTION + BC_RANGE) .OR. &
         (bcType >= BC_INFLOW    .AND. bcType <= BC_INFLOW    + BC_RANGE)      ) THEN 

!-------------------------------------------------------------------------------
!     Set tile pointers
!-------------------------------------------------------------------------------

      pTilePlag => pPatch%tilePlag

      pCvTile     => pTilePlag%cv
      pCvTileMass => pTilePlag%cvTileMass
      pDvTile     => pTilePlag%dv

!-------------------------------------------------------------------------------
!     Loop over face patches
!-------------------------------------------------------------------------------

      DO ifc = 1, pPatch%nBFaces
        iTile = ifc
        icg   = bf2c(ifc)

!-------------------------------------------------------------------------------
!       Set cell height, element type and vertex coordinates
!-------------------------------------------------------------------------------

        cellHeight = pVol(icg)/fn(XYZMAG,ifc)

        iv1 = pPatch%bv(bf2v(1,ifc))
        iv2 = pPatch%bv(bf2v(2,ifc))
        iv3 = pPatch%bv(bf2v(3,ifc))

        xyzVertex(XCOORD:ZCOORD,1) = pXyz(XCOORD:ZCOORD,iv1)
        xyzVertex(XCOORD:ZCOORD,2) = pXyz(XCOORD:ZCOORD,iv2)
        xyzVertex(XCOORD:ZCOORD,3) = pXyz(XCOORD:ZCOORD,iv3)

        IF ( bf2v(4,ifc) == VERT_NONE ) THEN
          nVert = 3
          iv4   = VERT_NONE
          xyzVertex(XCOORD:ZCOORD,4) = 0.0_RFREAL
        ELSE
          nVert = 4
          iv4 = pPatch%bv(bf2v(4,ifc))
          xyzVertex(XCOORD:ZCOORD,4) = pXyz(XCOORD:ZCOORD,iv4)
        ENDIF ! iv4

!-------------------------------------------------------------------------------
!       Invoke injection algorithm
!-------------------------------------------------------------------------------

        SELECT CASE(ejecModel)

!------------------------------------------------------------------------------
!         Ejection Model 1
!------------------------------------------------------------------------------

          CASE(PLAG_EJEC_MODEL1)
            CALL PLAG_InvokeEjecModel1( pRegion,pPlag,pTilePlag,           &
                                        iTile,icg,burnStat,nVert,          &
                                        normDirFlag,xyzVertex,volMeanPart, &
                                        fn(XCOORD:ZCOORD,ifc),             &
                                        fc(XCOORD:ZCOORD,ifc),             &
                                        cellHeight                         )

!------------------------------------------------------------------------------
!         Conservative Random Ejection (CRE) Model
!------------------------------------------------------------------------------

          CASE(PLAG_EJEC_CRE)
            CALL PLAG_InvokeConsRandEjec( pRegion,pPlag,pTilePlag,           &
                                          iTile,icg,burnStat,nVert,          &
                                          normDirFlag,xyzVertex,volMeanPart, &
                                          fn(XCOORD:ZCOORD,ifc),             &
                                          fc(XCOORD:ZCOORD,ifc),             &
                                          cellHeight                         )

!------------------------------------------------------------------------------
!         Default case: trap error
!------------------------------------------------------------------------------

          CASE DEFAULT
            CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )
        END SELECT ! ejecModel      
      END DO ! ifc

!-------------------------------------------------------------------------------
!     Reset injected particle size
!-------------------------------------------------------------------------------

      pTilePlag%nPclsInjc(:) = 0

    END IF ! bcType
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLU_InjectionDriver

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_InjectionDriver.F90,v $
! Revision 1.14  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.13  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.12  2006/09/18 20:36:20  fnajjar
! Activated tile datastructure for inflow bc
!
! Revision 1.11  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.10  2005/12/24 21:38:51  haselbac
! Cosmetics only
!
! Revision 1.9  2004/08/17 21:00:11  fnajjar
! Moved iv4 into IF statement for proper programming definition
!
! Revision 1.8  2004/08/15 20:36:41  fnajjar
! Bug fix for nVert, need to use bf2v instead of iv4
!
! Revision 1.7  2004/07/12 15:52:39  fnajjar
! Included interface and routine for conservative random ejection model
!
! Revision 1.6  2004/07/12 14:42:24  fnajjar
! Created proper interface to prepare for CRE ejection model
!
! Revision 1.5  2004/06/16 23:02:24  fnajjar
! Renamed variables for CRE kernel
!
! Revision 1.4  2004/05/21 23:28:45  fnajjar
! Bug fix for iv to properly access global face-to-vertex mapping and made io-cleanup
!
! Revision 1.3  2004/03/26 21:28:59  fnajjar
! Added aiv status flag at injection and properly activated burn status
!
! Revision 1.2  2004/03/25 21:16:06  jferry
! made initial BurnStatus depend on whether burning interaction is used
!
! Revision 1.1  2004/03/08 22:36:08  fnajjar
! Initial import of RFLU-specific injection routines
!
! ******************************************************************************







