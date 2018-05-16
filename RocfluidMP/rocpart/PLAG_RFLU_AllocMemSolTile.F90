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
! Purpose: Allocate memory for Lagrangian particle tile solution.
!
! Description: None.
!
! Input: 
!   pRegion     Region pointer
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_RFLU_AllocMemSolTile.F90,v 1.6 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_AllocMemSolTile(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModParameters  
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModPartLag, ONLY: t_tile_plag 
  USE ModMPI
   
  USE PLAG_ModParameters  
   
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
  INTEGER :: errorFlag,iCont,iPatch,iTile,iVar,nCont,nCv,nDv,nTiles
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_tile_plag), POINTER :: pTilePlag  

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_AllocMemSolTile.F90,v $ $Revision: 1.6 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_AllocMemSolTile',&
  'PLAG_RFLU_AllocMemSolTile.F90')

! ******************************************************************************
! Set pointers and variables 
! ******************************************************************************

  pGrid => pRegion%grid

  nCont = pRegion%plagInput%nCont

! ******************************************************************************
! Allocate memory
! ****************************************************************************** 

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
    
    IF ( (pPatch%bcType >= BC_INJECTION .AND. pPatch%bcType <= BC_INJECTION + BC_RANGE) .OR. &
         (pPatch%bcType >= BC_INFLOW    .AND. pPatch%bcType <= BC_INFLOW    + BC_RANGE)      ) THEN 
      pTilePlag => pPatch%tilePlag
            
      pPatch%tilePlag%nCv = CV_TILE_LAST + nCont
      pPatch%tilePlag%nDv = DV_TILE_LAST
            
      nCv = pTilePlag%nCv
      nDv = pTilePlag%nDv
            
      nTiles = pPatch%nBFaces     
            
      ALLOCATE(pTilePlag%cv(nCv,nTiles),STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global, ERR_ALLOCATE,__LINE__,'pTilePlag%cv') 
      END IF ! global%error

      ALLOCATE(pTilePlag%dv(nDv,nTiles),STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global, ERR_ALLOCATE,__LINE__,'pTilePlag%dv') 
      END IF ! global%error

      ALLOCATE(pTilePlag%cvTileMass(nCont),STAT=errorFlag)
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop(global, ERR_ALLOCATE,__LINE__,'pTilePlag%cvTileMass') 
      END IF ! global%error

      DO iCont = 1,nCont
        pTilePlag%cvTileMass(iCont) = CV_TILE_LAST + iCont
      END DO ! iCont            
    END IF ! pPatch%bcType   
  END DO ! iPatch

! ******************************************************************************
! Initialize memory
! ****************************************************************************** 

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
    
    IF ( (pPatch%bcType >= BC_INJECTION .AND. pPatch%bcType <= BC_INJECTION + BC_RANGE) .OR. &
         (pPatch%bcType >= BC_INFLOW    .AND. pPatch%bcType <= BC_INFLOW    + BC_RANGE)      ) THEN 
      pTilePlag => pPatch%tilePlag
      
      nCv = pTilePlag%nCv
      nDv = pTilePlag%nDv
            
      nTiles = pPatch%nBFaces     
            
      DO iTile = 1,nTiles
        DO iVar = 1,nCv
          pTilePlag%cv(iVar,iTile) = 0.0_RFREAL
        END DO ! iVar

        DO iVar = 1,nDv
          pTilePlag%dv(iVar,iTile) = 0.0_RFREAL
        END DO ! iVar     
      END DO ! iTile
    END IF ! pPatch%bcType   
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_AllocMemSolTile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_AllocMemSolTile.F90,v $
! Revision 1.6  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/09/18 20:32:49  fnajjar
! Activated tile datastructure for inflow bc
!
! Revision 1.3  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.2  2004/03/08 23:02:28  fnajjar
! Added initialization section within DO-loop construct
!
! Revision 1.1  2004/02/26 21:00:38  haselbac
! Initial revision
!
!******************************************************************************







