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
! Purpose: Initialize patch data for Rocpart.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes:
!   1. Only treat injection boundaries right now.  
!   2. Avoid error if randUnif is 0, and set timefactor such that 
!      EXP(-50) = 1.9E-22
!
!******************************************************************************
!
! $Id: PLAG_InitPatchData.F90,v 1.7 2009/03/02 00:19:36 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_InitPatchData(pRegion)

  USE ModDataTypes
  USE ModPartLag, ONLY: t_tile_plag
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModRandom, ONLY: Rand1Uniform
  USE ModError
  USE ModParameters
  USE ModMPI

  USE PLAG_ModParameters
  
  USE PLAG_ModInterfaces, ONLY: PLAG_injcMakeParticle

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion  
   
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: injcDiamDist,iPatch,iTile,nPatches,nTiles
#ifdef RFLO
  INTEGER :: iLev,n1,n2
#endif
  REAL(RFREAL) :: randUnif
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_tile_plag), POINTER :: pTilePlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InitPatchData.F90,v $ $Revision: 1.7 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_InitPatchData',&
  'PLAG_InitPatchData.F90')

#ifdef RFLO
  iLev     = pRegion%currLevel
  nPatches = pRegion%nPatches
#endif
#ifdef RFLU
  nPatches = pRegion%grid%nPatches
#endif

  injcDiamDist = pRegion%plagInput%injcDiamDist

! ******************************************************************************
! Loop over patches
! ******************************************************************************

  DO iPatch=1,nPatches
#ifdef RFLO
    pPatch => pRegion%levels(iLev)%patches(iPatch)
#endif
#ifdef RFLU
    pPatch => pRegion%patches(iPatch)
#endif

! ==============================================================================
!   Select patch type
! ==============================================================================

    SELECT CASE ( pPatch%bcType )
    
! ------------------------------------------------------------------------------    
!     Injection boundary
! ------------------------------------------------------------------------------    
    
#ifdef RFLU
      CASE ( BC_INJECTION:BC_INJECTION+BC_RANGE, &
             BC_INFLOW:BC_INFLOW+BC_RANGE        ) 
#else
      CASE ( BC_INJECTION:BC_INJECTION+BC_RANGE )
#endif
! ----- Get dimensions ---------------------------------------------------------

#ifdef RFLO
        n1     = ABS(pPatch%l1end-pPatch%l1beg) + 1    
        n2     = ABS(pPatch%l2end-pPatch%l2beg) + 1
        nTiles = n1*n2
#endif
#ifdef RFLU
        nTiles = pPatch%nBFaces
#endif

        pTilePlag => pPatch%tilePlag


! ----- Loop over tiles --------------------------------------------------------

        DO iTile = 1,nTiles
          CALL PLAG_injcMakeParticle(pRegion,injcDiamDist, &
                                     pTilePlag%dv(DV_TILE_DIAM,iTile), &
                                     pTilePlag%dv(DV_TILE_SPLOAD,iTile))

          randUnif = Rand1Uniform(pRegion%randData) 

          IF ( randUnif <= 0.0_RFREAL ) THEN          
            pTilePlag%dv(DV_TILE_COUNTDOWN,iTile) = 50.0_RFREAL 
          ELSE
            pTilePlag%dv(DV_TILE_COUNTDOWN,iTile) = -LOG(randUnif)
          END IF ! randUnif
        END DO ! iTile
    END SELECT ! pPatch%bcType
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_InitPatchData

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InitPatchData.F90,v $
! Revision 1.7  2009/03/02 00:19:36  mtcampbe
! Added some ifdefs around Rocflo to disable particle injection on INFLOW
! boundaries and added some checks around MPI tags utilizing a new global
! data item, global%mpiTagMax.
!
! Revision 1.6  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/09/18 20:29:07  fnajjar
! Activated tile infrastructure for inflow bc
!
! Revision 1.3  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.2  2004/06/16 22:58:13  fnajjar
! Renamed injcModel to injcDiamDist and dv(TIMEFCTR) to dv(COUNTDOWN) for CRE kernel
!
! Revision 1.1  2004/03/05 23:15:54  haselbac
! Initial revision
!
!******************************************************************************







