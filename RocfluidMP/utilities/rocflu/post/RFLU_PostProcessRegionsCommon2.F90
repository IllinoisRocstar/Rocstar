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
! Purpose: Wrapper for postprocessing results without merging the regions.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_PostProcessRegionsCommon2.F90,v 1.8 2008/12/06 08:45:06 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PostProcessRegionsCommon2(pRegion)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global 
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModParameters

  USE RFLU_ModBoundLists
  USE RFLU_ModCellMapping
  USE RFLU_ModDeallocateMemory
  USE RFLU_ModDimensions 
  USE RFLU_ModFaceList
  USE RFLU_ModGeometry
  USE RFLU_ModPatchCoeffs
  USE RFLU_ModPlottingVars  
  USE RFLU_ModStencilsCells
  USE RFLU_ModStencilsVert
  USE RFLU_ModVertexLists

#ifdef PLAG
  USE PLAG_ModSurfStats
#endif

  USE ModInterfaces, ONLY: RFLU_DeallocMemSolWrapper, &
                           RFLU_DeallocMemVertWrapper, &
                           RFLU_DecideBuildGeometry, & 
                           RFLU_DestroyGrid, &
                           RFLU_PrintFlowInfoWrapper, &
                           RFLU_SetVarInfoWrapper, &
                           RFLU_SetVarsWrapper
    
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
  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PostProcessRegionsCommon2.F90,v $ $Revision: 1.8 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_PostProcessRegionsCommon2',&
  'RFLU_PostProcessRegionsCommon2.F90')

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  IF ( global%postPlotType == PLOT_GRID_FLOW ) THEN
    IF ( RFLU_DecideComputePlottingVars(pRegion) .EQV. .TRUE. ) THEN
      CALL RFLU_DestroyPlottingVarMaps(pRegion)
      CALL RFLU_DestroyPlottingVars(pRegion)
    END IF ! RFLU_DecideComputePlottingVars
  END IF ! global%postPlotType

  IF ( RFLU_DecideBuildGeometry(pRegion%global) .EQV. .TRUE. ) THEN
    CALL RFLU_DestroyGeometry(pRegion)                            
  END IF ! RFLU_DecideBuildGeometry 

  CALL RFLU_DestroyFaceList(pRegion) 
  CALL RFLU_DestroyBVertexLists(pRegion)                   
  CALL RFLU_DestroyCellMapping(pRegion)             
  CALL RFLU_DeallocMemSolWrapper(pRegion)
  CALL RFLU_DeallocateMemoryGSpeeds(pRegion)

  IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN 
    CALL RFLU_DestroyPatchCoeffs(pRegion)
  END IF ! global%patchCoeffFlag     

#ifdef PLAG
  IF ( global%plagUsed .EQV. .TRUE. ) THEN 
    CALL PLAG_DestroySurfStats(pRegion)
  END IF ! global%plagUsed     
#endif

  CALL RFLU_DestroyGrid(pRegion)

  IF ( (global%postPlotType == PLOT_GRID_FLOW) .AND. & 
       (global%postInterpType /= INTERP_TYPE_NONE) ) THEN 
    CALL RFLU_DeallocMemVertWrapper(pRegion)
  END IF ! global%postPlotType      
                        
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PostProcessRegionsCommon2

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PostProcessRegionsCommon2.F90,v $
! Revision 1.8  2008/12/06 08:45:06  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:18:17  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2007/03/19 21:45:19  haselbac
! Adapted to changes related to plotting variables
!
! Revision 1.5  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.4  2005/12/10 16:57:46  haselbac
! Added use of RFLU_DecideComputePlottingVars
!
! Revision 1.3  2005/11/10 16:51:30  fnajjar
! Added plagUsed IF statement around PLAG routines
!
! Revision 1.2  2005/10/05 20:53:00  haselbac
! Fixed missing module bugs
!
! Revision 1.1  2005/10/05 20:23:34  haselbac
! Initial revision
!
! ******************************************************************************







