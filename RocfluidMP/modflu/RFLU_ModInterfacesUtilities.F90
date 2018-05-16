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
! Purpose: Set explicit interfaces to subroutines and functions for Rocflu
!   utilities.
!
! Description: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModInterfacesUtilities.F90,v 1.25 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************
  
MODULE RFLU_ModInterfacesUtilities

  IMPLICIT NONE

  INTERFACE

! ==============================================================================
! Rocflu partitioner
! ==============================================================================

  SUBROUTINE RFLU_ReadConvGridWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_ReadConvGridWrapper

  SUBROUTINE RFLU_ReadFormatsSection
  END SUBROUTINE RFLU_ReadFormatsSection

  SUBROUTINE RFLU_ReadInputFile(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:), POINTER :: regions       
  END SUBROUTINE RFLU_ReadInputFile

  SUBROUTINE RFLU_USER_EnforcePatchCoords(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_USER_EnforcePatchCoords

! ==============================================================================
! Rocflu post-processor
! ==============================================================================

  SUBROUTINE RFLU_AllocateMemoryVert(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_AllocateMemoryVert

  SUBROUTINE RFLU_AllocMemSolWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_AllocMemSolWrapper
  
  SUBROUTINE RFLU_AllocMemVertWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_AllocMemVertWrapper

  SUBROUTINE RFLU_ComputeExactFlowError(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_ComputeExactFlowError
  
  SUBROUTINE RFLU_ComputeExactFlowProbeError(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_ComputeExactFlowProbeError  
  
  SUBROUTINE RFLU_ComputeVertexVariables(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_ComputeVertexVariables
  
  SUBROUTINE RFLU_DeallocateMemoryVert(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_DeallocateMemoryVert

  SUBROUTINE RFLU_DeallocMemSolWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_DeallocMemSolWrapper
  
  SUBROUTINE RFLU_DeallocMemVertWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_DeallocMemVertWrapper

  LOGICAL FUNCTION RFLU_DecideBuildGeometry(global)
    USE ModGlobal, ONLY: t_global 
    TYPE(t_global), POINTER :: global  
  END FUNCTION RFLU_DecideBuildGeometry  

  LOGICAL FUNCTION RFLU_DecideBuildStencilsWeights(global)
    USE ModGlobal, ONLY: t_global 
    TYPE(t_global), POINTER :: global  
  END FUNCTION RFLU_DecideBuildStencilsWeights

  SUBROUTINE RFLU_InterpolateWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InterpolateWrapper
  
  SUBROUTINE RFLU_MergePostProcessRegions(levels)
    USE ModDataStruct, ONLY: t_level
    TYPE(t_level), POINTER :: levels(:)
  END SUBROUTINE RFLU_MergePostProcessRegions    
  
  SUBROUTINE RFLU_PostProcessRegions(levels)
    USE ModDataStruct, ONLY: t_level
    TYPE(t_level), POINTER :: levels(:)
  END SUBROUTINE RFLU_PostProcessRegions 
  
  SUBROUTINE RFLU_PostProcessRegionsCommon1(pRegion,postInfoFileExists)
    USE ModDataStruct, ONLY: t_region 
    LOGICAL, INTENT(IN) :: postInfoFileExists 
    TYPE(t_region), POINTER :: pRegion   
  END SUBROUTINE RFLU_PostProcessRegionsCommon1 
  
  SUBROUTINE RFLU_PostProcessRegionsCommon2(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion   
  END SUBROUTINE RFLU_PostProcessRegionsCommon2
     
  SUBROUTINE RFLU_PostProcessRegions_ENS(levels)
    USE ModDataStruct, ONLY: t_level
    TYPE(t_level), POINTER :: levels(:)
  END SUBROUTINE RFLU_PostProcessRegions_ENS       
  
  SUBROUTINE RFLU_ReadPostInfo(pRegion,readMode)
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: readMode
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_ReadPostInfo

  SUBROUTINE RFLU_SetPatchPlotFlags(pRegion)
    USE ModDataStruct, ONLY: t_region  
    TYPE(t_region), POINTER :: pRegion        
  END SUBROUTINE RFLU_SetPatchPlotFlags

! ==============================================================================
! Rocflu initializor
! ==============================================================================

  SUBROUTINE RFLU_InitBcDataHardCode(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InitBcDataHardCode
  
  SUBROUTINE RFLU_InitFlowHardCode(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InitFlowHardCode

  SUBROUTINE RFLU_InitFlowHardCodeLim(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InitFlowHardCodeLim
  
  SUBROUTINE RFLU_InitFlowHardCodeLimWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InitFlowHardCodeLimWrapper 

  SUBROUTINE RFLU_InitFlowHardCodeWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InitFlowHardCodeWrapper 

  SUBROUTINE RFLU_InitFlowScratch(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InitFlowScratch

  SUBROUTINE RFLU_InitFlowScratchWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_InitFlowScratchWrapper 

  SUBROUTINE RFLU_InitFlowSerialWrapper(pRegion,pRegionSerial)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion,pRegionSerial  
  END SUBROUTINE RFLU_InitFlowSerialWrapper 

! ==============================================================================
! Rocflu picker
! ==============================================================================

  SUBROUTINE RFLU_PickRegionsCoord(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE RFLU_PickRegionsCoord

  SUBROUTINE RFLU_PickRegionsManual(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:), POINTER :: regions   
  END SUBROUTINE RFLU_PickRegionsManual

  SUBROUTINE RFLU_PickSpecialCells(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_PickSpecialCells
  
  SUBROUTINE RFLU_PickSpecialFaces(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_PickSpecialFaces  

  SUBROUTINE RFLU_WritePostInfo(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_WritePostInfo

  END INTERFACE

END MODULE RFLU_ModInterfacesUtilities

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModInterfacesUtilities.F90,v $
! Revision 1.25  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.24  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.23  2006/01/06 22:22:11  haselbac
! Added if for RFLU_DecideBuildStencilsWeights
!
! Revision 1.22  2005/11/10 02:27:36  haselbac
! Modified interfaces for limited hard-coded init following name change
!
! Revision 1.21  2005/10/05 20:07:16  haselbac
! Added interfaces, removed some old unnecessary ones
!
! Revision 1.20  2005/08/09 00:58:43  haselbac
! Added entry for RFLU_SetPatchPlotFlags
!
! Revision 1.19  2005/04/29 12:56:19  haselbac
! Added interface for RFLU_ComputeExactFlowProbeError
!
! Revision 1.18  2005/04/15 15:07:00  haselbac
! Added section for rfluinit interfaces, changed section from rfluprep to rflupart
!
! Revision 1.17  2005/03/29 22:30:12  haselbac
! Added interface for RFLU_InitFlowHardCodeLimited
!
! Revision 1.16  2004/12/29 21:07:53  haselbac
! Added entries for new procedures
!
! Revision 1.15  2004/10/19 19:28:13  haselbac
! Added interface for RFLU_ConvGridWrapper
!
! Revision 1.14  2004/09/27 01:44:42  haselbac
! Clean-up and added interf for RFLU_PickSpecialFaces
!
! Revision 1.13  2004/07/21 14:59:28  haselbac
! Added interface for RFLU_DecideBuildGeometry
!
! Revision 1.12  2004/07/06 15:14:45  haselbac
! Removed many subroutines and functions bcos moved into modules
!
! Revision 1.11  2004/03/19 21:20:26  haselbac
! Removed RFLU_InitFlowScratch from prep, added RFLU_ReInitFlowWrapper for 
! rinit
!
! Revision 1.10  2004/02/26 21:02:05  haselbac
! Added interfaces for alloc/dealloc routines
!
! Revision 1.9  2004/01/29 22:57:39  haselbac
! Added ifs for RFLU_InitBcDataHardCode and RFLU_ComputeExactFlowError
!
! Revision 1.8  2003/11/25 21:03:27  haselbac
! Added interfaces for new routines
!
! Revision 1.7  2003/09/15 00:37:43  haselbac
! Added RFLU_InitFlowHardCode/Scratch
!
! Revision 1.6  2003/08/19 22:49:32  haselbac
! Added interfaces for COBALT conversion routines
!
! Revision 1.5  2003/08/07 15:33:03  haselbac
! Added RFLU_PickRegionsCoord and RFLU_PickRegionsManual
!
! Revision 1.4  2003/06/04 22:11:41  haselbac
! Added, deleted, and modified some interfaces
!
! Revision 1.3  2003/05/05 18:40:59  haselbac
! Added new merge routines
!
! Revision 1.2  2003/04/10 18:48:06  haselbac
! Added and deleted interf for reading CENTAUR grids
!
! Revision 1.1  2003/04/10 14:37:10  haselbac
! Initial revision
!
! ******************************************************************************






