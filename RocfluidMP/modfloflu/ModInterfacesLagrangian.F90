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
! Purpose: set explicit interfaces to subroutines and functions.
!
! Description: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModInterfacesLagrangian.F90,v 1.30 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfacesLagrangian

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! Common routines
! =============================================================================

  SUBROUTINE PLAG_AllocateMemory( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_AllocateMemory

  SUBROUTINE PLAG_AllocateDataBuffers( regions, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE PLAG_AllocateDataBuffers 

  SUBROUTINE PLAG_BuildVersionString( versionString )
    CHARACTER(*) :: versionString
  END SUBROUTINE PLAG_BuildVersionString

  SUBROUTINE PLAG_CECellsWrapper( regions )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region), POINTER  :: regions(:)
  END SUBROUTINE PLAG_CECellsWrapper
 
  SUBROUTINE PLAG_InitPatchData(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_InitPatchData
 
  SUBROUTINE PLAG_InitSolution( iReg,region )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iReg
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_InitSolution

  SUBROUTINE PLAG_NonCvUpdate( region )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region), POINTER :: region
  END SUBROUTINE PLAG_NonCvUpdate
  
  SUBROUTINE PLAG_PrintUserInput( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_PrintUserInput

  SUBROUTINE PLAG_ReadSolution( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_ReadSolution

  SUBROUTINE PLAG_StatMapping( global )
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE PLAG_StatMapping

  SUBROUTINE PLAG_WriteSolution( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_WriteSolution

  SUBROUTINE PLAG_UserInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_UserInput

  SUBROUTINE PLAG_PatchUpdateWrapper( regions )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region), POINTER  :: regions(:)
  END SUBROUTINE PLAG_PatchUpdateWrapper
 
  SUBROUTINE PLAG_RkInit( region,iStage )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region), INTENT(INOUT) :: region
    INTEGER,        INTENT(IN)    :: iStage
  END SUBROUTINE PLAG_RkInit
  
  SUBROUTINE PLAG_RkUpdateWrapper( region, iReg, iStage )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region)      :: region
    INTEGER, INTENT(IN) :: iReg, iStage
  END SUBROUTINE PLAG_RkUpdateWrapper

# ifdef RFLO
! =============================================================================
! Rocflo-specific routines
! =============================================================================

  SUBROUTINE PLAG_RFLO_SetMetrics( regions )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region), POINTER  :: regions(:)
  END SUBROUTINE PLAG_RFLO_SetMetrics
#endif

# ifdef RFLU
! =============================================================================
! Rocflu-specific routines
! =============================================================================

  SUBROUTINE PLAG_INRT_AllocMemTStep(pRegion,pPlag)
    USE ModDataStruct, ONLY: t_region
    USE ModPartLag, ONLY: t_plag
    TYPE(t_region), POINTER :: pRegion 
    TYPE(t_plag), POINTER :: pPlag     
  END SUBROUTINE PLAG_INRT_AllocMemTStep

  SUBROUTINE PLAG_INRT_DeallocMemTStep(pRegion,pPlag)
    USE ModDataStruct, ONLY: t_region
    USE ModPartLag, ONLY: t_plag
    TYPE(t_region), POINTER :: pRegion 
    TYPE(t_plag), POINTER :: pPlag     
  END SUBROUTINE PLAG_INRT_DeallocMemTStep

  SUBROUTINE PLAG_RFLU_AllocMemSol(pRegion,pPlag)
    USE ModDataStruct, ONLY: t_region
    USE ModPartLag, ONLY: t_plag
    TYPE(t_region), POINTER :: pRegion 
    TYPE(t_plag), POINTER :: pPlag     
  END SUBROUTINE PLAG_RFLU_AllocMemSol

  SUBROUTINE PLAG_RFLU_AllocMemSolTile(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_AllocMemSolTile

  SUBROUTINE PLAG_RFLU_AllocMemTStep(pRegion,pPlag)
    USE ModDataStruct, ONLY: t_region
    USE ModPartLag, ONLY: t_plag
    TYPE(t_region), POINTER :: pRegion 
    TYPE(t_plag), POINTER :: pPlag    
  END SUBROUTINE PLAG_RFLU_AllocMemTStep

  SUBROUTINE PLAG_RFLU_AllocMemTStepTile(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE PLAG_RFLU_AllocMemTStepTile  
  
  SUBROUTINE PLAG_RFLU_CorrectMixtProperties(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_CorrectMixtProperties
  
  SUBROUTINE PLAG_RFLU_DeallocMemSol(pRegion,pPlag)
    USE ModDataStruct, ONLY: t_region
    USE ModPartLag, ONLY: t_plag
    TYPE(t_region), POINTER :: pRegion 
    TYPE(t_plag), POINTER :: pPlag     
  END SUBROUTINE PLAG_RFLU_DeallocMemSol

  SUBROUTINE PLAG_RFLU_DeallocMemSolTile(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_DeallocMemSolTile

  SUBROUTINE PLAG_RFLU_DeallocMemTStep(pRegion,pPlag)
    USE ModDataStruct, ONLY: t_region
    USE ModPartLag, ONLY: t_plag
    TYPE(t_region), POINTER :: pRegion 
    TYPE(t_plag), POINTER :: pPlag    
  END SUBROUTINE PLAG_RFLU_DeallocMemTStep
  
  SUBROUTINE PLAG_RFLU_DeallocMemTStepTile(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE PLAG_RFLU_DeallocMemTStepTile  
  
  SUBROUTINE PLAG_RFLU_InitSolutionScratch(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE PLAG_RFLU_InitSolutionScratch
  
  SUBROUTINE PLAG_RFLU_InitSolutionRandom(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE PLAG_RFLU_InitSolutionRandom
    
  SUBROUTINE PLAG_RFLU_InitSolFromSerial(pRegion,pRegionSerial)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion,pRegionSerial
  END SUBROUTINE PLAG_RFLU_InitSolFromSerial
  
  SUBROUTINE PLAG_RFLU_InitSolSerialWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE PLAG_RFLU_InitSolSerialWrapper  
    
  SUBROUTINE PLAG_RFLU_ReadSolutionASCII(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_ReadSolutionASCII
  
  SUBROUTINE PLAG_RFLU_ReadSolutionBinary(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_ReadSolutionBinary

  SUBROUTINE PLAG_RFLU_WriteSolutionASCII(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_WriteSolutionASCII

  SUBROUTINE PLAG_RFLU_WriteSolutionBinary(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_WriteSolutionBinary
#endif
     
  END INTERFACE

END MODULE ModInterfacesLagrangian

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfacesLagrangian.F90,v $
! Revision 1.30  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.29  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.28  2007/03/20 17:32:33  fnajjar
! Deleted obsolete entries with new module PLAG_ModDimensions
!
! Revision 1.27  2007/03/15 21:59:32  haselbac
! Renamed IF for PLAG_RFLU_InitSolSerial
!
! Revision 1.26  2006/05/05 17:19:34  haselbac
! Added if for PLAG_RFLU_InitSolSerial
!
! Revision 1.25  2005/11/30 22:15:03  fnajjar
! Added IF for PLAG_RFLU_CorrectMixtProperties
!
! Revision 1.24  2005/05/18 22:07:16  fnajjar
! Added interface for new init routine
!
! Revision 1.23  2004/12/29 23:26:58  wasistho
! prepared statistics for PLAG and PEUL
!
! Revision 1.22  2004/12/01 00:09:07  wasistho
! added BuildVersionString
!
! Revision 1.21  2004/11/14 19:44:30  haselbac
! Changed interface
!
! Revision 1.20  2004/11/04 16:43:07  fnajjar
! Added Interface call to PLAG_SetDimensions
!
! Revision 1.19  2004/10/10 20:04:35  fnajjar
! Included interface for solution generated by random state
!
! Revision 1.18  2004/08/23 23:08:09  fnajjar
! Added interface calls for binary IO
!
! Revision 1.17  2004/07/28 18:55:32  fnajjar
! Included interface for PLAG_SetMaxDimensions
!
! Revision 1.16  2004/07/26 19:03:07  fnajjar
! Included interface call to PLAG_INRT_DeallocMemTStep routine
!
! Revision 1.15  2004/07/26 17:05:50  fnajjar
! moved allocation of inrtSources into Rocpart
!
! Revision 1.14  2004/03/05 23:21:50  haselbac
! Added interface for PLAG_InitPatchData
!
! Revision 1.13  2004/02/26 21:01:58  haselbac
! Added RFLU routines, modified PLAG_PatchUpdateWrapper entry
!
! Revision 1.12  2004/02/06 21:27:24  fnajjar
! Included proper INTENT to Interfaces
!
! Revision 1.11  2003/11/12 21:18:17  fnajjar
! Added Corner-Edge cells calls
!
! Revision 1.10  2003/04/14 21:12:25  fnajjar
! Bug fix to include POINTER attribute appropriately for PLAG_patchUpdateWrapper
!
! Revision 1.9  2003/04/14 18:56:07  fnajjar
! Included POINTER attribute to regions for PLAG_PatchUpdateWrapper Interface
!
! Revision 1.8  2003/04/14 18:13:28  fnajjar
! Removed iReg from Interface sequence for PLAG_PatchUpdateWrapper
!
! Revision 1.7  2003/03/28 19:39:31  fnajjar
! Aligned with routines pertinent to RocfluidMP
!
! Revision 1.6  2003/03/04 22:12:34  jferry
! Initial import of Rocinteract
!
! Revision 1.5  2003/02/04 19:32:57  f-najjar
! Added Interfaces to PLAG_InjcTileUpdate PLAG_NonCvUpdate
!
! Revision 1.4  2003/01/23 17:06:27  f-najjar
! Included Interface call to PLAG_PatchBufferSendRecv
!
! Revision 1.3  2003/01/13 18:59:19  f-najjar
! Added PLAG_allocateDataBuffers
!
! Revision 1.2  2003/01/10 19:07:55  f-najjar
! Added call to PLAG_PatchUpdate
!
! Revision 1.1  2002/12/27 22:07:14  jblazek
! Splitted up RFLO_ModInterfaces and ModInterfaces.
!
!******************************************************************************






