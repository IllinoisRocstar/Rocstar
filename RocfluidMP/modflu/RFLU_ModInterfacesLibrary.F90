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
! Purpose: Set explicit interfaces to subroutines and functions in Rocflu
!   library.
!
! Description: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModInterfacesLibrary.F90,v 1.41 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************
  
MODULE RFLU_ModInterfacesLibrary

  IMPLICIT NONE

  INTERFACE

  SUBROUTINE dgesdd(JOBZ,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,IWORK,INFO)
    CHARACTER :: JOBZ
    INTEGER :: INFO,LDA,LDU,LDVT,LWORK,M,N
    INTEGER :: IWORK(:)
    DOUBLE PRECISION :: A(:,:),S(:),U(:,:),VT(:,:),WORK(:)  
  END SUBROUTINE dgesdd

  SUBROUTINE RFLU_AllocateMemoryTbc(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_AllocateMemoryTbc

  SUBROUTINE RFLU_BuildDataStruct(global,levels)
    USE ModGlobal, ONLY: t_global
    USE ModDataStruct, ONLY: t_level
    TYPE(t_global), POINTER :: global
    TYPE(t_level), POINTER :: levels(:)
  END SUBROUTINE RFLU_BuildDataStruct  

  SUBROUTINE RFLU_BuildInterpStencilC2V(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_BuildInterpStencilC2V

  SUBROUTINE RFLU_CheckDerivedUserInput(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:), POINTER :: regions  
  END SUBROUTINE RFLU_CheckDerivedUserInput

  SUBROUTINE RFLU_CheckPositivity(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_CheckPositivity 
 
  SUBROUTINE RFLU_CheckPositivity_GL(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_CheckPositivity_GL
 
  SUBROUTINE RFLU_CheckPositivityWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_CheckPositivityWrapper     
 
  SUBROUTINE RFLU_CheckUserInput(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:), POINTER :: regions  
  END SUBROUTINE RFLU_CheckUserInput
 
  SUBROUTINE RFLU_CheckValidity(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_CheckValidity   
 
  SUBROUTINE RFLU_CheckValidity_GL(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_CheckValidity_GL

  SUBROUTINE RFLU_CheckValidityWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_CheckValidityWrapper

  SUBROUTINE RFLU_ClosePostInfo(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_ClosePostInfo
 
  SUBROUTINE RFLU_CloseRestartInfo(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_CloseRestartInfo 
 
  SUBROUTINE RFLU_CompInterpWeightsC2V(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_CompInterpWeightsC2V

  SUBROUTINE RFLU_ComputeDCUHREInfo(global,NDIM,NF,KEY,MAXCLS,NW)
    USE ModGlobal, ONLY: t_global
    INTEGER, INTENT(IN) :: KEY,NDIM,NF
    INTEGER, INTENT(OUT) :: NW
    INTEGER, INTENT(INOUT) :: MAXCLS
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_ComputeDCUHREInfo

  SUBROUTINE RFLU_ConvertCvCons2Prim(pRegion,cvStateFuture)
    USE ModDataStruct, ONLY: t_region  
    INTEGER, INTENT(IN) :: cvStateFuture
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_ConvertCvCons2Prim

  SUBROUTINE RFLU_ConvertCvPrim2Cons(pRegion,cvStateFuture)
    USE ModDataStruct, ONLY: t_region  
    INTEGER, INTENT(IN) :: cvStateFuture
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_ConvertCvPrim2Cons

  SUBROUTINE RFLU_CreateGrid(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_CreateGrid

  SUBROUTINE RFLU_CreateRegions(global,iLev,levels)
    USE ModGlobal, ONLY: t_global
    USE ModDataStruct, ONLY: t_level
    INTEGER, INTENT(IN) :: iLev
    TYPE(t_global), POINTER :: global
    TYPE(t_level), DIMENSION(:), POINTER :: levels
  END SUBROUTINE RFLU_CreateRegions

  LOGICAL FUNCTION RFLU_DecideNeedBGradFace(pRegion,pPatch)
!  LOGICAL FUNCTION RFLU_DecideNeedBGradFace(pRegion,iPatch)
    USE ModDataStruct, ONLY: t_region
    USE ModBndPatch, ONLY: t_patch
    TYPE(t_region), POINTER :: pRegion
    INTEGER :: iPatch
    TYPE(t_patch), POINTER :: pPatch
  END FUNCTION RFLU_DecideNeedBGradFace

  LOGICAL FUNCTION RFLU_DecidePrint(global)
    USE ModGlobal, ONLY: t_global    
    TYPE(t_global), POINTER :: global    
  END FUNCTION RFLU_DecidePrint
    
  LOGICAL FUNCTION RFLU_DecideWrite(global)
    USE ModGlobal, ONLY: t_global    
    TYPE(t_global), POINTER :: global    
  END FUNCTION RFLU_DecideWrite
  
  SUBROUTINE RFLU_DerivedInputValues(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:), POINTER :: regions  
  END SUBROUTINE RFLU_DerivedInputValues

  SUBROUTINE RFLU_DestroyGrid(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_DestroyGrid

  SUBROUTINE RFLU_EnforceBoundsWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_EnforceBoundsWrapper

  SUBROUTINE RFLU_GetUserInput(regions,inPrep)
    USE ModDataStruct, ONLY: t_region
    LOGICAL, OPTIONAL :: inPrep
    TYPE(t_region), DIMENSION(:), POINTER :: regions
  END SUBROUTINE RFLU_GetUserInput

  INTEGER FUNCTION RFLU_GetCvLoc(global,fluidModel,var)
    USE ModGlobal, ONLY: t_global
    INTEGER, INTENT(IN) :: fluidModel,var
    TYPE(t_global), POINTER :: global
  END FUNCTION RFLU_GetCvLoc
  
  SUBROUTINE RFLU_InitGlobal(casename,verbLevel,communicator,global)
    USE ModDataTypes
    USE ModGlobal, ONLY: t_global
    CHARACTER(*), INTENT(IN) :: casename
    INTEGER, INTENT(IN) :: communicator,verbLevel    
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_InitGlobal

  SUBROUTINE RFLU_InitInputValues(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:), POINTER :: regions  
  END SUBROUTINE RFLU_InitInputValues

  SUBROUTINE RFLU_InvertMatrixSVD(global,nRows,nCols,a,aInv,sCount)
    USE ModDataTypes  
    USE ModGlobal, ONLY: t_global
    INTEGER, INTENT(IN) :: nRows,nCols
    INTEGER, INTENT(OUT) :: sCount
    REAL(RFREAL) :: a(nRows,nCols),aInv(nCols,nRows)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_InvertMatrixSVD

  SUBROUTINE RFLU_OpenPostInfo(global,fileStatus,fileExists)
    USE ModGlobal, ONLY: t_global
    LOGICAL, INTENT(OUT) :: fileExists
    INTEGER, INTENT(IN) :: fileStatus
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_OpenPostInfo

  SUBROUTINE RFLU_OpenRestartInfo(global,filePosition,fileExists)
    USE ModGlobal, ONLY: t_global
    LOGICAL, INTENT(OUT) :: fileExists
    INTEGER, INTENT(IN) :: filePosition
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_OpenRestartInfo 

  SUBROUTINE RFLU_PrintChangeInfo(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE RFLU_PrintChangeInfo

  SUBROUTINE RFLU_PrintFlowInfo(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE RFLU_PrintFlowInfo

  SUBROUTINE RFLU_PrintFlowInfoWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE RFLU_PrintFlowInfoWrapper

  SUBROUTINE RFLU_PrintGridInfo(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE RFLU_PrintGridInfo

  SUBROUTINE RFLU_PrintLocInfo(pRegion,locUnsorted,nLocUnsorted,locInfoMode, & 
                               outputMode)
    USE ModParameters
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: locInfoMode,nLocUnsorted,outputMode
    INTEGER, INTENT(INOUT) :: locUnsorted(nLocUnsorted,MIN_VAL:MAX_VAL)
    TYPE(t_region), POINTER :: pRegion 
  END SUBROUTINE RFLU_PrintLocInfo

  SUBROUTINE RFLU_PrintUserInput(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLU_PrintUserInput

  SUBROUTINE RFLU_PrintWarnInfo(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global    
  END SUBROUTINE RFLU_PrintWarnInfo

  SUBROUTINE RFLU_RandomInit( regions )
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLU_RandomInit

  SUBROUTINE RFLU_ReadRestartInfo(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global  
  END SUBROUTINE RFLU_ReadRestartInfo

  SUBROUTINE RFLU_ReadTbcInputFile(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_ReadTbcInputFile

  SUBROUTINE RFLU_ReadTbcSection(pRegion,tbcType)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
    INTEGER, INTENT(IN)     :: tbcType
  END SUBROUTINE RFLU_ReadTbcSection

  SUBROUTINE RFLU_ScalarFirst(pRegion,nVarScal,cvScal,resScal)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: nVarScal
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cvScal
    REAL(RFREAL), DIMENSION(:,:), INTENT(INOUT) :: resScal
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_ScalarFirst
  
  SUBROUTINE RFLU_ScalarSecond(pRegion,nVarScal,cvScal,gradCellScal,resScal)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: nVarScal
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cvScal
    REAL(RFREAL), DIMENSION(:,:), INTENT(INOUT) :: resScal
    REAL(RFREAL), DIMENSION(:,:,:), INTENT(IN) :: gradCellScal    
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_ScalarSecond     
  
  SUBROUTINE RFLU_ScalarFirstPatch(pRegion,pPatch,nVarScal,cvScal,valScal, & 
                                   resScal)
    USE ModDataTypes
    USE ModBndPatch, ONLY: t_bcvalues,t_patch
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: nVarScal
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cvScal
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: resScal
    TYPE(t_bcvalues) :: valScal
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_ScalarFirstPatch 
  
  SUBROUTINE RFLU_ScalarSecondPatch(pRegion,pPatch,nVarScal,cvScal, & 
                                    gradCellScal,valScal,resScal)
    USE ModDataTypes
    USE ModBndPatch, ONLY: t_bcvalues,t_patch
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: nVarScal
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cvScal
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: resScal
    REAL(RFREAL), DIMENSION(:,:,:), INTENT(IN) :: gradCellScal    
    TYPE(t_bcvalues) :: valScal
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_ScalarSecondPatch   
  
  SUBROUTINE RFLU_ScalarCheckPositivity(pRegion,moduleType,nVarScal,cvScal)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: nVarScal,moduleType
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cvScal    
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE RFLU_ScalarCheckPositivity  
  
  SUBROUTINE RFLU_ScalarInitRhs(pRegion,nVarScal,dissScal,resScal)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: nVarScal
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: dissScal
    REAL(RFREAL), DIMENSION(:,:), INTENT(INOUT) :: resScal
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_ScalarInitRhs
  
  SUBROUTINE RFLU_ScalarViscousFluxes(pRegion,nVarScal,tvScal,gradScal,resScal)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: nVarScal
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: tvScal
    REAL(RFREAL), DIMENSION(:,:), INTENT(INOUT) :: resScal 
    REAL(RFREAL), DIMENSION(:,:,:), INTENT(IN) :: gradScal   
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_ScalarViscousFluxes 
  
  SUBROUTINE RFLU_SetDependentVars(pRegion,icgBeg,icgEnd)
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: icgBeg,icgEnd
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_SetDependentVars
    
  SUBROUTINE RFLU_SetDerivedUserInput(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:), POINTER :: regions  
  END SUBROUTINE RFLU_SetDerivedUserInput 
  
  SUBROUTINE RFLU_SetGasVars(pRegion,icgBeg,icgEnd)
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: icgBeg,icgEnd
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_SetGasVars
  
  SUBROUTINE RFLU_SetModuleType(global,moduleType)
    USE ModGlobal, ONLY: t_global    
    INTEGER, INTENT(IN) :: moduleType
    TYPE(t_global), POINTER :: global      
  END SUBROUTINE RFLU_SetModuleType
  
  SUBROUTINE RFLU_SetRestartTimeFlag(global)
    USE ModGlobal, ONLY: t_global    
    TYPE(t_global), POINTER :: global    
  END SUBROUTINE RFLU_SetRestartTimeFlag

  SUBROUTINE RFLU_SetTransportVars(pRegion,icgBeg,icgEnd)
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: icgBeg,icgEnd
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_SetTransportVars  

  SUBROUTINE RFLU_SetVarInfo(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE RFLU_SetVarInfo 
  
  SUBROUTINE RFLU_SetVarInfoWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE RFLU_SetVarInfoWrapper  

  SUBROUTINE RFLU_SetVars(pRegion,icgBeg,icgEnd)
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: icgBeg,icgEnd
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_SetVars
  
  SUBROUTINE RFLU_SetVarsContWrapper(pRegion,icgBeg,icgEnd)
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: icgBeg,icgEnd
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_SetVarsContWrapper
  
  SUBROUTINE RFLU_SetVarsDiscWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_SetVarsDiscWrapper    
  
  SUBROUTINE RFLU_SetVarsWrapper(pRegion,icgBeg,icgEnd)
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: icgBeg,icgEnd
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_SetVarsWrapper

  LOGICAL FUNCTION RFLU_TestIsFirstRegion(pRegion) 
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion                                       
  END FUNCTION RFLU_TestIsFirstRegion  

  SUBROUTINE RFLU_UpdateCommLists(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_UpdateCommLists

  SUBROUTINE RFLU_UserInput(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:), POINTER :: regions
  END SUBROUTINE RFLU_UserInput
  
  SUBROUTINE RFLU_WriteRestartInfo(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_WriteRestartInfo 

#ifdef STATS
  SUBROUTINE RFLU_ReadStat( region )
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region    
  END SUBROUTINE RFLU_ReadStat

  SUBROUTINE RFLU_WriteStat( region )
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region    
  END SUBROUTINE RFLU_WriteStat
#endif

  SUBROUTINE RFLU_WriteVersionString(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global    
  END SUBROUTINE RFLU_WriteVersionString

  SUBROUTINE RFLU_ZeroVirtualCellVars(pRegion,var)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region  
    REAL(RFREAL), POINTER :: var(:,:)
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_ZeroVirtualCellVars

  END INTERFACE

END MODULE RFLU_ModInterfacesLibrary

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModInterfacesLibrary.F90,v $
! Revision 1.41  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.40  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.39  2006/08/19 15:39:08  mparmar
! Added RFLU_DecideNeedBGradFace
!
! Revision 1.38  2006/08/04 03:05:58  haselbac
! Removed interface for RFLU_TestIsBoundaryCell
!
! Revision 1.37  2006/03/26 20:22:04  haselbac
! Added ifs for GL routines
!
! Revision 1.36  2006/02/06 23:56:55  haselbac
! Adapted to changes in arg list of RFLU_InitGlobal
!
! Revision 1.35  2005/12/24 21:30:10  haselbac
! Removed ifs for routines moved into modules
!
! Revision 1.34  2005/11/10 22:23:34  fnajjar
! ACH: Added ifs for new SetVars routines
!
! Revision 1.33  2005/06/09 20:21:40  haselbac
! Removed interface for RFLU_CheckMoveGridInput
!
! Revision 1.32  2005/05/16 20:42:43  haselbac
! Removed/added interfaces
!
! Revision 1.31  2005/04/29 12:55:30  haselbac
! Removed interface for RFLU_DecideWriteProbe
!
! Revision 1.30  2005/04/15 15:06:59  haselbac
! Removed intfcs for RFLU_InitFlowScratch, changed intfc for RFLU_SetVarsXyz
!
! Revision 1.29  2005/03/11 02:20:38  haselbac
! Added and removed interfaces
!
! Revision 1.28  2004/11/14 19:45:49  haselbac
! Added interfaces, some clean-up
!
! Revision 1.27  2004/11/06 03:18:53  haselbac
! Added interface for RFLU_GetCvLoc
!
! Revision 1.26  2004/11/05 21:48:22  fnajjar
! Added interface for RFLU_TestInCellFancy
!
! Revision 1.25  2004/11/05 20:25:43  haselbac
! Removed interface for RFLU_ComputeLineCellXSectDist
!
! Revision 1.24  2004/11/02 02:32:29  haselbac
! Added interfaces for RFLU_SetVarInfo and RFLU_SetVarInfoWrapper
!
! Revision 1.23  2004/10/19 19:28:07  haselbac
! Added and removed interfaces
!
! Revision 1.22  2004/07/30 22:47:36  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.21  2004/07/06 15:14:43  haselbac
! Removed many subroutines and functions bcos moved into modules
!
! Revision 1.20  2004/06/25 20:09:33  haselbac
! Added entry for RFLU_SetRestartTimeFlag
!
! Revision 1.19  2004/05/05 20:54:39  fnajjar
! Adapted entries for RFLU_ComputeLineCellXSect{Dist}
!
! Revision 1.18  2004/03/19 21:19:54  haselbac
! Added interface for RFLU_InitFlowScratch
!
! Revision 1.17  2004/03/17 04:26:53  haselbac
! Modified interface for RFLU_WriteDimensionsWrapper
!
! Revision 1.16  2004/03/01 23:54:53  haselbac
! Added interfaces for RFLU_ComputeLineCellXSect* routines
!
! Revision 1.15  2004/02/26 21:02:03  haselbac
! Added interfaces for read/write dimensions wrapper routines
!
! Revision 1.14  2004/02/23 23:04:05  haselbac
! Added interface for RFLU_ComputeExactFlowProudman
!
! Revision 1.13  2004/02/13 02:57:42  haselbac
! Added interface for RFLU_TestIsBoundaryCell
!
! Revision 1.12  2004/01/29 22:57:36  haselbac
! Added and deleted various interfaces
!
! Revision 1.11  2003/12/04 03:28:49  haselbac
! Added and removed several interfaces
!
! Revision 1.10  2003/11/25 21:03:25  haselbac
! Added interfaces for new routines
!
! Revision 1.9  2003/10/29 21:39:19  haselbac
! Added interfaces for RFLU_DecideXXX functions
!
! Revision 1.8  2003/08/07 15:32:26  haselbac
! Added RFLU_PrintWarnInfo
!
! Revision 1.7  2003/07/22 02:04:16  haselbac
! Removed RFLU_CreateInterpolant (in module)
!
! Revision 1.6  2003/06/20 22:34:33  haselbac
! Added interfaces for restart info file routines
!
! Revision 1.5  2003/06/09 14:03:46  haselbac
! Added RFLU_CreateRegions
!
! Revision 1.4  2003/06/04 22:09:28  haselbac
! Added, deleted, and modified some interfaces
!
! Revision 1.3  2003/06/04 20:05:54  jferry
! re-worked implementation of TBCs in unstructured code
!
! Revision 1.2  2003/04/28 22:45:15  haselbac
! Changed interface for RFLU_PrintLocInfo
!
! Revision 1.1  2003/04/10 14:37:10  haselbac
! Initial revision
!
! ******************************************************************************






