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
! Description: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_ModInterfaces.F90,v 1.48 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2001-2004 by the University of Illinois
!
!******************************************************************************
  
MODULE PLAG_ModInterfaces

  USE PLAG_ModInterfacesCECells
  
  IMPLICIT NONE
  
  INTERFACE

! =============================================================================
! PLAG specific library
! =============================================================================

  SUBROUTINE PLAG_AllocateMemoryTile( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_AllocateMemoryTile 

  SUBROUTINE PLAG_AppendDataFromBuffers( region, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)      :: region
    INTEGER, INTENT(IN) :: iReg
  END SUBROUTINE PLAG_AppendDataFromBuffers

  SUBROUTINE PLAG_BoundaryConditionsSet( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE PLAG_BoundaryConditionsSet

  SUBROUTINE PLAG_BufferDataRecv( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE PLAG_BufferDataRecv

  SUBROUTINE PLAG_BufferDataSend( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE PLAG_BufferDataSend

  SUBROUTINE PLAG_BufferSizeRecv( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE PLAG_BufferSizeRecv

  SUBROUTINE PLAG_BufferSizeSend( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE PLAG_BufferSizeSend

  SUBROUTINE PLAG_CalcBreakup( region, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
    INTEGER,        INTENT(IN)    :: iReg
  END SUBROUTINE PLAG_CalcBreakup
       
  SUBROUTINE PLAG_CalcDerivedVariables( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE PLAG_CalcDerivedVariables
  
  SUBROUTINE PLAG_CalcFaceCentroids( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_CalcFaceCentroids
        
  SUBROUTINE PLAG_CalcRhsPosition( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE PLAG_CalcRhsPosition

  SUBROUTINE PLAG_CheckUserInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_CheckUserInput

  SUBROUTINE PLAG_ClearRequests( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE PLAG_ClearRequests

  SUBROUTINE PLAG_ClearSizeSendRequests( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE PLAG_ClearSizeSendRequests

  SUBROUTINE PLAG_ClearDataSendRequests( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE PLAG_ClearDataSendRequests
  
  SUBROUTINE PLAG_CopyFaceVectors( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_CopyFaceVectors

  SUBROUTINE PLAG_DeallocateDataBuffers( regions, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE PLAG_DeallocateDataBuffers 
  
  SUBROUTINE PLAG_DerivedInputValues( regions )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region), POINTER  :: regions(:)
  END SUBROUTINE PLAG_DerivedInputValues 
  
  SUBROUTINE PLAG_GetCellIndices( region, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER :: iReg
  END SUBROUTINE PLAG_GetCellIndices 
  
  SUBROUTINE PLAG_GetCellIndicesOutflow( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_GetCellIndicesOutflow 

  SUBROUTINE PLAG_InCellTest( region, posPlag, indexSearch, &
                              ijkNR,ijkNRI,ijkNRJ,ijkNRK,   &
                              indexNew,cellLocate  )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER,                    INTENT(IN)    :: ijkNR,ijkNRI,ijkNRJ,ijkNRK
    INTEGER,      DIMENSION(4), INTENT(IN)    :: indexSearch
    INTEGER,      DIMENSION(4), INTENT(OUT)   :: indexNew
    LOGICAL,                    INTENT(OUT)   :: cellLocate  
    REAL(RFREAL), DIMENSION(3), INTENT(IN)    :: posPlag
  END SUBROUTINE PLAG_InCellTest

  SUBROUTINE PLAG_InCellTestRobust( region, posPlag, indexSearch, &
                                    indexNew,cellLocate  )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER,      DIMENSION(4), INTENT(IN)  :: indexSearch
    INTEGER,      DIMENSION(4), INTENT(OUT) :: indexNew
    LOGICAL,                    INTENT(OUT) :: cellLocate
    REAL(RFREAL), DIMENSION(3), INTENT(IN)  :: posPlag
  END SUBROUTINE PLAG_InCellTestRobust

  SUBROUTINE PLAG_InitInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_InitInputValues
  
  SUBROUTINE PLAG_InitMemory( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_InitMemory 

  SUBROUTINE PLAG_InjcEjectParticle( region, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)      :: region
    INTEGER, INTENT(IN) :: iReg
  END SUBROUTINE PLAG_InjcEjectParticle

  SUBROUTINE PLAG_InjcMakeParticle( region, injcDiamDist, diam, spLoad )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)            :: region
    INTEGER, INTENT(IN)       :: injcDiamDist
    REAL(RFREAL), INTENT(OUT) :: diam, spLoad
  END SUBROUTINE PLAG_InjcMakeParticle 

  SUBROUTINE PLAG_InjcSetInjection( region, pTilePlag, iTile, tCoeff, tSum, &
                                    poolVol, injectQ,  ratio )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE ModPartLag, ONLY    : t_tile_plag
    TYPE(t_region)              :: region
    TYPE(t_tile_plag), POINTER  :: pTilePlag
    INTEGER,      INTENT(IN)    :: iTile
    REAL(RFREAL), INTENT(IN)    :: tCoeff
    REAL(RFREAL), INTENT(OUT)   :: ratio
    REAL(RFREAL), INTENT(INOUT) :: poolVol,tSum
    LOGICAL, INTENT(OUT)        :: injectQ    
  END SUBROUTINE PLAG_InjcSetInjection 

  SUBROUTINE PLAG_InjcTileCalcRhs( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_InjcTileCalcRhs
     
  SUBROUTINE PLAG_InjcTileInitialize( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_InjcTileInitialize
  
  SUBROUTINE PLAG_InjcTileRKUpdate( region, iStage )
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region) :: region
    INTEGER        :: iStage
  END SUBROUTINE PLAG_InjcTileRKUpdate 
  
  SUBROUTINE PLAG_InjcTileUpdate( region, iReg, iStage )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region) :: region
    INTEGER        :: iReg, iStage
  END SUBROUTINE PLAG_InjcTileUpdate
  
  SUBROUTINE PLAG_InjcTileZeroRhs( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_InjcTileZeroRhs
    
  SUBROUTINE PLAG_IntrpMixtProperties( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE PLAG_IntrpMixtProperties
  
  SUBROUTINE PLAG_PatchBufferSendRecv( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_PatchBufferSendRecv
   
  SUBROUTINE PLAG_PatchExchangeConf( region,regionDes,patch,patchDes,iReg,iRegDes )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionDes
    TYPE(t_patch)  :: patch, patchDes
    INTEGER        :: iReg, iRegDes
  END SUBROUTINE PLAG_PatchExchangeConf

  SUBROUTINE PLAG_PatchGetBufferSize( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_PatchGetBufferSize 

  SUBROUTINE PLAG_PatchLoadDataBuffers( region, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: iReg
  END SUBROUTINE PLAG_PatchLoadDataBuffers 

  SUBROUTINE PLAG_PatchRemoveDataOutflow( region, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: iReg
  END SUBROUTINE PLAG_PatchRemoveDataOutflow 
  
  SUBROUTINE PLAG_PatchUpdate( regions, iReg )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region), POINTER :: regions(:)
    INTEGER        :: iReg
  END SUBROUTINE PLAG_PatchUpdate
  
  SUBROUTINE PLAG_ReadDisPartSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_ReadDisPartSection

  SUBROUTINE PLAG_ReadDisPartInitSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_ReadDisPartInitSection

  SUBROUTINE PLAG_ReadDisPartnContSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_ReadDisPartnContSection

  SUBROUTINE PLAG_ReadInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_ReadInputFile

  SUBROUTINE PLAG_ReallocateDataBuffers( regions, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE PLAG_ReallocateDataBuffers 
      
  SUBROUTINE PLAG_ReflectParticleData(pPatch,pPlag,ifl,iPcl,xLocOld,yLocOld, & 
                                      zLocOld,xLoc,yLoc,zLoc,xTraj,yTraj,zTraj)   
    USE ModDataTypes
    USE ModBndPatch, ONLY: t_patch
    USE ModPartLag, ONLY: t_plag    
    INTEGER, INTENT(IN) :: ifl,iPcl
    REAL(RFREAL), INTENT(INOUT) :: xLoc,xLocOld,yLoc,yLocOld,zLoc,zLocOld
    REAL(RFREAL), INTENT(INOUT), OPTIONAL :: xTraj,yTraj,zTraj                                   
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_plag), POINTER :: pPlag    
  END SUBROUTINE PLAG_ReflectParticleData 
  
  SUBROUTINE PLAG_ReadPdfFromFile(regions, brbeg, brend )
    USE ModDataStruct, ONLY : t_region 
    TYPE(t_region), POINTER :: regions(:)
    INTEGER,INTENT(IN)      :: brbeg,brend
  END SUBROUTINE PLAG_ReadPdfFromFile

  SUBROUTINE PLAG_RKUpdate( region, iStage )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region), INTENT(INOUT) :: region
    INTEGER,        INTENT(IN)    :: iStage
  END SUBROUTINE PLAG_RKUpdate 
  
  SUBROUTINE PLAG_RkUpdateGeneric(region,iStage,ipBeg,ipEnd,ivBeg,ivEnd, &
                                  cv,cvOld,rhs,rhsSum)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region
    INTEGER, INTENT(IN) :: ipBeg,ipEnd,iStage,ivBeg,ivEnd
    REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,cvOld,rhs,rhsSum    
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_RkUpdateGeneric
  
  SUBROUTINE PLAG_Update( region, iReg, iStage )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region), INTENT(INOUT) :: region
    INTEGER,        INTENT(IN)    :: iReg, iStage
  END SUBROUTINE PLAG_Update
  
#ifdef RFLO
  SUBROUTINE PLAG_UpdateDataStruct( region, iReg )
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region) :: region
    INTEGER, INTENT(IN)  :: iReg
  END SUBROUTINE PLAG_UpdateDataStruct
#endif 
  
  SUBROUTINE PLAG_WallBounce( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_WallBounce
  
  SUBROUTINE PLAG_ZeroRhs( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE PLAG_ZeroRhs

! -----------------------------------------------------------------------------
! RFLU-specific routines
! -----------------------------------------------------------------------------
  
#ifdef RFLU
  SUBROUTINE PLAG_RFLU_InitSolSerial_1D(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_InitSolSerial_1D

  SUBROUTINE PLAG_RFLU_InitSolSerial_3D(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_InitSolSerial_3D

  SUBROUTINE PLAG_RFLU_InjcTileCalcRhs(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_InjcTileCalcRhs

  SUBROUTINE PLAG_RFLU_InjectionDriver(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion  
  END SUBROUTINE PLAG_RFLU_InjectionDriver
  
  SUBROUTINE PLAG_RFLU_Update( pRegion, iStage )
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region), POINTER :: pRegion
    INTEGER,    INTENT(IN)  :: iStage
  END SUBROUTINE PLAG_RFLU_Update

  SUBROUTINE PLAG_UpdateDataStruct( pRegion )
    USE ModDataStruct, ONLY : t_region    
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE PLAG_UpdateDataStruct
#endif 
  
  END INTERFACE

END MODULE PLAG_ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModInterfaces.F90,v $
! Revision 1.48  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.47  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.46  2007/03/15 22:00:18  haselbac
! Added IFs for PLAG_RFLU_InitSolSerial_{1,3}D.F90
!
! Revision 1.45  2005/05/19 16:01:40  fnajjar
! Removed Interfaces for obsolete routine calls
!
! Revision 1.44  2005/04/27 15:00:08  fnajjar
! Removed interface calls to individual FindCells routines
!
! Revision 1.43  2005/04/25 18:39:10  luca1
! Imposed PDF from file option for random particle ejection
!
! Revision 1.42  2005/03/11 02:23:06  haselbac
! Adapted interfaces for new PLAG_RFLU_FindCellsTrajXYZ routines
!
! Revision 1.41  2005/01/01 21:33:12  haselbac
! Added interface for Apte method, modif interface for PLAG_ReflectParticleData
!
! Revision 1.40  2004/11/05 21:51:01  fnajjar
! Added interfaces to particle-cell search routines
!
! Revision 1.39  2004/11/05 21:09:14  haselbac
! Fixed bad check-in
!
! Revision 1.38  2004/11/05 20:33:57  haselbac
! Adapted interface for PLAG_ReflectParticleData
!
! Revision 1.37  2004/10/11 22:09:22  haselbac
! Renamed procedures
!
! Revision 1.36  2004/10/08 22:09:26  haselbac
! Added entry for PLAG_RFLU_FindParticleCellsBrut
!
! Revision 1.35  2004/08/20 23:27:13  fnajjar
! Added Infrastructure for Plag prep tool
!
! Revision 1.34  2004/06/16 22:58:48  fnajjar
! Renamed injcModel to injcDiamDist for CRE kernel
!
! Revision 1.33  2004/04/09 22:57:45  fnajjar
! Added Interfaces for RFLO-specific routines
!
! Revision 1.32  2004/04/08 01:33:00  haselbac
! Added interface for PLAG_ReflectParticleData
!
! Revision 1.31  2004/03/26 21:26:26  fnajjar
! Added new routines for RFLU-specific routines
!
! Revision 1.30  2004/03/15 21:05:54  haselbac
! Deleted/added interface
!
! Revision 1.29  2004/03/08 22:17:22  fnajjar
! Included interface calls for RFLU-specific injection routines
!
! Revision 1.28  2004/03/05 23:17:05  haselbac
! Added interface for PLAG_RFLU_FindParticleCellsTraj
!
! Revision 1.27  2004/02/26 21:02:12  haselbac
! Added common and RFLU-specific interfaces
!
! Revision 1.26  2004/02/25 21:56:45  fnajjar
! Included generic RKUpdate for PLAG
!
! Revision 1.25  2004/02/06 21:19:36  fnajjar
! Included proper INTENT to Interfaces
!
! Revision 1.24  2004/01/15 21:10:52  fnajjar
! Separated Interfaces for corner-edge cell routines
!
! Revision 1.23  2003/11/21 22:44:35  fnajjar
! Changed global to region for updated Random Number Generator
!
! Revision 1.22  2003/11/12 21:34:36  fnajjar
! Added Corner-Edge cells Interface calls
!
! Revision 1.21  2003/11/03 21:22:20  fnajjar
! Added PLAG_copyFaceVectors
!
! Revision 1.20  2003/09/13 20:14:21  fnajjar
! Added infrastructure for Breakup model
!
! Revision 1.19  2003/05/28 15:16:11  fnajjar
! Removed obsolete PLAG_mixt calls as embedded in Rocinteract
!
! Revision 1.18  2003/04/17 00:12:08  fnajjar
! Added INTENT(IN) for ijkNR, ijkNRI, ijkNRJ, ijkNRK
!
! Revision 1.17  2003/04/14 19:02:57  fnajjar
! Bug fix to include POINTER attribute for PLAG_PatchUpdate
!
! Revision 1.16  2003/04/14 14:32:20  fnajjar
! Added PLAG_initInputValues for proper initialization
!
! Revision 1.15  2003/03/28 19:50:39  fnajjar
! Moved interfaces for wrapper routines
!
! Revision 1.14  2003/02/21 17:06:55  fnajjar
! Added Interfaces to Data Send and Recv
!
! Revision 1.13  2003/01/24 22:07:57  f-najjar
! Include Interface call to PLAG_ClearRequests
!
! Revision 1.12  2003/01/23 17:12:14  f-najjar
! Included Interface for PLAG_BufferSizeRecv and PLAG_BufferSizeSend
!
! Revision 1.11  2003/01/23 17:05:06  f-najjar
! Moved Interface call to PLAG_patchBufferSendRecv into ModInterfacesLagrangian
!
! Revision 1.10  2003/01/17 19:32:20  f-najjar
! Included correct INTERFACE for PLAG_PatchRemoveDataOutflow and PLAG_PatchLoadDataBuffers
!
! Revision 1.9  2003/01/16 22:46:18  f-najjar
! Include iReg to calling sequence for PLAG_GetCellIndices
!
! Revision 1.8  2003/01/16 22:30:10  f-najjar
! Included INTERFACE call for PLAG_PatchBufferSendRecv
!
! Revision 1.7  2003/01/16 22:28:40  f-najjar
! Included INTERFACE call for PLAG_AppendDataFromBuffers
!
! Revision 1.6  2003/01/13 19:00:33  f-najjar
! Removed PLAG_allocateDataBuffers
!
! Revision 1.5  2003/01/10 19:37:40  f-najjar
! Included Interface call to PLAG_BoundaryConditionsSet
!
! Revision 1.4  2003/01/10 19:28:04  f-najjar
! Included Interface call for PLAG_PatchExchangeConf
!
! Revision 1.3  2003/01/10 19:09:10  f-najjar
! Included iReg in calling sequence
!
! Revision 1.2  2002/10/25 14:11:12  f-najjar
! Included interface calls to PLAG subroutines
!
! Revision 1.1.1.1  2001/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************






