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
! Purpose: set explicit interfaces to subroutines and functions
!          for corner and edge cells.
!
! Description: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_ModInterfacesCECells.F90,v 1.7 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************
  
MODULE PLAG_ModInterfacesCECells

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! PLAG specific library
! =============================================================================
 
  SUBROUTINE PLAG_CECellsAllocateData( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_CECellsAllocateData

  SUBROUTINE PLAG_CECellsClearRequestsData( regions, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_CECellsClearRequestsData

  SUBROUTINE PLAG_CECellsClearRequestsSize( regions, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_CECellsClearRequestsSize

  SUBROUTINE PLAG_CECellsDeallocateData( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_CECellsDeallocateData 

  SUBROUTINE PLAG_CECellsExchange( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_CECellsExchange

  SUBROUTINE PLAG_CECellsFaceCentroids( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_CECellsFaceCentroids
        
  SUBROUTINE PLAG_CECellsFaceVectors( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_CECellsFaceVectors
        
  SUBROUTINE PLAG_CECellsGetBufferSize( region,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)      :: region
    INTEGER, INTENT(IN) :: iReg
  END SUBROUTINE PLAG_CECellsGetBufferSize

  SUBROUTINE PLAG_CECellsLoadDataWrapper( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_CECellsLoadDataWrapper

  SUBROUTINE PLAG_CECellsRecvData( regions, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_CECellsRecvData

  SUBROUTINE PLAG_CECellsRecvSize( regions, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_CECellsRecvSize

  SUBROUTINE PLAG_CECellsSendData( regions, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_CECellsSendData

  SUBROUTINE PLAG_CECellsSendRecvWrapper( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_CECellsSendRecvWrapper

  SUBROUTINE PLAG_CECellsSendSize( regions, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_CECellsSendSize

  SUBROUTINE PLAG_CornCellsLoadData( regions, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_CornCellsLoadData

  SUBROUTINE PLAG_CornCellsLoadSendBuff( regions, iReg, ir, nBuffSizeEdge, &
                                         nBuffSizeCorn )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg, ir
    INTEGER, INTENT(IN)     :: nBuffSizeEdge
    INTEGER, INTENT(OUT)    :: nBuffSizeCorn
  END SUBROUTINE PLAG_CornCellsLoadSendBuff

  SUBROUTINE PLAG_EdgeCellsLoadData( regions, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_EdgeCellsLoadData

  SUBROUTINE PLAG_EdgeCellsLoadSendBuff( regions, iReg, ir, nBuffSizeEdge )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg, ir
    INTEGER, INTENT(OUT)    :: nBuffSizeEdge
  END SUBROUTINE PLAG_EdgeCellsLoadSendBuff

  SUBROUTINE PLAG_RFLO_ClearSendRequests( regions, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_RFLO_ClearSendRequests

  SUBROUTINE PLAG_RFLO_FindGridMapping( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_RFLO_FindGridMapping

  SUBROUTINE PLAG_RFLO_FindSourceCell( regions,iReg,iLev,ic,jc,kc,&
                                       found,iRegSrc,indexMapMat )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg, iLev
    INTEGER, INTENT(INOUT)  :: ic, jc, kc
    INTEGER, INTENT(OUT)    :: iRegSrc, indexMapMat(3,4)
    LOGICAL, INTENT(OUT)    :: found
  END SUBROUTINE PLAG_RFLO_FindSourceCell

  SUBROUTINE PLAG_RFLO_GetFaceMapping( mapMat, srcDir, srcFace )
    INTEGER, INTENT(IN)  :: mapMat(3,4)
    INTEGER, INTENT(OUT) :: srcDir(3)
    INTEGER, INTENT(OUT) :: srcFace(6)
  END SUBROUTINE PLAG_RFLO_GetFaceMapping

  SUBROUTINE PLAG_RFLO_RecvMetrics( regions, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_RFLO_RecvMetrics

  SUBROUTINE PLAG_RFLO_SendMetrics( regions, iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER, INTENT(IN)     :: iReg
  END SUBROUTINE PLAG_RFLO_SendMetrics
  
  SUBROUTINE PLAG_RFLO_SendRecvMetrics( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PLAG_RFLO_sendRecvMetrics

  SUBROUTINE PLAG_RFLO_SourceCell( region,regionSrc,patch,patchSrc, &
                                   iLev,ic,jc,kc,found,indexMapMat )
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch,   ONLY : t_patch
    TYPE(t_region)          :: region, regionSrc
    TYPE(t_patch), POINTER  :: patch, patchSrc
    INTEGER, INTENT(IN)     :: iLev
    INTEGER, INTENT(INOUT)  :: ic,jc,kc
    INTEGER, INTENT(OUT)    :: indexMapMat(3,4)
    LOGICAL, INTENT(INOUT)  :: found
  END SUBROUTINE PLAG_RFLO_SourceCell
 
  END INTERFACE

END MODULE PLAG_ModInterfacesCECells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModInterfacesCECells.F90,v $
! Revision 1.7  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2004/03/19 23:49:12  fnajjar
! Added interface call to loading buffers for MPI-based communications
!
! Revision 1.4  2004/03/10 23:12:54  fnajjar
! Included interfaces for MPI-based routines for corner-edge cells
!
! Revision 1.3  2004/02/10 21:23:08  fnajjar
! Added call and interfaces for index mapping between corner-edge regions
!
! Revision 1.2  2004/01/26 22:55:31  fnajjar
! Included calls to routines for corner-edge data loading
!
! Revision 1.1  2004/01/15 21:17:37  fnajjar
! Initial Import of Interface for corner-edge cell routines
!
!******************************************************************************






