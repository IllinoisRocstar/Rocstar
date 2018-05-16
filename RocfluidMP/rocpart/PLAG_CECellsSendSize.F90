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
! Purpose: send buffer size to edge and corner cells of an adjacent region.
!
! Description: kernel is pertinent when the other region is located
!              on a different processor.
!
! Input: 
!   regions = data of all regions
!   iReg    = current region.
!
! Output: buffer size sent.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_CECellsSendSize.F90,v 1.5 2009/03/02 00:19:36 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CECellsSendSize( regions,iReg )

  USE ModDataTypes
  USE ModError
  USE ModMPI
  USE ModParameters
  USE ModGlobal,     ONLY : t_global
  USE ModDataStruct, ONLY : t_dCell, t_dCellTransf, t_region, t_level 
  USE ModPartLag,    ONLY : t_plag, t_buffer_plag

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)

  INTEGER, INTENT(IN)     :: iReg

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: dest,i,j,k,icell,iCorner,iEdge,ijk,iLev,ir,iRegDes, iRegSrc,    &
             iRequestPlag,nBuffSizeCornSrc,nBuffSizeEdgeSrc,nBuffSizeTotSrc, &
             nCorners,nDimBuffSize,nEdges,tag

  TYPE(t_buffer_plag), POINTER :: pCornCellsXBuff, pEdgeCellsXBuff
  TYPE(t_dCellTransf), POINTER :: pSendEcCell
  TYPE(t_global),      POINTER :: global
  TYPE(t_level),       POINTER :: pLevelSrc
  TYPE(t_plag),        POINTER :: pPlag
  TYPE(t_region),      POINTER :: pRegionSrc, pRegionDes 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CECellsSendSize.F90,v $ $Revision: 1.5 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PLAG_CECellsSendSize',&
  'PLAG_CECellsSendSize.F90' )

! ******************************************************************************
! Get dimensions 
! ******************************************************************************

  iLev = regions(iReg)%currLevel
  nCorners  = 8
  nEdges    = 12
  
  nDimBuffSize = 1 

  nBuffSizeCornSrc = 0
  nBuffSizeEdgeSrc = 0
  nBuffSizeTotSrc  = 0

! ******************************************************************************
! Set pointers 
! ******************************************************************************
  
  pRegionSrc => regions(iReg)
  pLevelSrc  => pRegionSrc%levels(iLev)
  pPlag      => pLevelSrc%plag

! ******************************************************************************
! Load send buffer size
! ******************************************************************************

  DO ir=1,global%nRegions
    IF (regions(ir)%procid == global%myProcid) GOTO 1999

    IF ( pLevelSrc%sendEcCells(ir)%nCells > 0 ) THEN
      pSendEcCell => pLevelSrc%sendEcCells(ir)

! =============================================================================
!     Reset buffer size in active communicating regions to zero
!      since RFLO has possibility to have multiple active edges/corners
!      on various regions interacting; while for PLAG, edges/corners interact
!      dynamically when buffer size is non-null.
! =============================================================================

      nBuffSizeCornSrc = 0
      nBuffSizeEdgeSrc = 0     
      nBuffSizeTotSrc  = 0

! =============================================================================
!     Loop over edges of source region
!      Loading buffer size for edge
! =============================================================================

      DO iEdge=1,nEdges
        IF( .NOT. pLevelSrc%edgeCells(iEdge)%interact ) GOTO 2999
 
! -- Bypass for degenerate edge cells -----------------------------------------

        IF( pLevelSrc%edgeCells(iEdge)%degenrt /= DEGENERAT_NONE ) GOTO 2999
             
        DO ijk=1,UBOUND(pLevelSrc%edgeCells(iEdge)%cells,1)
          iRegDes = pLevelSrc%edgeCells(iEdge)%cells(ijk)%srcRegion
          pEdgeCellsXBuff => pLevelSrc%edgeCells(iEdge)%cells(ijk)%bufferExchPlag
          IF ( iRegDes == ir .AND. regions(iRegDes)%procid /= global%myProcid ) THEN               
            nBuffSizeEdgeSrc = nBuffSizeEdgeSrc +pEdgeCellsXBuff%nBuffSize 
          ENDIF   ! iRegDes
        ENDDO     ! ijk

#ifdef PLAG_CECELLS_MPI_DEBUG
  IF ( nBuffSizeEdgeSrc > 0 ) & 
    PRINT*,' PLAG_CECellsSendSize: procId, iReg, iR, procIdiR, iEdge, nBuffSizeEdgeSrc,iRegDes  = ',&
       global%myProcid, iReg, iR, regions(ir)%procid ,iEdge, nBuffSizeEdgeSrc,iRegDes
#endif

2999    CONTINUE
      ENDDO       ! iEdge

! =============================================================================
!     Loop over corners of source region
!      Loading buffer size for corner
! =============================================================================

      DO iCorner=1,nCorners
        IF (.NOT. pLevelSrc%cornerCells(iCorner)%interact) GOTO 3999

! -- Bypass for degenerate corner cells ---------------------------------------

        IF( pLevelSrc%cornerCells(iCorner)%degenrt /= DEGENERAT_NONE ) GOTO 3999
        
        DO ijk=1,UBOUND(pLevelSrc%cornerCells(iCorner)%cells,1)
          iRegDes = pLevelSrc%cornerCells(iCorner)%cells(ijk)%srcRegion
          pCornCellsXBuff => pLevelSrc%cornerCells(iCorner)%cells(ijk)%bufferExchPlag
          IF ( iRegDes == ir .AND. regions(iRegDes)%procid /= global%myProcid ) THEN               
            nBuffSizeCornSrc = nBuffSizeCornSrc +pCornCellsXBuff%nBuffSize 
          ENDIF   ! iRegDes 
        ENDDO     ! ijk

#ifdef PLAG_CECELLS_MPI_DEBUG
  IF ( nBuffSizeCornSrc > 0 ) &
    PRINT*,' PLAG_CECellsSendSize : procId, iReg, iR, iCorner, nBuffSizeCornSrc,iRegDes   = ',&
     global%myProcid, iReg, iR, regions(ir)%procid ,iCorner, nBuffSizeCornSrc,iRegDes
#endif

3999    CONTINUE
      ENDDO         ! iCorner

! =============================================================================
!     Send buffers to destination processor 
! =============================================================================

      nBuffSizeTotSrc           = nBuffSizeCornSrc +nBuffSizeEdgeSrc
      pSendEcCell%nBuffSizePlag = nBuffSizeTotSrc 
      iRequestPlag              = pSendEcCell%iRequestPlag

#ifdef MPI
      dest = regions(ir)%procid
      tag  = regions(ir)%localNumber +PLAG_TAG_SHIFT +MPI_PATCHOFF +1000

      IF(tag .gt. global%mpiTagMax) tag = MOD(tag,global%mpiTagMax)
#ifdef PLAG_CECELLS_MPI_DEBUG 
  IF ( nBuffSizeTotSrc > 0 ) &             
   PRINT*,'  PLAG_CECellsSendSize: iRegDes, iRegSrc, procDes, procSrc, tagSrc, nBuffSizePlag  = ',&
    ir, iReg, dest, global%myProcid,tag, pSendEcCell%nBuffSizePlag
#endif

      CALL MPI_Isend( pSendEcCell%nBuffSizePlag,nDimBuffSize,MPI_INTEGER, &
                      dest,tag,global%mpiComm,                            &
                      pPlag%requestsCECells(iRequestPlag),global%mpierr   )
      IF ( global%mpierr /= ERR_NONE ) &
        CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

    ENDIF      ! some cells to send

1999 CONTINUE
  ENDDO        ! ir

! ******************************************************************************
! finalize
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CECellsSendSize

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CECellsSendSize.F90,v $
! Revision 1.5  2009/03/02 00:19:36  mtcampbe
! Added some ifdefs around Rocflo to disable particle injection on INFLOW
! boundaries and added some checks around MPI tags utilizing a new global
! data item, global%mpiTagMax.
!
! Revision 1.4  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 20:57:20  fnajjar
! Initial revision after changing case
!
! Revision 1.4  2004/11/29 19:26:25  fnajjar
! Added bypass statement for dengerate cells
!
! Revision 1.3  2004/03/20 21:33:20  fnajjar
! Reset buffer sizes in active communicating regions to zero
!
! Revision 1.2  2004/03/10 23:42:31  fnajjar
! Activated IF-statement based on buffer size in ifdef construct
!
! Revision 1.1  2004/03/10 23:16:09  fnajjar
! Initial import of routines to MPI-communicate buffer sizes
!
!******************************************************************************







