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
! Purpose: send RFLO metrics to edge and corner cells of an adjacent region.
!
! Description: kernel is pertinent when the other region is located
!              on a different processor.
!
! Input: regions = data of all regions
!        iReg    = current region.
!
! Output: new values of RFLO metrics variables.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_RFLO_SendMetrics.F90,v 1.5 2009/03/02 00:19:36 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLO_SendMetrics( regions,iReg )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region, t_level, t_dCellTransf
  USE ModPartLag, ONLY    : t_plag
  USE ModError
  USE ModMPI
  USE ModParameters
  USE ModIndexing, ONLY   : GetIJK
  USE ModInterfaces, ONLY : RFLO_GetCellOffset,         &
                            RFLO_GetCornerCellsIndices, &
                            RFLO_GetEdgeCellsIndices,   &
                            RFLO_GetNodeOffset
  USE PLAG_ModInterfaces, ONLY : PLAG_RFLO_GetFaceMapping
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER, INTENT(IN)     :: iReg

! ... loop variables
  INTEGER :: i,j,k, iCorner, iEdge, iFace, ijk, ir, ldir 

! ... local variables
  INTEGER :: iLev, iRegSrc, icell, ibuff, nDim, dest, tag
  INTEGER :: nBuffSize, nCorners, nEdges, nDir
  INTEGER :: nFaces, nFaceCentroidSize, nFaceNormalSize
  INTEGER :: iCOffSrc, ijCOffSrc
  INTEGER :: ijkCCSrc, ijkECSrc
  INTEGER :: iCCSrc, jCCSrc, kCCSrc
  INTEGER :: iECSrc, jECSrc, kECSrc
  INTEGER :: nDumCellsSrc
  INTEGER :: iNOff, ijNOff
  INTEGER :: iNOffSrc, ijNOffSrc, ijkN, ijkCNSrc, ijkENSrc
  INTEGER :: iCNSrc, jCNSrc, kCNSrc
  INTEGER :: iENSrc, jENSrc, kENSrc
  INTEGER :: nodeCornSrc, nodeEdgeSrc
  INTEGER :: idirSrc, jdirSrc, kdirSrc, ldirSrc, iFaceSrc
  INTEGER :: srcDir(3),srcFace(6)
  INTEGER :: srcIndexMapMat(3,4)
 
  REAL(RFREAL), DIMENSION(3,3)            :: sFaceCorn, sFaceEdge 
  REAL(RFREAL), POINTER, DIMENSION(:,:)   :: pSi, pSj, pSk
  REAL(RFREAL), POINTER, DIMENSION(:,:,:) :: pFc
 
  TYPE(t_region),      POINTER :: pRegion 
  TYPE(t_level),       POINTER :: pLevel
  TYPE(t_dCellTransf), POINTER :: pSendEcCell
  TYPE(t_plag),        POINTER :: pPlag
  TYPE(t_global),      POINTER :: global
  
!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PLAG_RFLO_SendMetrics',&
  'PLAG_RFLO_SendMetrics.F90' )

! Get dimensions --------------------------------------------------------------

  iLev = regions(iReg)%currLevel
  nCorners  = 8
  nEdges    = 12
  nFaces    = 6
  nDir      = 3
  
! Compute buffer size has to store information --------------------------------
!  on the total number of faces for each cell .
  
  nFaceCentroidSize = ZCOORD*KCOORD
  nFaceNormalSize   =      3*KCOORD
  nBuffSize = (nFaceCentroidSize + nFaceNormalSize)

! Set pointers ----------------------------------------------------------------
  
  pRegion => regions(iReg)
  pPlag   => regions(iReg)%levels(iLev)%plag
  pFc     => regions(iReg)%levels(iLev)%plag%fc
  pSi     => regions(iReg)%levels(iLev)%plag%si
  pSj     => regions(iReg)%levels(iLev)%plag%sj
  pSk     => regions(iReg)%levels(iLev)%plag%sk

! Get node offset -------------------------------------------------------------

  CALL RFLO_GetNodeOffset( pRegion,iLev,iNOff,ijNOff ) 

! Fill send buffers -----------------------------------------------------------

  DO ir=1,global%nRegions
    IF (regions(ir)%procid /= global%myProcid) THEN
    IF (regions(iReg)%levels(iLev)%sendEcCells(ir)%nCells > 0) THEN

      pSendEcCell => regions(iReg)%levels(iLev)%sendEcCells(ir)
      pLevel      => regions(ir)%levels(iLev)
      nDim        =  pSendEcCell%nCells*nFaces
      ibuff       =  0

! - Load edges ----------------------------------------------------------------

      DO iEdge=1,nEdges

! -- Bypass for noninteracting regions ----------------------------------------

        IF( .NOT. pLevel%edgeCells(iEdge)%interact ) GOTO 1999

! -- Bypass for degenerate edge cells -----------------------------------------

        IF( pLevel%edgeCells(iEdge)%degenrt /= DEGENERAT_NONE ) GOTO 1999

! -- Source region infrastructure ---------------------------------------------

        DO ijk=1,UBOUND(pLevel%edgeCells(iEdge)%cells,1)
          iRegSrc      = pLevel%edgeCells(iEdge)%cells(ijk)%srcRegion

          IF ( iRegSrc == iReg ) THEN
            ijkECSrc     = pLevel%edgeCells(iEdge)%cells(ijk)%srcCell          
            nDumCellsSrc = regions(iRegSrc)%nDumCells

            CALL RFLO_GetCellOffset( regions(iRegSrc),iLev,iCOffSrc,ijCOffSrc )
            CALL GetIJK( ijkECSrc,iCOffSrc,ijCOffSrc,nDumCellsSrc, &
                         iECSrc,jECSrc,kECSrc )

            CALL RFLO_GetNodeOffset( regions(iRegSrc),iLev,iNOffSrc,ijNOffSrc )
            ijkENSrc = IndIJK(iECSrc,jECSrc,kECSrc,iNOffSrc,ijNOffSrc)

! -- Find face mapping for source region --------------------------------------
        
            srcIndexMapMat= pLevel%edgeCells(iEdge)%cells(ijk)%srcIndexMapMat 

            CALL PLAG_RFLO_GetFaceMapping(srcIndexMapMat,srcDir,srcFace)

            idirSrc = srcDir(1); jdirSrc = srcDir(2); kdirSrc = srcDir(3);
         
! --- Loop over all the Faces -------------------------------------------------
! --- Note: assumption of coordinate alignment between regions removed --------

            DO iFace = 1, nFaces
              ibuff = ibuff + 1

              iFaceSrc = srcFace(iFace)
              SELECT CASE (iFaceSrc)
                CASE(1,3,5)
                  nodeEdgeSrc = IndIJK(iECSrc,jECSrc,kECSrc,iNOffSrc,ijNOffSrc)

                CASE(2)
                  nodeEdgeSrc = IndIJK(iECSrc+1,jECSrc,kECSrc,iNOffSrc,ijNOffSrc)

                CASE(4)
                  nodeEdgeSrc = IndIJK(iECSrc,jECSrc+1,kECSrc,iNOffSrc,ijNOffSrc)

                CASE(6)
                  nodeEdgeSrc = IndIJK(iECSrc,jECSrc,kECSrc+1,iNOffSrc,ijNOffSrc)
              END SELECT ! iFaceSrc

! --- Determine direction mapping for face normals ----------------------------

              DO ldir = 1,nDir
                ldirSrc=srcDir(ldir)
                SELECT CASE (ldirSrc)
                  CASE(ICOORD)
                    sFaceEdge(XCOORD:ZCOORD,ldirSrc) = pSi(XCOORD:ZCOORD,nodeEdgeSrc)
                
                  CASE(JCOORD)
                    sFaceEdge(XCOORD:ZCOORD,ldirSrc) = pSj(XCOORD:ZCOORD,nodeEdgeSrc)
                
                  CASE(KCOORD)
                    sFaceEdge(XCOORD:ZCOORD,ldirSrc) = pSk(XCOORD:ZCOORD,nodeEdgeSrc)
                END SELECT ! ldirSrc
              ENDDO ! ldir

! ---- Load face centroids -----------------------------------------------------

              pSendEcCell%buffMetrics(ibuff        ) = pFc(XCOORD,idirSrc,nodeEdgeSrc)
              pSendEcCell%buffMetrics(ibuff+   nDim) = pFc(XCOORD,jdirSrc,nodeEdgeSrc)
              pSendEcCell%buffMetrics(ibuff+ 2*nDim) = pFc(XCOORD,kdirSrc,nodeEdgeSrc)

              pSendEcCell%buffMetrics(ibuff+ 3*nDim) = pFc(YCOORD,idirSrc,nodeEdgeSrc)
              pSendEcCell%buffMetrics(ibuff+ 4*nDim) = pFc(YCOORD,jdirSrc,nodeEdgeSrc)
              pSendEcCell%buffMetrics(ibuff+ 5*nDim) = pFc(YCOORD,kdirSrc,nodeEdgeSrc)

              pSendEcCell%buffMetrics(ibuff+ 6*nDim) = pFc(ZCOORD,idirSrc,nodeEdgeSrc)
              pSendEcCell%buffMetrics(ibuff+ 7*nDim) = pFc(ZCOORD,jdirSrc,nodeEdgeSrc)
              pSendEcCell%buffMetrics(ibuff+ 8*nDim) = pFc(ZCOORD,kdirSrc,nodeEdgeSrc)

! --- Load face normals -------------------------------------------------------

              pSendEcCell%buffMetrics(ibuff+ 9*nDim) = sFaceEdge(XCOORD,ICOORD)
              pSendEcCell%buffMetrics(ibuff+10*nDim) = sFaceEdge(YCOORD,ICOORD)
              pSendEcCell%buffMetrics(ibuff+11*nDim) = sFaceEdge(ZCOORD,ICOORD)

              pSendEcCell%buffMetrics(ibuff+12*nDim) = sFaceEdge(XCOORD,JCOORD)
              pSendEcCell%buffMetrics(ibuff+13*nDim) = sFaceEdge(YCOORD,JCOORD)
              pSendEcCell%buffMetrics(ibuff+14*nDim) = sFaceEdge(ZCOORD,JCOORD)

              pSendEcCell%buffMetrics(ibuff+15*nDim) = sFaceEdge(XCOORD,KCOORD)
              pSendEcCell%buffMetrics(ibuff+16*nDim) = sFaceEdge(YCOORD,KCOORD)
              pSendEcCell%buffMetrics(ibuff+17*nDim) = sFaceEdge(ZCOORD,KCOORD)

            ENDDO ! iFace  

          ENDIF   ! iRegSrc
        ENDDO     ! ijk

1999    CONTINUE
      ENDDO       ! iEdge

! - Load corners --------------------------------------------------------------

      DO iCorner=1,nCorners
        IF (.NOT. pLevel%cornerCells(iCorner)%interact) GOTO 2999

! -- Bypass for degenerate corner cells ---------------------------------------

        IF( pLevel%cornerCells(iCorner)%degenrt /= DEGENERAT_NONE ) GOTO 2999

! -- Source region infrastructure ---------------------------------------------
        
          DO ijk=1,UBOUND(pLevel%cornerCells(iCorner)%cells,1)
            iRegSrc      = pLevel%cornerCells(iCorner)%cells(ijk)%srcRegion
            IF ( iRegSrc == iReg ) THEN
              ijkCCSrc     = pLevel%cornerCells(iCorner)%cells(ijk)%srcCell
              nDumCellsSrc = regions(iRegSrc)%nDumCells

              CALL RFLO_GetCellOffset( regions(iRegSrc),iLev,iCOffSrc,ijCOffSrc )
              CALL GetIJK( ijkCCSrc,iCOffSrc,ijCOffSrc,nDumCellsSrc, &
                           iCCSrc,jCCSrc,kCCSrc )

              CALL RFLO_GetNodeOffset( regions(iRegSrc),iLev,iNOffSrc,ijNOffSrc )
              ijkCNSrc = IndIJK(iCCSrc,jCCSrc,kCCSrc,iNOffSrc,ijNOffSrc)

! -- Find face mapping for source region --------------------------------------
        
              srcIndexMapMat= pLevel%cornerCells(iCorner)%cells(ijk)%srcIndexMapMat 

              CALL PLAG_RFLO_GetFaceMapping(srcIndexMapMat,srcDir,srcFace)

              idirSrc = srcDir(1); jdirSrc = srcDir(2); kdirSrc = srcDir(3);

! --- Loop over all the Faces -------------------------------------------------
! --- Note: assumption of coordinate alignment between regions removed --------

              DO iFace = 1, nFaces
                ibuff = ibuff + 1 

                iFaceSrc = srcFace(iFace)
                SELECT CASE (iFaceSrc)
                  CASE(1,3,5)
                    nodeCornSrc = IndIJK(iCCSrc,jCCSrc,kCCSrc,iNOffSrc,ijNOffSrc)

                  CASE(2)
                    nodeCornSrc = IndIJK(iCCSrc+1,jCCSrc,kCCSrc,iNOffSrc,ijNOffSrc)

                  CASE(4)
                     nodeCornSrc = IndIJK(iCCSrc,jCCSrc+1,kCCSrc,iNOffSrc,ijNOffSrc)

                  CASE(6)
                    nodeCornSrc = IndIJK(iCCSrc,jCCSrc,kCCSrc+1,iNOffSrc,ijNOffSrc)
                END SELECT  ! iFaceSrc            

! --- Determine direction mapping for face normals ----------------------------

                DO ldir = 1,nDir
                  ldirSrc=srcDir(ldir)
                  SELECT CASE (ldirSrc)
                    CASE(ICOORD)
                      sFaceCorn(XCOORD:ZCOORD,ldirSrc) = pSi(XCOORD:ZCOORD,nodeCornSrc)
                
                    CASE(JCOORD)
                      sFaceCorn(XCOORD:ZCOORD,ldirSrc) = pSj(XCOORD:ZCOORD,nodeCornSrc)
                
                    CASE(KCOORD)
                      sFaceCorn(XCOORD:ZCOORD,ldirSrc) = pSk(XCOORD:ZCOORD,nodeCornSrc)
                  END SELECT ! ldirSrc
                ENDDO ! ldir
            
! --- Load face centroids -----------------------------------------------------

                pSendEcCell%buffMetrics(ibuff        ) = pFc(XCOORD,idirSrc,nodeCornSrc)
                pSendEcCell%buffMetrics(ibuff+   nDim) = pFc(XCOORD,jdirSrc,nodeCornSrc)
                pSendEcCell%buffMetrics(ibuff+ 2*nDim) = pFc(XCOORD,kdirSrc,nodeCornSrc)

                pSendEcCell%buffMetrics(ibuff+ 3*nDim) = pFc(YCOORD,idirSrc,nodeCornSrc)
                pSendEcCell%buffMetrics(ibuff+ 4*nDim) = pFc(YCOORD,jdirSrc,nodeCornSrc)
                pSendEcCell%buffMetrics(ibuff+ 5*nDim) = pFc(YCOORD,kdirSrc,nodeCornSrc)

                pSendEcCell%buffMetrics(ibuff+ 6*nDim) = pFc(ZCOORD,idirSrc,nodeCornSrc)
                pSendEcCell%buffMetrics(ibuff+ 7*nDim) = pFc(ZCOORD,jdirSrc,nodeCornSrc)
                pSendEcCell%buffMetrics(ibuff+ 8*nDim) = pFc(ZCOORD,kdirSrc,nodeCornSrc)

! --- Load face normals -------------------------------------------------------

                pSendEcCell%buffMetrics(ibuff+ 9*nDim) = sFaceCorn(XCOORD,ICOORD)
                pSendEcCell%buffMetrics(ibuff+10*nDim) = sFaceCorn(YCOORD,ICOORD)
                pSendEcCell%buffMetrics(ibuff+11*nDim) = sFaceCorn(ZCOORD,ICOORD)

                pSendEcCell%buffMetrics(ibuff+12*nDim) = sFaceCorn(XCOORD,JCOORD)
                pSendEcCell%buffMetrics(ibuff+13*nDim) = sFaceCorn(YCOORD,JCOORD)
                pSendEcCell%buffMetrics(ibuff+14*nDim) = sFaceCorn(ZCOORD,JCOORD)

                pSendEcCell%buffMetrics(ibuff+15*nDim) = sFaceCorn(XCOORD,KCOORD)
                pSendEcCell%buffMetrics(ibuff+16*nDim) = sFaceCorn(YCOORD,KCOORD)
                pSendEcCell%buffMetrics(ibuff+17*nDim) = sFaceCorn(ZCOORD,KCOORD)

              ENDDO ! iFace

            ENDIF   ! iRegSrc 
          ENDDO     ! ijk

2999    CONTINUE
      ENDDO         ! iCorner

! Send buffers to destination processor ---------------------------------------

#ifdef MPI
      dest = regions(ir)%procid
      tag  = regions(ir)%localNumber+ PLAG_TAG_SHIFT +MPI_PATCHOFF +10
      IF(tag .gt. global%mpiTagMax) tag = MOD(tag,global%mpiTagMax)
      CALL MPI_Isend( pSendEcCell%buffMetrics,nBuffSize*nDim,MPI_RFREAL, &
                      dest,tag,global%mpiComm, &
                      pPlag%requestsMetrics(pSendEcCell%iRequestMetrics),&
                      global%mpierr )
      IF ( global%mpierr /= ERR_NONE ) &
        CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

    ENDIF      ! some cells to send
    ENDIF      ! not my processor
  ENDDO        ! ir

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLO_SendMetrics

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLO_SendMetrics.F90,v $
! Revision 1.5  2009/03/02 00:19:36  mtcampbe
! Added some ifdefs around Rocflo to disable particle injection on INFLOW
! boundaries and added some checks around MPI tags utilizing a new global
! data item, global%mpiTagMax.
!
! Revision 1.4  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 20:58:12  fnajjar
! Initial revision after changing case
!
! Revision 1.8  2004/11/29 19:22:39  fnajjar
! Added bypass statement for dengerate cells
!
! Revision 1.7  2004/03/21 00:43:32  fnajjar
! Fixed tags to be smaller number since Frost run-time system complains about size
!
! Revision 1.6  2004/03/06 21:25:05  fnajjar
! Added PLAG_TAG_SHIFT to MPI-based communication tags
!
! Revision 1.5  2004/02/13 01:24:49  fnajjar
! Included missing comma in ModInterfaces calling
!
! Revision 1.4  2004/02/11 23:12:35  fnajjar
! Included RFLO_GetNodeOffset in ModInterfaces call
!
! Revision 1.3  2004/02/10 21:46:37  fnajjar
! Added capability to remove coordinate alignment between corner-edge regions
!
! Revision 1.2  2004/01/28 21:22:07  fnajjar
! Moved statements inside iRegSrc IF loop to fix null state of iRegSrc
!
! Revision 1.1  2004/01/15 21:16:48  fnajjar
! Initial import for corner-edge cell metrics
!
!******************************************************************************







