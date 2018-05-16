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
! Purpose: receives RFLO metrics to edge and corner cells of an adjacent 
!          region.
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
! $Id: PLAG_RFLO_RecvMetrics.F90,v 1.3 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLO_RecvMetrics( regions,iReg )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region, t_level, t_dCellTransf
  USE ModInterfaces, ONLY : RFLO_GetCellOffset, RFLO_GetEdgeCellsIndices, &
                            RFLO_GetCornerCellsIndices
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER, INTENT(IN)     :: iReg

! ... loop variables
  INTEGER :: ir, iCorner, iEdge, iFace, i, j, k, ijk

! ... local variables
#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: ibuff, icell, iLev, iRegSrc, nDim, source, tag
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iCOff, ijCOff
  INTEGER :: nBuffSize, nCorners, nEdges
  INTEGER :: nFaces, nFaceCentroidSize, nFaceNormalSize
  INTEGER :: iNOff, ijNOff, ijkN
  INTEGER :: nDumCellsSrc
  INTEGER :: nodeCorn, nodeEdge

  REAL(RFREAL), POINTER, DIMENSION(:,:)   :: pSi, pSj, pSk
  REAL(RFREAL), POINTER, DIMENSION(:,:,:) :: pFc
 
  TYPE(t_region),      POINTER :: pRegion 
  TYPE(t_level),       POINTER :: pLevel
  TYPE(t_dCellTransf), POINTER :: pRecvEcCell
  TYPE(t_global),      POINTER :: global

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PLAG_RFLO_RecvMetrics',&
  'PLAG_RFLO_RecvMetrics.F90' )

! Get dimensions --------------------------------------------------------------

  iLev = regions(iReg)%currLevel
  nCorners  = 8
  nEdges    = 12
  nFaces    = 6
  
! Compute buffer size has to store information --------------------------------
!  on the total number of faces for each cell .
  
  nFaceCentroidSize = ZCOORD*KCOORD
  nFaceNormalSize   =      3*KCOORD
  nBuffSize = (nFaceCentroidSize + nFaceNormalSize)

! Set pointers ----------------------------------------------------------------
  
  pRegion => regions(iReg)
  pLevel  => regions(iReg)%levels(iLev) 
  pFc     => regions(iReg)%levels(iLev)%plag%fc
  pSi     => regions(iReg)%levels(iLev)%plag%si
  pSj     => regions(iReg)%levels(iLev)%plag%sj
  pSk     => regions(iReg)%levels(iLev)%plag%sk

! Get node offset -------------------------------------------------------------

  CALL RFLO_GetNodeOffset( pRegion,iLev,iNOff,ijNOff ) 

! copy data from buffer to dummy cells

  DO ir=1,global%nRegions
    IF (regions(ir)%procid /= global%myProcid) THEN
    IF (pLevel%recvEcCells(ir)%nCells > 0) THEN

      pRecvEcCell => pLevel%recvEcCells(ir)
      nDim        =  pRecvEcCell%nCells*nFaces
      ibuff       =  0

! - Receive buffers from source processor ---------------------------------------

#ifdef MPI
      source = regions(ir)%procid
      tag    = regions(iReg)%localNumber+ PLAG_TAG_SHIFT +MPI_PATCHOFF +10

      CALL MPI_Recv( pRecvEcCell%buffMetrics,nBuffSize*nDim,MPI_RFREAL, &
                     source,tag,global%mpiComm,status,global%mpierr )
      IF ( global%mpierr /= ERR_NONE ) &
        CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! -- Load edges ---------------------------------------------------------------

      DO iEdge=1,nEdges

! --- Bypass for noninteracting regions ---------------------------------------

        IF( .NOT. pLevel%edgeCells(iEdge)%interact ) GOTO 1999

! -- Bypass for degenerate edge cells -----------------------------------------

        IF( pLevel%edgeCells(iEdge)%degenrt /= DEGENERAT_NONE ) GOTO 1999

! --- Current region infrastructure -------------------------------------------

        CALL RFLO_GetEdgeCellsIndices( pRegion,iLev,iEdge, &
                                       ibeg,iend,jbeg,jend,kbeg,kend )

        ijk = 0
        DO k=kbeg,kend
        DO j=jbeg,jend
        DO i=ibeg,iend 
          ijk     =  ijk + 1
          ijkN    = IndIJK(i,j,k,iNOff, ijNOff)
          iRegSrc = pLevel%edgeCells(iEdge)%cells(ijk)%srcRegion
          IF ( iRegSrc == ir ) THEN

            IF ( pLevel%edgeCells(iEdge)%cells(ijk)%rotate ) THEN
                    ! rotational periodicity
            ELSE

! --- Loop over all the Faces -------------------------------------------------
! --- Note: assumes coordinate alignment between regions ----------------------

              DO iFace = 1, nFaces
                ibuff = ibuff + 1

                SELECT CASE (iFace)
                  CASE(1,3,5)
                    nodeEdge    = IndIJK(i,j,k,iNOff, ijNOff)

                  CASE(2)
                    nodeEdge    = IndIJK(i+1,j,k,iNOff, ijNOff)

                  CASE(4)
                    nodeEdge    = IndIJK(i,j+1,k,iNOff, ijNOff)

                  CASE(6)
                    nodeEdge    = IndIJK(i,j,k+1,iNOff, ijNOff)
                END SELECT  

! ---- Load face centroids -----------------------------------------------------

                pFc(XCOORD,ICOORD,nodeEdge) = pRecvEcCell%buffMetrics(ibuff        )
                pFc(XCOORD,JCOORD,nodeEdge) = pRecvEcCell%buffMetrics(ibuff+   nDim)
                pFc(XCOORD,ZCOORD,nodeEdge) = pRecvEcCell%buffMetrics(ibuff+ 2*nDim)

                pFc(YCOORD,ICOORD,nodeEdge) = pRecvEcCell%buffMetrics(ibuff+ 3*nDim)
                pFc(YCOORD,JCOORD,nodeEdge) = pRecvEcCell%buffMetrics(ibuff+ 4*nDim)
                pFc(YCOORD,ZCOORD,nodeEdge) = pRecvEcCell%buffMetrics(ibuff+ 5*nDim)
                
                pFc(ZCOORD,ICOORD,nodeEdge) = pRecvEcCell%buffMetrics(ibuff+ 6*nDim)
                pFc(ZCOORD,JCOORD,nodeEdge) = pRecvEcCell%buffMetrics(ibuff+ 7*nDim)
                pFc(ZCOORD,ZCOORD,nodeEdge) = pRecvEcCell%buffMetrics(ibuff+ 8*nDim)

! --- Load face normals -------------------------------------------------------

                pSi(XCOORD,nodeEdge)        = pRecvEcCell%buffMetrics(ibuff+ 9*nDim)
                pSi(YCOORD,nodeEdge)        = pRecvEcCell%buffMetrics(ibuff+10*nDim)
                pSi(ZCOORD,nodeEdge)        = pRecvEcCell%buffMetrics(ibuff+11*nDim)

                pSj(XCOORD,nodeEdge)        = pRecvEcCell%buffMetrics(ibuff+12*nDim)
                pSj(YCOORD,nodeEdge)        = pRecvEcCell%buffMetrics(ibuff+13*nDim)
                pSj(ZCOORD,nodeEdge)        = pRecvEcCell%buffMetrics(ibuff+14*nDim)

                pSk(XCOORD,nodeEdge)        = pRecvEcCell%buffMetrics(ibuff+15*nDim)
                pSk(YCOORD,nodeEdge)        = pRecvEcCell%buffMetrics(ibuff+16*nDim)
                pSk(ZCOORD,nodeEdge)        = pRecvEcCell%buffMetrics(ibuff+17*nDim)
              ENDDO ! iFace 

             ENDIF  ! rotate
            ENDIF   !iRegSrc

          ENDDO     ! i
          ENDDO     ! j
          ENDDO     ! k

1999    CONTINUE
      ENDDO    ! iEdge

! --- Load corners ------------------------------------------------------------

      DO iCorner=1,nCorners

! -- Bypass for noninteracting regions ----------------------------------------

        IF (.NOT. pLevel%cornerCells(iCorner)%interact) GOTO 2999

! -- Bypass for degenerate corner cells ---------------------------------------

        IF( pLevel%cornerCells(iCorner)%degenrt /= DEGENERAT_NONE ) GOTO 2999

! -- Current region infrastructure ---------------------------------------------

        CALL RFLO_GetCornerCellsIndices( pRegion,iLev,iCorner, &
                                         ibeg,iend,jbeg,jend,kbeg,kend )

        ijk = 0
        DO k=kbeg,kend
        DO j=jbeg,jend
        DO i=ibeg,iend
          ijk     = ijk + 1
          ijkN    = IndIJK(i,j,k,iNOff,ijNOff)
          iRegSrc = pLevel%cornerCells(iCorner)%cells(ijk)%srcRegion
          IF ( iRegSrc == ir ) THEN

            IF ( pLevel%cornerCells(iCorner)%cells(ijk)%rotate ) THEN
                   ! rotational periodicity
            ELSE

! --- Loop over all the Faces -------------------------------------------------
! --- Note: assumes coordinate alignment between regions ----------------------

              DO iFace = 1, nFaces
                ibuff = ibuff + 1

                SELECT CASE (iFace)
                  CASE(1,3,5)
                    nodeCorn = IndIJK(i,j,k,iNOff, ijNOff)

                  CASE(2)
                    nodeCorn = IndIJK(i+1,j,k,iNOff, ijNOff)

                  CASE(4)
                    nodeCorn = IndIJK(i,j+1,k,iNOff, ijNOff)

                  CASE(6)
                    nodeCorn = IndIJK(i,j,k+1,iNOff, ijNOff)
                END SELECT  

! ---- Load face centroids -----------------------------------------------------

                pFc(XCOORD,ICOORD,nodeCorn) = pRecvEcCell%buffMetrics(ibuff        )
                pFc(XCOORD,JCOORD,nodeCorn) = pRecvEcCell%buffMetrics(ibuff+   nDim)
                pFc(XCOORD,ZCOORD,nodeCorn) = pRecvEcCell%buffMetrics(ibuff+ 2*nDim)

                pFc(YCOORD,ICOORD,nodeCorn) = pRecvEcCell%buffMetrics(ibuff+ 3*nDim)
                pFc(YCOORD,JCOORD,nodeCorn) = pRecvEcCell%buffMetrics(ibuff+ 4*nDim)
                pFc(YCOORD,ZCOORD,nodeCorn) = pRecvEcCell%buffMetrics(ibuff+ 5*nDim)
                
                pFc(ZCOORD,ICOORD,nodeCorn) = pRecvEcCell%buffMetrics(ibuff+ 6*nDim)
                pFc(ZCOORD,JCOORD,nodeCorn) = pRecvEcCell%buffMetrics(ibuff+ 7*nDim)
                pFc(ZCOORD,ZCOORD,nodeCorn) = pRecvEcCell%buffMetrics(ibuff+ 8*nDim)

! --- Load face normals -------------------------------------------------------

                pSi(XCOORD,nodeCorn)        = pRecvEcCell%buffMetrics(ibuff+ 9*nDim)
                pSi(YCOORD,nodeCorn)        = pRecvEcCell%buffMetrics(ibuff+10*nDim)
                pSi(ZCOORD,nodeCorn)        = pRecvEcCell%buffMetrics(ibuff+11*nDim)
                pSj(XCOORD,nodeCorn)        = pRecvEcCell%buffMetrics(ibuff+12*nDim)
                pSj(YCOORD,nodeCorn)        = pRecvEcCell%buffMetrics(ibuff+13*nDim)
                pSj(ZCOORD,nodeCorn)        = pRecvEcCell%buffMetrics(ibuff+14*nDim)
                pSk(XCOORD,nodeCorn)        = pRecvEcCell%buffMetrics(ibuff+15*nDim)
                pSk(YCOORD,nodeCorn)        = pRecvEcCell%buffMetrics(ibuff+16*nDim)
                pSk(ZCOORD,nodeCorn)        = pRecvEcCell%buffMetrics(ibuff+17*nDim)
              ENDDO ! iFace 

             ENDIF  ! rotate
            ENDIF   !iRegSrc

          ENDDO     ! i
          ENDDO     ! j
          ENDDO     ! k

2999    CONTINUE 
      ENDDO    ! icorner

    ENDIF      ! some cells to receive
    ENDIF      ! not my processor
  ENDDO        ! ir

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLO_RecvMetrics

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLO_RecvMetrics.F90,v $
! Revision 1.3  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:58:11  fnajjar
! Initial revision after changing case
!
! Revision 1.4  2004/11/29 19:22:39  fnajjar
! Added bypass statement for dengerate cells
!
! Revision 1.3  2004/03/21 00:43:32  fnajjar
! Fixed tags to be smaller number since Frost run-time system complains about size
!
! Revision 1.2  2004/03/06 21:25:05  fnajjar
! Added PLAG_TAG_SHIFT to MPI-based communication tags
!
! Revision 1.1  2004/01/15 21:16:48  fnajjar
! Initial import for corner-edge cell metrics
!
!******************************************************************************







