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
! Purpose: copy face vectors in corner and edge cells from an adjacent 
!          region on the same processor
!
! Description: none.
!
! Input: regions = data of all regions,
!        iReg    = current region number.
!
! Output: plag%si,sj,sk = face vectors.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_CECellsFaceVectors.F90,v 1.4 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CECellsFaceVectors( regions, iReg )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region, t_level, t_dCell
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag 
  USE ModIndexing, ONLY   : GetIJK
  USE ModInterfaces, ONLY : RFLO_GetCellOffset,         &
                            RFLO_GetCornerCellsIndices, &
                            RFLO_GetEdgeCellsIndices,   &
                            RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  USE PLAG_ModInterfaces, ONLY : PLAG_RFLO_GetFaceMapping
  IMPLICIT NONE
  
#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER, INTENT(IN)     :: iReg

! ... loop variables
  INTEGER :: i, iCorner, iEdge, iFace, ijk, j, k, ldir

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iLev, nCorners, nEdges, nFaces, nDir
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff
  INTEGER :: iCOffSrc, ijCOffSrc, iRegSrc
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: ijkC, ijkCCSrc, ijkECSrc

  INTEGER :: iCCSrc, jCCSrc, kCCSrc
  INTEGER :: iECSrc, jECSrc, kECSrc
  INTEGER :: nDumCellsSrc

  INTEGER :: iNOffSrc, ijNOffSrc, ijkN, ijkCNSrc, ijkENSrc
  INTEGER :: iCNSrc, jCNSrc, kCNSrc
  INTEGER :: iENSrc, jENSrc, kENSrc
  INTEGER :: nodeCorn, nodeCornSrc, nodeEdge, nodeEdgeSrc 
  INTEGER :: idirSrc, jdirSrc, kdirSrc, ldirSrc, iFaceSrc
  INTEGER :: srcDir(3),srcFace(6)
  INTEGER :: srcIndexMapMat(3,4)
  
  REAL(RFREAL), DIMENSION(3,3)          :: sFaceCorn, sFaceEdge
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pSi, pSj, pSk, &
                                           pSiSrc, pSjSrc, pSkSrc

  TYPE(t_region),      POINTER :: pRegion
  TYPE(t_level),       POINTER :: pLevel
  TYPE(t_global),      POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PLAG_CECellsFaceVectors.F90,v $ $Revision: 1.4 $'

  global => regions(iReg)%global
    
  CALL RegisterFunction( global, 'PLAG_CECellsFaceVectors',&
  'PLAG_CECellsFaceVectors.F90' )

! Get dimensions --------------------------------------------------------------

  iLev  = regions(iReg)%currLevel
  nCorners = 8
  nEdges   = 12
  nFaces   = 6
  nDir     = 3

! Set pointers ----------------------------------------------------------------
  
  pRegion => regions(iReg)
  pLevel  => regions(iReg)%levels(iLev) 

  pSi     => regions(iReg)%levels(iLev)%plag%si
  pSj     => regions(iReg)%levels(iLev)%plag%sj
  pSk     => regions(iReg)%levels(iLev)%plag%sk

! Get node offset ------------------------------------------------------------

  CALL RFLO_GetNodeOffset( pRegion,iLev,iNOff,ijNOff )

! Loop over edges -------------------------------------------------------------

  DO iEdge=1,nEdges 

! - Bypass for noninteracting regions -----------------------------------------

    IF( .NOT. pLevel%edgeCells(iEdge)%interact ) GOTO 1999

! -- Bypass for degenerate edge cells -----------------------------------------

    IF( pLevel%edgeCells(iEdge)%degenrt /= DEGENERAT_NONE ) GOTO 1999

! - Loop over edge cell indices -----------------------------------------------
        
    CALL RFLO_GetEdgeCellsIndices( pRegion,iLev,iEdge, &
                                   ibeg,iend,jbeg,jend,kbeg,kend )

    ijk = 0 
    DO k=kbeg,kend
    DO j=jbeg,jend
    DO i=ibeg,iend
      ijk     =  ijk + 1
      ijkN = IndIJK(i,j,k,iNOff, ijNOff)

! -- Source region infrastructure ----------------------------------------
            
      iRegSrc      = pLevel%edgeCells(iEdge)%cells(ijk)%srcRegion            

      IF ( iRegSrc > 0 ) THEN
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

! -- Check if both regions are on the same processor --------------------------

        IF ( regions(iRegSrc)%procid == global%myProcid ) THEN

! --- Set pointers for edge cells ---------------------------------------------

          pSiSrc => regions(iRegSrc)%levels(iLev)%plag%si
          pSjSrc => regions(iRegSrc)%levels(iLev)%plag%sj
          pSkSrc => regions(iRegSrc)%levels(iLev)%plag%sk

! --- Loop over all the Faces -------------------------------------------------
! --- Note: assumption of coordinate alignment between regions removed --------

          DO iFace = 1, nFaces
            SELECT CASE (iFace)
              CASE(1,3,5)
                nodeEdge    = IndIJK(i,j,k,iNOff, ijNOff)

              CASE(2)
                nodeEdge    = IndIJK(i+1,j,k,iNOff, ijNOff)

              CASE(4)
                nodeEdge    = IndIJK(i,j+1,k,iNOff, ijNOff)

              CASE(6)
                nodeEdge    = IndIJK(i,j,k+1,iNOff, ijNOff)
            END SELECT ! iFace  

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
                  sFaceEdge(XCOORD:ZCOORD,ldirSrc) = pSiSrc(XCOORD:ZCOORD,nodeEdgeSrc)
                
                CASE(JCOORD)
                  sFaceEdge(XCOORD:ZCOORD,ldirSrc) = pSjSrc(XCOORD:ZCOORD,nodeEdgeSrc)
                
                CASE(KCOORD)
                  sFaceEdge(XCOORD:ZCOORD,ldirSrc) = pSkSrc(XCOORD:ZCOORD,nodeEdgeSrc)
              END SELECT ! ldirSrc
            ENDDO ! ldir

! --- Load face normal vectors ------------------------------------------------
            
            pSi(XCOORD:ZCOORD,nodeEdge) = sFaceEdge(XCOORD:ZCOORD,ICOORD)
            pSj(XCOORD:ZCOORD,nodeEdge) = sFaceEdge(XCOORD:ZCOORD,JCOORD)
            pSk(XCOORD:ZCOORD,nodeEdge) = sFaceEdge(XCOORD:ZCOORD,KCOORD)
          ENDDO ! iFace  

        END IF ! procid 
      END IF   ! iRegSrc 

    ENDDO   ! i
    ENDDO   ! j
    ENDDO   ! k

1999 CONTINUE 
      
  ENDDO      ! iEdge
      
! Loop over corners -----------------------------------------------------------

  DO iCorner=1,nCorners

! - Bypass for noninteracting regions -----------------------------------------

    IF( .NOT. pLevel%cornerCells(iCorner)%interact ) GOTO 2999

! -- Bypass for degenerate corner cells ---------------------------------------

    IF( pLevel%cornerCells(iCorner)%degenrt /= DEGENERAT_NONE ) GOTO 2999

! - Loop over corner cell indices ---------------------------------------------
        
    CALL RFLO_GetCornerCellsIndices( pRegion,iLev,iCorner, &
                                     ibeg,iend,jbeg,jend,kbeg,kend )
        
    ijk = 0 
    DO k=kbeg,kend
    DO j=jbeg,jend
    DO i=ibeg,iend
      ijk  =  ijk + 1
      ijkN = IndIJK(i,j,k,iNOff, ijNOff)
  
! - Source region infrastructure -----------------------------------------
            
      iRegSrc  =  pLevel%cornerCells(iCorner)%cells(ijk)%srcRegion 

      IF ( iRegSrc > 0 ) THEN
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

! -- Check if both regions are on the same processor --------------------------

        IF ( regions(iRegSrc)%procid == global%myProcid ) THEN

! --- Set pointers for corner cells -------------------------------------------

          pSiSrc => regions(iRegSrc)%levels(iLev)%plag%si
          pSjSrc => regions(iRegSrc)%levels(iLev)%plag%sj
          pSkSrc => regions(iRegSrc)%levels(iLev)%plag%sk

! --- Loop over all the Faces -------------------------------------------------
! --- Note: assumption of coordinate alignment between regions removed --------

          DO iFace = 1, nFaces
            SELECT CASE (iFace)
              CASE(1,3,5)
                nodeCorn    = IndIJK(i,j,k,iNOff, ijNOff)

              CASE(2)
                nodeCorn    = IndIJK(i+1,j,k,iNOff, ijNOff)

              CASE(4)
                nodeCorn    = IndIJK(i,j+1,k,iNOff, ijNOff)

              CASE(6)
                nodeCorn    = IndIJK(i,j,k+1,iNOff, ijNOff)
            END SELECT  ! iFace 

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
                  sFaceCorn(XCOORD:ZCOORD,ldirSrc) = pSiSrc(XCOORD:ZCOORD,nodeCornSrc)
                
                CASE(JCOORD)
                  sFaceCorn(XCOORD:ZCOORD,ldirSrc) = pSjSrc(XCOORD:ZCOORD,nodeCornSrc)
                
                CASE(KCOORD)
                  sFaceCorn(XCOORD:ZCOORD,ldirSrc) = pSkSrc(XCOORD:ZCOORD,nodeCornSrc)
              END SELECT ! ldirSrc
            ENDDO ! ldir

! --- Load face normal vectors ------------------------------------------------

            pSi(XCOORD:ZCOORD,nodeCorn) = sFaceCorn(XCOORD:ZCOORD,ICOORD)
            pSj(XCOORD:ZCOORD,nodeCorn) = sFaceCorn(XCOORD:ZCOORD,JCOORD)
            pSk(XCOORD:ZCOORD,nodeCorn) = sFaceCorn(XCOORD:ZCOORD,KCOORD)
          ENDDO ! iFace 

        END IF ! procid 
      END IF   ! iRegSrc 
 
    ENDDO   ! i
    ENDDO   ! j
    ENDDO   ! k

2999 CONTINUE

  ENDDO       ! iCorner

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CECellsFaceVectors

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CECellsFaceVectors.F90,v $
! Revision 1.4  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 20:57:12  fnajjar
! Initial revision after changing case
!
! Revision 1.7  2004/11/29 19:25:07  fnajjar
! Added bypass statement for dengerate cells
!
! Revision 1.6  2004/02/11 23:12:35  fnajjar
! Included RFLO_GetNodeOffset in ModInterfaces call
!
! Revision 1.5  2004/02/10 21:46:37  fnajjar
! Added capability to remove coordinate alignment between corner-edge regions
!
! Revision 1.4  2004/01/28 16:10:28  fnajjar
! Moved statements inside IF iReg for correct syntax
!
! Revision 1.3  2004/01/14 21:27:24  fnajjar
! Included definition of nFaces
!
! Revision 1.2  2003/12/05 23:45:30  fnajjar
! Fixed metrics for edge and corner cells based on face looping
!
! Revision 1.1  2003/11/12 21:37:59  fnajjar
! Initial import of Corner-Edge cells Infrastructure
!
!******************************************************************************







