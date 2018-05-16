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
! Purpose: copy values to edge and corner cells of an adjacent region
!          and append Lagrangian particle datastructure from buffers
!          communicated via the corner and edge cell infrastructure.
!
! Description: this is for the case if the other region is located
!              on the same processor.
!
! Input: regions = data of all regions,
!        iReg    = current region.
!
! Output: regionDes%level%plag%aiv,arv,cv,dv,tv = Lagrangian particles data.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_CECellsExchange.F90,v 1.4 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CECellsExchange( regions,iReg )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag, t_plag_input, t_buffer_plag
  USE ModDataStruct, ONLY : t_region, t_level, t_dCell
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetCellOffset,         &
                            RFLO_GetCornerCellsIndices, &
                            RFLO_GetEdgeCellsIndices
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER, INTENT(IN)     :: iReg

! ... loop variables
  INTEGER :: i, j, k, ijk, iEdge, iCorner, iPcls

! ... local variables
  INTEGER :: iBuff, icell, iRegDes, iLev, iCOff, ijCOff
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: nCornBuffSize, nEdgeBuffSize
  INTEGER :: nCorners, nEdges
  INTEGER :: nPclsDes, nPclsEnd, nPclsPrev, nPclsSrc, nPclsStart
  INTEGER :: iCornerBuffLoaded, iEdgeBuffLoaded
  INTEGER :: errorFlag,iCCMax, iECMax

  INTEGER, POINTER, DIMENSION(:,:) :: pAivBuff, pAivOldBuff
  INTEGER, POINTER, DIMENSION(:,:) :: pAivDes, pAivOldDes
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: cornCellCounter, edgeCellCounter

  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArvBuff, pArvOldBuff,         &
                                           pCvBuff, pCvOldBuff, pDvBuff,  &
                                           pTvBuff, pRhsBuff, pRhsSumBuff

  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArvDes, pArvOldDes,           &
                                           pCvDes, pCvOldDes, pDvDes,     &
                                           pTvDes, pRhsDes, pRhsSumDes

  TYPE(t_region),      POINTER :: pRegion, pRegionDes
  TYPE(t_level),       POINTER :: pLevel
  TYPE(t_buffer_plag), POINTER :: pCornCellsXBuff, pEdgeCellsXBuff
  TYPE(t_plag),        POINTER :: pPlagDes
  TYPE(t_global),      POINTER :: global

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PLAG_CECellsExchange',&
  'PLAG_CECellsExchange.F90' )

! Get dimensions --------------------------------------------------------------

  iLev = regions(iReg)%currLevel

  CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )

! Initialize values -----------------------------------------------------------

  nCorners = 8
  nEdges   = 12

  nPclsStart = 0
  nPclsEnd   = 0
  nPclsPrev  = 0
  nPclsDes   = 0

  nCornBuffSize = 0
  nEdgeBuffSize = 0

  iCornerBuffLoaded = 0
  iEdgeBuffLoaded   = 0

! Set pointers ----------------------------------------------------------------

  pRegion => regions(iReg)
  pLevel  => regions(iReg)%levels(iLev)

  nPclsSrc = pLevel%plag%nPcls

! Allocate corner cell buffer counter -----------------------------------------

  iCCMax = 0
  DO iCorner=1,nCorners
    IF( pLevel%cornerCells(iCorner)%interact ) &
      iCCMax = MAX(iCCMax,UBOUND(pLevel%cornerCells(iCorner)%cells,1))
  ENDDO ! iCorner

  ALLOCATE( cornCellCounter(nCorners,iCCMax),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= ERR_NONE) &
    CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  cornCellCounter = 0

#ifdef PLAG_CECELLS_DEBUG
    PRINT*,'PLAG_CECellsExchange: iReg,iCCMax = ', iReg, iCCMax
#endif

! Allocate edge cell buffer counter -------------------------------------------

  iECMax = 0
  DO iEdge=1,nEdges
    IF( pLevel%edgeCells(iEdge)%interact ) &
      iECMax = MAX(iECMax,UBOUND(pLevel%edgeCells(iEdge)%cells,1))
  ENDDO ! iEdge

  ALLOCATE( edgeCellCounter(nEdges,iECMax),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= ERR_NONE) &
    CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  edgeCellCounter = 0

#ifdef PLAG_CECELLS_DEBUG
    PRINT*,'PLAG_CECellsExchange: iReg,iECMax = ', iReg, iECMax
#endif

! Corner cells ----------------------------------------------------------------

  DO iCorner=1,nCorners

! - Bypass for noninteracting regions -----------------------------------------

    IF( .NOT. pLevel%cornerCells(iCorner)%interact ) GOTO 1999

! - Bypass for degenerate corner cells ----------------------------------------

      IF( pLevel%cornerCells(iCorner)%degenrt /= DEGENERAT_NONE ) GOTO 1999

! - Get corner indices --------------------------------------------------------

    CALL RFLO_GetCornerCellsIndices( pRegion,iLev,iCorner, &
                                     ibeg,iend,jbeg,jend,kbeg,kend )

    ijk = 0
    iCornerBuffLoaded = 0
    DO k=kbeg,kend
    DO j=jbeg,jend
    DO i=ibeg,iend
      ijk     =  ijk + 1

! -- Set pointers -------------------------------------------------------------

      pCornCellsXBuff => pLevel%cornerCells(iCorner)%cells(ijk)%bufferExchPlag

      pAivBuff    => pCornCellsXBuff%aiv
      pArvBuff    => pCornCellsXBuff%arv
      pCvBuff     => pCornCellsXBuff%cv
      pDvBuff     => pCornCellsXBuff%dv
      pTvBuff     => pCornCellsXBuff%tv
      pRhsBuff    => pCornCellsXBuff%rhs
      pRhsSumBuff => pCornCellsXBuff%rhsSum

      pAivOldBuff => pCornCellsXBuff%aivOld
      pArvOldBuff => pCornCellsXBuff%arvOld
      pCvOldBuff  => pCornCellsXBuff%cvOld

! -- Destination region infrastructure ----------------------------------------

      icell   =  pLevel%cornerCells(iCorner)%cells(ijk)%srcCell
      iRegDes =  pLevel%cornerCells(iCorner)%cells(ijk)%srcRegion

      IF ( iRegDes > 0 ) THEN
        pRegionDes => regions(iRegDes)
        pPlagDes   => pRegionDes%levels(iLev)%plag

! -- Get buffer size and start appending for non-null size --------------------

        nCornBuffSize = pCornCellsXBuff%nBuffSize
        IF ( cornCellCounter(iCorner,ijk) >= nCornBuffSize ) CYCLE

        nPclsStart = 0
        nPclsEnd   = 0
        nPclsPrev  = 0
        nPclsDes   = 0

        IF ( pRegionDes%procid == global%myProcid .AND. &
             nCornBuffSize > 0                          ) THEN

#ifdef PLAG_CECELLS_DEBUG
  WRITE(STDOUT,'(A,1PE12.5,3(2X,I5),6(3X,I4))') &
  '     PLAG_CECellsExchange: time, iReg, iRegDes, iCorner, nCornBuffSize,iCornerBuffLoaded, ijk, i ,j, k  = ',&
    global%currentTime+global%dtMin,iReg, iRegDes,iCorner, nCornBuffSize,iCornerBuffLoaded, ijk, i ,j, k
#endif

! -- Set pointers on destination region ---------------------------------------

          pAivDes    => pPlagDes%aiv

          pArvDes    => pPlagDes%arv
          pCvDes     => pPlagDes%cv
          pDvDes     => pPlagDes%dv
          pTvDes     => pPlagDes%tv
          pRhsDes    => pPlagDes%rhs
          pRhsSumDes => pPlagDes%rhsSum

          pAivOldDes => pPlagDes%aivOld
          pArvOldDes => pPlagDes%arvOld
          pCvOldDes  => pPlagDes%cvOld

! -- Set loop extent-----------------------------------------------------------

          nPclsDes   = pPlagDes%nPcls
          nPclsPrev  = nPclsDes

          nPclsStart = nPclsDes+1
          nPclsEnd   = nPclsStart + (nCornBuffSize-1)

! -- Append to PLAG datastructure with buffer arrays --------------------------

          DO iPcls = nPclsStart,nPclsEnd
            iBuff = iPcls-nPclsStart+1
            cornCellCounter(iCorner,ijk) = cornCellCounter(iCorner,ijk)+1
            iCornerBuffLoaded = cornCellCounter(iCorner,ijk)

            pAivDes(:,iPcls) = pAivBuff(:,iBuff)
            pArvDes(:,iPcls) = pArvBuff(:,iBuff)
            pCvDes(:,iPcls)  = pCvBuff( :,iBuff)
            pDvDes(:,iPcls)  = pDvBuff( :,iBuff)
            pTvDes(:,iPcls)  = pTvBuff( :,iBuff)
            pRhsDes(:,iPcls) = pRhsBuff(:,iBuff)
            pRhsSumDes(:,iPcls) = pRhsSumBuff(:,iBuff)

            pAivOldDes(:,iPcls) = pAivOldBuff(:,iBuff)
            pArvOldDes(:,iPcls) = pArvOldBuff(:,iBuff)
            pCvOldDes(:,iPcls)  = pCvOldBuff( :,iBuff)
          END DO ! iPcls

! -- Get new particle datasize ------------------------------------------------

          nPclsDes = nPclsDes+nCornBuffSize
          pPlagDes%nPcls = nPclsDes

#ifdef PLAG_CECELLS_DEBUG
    WRITE(STDOUT,'(A,A,2X,1PE15.7,2X,3(I3,2X),3(I4,3X))') &
  '     PLAG_CECellsExchange: time, iReg, iRegDes, iCorner, nPclsDes, nPclsSrc,iCornerBuffLoaded, ',&
  ' pAivDes(PIdini,Regini,RegC,ICells,IndexIJK,)', &
  global%currentTime+global%dtMin,iReg, iRegDes, iCorner, pPlagDes%nPcls, nPclsSrc,iCornerBuffLoaded
        DO iPcls = 1, pPlagDes%nPcls
         WRITE(STDOUT,'(9(I4,3X),8(1PE12.5,3X))') &
         iPcls, &
         pAivDes(AIV_PLAG_PIDINI,iPcls),&
         pAivDes(AIV_PLAG_REGINI,iPcls),&
         pAivDes(AIV_PLAG_REGCRT,iPcls),&
         pAivDes(AIV_PLAG_ICELLS,iPcls),&
         pAivDes(AIV_PLAG_INDEXI,iPcls),&
         pAivDes(AIV_PLAG_INDEXJ,iPcls),&
         pAivDes(AIV_PLAG_INDEXK,iPcls),&
         pAivDes(AIV_PLAG_BURNSTAT,iPcls),&
         pCvDes(CV_PLAG_XPOS:CV_PLAG_ZPOS,iPcls),&
         pCvDes(CV_PLAG_XMOM:CV_PLAG_ENER,iPcls),&
         pCvDes(CV_PLAG_ENERVAPOR,iPcls)
        ENDDO ! iPcls
#endif

        ENDIF ! procid
      ENDIF   ! iRegDes

    ENDDO     ! i
    ENDDO     ! j
    ENDDO     ! k

#ifdef PLAG_CECELLS_DEBUG
 IF ( nCornBuffSize > 0 ) THEN

    WRITE(STDOUT,'(A,A,2X,1PE15.7,2X,3(I3,2X),5(I4,3X))') &
  '     PLAG_CECellsExchange-iReg: iReg, iCorner, iRegDes, nCornBuffSize, ',&
  ' nPclsSrc,nPclsDesPrev,nPclsDes,nPclsStart,nPclsEnd = ',&
    global%currentTime+global%dtMin,iReg,iCorner,iRegDes, &
    nCornBuffSize,nPclsSrc,nPclsPrev,nPclsDes,nPclsStart,nPclsEnd

 ENDIF ! nCornBuffSize

#endif

! -- Reset buffer size to null ------------------------------------------------
! --   This insures no data clobbering occur for null particle size -----------

    IF ( pRegionDes%procid == global%myProcid ) THEN
      DO ijk = 1, UBOUND(pLevel%cornerCells(iCorner)%cells,1)
        pLevel%cornerCells(iCorner)%cells(ijk)%bufferExchPlag%nBuffSize = 0
      ENDDO ! ijk
    ENDIF ! pRegionDes%procid

1999 CONTINUE

  ENDDO         ! iCorner

! Edge cells ------------------------------------------------------------------

  DO iEdge=1,nEdges

! - Bypass for noninteracting regions -----------------------------------------

    IF( .NOT. pLevel%edgeCells(iEdge)%interact ) GOTO 2999

! - Bypass for degenerate edge cells ------------------------------------------

      IF( pLevel%edgeCells(iEdge)%degenrt /= DEGENERAT_NONE ) GOTO 2999

! - Get edge indices ----------------------------------------------------------

    CALL RFLO_GetEdgeCellsIndices( pRegion,iLev,iEdge, &
                                   ibeg,iend,jbeg,jend,kbeg,kend )

    ijk = 0
!    iEdgeBuffLoaded = 0
    DO k=kbeg,kend
    DO j=jbeg,jend
    DO i=ibeg,iend
      ijk     =  ijk + 1

! -- Set pointers -------------------------------------------------------------

      pEdgeCellsXBuff => pLevel%edgeCells(iEdge)%cells(ijk)%bufferExchPlag

      pAivBuff    => pEdgeCellsXBuff%aiv
      pArvBuff    => pEdgeCellsXBuff%arv
      pCvBuff     => pEdgeCellsXBuff%cv
      pDvBuff     => pEdgeCellsXBuff%dv
      pTvBuff     => pEdgeCellsXBuff%tv
      pRhsBuff    => pEdgeCellsXBuff%rhs
      pRhsSumBuff => pEdgeCellsXBuff%rhsSum

      pAivOldBuff => pEdgeCellsXBuff%aivOld
      pArvOldBuff => pEdgeCellsXBuff%arvOld
      pCvOldBuff  => pEdgeCellsXBuff%cvOld

! -- Destination region infrastructure ----------------------------------------

      icell   =  pLevel%edgeCells(iEdge)%cells(ijk)%srcCell
      iRegDes =  pLevel%edgeCells(iEdge)%cells(ijk)%srcRegion

      IF ( iRegDes > 0 ) THEN
        pRegionDes => regions(iRegDes)
        pPlagDes   => pRegionDes%levels(iLev)%plag

! -- Get buffer size and start appending for non-null size --------------------

        nEdgeBuffSize = pEdgeCellsXBuff%nBuffSize
        IF ( edgeCellCounter(iEdge,ijk) >= nEdgeBuffSize ) CYCLE

        nPclsStart = 0
        nPclsEnd   = 0
        nPclsPrev  = 0
        nPclsDes   = 0

        IF ( pRegionDes%procid == global%myProcid .AND. &
             nEdgeBuffSize > 0                          ) THEN

#ifdef PLAG_CECELLS_DEBUG
  WRITE(STDOUT,'(A,1PE12.5,3(2X,I3),5(3X,I4))') &
  '     PLAG_CECellsExchange: time, iReg, iRegDes, iEdge, nEdgeBuffSize,iEdgeBuffLoaded,  ijk, i ,j, k  = ',&
    global%currentTime+global%dtMin,iReg, iRegDes,iEdge, nEdgeBuffSize, iEdgeBuffLoaded ,ijk, i ,j, k
#endif

! -- Set pointers on destination region ---------------------------------------

          pAivDes    => pPlagDes%aiv

          pArvDes    => pPlagDes%arv
          pCvDes     => pPlagDes%cv
          pDvDes     => pPlagDes%dv
          pTvDes     => pPlagDes%tv
          pRhsDes    => pPlagDes%rhs
          pRhsSumDes => pPlagDes%rhsSum

          pAivOldDes => pPlagDes%aivOld
          pArvOldDes => pPlagDes%arvOld
          pCvOldDes  => pPlagDes%cvOld

! -- Set loop extent-----------------------------------------------------------

          nPclsDes   = pPlagDes%nPcls
          nPclsPrev  = nPclsDes

          nPclsStart = nPclsDes+1
          nPclsEnd   = nPclsStart + (nEdgeBuffSize-1)

! -- Append to PLAG datastructure with buffer arrays --------------------------

          DO iPcls = nPclsStart,nPclsEnd
            iBuff = iPcls-nPclsStart+1
            edgeCellCounter(iEdge,ijk) = edgeCellCounter(iEdge,ijk)+1
            iEdgeBuffLoaded = edgeCellCounter(iEdge,ijk)

            pAivDes(:,iPcls) = pAivBuff(:,iBuff)
            pArvDes(:,iPcls) = pArvBuff(:,iBuff)
            pCvDes(:,iPcls)  = pCvBuff( :,iBuff)
            pDvDes(:,iPcls)  = pDvBuff( :,iBuff)
            pTvDes(:,iPcls)  = pTvBuff( :,iBuff)
            pRhsDes(:,iPcls) = pRhsBuff(:,iBuff)
            pRhsSumDes(:,iPcls) = pRhsSumBuff(:,iBuff)

            pAivOldDes(:,iPcls) = pAivOldBuff(:,iBuff)
            pArvOldDes(:,iPcls) = pArvOldBuff(:,iBuff)
            pCvOldDes(:,iPcls)  = pCvOldBuff( :,iBuff)
          END DO ! iPcls

! -- Get new particle datasize ------------------------------------------------

          nPclsDes = nPclsDes+nEdgeBuffSize
          pPlagDes%nPcls = nPclsDes

#ifdef PLAG_CECELLS_DEBUG
    WRITE(STDOUT,'(A,A,2X,1PE15.7,2X,3(I3,2X),3(I4,3X))') &
  '     PLAG_CECellsExchange: time, iReg, iRegDes, iEdge, nPclsDes, nPclsSrc,iEdgeBuffLoaded, ',&
  ' pAivDes(PIdini,Regini,RegC,ICells,IndexIJK,)', &
  global%currentTime+global%dtMin,iReg, iRegDes, iEdge, pPlagDes%nPcls, nPclsSrc,iEdgeBuffLoaded
        DO iPcls = 1, pPlagDes%nPcls
         WRITE(STDOUT,'(9(I4,3X),8(1PE12.5,3X))') &
         iPcls, &
         pAivDes(AIV_PLAG_PIDINI,iPcls),&
         pAivDes(AIV_PLAG_REGINI,iPcls),&
         pAivDes(AIV_PLAG_REGCRT,iPcls),&
         pAivDes(AIV_PLAG_ICELLS,iPcls),&
         pAivDes(AIV_PLAG_INDEXI,iPcls),&
         pAivDes(AIV_PLAG_INDEXJ,iPcls),&
         pAivDes(AIV_PLAG_INDEXK,iPcls),&
         pAivDes(AIV_PLAG_BURNSTAT,iPcls),&
         pCvDes(CV_PLAG_XPOS:CV_PLAG_ZPOS,iPcls),&
         pCvDes(CV_PLAG_XMOM:CV_PLAG_ENER,iPcls),&
         pCvDes(CV_PLAG_ENERVAPOR,iPcls)
        ENDDO ! iPcls
#endif

        ENDIF ! procid
      ENDIF   ! iRegDes

    ENDDO     ! i
    ENDDO     ! j
    ENDDO     ! k

#ifdef PLAG_CECELLS_DEBUG
 IF ( nEdgeBuffSize > 0 ) THEN

    WRITE(STDOUT,'(A,A,2X,1PE15.7,2X,3(I3,2X),5(I4,3X))') &
  '     PLAG_CECellsExchange-iReg: iReg, iEdge, iRegDes, nEdgeBuffSize, ',&
  ' nPclsSrc,nPclsDesPrev,nPclsDes,nPclsStart,nPclsEnd = ',&
    global%currentTime+global%dtMin,iReg,iEdge,iRegDes, &
    nEdgeBuffSize,nPclsSrc,nPclsPrev,nPclsDes,nPclsStart,nPclsEnd

 ENDIF ! nEdgeBuffSize

#endif

! -- Reset buffer size to null ------------------------------------------------
! --   This insures no data clobbering occur for null particle size -----------

    IF ( pRegionDes%procid == global%myProcid ) THEN
      DO ijk = 1, UBOUND(pLevel%edgeCells(iEdge)%cells,1)
        pLevel%edgeCells(iEdge)%cells(ijk)%bufferExchPlag%nBuffSize = 0
      ENDDO ! ijk
    ENDIF ! pRegionDes%procid

2999 CONTINUE

  ENDDO         ! iEdge

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CECellsExchange

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CECellsExchange.F90,v $
! Revision 1.4  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 20:57:09  fnajjar
! Initial revision after changing case
!
! Revision 1.10  2004/11/29 19:24:34  fnajjar
! Added bypass statement for dengerate cells
!
! Revision 1.9  2004/03/25 21:16:43  jferry
! fixed Vapor Energy bug
!
! Revision 1.8  2004/03/20 21:59:58  fnajjar
! Included reset buffer size to zero in IF statement of on-processor communication
!
! Revision 1.7  2004/03/20 21:51:51  fnajjar
! Reset buffer size to null to insure no data clobbering
!
! Revision 1.6  2004/03/20 00:17:58  fnajjar
! Moved and added IFDEF with WRITE statement for corner cell section
!
! Revision 1.5  2004/03/19 23:49:50  fnajjar
! Fixed WRITE formating
!
! Revision 1.4  2004/03/18 21:42:14  fnajjar
! Various bug fixed for proper buffer loading
!
! Revision 1.3  2004/01/26 22:54:01  fnajjar
! Renamed ifdef PLAG_DEBUG to PLAG_CECELLS_DEBUG
!
! Revision 1.2  2003/11/21 22:41:18  fnajjar
! Activated IFDEF PLAG_DEBUG
!
! Revision 1.1  2003/11/12 21:37:59  fnajjar
! Initial import of Corner-Edge cells Infrastructure
!
!******************************************************************************







