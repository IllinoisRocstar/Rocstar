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
! Purpose: obtains appropriate buffer sizes for corner and edge cells. 
!
! Description: none.
!
! Input: region = current region.
!        iReg   = region number
!
! Output: region%level%cornerCells(:)%bufferExchPlag%nBuffSize = buffer size.
!         region%level%edgeCells(:)%bufferExchPlag%nBuffSize   = buffer size.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_CECellsGetBufferSize.F90,v 1.3 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CECellsGetBufferSize( region, iReg )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag, t_plag_input, t_buffer_plag
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global  
  USE ModInterfaces, ONLY : RFLO_GetCellOffset, RFLO_GetCornerCellsIndices, &
                            RFLO_GetEdgeCellsIndices, RFLO_GetDimensDummy,  &
                            RFLO_GetDimensPhys

  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region)      :: region

  INTEGER, INTENT(IN) :: iReg

! ... loop variables
  INTEGER :: i, j, k, ijk, iCorner, iEdge, iPcls 

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iCPlag, jCPlag, kCPlag
  INTEGER :: iCOff, ijCOff, ijkCCell, ijkCPlag, ijkECell, iLev
  INTEGER :: ibegDummCell, iendDummCell, ibegPhysCell, iendPhysCell
  INTEGER :: ibegCornCell, iendCornCell, ibegEdgeCell, iendEdgeCell
  INTEGER :: lboundSkipSum, nBuffSizeTot, nDumCells, nPatches, nPcls
  INTEGER :: nCorners, nEdges
   
  INTEGER,          DIMENSION(6)   :: lboundSkip
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv
  
  LOGICAL :: pclsFound

  TYPE(t_level),       POINTER :: pLevel
  TYPE(t_plag),        POINTER :: pPlag
  TYPE(t_buffer_plag), POINTER :: pCornCellsXBuff, pEdgeCellsXBuff
  TYPE(t_global),      POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CECellsGetBufferSize.F90,v $ $Revision: 1.3 $'

  global => region%global
    
  CALL RegisterFunction( global, 'PLAG_CECellsGetBufferSize',&
  'PLAG_CECellsGetBufferSize.F90' )

! Get dimensions --------------------------------------------------------------

  iLev         = region%currLevel
  nDumCells    = region%nDumCells
  nCorners     = 8
  nEdges       = 12
  
  nPcls        = region%levels(iLev)%plag%nPcls 
  nBuffSizeTot = region%plagInput%nPclsBuffTot

! Set pointers ----------------------------------------------------------------
  
  pLevel => region%levels(iLev)
  pPlag  => pLevel%plag   
  pAiv   => pPlag%aiv

! Initialize buffer size ----------------------------------------------------

  DO iCorner = 1, nCorners
    IF( pLevel%cornerCells(iCorner)%interact .EQV. .TRUE. .AND.    &
        pLevel%cornerCells(iCorner)%degenrt   ==   DEGENERAT_NONE ) THEN 
      DO ijk = 1, UBOUND(pLevel%cornerCells(iCorner)%cells,1)
        pLevel%cornerCells(iCorner)%cells(ijk)%bufferExchPlag%nBuffSize  = 0
      ENDDO ! ijk
    ENDIF   ! interact
  ENDDO     ! iCorner

  DO iEdge = 1, nEdges
    IF( pLevel%edgeCells(iEdge)%interact .EQV. .TRUE. .AND.    &
        pLevel%edgeCells(iEdge)%degenrt   ==    DEGENERAT_NONE ) THEN
      DO ijk = 1, UBOUND(pLevel%edgeCells(iedge)%cells,1)   
        pLevel%edgeCells(iEdge)%cells(ijk)%bufferExchPlag%nBuffSize  = 0
      ENDDO ! ijk
    ENDIF   ! interact
  ENDDO     ! iEdge

! Exit for null number of particles -------------------------------------------
!   after initializing buffer size to zero ------------------------------------

  IF ( nPcls == 0 ) GOTO 999

! Get cell dimensions ---------------------------------------------------------
  
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  
  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset(  region,iLev,iCOff,ijCOff )
  ibegDummCell = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iendDummCell = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

#ifdef PLAG_CECELLS_DEBUG
  WRITE(*,*) '  Inside PLAG_CECellsGetBufferSize: time,iReg,ijkDummBeg-End = ',&
 global%currentTime+ global%dtMin,iReg,idcbeg,idcend,jdcbeg,jdcend,kdcbeg,kdcend,ibegDummCell,iendDummCell
#endif 
       
! Get physical dimensions  ----------------------------------------------------
      
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibegPhysCell = IndIJK(ipcbeg,jpcbeg,kpcbeg,iCOff,ijCOff)
  iendPhysCell = IndIJK(ipcend,jpcend,kpcend,iCOff,ijCOff)

#ifdef PLAG_CECELLS_DEBUG
  WRITE(*,*) '  Inside PLAG_CECellsGetBufferSize: time,iReg,iPcls,ijkPhysBeg-End = ',&
      global%currentTime+ global%dtMin,iReg,ipcbeg,ipcend,jpcbeg,jpcend,kpcbeg,kpcend,ibegPhysCell,iendPhysCell
#endif 

! Loop over particles ---------------------------------------------------------

  DO iPcls = 1, nPcls
    ijkCPlag = pAiv(AIV_PLAG_ICELLS,iPcls)
    iCPlag   = pAiv(AIV_PLAG_INDEXI,iPcls)
    jCPlag   = pAiv(AIV_PLAG_INDEXJ,iPcls)
    kCPlag   = pAiv(AIV_PLAG_INDEXK,iPcls)    

    pclsFound = .FALSE.

! - Set lboundSkip(:) to 0 for values of lbound adjacent to a boundary -------

    lboundSkip(1:6) = 0 

    IF (iCPlag < ipcbeg) lboundSkip(1) = 1       ! to include ipcbeg = ipcend
    IF (iCPlag > ipcend) lboundSkip(2) = 1       ! case, ELSE IF not used here

    IF (jCPlag < jpcbeg) lboundSkip(3) = 1 
    IF (jCPlag > jpcend) lboundSkip(4) = 1

    IF (kCPlag < kpcbeg) lboundSkip(5) = 1
    IF (kCPlag > kpcend) lboundSkip(6) = 1

    lboundSkipSum = SUM(lboundSkip)

#ifdef PLAG_CECELLS_DEBUG
!   IF ( iPcls == 6 ) THEN
    WRITE(*,*) '  Inside PLAG_CECellsGetBufferSize: time,iReg,iPcls,iCPlag,jCPlag,kCPlag,',&
               ' ijkCPlag,lboundSkip,lboundSkipSum = ',&
        global%currentTime+ global%dtMin,iReg,iPcls,iCPlag,jCPlag,kCPlag, &
        ijkCPlag,lboundSkip,lboundSkipSum
!   ENDIF
#endif 
 
    IF ( lboundSkipSum == 0  ) GOTO 1999

#ifdef PLAG_CECELLS_DEBUG
!   IF ( iPcls == 6 ) THEN
    WRITE(*,*) '  Inside PLAG_CECellsGetBufferSize: Particle Fails Skip Test-',&
               'iReg,iPcls,iCPlag,jCPlag,kCPlag,ijkCPlag = ',&
        iReg,iPcls,iCPlag,jCPlag,kCPlag,ijkCPlag
!   ENDIF
#endif 

! - Corner Cells --------------------------------------------------------------

    DO iCorner=1,nCorners

! - Bypass for noninteracting regions -----------------------------------------

      IF( .NOT. pLevel%cornerCells(iCorner)%interact ) GOTO 2999

! - Bypass for degenerate corner cells ----------------------------------------

      IF( pLevel%cornerCells(iCorner)%degenrt /= DEGENERAT_NONE ) GOTO 2999

#ifdef PLAG_CECELLS_DEBUG
!       IF ( iPcls == 6 ) THEN
       WRITE(*,*) '  Inside PLAG_CECellsGetBufferSize: Entering Corner Test -iReg,iPcls,iCorner = ',&
        iReg,iPcls,iCorner
!       ENDIF
#endif 

! - Loop over corner cell indices ---------------------------------------------

      CALL RFLO_GetCornerCellsIndices( region,iLev,iCorner, &
                                       ibeg,iend,jbeg,jend,kbeg,kend )

      ibegCornCell = IndIJK(ibeg,jbeg,kbeg,iCOff,ijCOff)
      iendCornCell = IndIJK(iend,jend,kend,iCOff,ijCOff)

! - Increment buffer size 

      ijk = 0
      DO k = kbeg, kend
      DO j = jbeg, jend
      DO i = ibeg, iend
        ijk = ijk+1
        pCornCellsXBuff => pLevel%cornerCells(iCorner)%cells(ijk)%bufferExchPlag
  
        ijkCCell = IndIJK(i,j,k,iCOff,ijCOff)
        IF ( ijkCPlag == ijkCCell ) THEN
          pclsFound = .TRUE.
          pCornCellsXBuff%nBuffSize = pCornCellsXBuff%nBuffSize + 1
        ENDIF ! ijkCPlag
      ENDDO   ! i
      ENDDO   ! j
      ENDDO   ! k

#ifdef PLAG_CECELLS_DEBUG
!      IF ( iPcls == 6 ) THEN
        WRITE(*,*) & 
     '  Inside PLAG_CECellsGetBufferSize: iReg,iPcls,iCorner,ijkbeg-end,ibeg-endCornCell ',&
      ' pCornCellsXBuff%nBuffSize=',&
        iReg,iPcls,iEdge,ibeg,iend,jbeg,jend,kbeg,kend,ibegCornCell, iendCornCell, &
        pCornCellsXBuff%nBuffSize
!       ENDIF
#endif 

! -- Trap error ---------------------------------------------------------------

!      IF ( pCornCellsXBuff%nBuffSize > nBuffSizeTot ) & 
!        CALL ErrorStop( global, ERR_PLAG_ARRAYSIZE,__LINE__ )

2999  CONTINUE
        
    ENDDO       ! iCorner    

! If particle is found in corner, bypass edge search --------------------------

    IF (pclsFound) GOTO 1999
     
! Edge Cells ------------------------------------------------------------------

    DO iEdge=1,nEdges

! - Bypass for noninteracting regions -----------------------------------------

      IF( .NOT. pLevel%edgeCells(iEdge)%interact ) GOTO 3999

! - Bypass for degenerate edge cells ------------------------------------------

      IF( pLevel%edgeCells(iEdge)%degenrt /= DEGENERAT_NONE ) GOTO 3999
      
#ifdef PLAG_CECELLS_DEBUG
!      IF ( iPcls == 6 ) THEN
       WRITE(STDOUT,'(A,7(3X,I4))') &
       '  Inside PLAG_CECellsGetBufferSize: Entering Edge Test-iReg,iPcls,iEdge,ijkCPlag,i-j-kCPlag',&     
        iReg,iPcls,iEdge,ijkCPlag,iCPlag,jCPlag,kCPlag

!      END IF ! iPcls 
#endif 

! - Loop over edge cell indices -----------------------------------------------

      CALL RFLO_GetEdgeCellsIndices( region,iLev,iEdge, &
                                     ibeg,iend,jbeg,jend,kbeg,kend )
      ibegEdgeCell = IndIJK(ibeg,jbeg,kbeg,iCOff,ijCOff)
      iendEdgeCell = IndIJK(iend,jend,kend,iCOff,ijCOff)

#ifdef PLAG_CECELLS_DEBUG
  WRITE(STDOUT,'(A,12(2X,I4))') &
  '  Inside PLAG_CECellsGetBufferSize: iReg,iEdge,ijkbeg-end,ibeg-endEdgeCell =',&
        iReg,iEdge,ibeg,iend,jbeg,jend,kbeg,kend,ibegEdgeCell,iendEdgeCell
#endif 

! - Increment buffer size 

      ijk = 0
      DO k = kbeg, kend
      DO j = jbeg, jend
      DO i = ibeg, iend
        ijk = ijk+1    
        pEdgeCellsXBuff => pLevel%edgeCells(iEdge)%cells(ijk)%bufferExchPlag     

        ijkECell = IndIJK(i,j,k,iCOff,ijCOff)
        IF ( ijkCPlag == ijkECell ) THEN
          pclsFound = .TRUE.
          pEdgeCellsXBuff%nBuffSize = pEdgeCellsXBuff%nBuffSize + 1
        ENDIF ! ijkCPlag

#ifdef PLAG_CECELLS_DEBUG
!      IF ( iPcls == 6 ) THEN
!  WRITE(STDOUT,'(A,12(2X,I4))') &
!  '  Inside PLAG_CECellsGetBufferSize: iReg,iPcls,iEdge,ijkbeg-end,ibeg-endEdgeCell,pEdgeCellsXBuff%nBuffSize=',&
!        iReg,iPcls,iEdge,ibeg,iend,jbeg,jend,kbeg,kend,ibegEdgeCell, iendEdgeCell, &
!        pEdgeCellsXBuff%nBuffSize

  WRITE(STDOUT,'(A,12(2X,I4))') &
  '  Inside PLAG_CECellsGetBufferSize: iReg,iPcls,iEdge,ijk,ijkCPlag,ijkECell,pEdgeCellsXBuff%nBuffSize=',&
        iReg,iPcls,iEdge,i,j,k, &
        ijkCPlag, ijkECell,&
        pEdgeCellsXBuff%nBuffSize
!      ENDIF ! iPcls  
#endif

      ENDDO   ! i
      ENDDO   ! j
      ENDDO   ! k


! -- Trap error ---------------------------------------------------------------

!      IF ( pEdgeCellsXBuff%nBuffSize > nBuffSizeTot ) & 
!        CALL ErrorStop( global, ERR_PLAG_ARRAYSIZE,__LINE__ )

3999  CONTINUE
        
   ENDDO      ! iEdge

#ifdef PLAG_CECELLS_DEBUG
!   IF ( iPcls == 6 ) THEN
    WRITE(*,*) '  Inside PLAG_CECellsGetBufferSize: End of Particle Fails Skip Test-iReg,iPcls = ',&
      iReg,iPcls
!   ENDIF ! iPcls  
#endif

1999 CONTINUE
        
  ENDDO          ! iPcls

! finalize --------------------------------------------------------------------

999 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CECellsGetBufferSize

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CECellsGetBufferSize.F90,v $
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:13  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2004/11/29 19:26:06  fnajjar
! Added bypass statement for dengerate cells and debugging IO
!
! Revision 1.4  2004/03/20 23:36:43  fnajjar
! Initialized buffer sizes before checking particle size to alleviate data clobbering for zero nPcls
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







