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
! Purpose: loads data buffer size for edge cells
!          and shrinks particle datastructure. 
!
! Description: none.
!
! Input: regions = data of all regions,
!        iReg    = current region number.
!
! Output: region%level%edgeCells%buffExchPlag%aiv,arv,cv,dv,tv   = buffer data.
!         region%level%cornerCells%buffExchPlag%aiv,arv,cv,dv,tv = buffer data.
!         region%level%plag%aiv,arv,cv,dv,tv = Plag data.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_EdgeCellsLoadData.F90,v 1.5 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_EdgeCellsLoadData( regions, iReg )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag, t_plag_input, t_buffer_plag
  USE ModDataStruct, ONLY : t_region, t_level, t_dCell
  USE ModGlobal, ONLY     : t_global  
  USE ModIndexing, ONLY   : GetIJK
  USE ModInterfaces, ONLY : RFLO_GetCellOffset,         &
                            RFLO_GetEdgeCellsIndices,   &
                            RFLO_GetDimensPhys
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER, INTENT(IN)     :: iReg

! ... loop variables
  INTEGER :: i, j, k, ijk, iEdge, iEdgeCellBuff, iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iLev, iPclsRegIn, nEdges, nPcls, nPclsPrev
  INTEGER :: iCOff, ijCOff
  INTEGER :: ijkCSrc, ijkESrc
  INTEGER :: iECDes, jECDes, kECDes, ijkECDes
  INTEGER :: iCPlag, ipcbeg, ipcend, ibeg, iend, ijkCPlag
  INTEGER :: jCPlag, jpcbeg, jpcend, jbeg, jend
  INTEGER :: kCPlag, kpcbeg, kpcend, kbeg, kend
  INTEGER :: iCOffDes, ijCOffDes, nDumCellsDes, iRegDes
  INTEGER :: ibegEdgeCell,iendEdgeCell 
  INTEGER :: errorFlag,iECMax
  INTEGER :: lPclsFoundInEdgeCellSum
  INTEGER :: nPclsBeg, nPclsEnd

  INTEGER,          DIMENSION(12)  :: lPclsFoundInEdgeCell
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv, pAivOld
  INTEGER, POINTER, DIMENSION(:,:) :: pAivC, pAivOldC
  INTEGER, POINTER, DIMENSION(:,:) :: pAivE, pAivOldE
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: edgeCellCounter

  LOGICAL :: pclsFoundInEdgeCell
  
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv, pArvOld, pCv, pCvOld, &
                                           pDv, pTv, pRhs, pRhsSum
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArvE, pArvOldE,pCvE, pCvOldE, &
                                           pDvE, pTvE, pRhsE, pRhsSumE
     
  TYPE(t_region),      POINTER :: pRegion
  TYPE(t_level),       POINTER :: pLevel
  TYPE(t_plag),        POINTER :: pPlag 
  TYPE(t_buffer_plag), POINTER :: pEdgeCellsXBuff
  TYPE(t_global),      POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PLAG_EdgeCellsLoadData.F90,v $ $Revision: 1.5 $'

  global => regions(iReg)%global
    
  CALL RegisterFunction( global, 'PLAG_EdgeCellsLoadData',&
  'PLAG_EdgeCellsLoadData.F90' )

! Get dimensions --------------------------------------------------------------

  iLev   = regions(iReg)%currLevel
  nPcls  = regions(iReg)%levels(iLev)%plag%nPcls 
  nEdges = 12

! Set pointers ----------------------------------------------------------------
  
  pRegion => regions(iReg)
  pLevel  => regions(iReg)%levels(iLev) 
  pPlag   => pLevel%plag   
  pAiv    => pPlag%aiv
  pArv    => pPlag%arv
  pCv     => pPlag%cv
  pDv     => pPlag%dv
  pTv     => pPlag%tv
  pRhs    => pPlag%rhs
  pRhsSum => pPlag%rhsSum

  pAivOld => pPlag%aivOld
  pArvOld => pPlag%arvOld
  pCvOld  => pPlag%cvOld

! Exit for null number of particles -------------------------------------------

  IF ( pPlag%nPcls == 0 ) GOTO 3999  
  
! Get grid dimensions ---------------------------------------------------------
  
  CALL RFLO_GetDimensPhys( pRegion,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )

! Get cell offset -------------------------------------------------------------

  CALL RFLO_GetCellOffset( pRegion,iLev,iCOff,ijCOff )
  
! Initialize counters for particles in inside and buffer regions --------------
  
  iPclsRegIn    = 0
  iEdgeCellBuff = 0

  pclsFoundInEdgeCell = .FALSE.
  lPclsFoundInEdgeCell = 0

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
    PRINT*,'PLAG_EdgeCellsLoadData: iReg,iECMax = ', iReg, iECMax
#endif

! Loop over particles ---------------------------------------------------------

  DO iPcls=1,nPcls

    iCPlag   = pAiv(AIV_PLAG_INDEXI,iPcls)
    jCPlag   = pAiv(AIV_PLAG_INDEXJ,iPcls)
    kCPlag   = pAiv(AIV_PLAG_INDEXK,iPcls)
    ijkCPlag = pAiv(AIV_PLAG_ICELLS,iPcls)
    
    lPclsFoundInEdgeCell = 0

! - Loop over edges -----------------------------------------------------------
   
    DO iEdge=1,nEdges 

! -- Bypass for noninteracting regions ----------------------------------------

      IF( .NOT. pLevel%edgeCells(iEdge)%interact ) GOTO 999

! -- Bypass for degenerate edge cells -----------------------------------------

      IF( pLevel%edgeCells(iEdge)%degenrt /= DEGENERAT_NONE ) GOTO 999

! -- Get edge cell indices ----------------------------------------------------

      CALL RFLO_GetEdgeCellsIndices( pRegion,iLev,iEdge, &
                                     ibeg,iend,jbeg,jend,kbeg,kend )

      ibegEdgeCell = IndIJK(ibeg,jbeg,kbeg,iCOff,ijCOff)
      iendEdgeCell = IndIJK(iend,jend,kend,iCOff,ijCOff)

! -- Check if particle cell is within edge cells extent -----------------------
! --  set flag to take care of multiple active edges --------------------------    

      IF ( iCPlag >= ibeg .AND. iCPlag <= iend .AND. &
           jCplag >= jbeg .AND. jCPlag <= jend .AND. &
           kCPlag >= kbeg .AND. kCPlag <= kend       ) THEN 
        lPclsFoundInEdgeCell(iEdge) = 1
        pclsFoundInEdgeCell = .TRUE. 
      ENDIF ! iCPlag

#ifdef PLAG_CECELLS_DEBUG
    PRINT*,'PLAG_EdgeCellsLoadData: iReg,iEdge,iPcls,ibeg,iend,jbeg,jend,kbeg,kend,',&
           'ibegEdgeCell, iendEdgeCell,i-j-kCPlag,ijkCPlag,nPcls,lPclsFoundInEdgeCell',&
            iReg,iEdge,iPcls,ibeg,iend,jbeg,jend,kbeg,kend, &
            ibegEdgeCell, iendEdgeCell,iCPlag,jCPlag,kCPlag,ijkCPlag,nPcls, &
            lPclsFoundInEdgeCell(iEdge) 
#endif

999 CONTINUE
    END DO ! iEdge
    
! - Loop over edges -----------------------------------------------------------

    DO iEdge=1,nEdges 

! -- Bypass for noninteracting regions ----------------------------------------

      IF( .NOT. pLevel%edgeCells(iEdge)%interact ) GOTO 2999

! -- Bypass for degenerate edge cells -----------------------------------------

      IF( pLevel%edgeCells(iEdge)%degenrt /= DEGENERAT_NONE ) GOTO 2999

! -- Get edge cell indices ----------------------------------------------------

      CALL RFLO_GetEdgeCellsIndices( pRegion,iLev,iEdge, &
                                     ibeg,iend,jbeg,jend,kbeg,kend )

! -- Determine sum of lPclsFoundInEdgeCell ------------------------------------
!      Such particle is in the physical domain --------------------------------

      lPclsFoundInEdgeCellSum = SUM( lPclsFoundInEdgeCell )     

! -- Particle is in physical domain -------------------------------------------
! -- Shift particle datastructure only if particle is not in its proper spot --
       
      IF ( lPclsFoundInEdgeCellSum == 0 ) THEN  

        iPclsRegIn = iPclsRegIn + 1 

#ifdef PLAG_CECELLS_DEBUG
    PRINT*,' PLAG_EdgeCellsLoadData: iReg, iEdge, iPcls,iPclsRegIn =',&
                                     iReg, iEdge, iPcls,iPclsRegIn
#endif

        IF ( iPclsRegIn /= iPcls ) THEN
          pAiv(   :,iPclsRegIn) = pAiv(   :,iPcls)
          pArv(   :,iPclsRegIn) = pArv(   :,iPcls)
          pCv(    :,iPclsRegIn) = pCv(    :,iPcls)
          pDv(    :,iPclsRegIn) = pDv(    :,iPcls)
          pTv(    :,iPclsRegIn) = pTv(    :,iPcls)
          pRhs(   :,iPclsRegIn) = pRhs(   :,iPcls)
          pRhsSum(:,iPclsRegIn) = pRhsSum(:,iPcls)

          pAivOld(:,iPclsRegIn) = pAivOld(:,iPcls)
          pArvOld(:,iPclsRegIn) = pArvOld(:,iPcls)
          pCvOld( :,iPclsRegIn) = pCvOld( :,iPcls)
        ENDIF ! iPclsRegIn

! -- Remove particle from active datastructure for specific edge --------------
                
      ENDIF ! lPclsFoundInEdgeCellSum
      
      IF ( lPclsFoundInEdgeCell(iEdge) == 1 ) THEN  

! -- Loop over edge cell indices ----------------------------------------------

        ijk = 0
        DO k=kbeg,kend
        DO j=jbeg,jend
        DO i=ibeg,iend
          ijk     = ijk + 1
          ijkESrc = IndIJK(i,j,k,iCOff, ijCOff)

! --- Set pointers ------------------------------------------------------------

          pEdgeCellsXBuff => pLevel%edgeCells(iEdge)%cells(ijk)%bufferExchPlag     
           
          pAivE    => pEdgeCellsXBuff%aiv
          pArvE    => pEdgeCellsXBuff%arv    
          pCvE     => pEdgeCellsXBuff%cv
          pDvE     => pEdgeCellsXBuff%dv
          pTvE     => pEdgeCellsXBuff%tv
          pRhsE    => pEdgeCellsXBuff%rhs
          pRhsSumE => pEdgeCellsXBuff%rhsSum
    
          pAivOldE => pEdgeCellsXBuff%aivOld    
          pArvOldE => pEdgeCellsXBuff%arvOld
          pCvOldE  => pEdgeCellsXBuff%cvOld
            
! --- Destination region infrastructure ---------------------------------------
            
          iRegDes      =  pLevel%edgeCells(iEdge)%cells(ijk)%srcRegion 

          IF ( iRegDes > 0 .AND. ijkCPlag == ijkESrc .AND. &
               pEdgeCellsXBuff%nBuffSize /= 0              ) THEN
            ijkECDes     = pLevel%edgeCells(iEdge)%cells(ijk)%srcCell           
            nDumCellsDes = regions(iRegDes)%nDumCells

            CALL RFLO_GetCellOffset( regions(iRegDes),iLev,iCOffDes,ijCOffDes )
            CALL GetIJK( ijkECDes,iCOffDes,ijCOffDes,nDumCellsDes, &
                         iECDes,jECDes,kECDes )

            edgeCellCounter(iEdge,ijk) = edgeCellCounter(iEdge,ijk)+1
            iEdgeCellBuff = edgeCellCounter(iEdge,ijk)

#ifdef PLAG_CECELLS_DEBUG
  WRITE(STDOUT,*) '      PLAG_EdgeCellsLoadData: iReg,iPcls,iEdge,iRegDes,iEdgeCellBuff ', &
       iReg, iPcls, iEdge, iRegDes, iEdgeCellBuff
#endif 
                    
! --- Update aiv field --------------------------------------------------------
    
            pAivE(AIV_PLAG_ICELLS,iEdgeCellBuff) = ijkECDes   
            pAivE(AIV_PLAG_INDEXI,iEdgeCellBuff) = iECDes
            pAivE(AIV_PLAG_INDEXJ,iEdgeCellBuff) = jECDes
            pAivE(AIV_PLAG_INDEXK,iEdgeCellBuff) = kECDes
            pAivE(AIV_PLAG_PIDINI,iEdgeCellBuff) = pAiv(AIV_PLAG_PIDINI,iPcls)
            pAivE(AIV_PLAG_REGINI,iEdgeCellBuff) = pAiv(AIV_PLAG_REGINI,iPcls)
            pAivE(AIV_PLAG_REGCRT,iEdgeCellBuff) = iRegDes
            pAivE(AIV_PLAG_BURNSTAT,iEdgeCellBuff) = pAiv(AIV_PLAG_BURNSTAT,iPcls)
            pAivE(AIV_PLAG_STATUS,iEdgeCellBuff) = pAiv(AIV_PLAG_STATUS,iPcls)
    
            pAivOldE(:           ,iEdgeCellBuff) = pAivE(:,iEdgeCellBuff) 
                       
! --- Load communication buffer arrays for corner cells -----------------------

            pArvE(   :,iEdgeCellBuff) = pArv(   :,iPcls)
            pCvE(    :,iEdgeCellBuff) = pCv(    :,iPcls)
            pDvE(    :,iEdgeCellBuff) = pDv(    :,iPcls)
            pTvE(    :,iEdgeCellBuff) = pTv(    :,iPcls)
            pRhsE(   :,iEdgeCellBuff) = pRhs(   :,iPcls)
            pRhsSumE(:,iEdgeCellBuff) = pRhsSum(:,iPcls)

            pArvOldE(:,iEdgeCellBuff) = pArvOld(:,iPcls)
            pCvOldE( :,iEdgeCellBuff) = pCvOld( :,iPcls)

#ifdef PLAG_CECELLS_DEBUG
  WRITE(STDOUT,*) &
  '      PLAG_EdgeCellsLoadData: iReg, iEdge, iEdgeCellBuff, nEdgeBuffSize, pAiv', &
  iReg, iEdge, iEdgeCellBuff, pEdgeCellsXBuff%nBuffSize,&
  pAivE(AIV_PLAG_ICELLS,iEdgeCellBuff),&
  pAivE(AIV_PLAG_INDEXI,iEdgeCellBuff),&
  pAivE(AIV_PLAG_INDEXJ,iEdgeCellBuff),&
  pAivE(AIV_PLAG_INDEXK,iEdgeCellBuff),&
  pAivE(AIV_PLAG_PIDINI,iEdgeCellBuff),&
  pAivE(AIV_PLAG_REGINI,iEdgeCellBuff),&
  pAivE(AIV_PLAG_REGCRT,iEdgeCellBuff),&
  pAivE(AIV_PLAG_BURNSTAT,iEdgeCellBuff)
#endif  

          ENDIF ! iRegDes

        ENDDO   ! i
        ENDDO   ! j
        ENDDO   ! k

      ENDIF   ! lPclsFoundInEdgeCell       

! ---- Exit edge search if particle is in physical domain ------------------

      IF ( lPclsFoundInEdgeCellSum == 0 ) GOTO 1999

#ifdef PLAG_CECELLS_DEBUG
  WRITE(STDOUT,'(A,2(2X,I3),2(2X,I4))') &
  '      PLAG_EdgeCellsLoadData: iReg, iEdge, iEdgeCellBuff', &
      iReg, iEdge, iEdgeCellBuff
#endif

2999 CONTINUE

    ENDDO    ! iEdge

1999 CONTINUE
  ENDDO      ! iPcls

! Get new particle datasize --------------------------------------------------

  nPclsPrev = pPlag%nPcls 
  IF ( pclsFoundInEdgeCell ) pPlag%nPcls = iPclsRegIn

#ifdef PLAG_CECELLS_DEBUG    
  WRITE(STDOUT,'(A,I4,2I8,2X,L1)') &
  '      PLAG_EdgeCellsLoadData: iReg,nPclsPrev,nPclsCurr,pclsFoundInEdgeCell = ',&
         iReg,nPclsPrev,pPlag%nPcls,pclsFoundInEdgeCell
#endif  

! reinitialize reshuffled particle datastructure to account for ---------------
!   region with null size particle --------------------------------------------
!   perform if data reshuffled and particle size null -------------------------

  IF ( pclsFoundInEdgeCell .AND. pPlag%nPcls == 0) THEN 
    nPclsBeg = MAX(1,pPlag%nPcls+1)
    nPclsEnd = nPclsPrev
  
    pPlag%aiv(:   ,nPclsBeg:nPclsEnd) = 0
    pPlag%aivOld(:,nPclsBeg:nPclsEnd) = 0
    pPlag%arv(:   ,nPclsBeg:nPclsEnd) = 0.0_RFREAL
    pPlag%arvOld(:,nPclsBeg:nPclsEnd) = 0.0_RFREAL
    pPlag%cv(:    ,nPclsBeg:nPclsEnd) = 0.0_RFREAL
    pPlag%cvOld(: ,nPclsBeg:nPclsEnd) = 0.0_RFREAL
    pPlag%rhs(:   ,nPclsBeg:nPclsEnd) = 0.0_RFREAL
    pPlag%rhsSum(:,nPclsBeg:nPclsEnd) = 0.0_RFREAL

#ifdef PLAG_CECELLS_DEBUG    
    WRITE(STDOUT,'(A,I4,2I8,2X,L1)') &
    '      PLAG_EdgeCellsLoadData: iReg,nPclsBeg,nPclsEnd = ',&
           iReg,nPclsBeg,nPclsEnd
#endif
  ENDIF ! pclsFoundInEdgeCell

! Deallocate edge cell buffer counter -----------------------------------------

  DEALLOCATE( edgeCellCounter,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= ERR_NONE) &
    CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

! finalize --------------------------------------------------------------------

3999 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_EdgeCellsLoadData

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_EdgeCellsLoadData.F90,v $
! Revision 1.5  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.2  2004/12/01 21:09:35  fnajjar
! Initial revision after changing case
!
! Revision 1.10  2004/11/29 19:27:08  fnajjar
! Added bypass statement for dengerate cells
!
! Revision 1.9  2004/04/09 23:09:16  fnajjar
! Added AIV_PLAG_STATUS to receive buffer and fixed nPclsBeg
!
! Revision 1.8  2004/03/20 21:53:16  fnajjar
! Included data reinitilization in an IF statement to alleviate data shrinkage
!
! Revision 1.7  2004/03/20 21:34:51  fnajjar
! Exit routine for null nPcls and reinitialized reshuffled data
!
! Revision 1.6  2004/03/19 23:51:30  fnajjar
! Reworked kernel to handle multiple active edges
!
! Revision 1.5  2004/03/18 21:41:52  fnajjar
! Various bug fixed for proper buffer loading
!
! Revision 1.4  2004/02/13 23:22:07  fnajjar
! Included new cv and aiv definitions for particle burning module
!
! Revision 1.3  2004/01/29 16:52:48  fnajjar
! Included search bypass for particles in physical domain
!
! Revision 1.2  2004/01/28 16:10:28  fnajjar
! Moved statements inside IF iReg for correct syntax
!
! Revision 1.1  2004/01/26 22:56:28  fnajjar
! Initial import for corner-edge cells to load buffer data
!
!******************************************************************************







