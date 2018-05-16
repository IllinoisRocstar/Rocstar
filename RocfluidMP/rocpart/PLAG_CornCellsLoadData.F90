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
! Purpose: loads data buffer size for corner cells
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
! $Id: PLAG_CornCellsLoadData.F90,v 1.4 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CornCellsLoadData( regions, iReg )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag, t_plag_input, t_buffer_plag
  USE ModDataStruct, ONLY : t_region, t_level, t_dCell
  USE ModGlobal, ONLY     : t_global  
  USE ModIndexing, ONLY   : GetIJK
  USE ModInterfaces, ONLY : RFLO_GetCellOffset,         &
                            RFLO_GetCornerCellsIndices, &
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
  INTEGER :: i, j, k, ijk, iCorner, iCornCellBuff, iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iLev, iPclsRegIn, nCorners, nPcls, nPclsPrev
  INTEGER :: iCOff, ijCOff
  INTEGER :: ijkCSrc
  INTEGER :: iCCDes, jCCDes, kCCDes, ijkCCDes
  INTEGER :: iCPlag, ipcbeg, ipcend, ibeg, iend, ijkCPlag
  INTEGER :: jCPlag, jpcbeg, jpcend, jbeg, jend
  INTEGER :: kCPlag, kpcbeg, kpcend, kbeg, kend
  INTEGER :: iCOffDes, ijCOffDes, nDumCellsDes, iRegDes
  INTEGER :: ibegCornCell, iendCornCell
  INTEGER :: errorFlag,iCCMax 
  INTEGER :: lPclsFoundInCornCellSum
  INTEGER :: nPclsBeg, nPclsEnd

  INTEGER,          DIMENSION(9)   :: lPclsFoundInCornCell
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv, pAivOld
  INTEGER, POINTER, DIMENSION(:,:) :: pAivC, pAivOldC
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: cornCellCounter
  
  LOGICAL :: pclsFoundInCornCell

  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv, pArvOld, pCv, pCvOld, &
                                           pDv, pTv, pRhs, pRhsSum
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArvC, pArvOldC,pCvC, pCvOldC, &
                                           pDvC, pTvC, pRhsC, pRhsSumC
     
  TYPE(t_region),      POINTER :: pRegion
  TYPE(t_level),       POINTER :: pLevel
  TYPE(t_plag),        POINTER :: pPlag 
  TYPE(t_buffer_plag), POINTER :: pCornCellsXBuff
  TYPE(t_global),      POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PLAG_CornCellsLoadData.F90,v $ $Revision: 1.4 $'

  global => regions(iReg)%global
    
  CALL RegisterFunction( global, 'PLAG_CornCellsLoadData',&
  'PLAG_CornCellsLoadData.F90' )

! Get dimensions --------------------------------------------------------------

  iLev     = regions(iReg)%currLevel
  nPcls    = regions(iReg)%levels(iLev)%plag%nPcls 
  nCorners = 8

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
  iCornCellBuff = 0

  pclsFoundInCornCell = .FALSE.
  lPclsFoundInCornCell = 0

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
    PRINT*,'PLAG_CornCellsLoadData: iReg,iCCMax = ', iReg, iCCMax
#endif

! Loop over particles ---------------------------------------------------------

  DO iPcls=1,nPcls

    iCPlag   = pAiv(AIV_PLAG_INDEXI,iPcls)
    jCPlag   = pAiv(AIV_PLAG_INDEXJ,iPcls)
    kCPlag   = pAiv(AIV_PLAG_INDEXK,iPcls)
    ijkCPlag = pAiv(AIV_PLAG_ICELLS,iPcls)
    
    lPclsFoundInCornCell = 0

! - Loop over corners ---------------------------------------------------------

    DO iCorner=1,nCorners

! -- Bypass for noninteracting regions ----------------------------------------

      IF ( .NOT. pLevel%cornerCells(iCorner)%interact ) GOTO 999

! -- Bypass for degenerate corner cells ---------------------------------------

      IF( pLevel%cornerCells(iCorner)%degenrt /= DEGENERAT_NONE ) GOTO 999

! -- Get corner cell indices --------------------------------------------------

      CALL RFLO_GetCornerCellsIndices( pRegion,iLev,iCorner, &
                                       ibeg,iend,jbeg,jend,kbeg,kend )

      ibegCornCell = IndIJK(ibeg,jbeg,kbeg,iCOff,ijCOff)
      iendCornCell = IndIJK(iend,jend,kend,iCOff,ijCOff)

! -- Check if particle cell is within corner cells extent ---------------------
! --  set flag to take care of multiple active edges --------------------------    

      IF ( iCPlag >= ibeg .AND. iCPlag <= iend .AND. &
           jCplag >= jbeg .AND. jCPlag <= jend .AND. &
           kCPlag >= kbeg .AND. kCPlag <= kend       ) THEN 
        lPclsFoundInCornCell(iCorner) = 1
        pclsFoundInCornCell = .TRUE. 
      ENDIF ! iCPlag

#ifdef PLAG_CECELLS_DEBUG
    PRINT*,'PLAG_CornCellsLoadData: iReg,iCorner,iPcls,ibeg,iend,jbeg,jend,kbeg,kend,',&
           'ibegEdgeCell, iendEdgeCell,i-j-kCPlag,ijkCPlag,nPcls,lPclsFoundInEdgeCell',&
            iReg,iCorner,iPcls,ibeg,iend,jbeg,jend,kbeg,kend, &
            ibegCornCell, iendCornCell,iCPlag,jCPlag,kCPlag,ijkCPlag,nPcls, &
            lPclsFoundInCornCell(iCorner) 
#endif

999 CONTINUE
    END DO ! iCorner

! - Loop over corners ---------------------------------------------------------

    DO iCorner=1,nCorners

! -- Bypass for noninteracting regions ----------------------------------------

      IF ( .NOT. pLevel%cornerCells(iCorner)%interact ) GOTO 2999

! -- Bypass for degenerate corner cells ---------------------------------------

      IF( pLevel%cornerCells(iCorner)%degenrt /= DEGENERAT_NONE ) GOTO 2999

! -- Get corner cell indices --------------------------------------------------

      CALL RFLO_GetCornerCellsIndices( pRegion,iLev,iCorner, &
                                       ibeg,iend,jbeg,jend,kbeg,kend )

! -- Determine sum of lPclsFoundInCornCell ------------------------------------
!      Such particle is in the physical domain --------------------------------

      lPclsFoundInCornCellSum = SUM( lPclsFoundInCornCell )   

! -- Check if particle cell is within edge cells extent -----------------------
     
      IF ( iCPlag >= ibeg .AND. iCPlag <= iend .AND. &
           jCplag >= jbeg .AND. jCPlag <= jend .AND. &
           kCPlag >= kbeg .AND. kCPlag <= kend       ) & 
        pclsFoundInCornCell = .TRUE.

      IF ( lPclsFoundInCornCellSum == 0 ) THEN  

! --- Particle is in physical domain ------------------------------------------
! --- Shift particle datastructure only if particle is not in its proper spot -

        iPclsRegIn = iPclsRegIn + 1 

#ifdef PLAG_CECELLS_DEBUG
    PRINT*,' PLAG_CornCellsLoadData: iReg, iCorner, iPcls, iPclsRegIn =',&
                                     iReg, iCorner, iPcls, iPclsRegIn
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

      ENDIF ! lPclsFoundInCornCellSum

! -- Remove particle from active datastructure for specific corner ------------
     
      IF ( lPclsFoundInCornCell(iCorner) == 1 ) THEN

! -- Loop over corner cell indices --------------------------------------------

        ijk  = 0 
        DO k=kbeg,kend
        DO j=jbeg,jend
        DO i=ibeg,iend
          ijk     = ijk + 1
          ijkCSrc = IndIJK(i,j,k,iCOff, ijCOff)
 
! --- Set pointers ------------------------------------------------------------
            
          pCornCellsXBuff => pLevel%cornerCells(iCorner)%cells(ijk)%bufferExchPlag  
           
          pAivC    => pCornCellsXBuff%aiv
          pArvC    => pCornCellsXBuff%arv    
          pCvC     => pCornCellsXBuff%cv
          pDvC     => pCornCellsXBuff%dv
          pTvC     => pCornCellsXBuff%tv
          pRhsC    => pCornCellsXBuff%rhs
          pRhsSumC => pCornCellsXBuff%rhsSum
    
          pAivOldC => pCornCellsXBuff%aivOld    
          pArvOldC => pCornCellsXBuff%arvOld
          pCvOldC  => pCornCellsXBuff%cvOld
           
! --- Destination region infrastructure --------------------------------------
            
          iRegDes      =  pLevel%cornerCells(iCorner)%cells(ijk)%srcRegion 

          IF ( iRegDes > 0 .AND. ijkCPlag == ijkCSrc .AND. &
               pCornCellsXBuff%nBuffSize /= 0              ) THEN
            ijkCCDes     = pLevel%cornerCells(iCorner)%cells(ijk)%srcCell           
            nDumCellsDes = regions(iRegDes)%nDumCells

            CALL RFLO_GetCellOffset( regions(iRegDes),iLev,iCOffDes,ijCOffDes )
            CALL GetIJK( ijkCCDes,iCOffDes,ijCOffDes,nDumCellsDes, &
                         iCCDes,jCCDes,kCCDes )

            cornCellCounter(iCorner,ijk) = cornCellCounter(iCorner,ijk)+1
            iCornCellBuff = cornCellCounter(iCorner,ijk)

#ifdef PLAG_CECELLS_DEBUG
  WRITE(STDOUT,*) '      PLAG_CornCellsLoadData: iReg,iPcls,iCorner,iRegDes,iCornCellBuff ', &
       iReg,iPcls,iCorner,iRegDes,iCornCellBuff
#endif 

! --- Update aiv field --------------------------------------------------------
    
            pAivC(AIV_PLAG_ICELLS,iCornCellBuff) = ijkCCDes   
            pAivC(AIV_PLAG_INDEXI,iCornCellBuff) = iCCDes
            pAivC(AIV_PLAG_INDEXJ,iCornCellBuff) = jCCDes
            pAivC(AIV_PLAG_INDEXK,iCornCellBuff) = kCCDes
            pAivC(AIV_PLAG_PIDINI,iCornCellBuff) = pAiv(AIV_PLAG_PIDINI,iPcls)
            pAivC(AIV_PLAG_REGINI,iCornCellBuff) = pAiv(AIV_PLAG_REGINI,iPcls)
            pAivC(AIV_PLAG_REGCRT,iCornCellBuff) = iRegDes
            pAivC(AIV_PLAG_BURNSTAT,iCornCellBuff) = pAiv(AIV_PLAG_BURNSTAT,iPcls)
            pAivC(AIV_PLAG_STATUS,iCornCellBuff) = pAiv(AIV_PLAG_STATUS,iPcls)
    
            pAivOldC(:           ,iCornCellBuff) = pAivC(:,iCornCellBuff) 
                       
! --- Load communication buffer arrays for corner cells -----------------------

            pArvC(   :,iCornCellBuff) = pArv(   :,iPcls)
            pCvC(    :,iCornCellBuff) = pCv(    :,iPcls)
            pDvC(    :,iCornCellBuff) = pDv(    :,iPcls)
            pTvC(    :,iCornCellBuff) = pTv(    :,iPcls)
            pRhsC(   :,iCornCellBuff) = pRhs(   :,iPcls)
            pRhsSumC(:,iCornCellBuff) = pRhsSum(:,iPcls)

            pArvOldC(:,iCornCellBuff) = pArvOld(:,iPcls)
            pCvOldC( :,iCornCellBuff) = pCvOld( :,iPcls)

#ifdef PLAG_CECELLS_DEBUG
  WRITE(STDOUT,*) &
  '      PLAG_CornCellsLoadData: iReg, iCorner, iCornCellBuff, nCornBuffSize, pAiv', &
  iReg,iCorner, iCornCellBuff, pCornCellsXBuff%nBuffSize,&
  pAivC(AIV_PLAG_ICELLS,iCornCellBuff),&
  pAivC(AIV_PLAG_INDEXI,iCornCellBuff),&
  pAivC(AIV_PLAG_INDEXJ,iCornCellBuff),&
  pAivC(AIV_PLAG_INDEXK,iCornCellBuff),&
  pAivC(AIV_PLAG_PIDINI,iCornCellBuff),&
  pAivC(AIV_PLAG_REGINI,iCornCellBuff),&
  pAivC(AIV_PLAG_REGCRT,iCornCellBuff),&
  pAivC(AIV_PLAG_BURNSTAT,iCornCellBuff)
#endif 

          ENDIF ! iRegDes

        ENDDO   ! i
        ENDDO   ! j
        ENDDO   ! k

      ENDIF   ! pclsFoundInCornCell         

! ---- Exit edge search if particle is in physical domain ------------------

      IF ( lPclsFoundInCornCellSum == 0 ) GOTO 1999

#ifdef PLAG_CECELLS_DEBUG
  WRITE(STDOUT,'(A,2(2X,I3),2(2X,I4))') &
  '      PLAG_CornCellsLoadData: iReg, iCorner, iCornCellBuff', &
      iReg, iCorner, iCornCellBuff
#endif

2999 CONTINUE

    ENDDO    ! iCorner

1999 CONTINUE
  ENDDO      ! iPcls

! Get new particle datasize --------------------------------------------------

  nPclsPrev = pPlag%nPcls 
  IF ( pclsFoundInCornCell ) pPlag%nPcls = iPclsRegIn

#ifdef PLAG_CECELLS_DEBUG    
  WRITE(STDOUT,'(A,I4,2I8,2X,L1)') &
  '      PLAG_CornCellsLoadData: iReg,nPclsPrev,nPclsCurr,pclsFoundInCornCell = ',&
         iReg,nPclsPrev,pPlag%nPcls,pclsFoundInCornCell
#endif

! reinitialize reshuffled particle datastructure to account for ---------------
!   region with null size particle --------------------------------------------
!   perform if data reshuffled and particle size null -------------------------
  
  IF ( pclsFoundInCornCell .AND. pPlag%nPcls == 0) THEN 
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
    '      PLAG_CornCellsLoadData: iReg,nPclsBeg,nPclsEnd = ',&
           iReg,nPclsBeg,nPclsEnd
#endif
  ENDIF ! pclsFoundInCornCell

! Deallocate corner cell buffer counter ---------------------------------------

  DEALLOCATE( cornCellCounter,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= ERR_NONE) &
    CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

! finalize --------------------------------------------------------------------
3999 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CornCellsLoadData

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CornCellsLoadData.F90,v $
! Revision 1.4  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 20:57:27  fnajjar
! Initial revision after changing case
!
! Revision 1.10  2004/11/29 19:27:08  fnajjar
! Added bypass statement for dengerate cells
!
! Revision 1.9  2004/04/09 23:07:36  fnajjar
! Added AIV_PLAG_STATUS to receive buffer and fixed nPclsBeg
!
! Revision 1.8  2004/03/20 21:53:15  fnajjar
! Included data reinitilization in an IF statement to alleviate data shrinkage
!
! Revision 1.7  2004/03/20 21:34:23  fnajjar
! Exit routine for null nPcls and reinitialized reshuffled data
!
! Revision 1.6  2004/03/19 23:51:08  fnajjar
! Reworked kernel to handle multiple active corners
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







