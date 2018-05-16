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
! Purpose: Load buffer data from edge cells.
!
! Description: kernel is pertinent when the other region is located
!              on a different processor.
!
! Input: 
!   regions = data of all regions
!   iReg    = current region.
!   ir      = adjacent region
!
! Output: 
!   buffer data sent for all edges
!   nBuffSizeEdge = buffer size for all edges.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_EdgeCellsLoadSendBuff.F90,v 1.3 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_EdgeCellsLoadSendBuff( regions,iReg,ir,nBuffSizeEdge )

  USE ModDataTypes
  USE ModError
  USE ModMPI
  USE ModParameters
  USE ModGlobal,     ONLY : t_global
  USE ModDataStruct, ONLY : t_dCell, t_dCellTransf, t_region, t_level 
  USE ModPartLag,    ONLY : t_plag, t_buffer_plag

  USE PLAG_ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)

  INTEGER, INTENT(IN)     :: ir,iReg
  INTEGER, INTENT(OUT)    :: nBuffSizeEdge
 
! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: i,j,k,icell,iEdge,ijk,iLev,iRegDes,iRegSrc,nEdges
  INTEGER :: iAiv,iArv,iArvOld,iBuff,iBuffSendI, iBuffSendR,iCont, &
             iCv,iCvMass,iCvOld,iRhs,iRhsSum,iShiftI,iShiftR,      &
             nArv,nAiv,nCont,nCv,nDimI,nDimR,nSendBuffI,nSendBuffR

  INTEGER,      POINTER, DIMENSION(:)   :: pCvPlagMass,pSendBuffI
  INTEGER,      POINTER, DIMENSION(:,:) :: pAivE, pAivOldE

  REAL(RFREAL), POINTER, DIMENSION(:)   :: pSendBuffR
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArvE,pArvOldE,pCvE,pCvOldE, &
                                           pRhsE,pRhsSumE  

  TYPE(t_buffer_plag), POINTER :: pEdgeCellsXBuff
  TYPE(t_dCellTransf), POINTER :: pSendEcCell
  TYPE(t_global),      POINTER :: global
  TYPE(t_level),       POINTER :: pLevelSrc
  TYPE(t_plag),        POINTER :: pPlag
  TYPE(t_region),      POINTER :: pRegionSrc 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_EdgeCellsLoadSendBuff.F90,v $ $Revision: 1.3 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PLAG_EdgeCellsLoadSendBuff',&
  'PLAG_EdgeCellsLoadSendBuff.F90' )

! ******************************************************************************
! Get dimensions 
! ******************************************************************************

  iLev   = regions(iReg)%currLevel
  nEdges = 12

  nBuffSizeEdge = 0

! ******************************************************************************
! Set pointers 
! ******************************************************************************
  
  pRegionSrc  => regions(iReg)
  pLevelSrc   => pRegionSrc%levels(iLev)
  pPlag       => pLevelSrc%plag
  pCvPlagMass => pPlag%cvPlagMass

! ******************************************************************************
! Set send buffer dimensions
! ******************************************************************************

  nAiv  = pPlag%nAiv
  nArv  = pPlag%nArv
  nCont = pRegionSrc%plagInput%nCont
  nCv   = pPlag%nCv
       
  nDimI = nAiv
  nDimR = 2*nArv +4*nCv 

! ******************************************************************************
! Load send buffer data
! ******************************************************************************

  IF ( pLevelSrc%sendEcCells(ir)%nCells > 0 ) THEN
    pSendEcCell => pLevelSrc%sendEcCells(ir)
    pSendBuffI  => pSendEcCell%buffplagI
    pSendBuffR  => pSendEcCell%buffplagR 

! =============================================================================
!   Loop over edges of source region
!     Loading buffer data for edge
! =============================================================================

    DO iEdge=1,nEdges
      IF( .NOT. pLevelSrc%edgeCells(iEdge)%interact ) GOTO 2999

! -- Bypass for degenerate edge cells -----------------------------------------

      IF( pLevelSrc%edgeCells(iEdge)%degenrt /= DEGENERAT_NONE ) GOTO 2999
      
      iBuffSendI = 0; iBuffSendR = 0;
      iShiftI = 0; iShiftR = 0;

      DO ijk=1,UBOUND(pLevelSrc%edgeCells(iEdge)%cells,1)
        iRegDes = pLevelSrc%edgeCells(iEdge)%cells(ijk)%srcRegion

!------------------------------------------------------------------------------
!       Set pointers
!------------------------------------------------------------------------------

        pEdgeCellsXBuff => pLevelSrc%edgeCells(iEdge)%cells(ijk)%bufferExchPlag

        pAivE    => pEdgeCellsXBuff%aiv
        pArvE    => pEdgeCellsXBuff%arv    
        pCvE     => pEdgeCellsXBuff%cv
        pRhsE    => pEdgeCellsXBuff%rhs
        pRhsSumE => pEdgeCellsXBuff%rhsSum
    
        pAivOldE => pEdgeCellsXBuff%aivOld    
        pArvOldE => pEdgeCellsXBuff%arvOld
        pCvOldE  => pEdgeCellsXBuff%cvOld

        IF ( iRegDes == ir .AND. pEdgeCellsXBuff%nBuffSize /= 0 .AND. &
             regions(iRegDes)%procid /= global%myProcid               ) THEN 
          nBuffSizeEdge = nBuffSizeEdge +pEdgeCellsXBuff%nBuffSize

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         Load integer data buffers
!           compute shift and accumulate for various edge cells
!           provide initial shift to account for multiple particles 
!           in one cell
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          DO iBuff = 1, pEdgeCellsXBuff%nBuffSize
            iBuffSendI = iShiftI +nDimI*(iBuff-1) +1
            iAiv       = iBuffSendI

            pSendBuffI(iAiv  ) = pAivE(AIV_PLAG_PIDINI,iBuff)
            pSendBuffI(iAiv+1) = pAivE(AIV_PLAG_REGINI,iBuff)
            pSendBuffI(iAiv+2) = pAivE(AIV_PLAG_REGCRT,iBuff) 
            pSendBuffI(iAiv+3) = pAivE(AIV_PLAG_ICELLS,iBuff) 
            pSendBuffI(iAiv+4) = pAivE(AIV_PLAG_INDEXI,iBuff) 
            pSendBuffI(iAiv+5) = pAivE(AIV_PLAG_INDEXJ,iBuff)
            pSendBuffI(iAiv+6) = pAivE(AIV_PLAG_INDEXK,iBuff)
            pSendBuffI(iAiv+7) = pAivE(AIV_PLAG_BURNSTAT,iBuff)
            pSendBuffI(iAiv+8) = pAivE(AIV_PLAG_STATUS,iBuff)

#ifdef PLAG_CECELLS_MPI_DEBUG 
  IF ( iReg==1 ) WRITE(STDOUT,'(A,A,15(2X,I5))') &
    '  PLAG_EdgeCellsLoadSendBuff-INT: procSrc, iReg, procDes,',&
   '  iRegDes,iEdge,ijk,iBuff,iBuffSendI,iShiftI,iAiv = ',&
    global%myProcid,iReg,iReg,regions(iRegDes)%procid,iRegDes,iEdge, &
    ijk,iBuff,iBuffSendI,iShiftI,iAiv
  
  IF(iReg==1) WRITE(STDOUT,'(A,10(2X,I5))')'iBuff, iAiv,pSendBuffI =',&
             iBuff,iAiv, &
             pSendBuffI(iAiv  ),&
             pSendBuffI(iAiv+1),&
             pSendBuffI(iAiv+2),& 
             pSendBuffI(iAiv+3),& 
             pSendBuffI(iAiv+4),& 
             pSendBuffI(iAiv+5),&
             pSendBuffI(iAiv+6),&
             pSendBuffI(iAiv+7),&
             pSendBuffI(iAiv+8)
#endif
 
          ENDDO ! iBuff

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!          Load real data buffers
!            compute shift and accumulate for various edge cells
!           provide initial shift to account for multiple particles 
!           in one cell
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          DO iBuff = 1, pEdgeCellsXBuff%nBuffSize 
            iBuffSendR = iShiftR +nDimR*(iBuff-1) +1
            iCv        = iBuffSendR
            iRhs       = iBuffSendR +nCv 
            iRhsSum    = iBuffSendR +2*nCv 
            iCvOld     = iBuffSendR +3*nCv  
            iArv       = iBuffSendR +4*nCv 
            iArvOld    = iBuffSendR +4*nCv +nArv

#ifdef PLAG_CECELLS_MPI_DEBUG 
  IF ( iReg==1 ) WRITE(STDOUT,'(A,10(2X,I5))') &       
  '  PLAG_EdgeCellsLoadSendBuff--REAL:procSrc,iReg,procDes,iBuff,iBuffSendR,iShiftR',&
   global%myProcid, iReg, regions(iRegDes)%procid, iRegDes, iBuff, iBuffSendR,iShiftR
  IF ( iReg==1 ) WRITE(STDOUT,'(A,A,10(2X,I5))') &        
  '  PLAG_EdgeCellsLoadSendBuff--REAL:procSrc,iReg,procDes,iCv,',&
  'iRhs,iRhsSum, iCvOld, iArv, iArvOld = ',&
   global%myProcid, iReg,regions(iRegDes)%procid,iRegDes,iCv, &
   iRhs, iRhsSum, iCvOld, iArv, iArvOld 
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~           
!           Load real data buffers: cv
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            pSendBuffR(iCv  ) =  pCvE(CV_PLAG_XMOM,iBuff)
            pSendBuffR(iCv+1) =  pCvE(CV_PLAG_YMOM,iBuff) 
            pSendBuffR(iCv+2) =  pCvE(CV_PLAG_ZMOM,iBuff)
            pSendBuffR(iCv+3) =  pCvE(CV_PLAG_ENER,iBuff)
            pSendBuffR(iCv+4) =  pCvE(CV_PLAG_XPOS,iBuff)
            pSendBuffR(iCv+5) =  pCvE(CV_PLAG_YPOS,iBuff)
            pSendBuffR(iCv+6) =  pCvE(CV_PLAG_ZPOS,iBuff)
            pSendBuffR(iCv+7) =  pCvE(CV_PLAG_ENERVAPOR,iBuff)
            DO iCont = 1, nCont
              iCvMass = pCvPlagMass(iCont)
              pSendBuffR(iCv+(CV_PLAG_LAST-1)+iCont) = pCvE(iCvMass,iBuff)
            ENDDO ! iCont

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!           Load real data buffers: rhs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            pSendBuffR(iRhs  ) =  pRhsE(CV_PLAG_XMOM,iBuff)
            pSendBuffR(iRhs+1) =  pRhsE(CV_PLAG_YMOM,iBuff) 
            pSendBuffR(iRhs+2) =  pRhsE(CV_PLAG_ZMOM,iBuff)
            pSendBuffR(iRhs+3) =  pRhsE(CV_PLAG_ENER,iBuff)
            pSendBuffR(iRhs+4) =  pRhsE(CV_PLAG_XPOS,iBuff)
            pSendBuffR(iRhs+5) =  pRhsE(CV_PLAG_YPOS,iBuff)
            pSendBuffR(iRhs+6) =  pRhsE(CV_PLAG_ZPOS,iBuff)
            pSendBuffR(iRhs+7) =  pRhsE(CV_PLAG_ENERVAPOR,iBuff)
            DO iCont = 1, nCont
              iCvMass = pCvPlagMass(iCont)
              pSendBuffR(iRhs+(CV_PLAG_LAST-1)+iCont) = pRhsE(iCvMass,iBuff)
            ENDDO ! iCont

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!           Load real data buffers: rhsSum
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            pSendBuffR(iRhsSum  ) =  pRhsSumE(CV_PLAG_XMOM,iBuff)
            pSendBuffR(iRhsSum+1) =  pRhsSumE(CV_PLAG_YMOM,iBuff) 
            pSendBuffR(iRhsSum+2) =  pRhsSumE(CV_PLAG_ZMOM,iBuff)
            pSendBuffR(iRhsSum+3) =  pRhsSumE(CV_PLAG_ENER,iBuff)
            pSendBuffR(iRhsSum+4) =  pRhsSumE(CV_PLAG_XPOS,iBuff)
            pSendBuffR(iRhsSum+5) =  pRhsSumE(CV_PLAG_YPOS,iBuff)
            pSendBuffR(iRhsSum+6) =  pRhsSumE(CV_PLAG_ZPOS,iBuff)
            pSendBuffR(iRhsSum+7) =  pRhsSumE(CV_PLAG_ENERVAPOR,iBuff)
            DO iCont = 1, nCont
              iCvMass = pCvPlagMass(iCont)
              pSendBuffR(iRhsSum+(CV_PLAG_LAST-1)+iCont) = pRhsSumE(iCvMass,iBuff)
            ENDDO ! iCont          

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!           Load real data buffers: cvOld
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            pSendBuffR(iCvOld  ) =  pCvOldE(CV_PLAG_XMOM,iBuff)
            pSendBuffR(iCvOld+1) =  pCvOldE(CV_PLAG_YMOM,iBuff) 
            pSendBuffR(iCvOld+2) =  pCvOldE(CV_PLAG_ZMOM,iBuff)
            pSendBuffR(iCvOld+3) =  pCvOldE(CV_PLAG_ENER,iBuff)
            pSendBuffR(iCvOld+4) =  pCvOldE(CV_PLAG_XPOS,iBuff)
            pSendBuffR(iCvOld+5) =  pCvOldE(CV_PLAG_YPOS,iBuff)
            pSendBuffR(iCvOld+6) =  pCvOldE(CV_PLAG_ZPOS,iBuff)
            pSendBuffR(iCvOld+7) =  pCvOldE(CV_PLAG_ENERVAPOR,iBuff)
            DO iCont = 1, nCont
              iCvMass = pCvPlagMass(iCont)
              pSendBuffR(iCvOld+(CV_PLAG_LAST-1)+iCont) = pCvOldE(iCvMass,iBuff)
            ENDDO ! iCont

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!           Load real data buffers: arv, arvold  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            pSendBuffR(iArv)    =  pArvE(ARV_PLAG_SPLOAD,iBuff)
            pSendBuffR(iArvOld) =  pArvOldE(ARV_PLAG_SPLOAD,iBuff)

#ifdef PLAG_CECELLS_MPI_DEBUG 
  IF ( iReg==1 ) WRITE(STDOUT,'(A,A,10(2X,I5))') &                      
  '  PLAG_EdgeCellsLoadSendBuff-REAL: procSrc, iReg, procDes,iRegDes,',&
  'iEdge,iBuff,iBuffSendR,iCv = ',&
   global%myProcid,iReg,iReg,regions(iRegDes)%procid,iRegDes,&
   iEdge,iBuff,iBuffSendR,iCv
  
  IF  (iReg==1 ) WRITE(STDOUT,'(A,2(2X,I5),15(2X,1PE12.5))')'iBuff, iCv,pSendBuffR =',&
    iBuff,iCv, &
    pSendBuffR(iCv  ),&
    pSendBuffR(iCv+1),&
    pSendBuffR(iCv+2),& 
    pSendBuffR(iCv+3),& 
    pSendBuffR(iCv+4),& 
    pSendBuffR(iCv+5),&
    pSendBuffR(iCv+6),&
    pSendBuffR(iCv+7),&
    pSendBuffR(iCv+(CV_PLAG_LAST-1)+1:iCv+(CV_PLAG_LAST-1)+nCont)
#endif

          ENDDO ! iBuff

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         Set appropriate shifts for integer and real send buffers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          iShiftI = iBuffSendI +nDimI -1
          iShiftR = iBuffSendR +nDimR -1

        ENDIF   ! iRegDes
      ENDDO     ! ijk

#ifdef PLAG_CECELLS_MPI_DEBUG
  IF ( nBuffSizeEdge > 0 ) & 
    PRINT*,' PLAG_EdgeCellsLoadSendBuff: procId, iReg, iR, procIdiR, iEdge, nBuffSizeEdge,iRegDes  = ',&
       global%myProcid, iReg, iR, regions(ir)%procid ,iEdge, nBuffSizeEdge,iRegDes
#endif

2999  CONTINUE
    ENDDO       ! iEdge
  ENDIF      ! some cells to send

! ******************************************************************************
! finalize
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_EdgeCellsLoadSendBuff

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_EdgeCellsLoadSendBuff.F90,v $
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:02:42  fnajjar
! Initial revision after changing case
!
! Revision 1.3  2004/11/29 19:27:08  fnajjar
! Added bypass statement for dengerate cells
!
! Revision 1.2  2004/04/09 23:08:29  fnajjar
! Added AIV_PLAG_STATUS to send buffer
!
! Revision 1.1  2004/03/18 21:43:27  fnajjar
! Initial import for MPI-based data buffer communication
!
!******************************************************************************







