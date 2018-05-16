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
! Purpose: Load buffer data from corner cells.
!
! Description: kernel is pertinent when the other region is located
!              on a different processor.
!
! Input: 
!   regions = data of all regions
!   iReg    = current region.
!   ir      = adjacent region
!   nBuffSizeEdge = buffer size for all edges
!
! Output: 
!   buffer data sent for all edges.
!   nBuffSizeCorn = buffer size for all corners
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_CornCellsLoadSendBuff.F90,v 1.3 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CornCellsLoadSendBuff( regions,iReg,ir,nBuffSizeEdge, &
                                       nBuffSizeCorn )

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
  INTEGER, INTENT(IN)     :: nBuffSizeEdge
  INTEGER, INTENT(OUT)    :: nBuffSizeCorn
 
! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: i,j,k,icell,iCorner,ijk,iLev,iRegDes,iRegSrc,nCorners
  INTEGER :: iAiv,iArv,iArvOld,iBuff,iBuffSendI,iBuffSendR,iCont,  &
             iCv,iCvMass,iCvOld,iRhs,iRhsSum,iShiftI,iShiftR,      &
             nArv,nAiv,nCont,nCv,nDimI,nDimR,nSendBuffI,nSendBuffR

  INTEGER,      POINTER, DIMENSION(:)   :: pCvPlagMass,pSendBuffI
  INTEGER,      POINTER, DIMENSION(:,:) :: pAivC, pAivOldC
  
  REAL(RFREAL), POINTER, DIMENSION(:)   :: pSendBuffR
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArvC,pArvOldC,pCvC,pCvOldC, &
                                           pRhsC,pRhsSumC 

  TYPE(t_buffer_plag), POINTER :: pCornCellsXBuff
  TYPE(t_dCellTransf), POINTER :: pSendEcCell
  TYPE(t_global),      POINTER :: global
  TYPE(t_level),       POINTER :: pLevelSrc
  TYPE(t_plag),        POINTER :: pPlag
  TYPE(t_region),      POINTER :: pRegionSrc

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CornCellsLoadSendBuff.F90,v $ $Revision: 1.3 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PLAG_CornCellsLoadSendBuff',&
  'PLAG_CornCellsLoadSendBuff.F90' )

! ******************************************************************************
! Get dimensions 
! ******************************************************************************

  iLev     = regions(iReg)%currLevel
  nCorners = 8

  nBuffSizeCorn = 0

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
!     Loading buffer data for corner
! =============================================================================

    DO iCorner=1,nCorners
      IF( .NOT. pLevelSrc%cornerCells(iCorner)%interact ) GOTO 2999

! -- Bypass for degenerate corner cells ---------------------------------------

      IF( pLevelSrc%cornerCells(iCorner)%degenrt /= DEGENERAT_NONE ) GOTO 2999

      iBuffSendI = nBuffSizeEdge; iBuffSendR = nBuffSizeEdge;
      iShiftI = nBuffSizeEdge; iShiftR = nBuffSizeEdge;

      DO ijk=1,UBOUND(pLevelSrc%cornerCells(iCorner)%cells,1)
        iRegDes = pLevelSrc%cornerCells(iCorner)%cells(ijk)%srcRegion

!------------------------------------------------------------------------------
!       Set pointers
!------------------------------------------------------------------------------

        pCornCellsXBuff => pLevelSrc%cornerCells(iCorner)%cells(ijk)%bufferExchPlag

        pAivC    => pCornCellsXBuff%aiv
        pArvC    => pCornCellsXBuff%arv    
        pCvC     => pCornCellsXBuff%cv
        pRhsC    => pCornCellsXBuff%rhs
        pRhsSumC => pCornCellsXBuff%rhsSum
    
        pAivOldC => pCornCellsXBuff%aivOld    
        pArvOldC => pCornCellsXBuff%arvOld
        pCvOldC  => pCornCellsXBuff%cvOld

        IF ( iRegDes == ir .AND. pCornCellsXBuff%nBuffSize /= 0 .AND. &
             regions(iRegDes)%procid /= global%myProcid               ) THEN 
          nBuffSizeCorn = nBuffSizeCorn +pCornCellsXBuff%nBuffSize

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         Load integer data buffers
!           compute shift and accumulate for various corner cells
!           indices of send buffers have to start where edges left off
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
          DO iBuff = 1, pCornCellsXBuff%nBuffSize
            iBuffSendI = iShiftI +nDimI*(iBuff-1) +1
            iAiv       = iBuffSendI

            pSendBuffI(iAiv  ) = pAivC(AIV_PLAG_PIDINI,iBuff)
            pSendBuffI(iAiv+1) = pAivC(AIV_PLAG_REGINI,iBuff)
            pSendBuffI(iAiv+2) = pAivC(AIV_PLAG_REGCRT,iBuff) 
            pSendBuffI(iAiv+3) = pAivC(AIV_PLAG_ICELLS,iBuff) 
            pSendBuffI(iAiv+4) = pAivC(AIV_PLAG_INDEXI,iBuff) 
            pSendBuffI(iAiv+5) = pAivC(AIV_PLAG_INDEXJ,iBuff)
            pSendBuffI(iAiv+6) = pAivC(AIV_PLAG_INDEXK,iBuff)
            pSendBuffI(iAiv+7) = pAivC(AIV_PLAG_BURNSTAT,iBuff)
            pSendBuffI(iAiv+8) = pAivC(AIV_PLAG_STATUS,iBuff)

#ifdef PLAG_CECELLS_MPI_DEBUG 
  IF(iReg==1)&                      
  WRITE(STDOUT,*) '  PLAG_CornCellsLoadSendBuff-INT: procSrc, iReg, procDes,iRegDes,iCorner,iBuff,iBuffSendI,iAiv = ',&
    global%myProcid,iReg,iReg,regions(iRegDes)%procid,iRegDes,iCorner,iBuff,iBuffSendI,iAiv
#endif

          ENDDO ! iBuff

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!          Load real data buffers
!            compute shift and accumulate for various edge cells
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          DO iBuff = 1, pCornCellsXBuff%nBuffSize
            iBuffSendR = iShiftR +nDimR*(iBuff-1) +1
            iCv       = iBuffSendR
            iRhs      = iBuffSendR +nCv 
            iRhsSum   = iBuffSendR +2*nCv 
            iCvOld    = iBuffSendR +3*nCv  
            iArv      = iBuffSendR +4*nCv 
            iArvOld   = iBuffSendR +4*nCv +nArv

#ifdef PLAG_CECELLS_MPI_DEBUG 
  IF(iReg==1)&       
  WRITE(STDOUT,*) '  PLAG_EdgeCellsLoadSendBuff--REAL:procSrc,iReg,procDes,iCorner,iBuff,iBuffSendR',&
  global%myProcid, iReg, regions(iRegDes)%procid, iRegDes, iCorner,iBuff, iBuffSendR
  IF(iReg==1)&       
  WRITE(STDOUT,*) '  PLAG_EdgeCellsLoadSendBuff--REAL:procSrc,iReg,procDes,iRegDes,iCv,iRhs,iRhsSum, iCvOld, iArv, iArvOld = ',&
     global%myProcid, iReg,regions(iRegDes)%procid,iRegDes,iCv,iRhs, iRhsSum, iCvOld, iArv, iArvOld 
#endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~           
!           Load real data buffers: cv
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            pSendBuffR(iCv  ) =  pCvC(CV_PLAG_XMOM,iBuff)
            pSendBuffR(iCv+1) =  pCvC(CV_PLAG_YMOM,iBuff) 
            pSendBuffR(iCv+2) =  pCvC(CV_PLAG_ZMOM,iBuff)
            pSendBuffR(iCv+3) =  pCvC(CV_PLAG_ENER,iBuff)
            pSendBuffR(iCv+4) =  pCvC(CV_PLAG_XPOS,iBuff)
            pSendBuffR(iCv+5) =  pCvC(CV_PLAG_YPOS,iBuff)
            pSendBuffR(iCv+6) =  pCvC(CV_PLAG_ZPOS,iBuff)
            pSendBuffR(iCv+7) =  pCvC(CV_PLAG_ENERVAPOR,iBuff)
            DO iCont = 1, nCont
              iCvMass = pCvPlagMass(iCont)
              pSendBuffR(iCv+(CV_PLAG_LAST-1)+iCont) = pCvC(iCvMass,iBuff)
            ENDDO ! iCont

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!           Load real data buffers: rhs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            pSendBuffR(iRhs  ) =  pRhsC(CV_PLAG_XMOM,iBuff)
            pSendBuffR(iRhs+1) =  pRhsC(CV_PLAG_YMOM,iBuff) 
            pSendBuffR(iRhs+2) =  pRhsC(CV_PLAG_ZMOM,iBuff)
            pSendBuffR(iRhs+3) =  pRhsC(CV_PLAG_ENER,iBuff)
            pSendBuffR(iRhs+4) =  pRhsC(CV_PLAG_XPOS,iBuff)
            pSendBuffR(iRhs+5) =  pRhsC(CV_PLAG_YPOS,iBuff)
            pSendBuffR(iRhs+6) =  pRhsC(CV_PLAG_ZPOS,iBuff)
            pSendBuffR(iRhs+7) =  pRhsC(CV_PLAG_ENERVAPOR,iBuff)
            DO iCont = 1, nCont
              iCvMass = pCvPlagMass(iCont)
              pSendBuffR(iRhs+(CV_PLAG_LAST-1)+iCont) = pRhsC(iCvMass,iBuff)
            ENDDO ! iCont

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!           Load real data buffers: rhsSum
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            pSendBuffR(iRhsSum  ) =  pRhsSumC(CV_PLAG_XMOM,iBuff)
            pSendBuffR(iRhsSum+1) =  pRhsSumC(CV_PLAG_YMOM,iBuff) 
            pSendBuffR(iRhsSum+2) =  pRhsSumC(CV_PLAG_ZMOM,iBuff)
            pSendBuffR(iRhsSum+3) =  pRhsSumC(CV_PLAG_ENER,iBuff)
            pSendBuffR(iRhsSum+4) =  pRhsSumC(CV_PLAG_XPOS,iBuff)
            pSendBuffR(iRhsSum+5) =  pRhsSumC(CV_PLAG_YPOS,iBuff)
            pSendBuffR(iRhsSum+6) =  pRhsSumC(CV_PLAG_ZPOS,iBuff)
            pSendBuffR(iRhsSum+7) =  pRhsSumC(CV_PLAG_ENERVAPOR,iBuff)
            DO iCont = 1, nCont
              iCvMass = pCvPlagMass(iCont)
              pSendBuffR(iRhsSum+(CV_PLAG_LAST-1)+iCont) = pRhsSumC(iCvMass,iBuff)
            ENDDO ! iCont          

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!           Load real data buffers: cvOld
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            pSendBuffR(iCvOld  ) =  pCvOldC(CV_PLAG_XMOM,iBuff)
            pSendBuffR(iCvOld+1) =  pCvOldC(CV_PLAG_YMOM,iBuff) 
            pSendBuffR(iCvOld+2) =  pCvOldC(CV_PLAG_ZMOM,iBuff)
            pSendBuffR(iCvOld+3) =  pCvOldC(CV_PLAG_ENER,iBuff)
            pSendBuffR(iCvOld+4) =  pCvOldC(CV_PLAG_XPOS,iBuff)
            pSendBuffR(iCvOld+5) =  pCvOldC(CV_PLAG_YPOS,iBuff)
            pSendBuffR(iCvOld+6) =  pCvOldC(CV_PLAG_ZPOS,iBuff)
            pSendBuffR(iCvOld+7) =  pCvOldC(CV_PLAG_ENERVAPOR,iBuff)
            DO iCont = 1, nCont
              iCvMass = pCvPlagMass(iCont)
              pSendBuffR(iCvOld+(CV_PLAG_LAST-1)+iCont) = pCvOldC(iCvMass,iBuff)
            ENDDO ! iCont

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!           Load real data buffers: arv, arvold  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            pSendBuffR(iArv)    =  pArvC(ARV_PLAG_SPLOAD,iBuff)
            pSendBuffR(iArvOld) =  pArvOldC(ARV_PLAG_SPLOAD,iBuff)

          ENDDO ! iBuff

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         Set appropriate shifts for integer and real send buffers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          iShiftI = iBuffSendI +nDimI -1
          iShiftR = iBuffSendR +nDimR -1

        ENDIF   ! iRegDes
      ENDDO     ! ijk

#ifdef PLAG_CECELLS_MPI_DEBUG
  IF ( nBuffSizeCorn > 0 ) & 
    PRINT*,' PLAG_CornCellsLoadSendBuff: procId, iReg, iR, procIdiR, iCorner, nBuffSizeEdge,iRegDes  = ',&
       global%myProcid, iReg, iR, regions(ir)%procid ,iCorner, nBuffSizeEdge,iRegDes
#endif

2999  CONTINUE
    ENDDO       ! iCorner
  ENDIF      ! some cells to send

! ******************************************************************************
! finalize
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CornCellsLoadSendBuff

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CornCellsLoadSendBuff.F90,v $
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:28  fnajjar
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







