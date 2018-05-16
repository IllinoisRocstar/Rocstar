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
! Purpose: receive buffer data to edge and corner cells of an adjacent region. 
!          and append Lagrangian particle datastructure from buffers 
!          communicated via the corner and edge cell infrastructure.
!
! Description: kernel is pertinent when the other region is located
!              on a different processor.
!
! Input: 
!   regions = data of all regions
!   iReg    = current region.
!
! Output: buffer size received.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_CECellsRecvData.F90,v 1.5 2009/03/02 00:19:36 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CECellsRecvData( regions,iReg )

  USE ModDataTypes
  USE ModError
  USE ModMPI
  USE ModParameters
  USE ModGlobal,     ONLY : t_global
  USE ModDataStruct, ONLY : t_dCellTransf, t_level, t_region 
  USE ModPartLag,    ONLY : t_plag,t_buffer_plag

  USE PLAG_ModParameters

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

#ifdef MPI
  INTEGER :: statusPlag(MPI_STATUS_SIZE)
#endif

  INTEGER :: iAiv,iArv,iArvOld,iBuff,iBuffSendI,iBuffSendR,iCont,iCv, &
             iCvMass,iCvOld,iLev,iPcl,iPclBeg,iPclEnd,ir,iRhs,iRhsSum,&
             source,tagI,tagR
  INTEGER :: nArv,nAiv,nBuffSizeRecv,nCont,nCv,nDimI,nDimR,nPcls,nPclsPrev, &
             nRecvBuffI,nRecvBuffR
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass, pRecvBuffI 

  INTEGER, POINTER, DIMENSION(:,:) :: pAiv, pAivOld

  REAL(RFREAL), POINTER, DIMENSION(:)   :: pRecvBuffR
  
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv,pArvOld,pCv,pCvOld ,&
                                           pRhs,pRhsSum

  TYPE(t_dCellTransf), POINTER :: pRecvEcCell
  TYPE(t_global),      POINTER :: global
  TYPE(t_level),       POINTER :: pLevel
  TYPE(t_plag),        POINTER :: pPlag
  TYPE(t_region),      POINTER :: pRegion

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CECellsRecvData.F90,v $ $Revision: 1.5 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PLAG_CECellsRecvData',&
  'PLAG_CECellsRecvData.F90' )

! ******************************************************************************
! Set pointers 
! ******************************************************************************
  
  iLev = regions(iReg)%currLevel

  pRegion => regions(iReg)
  pLevel  => regions(iReg)%levels(iLev) 
  pPlag   => pLevel%plag
  
  pCvPlagMass => pPlag%cvPlagMass

  pAiv    => pPlag%aiv

  pArv    => pPlag%arv
  pCv     => pPlag%cv
  pRhs    => pPlag%rhs
  pRhsSum => pPlag%rhsSum
    
  pAivOld => pPlag%aivOld    
  pArvOld => pPlag%arvOld
  pCvOld  => pPlag%cvOld  

! ******************************************************************************
! Get dimensions 
! ******************************************************************************
  
  nCont = regions(iReg)%plagInput%nCont

  nAiv  = pPlag%nAiv
  nArv  = pPlag%nArv
  nCv   = pPlag%nCv

  nDimI = nAiv
  nDimR = 2*nArv +4*nCv

! ******************************************************************************
! Receive buffer data from source processor 
! ******************************************************************************

  DO ir=1,global%nRegions
    IF (regions(ir)%procid == global%myProcid) GOTO 999

    IF (pLevel%recvEcCells(ir)%nCells > 0) THEN
      pRecvEcCell => pLevel%recvEcCells(ir)
      pRecvBuffI  => pRecvEcCell%buffplagI
      pRecvBuffR  => pRecvEcCell%buffplagR 

! =============================================================================
!     Bypass MPI communication for null buffer size 
! =============================================================================       

      nBuffSizeRecv =  pRecvEcCell%nBuffSizePlag 

      IF ( nBuffSizeRecv == 0 ) GOTO 1999

! =============================================================================
!     Set receive buffer sizes 
! =============================================================================

      nRecvBuffI = nDimI * nBuffSizeRecv
      nRecvBuffR = nDimR * nBuffSizeRecv
      iBuffSendI = 0
      iBuffSendR = 0

#ifdef MPI
      source = regions(ir)%procid

! =============================================================================
!     Integer buffer data
! =============================================================================

      tagI    = regions(iReg)%localNumber +PLAG_TAG_SHIFT +MPI_PATCHOFF +2000
        IF(tagI .gt. global%mpiTagMax) tagI = MOD(tagI,global%mpiTagMax)

      CALL MPI_Recv( pRecvBuffI,nRecvBuffI,MPI_INTEGER, &
                     source,tagI,global%mpiComm,statusPlag,global%mpierr )
      IF ( global%mpierr /= ERR_NONE ) &
        CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

#ifdef PLAG_CECELLS_MPI_DEBUG 
   WRITE(STDOUT,'(A,6(2X,I5))') &
   '  PLAG_CECellsRecvData-INT: iRegDes, iRegSrc, procSrc, tagSrc, nBuffSizeRecv, nRecvBuffI  = ',&
      iReg, ir,source,tagI,nBuffSizeRecv,nRecvBuffI
#endif

! =============================================================================
!     Real buffer data
! =============================================================================

      tagR    = regions(iReg)%localNumber +PLAG_TAG_SHIFT +MPI_PATCHOFF +3000

        IF(tagR .gt. global%mpiTagMax) tagR = MOD(tagR,global%mpiTagMax)

      CALL MPI_Recv( pRecvBuffR,nRecvBuffR,MPI_RFREAL, &
                     source,tagR,global%mpiComm,statusPlag,global%mpierr )
      IF ( global%mpierr /= ERR_NONE ) &
        CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

#ifdef PLAG_CECELLS_MPI_DEBUG 
   WRITE(STDOUT,'(A,6(2X,I5))') & 
   '  PLAG_CECellsRecvSize-REAL: iRegDes, iRegSrc, procSrc, tagSrc, nBuffSizeRecv, nRecvBuffR  = ',&
      iReg, ir,source,tagR,nBuffSizeRecv,nRecvBuffR
#endif
#endif

#ifdef PLAG_CECELLS_MPI_DEBUG
 IF ( nBuffSizeRecv > 0 ) THEN
    DO iBuff=1,nBuffSizeRecv
      iBuffSendI = nDimI*(iBuff-1) +1
      iAiv       = iBuffSendI

     PRINT*,'iBuff,iBuffSendI,iAiv,pRecvBuffI = ',&
     iBuff,iBuffSendI,iAiv,&
     pRecvBuffI(iAiv  ),pRecvBuffI(iAiv+1),pRecvBuffI(iAiv+2),&
     pRecvBuffI(iAiv+3),pRecvBuffI(iAiv+4),pRecvBuffI(iAiv+5),&
     pRecvBuffI(iAiv+6),pRecvBuffI(iAiv+7)
     
     iBuffSendR = nDimR*(iBuff-1) +1
     iCv        = iBuffSendR
     iRhs       = iBuffSendR +nCv 
     iRhsSum    = iBuffSendR +2*nCv 
     iCvOld     = iBuffSendR +3*nCv  
     iArv       = iBuffSendR +4*nCv 
     iArvOld    = iBuffSendR +4*nCv +nArv

     PRINT*,'iBuff,iBuffSendR,iCv,pRecvBuffR = ',&
     iBuff,iCv,&
     pRecvBuffR(iCv  ),pRecvBuffR(iCv+1),pRecvBuffR(iCv+2),& 
     pRecvBuffR(iCv+3),pRecvBuffR(iCv+4),pRecvBuffR(iCv+5),& 
     pRecvBuffR(iCv+6),pRecvBuffR(iCv+7),pRecvBuffR(iCv+(CV_PLAG_LAST-1)+1:iCv+(CV_PLAG_LAST-1)+nCont)

    END DO ! iBuff   
 ENDIF ! nBuffSizeRecv
#endif

! =============================================================================
!   Append data from buffers to PLAG datastructure
! =============================================================================

!------------------------------------------------------------------------------
!   Set loop extent
!------------------------------------------------------------------------------
    
    nPcls = pPlag%nPcls 
    nPclsPrev = nPcls
          
    iPclBeg = nPcls +1
    iPclEnd = iPclBeg +(nBuffSizeRecv-1)

    iBuffSendI = 0
    iBuffSendR = 0

!------------------------------------------------------------------------------
!   Append to PLAG datastructure with receive buffer arrays 
!------------------------------------------------------------------------------
          
    DO iPcl = iPclBeg,iPclEnd 
      iBuff = iPcl-iPclBeg+1

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Load aiv from integer data buffers
!       compute shift and accumulate for various edge cells
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      iBuffSendI = nDimI*(iBuff-1) +1
      iAiv       = iBuffSendI

!     pAiv(:,iPcl) = pRecvBuffI(:  ) 

      pAiv(AIV_PLAG_PIDINI,iPcl) = pRecvBuffI(iAiv  )
      pAiv(AIV_PLAG_REGINI,iPcl) = pRecvBuffI(iAiv+1)
      pAiv(AIV_PLAG_REGCRT,iPcl) = pRecvBuffI(iAiv+2)
      pAiv(AIV_PLAG_ICELLS,iPcl) = pRecvBuffI(iAiv+3)
      pAiv(AIV_PLAG_INDEXI,iPcl) = pRecvBuffI(iAiv+4)
      pAiv(AIV_PLAG_INDEXJ,iPcl) = pRecvBuffI(iAiv+5)
      pAiv(AIV_PLAG_INDEXK,iPcl) = pRecvBuffI(iAiv+6) 
      pAiv(AIV_PLAG_BURNSTAT,iPcl) = pRecvBuffI(iAiv+7)
      pAiv(AIV_PLAG_STATUS,iPcl) = pRecvBuffI(iAiv+8)

      pAivOld(1:nAiv,iPcl) = pAiv(1:nAiv,iPcl)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Load from real data buffers
!       compute shift and accumulate for various edge cells
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      iBuffSendR = nDimR*(iBuff-1) +1
      iCv        = iBuffSendR
      iRhs       = iBuffSendR +nCv 
      iRhsSum    = iBuffSendR +2*nCv 
      iCvOld     = iBuffSendR +3*nCv  
      iArv       = iBuffSendR +4*nCv 
      iArvOld    = iBuffSendR +4*nCv +nArv

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~           
!     Load real data buffers: cv
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      pCv(CV_PLAG_XMOM,iPcl) = pRecvBuffR(iCv  )
      pCv(CV_PLAG_YMOM,iPcl) = pRecvBuffR(iCv+1)   
      pCv(CV_PLAG_ZMOM,iPcl) = pRecvBuffR(iCv+2) 
      pCv(CV_PLAG_ENER,iPcl) = pRecvBuffR(iCv+3)
      pCv(CV_PLAG_XPOS,iPcl) = pRecvBuffR(iCv+4)
      pCv(CV_PLAG_YPOS,iPcl) = pRecvBuffR(iCv+5) 
      pCv(CV_PLAG_ZPOS,iPcl) = pRecvBuffR(iCv+6)
      pCv(CV_PLAG_ENERVAPOR,iPcl) = pRecvBuffR(iCv+7) 
      DO iCont = 1, nCont
        iCvMass = pCvPlagMass(iCont)
        pCv(iCvMass,iPcl) = pRecvBuffR(iCv+(CV_PLAG_LAST-1)+iCont)
      ENDDO ! iCont

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Load real data buffers: rhs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      pRhs(CV_PLAG_XMOM,iPcl) = pRecvBuffR(iRhs  ) 
      pRhs(CV_PLAG_YMOM,iPcl) = pRecvBuffR(iRhs+1)   
      pRhs(CV_PLAG_ZMOM,iPcl) = pRecvBuffR(iRhs+2) 
      pRhs(CV_PLAG_ENER,iPcl) = pRecvBuffR(iRhs+3)  
      pRhs(CV_PLAG_XPOS,iPcl) = pRecvBuffR(iRhs+4) 
      pRhs(CV_PLAG_YPOS,iPcl) = pRecvBuffR(iRhs+5) 
      pRhs(CV_PLAG_ZPOS,iPcl) = pRecvBuffR(iRhs+6)
      pRhs(CV_PLAG_ENERVAPOR,iPcl)= pRecvBuffR(iRhs+7) 
      DO iCont = 1, nCont
        iCvMass = pCvPlagMass(iCont)
         pRhs(iCvMass,iPcl) = pRecvBuffR(iRhs+(CV_PLAG_LAST-1)+iCont) 
      ENDDO ! iCont

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      Load real data buffers: rhsSum
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      pRhsSum(CV_PLAG_XMOM,iPcl) = pRecvBuffR(iRhsSum  )
      pRhsSum(CV_PLAG_YMOM,iPcl) = pRecvBuffR(iRhsSum+1)   
      pRhsSum(CV_PLAG_ZMOM,iPcl) = pRecvBuffR(iRhsSum+2)  
      pRhsSum(CV_PLAG_ENER,iPcl) = pRecvBuffR(iRhsSum+3)  
      pRhsSum(CV_PLAG_XPOS,iPcl) = pRecvBuffR(iRhsSum+4) 
      pRhsSum(CV_PLAG_YPOS,iPcl) = pRecvBuffR(iRhsSum+5)  
      pRhsSum(CV_PLAG_ZPOS,iPcl) = pRecvBuffR(iRhsSum+6) 
      pRhsSum(CV_PLAG_ENERVAPOR,iPcl) = pRecvBuffR(iRhsSum+7)  
      DO iCont = 1, nCont
        iCvMass = pCvPlagMass(iCont)
        pRhsSum(iCvMass,iPcl) = pRecvBuffR(iRhsSum+(CV_PLAG_LAST-1)+iCont) 
      ENDDO ! iCont          

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      Load real data buffers: cvOld
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       pCvOld(CV_PLAG_XMOM,iPcl) = pRecvBuffR(iCvOld  ) 
       pCvOld(CV_PLAG_YMOM,iPcl) = pRecvBuffR(iCvOld+1) 
       pCvOld(CV_PLAG_ZMOM,iPcl) = pRecvBuffR(iCvOld+2) 
       pCvOld(CV_PLAG_ENER,iPcl) = pRecvBuffR(iCvOld+3) 
       pCvOld(CV_PLAG_XPOS,iPcl) = pRecvBuffR(iCvOld+4) 
       pCvOld(CV_PLAG_YPOS,iPcl) = pRecvBuffR(iCvOld+5)  
       pCvOld(CV_PLAG_ZPOS,iPcl) = pRecvBuffR(iCvOld+6) 
       pCvOld(CV_PLAG_ENERVAPOR,iPcl) = pRecvBuffR(iCvOld+7)
       DO iCont = 1, nCont
         iCvMass = pCvPlagMass(iCont)
         pCvOld(iCvMass,iPcl) = pRecvBuffR(iCvOld+(CV_PLAG_LAST-1)+iCont)
       ENDDO ! iCont

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Load real data buffers: arv, arvold  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

       pArv(ARV_PLAG_SPLOAD,iPcl) = pRecvBuffR(iArv) 
       pArvOld(ARV_PLAG_SPLOAD,iPcl) = pRecvBuffR(iArvOld) 
     END DO ! iPcl

! =============================================================================
!   Get new particle datasize
! =============================================================================
          
     nPcls = nPcls +nBuffSizeRecv
     pPlag%nPcls = nPcls

#ifdef PLAG_CECELLS_MPI_DEBUG
 IF ( nBuffSizeRecv > 0 ) THEN       
    WRITE(STDOUT,'(A,A,2X,1PE15.7,2X,3(I3,2X),5(I4,3X))') &
  '     PLAG_CECellsRecvData: time, procId,iReg, nBuffSizeRecv, ',&
  ' nPcls,nPclsPrev,iPclBeg,iPclEnd = ',&
    global%currentTime+global%dtMin,iReg,global%myProcid, &
    nBuffSizeRecv,nPcls,nPclsPrev,iPclBeg,iPclEnd 
    DO iPcl=1,nPcls
     PRINT*,'iPcl,aiv = ',&
     iPcl,&
     pAiv(AIV_PLAG_PIDINI,iPcl),&
     pAiv(AIV_PLAG_REGINI,iPcl),&
     pAiv(AIV_PLAG_REGCRT,iPcl),&
     pAiv(AIV_PLAG_ICELLS,iPcl),&
     pAiv(AIV_PLAG_INDEXI,iPcl),&
     pAiv(AIV_PLAG_INDEXJ,iPcl),&
     pAiv(AIV_PLAG_INDEXK,iPcl),& 
     pAiv(AIV_PLAG_BURNSTAT,iPcl),&
     pAiv(AIV_PLAG_STATUS,iPcl)
    END DO ! iPcl
    DO iPcl=1,nPcls
     PRINT*,'iPcl,cv = ',&
     iPcl,&
      pCv(CV_PLAG_XMOM,iPcl),&
      pCv(CV_PLAG_YMOM,iPcl),&   
      pCv(CV_PLAG_ZMOM,iPcl),& 
      pCv(CV_PLAG_ENER,iPcl),&
      pCv(CV_PLAG_XPOS,iPcl),&
      pCv(CV_PLAG_YPOS,iPcl),&
      pCv(CV_PLAG_ZPOS,iPcl),&
      pCv(CV_PLAG_ENERVAPOR,iPcl),&
      pCv(CV_PLAG_LAST+1:CV_PLAG_LAST+nCont,iPcl)
    END DO ! iPcl    
 ENDIF ! nBuffSizeRecv
#endif

1999 CONTINUE
    ENDIF      ! some cells to receive

999 CONTINUE
  ENDDO        ! ir

! ******************************************************************************
! finalize
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CECellsRecvData

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CECellsRecvData.F90,v $
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
! Revision 1.1  2004/12/01 20:57:15  fnajjar
! Initial revision after changing case
!
! Revision 1.2  2004/04/09 23:06:50  fnajjar
! Added AIV_PLAG_STATUS to receive buffer
!
! Revision 1.1  2004/03/18 21:43:27  fnajjar
! Initial import for MPI-based data buffer communication
!
!******************************************************************************







