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
! Purpose: send buffer data from edge and corner cells to an adjacent region.
!
! Description: kernel is pertinent when the other region is located
!              on a different processor.
!
! Input: 
!   regions = data of all regions
!   iReg    = current region.
!
! Output: 
!   buffer data sent.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_CECellsSendData.F90,v 1.4 2009/03/02 00:19:36 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CECellsSendData( regions,iReg )

  USE ModDataTypes
  USE ModError
  USE ModMPI
  USE ModParameters
  USE ModGlobal,     ONLY : t_global
  USE ModDataStruct, ONLY : t_dCell, t_dCellTransf, t_region, t_level 
  USE ModPartLag,    ONLY : t_plag, t_buffer_plag

  USE PLAG_ModParameters
  
  USE PLAG_ModInterfaces, ONLY : PLAG_EdgeCellsLoadSendBuff, &
                                 PLAG_CornCellsLoadSendBuff

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

  INTEGER :: dest,iLev,ir,iRegDes,iRegSrc,iRequestPlag,tagI,tagR
  INTEGER :: nBuffSizeCornSrc,nBuffSizeEdgeSrc,nBuffSizeTotSrc 
  INTEGER :: nArv,nAiv,nCont,nCv,nDimI,nDimR,nSendBuffI,nSendBuffR
  INTEGER,      POINTER, DIMENSION(:)   :: pSendBuffI
  
  REAL(RFREAL), POINTER, DIMENSION(:)   :: pSendBuffR

  TYPE(t_dCellTransf), POINTER :: pSendEcCell
  TYPE(t_global),      POINTER :: global
  TYPE(t_level),       POINTER :: pLevelSrc
  TYPE(t_plag),        POINTER :: pPlag
  TYPE(t_region),      POINTER :: pRegionSrc

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CECellsSendData.F90,v $ $Revision: 1.4 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PLAG_CECellsSendData',&
  'PLAG_CECellsSendData.F90' )

! ******************************************************************************
! Get dimensions 
! ******************************************************************************

  iLev = regions(iReg)%currLevel

  nBuffSizeEdgeSrc = 0
  nBuffSizeCornSrc = 0
  nBuffSizeTotSrc  = 0

! ******************************************************************************
! Set pointers 
! ******************************************************************************
  
  pRegionSrc  => regions(iReg)
  pLevelSrc   => pRegionSrc%levels(iLev)
  pPlag       => pLevelSrc%plag

! ******************************************************************************
! Set send buffer dimensions
! ******************************************************************************

  nAiv  = pPlag%nAiv
  nArv  = pPlag%nArv
  nCv   = pPlag%nCv

  nDimI = nAiv
  nDimR = 2*nArv +4*nCv     

! ******************************************************************************
! Load send buffer data
! ******************************************************************************

  DO ir=1,global%nRegions
    IF (regions(ir)%procid == global%myProcid) GOTO 999

    IF ( pLevelSrc%sendEcCells(ir)%nCells > 0 ) THEN
      pSendEcCell => pLevelSrc%sendEcCells(ir)
      pSendBuffI  => pSendEcCell%buffplagI
      pSendBuffR  => pSendEcCell%buffplagR 

! =============================================================================
!     Bypass MPI communication for null buffer size 
! =============================================================================       

      IF ( pSendEcCell%nBuffSizePlag == 0 ) GOTO 1999

! =============================================================================
!     Load buffer data for edges
! =============================================================================

      CALL PLAG_EdgeCellsLoadSendBuff( regions,iReg,ir,nBuffSizeEdgeSrc )

! =============================================================================
!     Load buffer data for corners
! =============================================================================

      CALL PLAG_CornCellsLoadSendBuff( regions,iReg,ir,nBuffSizeEdgeSrc, &
                                       nBuffSizeCornSrc )

! =============================================================================
!     Trap error for inconsistent buffer sizes
! =============================================================================

      nBuffSizeTotSrc=  nBuffSizeCornSrc+nBuffSizeEdgeSrc
      IF ( nBuffSizeTotSrc /= pSendEcCell%nBuffSizePlag ) THEN
        WRITE(STDOUT,*) 'PLAG_CECellsSendData: Error inconsistent buffer sizes'
        WRITE(STDOUT,*) '                      nBuffSizeEdgeSrc = ', nBuffSizeEdgeSrc
        WRITE(STDOUT,*) '                      nBuffSizeCornSrc = ', nBuffSizeCornSrc
        WRITE(STDOUT,*) '                      nBuffSizeTotSrc = ', nBuffSizeTotSrc
        WRITE(STDOUT,*) '            pSendEcCell%nBuffSizePlag = ', pSendEcCell%nBuffSizePlag

! TEMPORARY
!        CALL ErrorStop(global,ERR_PLAG_BUFFSIZE,__LINE__,errorString)
#ifdef MPI
         CALL MPI_Finalize(global%mpierr)
#endif
         STOP
! END TEMPORARY
      ENDIF ! nBuffSizeTotSrc

! =============================================================================
!     Send buffer data to destination processor 
! =============================================================================

      nSendBuffI = nDimI * pSendEcCell%nBuffSizePlag
      nSendBuffR = nDimR * pSendEcCell%nBuffSizePlag

      iRequestPlag = pSendEcCell%iRequestPlag

#ifdef MPI
      dest = regions(ir)%procid

!------------------------------------------------------------------------------
!     Integer buffer data
!------------------------------------------------------------------------------

      tagI  = regions(ir)%localNumber +PLAG_TAG_SHIFT +MPI_PATCHOFF +2000

      IF(tagI .gt. global%mpiTagMax) tagI = MOD(tagI,global%mpiTagMax)

#ifdef PLAG_CECELLS_MPI_DEBUG 
  WRITE(STDOUT,'(A,A,7(2X,I5))') &
  '  PLAG_CECellsSendData-INT: iRegDes, iRegSrc, procDes, procSrc,',&
  'tagSrc, nBuffSizePlag,nSendBuffI  = ',&
    ir, iReg, dest, global%myProcid,tagI, pSendEcCell%nBuffSizePlag,nSendBuffI
#endif

      CALL MPI_Isend( pSendBuffI,nSendBuffI,MPI_INTEGER, &
                      dest,tagI,global%mpiComm,                            &
                      pPlag%requestsCECellsI(iRequestPlag),global%mpierr   )
      IF ( global%mpierr /= ERR_NONE ) &
        CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

!------------------------------------------------------------------------------
!     Real buffer data
!------------------------------------------------------------------------------

      tagR  = regions(ir)%localNumber +PLAG_TAG_SHIFT +MPI_PATCHOFF +3000
      IF(tagR .gt. global%mpiTagMax) tagR = MOD(tagR,global%mpiTagMax)

#ifdef PLAG_CECELLS_MPI_DEBUG 
  WRITE(STDOUT,'(A,A,7(2X,I5))') &
  '  PLAG_CECellsSendData-REAL: iRegDes, iRegSrc, procDes, procSrc,',&
  'tagSrc, nBuffSizePlag,nSendBuffR  = ',&
    ir, iReg, dest, global%myProcid,tagR, pSendEcCell%nBuffSizePlag,nSendBuffR
#endif

      CALL MPI_Isend( pSendBuffR,nSendBuffR,MPI_RFREAL, &
                      dest,tagR,global%mpiComm,                            &
                      pPlag%requestsCECellsR(iRequestPlag),global%mpierr   )
      IF ( global%mpierr /= ERR_NONE ) &
        CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

#endif
1999 CONTINUE
    ENDIF      ! some cells to send

999 CONTINUE
  ENDDO        ! ir

! ******************************************************************************
! finalize
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CECellsSendData

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CECellsSendData.F90,v $
! Revision 1.4  2009/03/02 00:19:36  mtcampbe
! Added some ifdefs around Rocflo to disable particle injection on INFLOW
! boundaries and added some checks around MPI tags utilizing a new global
! data item, global%mpiTagMax.
!
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:17  fnajjar
! Initial revision after changing case
!
! Revision 1.2  2004/03/20 21:32:34  fnajjar
! Added more writing statement for error trapping
!
! Revision 1.1  2004/03/18 21:43:27  fnajjar
! Initial import for MPI-based data buffer communication
!
!******************************************************************************







