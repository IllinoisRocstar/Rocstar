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
! Purpose: receive buffer size to edge and corner cells of an adjacent region. 
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
! $Id: PLAG_CECellsRecvSize.F90,v 1.4 2009/03/02 00:19:36 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CECellsRecvSize( regions,iReg )

  USE ModDataTypes
  USE ModError
  USE ModMPI
  USE ModParameters
  USE ModGlobal,     ONLY : t_global
  USE ModDataStruct, ONLY : t_dCellTransf, t_level, t_region 

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

  INTEGER ::  iLev,ir,nDimBuffSize,source,tag

  TYPE(t_dCellTransf), POINTER :: pRecvEcCell
  TYPE(t_global),      POINTER :: global
  TYPE(t_level),       POINTER :: pLevel
  TYPE(t_region),      POINTER :: pRegion 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CECellsRecvSize.F90,v $ $Revision: 1.4 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PLAG_CECellsRecvSize',&
  'PLAG_CECellsRecvSize.F90' )

! ******************************************************************************
! Get dimensions 
! ******************************************************************************

  iLev = regions(iReg)%currLevel
  
  nDimBuffSize = 1 

! ******************************************************************************
! Set pointers 
! ******************************************************************************
 
  pRegion => regions(iReg)
  pLevel  => regions(iReg)%levels(iLev) 

! ******************************************************************************
! Receive buffer size from source processor 
! ******************************************************************************

  DO ir=1,global%nRegions
    IF (regions(ir)%procid /= global%myProcid) THEN
    IF (pLevel%recvEcCells(ir)%nCells > 0) THEN
      pRecvEcCell => pLevel%recvEcCells(ir)

#ifdef MPI
      source = regions(ir)%procid
      tag    = regions(iReg)%localNumber +PLAG_TAG_SHIFT +MPI_PATCHOFF +1000

        IF(tag .gt. global%mpiTagMax) tag = MOD(tag,global%mpiTagMax)
      CALL MPI_Recv( pRecvEcCell%nBuffSizePlag,nDimBuffSize,MPI_INTEGER, &
                     source,tag,global%mpiComm,statusPlag,global%mpierr )
      IF ( global%mpierr /= ERR_NONE ) &
        CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

#ifdef PLAG_CECELLS_MPI_DEBUG 
   IF ( pRecvEcCell%nBuffSizePlag > 0 ) &                 
     WRITE(STDOUT,*) '  PLAG_CECellsRecvSize: iRegDes, iRegSrc, procSrc, tagSrc, nBuffSizePlag  = ',&
      iReg, ir,source,tag,pRecvEcCell%nBuffSizePlag
#endif

    ENDIF      ! some cells to receive
    ENDIF      ! not my processor
  ENDDO        ! ir

! ******************************************************************************
! finalize
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CECellsRecvSize

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CECellsRecvSize.F90,v $
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
! Revision 1.1  2004/12/01 20:57:16  fnajjar
! Initial revision after changing case
!
! Revision 1.1  2004/03/10 23:16:09  fnajjar
! Initial import of routines to MPI-communicate buffer sizes
!
!******************************************************************************







