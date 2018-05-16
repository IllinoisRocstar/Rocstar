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
! Purpose: wait until data is received by other processors communicating
!          with the current region specific for PEUL Module.
!
! Description: none.
!
! Input: regions  = all regiona
!        iReg     = current region
!        geometry = true if geometrical data were send, false otherwise.
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_ClearSendRequests.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_ClearSendRequests( regions,iReg ) ! PUBLIC

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartEul, ONLY    : t_peul
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER    :: regions(:)
  INTEGER,        INTENT(IN) :: iReg

! ... loop variables
  INTEGER :: iPatch, ir

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iLev, nPatches, bcType, iRegSrc, iRequest
#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif

  LOGICAL doWait

  TYPE(t_patch),  POINTER :: patch
  TYPE(t_peul),   POINTER :: pPeul
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_ClearSendRequests.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PEUL_ClearSendRequests',&
  'PEUL_ClearSendRequests.F90' )

#ifdef MPI
! get dimensions --------------------------------------------------------------

  iLev     = regions(iReg)%currLevel
  nPatches = regions(iReg)%nPatches

  pPeul => regions(iReg)%levels(iLev)%peul

! wait for patch data being received by other processors ----------------------

  DO iPatch=1,nPatches

    patch => regions(iReg)%levels(iLev)%patches(iPatch)

    bcType   = patch%bcType
    iRegSrc  = patch%srcRegion
    iRequest = patch%bufferPeul%iRequest

! - region interface, periodic boundary

    doWait = ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
              (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
              (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE))

    IF (iRegSrc > 0) THEN
      IF (doWait .AND. (regions(iRegSrc)%procid /= global%myProcid)) THEN
        CALL MPI_Wait( pPeul%requests(iRequest), status, global%mpierr )
        IF (global%mpierr /= ERR_NONE) &
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
    ENDIF

  ENDDO   ! iPatch

! wait for edges & corners being received by other processors -----------------

  IF (global%nProcAlloc>1) THEN
    DO ir=1,global%nRegions
      IF (regions(iReg)%levels(iLev)%sndPeulEcCells(ir)%nCells > 0) THEN
        iRequest = regions(iReg)%levels(iLev)%sndPeulEcCells(ir)%iRequest
        CALL MPI_Wait( global%requests(iRequest),status,global%mpierr )
        IF (global%mpierr /= 0) &
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
    ENDDO
  ENDIF

#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_ClearSendRequests

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_ClearSendRequests.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:28  haselbac
! Initial revision after changing case
!
! Revision 1.3  2004/04/15 16:04:03  jferry
! minor formatting (removed trailing spaces)
!
! Revision 1.2  2004/03/02 21:44:07  jferry
! Added corner and edge cell data structures and routines
!
! Revision 1.1  2003/04/09 14:33:38  fnajjar
! Initial Import of MPI-based rocsmoke
!
!******************************************************************************







