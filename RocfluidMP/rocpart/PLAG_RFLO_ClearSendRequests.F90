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
!          with the current region.
!
! Description: none.
!
! Input: regions  = all regiona
!        iReg     = current region
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_RFLO_ClearSendRequests.F90,v 1.3 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLO_ClearSendRequests( regions,iReg )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER, INTENT(IN)     :: iReg

! ... loop variables
  INTEGER :: ir

! ... local variables
  INTEGER :: iLev, iRegSrc, iRequest
#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif

  LOGICAL doWait

  TYPE(t_plag),   POINTER :: pPlag
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'PLAG_RFLO_ClearSendRequests',&
  'PLAG_RFLO_ClearSendRequests.F90' )

#ifdef MPI

! get dimension & set pointer -------------------------------------------------

  iLev = 1
  pPlag => regions(iReg)%levels(iLev)%plag

! wait for edges & corners being received by other processors -----------------

  IF ( global%nProcAlloc>1 ) THEN
    DO ir=1,global%nRegions
      IF (regions(iReg)%levels(iLev)%sendEcCells(ir)%nCells > 0) THEN
        iRequest = regions(iReg)%levels(iLev)%sendEcCells(ir)%iRequestMetrics
        CALL MPI_Wait( pPlag%requestsMetrics(iRequest),status,global%mpierr )
        IF (global%mpierr /= ERR_NONE) &
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
    ENDDO
  ENDIF
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLO_ClearSendRequests

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLO_ClearSendRequests.F90,v $
! Revision 1.3  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:58:09  fnajjar
! Initial revision after changing case
!
! Revision 1.1  2004/01/15 21:16:48  fnajjar
! Initial import for corner-edge cell metrics
!
!******************************************************************************







