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
!        geometry = true if geometrical data were send, false otherwise.
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ClearSendRequests.F90,v 1.4 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ClearSendRequests( regions,iReg,geometry )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

  LOGICAL :: geometry

! ... loop variables
  INTEGER :: iPatch, ir

! ... local variables
  INTEGER :: iLev, nPatches, bcType, iRegSrc, iRequest
#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif

  LOGICAL doWait

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ClearSendRequests',&
  'RFLO_ClearSendRequests.F90' )

#ifdef MPI
! get dimensions --------------------------------------------------------------

  IF (geometry) THEN
    iLev = 1
  ELSE
    iLev = regions(iReg)%currLevel
  ENDIF
  nPatches = regions(iReg)%nPatches

! wait for patch data being received by other processors ----------------------

  DO iPatch=1,nPatches

    patch => regions(iReg)%levels(iLev)%patches(iPatch)

    bcType  = patch%bcType
    iRegSrc = patch%srcRegion

! - region interface, periodic boundary

    IF (geometry) THEN
      doWait = ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) &
                .OR. &
                (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) &
                .OR. &
                (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE))
    ELSE
      doWait = ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) &
                .OR. &
                (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) &
                .OR. &
                (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) &
                .OR. &
                (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) &
                .OR. &
                (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE))
    ENDIF

    IF (iRegSrc > 0) THEN
      IF (doWait .AND. (regions(iRegSrc)%procid /= global%myProcid)) THEN
        CALL MPI_Wait( global%requests(patch%mixt%iRequest),status, &
                       global%mpierr )
        IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
    ENDIF

  ENDDO   ! iPatch

! wait for edges & corners being received by other processors -----------------

  IF (.NOT. geometry .AND. global%nProcAlloc>1) THEN
    DO ir=1,global%nRegions
      IF (regions(iReg)%levels(iLev)%sendEcCells(ir)%nCells > 0) THEN
        iRequest = regions(iReg)%levels(iLev)%sendEcCells(ir)%iRequest
        CALL MPI_Wait( global%requests(iRequest),status,global%mpierr )
        IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
    ENDDO
  ENDIF
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ClearSendRequests

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ClearSendRequests.F90,v $
! Revision 1.4  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:39:32  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.8  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.7  2003/02/24 16:57:43  jferry
! Moved #endif for #ifdef MPI to encompass final section of code
!
! Revision 1.6  2003/02/14 22:32:37  jblazek
! Finished implementation of corener and edge cells.
!
! Revision 1.5  2002/09/30 19:58:17  jblazek
! Added check of iRegSrc before IF (dowait ...) THEN statement.
!
! Revision 1.4  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.3  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/07/12 21:50:07  jblazek
! Added tool to split single grid into multiple regions.
!
! Revision 1.1  2002/04/01 19:36:08  jblazek
! Added routine to clear send requests.
!
!******************************************************************************







