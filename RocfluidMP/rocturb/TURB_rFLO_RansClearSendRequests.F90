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
! Purpose: wait until RaNS data is received by other processors communicating
!          with the current region
!
! Description: none.
!
! Input: regions  = all regions
!        iReg     = current region
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_rFLO_RansClearSendRequests.F90,v 1.4 2008/12/06 08:44:45 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RFLO_RansClearSendRequests( regions,iReg ) ! PUBLIC

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE TURB_ModInterfaces, ONLY : TURB_RansSAGetEddyVis
  USE ModError
  USE ModMPI
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER, INTENT(IN)     :: iReg

! ... loop variables
  INTEGER :: iPatch, ir

! ... local variables
  INTEGER :: iLev, nPatches, bcType, iRegSrc, iRequest
#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif

  LOGICAL doWait

  TYPE(t_patch),  POINTER :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'TURB_RFLO_RansClearSendRequests',&
  'TURB_rFLO_RansClearSendRequests.F90' )

  IF (regions(iReg)%turbInput%modelClass /= MODEL_RANS) GOTO 999

#ifdef MPI
! get dimensions --------------------------------------------------------------

  iLev     = regions(iReg)%currLevel
  nPatches = regions(iReg)%nPatches
  
! wait for patch data being received by other processors ----------------------

  DO iPatch=1,nPatches

    patch => regions(iReg)%levels(iLev)%patches(iPatch)

    bcType   = patch%bcType
    iRegSrc  = patch%srcRegion
    iRequest = patch%turb%iRequest
    
! - region interface, periodic boundary

    doWait = ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
              (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
              (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
              (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
              (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE))

    IF (iRegSrc > 0) THEN
      IF (doWait .AND. (regions(iRegSrc)%procid /= global%myProcid)) THEN
        CALL MPI_Wait( global%requests(iRequest), status, global%mpierr )
        IF (global%mpierr /= ERR_NONE) &
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
    ENDIF

  ENDDO   ! iPatch

! wait for edges & corners being received by other processors -----------------

  IF (global%nProcAlloc>1) THEN
    DO ir=1,global%nRegions
      IF (regions(iReg)%levels(iLev)%sndTurbEcCells(ir)%nCells > 0) THEN
        iRequest = regions(iReg)%levels(iLev)%sndTurbEcCells(ir)%iRequest
        CALL MPI_Wait( global%requests(iRequest),status,global%mpierr )
        IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF
    ENDDO
  ENDIF
#endif

! get mu_t and kappa_t from tilde[nu]

  CALL TURB_RansSAGetEddyVis( regions(iReg) )

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_RFLO_RansClearSendRequests

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_rFLO_RansClearSendRequests.F90,v $
! Revision 1.4  2008/12/06 08:44:45  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:40:54  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/03/11 03:26:33  wasistho
! changed rocturb nomenclature
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.2  2004/01/23 00:33:41  wasistho
! added clear-send-request of RaNS edge/corner variables
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







