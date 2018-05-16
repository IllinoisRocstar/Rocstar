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
! Purpose: allocate data buffers (send & receive) for inter-region
!          communication for PEUL module.
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = region number
!
! Output: regions%levels%patches%bufferPeul%... = send & receive buffers
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_AllocateDataBuffers.F90,v 1.3 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_AllocateDataBuffers( regions,iReg ) ! PUBLIC

  USE ModDataTypes
  USE ModPartEul,    ONLY : t_peul, t_buffer_peul
  USE ModBndPatch,   ONLY : t_patch
  USE ModDataStruct, ONLY : t_region, t_level, t_dCellTransf
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER, INTENT(IN)     :: iReg

! ... loop variables
  INTEGER :: iPatch, iLev, ir

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: bcType, iRegSrc, n1, n2, n1Src, n2Src, nCv, nEqs, nEqsSrc
  INTEGER :: ndc, ndcSrc, ndim, ndimSrc, errorFlag

  TYPE(t_global),      POINTER :: global
  TYPE(t_level),       POINTER :: level
  TYPE(t_patch),       POINTER :: pPatch
  TYPE(t_buffer_peul), POINTER :: pBuffPeul
  TYPE(t_peul),        POINTER :: pPeul
  TYPE(t_dCellTransf), POINTER :: sendEcCell,    recvEcCell
  TYPE(t_dCellTransf), POINTER :: sndPeulEcCell, rcvPeulEcCell

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_AllocateDataBuffers.F90,v $ $Revision: 1.3 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'RFLO_AllocateDataBuffers',&
  'PEUL_AllocateDataBuffers.F90' )

! data buffers for patches

  DO iLev=1,regions(iReg)%nGridLevels

    pPeul => regions(iReg)%levels(iLev)%peul

    pPeul%nRequests = 0
    nCv  =  pPeul%nCv

    DO iPatch=1,regions(iReg)%nPatches

      pPatch => regions(iReg)%levels(iLev)%patches(iPatch)
      bcType = pPatch%bcType

      pBuffPeul => pPatch%bufferPeul

      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
          (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
          (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
        iRegSrc = pPatch%srcRegion
        IF (regions(iRegSrc)%procid /= global%myProcid) THEN  ! other processor
          n1      = ABS(pPatch%l1end   -pPatch%l1beg   ) + 2    ! large enough
          n2      = ABS(pPatch%l2end   -pPatch%l2beg   ) + 2    ! for NODES!
          n1Src   = ABS(pPatch%srcL1end-pPatch%srcL1beg) + 2
          n2Src   = ABS(pPatch%srcL2end-pPatch%srcL2beg) + 2
          nEqs    = nCv
          nEqsSrc = nCv
          ndc     = regions(iReg   )%nDumCells
          ndcSrc  = regions(iRegSrc)%nDumCells
          ndim    = n1*n2*nEqs*ndc
          ndimSrc = n1Src*n2Src*nEqsSrc*ndcSrc

          ALLOCATE( pBuffPeul%sendBuff(ndimSrc),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) &
            CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

          ALLOCATE( pBuffPeul%recvBuff(ndim   ),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= ERR_NONE) &
            CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

          pBuffPeul%nSendBuff = ndimSrc
          pBuffPeul%nRecvBuff = ndim
          pPeul%nRequests     = pPeul%nRequests + 1
          pBuffPeul%iRequest  = pPeul%nRequests
        ENDIF
      ELSE IF ((bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
               (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE)) THEN
        CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__ )  ! #### TEMPORARY ####
      ENDIF  ! bcType

    ENDDO    ! iPatch

! Allocate array for send requests --------------------------------------------

#ifdef MPI
    ALLOCATE( pPeul%requests(pPeul%nRequests),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN
      CALL ErrorStop( global, ERR_ALLOCATE,__LINE__,'pPeul%requests' )
    END IF ! global%error
#endif

  ENDDO      ! iLev

! data buffers for edge & corner cells

  IF (global%nProcAlloc > 1) THEN    ! only if multiple processors

    DO iLev=1,regions(iReg)%nGridLevels

! --- allocate send and receive region data

      level => regions(iReg)%levels(iLev)
      nCv   =  level%peul%nCv

      ALLOCATE( level%sndPeulEcCells(global%nRegions),stat=errorFlag )
      ALLOCATE( level%rcvPeulEcCells(global%nRegions),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      DO ir=1,global%nRegions
        sendEcCell    => regions(iReg)%levels(iLev)%sendEcCells(ir)
        recvEcCell    => regions(iReg)%levels(iLev)%recvEcCells(ir)
        sndPeulEcCell => regions(iReg)%levels(iLev)%sndPeulEcCells(ir)
        rcvPeulEcCell => regions(iReg)%levels(iLev)%rcvPeulEcCells(ir)
        sndPeulEcCell%nCells = sendEcCell%nCells
        rcvPeulEcCell%nCells = recvEcCell%nCells

        IF (sndPeulEcCell%nCells > 0) THEN
          global%nRequests       = global%nRequests + 1
          sndPeulEcCell%iRequest = global%nRequests
          ALLOCATE( sndPeulEcCell%buff(sndPeulEcCell%nCells*nCv), &
                    stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
        ENDIF
        IF (rcvPeulEcCell%nCells > 0) THEN
          rcvPeulEcCell%iRequest = -999999
          ALLOCATE( rcvPeulEcCell%buff(rcvPeulEcCell%nCells*nCv), &
                    stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
        ENDIF
      ENDDO    ! ir
    ENDDO      ! iLev

  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_AllocateDataBuffers

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_AllocateDataBuffers.F90,v $
! Revision 1.3  2008/12/06 08:44:38  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:12  haselbac
! Initial revision after changing case
!
! Revision 1.5  2004/04/15 16:04:03  jferry
! minor formatting (removed trailing spaces)
!
! Revision 1.4  2004/03/02 21:44:07  jferry
! Added corner and edge cell data structures and routines
!
! Revision 1.3  2003/05/05 21:58:29  fnajjar
! Moved nRequests and nCv outside iPatch loop
!
! Revision 1.2  2003/05/05 21:49:50  fnajjar
! Moved pPeul pointer outside iPatch loop
!
! Revision 1.1  2003/04/09 14:32:06  fnajjar
! Initial Import of MPI-based rocsmoke
!
!******************************************************************************







