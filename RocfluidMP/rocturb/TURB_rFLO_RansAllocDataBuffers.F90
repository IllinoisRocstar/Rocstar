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
! Purpose: allocate RaNS data buffers (send & receive) for inter-region
!          communication (mixture only).
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = region number
!
! Output: regions%levels%patches%turb%... = send & receive buffers
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_rFLO_RansAllocDataBuffers.F90,v 1.4 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RFLO_RansAllocDataBuffers( regions,iReg ) ! PUBLIC

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region, t_level, t_dCellTransf
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

! ... loop variables
  INTEGER :: iPatch, iLev, ir

! ... local variables
  INTEGER :: bcType, iRegSrc, n1, n2, n1Src, n2Src, nEqs, nEqsSrc, &
             ndc, ndcSrc, ndim, ndimSrc, nCv, errorFlag

  TYPE(t_patch), POINTER       :: patch
  TYPE(t_global), POINTER      :: global
  TYPE(t_level), POINTER       :: level
  TYPE(t_dCellTransf), POINTER :: sendEcCell, recvEcCell
  TYPE(t_dCellTransf), POINTER :: sndTurbEcCell, rcvTurbEcCell

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'TURB_RFLO_RansAllocDataBuffers',&
  'TURB_rFLO_RansAllocDataBuffers.F90' )

  IF (regions(iReg)%turbInput%modelClass /= MODEL_RANS) GOTO 999

! data buffers for patches ----------------------------------------------------

  nCv = regions(iReg)%turbInput%nCv

  DO iLev=1,regions(iReg)%nGridLevels
    DO iPatch=1,regions(iReg)%nPatches
  
      patch => regions(iReg)%levels(iLev)%patches(iPatch)
      bcType = patch%bcType

      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
          (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
          (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
        iRegSrc = patch%srcRegion
        IF (regions(iRegSrc)%procid /= global%myProcid) THEN  ! other processor
          n1      = ABS(patch%l1end   -patch%l1beg   ) + 2    ! large enough
          n2      = ABS(patch%l2end   -patch%l2beg   ) + 2    ! for NODES!
          n1Src   = ABS(patch%srcL1end-patch%srcL1beg) + 2
          n2Src   = ABS(patch%srcL2end-patch%srcL2beg) + 2
          nEqs    = nCv
          nEqsSrc = nCv
          ndc     = regions(iReg   )%nDumCells
          ndcSrc  = regions(iRegSrc)%nDumCells
          ndim    = n1*n2*nEqs*ndc
          ndimSrc = n1Src*n2Src*nEqsSrc*ndcSrc
          ALLOCATE( patch%turb%sendBuff(ndimSrc),stat=errorFlag )
          ALLOCATE( patch%turb%recvBuff(ndim   ),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
          patch%turb%nSendBuff = ndimSrc
          patch%turb%nRecvBuff = ndim
          global%nRequests        = global%nRequests + 1
          patch%turb%iRequest  = global%nRequests
        ENDIF
      ELSE IF ((bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
               (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE)) THEN
        CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__ )  ! #### TEMPORARY ####
      ENDIF  ! bcType

    ENDDO    ! iPatch
  ENDDO      ! iLev

! data buffers for edge & corner cells

  IF (global%nProcAlloc > 1) THEN    ! only if multiple processors

    DO iLev=1,regions(iReg)%nGridLevels

! --- allocate send and receive region data

      level => regions(iReg)%levels(iLev)

      ALLOCATE( level%sndTurbEcCells(global%nRegions),stat=errorFlag )
      ALLOCATE( level%rcvTurbEcCells(global%nRegions),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      DO ir=1,global%nRegions
        sendEcCell    => regions(iReg)%levels(iLev)%sendEcCells(ir)
        recvEcCell    => regions(iReg)%levels(iLev)%recvEcCells(ir)
        sndTurbEcCell => regions(iReg)%levels(iLev)%sndTurbEcCells(ir)
        rcvTurbEcCell => regions(iReg)%levels(iLev)%rcvTurbEcCells(ir)
        sndTurbEcCell%nCells = sendEcCell%nCells
        rcvTurbEcCell%nCells = recvEcCell%nCells

        IF (sndTurbEcCell%nCells > 0) THEN
          global%nRequests       = global%nRequests + 1
          sndTurbEcCell%iRequest = global%nRequests
          ALLOCATE( sndTurbEcCell%buff(sndTurbEcCell%nCells*nCv), &
                    stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
        ENDIF
        IF (rcvTurbEcCell%nCells > 0) THEN
          rcvTurbEcCell%iRequest = -999999
          ALLOCATE( rcvTurbEcCell%buff(rcvTurbEcCell%nCells*nCv), &
                    stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
        ENDIF
      ENDDO    ! ir
    ENDDO      ! iLev

  ENDIF

! finalize

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_RFLO_RansAllocDataBuffers

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_rFLO_RansAllocDataBuffers.F90,v $
! Revision 1.4  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:40:53  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/03/11 03:26:33  wasistho
! changed rocturb nomenclature
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.2  2004/01/23 00:27:58  wasistho
! add data allocations for RaNS variables at edge and corners
!
! Revision 1.1  2003/10/07 02:17:02  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







