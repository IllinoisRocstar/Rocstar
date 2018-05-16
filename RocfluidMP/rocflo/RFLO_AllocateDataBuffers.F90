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
!          communication (mixture only).
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = region number
!
! Output: regions%levels%patches%mixt%... = send & receive buffers
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_AllocateDataBuffers.F90,v 1.4 2008/12/06 08:44:25 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_AllocateDataBuffers( regions,iReg )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region, t_dCellTransf
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

! ... loop variables
  INTEGER :: iPatch, iLev, ir

! ... local variables
  INTEGER :: bcType, iRegSrc, n1, n2, n1Src, n2Src, nEqs, nEqsSrc, &
             ndc, ndcSrc, ndim, ndimSrc, errorFlag

  TYPE(t_patch), POINTER       :: patch
  TYPE(t_global), POINTER      :: global
  TYPE(t_dCellTransf), POINTER :: sendEcCell, recvEcCell

!******************************************************************************

  global => regions(iReg)%global

  CALL RegisterFunction( global,'RFLO_AllocateDataBuffers',&
  'RFLO_AllocateDataBuffers.F90' )

! data buffers for patches

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
          nEqs    = CV_MIXT_NEQS      ! valid for mixture
          nEqsSrc = CV_MIXT_NEQS      ! valid for mixture
          ndc     = regions(iReg   )%nDumCells
          ndcSrc  = regions(iRegSrc)%nDumCells
          ndim    = n1*n2*nEqs*ndc
          ndimSrc = n1Src*n2Src*nEqsSrc*ndcSrc
          ALLOCATE( patch%mixt%sendBuff(ndimSrc),stat=errorFlag )
          ALLOCATE( patch%mixt%recvBuff(ndim   ),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
          patch%mixt%nSendBuff = ndimSrc
          patch%mixt%nRecvBuff = ndim
          global%nRequests        = global%nRequests + 1
          patch%mixt%iRequest  = global%nRequests
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
      DO ir=1,global%nRegions
        sendEcCell => regions(iReg)%levels(iLev)%sendEcCells(ir)
        recvEcCell => regions(iReg)%levels(iLev)%recvEcCells(ir)
        IF (sendEcCell%nCells > 0) THEN
          global%nRequests    = global%nRequests + 1
          sendEcCell%iRequest = global%nRequests
          ALLOCATE( sendEcCell%buff(sendEcCell%nCells*CV_MIXT_NEQS), &
                    stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
        ENDIF
        IF (recvEcCell%nCells > 0) THEN
          recvEcCell%iRequest = -999999
          ALLOCATE( recvEcCell%buff(recvEcCell%nCells*CV_MIXT_NEQS), &
                    stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
        ENDIF
      ENDDO    ! ir
    ENDDO      ! iLev

  ENDIF

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_AllocateDataBuffers

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_AllocateDataBuffers.F90,v $
! Revision 1.4  2008/12/06 08:44:25  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:36  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:39:25  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.12  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.11  2003/02/14 22:32:36  jblazek
! Finished implementation of corener and edge cells.
!
! Revision 1.10  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.9  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.8  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.7  2002/03/29 23:15:22  jblazek
! Corrected bug in MPI send.
!
! Revision 1.6  2002/03/18 23:11:32  jblazek
! Finished multiblock and MPI.
!
! Revision 1.5  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.3  2002/01/08 22:09:17  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.2  2001/12/22 00:09:39  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1  2001/12/11 21:59:29  jblazek
! memory allocation added.
!
!******************************************************************************







