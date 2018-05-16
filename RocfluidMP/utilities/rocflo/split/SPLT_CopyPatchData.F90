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
! Purpose: copy topology data associated with old patch onto a new patch
!          created by splitting the domain.
!
! Description: none.
!
! Input: regionNew = new region
!        patchOld  = patch of old grid
!        l1Off     = offset of the old patch in l1-direction
!        l2Off     = offset of the old patch in l2-direction
!        l1Cells   = max. number of cells in l1-direction
!        l2Cells   = max. number of cells in l2-direction.
!
! Output: regionNew%levels(1)%patch = new topology.
!
! Notes: routine can handle periodic boundaries or interior cuts, but NOT
!        boundaries between regions.
!
!******************************************************************************
!
! $Id: SPLT_CopyPatchData.F90,v 1.3 2008/12/06 08:44:51 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE CopyPatchData( regionsNew,patchOld,splitDirection,iReg, &
                          l1Off,l2Off,l1Cells,l2Cells )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regionsNew(:)
  TYPE(t_patch), POINTER  :: patchOld

  INTEGER :: splitDirection, iReg, l1Off, l2Off, l1Cells, l2Cells

! ... loop variables
  INTEGER :: ir

! ... local variables
  INTEGER :: bcType, lbound, lboundSrc
  INTEGER :: ibeg, jbeg, kbeg, l1OffSrc, l2OffSrc, l1CellsSrc, l2CellsSrc

  TYPE(t_patch), POINTER :: patchNew

!******************************************************************************

  patchNew => regionsNew(iReg)%levels(1)%patches(regionsNew(iReg)%nPatches)

  patchNew%bcType    = patchOld%bcType
  patchNew%lbound    = patchOld%lbound
  patchNew%bcCoupled = patchOld%bcCoupled
  patchNew%l1beg     = MAX(patchOld%l1beg-l1Off+1,1)
  patchNew%l1end     = MIN(patchOld%l1end-l1Off+1,l1Cells)
  patchNew%l2beg     = MAX(patchOld%l2beg-l2Off+1,1)
  patchNew%l2end     = MIN(patchOld%l2end-l2Off+1,l2Cells)

! zero out data of source patch

  patchNew%align     = patchOld%align
  patchNew%srcLbound = 0
  patchNew%srcRegion = 0
  patchNew%srcL1beg  = 0
  patchNew%srcL1end  = 0
  patchNew%srcL2beg  = 0
  patchNew%srcL2end  = 0

  bcType = patchNew%bcType

! source patch is periodic boundary -------------------------------------------

  IF ((bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) .OR. &
      (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE)) THEN

    patchNew%srcLbound = patchOld%srcLbound
    lboundSrc          = patchOld%srcLbound

! - split between the boundaries

    IF (((lboundSrc==1 .OR. lboundSrc==2) .AND. splitDirection==1) .OR. &
        ((lboundSrc==3 .OR. lboundSrc==4) .AND. splitDirection==2) .OR. &
        ((lboundSrc==5 .OR. lboundSrc==6) .AND. splitDirection==3)) THEN
      IF (lboundSrc==1 .OR. lboundSrc==3 .OR. lboundSrc==5) THEN
        patchNew%srcRegion = 1
      ELSE
        patchNew%srcRegion = regionsNew(iReg)%global%nRegions
      ENDIF
      patchNew%srcL1beg = patchOld%srcL1beg
      patchNew%srcL1end = patchOld%srcL1end
      patchNew%srcL2beg = patchOld%srcL2beg
      patchNew%srcL2end = patchOld%srcL2end

! - split along the boundaries

    ELSE
      patchNew%srcRegion = iReg
      patchNew%srcL1beg  = MAX(patchOld%srcL1beg-l1Off+1,1)
      patchNew%srcL1end  = MIN(patchOld%srcL1end-l1Off+1,l1Cells)
      patchNew%srcL2beg  = MAX(patchOld%srcL2beg-l2Off+1,1)
      patchNew%srcL2end  = MIN(patchOld%srcL2end-l2Off+1,l2Cells)
    ENDIF

! source patch is interior cut ------------------------------------------------

  ELSE IF (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) THEN

    patchNew%srcLbound = patchOld%srcLbound
    lbound             = patchOld%lbound
    lboundSrc          = patchOld%srcLbound

! - source patch on the same boundary (C-grid)

    IF (lbound == lboundSrc) THEN

! --- boundary cannot have been splitted

      IF (((lbound==1 .OR. lbound==2) .AND. splitDirection==1) .OR. &
          ((lbound==3 .OR. lbound==4) .AND. splitDirection==2) .OR. &
          ((lbound==5 .OR. lbound==6) .AND. splitDirection==3)) THEN
        patchNew%srcRegion = iReg
        patchNew%srcL1beg  = patchOld%srcL1beg
        patchNew%srcL1end  = patchOld%srcL1end
        patchNew%srcL2beg  = patchOld%srcL2beg
        patchNew%srcL2end  = patchOld%srcL2end

! --- boundary possibly splitted

      ELSE
        patchNew%srcRegion = regionsNew(iReg)%global%nRegions - iReg + 1
        IF (patchOld%srcL1beg > patchOld%srcL1End) THEN    ! l1 reversed
          IF (((lbound==1 .OR. lbound==2) .AND. splitDirection==3) .OR. &
              ((lbound==3 .OR. lbound==4) .AND. splitDirection==1) .OR. &
              ((lbound==5 .OR. lbound==6) .AND. splitDirection==2)) THEN
            patchNew%srcRegion = iReg
          ENDIF
        ELSE                                               ! l2 reversed
          IF (((lbound==1 .OR. lbound==2) .AND. splitDirection==2) .OR. &
              ((lbound==3 .OR. lbound==4) .AND. splitDirection==3) .OR. &
              ((lbound==5 .OR. lbound==6) .AND. splitDirection==1)) THEN
            patchNew%srcRegion = iReg
          ENDIF
        ENDIF
        ibeg = 1
        jbeg = 1
        kbeg = 1
        DO ir=1,patchNew%srcRegion-1
          IF (splitDirection == 1) THEN
            ibeg = ibeg + regionsNew(ir)%levels(1)%grid%ipc
          ELSE IF (splitDirection == 2) THEN
            jbeg = jbeg + regionsNew(ir)%levels(1)%grid%jpc
          ELSE
            kbeg = kbeg + regionsNew(ir)%levels(1)%grid%kpc
          ENDIF
        ENDDO
        IF (lbound==1 .OR. lbound==2) THEN
          l1OffSrc   = jbeg
          l2OffSrc   = kbeg
          l1CellsSrc = regionsNew(patchNew%srcRegion)%levels(1)%grid%jpc
          l2CellsSrc = regionsNew(patchNew%srcRegion)%levels(1)%grid%kpc
        ELSE IF (lbound==3 .OR. lbound==4) THEN
          l1OffSrc   = kbeg
          l2OffSrc   = ibeg
          l1CellsSrc = regionsNew(patchNew%srcRegion)%levels(1)%grid%kpc
          l2CellsSrc = regionsNew(patchNew%srcRegion)%levels(1)%grid%ipc
        ELSE IF (lbound==5 .OR. lbound==6) THEN
          l1OffSrc   = ibeg
          l2OffSrc   = jbeg
          l1CellsSrc = regionsNew(patchNew%srcRegion)%levels(1)%grid%ipc
          l2CellsSrc = regionsNew(patchNew%srcRegion)%levels(1)%grid%jpc
        ENDIF
        patchNew%srcL1beg = MAX(patchOld%srcL1beg-l1OffSrc+1,1)
        patchNew%srcL1beg = MIN(patchNew%srcL1beg,l1CellsSrc)
        patchNew%srcL1end = MIN(patchOld%srcL1end-l1OffSrc+1,l1CellsSrc)
        patchNew%srcL1end = MAX(patchNew%srcL1end,1)
        patchNew%srcL2beg = MAX(patchOld%srcL2beg-l2OffSrc+1,1)
        patchNew%srcL2beg = MIN(patchNew%srcL2beg,l2CellsSrc)
        patchNew%srcL2end = MIN(patchOld%srcL2end-l2OffSrc+1,l2CellsSrc)
        patchNew%srcL2end = MAX(patchNew%srcL2end,1)
      ENDIF

! - source patch on different boundary (O-grid)

    ELSE

! --- boundary cannot have been splitted

      IF (((lbound==1 .OR. lbound==2) .AND. splitDirection==1) .OR. &
          ((lbound==3 .OR. lbound==4) .AND. splitDirection==2) .OR. &
          ((lbound==5 .OR. lbound==6) .AND. splitDirection==3)) THEN
        IF (lboundSrc==1 .OR. lboundSrc==3 .OR. lboundSrc==5) THEN
          patchNew%srcRegion = 1
        ELSE
          patchNew%srcRegion = regionsNew(iReg)%global%nRegions
        ENDIF
        patchNew%srcL1beg  = patchOld%srcL1beg
        patchNew%srcL1end  = patchOld%srcL1end
        patchNew%srcL2beg  = patchOld%srcL2beg
        patchNew%srcL2end  = patchOld%srcL2end

! --- boundary splitted

      ELSE
        patchNew%srcRegion = iReg
        patchNew%srcL1beg  = MAX(patchOld%srcL1beg-l1Off+1,1)
        patchNew%srcL1end  = MIN(patchOld%srcL1end-l1Off+1,l1Cells)
        patchNew%srcL2beg  = MAX(patchOld%srcL2beg-l2Off+1,1)
        patchNew%srcL2end  = MIN(patchOld%srcL2end-l2Off+1,l2Cells)
      ENDIF

    ENDIF  ! lbound == lboundSrc

  ENDIF    ! bcType

END SUBROUTINE CopyPatchData

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SPLT_CopyPatchData.F90,v $
! Revision 1.3  2008/12/06 08:44:51  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:01  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/03 02:41:15  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:45:30  wasistho
! lower to upper case
!
! Revision 1.4  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.3  2002/09/27 00:57:11  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.2  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/07/12 21:50:07  jblazek
! Added tool to split single grid into multiple regions.
!
!******************************************************************************






