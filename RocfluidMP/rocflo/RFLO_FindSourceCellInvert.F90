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
! Purpose: given the indices of a dummy cell, find the source region
!          and the source cell by transformation between patches.
!
! Description: the only difference with RFLO_findSourceCell is that here
!              the loop over region patches is inverted from the highest to
!              the lowest patch number.
!
! Input: regions  = dimensions and topology of all regions
!        iReg     = current region
!        iLev     = current grid level
!        i/j/kc   = indices of the dummy cell
!
! Output: found   = flag if a source cell was found and is within the physical
!                   domain of the source region (true/false)
!         rotate  = rotational periodicity involved (true/false)
!         i/j/kc  = indices of the source cell
!         icell   = index of the source cell
!         iRegSrc = index of the source region
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_FindSourceCellInvert.F90,v 1.3 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_FindSourceCellInvert( regions,iReg,iLev,ic,jc,kc,icell, &
                                      found,rotate,iRegSrc )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetPatchIndices, &
                            RFLO_SourceCell
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  INTEGER :: iReg, iLev, ic, jc, kc, icell, iRegSrc

  LOGICAL :: found, rotate

  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: ilb, iPatch

! ... local variables
  LOGICAL :: hit, debug

  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend, iPatchSrc
  INTEGER :: bcType, lbound, ibeg, iend, jbeg, jend, kbeg, kend

  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch, patchSrc

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_FindSourceCellInvert',&
  'RFLO_FindSourceCellInvert.F90' )

  found = .false.

! get dimensions --------------------------------------------------------------

  CALL RFLO_GetDimensPhys( regions(iReg),iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )

! find a suitable patch to start from -----------------------------------------

  hit = .false.     ! no patch found yet

! cell within the patch?

  DO iPatch=regions(iReg)%nPatches,1,-1
    patch => regions(iReg)%levels(iLev)%patches(iPatch)
    bcType = patch%bcType
    lbound = patch%lbound
    IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
        (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
        (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
      CALL RFLO_GetPatchIndices( regions(iReg),patch,iLev, &
                                 ibeg,iend,jbeg,jend,kbeg,kend )
      IF      ((lbound==1 .AND. ic<ipcbeg) .OR. &
               (lbound==2 .AND. ic>ipcend)) THEN
        IF ((jc>=jbeg .AND. jc<=jend) .AND. &
            (kc>=kbeg .AND. kc<=kend)) THEN
          hit = .true.
          EXIT
        ENDIF
      ELSE IF ((lbound==3 .AND. jc<jpcbeg) .OR. &
               (lbound==4 .AND. jc>jpcend)) THEN 
        IF ((ic>=ibeg .AND. ic<=iend) .AND. &
            (kc>=kbeg .AND. kc<=kend)) THEN
          hit = .true.
          EXIT
        ENDIF
      ELSE IF ((lbound==5 .AND. kc<kpcbeg) .OR. &
               (lbound==6 .AND. kc>kpcend)) THEN
        IF ((jc>=jbeg .AND. jc<=jend) .AND. &
            (ic>=ibeg .AND. ic<=iend)) THEN
          hit = .true.
          EXIT
        ENDIF
      ENDIF
    ENDIF  ! bcType
  ENDDO    ! iPatch

! cell just outside the patch?

  IF (.NOT. hit) THEN
    DO iPatch=regions(iReg)%nPatches,1,-1
      patch => regions(iReg)%levels(iLev)%patches(iPatch)
      bcType = patch%bcType
      lbound = patch%lbound
      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
          (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
          (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
        CALL RFLO_GetPatchIndices( regions(iReg),patch,iLev, &
                                   ibeg,iend,jbeg,jend,kbeg,kend )
        IF      ((lbound==1 .AND. ic<ipcbeg) .OR. &
                 (lbound==2 .AND. ic>ipcend)) THEN       ! face 1, 2
          IF (kc<kpcbeg .AND. kbeg==kpcbeg .AND. &
              (jc>=jbeg .AND. jc<=jend)) THEN
            hit = .true.
            EXIT
          ENDIF
          IF (kc>kpcend .AND. kend==kpcend .AND. &
              (jc>=jbeg .AND. jc<=jend)) THEN
            hit = .true.
            EXIT
          ENDIF
          IF (jc<jpcbeg .AND. jbeg==jpcbeg .AND. &
              (kc>=kbeg .AND. kc<=kend)) THEN
            hit = .true.
            EXIT
          ENDIF
          IF (jc>jpcend .AND. jend==jpcend .AND. &
              (kc>=kbeg .AND. kc<=kend)) THEN
            hit = .true.
            EXIT
          ENDIF
        ELSE IF ((lbound==3 .AND. jc<jpcbeg) .OR. &
                 (lbound==4 .AND. jc>jpcend)) THEN       ! face 3, 4
          IF (kc<kpcbeg .AND. kbeg==kpcbeg .AND. &
              (ic>=ibeg .AND. ic<=iend)) THEN
            hit = .true.
            EXIT
          ENDIF
          IF (kc>kpcend .AND. kend==kpcend .AND. &
              (ic>=ibeg .AND. ic<=iend)) THEN
            hit = .true.
            EXIT
          ENDIF
          IF (ic<ipcbeg .AND. ibeg==ipcbeg .AND. &
              (kc>=kbeg .AND. kc<=kend)) THEN
            hit = .true.
            EXIT
          ENDIF
          IF (ic>ipcend .AND. iend==ipcend .AND. &
              (kc>=kbeg .AND. kc<=kend)) THEN
            hit = .true.
            EXIT
          ENDIF
        ELSE IF ((lbound==5 .AND. kc<kpcbeg) .OR. &
                 (lbound==6 .AND. kc>kpcend)) THEN       ! face 5, 6
          IF (jc<jpcbeg .AND. jbeg==jpcbeg .AND. &
              (ic>=ibeg .AND. ic<=iend)) THEN
            hit = .true.
            EXIT
          ENDIF
          IF (jc>jpcend .AND. jend==jpcend .AND. &
              (ic>=ibeg .AND. ic<=iend)) THEN
            hit = .true.
            EXIT
          ENDIF
          IF (ic<ipcbeg .AND. ibeg==ipcbeg .AND. &
              (jc>=jbeg .AND. jc<=jend)) THEN
            hit = .true.
            EXIT
          ENDIF
          IF (ic>ipcend .AND. iend==ipcend .AND. &
              (jc>=jbeg .AND. jc<=jend)) THEN
            hit = .true.
            EXIT
          ENDIF
        ENDIF
      ENDIF  ! bcType
    ENDDO    ! iPatch
  ENDIF      ! .NOT. hit

! cell at some corner?

  IF (.NOT. hit) THEN
    DO iPatch=regions(iReg)%nPatches,1,-1
      patch => regions(iReg)%levels(iLev)%patches(iPatch)
      bcType = patch%bcType
      lbound = patch%lbound
      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
          (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
          (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN
        CALL RFLO_GetPatchIndices( regions(iReg),patch,iLev, &
                                   ibeg,iend,jbeg,jend,kbeg,kend )
        IF      (ic<ipcbeg .AND. jc<jpcbeg .AND. kc<kpcbeg) THEN   ! corner 1
          IF ((lbound==1 .AND. jbeg==jpcbeg .AND. kbeg==kpcbeg) .OR. &
              (lbound==3 .AND. ibeg==ipcbeg .AND. kbeg==kpcbeg) .OR. &
              (lbound==5 .AND. ibeg==ipcbeg .AND. jbeg==jpcbeg)) THEN
            hit = .true.
            EXIT
          ENDIF
        ELSE IF (ic<ipcbeg .AND. jc<jpcbeg .AND. kc>kpcend) THEN   ! corner 2
          IF ((lbound==1 .AND. jbeg==jpcbeg .AND. kend==kpcend) .OR. &
              (lbound==3 .AND. ibeg==ipcbeg .AND. kend==kpcend) .OR. &
              (lbound==6 .AND. ibeg==ipcbeg .AND. jbeg==jpcbeg)) THEN
            hit = .true.
            EXIT
          ENDIF
        ELSE IF (ic<ipcbeg .AND. jc>jpcend .AND. kc>kpcend) THEN   ! corner 3
          IF ((lbound==1 .AND. jend==jpcend .AND. kend==kpcend) .OR. &
              (lbound==4 .AND. ibeg==ipcbeg .AND. kend==kpcend) .OR. &
              (lbound==6 .AND. ibeg==ipcbeg .AND. jend==jpcend)) THEN
            hit = .true.
            EXIT
          ENDIF
        ELSE IF (ic<ipcbeg .AND. jc>jpcend .AND. kc<kpcbeg) THEN   ! corner 4
          IF ((lbound==1 .AND. jend==jpcend .AND. kbeg==kpcbeg) .OR. &
              (lbound==4 .AND. ibeg==ipcbeg .AND. kbeg==kpcbeg) .OR. &
              (lbound==5 .AND. ibeg==ipcbeg .AND. jend==jpcend)) THEN
            hit = .true.
            EXIT
          ENDIF
        ELSE IF (ic>ipcend .AND. jc<jpcbeg .AND. kc<kpcbeg) THEN   ! corner 5
          IF ((lbound==2 .AND. jbeg==jpcbeg .AND. kbeg==kpcbeg) .OR. &
              (lbound==3 .AND. iend==ipcend .AND. kbeg==kpcbeg) .OR. &
              (lbound==5 .AND. iend==ipcend .AND. jbeg==jpcbeg)) THEN
            hit = .true.
            EXIT
          ENDIF
        ELSE IF (ic>ipcend .AND. jc<jpcbeg .AND. kc>kpcend) THEN   ! corner 6
          IF ((lbound==2 .AND. jbeg==jpcbeg .AND. kend==kpcend) .OR. &
              (lbound==3 .AND. iend==ipcend .AND. kend==kpcend) .OR. &
              (lbound==6 .AND. iend==ipcend .AND. jbeg==jpcbeg)) THEN
            hit = .true.
            EXIT
          ENDIF
        ELSE IF (ic>ipcend .AND. jc>jpcend .AND. kc>kpcend) THEN   ! corner 7
          IF ((lbound==2 .AND. jend==jpcend .AND. kend==kpcend) .OR. &
              (lbound==4 .AND. iend==ipcend .AND. kend==kpcend) .OR. &
              (lbound==6 .AND. iend==ipcend .AND. jend==jpcend)) THEN
            hit = .true.
            EXIT
          ENDIF
        ELSE IF (ic>ipcend .AND. jc>jpcend .AND. kc<kpcbeg) THEN   ! corner 8
          IF ((lbound==2 .AND. jend==jpcend .AND. kbeg==kpcbeg) .OR. &
              (lbound==4 .AND. iend==ipcend .AND. kbeg==kpcbeg) .OR. &
              (lbound==5 .AND. iend==ipcend .AND. jend==jpcend)) THEN
            hit = .true.
            EXIT
          ENDIF
        ENDIF
      ENDIF  ! bcType
    ENDDO    ! iPatch
  ENDIF      ! .NOT. hit

! if patch was found, do the transformation -----------------------------------

  IF (hit) THEN
    bcType    = patch%bcType
    iRegSrc   = patch%srcRegion
    iPatchSrc = patch%srcPatch
    patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)
    CALL RFLO_SourceCell( regions(iReg),regions(iRegSrc),patch,patchSrc, &
                          iLev,ic,jc,kc,icell,found )
    IF (found) THEN
      IF (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE) THEN
        rotate = .true.
      ELSE
        rotate = .false.
      ENDIF
    ENDIF
  ELSE
    iRegSrc = iReg
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_FindSourceCellInvert

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_FindSourceCellInvert.F90,v $
! Revision 1.3  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.1  2004/08/21 00:37:16  wasistho
! initial checkin to search for degenerated edge/corners
!
! Revision 1.6  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.2  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.1  2003/02/03 19:20:47  jblazek
! Added treatment of edge and corner cells for one processor.
!
!******************************************************************************







