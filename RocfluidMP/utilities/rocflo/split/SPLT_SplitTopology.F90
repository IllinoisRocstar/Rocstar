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
! Purpose: generate new topology for the splitted grid.
!
! Description: none.
!
! Input: splitDirection = direction (i,j,k) in which the grid is to be splitted
!        regionsOld     = old topology.
!
! Output: regionsNew = new topology.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: SPLT_SplitTopology.F90,v 1.4 2008/12/06 08:44:51 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE SplitTopology( splitDirection,regionsOld,regionsNew )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE SPLT_ModInterfaces, ONLY : CopyPatchData
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  INTEGER :: splitDirection

  TYPE(t_region), POINTER :: regionsOld(:), regionsNew(:)

! ... loop variables
  INTEGER :: iReg, iPatchOld

! ... local variables
  INTEGER :: lbound, l1beg, l1end, l2beg, l2end
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, iPatch

  TYPE(t_global), POINTER :: globalNew
  TYPE(t_grid) , POINTER  :: gridOld, gridNew
  TYPE(t_patch), POINTER  :: patchOld, patchNew

!******************************************************************************

  globalNew => regionsNew(1)%global

  CALL RegisterFunction( globalNew,'SplitTopology',&
  'SPLT_SplitTopology.F90' )

! copy old patches ------------------------------------------------------------

  gridOld => regionsOld(1)%levels(1)%grid

  DO iPatchOld=1,regionsOld(1)%nPatches
    patchOld => regionsOld(1)%levels(1)%patches(iPatchOld)
    lbound   = patchOld%lbound
    l1beg    = patchOld%l1beg
    l1end    = patchOld%l1end
    l2beg    = patchOld%l2beg
    l2end    = patchOld%l2end
    ibeg     = 1                ! offsets of new regions wrp. to old grid
    jbeg     = 1
    kbeg     = 1

! - loop over all new regions

    DO iReg=1,globalNew%nRegions

      gridNew => regionsNew(iReg)%levels(1)%grid
      iend    = ibeg + gridNew%ipc - 1
      jend    = jbeg + gridNew%jpc - 1
      kend    = kbeg + gridNew%kpc - 1

! --- i=const faces

      IF (lbound==1 .AND. ibeg==1) THEN
        IF ((l1beg<=jend .AND. l1end>=jbeg) .AND. &
            (l2beg<=kend .AND. l2end>=kbeg)) THEN
          regionsNew(iReg)%nPatches = regionsNew(iReg)%nPatches + 1
          CALL CopyPatchData( regionsNew,patchOld,splitDirection,iReg, &
                              jbeg,kbeg,gridNew%jpc,gridNew%kpc )
        ENDIF

      ELSE IF (lbound==2 .AND. iend==gridOld%ipc) THEN
        IF ((l1beg<=jend .AND. l1end>=jbeg) .AND. &
            (l2beg<=kend .AND. l2end>=kbeg)) THEN
          regionsNew(iReg)%nPatches = regionsNew(iReg)%nPatches + 1
          CALL CopyPatchData( regionsNew,patchOld,splitDirection,iReg, &
                              jbeg,kbeg,gridNew%jpc,gridNew%kpc )
        ENDIF

! --- j=const. faces

      ELSE IF (lbound==3 .AND. jbeg==1) THEN
        IF ((l1beg<=kend .AND. l1end>=kbeg) .AND. &
            (l2beg<=iend .AND. l2end>=ibeg)) THEN
          regionsNew(iReg)%nPatches = regionsNew(iReg)%nPatches + 1
          CALL CopyPatchData( regionsNew,patchOld,splitDirection,iReg, &
                              kbeg,ibeg,gridNew%kpc,gridNew%ipc )
        ENDIF

      ELSE IF (lbound==4 .AND. jend==gridOld%jpc) THEN
        IF ((l1beg<=kend .AND. l1end>=kbeg) .AND. &
            (l2beg<=iend .AND. l2end>=ibeg)) THEN
          regionsNew(iReg)%nPatches = regionsNew(iReg)%nPatches + 1
          CALL CopyPatchData( regionsNew,patchOld,splitDirection,iReg, &
                              kbeg,ibeg,gridNew%kpc,gridNew%ipc )
        ENDIF

! --- k=const. faces

      ELSE IF (lbound==5 .AND. kbeg==1) THEN
        IF ((l1beg<=iend .AND. l1end>=ibeg) .AND. &
            (l2beg<=jend .AND. l2end>=jbeg)) THEN
          regionsNew(iReg)%nPatches = regionsNew(iReg)%nPatches + 1
          CALL CopyPatchData( regionsNew,patchOld,splitDirection,iReg, &
                              ibeg,jbeg,gridNew%ipc,gridNew%jpc )
        ENDIF

      ELSE IF (lbound==6 .AND. kend==gridOld%kpc) THEN
        IF ((l1beg<=iend .AND. l1end>=ibeg) .AND. &
            (l2beg<=jend .AND. l2end>=jbeg)) THEN
          regionsNew(iReg)%nPatches = regionsNew(iReg)%nPatches + 1
          CALL CopyPatchData( regionsNew,patchOld,splitDirection,iReg, &
                              ibeg,jbeg,gridNew%ipc,gridNew%jpc )
        ENDIF
      ENDIF    ! lbound

! --- increase offsets

      IF (splitDirection == 1) THEN
        ibeg = ibeg + gridNew%ipc
      ELSE IF (splitDirection == 2) THEN
        jbeg = jbeg + gridNew%jpc
      ELSE
        kbeg = kbeg + gridNew%kpc
      ENDIF

    ENDDO  ! iReg

  ENDDO    ! iPatchOld

! generate connections between regions ----------------------------------------
! faces 2, 4, 6

  DO iReg=1,globalNew%nRegions-1
    regionsNew(iReg)%nPatches = regionsNew(iReg)%nPatches + 1
    iPatch             = regionsNew(iReg)%nPatches
    patchNew           => regionsNew(iReg)%levels(1)%patches(iPatch)
    gridNew            => regionsNew(iReg)%levels(1)%grid
    patchNew%bcType    = BC_REGIONCONF
    patchNew%bcCoupled = 0
    patchNew%align     = .true.
    patchNew%srcRegion = iReg + 1

    IF (splitDirection == 1) THEN
      patchNew%lbound    = 2
      patchNew%l1beg     = 1
      patchNew%l1end     = gridNew%jpc
      patchNew%l2beg     = 1
      patchNew%l2end     = gridNew%kpc
      patchNew%srcLbound = 1
      patchNew%srcL1beg  = 1
      patchNew%srcL1end  = gridNew%jpc
      patchNew%srcL2beg  = 1
      patchNew%srcL2end  = gridNew%kpc
    ELSE IF (splitDirection == 2) THEN
      patchNew%lbound    = 4
      patchNew%l1beg     = 1
      patchNew%l1end     = gridNew%kpc
      patchNew%l2beg     = 1
      patchNew%l2end     = gridNew%ipc
      patchNew%srcLbound = 3
      patchNew%srcL1beg  = 1
      patchNew%srcL1end  = gridNew%kpc
      patchNew%srcL2beg  = 1
      patchNew%srcL2end  = gridNew%ipc
    ELSE
      patchNew%lbound    = 6
      patchNew%l1beg     = 1
      patchNew%l1end     = gridNew%ipc
      patchNew%l2beg     = 1
      patchNew%l2end     = gridNew%jpc
      patchNew%srcLbound = 5
      patchNew%srcL1beg  = 1
      patchNew%srcL1end  = gridNew%ipc
      patchNew%srcL2beg  = 1
      patchNew%srcL2end  = gridNew%jpc
    ENDIF
  ENDDO  ! iReg

! faces 1, 3, 5

  DO iReg=2,globalNew%nRegions
    regionsNew(iReg)%nPatches = regionsNew(iReg)%nPatches + 1
    iPatch             = regionsNew(iReg)%nPatches
    patchNew           => regionsNew(iReg)%levels(1)%patches(iPatch)
    gridNew            => regionsNew(iReg)%levels(1)%grid
    patchNew%bcType    = BC_REGIONCONF
    patchNew%bcCoupled = 0
    patchNew%align     = .true.
    patchNew%srcRegion = iReg - 1

    IF (splitDirection == 1) THEN
      patchNew%lbound    = 1
      patchNew%l1beg     = 1
      patchNew%l1end     = gridNew%jpc
      patchNew%l2beg     = 1
      patchNew%l2end     = gridNew%kpc
      patchNew%srcLbound = 2
      patchNew%srcL1beg  = 1
      patchNew%srcL1end  = gridNew%jpc
      patchNew%srcL2beg  = 1
      patchNew%srcL2end  = gridNew%kpc
    ELSE IF (splitDirection == 2) THEN
      patchNew%lbound    = 3
      patchNew%l1beg     = 1
      patchNew%l1end     = gridNew%kpc
      patchNew%l2beg     = 1
      patchNew%l2end     = gridNew%ipc
      patchNew%srcLbound = 4
      patchNew%srcL1beg  = 1
      patchNew%srcL1end  = gridNew%kpc
      patchNew%srcL2beg  = 1
      patchNew%srcL2end  = gridNew%ipc
    ELSE
      patchNew%lbound    = 5
      patchNew%l1beg     = 1
      patchNew%l1end     = gridNew%ipc
      patchNew%l2beg     = 1
      patchNew%l2end     = gridNew%jpc
      patchNew%srcLbound = 6
      patchNew%srcL1beg  = 1
      patchNew%srcL1end  = gridNew%ipc
      patchNew%srcL2beg  = 1
      patchNew%srcL2end  = gridNew%jpc
    ENDIF
  ENDDO  ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( globalNew )

END SUBROUTINE SplitTopology

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SPLT_SplitTopology.F90,v $
! Revision 1.4  2008/12/06 08:44:51  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:01  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/03 03:33:38  wasistho
! rflo_modinterfacessplit to splt_modinterfaces
!
! Revision 1.1  2004/12/03 02:41:15  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:45:30  wasistho
! lower to upper case
!
! Revision 1.6  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.5  2003/03/20 22:32:37  haselbac
! Renamed ModInterfaces
!
! Revision 1.4  2003/03/20 19:45:54  haselbac
! Corrected mistake in phased check-in
!
! Revision 1.3  2002/09/27 00:57:11  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.2  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/07/12 21:50:08  jblazek
! Added tool to split single grid into multiple regions.
!
!******************************************************************************







