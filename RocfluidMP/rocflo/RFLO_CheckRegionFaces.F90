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
! Purpose: check if boundary conditions are specified for all boundary
!          cell faces.
!
! Description: (1) mark all cells with 0
!              (2) mark boundary cells of all patches (add values to find
!                  overlaps and to check later edges and corners of regions)
!              (3) check if if all boundary cells are set to 1 (or =2 at edges,
!                  =3 at corners).
!
! Input: regions = dimensions and topology of all regions.
!
! Output: none.
!
! Notes: patch indices are checked to be within region dimensions.
!
!******************************************************************************
!
! $Id: RFLO_CheckRegionFaces.F90,v 1.3 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CheckRegionFaces( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
                            RFLO_GetPatchIndices
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iLev, iPatch, i, j, k

! ... local variables
  CHARACTER(CHRLEN) :: msg

  INTEGER          :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend, errorFlag
  INTEGER          :: lbound, ibeg, iend, jbeg, jend, kbeg, kend, errcode
  INTEGER, POINTER :: marker(:,:,:)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_CheckRegionFaces',&
  'RFLO_CheckRegionFaces.F90' )

  IF (global%myProcid > global%nRegions-1) RETURN

! start -----------------------------------------------------------------------

  NULLIFY( marker )

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

      DO iLev=1,regions(iReg)%nGridLevels
        CALL RFLO_GetDimensPhys( regions(iReg),iLev,ipcbeg,ipcend, &
                                 jpcbeg,jpcend,kpcbeg,kpcend )

! ----- allocate marker and set = 0

        IF (ASSOCIATED(marker)) THEN
          DEALLOCATE( marker,stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
        ENDIF
        ALLOCATE( marker(ipcbeg:ipcend,jpcbeg:jpcend,kpcbeg:kpcend),stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

        marker(:,:,:) = 0

! ----- loop over patches and mark boundary cells

        DO iPatch=1,regions(iReg)%nPatches

          patch  => regions(iReg)%levels(iLev)%patches(iPatch)
          lbound =  patch%lbound
          CALL RFLO_GetPatchIndices( regions(iReg), &
                    regions(iReg)%levels(iLev)%patches(iPatch), &
                    iLev,ibeg,iend,jbeg,jend,kbeg,kend )

          IF ((ibeg<ipcbeg .OR. ibeg>ipcend) .OR. &   ! check dimensions
              (iend<ipcbeg .OR. iend>ipcend) .OR. &
              (jbeg<jpcbeg .OR. jbeg>jpcend) .OR. &
              (jend<jpcbeg .OR. jend>jpcend) .OR. &
              (kbeg<kpcbeg .OR. kbeg>kpcend) .OR. &
              (kend<kpcbeg .OR. kend>kpcend)) THEN
            WRITE(msg,1000) iReg,iLev,iPatch,patch%bcType
            CALL ErrorStop( global,ERR_PATCH_RANGE,__LINE__,msg )
          ENDIF

          DO k=kbeg,kend
            DO j=jbeg,jend
              DO i=ibeg,iend
                marker(i,j,k) = marker(i,j,k) + 1
              ENDDO  ! i
            ENDDO    ! j
          ENDDO      ! k

        ENDDO   ! iPatch

! ----- check faces 1 and 2

        DO k=kpcbeg,kpcend
          DO j=jpcbeg,jpcend

            IF (marker(ipcbeg,j,k) == 0) THEN
              WRITE(msg,1005) iReg,iLev,1,ipcbeg,j,k
              CALL ErrorStop( global,ERR_PATCH_NOTCOVERED,__LINE__,msg )
            ENDIF
            IF (k>kpcbeg .AND. k<kpcend .AND. &
                j>jpcbeg .AND. j<jpcend .AND. marker(ipcbeg,j,k) > 1) THEN
              WRITE(msg,1005) iReg,iLev,1,ipcbeg,j,k
              CALL ErrorStop( global,ERR_PATCH_OVERLAP,__LINE__,msg )
            ENDIF

            IF (marker(ipcend,j,k) == 0) THEN
              WRITE(msg,1005) iReg,iLev,2,ipcend,j,k
              CALL ErrorStop( global,ERR_PATCH_NOTCOVERED,__LINE__,msg )
            ENDIF
            IF (k>kpcbeg .AND. k<kpcend .AND. &
                j>jpcbeg .AND. j<jpcend .AND. marker(ipcend,j,k) > 1) THEN
              WRITE(msg,1005) iReg,iLev,1,ipcend,j,k
              CALL ErrorStop( global,ERR_PATCH_OVERLAP,__LINE__,msg )
            ENDIF

            IF (k==kpcbeg .OR. k==kpcend) THEN
              IF (j==jpcbeg .OR. j==jpcend) THEN
                IF (marker(ipcbeg,j,k) /= 3) THEN
                  WRITE(msg,1005) iReg,iLev,1,ipcbeg,j,k
                  IF (marker(ipcbeg,j,k) < 3) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(ipcbeg,j,k) > 3) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
                IF (marker(ipcend,j,k) /= 3) THEN
                  WRITE(msg,1005) iReg,iLev,2,ipcend,j,k
                  IF (marker(ipcend,j,k) < 3) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(ipcend,j,k) > 3) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
              ELSE
                IF (marker(ipcbeg,j,k) /= 2) THEN
                  WRITE(msg,1005) iReg,iLev,1,ipcbeg,j,k
                  IF (marker(ipcbeg,j,k) < 2) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(ipcbeg,j,k) > 2) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
                IF (marker(ipcend,j,k) /= 2) THEN
                  WRITE(msg,1005) iReg,iLev,2,ipcend,j,k
                  IF (marker(ipcend,j,k) < 2) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(ipcend,j,k) > 2) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
              ENDIF
            ENDIF

            IF (j==jpcbeg .OR. j==jpcend) THEN
              IF (k/=kpcbeg .AND. k/=kpcend) THEN
                IF (marker(ipcbeg,j,k) /= 2) THEN
                  WRITE(msg,1005) iReg,iLev,1,ipcbeg,j,k
                  IF (marker(ipcbeg,j,k) < 2) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(ipcbeg,j,k) > 2) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
                IF (marker(ipcend,j,k) /= 2) THEN
                  WRITE(msg,1005) iReg,iLev,2,ipcend,j,k
                  IF (marker(ipcend,j,k) < 2) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(ipcend,j,k) > 2) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
              ENDIF
            ENDIF

          ENDDO   ! j
        ENDDO     ! k

! ----- check faces 3 and 4

        DO k=kpcbeg,kpcend
          DO i=ipcbeg,ipcend

            IF (marker(i,jpcbeg,k) == 0) THEN
              WRITE(msg,1005) iReg,iLev,3,i,jpcbeg,k
              CALL ErrorStop( global,ERR_PATCH_NOTCOVERED,__LINE__,msg )
            ENDIF
            IF (k>kpcbeg .AND. k<kpcend .AND. &
                i>ipcbeg .AND. i<ipcend .AND. marker(i,jpcbeg,k) > 1) THEN
              WRITE(msg,1005) iReg,iLev,3,i,jpcbeg,k
              CALL ErrorStop( global,ERR_PATCH_OVERLAP,__LINE__,msg )
            ENDIF

            IF (marker(i,jpcend,k) == 0) THEN
              WRITE(msg,1005) iReg,iLev,4,i,jpcend,k
              CALL ErrorStop( global,ERR_PATCH_NOTCOVERED,__LINE__,msg )
            ENDIF
            IF (k>kpcbeg .AND. k<kpcend .AND. &
                i>ipcbeg .AND. i<ipcend .AND. marker(i,jpcend,k) > 1) THEN
              WRITE(msg,1005) iReg,iLev,4,i,jpcend,k
              CALL ErrorStop( global,ERR_PATCH_OVERLAP,__LINE__,msg )
            ENDIF

            IF (k==kpcbeg .OR. k==kpcend) THEN
              IF (i==ipcbeg .OR. i==ipcend) THEN
                IF (marker(i,jpcbeg,k) /= 3) THEN
                  WRITE(msg,1005) iReg,iLev,3,i,jpcbeg,k
                  IF (marker(i,jpcbeg,k) < 3) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(i,jpcbeg,k) > 3) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
                IF (marker(i,jpcend,k) /= 3) THEN
                  WRITE(msg,1005) iReg,iLev,4,i,jpcend,k
                  IF (marker(i,jpcend,k) < 3) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(i,jpcend,k) > 3) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
              ELSE
                IF (marker(i,jpcbeg,k) /= 2) THEN
                  WRITE(msg,1005) iReg,iLev,3,i,jpcbeg,k
                  IF (marker(i,jpcbeg,k) < 2) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(i,jpcbeg,k) > 2) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
                IF (marker(i,jpcend,k) /= 2) THEN
                  WRITE(msg,1005) iReg,iLev,4,i,jpcend,k
                  IF (marker(i,jpcend,k) < 2) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(i,jpcend,k) > 2) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
              ENDIF
            ENDIF

            IF (i==ipcbeg .OR. i==ipcend) THEN
              IF (k/=kpcbeg .AND. k/=kpcend) THEN
                IF (marker(i,jpcbeg,k) /= 2) THEN
                  WRITE(msg,1005) iReg,iLev,3,i,jpcbeg,k
                  IF (marker(i,jpcbeg,k) < 2) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(i,jpcbeg,k) > 2) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
                IF (marker(i,jpcend,k) /= 2) THEN
                  WRITE(msg,1005) iReg,iLev,4,i,jpcend,k
                  IF (marker(i,jpcend,k) < 2) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(i,jpcend,k) > 2) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
              ENDIF
            ENDIF

          ENDDO   ! i
        ENDDO     ! k

! ----- check faces 5 and 6

        DO j=jpcbeg,jpcend
          DO i=ipcbeg,ipcend

            IF (marker(i,j,kpcbeg) == 0) THEN
              WRITE(msg,1005) iReg,iLev,5,i,j,kpcbeg
              CALL ErrorStop( global,ERR_PATCH_NOTCOVERED,__LINE__,msg )
            ENDIF
            IF (i>ipcbeg .AND. i<ipcend .AND. &
                j>jpcbeg .AND. j<jpcend .AND. marker(i,j,kpcbeg) > 1) THEN
              WRITE(msg,1005) iReg,iLev,5,i,j,kpcbeg
              CALL ErrorStop( global,ERR_PATCH_OVERLAP,__LINE__,msg )
            ENDIF
            IF (marker(i,j,kpcend) == 0) THEN
              WRITE(msg,1005) iReg,iLev,6,i,j,kpcend
              CALL ErrorStop( global,ERR_PATCH_NOTCOVERED,__LINE__,msg )
            ENDIF
            IF (i>ipcbeg .AND. i<ipcend .AND. &
                j>jpcbeg .AND. j<jpcend .AND. marker(i,j,kpcend) > 1) THEN
              WRITE(msg,1005) iReg,iLev,6,i,j,kpcend
              CALL ErrorStop( global,ERR_PATCH_OVERLAP,__LINE__,msg )
            ENDIF

            IF (j==jpcbeg .OR. j==jpcend) THEN
              IF (i==ipcbeg .OR. i==ipcend) THEN
                IF (marker(i,j,kpcbeg) /= 3) THEN
                  WRITE(msg,1005) iReg,iLev,5,i,j,kpcbeg
                  IF (marker(i,j,kpcbeg) < 3) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(i,j,kpcbeg) > 3) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
                IF (marker(i,j,kpcend) /= 3) THEN
                  WRITE(msg,1005) iReg,iLev,6,i,j,kpcend
                  IF (marker(i,j,kpcend) < 3) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(i,j,kpcend) > 3) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
              ELSE
                IF (marker(i,j,kpcbeg) /= 2) THEN
                  WRITE(msg,1005) iReg,iLev,5,i,j,kpcbeg
                  IF (marker(i,j,kpcbeg) < 2) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(i,j,kpcbeg) > 2) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
                IF (marker(i,j,kpcend) /= 2) THEN
                  WRITE(msg,1005) iReg,iLev,6,i,j,kpcend
                  IF (marker(i,j,kpcend) < 2) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(i,j,kpcend) > 2) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
              ENDIF
            ENDIF

            IF (i==ipcbeg .OR. i==ipcend) THEN
              IF (j/=jpcbeg .AND. j/=jpcend) THEN
                IF (marker(i,j,kpcbeg) /= 2) THEN
                  WRITE(msg,1005) iReg,iLev,5,i,j,kpcbeg
                  IF (marker(i,j,kpcbeg) < 2) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(i,j,kpcbeg) > 2) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
                IF (marker(i,j,kpcend) /= 2) THEN
                  WRITE(msg,1005) iReg,iLev,6,i,j,kpcend
                  IF (marker(i,j,kpcend) < 2) errcode = ERR_PATCH_NOTCOVERED
                  IF (marker(i,j,kpcend) > 2) errcode = ERR_PATCH_OVERLAP
                  CALL ErrorStop( global,errcode,__LINE__,msg )
                ENDIF
              ENDIF
            ENDIF

          ENDDO   ! i
        ENDDO     ! j

      ENDDO       ! iLev

    ENDIF         ! region on this processor and active
  ENDDO           ! iReg

! finalize --------------------------------------------------------------------

  IF (ASSOCIATED(marker)) DEALLOCATE( marker,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  CALL DeregisterFunction( global )

1000 FORMAT('Region ',I5,', level ',I1,', patch ',I3,', BC type ',I3)
1005 FORMAT('Region ',I5,', level ',I1,', face ',I1,', i,j,k= ',3I6)

END SUBROUTINE RFLO_CheckRegionFaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CheckRegionFaces.F90,v $
! Revision 1.3  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.13  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.9  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.8  2003/02/14 22:32:36  jblazek
! Finished implementation of corener and edge cells.
!
! Revision 1.7  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.6  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.5  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/03/30 00:50:49  jblazek
! Cleaned up with flint.
!
! Revision 1.3  2002/03/22 19:26:24  jblazek
! Unused pointers are nullified.
!
! Revision 1.2  2002/03/18 23:11:32  jblazek
! Finished multiblock and MPI.
!
! Revision 1.1  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.5  2002/02/04 15:30:25  jblazek
! Added injection boundary condition.
!
! Revision 1.4  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.3  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.2  2001/12/22 00:09:39  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1  2001/12/19 23:09:21  jblazek
! Added routines to read grid and solution.
!
!******************************************************************************







