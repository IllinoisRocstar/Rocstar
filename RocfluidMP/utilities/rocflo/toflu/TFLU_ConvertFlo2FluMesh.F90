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
! Purpose: collect sizes of required arrays and compute unstructured grid
!          and its connectivity from structured grid in the current region
!
! Description: collect array sizes at first call, i.e. when iFlag=0.
!              At second call (iFlag=1) compute the x,y,z coordinates of the
!              the unstructured vertices and obtain connectivity for all
!              all hexes in the whole region. The connectivity of quads and
!              vertices at boundary patches are computed in subroutine
!              ConvertFlo2FluPatch.
!
! Input: iFlag      = flag to determine the type of processing:
!                     0 for collection of array sizes
!                     1 for computation of unstructured grid and connectivity
!        iReg       = region number
!        regions    = data for all regions
!
! Output: global%tofluNPatches, NVerts, NHexs, NbfMax, NbnMax (iFlag=0)
!         grid%tofluLoc2g, global%tofluLoc2g, tofluXyz, tofluHex2v (iFlag=1)
!         global%tofluNFaces (iFlag=0), NEdges (iFlag=0 and 1)
!
! Notes: grid%tofluLoc2g is not deallocated in loop over regions in the main
!        main code ROCFLO_toFlu (unlike grid%xyz), as connectivity in previous 
!        regions are needed by the current region. This is done for simplicity
!        in coding, otherwise additional surface arrays for connectivity at
!        patches of previous regions are required including a copy procedure
!        to these arrays. The memory penalty due to this unallocated arrays
!        is believed to be insubstantial, as the values stored are integer.
!
!******************************************************************************
!
! $Id: TFLU_ConvertFlo2FluMesh.F90,v 1.5 2008/12/06 08:44:53 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ConvertFlo2FluMesh( iFlag,iReg,regions )

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE ModBndPatch, ONLY   : t_patch
  USE TFLU_ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, &
         RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset, RFLO_GetPatchIndices, &
         RFLO_GetPatchIndicesNodes, RFLO_GetPatchDirection, &
         RFLO_GetPatchMapping, RFLO_ReadGridRegion 
                               
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: iFlag, iReg
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iPatch, i, j, k, iPatchSrc

! ... local variables
  TYPE(t_global), POINTER :: global
  TYPE(t_grid) , POINTER  :: grid, gridSrc
  TYPE(t_patch) , POINTER :: patch, patchSrc

  INTEGER :: iLev, lbound, bcType, iNOff, ijNOff, ijkN, ibn, ien, found(6)
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: iaddb, jaddb, kaddb, iadde, jadde, kadde, madd
  INTEGER :: idir, jdir, kdir, idirSrc, jdirSrc, kdirSrc
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, dims(2), errorFlag
  INTEGER :: ibegSrc, iendSrc, jbegSrc, jendSrc, kbegSrc, kendSrc 
  INTEGER :: is, js, ks, l1SrcDir, l2SrcDir, mapMat(3,4)
  INTEGER :: h, iRegSrc, iLevSrc, lbs, l1bs, l1es, l2bs, l2es
  LOGICAL :: align

!******************************************************************************

  global => regions(1)%global
  CALL RegisterFunction( global,'ConvertFlo2FluMesh',&
  'TFLU_ConvertFlo2FluMesh.F90' )

! collect sizes and check bnd patches (iFlag=0) -------------------------------
! memory allocations and read grid (iFlag=1)
  
! allocate memory for grid

  iLev = regions(iReg)%currLevel

  CALL RFLO_GetDimensPhysNodes( regions(iReg),iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetDimensDummyNodes( regions(iReg),iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )

! detect any region face with mixed connecting and non-conn. bc.

  iaddb = 0
  jaddb = 0
  kaddb = 0
  iadde = 0
  jadde = 0
  kadde = 0

  found(:) = 0

  DO iPatch=1,regions(iReg)%nPatches

    patch   => regions(iReg)%levels(iLev)%patches(iPatch)
    lbound  =  patch%lbound
    bcType  =  patch%bcType
    iRegSrc =  patch%srcRegion

! - note, current algorithm only holds for connecting bc applied on whole face
!   of a region, not partially or mixed with other bcType

    CALL RFLO_GetPatchIndicesNodes( regions(iReg),patch,iLev, &
                                    ibeg,iend,jbeg,jend,kbeg,kend )

    IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .AND. &
        iRegSrc < iReg) THEN
      IF (lbound==1) iaddb =  1
      IF (lbound==2) iadde = -1
      IF (lbound==3) jaddb =  1
      IF (lbound==4) jadde = -1
      IF (lbound==5) kaddb =  1
      IF (lbound==6) kadde = -1

      IF (lbound==1) THEN
        IF (found(1)<0) GOTO 666
        found(1)=1
      ELSEIF (lbound==2) THEN
        IF (found(2)<0) GOTO 666
        found(2)=1
      ELSEIF (lbound==3) THEN
        IF (found(3)<0) GOTO 666
        found(3)=1
      ELSEIF (lbound==4) THEN
        IF (found(4)<0) GOTO 666
        found(4)=1
      ELSEIF (lbound==5) THEN
        IF (found(5)<0) GOTO 666
        found(5)=1
      ELSEIF (lbound==6) THEN
        IF (found(6)<0) GOTO 666
        found(6)=1
      ENDIF

    ELSE  ! bcType

      IF (lbound==1) THEN
        IF (found(1)>0) GOTO 666
        found(1)=-1
      ELSEIF (lbound==2) THEN
        IF (found(2)>0) GOTO 666
        found(2)=-1
      ELSEIF (lbound==3) THEN
        IF (found(3)>0) GOTO 666
        found(3)=-1
      ELSEIF (lbound==4) THEN
        IF (found(4)>0) GOTO 666
        found(4)=-1
      ELSEIF (lbound==5) THEN
        IF (found(5)>0) GOTO 666
        found(5)=-1
      ELSEIF (lbound==6) THEN
        IF (found(6)>0) GOTO 666
        found(6)=-1
      ENDIF
    ENDIF   ! bcType
  ENDDO     ! iPatch

  GOTO 777

666 CONTINUE

  WRITE(STDOUT,*)'ERROR: Region', iReg, &
                 ' has mixed conn./phys. bc on same block side'
!  CALL ErrorStop( global,ERR_PATCH_DIMENS,__LINE__, &
!          'Only works if connecting BC covers whole, not partial, patch' )

777 CONTINUE

  IF (iFlag == 0) THEN

! - increment global number of vertices and hexes for memory allocations

    DO k=kpnbeg+kaddb,kpnend+kadde
      DO j=jpnbeg+jaddb,jpnend+jadde
        DO i=ipnbeg+iaddb,ipnend+iadde
          global%tofluNVerts = global%tofluNVerts + 1  
        ENDDO
      ENDDO
    ENDDO

    DO k=kpnbeg,kpnend-1
      DO j=jpnbeg,jpnend-1
        DO i=ipnbeg,ipnend-1
          global%tofluNHexs = global%tofluNHexs + 1  
        ENDDO
      ENDDO
    ENDDO

! - similar for MAX number of vertices and faces (quads) at non-conn. patches
!   global%tofluNFaces rides along

    DO iPatch=1,regions(iReg)%nPatches
      patch   => regions(iReg)%levels(iLev)%patches(iPatch)
      bcType  =  patch%bcType

      IF (bcType<BC_REGIONCONF .OR. bcType>BC_REGIONCONF+BC_RANGE) THEN

        global%tofluNPatches = global%tofluNPatches + 1

        dims(1) = ABS(patch%l1end-patch%l1beg) + 1    ! faces values
        dims(2) = ABS(patch%l2end-patch%l2beg) + 1
        global%tofluNbfMax   = MAX( global%tofluNbfMax,dims(1)*dims(2) )
        global%tofluNFaces   = global%tofluNFaces + dims(1)*dims(2)

        dims(1) = ABS(patch%l1end-patch%l1beg) + 2    ! nodes values
        dims(2) = ABS(patch%l2end-patch%l2beg) + 2
        global%tofluNbnMax   = MAX( global%tofluNbnMax,dims(1)*dims(2) )
      ENDIF
    ENDDO

! - increment global number of faces and edges for dimension file

! - I-faces
    DO k=kpnbeg  ,kpnend-1
      DO j=jpnbeg  ,jpnend-1
        DO i=ipnbeg+1,ipnend-1
          global%tofluNFaces = global%tofluNFaces + 1  
        ENDDO
      ENDDO
    ENDDO

! - J-faces
    DO k=kpnbeg  ,kpnend-1
      DO j=jpnbeg+1,jpnend-1
        DO i=ipnbeg  ,ipnend-1
          global%tofluNFaces = global%tofluNFaces + 1  
        ENDDO
      ENDDO
    ENDDO

! - K-faces
    DO k=kpnbeg+1,kpnend-1
      DO j=jpnbeg  ,jpnend-1
        DO i=ipnbeg  ,ipnend-1
          global%tofluNFaces = global%tofluNFaces + 1  
        ENDDO
      ENDDO
    ENDDO

! - Edges
    DO k=kpnbeg,kpnend
      DO j=jpnbeg,jpnend
        DO i=ipnbeg,ipnend
          global%tofluNEdges = global%tofluNEdges + 3  
        ENDDO
      ENDDO
    ENDDO

! - Patch faces and edges
    DO iPatch=1,regions(iReg)%nPatches
      patch   => regions(iReg)%levels(iLev)%patches(iPatch)
      bcType  =  patch%bcType
      iRegSrc =  patch%srcRegion

      IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) &
           .AND. iRegSrc < iReg) THEN

        dims(1) = ABS(patch%l1end-patch%l1beg) + 1    ! faces values
        dims(2) = ABS(patch%l2end-patch%l2beg) + 1
        global%tofluNFaces = global%tofluNFaces + dims(1)*dims(2)

        dims(1) = ABS(patch%l1end-patch%l1beg) + 2    ! nodes values
        dims(2) = ABS(patch%l2end-patch%l2beg) + 2
        global%tofluNEdges = global%tofluNEdges - dims(1)*dims(2)*3
      ENDIF  ! bcType
    ENDDO

  ELSE   ! iFlag 

    ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
    ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

    grid => regions(iReg)%levels(iLev)%grid

    ALLOCATE( grid%xyz(3,ibn:ien),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    ALLOCATE( grid%tofluLoc2g(ipnbeg:ipnend,jpnbeg:jpnend,kpnbeg:kpnend), &
              stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

! - read volume grid

    CALL RFLO_ReadGridRegion( iReg,regions )

! - obtain current region local-to-global (l2g) mapping -------------------------

    DO k=kpnbeg+kaddb,kpnend+kadde
      DO j=jpnbeg+jaddb,jpnend+jadde
        DO i=ipnbeg+iaddb,ipnend+iadde
          ijkN = IndIJK(i,j,k,iNOff,ijNOff)

          global%tofluNVerts = global%tofluNVerts + 1  
          global%tofluXyz(1,global%tofluNVerts) = grid%xyz(XCOORD,ijkN)  
          global%tofluXyz(2,global%tofluNVerts) = grid%xyz(YCOORD,ijkN)  
          global%tofluXyz(3,global%tofluNVerts) = grid%xyz(ZCOORD,ijkN)
          grid%tofluLoc2g(i,j,k) = global%tofluNVerts
        ENDDO
      ENDDO
    ENDDO
  ENDIF  ! iFlag

  IF (iFlag == 0) GOTO 888

! l2g mapping at connecting patches --------------------------------------------

  DO iPatch=1,regions(iReg)%nPatches

    patch     => regions(iReg)%levels(iLev)%patches(iPatch)
    lbound    =  patch%lbound
    bcType    =  patch%bcType
    iRegSrc   =  patch%srcRegion

! - only at connecting patches and previous source regions

    IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF + BC_RANGE) .AND. &
        iRegSrc < iReg) THEN

      align     = patch%align 
      iLevSrc   = regions(iRegSrc)%currLevel
      lbs       = patch%srcLbound
      l1bs      = MIN(patch%srcL1beg,patch%srcL1end)
      l1es      = MAX(patch%srcL1beg,patch%srcL1end)
      l2bs      = MIN(patch%srcL2beg,patch%srcL2end)
      l2es      = MAX(patch%srcL2beg,patch%srcL2end)

! --- search for iPatchSrc

      DO iPatchSrc=1,regions(iRegSrc)%nPatches
        patchSrc => regions(iRegSrc)%levels(iLevSrc)%patches(iPatchSrc)

        IF (patchSrc%bcType==bcType .AND. &
            patchSrc%lbound==lbs    .AND. &
            patchSrc%l1beg ==l1bs   .AND. &
            patchSrc%l1end ==l1es   .AND. &
            patchSrc%l2beg ==l2bs   .AND. &
            patchSrc%l2end ==l2es) THEN       ! OK, iPatchSrc found
          patch%srcPatch = iPatchSrc
        ENDIF
      ENDDO    ! iPatchSrc

      iPatchSrc =  patch%srcPatch
      patchSrc  => regions(iRegSrc)%levels(iLevSrc)%patches(iPatchSrc)
      gridSrc   => regions(iRegSrc)%levels(iLevSrc)%grid

! --- copy l2g mapping from the source patch already processed to current patch

      CALL RFLO_GetPatchIndicesNodes( regions(iReg),patch,iLev, &
                             ibeg,iend,jbeg,jend,kbeg,kend )
      CALL RFLO_GetPatchIndicesNodes( regions(iRegSrc),patchSrc,iLev, &
                             ibegSrc,iendSrc,jbegSrc,jendSrc,kbegSrc,kendSrc )
      CALL RFLO_GetPatchDirection( patch   ,idir   ,jdir,   kdir    )
      CALL RFLO_GetPatchDirection( patchSrc,idirSrc,jdirSrc,kdirSrc )

                                           l1SrcDir =  1
      IF (patch%srcL1beg > patch%srcL1end) l1SrcDir = -1
                                           l2SrcDir =  1
      IF (patch%srcL2beg > patch%srcL2end) l2SrcDir = -1

      CALL RFLO_GetPatchMapping( lbound,lbs,l1SrcDir,l2SrcDir,align, &
                                 idir,jdir,kdir,idirSrc,jdirSrc,kdirSrc, &
                                 ibeg,iend,jbeg,jend,kbeg,kend, &
                                 ibegSrc,iendSrc,jbegSrc,jendSrc, &
                                 kbegSrc,kendSrc,mapMat )

      IF (lbs==1) mapMat(1,4) = mapMat(1,4)+1
      IF (lbs==2) mapMat(1,4) = mapMat(1,4)-1
      IF (lbs==3) mapMat(2,4) = mapMat(2,4)+1
      IF (lbs==4) mapMat(2,4) = mapMat(2,4)-1
      IF (lbs==5) mapMat(3,4) = mapMat(3,4)+1
      IF (lbs==6) mapMat(3,4) = mapMat(3,4)-1
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            is = i*mapMat(1,1)+j*mapMat(1,2)+k*mapMat(1,3)+mapMat(1,4)
            js = i*mapMat(2,1)+j*mapMat(2,2)+k*mapMat(2,3)+mapMat(2,4)
            ks = i*mapMat(3,1)+j*mapMat(3,2)+k*mapMat(3,3)+mapMat(3,4)
            grid%tofluLoc2g(i,j,k) = gridSrc%tofluLoc2g(is,js,ks)
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k

    ENDIF ! bcType
  ENDDO ! iPatch

! obtain global hex connectivity ----------------------------------------------

  DO k=kpnbeg,kpnend-1
    DO j=jpnbeg,jpnend-1
      DO i=ipnbeg,ipnend-1
        global%tofluNHexs = global%tofluNHexs + 1  
        h = global%tofluNHexs
        global%tofluHex2v(1,h) = grid%tofluLoc2g(i  ,j  ,k  )
        global%tofluHex2v(2,h) = grid%tofluLoc2g(i  ,j  ,k+1)
        global%tofluHex2v(3,h) = grid%tofluLoc2g(i+1,j  ,k+1)
        global%tofluHex2v(4,h) = grid%tofluLoc2g(i+1,j  ,k  )
        global%tofluHex2v(5,h) = grid%tofluLoc2g(i  ,j+1,k  )
        global%tofluHex2v(6,h) = grid%tofluLoc2g(i  ,j+1,k+1)
        global%tofluHex2v(7,h) = grid%tofluLoc2g(i+1,j+1,k+1)
        global%tofluHex2v(8,h) = grid%tofluLoc2g(i+1,j+1,k  )
      ENDDO  ! i
    ENDDO  ! j
  ENDDO  ! k

888 CONTINUE

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE ConvertFlo2FluMesh

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TFLU_ConvertFlo2FluMesh.F90,v $
! Revision 1.5  2008/12/06 08:44:53  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:18:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/03/07 03:23:48  wasistho
! enabled serial execution
!
! Revision 1.2  2004/12/03 03:43:44  wasistho
! rflo_modinterfacestoflu to tflu_modinterfaces
!
! Revision 1.1  2004/12/03 02:59:30  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:58:20  wasistho
! lower to upper case
!
! Revision 1.2  2004/08/18 02:10:33  wasistho
! added new routines to create dimension file
!
! Revision 1.1.1.1  2004/08/17 01:41:39  wasistho
! initial checkin
!
!
!******************************************************************************







