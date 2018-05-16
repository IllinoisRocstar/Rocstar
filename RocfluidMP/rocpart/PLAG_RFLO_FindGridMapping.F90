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
! Purpose: find the grid transformation between regions connecting at 
!          corner-edge cells for proper metrics.
!
! Description: none.
!
! Input: regions = data of all regions,
!        iReg    = current region number.
!
! Output: indexMapMat   = index mapping matrix.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_RFLO_FindGridMapping.F90,v 1.4 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLO_FindGridMapping( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region, t_level, t_dCell
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag 
  USE ModIndexing, ONLY   : GetIJK
  USE ModInterfaces, ONLY : RFLO_GetCellOffset,         &
                            RFLO_GetCornerCellsIndices, &
                            RFLO_GetEdgeCellsIndices
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  USE PLAG_ModInterfaces, ONLY : PLAG_RFLO_FindSourceCell, &
                                 PLAG_RFLO_GetFaceMapping
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)


! ... loop variables
  INTEGER :: i, iCorner, iEdge, iReg, ijk, j, k, nt

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iLev, nCorners, nEdges, ntEdge, ntCorn
  INTEGER :: iRegSrc, ijkN, ic, jc, kc
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: iNOff, ijNOff

  INTEGER :: nDumCellsSrc
  INTEGER :: iRegTemp, iRegSrcTemp, iGridMapMat
  INTEGER :: srcDir(3),srcFace(6)
  INTEGER :: srcIndexMapMat(3,4)
  
  LOGICAL :: found

  TYPE(t_region),      POINTER :: pRegion
  TYPE(t_level),       POINTER :: pLevel
  TYPE(t_global),      POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PLAG_RFLO_FindGridMapping.F90,v $ $Revision: 1.4 $'

  global => regions(1)%global
    
  CALL RegisterFunction( global, 'PLAG_RFLO_FindGridMapping',&
  'PLAG_RFLO_FindGridMapping.F90' )

! Get dimensions --------------------------------------------------------------

  nCorners = 8
  nEdges   = 12
  ntEdge   = 2
  ntCorn   = 3

! Search for source regions ===================================================

  DO iReg=1,global%nRegions
    iLev=1

! - Set pointers --------------------------------------------------------------
  
      pLevel  => regions(iReg)%levels(iLev)

      pRegion => regions(iReg)

! - Get node offset -----------------------------------------------------------

      CALL RFLO_GetNodeOffset( pRegion,iLev,iNOff,ijNOff )

! - Loop over edges -----------------------------------------------------------

      DO iEdge=1,nEdges 

! -- Bypass for noninteracting regions ----------------------------------------

        IF( .NOT. pLevel%edgeCells(iEdge)%interact ) GOTO 1999

! -- Loop over edge cell indices ----------------------------------------------
        
        CALL RFLO_GetEdgeCellsIndices( pRegion,iLev,iEdge, &
                                       ibeg,iend,jbeg,jend,kbeg,kend )
        ijk = 0 
        DO k=kbeg,kend
        DO j=jbeg,jend
        DO i=ibeg,iend
          ijk = ijk + 1

! --- Source region infrastructure --------------------------------------------
            
          iRegSrc = pLevel%edgeCells(iEdge)%cells(ijk)%srcRegion
     
          ic = i; jc = j; kc =k;
          iRegTemp=iReg 
          DO nt=1,ntEdge
            CALL PLAG_RFLO_FindSourceCell( regions,iRegTemp,iLev,ic,jc,kc, &
                                           found,iRegSrcTemp,srcIndexMapMat )
            IF(found) EXIT
            iRegTemp = iRegSrcTemp
          END DO ! nt

! --- Extract mapping between regions sharing an edge -------------------------
                           
          pLevel%edgeCells(iEdge)%cells(ijk)%srcIndexMapMat = srcIndexMapMat
        ENDDO   ! i
        ENDDO   ! j
        ENDDO   ! k

1999    CONTINUE 
      
      ENDDO      ! iEdge
      
! - Loop over corners ---------------------------------------------------------

      DO iCorner=1,nCorners

! -- Bypass for noninteracting regions ----------------------------------------

        IF( .NOT. pLevel%cornerCells(iCorner)%interact ) GOTO 2999

! -- Loop over corner cell indices --------------------------------------------
        
        CALL RFLO_GetCornerCellsIndices( pRegion,iLev,iCorner, &
                                         ibeg,iend,jbeg,jend,kbeg,kend )
        
        ijk = 0 
        DO k=kbeg,kend
        DO j=jbeg,jend
        DO i=ibeg,iend
          ijk = ijk + 1
  
! -- Source region infrastructure ---------------------------------------------
            
          iRegSrc  =  pLevel%cornerCells(iCorner)%cells(ijk)%srcRegion 
          
          ic = i; jc = j; kc =k;
          iRegTemp=iReg
          DO nt=1,ntCorn
            CALL PLAG_RFLO_FindSourceCell( regions,iRegTemp,iLev,ic,jc,kc, &
                                           found,iRegSrcTemp,srcIndexMapMat )
            IF(found) EXIT
            iRegTemp = iRegSrcTemp
          END DO ! nt

! -- Extract mapping between regions sharing a corner -------------------------

          pLevel%cornerCells(iCorner)%cells(ijk)%srcIndexMapMat = srcIndexMapMat
        ENDDO   ! i
        ENDDO   ! j
        ENDDO   ! k

2999    CONTINUE

      ENDDO     ! iCorner
  ENDDO         ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )
  
END SUBROUTINE PLAG_RFLO_FindGridMapping

!******************************************************************************
SUBROUTINE PLAG_RFLO_FindSourceCell( regions,iReg,iLev,ic,jc,kc, &
                                     found,iRegSrc,indexMapMat )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetPatchIndices
  USE ModError
  USE ModParameters
  USE PLAG_ModInterfaces, ONLY : PLAG_RFLO_SourceCell
  IMPLICIT NONE

! ... parameters
  INTEGER, INTENT(IN)    :: iReg, iLev
  INTEGER, INTENT(INOUT) :: ic, jc, kc
  INTEGER, INTENT(OUT)   :: iRegSrc
  INTEGER, INTENT(OUT)   :: indexMapMat(3,4)
 
  LOGICAL, INTENT(OUT)    :: found

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

  CALL RegisterFunction( global,'PLAG_RFLO_FindSourceCell',&
  'PLAG_RFLO_FindGridMapping.F90' )

  found = .false.

! Get dimensions --------------------------------------------------------------

  CALL RFLO_GetDimensPhys( regions(iReg),iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )

! Find a suitable patch to start from -----------------------------------------

  hit = .false.     ! no patch found yet

! Check if cell is within the patch -------------------------------------------

  DO iPatch=1,regions(iReg)%nPatches
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

! cell just outside the patch? ------------------------------------------------

  IF (.NOT. hit) THEN
    DO iPatch=1,regions(iReg)%nPatches
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

! cell at some corner? --------------------------------------------------------

  IF (.NOT. hit) THEN
    DO iPatch=1,regions(iReg)%nPatches
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
    CALL PLAG_RFLO_SourceCell( regions(iReg),regions(iRegSrc),patch,patchSrc, &
                               iLev,ic,jc,kc,found,indexMapMat )
  ELSE
    iRegSrc = iReg
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLO_FindSourceCell

! #############################################################################
! #############################################################################

SUBROUTINE PLAG_RFLO_SourceCell( region,regionSrc,patch,patchSrc, &
                                  iLev,ic,jc,kc,found,indexMapMat )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModIndexing, ONLY   : IndIJKMap, GetIJK
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
        RFLO_GetCellOffset, RFLO_GetPatchMapping, RFLO_GetDimensPhys
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  INTEGER, INTENT(IN)    :: iLev
  INTEGER, INTENT(INOUT) :: ic,jc,kc
  INTEGER, INTENT(OUT)   :: indexMapMat(3,4)

  LOGICAL, INTENT(INOUT) :: found

  TYPE(t_region)          :: region, regionSrc
  TYPE(t_patch), POINTER  :: patch, patchSrc

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: ibegSrc, iendSrc, jbegSrc, jendSrc, kbegSrc, kendSrc, &
             idirSrc, jdirSrc, kdirSrc, iCOffSrc, ijCOffSrc, ijkCSrc
  INTEGER :: lb, lbs, l1SrcDir, l2SrcDir, mapMat(3,4)
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: icell
  INTEGER :: icSrc, jcSrc, kcSrc

  LOGICAL :: align

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_SourceCell',&
  'PLAG_RFLO_FindGridMapping.F90' )

! Get dimensions --------------------------------------------------------------

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend, &
                             jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchIndices( regionSrc,patchSrc,iLev,ibegSrc,iendSrc, &
                             jbegSrc,jendSrc,kbegSrc,kendSrc )
  CALL RFLO_GetPatchDirection( patch   ,idir   ,jdir,   kdir    )
  CALL RFLO_GetPatchDirection( patchSrc,idirSrc,jdirSrc,kdirSrc )
  CALL RFLO_GetCellOffset( regionSrc,iLev,iCOffSrc,ijCOffSrc )

! Determine mapping between patches -------------------------------------------

                                       l1SrcDir =  1
  IF (patch%srcL1beg > patch%srcL1end) l1SrcDir = -1
                                       l2SrcDir =  1
  IF (patch%srcL2beg > patch%srcL2end) l2SrcDir = -1

  lb    = patch%lbound
  lbs   = patch%srcLbound
  align = patch%align

  CALL RFLO_GetPatchMapping( lb,lbs,l1SrcDir,l2SrcDir,align, &
                             idir,jdir,kdir,idirSrc,jdirSrc,kdirSrc, &
                             ibeg,iend,jbeg,jend,kbeg,kend, &
                             ibegSrc,iendSrc,jbegSrc,jendSrc,kbegSrc,kendSrc, &
                             mapMat )

! Get cell indices on source patch --------------------------------------------

  icell = IndIJKMap( ic,jc,kc,mapMat,iCOffSrc,ijCOffSrc )
  CALL GetIJK( icell,iCOffSrc,ijCOffSrc,regionSrc%nDumCells,ic,jc,kc )

! Check if ic,jc,kc within physical domain ------------------------------------

  CALL RFLO_GetDimensPhys( regionSrc,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )

  IF ((ic>=ipcbeg .AND. ic<=ipcend) .AND. &
      (jc>=jpcbeg .AND. jc<=jpcend) .AND. &
      (kc>=kpcbeg .AND. kc<=kpcend)) THEN
    found = .true.
    indexMapMat = mapMat
  ELSE
    found = .false.
    indexMapMat = -999999
  ENDIF

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLO_SourceCell

! #############################################################################
! #############################################################################

SUBROUTINE PLAG_RFLO_GetFaceMapping( mapMat,srcDir,srcFace )

  USE ModDataTypes
  USE ModError
  IMPLICIT NONE

! ... parameters
  INTEGER, INTENT(IN)  :: mapMat(3,4)
  INTEGER, INTENT(OUT) :: srcDir(3)
  INTEGER, INTENT(OUT) :: srcFace(6)

! ... local variables
  INTEGER :: idir, jdir, kdir
  INTEGER :: idirSrc, jdirSrc,kdirSrc 
  INTEGER :: iface,jface,kface
  INTEGER :: ifaceSrc,jfaceSrc,kfaceSrc
  INTEGER :: ifaceSrcFlip,jfaceSrcFlip,kfaceSrcFlip
  INTEGER :: ifacePSrcFlip,jfacePSrcFlip,kfacePSrcFlip

!******************************************************************************

! Set dimensions --------------------------------------------------------------
  
  idir=1;jdir=2;kdir=3;
  iface =1; jface=3; kface=5;
  
! Determine grid mapping of source region -------------------------------------
!  use absolute value since negative value is meaningless ---------------------

  idirSrc = ABS( idir*mapMat(1,1) + jdir*mapMat(1,2) + kdir*mapMat(1,3) )
  jdirSrc = ABS( idir*mapMat(2,1) + jdir*mapMat(2,2) + kdir*mapMat(2,3) )
  kdirSrc = ABS( idir*mapMat(3,1) + jdir*mapMat(3,2) + kdir*mapMat(3,3) ) 

! determine face mapping of source region -------------------------------------
           
  ifaceSrc = iface*mapMat(1,1) + jface*mapMat(1,2) + kface*mapMat(1,3) 
  jfaceSrc = iface*mapMat(2,1) + jface*mapMat(2,2) + kface*mapMat(2,3) 
  kfaceSrc = iface*mapMat(3,1) + jface*mapMat(3,2) + kface*mapMat(3,3)     

! Swap face indices for negative values ---------------------------------------

  ifaceSrcFlip  = ifaceSrc
  ifacePSrcFlip = ifaceSrc+1

  jfaceSrcFlip  = jfaceSrc
  jfacePSrcFlip = jfaceSrc+1
           
  kfaceSrcFlip =  kfaceSrc
  kfacePSrcFlip = kfaceSrc+1       
           
  IF ( ifaceSrc < 0 ) THEN
    ifaceSrcFlip  = ABS(ifaceSrc)+1 
    ifacePSrcFlip = ABS(ifaceSrc)
  ENDIF
           
  IF ( jfaceSrc < 0 ) THEN
    jfaceSrcFlip  = ABS(jfaceSrc)+1 
    jfacePSrcFlip = ABS(jfaceSrc)
  ENDIF
           
  IF ( kfaceSrc < 0 ) THEN
    kfaceSrcFlip  = ABS(kfaceSrc)+1 
    kfacePSrcFlip = ABS(kfaceSrc)
  ENDIF

! Load face indices of source region ------------------------------------------

  srcDir(1:3) = (/idirSrc,jdirSrc,kdirSrc/)
           
  srcFace(1:6)=(/ifaceSrcFlip,ifacePSrcFlip,jfaceSrcFlip,jfacePSrcFlip, &
                 kfaceSrcFlip,kfacePSrcFlip/)

END SUBROUTINE PLAG_RFLO_GetFaceMapping

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLO_FindGridMapping.F90,v $
! Revision 1.4  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 20:58:10  fnajjar
! Initial revision after changing case
!
! Revision 1.1  2004/02/10 21:24:06  fnajjar
! Initial import of index mapping capability for corner-edge regions
!
!******************************************************************************









