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
! ******************************************************************************
!
! Purpose: Suite pertaining degenerate edges/corners routines.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModDegenerateCornEdge.F90,v 1.3 2008/12/06 08:44:15 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModDegenerateCornEdge

  USE ModGlobal, ONLY    : t_global 
  USE ModDataStruct, ONLY: t_region, t_level
  USE ModBndPatch, ONLY  : t_patch
  USE ModGrid, ONLY      : t_grid
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLO_FindDegeneratCell, &
            RFLO_MarkDegeneratVert, &
            RFLO_WriteDegeneratEC
 
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModDegenerateCornEdge.F90,v $ $Revision: 1.3 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS

!******************************************************************************
!
! Purpose: given the indices of a dummy cell, find out if this cell is
!          degenerative and what level of degeration.
!
! Description: at first call (handleCorn=0), find if the test edge lays within
!              a patch of the same region, at second call (handleCorn=1) find
!              if the test corner lays within an edge or a patch of the same 
!              region.
!
! Input: regions    = dimensions and topology of all regions
!        handleCorn = treat corner (1) or edge (0)
!        icount     = i^th component of corner (handle=1) or edge (handle=0)
!        iReg       = current region
!        iLev       = current grid level
!        iRegTest   = region number of the test dimmy cell
!        icellTest  = index of the test dimmy cell
!        i/j/kTest  = indices of the test dummy cell
!
! Output: global%degenrtEc = .true. if the test cell is degenerative
!         level%edgeCell(icount)%degenrt   = level of edge degeneration
!         level%cornerCell(icount)%degenrt = level of level degeneration
!
! Notes: type of degenerative edge   = 0 (none)
!                                    = 1 (edge within patch)
!                                         when < 4 int. cells meet at the edge 
!                                    < 0 (detached from adjacent patches)
!                                         when > 4 int. cells meet at the edge 
!        type of degenerative corner = 0 (none)
!                                    = 1 (corner within edge)
!                                         when 6 int. cells meet at the corner 
!                                    = 2 (corner within patch)
!                                         when 4 int. cells meet at the corner 
!                                    < 0 (detached from adjacent edges)
!                                         when >8 int.cells meet at the corner 
!
!        Assignment of negative level (disconnected edge/corners) happens 
!        outside this routine.
!
!******************************************************************************

SUBROUTINE RFLO_FindDegeneratCell( regions,handleCorn,icount,iReg,iLev, &
                                   iRegTest,icellTest,iTest,jTest,kTest )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetEdgeCellsIndices
  IMPLICIT NONE

! ... parameters
  INTEGER :: handleCorn, icount, iReg, iLev
  INTEGER :: iRegTest, icellTest, iTest, jTest, kTest

  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iPatch, iedge, i, j, k

! ... local variables
  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER  :: level
  TYPE(t_patch), POINTER  :: patch, patchSrc

  INTEGER :: iRegSrc, iPatchSrc, bcType, lbound, jReg, jcell
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, ijk, nshift
  INTEGER :: ibegSrc, iendSrc, jbegSrc, jendSrc, kbegSrc, kendSrc


!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_FindDegeratCell',&
       'RFLO_ModDegenerateCornEdge.F90' )

! get pointers ----------------------------------------------------------------

  level => regions(iReg)%levels(iLev)

! check if input edge or corner within patch dummy layers ---------------------

  DO iPatch=1,regions(iReg)%nPatches
    patch => regions(iReg)%levels(iLev)%patches(iPatch)
    bcType = patch%bcType
    lbound = patch%lbound

    IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
        (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
        (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN

      iRegSrc   = patch%srcRegion

      IF (iRegTest == iRegSrc) THEN
        iPatchSrc = patch%srcPatch
        patchSrc => regions(iRegSrc)%levels(iLev)%patches(iPatchSrc)

        CALL RFLO_GetPatchIndices( regions(iRegSrc),patchSrc,iLev,ibegSrc, &
                                   iendSrc,jbegSrc,jendSrc,kbegSrc,kendSrc )
        nshift = regions(iRegSrc)%nDumCells-1
        IF (patchSrc%lbound==1) iendSrc = iendSrc+nshift
        IF (patchSrc%lbound==2) ibegSrc = ibegSrc-nshift
        IF (patchSrc%lbound==3) jendSrc = jendSrc+nshift
        IF (patchSrc%lbound==4) jbegSrc = jbegSrc-nshift
        IF (patchSrc%lbound==5) kendSrc = kendSrc+nshift
        IF (patchSrc%lbound==6) kbegSrc = kbegSrc-nshift

        IF ((iTest>=ibegSrc .AND. iTest<=iendSrc) .AND. &
            (jTest>=jbegSrc .AND. jTest<=jendSrc) .AND. &
            (kTest>=kbegSrc .AND. kTest<=kendSrc)) THEN
          IF (handleCorn==0) THEN

! --------- edge within patch
            level%edgeCells(icount)%degenrt = DEGENERAT_EDGE_IN_PATCH
          ELSE

! --------- corner within patch
            level%cornerCells(icount)%degenrt = DEGENERAT_CORN_IN_PATCH
          ENDIF  ! handle
          global%degenrtEc = .TRUE.

        ENDIF    ! i,j,kTest
      ENDIF      ! iRegTest
    ENDIF        ! bcType
  ENDDO          ! iPatch

! check if input corner within an edge ----------------------------------------

  IF (handleCorn==1) THEN

    DO iedge=1,12
      IF (level%edgeCells(iedge)%interact) THEN
        CALL RFLO_GetEdgeCellsIndices( regions(iReg),iLev,iedge, &
                                       ibeg,iend,jbeg,jend,kbeg,kend )
        ijk = 0
        DO k=kbeg,kend
          DO j=jbeg,jend
            DO i=ibeg,iend
              ijk = ijk + 1
              jcell = level%edgeCells(iedge)%cells(ijk)%srcCell
              jReg  = level%edgeCells(iedge)%cells(ijk)%srcRegion
              IF (iRegTest==jReg .AND. icellTest==jcell) THEN

! ------------- corner within edge
                level%cornerCells(icount)%degenrt = DEGENERAT_CORN_IN_EDGE
                global%degenrtEc = .TRUE.
                GOTO 888
              ENDIF
            ENDDO  ! i
          ENDDO    ! j
        ENDDO      ! k
888     CONTINUE

      ENDIF      ! interact
    ENDDO        ! iedge
  ENDIF          ! handleCorn

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_FindDegeneratCell


!******************************************************************************
!
! Purpose: flag degenerate vertices.
!
! Description: grid%ijkDgen = 0 for normal vertices and 1 for degenerate 
!
! Input: regions = dimensions and topology of all regions.
!
! Output: regions = grid%ijkDgen
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_MarkDegeneratVert( regions )

  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes, RFLO_GetNodeOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iLev, iedge, icorner, i, j, k

! ... local variables
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, iNOff, ijNOff
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, ijk, errorFlag

  INTEGER, POINTER :: idgen(:)
  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER  :: level

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MarkDegeneratVert',&
       'RFLO_ModDegenerateCornEdge.F90' )

! find degenerate vertices ====================================================

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

      DO iLev=1,regions(iReg)%nGridLevels
        level => regions(iReg)%levels(iLev)
        idgen => regions(iReg)%levels(iLev)%grid%ijkDgen
        idgen =  0

        CALL RFLO_GetDimensPhysNodes( regions(iReg),iLev,ipnbeg,ipnend, &
                                    jpnbeg,jpnend,kpnbeg,kpnend )
        CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )

! ----- edge vertices ---------------------------------------------------------

        DO iedge=1,12

          SELECT CASE (iedge)
          CASE (1)
            ibeg = ipnbeg
            iend = ipnbeg
            jbeg = jpnbeg
            jend = jpnbeg
            kbeg = kpnbeg
            kend = kpnend
          CASE (2)
            ibeg = ipnbeg
            iend = ipnbeg
            jbeg = jpnbeg
            jend = jpnend
            kbeg = kpnend
            kend = kpnend
          CASE (3)
            ibeg = ipnbeg
            iend = ipnbeg
            jbeg = jpnend
            jend = jpnend
            kbeg = kpnbeg
            kend = kpnend
          CASE (4)
            ibeg = ipnbeg
            iend = ipnbeg
            jbeg = jpnbeg
            jend = jpnend
            kbeg = kpnbeg
            kend = kpnbeg
          CASE (5)
            ibeg = ipnend
            iend = ipnend
            jbeg = jpnbeg
            jend = jpnbeg
            kbeg = kpnbeg
            kend = kpnend
          CASE (6)
            ibeg = ipnend
            iend = ipnend
            jbeg = jpnbeg
            jend = jpnend
            kbeg = kpnend
            kend = kpnend
          CASE (7)
            ibeg = ipnend
            iend = ipnend
            jbeg = jpnend
            jend = jpnend
            kbeg = kpnbeg
            kend = kpnend
          CASE (8)
            ibeg = ipnend
            iend = ipnend
            jbeg = jpnbeg
            jend = jpnend
            kbeg = kpnbeg
            kend = kpnbeg
          CASE (9)
            ibeg = ipnbeg
            iend = ipnend
            jbeg = jpnbeg
            jend = jpnbeg
            kbeg = kpnbeg
            kend = kpnbeg
          CASE (10)
            ibeg = ipnbeg
            iend = ipnend
            jbeg = jpnbeg
            jend = jpnbeg
            kbeg = kpnend
            kend = kpnend
          CASE (11)
            ibeg = ipnbeg
            iend = ipnend
            jbeg = jpnend
            jend = jpnend
            kbeg = kpnend
            kend = kpnend
          CASE (12)
            ibeg = ipnbeg
            iend = ipnend
            jbeg = jpnend
            jend = jpnend
            kbeg = kpnbeg
            kend = kpnbeg
          END SELECT

          IF (level%edgeCells(iedge)%degenrt /= DEGENERAT_NONE) THEN
            DO k=kbeg,kend
              DO j=jbeg,jend
                DO i=ibeg,iend
                  ijk = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
                  idgen(ijk) = 1
                ENDDO
              ENDDO
            ENDDO
          ENDIF  ! degenrt
        ENDDO    ! iedge

! ----- corner cells ----------------------------------------------------------

        DO icorner=1,8
          SELECT CASE (icorner)
          CASE (1)
            ibeg = ipnbeg
            iend = ipnbeg
            jbeg = jpnbeg
            jend = jpnbeg
            kbeg = kpnbeg
            kend = kpnbeg
          CASE (2)
            ibeg = ipnbeg
            iend = ipnbeg
            jbeg = jpnbeg
            jend = jpnbeg
            kbeg = kpnend
            kend = kpnend
          CASE (3)
            ibeg = ipnbeg
            iend = ipnbeg
            jbeg = jpnend
            jend = jpnend
            kbeg = kpnend
            kend = kpnend
          CASE (4)
            ibeg = ipnbeg
            iend = ipnbeg
            jbeg = jpnend
            jend = jpnend
            kbeg = kpnbeg
            kend = kpnbeg
          CASE (5)
            ibeg = ipnend
            iend = ipnend
            jbeg = jpnbeg
            jend = jpnbeg
            kbeg = kpnbeg
            kend = kpnbeg
          CASE (6)
            ibeg = ipnend
            iend = ipnend
            jbeg = jpnbeg
            jend = jpnbeg
            kbeg = kpnend
            kend = kpnend
          CASE (7)
            ibeg = ipnend
            iend = ipnend
            jbeg = jpnend
            jend = jpnend
            kbeg = kpnend
            kend = kpnend
          CASE (8)
            ibeg = ipnend
            iend = ipnend
            jbeg = jpnend
            jend = jpnend
            kbeg = kpnbeg
            kend = kpnbeg
          END SELECT

          IF (level%cornerCells(icorner)%degenrt /= DEGENERAT_NONE) THEN
            DO k=kbeg,kend
              DO j=jbeg,jend
                DO i=ibeg,iend
                  ijk = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
                  idgen(ijk) = 1
                ENDDO
              ENDDO
            ENDDO
          ENDIF  ! degenrt
        ENDDO    ! icorner

      ENDDO      ! iLev
    ENDIF        ! myProcid
  ENDDO          ! iReg

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_MarkDegeneratVert


! ******************************************************************************
!
! Purpose: write degenerated edges and corners (if exist) in each region 
!          into a file
!
! Description: only MASTERPROC write the information as all edge and corner
!              data are available in this processor.
!
! Notes: none.
!
! ******************************************************************************

SUBROUTINE RFLO_WriteDegeneratEC( regions )

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iedge, icorner
   
! ... local variables
  CHARACTER(CHRLEN) :: iFileName
  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER  :: level

  INTEGER :: iLev, errorFlag, iFile
  
! ******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLO_WriteDegeneratEC',&
       'RFLO_ModDegenerateCornEdge.F90')

! start ------------------------------------------------------------------------

  iFile = IF_DEGENRT

  WRITE(iFileName,'(A)') &
        TRIM(global%outDir)//TRIM(global%casename)//'.degec'

  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN",IOSTAT=errorFlag)   
  global%error = errorFlag        

  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,&
         __LINE__,iFileName)
  END IF ! global%error

! header and general information ---------------------------------------------

  WRITE(iFile,1005)' ROCFLO degenerated edges and corners'  

  WRITE(iFile,1000)

! print degenerated edge and corner per region (if exist) --------------------

  DO iReg=1,global%nRegions

    iLev  =  regions(iReg)%currLevel
    level => regions(iReg)%levels(iLev)

    WRITE(iFile,1025) iReg
    WRITE(iFile,1005) '  Edge #       Degeneration-type' 
      
    DO iedge = 1,12
      IF (level%edgeCells(iedge)%degenrt /= DEGENERAT_NONE) &
        WRITE(iFile,1015) iedge, level%edgeCells(iedge)%degenrt
    ENDDO          

    WRITE(iFile,1005) '  Corner #     Degeneration-type' 
      
    DO icorner = 1,8
      IF (level%cornerCells(icorner)%degenrt /= DEGENERAT_NONE) &
        WRITE(iFile,1015) icorner, level%cornerCells(icorner)%degenrt
    ENDDO          

    WRITE(iFile,1000)
  ENDDO  ! iReg

! print legend explaining meaning of degeration values ----------------------

  WRITE(iFile,1005)' Edge degeneration type:'  
  WRITE(iFile,1030)'    1 = EDGE IN PATCH (< 4 interior cells meet at the edge)'  
  WRITE(iFile,1030)'   -1 = EDGE DETACHED (> 4 interior cells meet at the edge)'  

  WRITE(iFile,1005)' Corner degeneration type:'
  WRITE(iFile,1030)'    1 = CORNER IN EDGE  (  6 interior cells meet at the corner)'  
  WRITE(iFile,1030)'    2 = CORNER IN PATCH (  4 interior cells meet at the corner)'  
  WRITE(iFile,1030)'   -1 = CORNER DETACHED (> 8 interior cells meet at the corner)'  

! close file -------------------------------------------------------------------

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag      
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,&
         __LINE__,iFileName)
  END IF

! finalize ---------------------------------------------------------------------

  CALL DeregisterFunction( global )

1000 FORMAT(/,1X,40('-'))
1005 FORMAT(/,A)
1015 FORMAT(I8,8X,I8)
1025 FORMAT('   Region ',I6,':')
1030 FORMAT(A)
1035 FORMAT(/,1X,40('-'),/)
  
END SUBROUTINE RFLO_WriteDegeneratEC

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLO_ModDegenerateCornEdge

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModDegenerateCornEdge.F90,v $
! Revision 1.3  2008/12/06 08:44:15  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2005/11/11 07:24:58  wasistho
! initial import RFLO_ModDegenerateCornEdge
!
!
! ******************************************************************************









