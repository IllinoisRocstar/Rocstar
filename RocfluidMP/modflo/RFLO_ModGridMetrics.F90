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
! Purpose: Suite for grid metrics related routines.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModGridMetrics.F90,v 1.7 2008/12/06 08:44:16 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModGridMetrics

  USE ModGlobal, ONLY    : t_global 
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY      : t_grid
  USE ModBndPatch, ONLY  : t_patch
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLO_ArcLengthPatch, &
            RFLO_CalcGridMetrics, &
            RFLO_GridQualityGlobal

! private :
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModGridMetrics.F90,v $ $Revision: 1.7 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS


!******************************************************************************
!
! Purpose: calculate approximate arclengths for every grid line
!          between two oposite patch boundaries (on the finest
!          grid only).
!
! Description: none.
!
! Input: region = grid dimensions
!        patch  = current patch
!        xyzRef = reference coordinates.
!
! Output: arcLen1 = arclength in first  coordinate direction of patch
!         arcLen2 = arclength in second coordinate direction of patch
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_ArcLengthPatch( region,patch,xyzRef )

  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetPatchIndicesNodes

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters

  TYPE(t_region)         :: region
  TYPE(t_patch), POINTER :: patch
  REAL(RFREAL), POINTER  :: xyzRef(:,:)

! ... loop variables
  INTEGER :: l1, l2

! ... local variables
  INTEGER :: iLev, lbound, ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: l1b, l1e, l2b, l2e, lc, ijkN, ijkN1, ijkN2, iNOff, ijNOff
  INTEGER :: k1, k2, switch(6,5)

  REAL(RFREAL), POINTER :: arcLen1(:), arclen2(:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_ArcLengthPatch',&
       'RFLO_ModGridMetrics.F90' )

! get dimensions and pointers -------------------------------------------------

  lbound = patch%lbound
  iLev   = 1

  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                  ibeg,iend,jbeg,jend,kbeg,kend )

  arclen1 => patch%arclen1
  arclen2 => patch%arclen2

! set boundary switch ---------------------------------------------------------
! switch(:,1-2) = first/last index in l1-direction
! switch(:,3-4) = first/last index in l2-direction
! switch(:,  5) = constant index

  switch(1,:) = (/jbeg, jend, kbeg, kend, ibeg/)
  switch(2,:) = (/jbeg, jend, kbeg, kend, iend/)
  switch(3,:) = (/kbeg, kend, ibeg, iend, jbeg/)
  switch(4,:) = (/kbeg, kend, ibeg, iend, jend/)
  switch(5,:) = (/ibeg, iend, jbeg, jend, kbeg/)
  switch(6,:) = (/ibeg, iend, jbeg, jend, kend/)

  l1b = switch(lbound,1)
  l1e = switch(lbound,2)
  l2b = switch(lbound,3)
  l2e = switch(lbound,4)
  lc  = switch(lbound,5)

! compute arclengths of current patch ------------------------------------------

  arclen1(:) = 0._RFREAL
  arclen2(:) = 0._RFREAL

  DO l2=l2b+1,l2e
    k2 = l2-l2b+1
    DO l1=l1b+1,l1e
      k1 = l1-l1b+1
      IF (lbound==1 .OR. lbound==2) THEN
        ijkN      = IndIJK(lc,l1    ,l2    ,iNOff,ijNOff)
        ijkN1     = IndIJK(lc,l1-1  ,l2    ,iNOff,ijNOff)
        ijkN2     = IndIJK(lc,l1    ,l2-1  ,iNOff,ijNOff)
      ELSEIF (lbound==3 .OR. lbound==4) THEN
        ijkN      = IndIJK(l2    ,lc,l1    ,iNOff,ijNOff)
        ijkN1     = IndIJK(l2    ,lc,l1-1  ,iNOff,ijNOff)
        ijkN2     = IndIJK(l2-1  ,lc,l1    ,iNOff,ijNOff)
      ELSEIF (lbound==5 .OR. lbound==6) THEN
        ijkN      = IndIJK(l1    ,l2    ,lc,iNOff,ijNOff)
        ijkN1     = IndIJK(l1-1  ,l2    ,lc,iNOff,ijNOff)
        ijkN2     = IndIJK(l1    ,l2-1  ,lc,iNOff,ijNOff)
      ENDIF
      arclen1(k2) = arclen1(k2) + &
            SQRT((xyzRef(XCOORD,ijkN)-xyzRef(XCOORD,ijkN1))**2 + &
                 (xyzRef(YCOORD,ijkN)-xyzRef(YCOORD,ijkN1))**2 + &
                 (xyzRef(ZCOORD,ijkN)-xyzRef(ZCOORD,ijkN1))**2)
      arclen2(k1) = arclen2(k1) + &
            SQRT((xyzRef(XCOORD,ijkN)-xyzRef(XCOORD,ijkN2))**2 + &
                 (xyzRef(YCOORD,ijkN)-xyzRef(YCOORD,ijkN2))**2 + &
                 (xyzRef(ZCOORD,ijkN)-xyzRef(ZCOORD,ijkN2))**2)
    ENDDO ! l1
  ENDDO   ! l2

  DO l2=l2b,l2b
    DO l1=l1b+1,l1e
      IF (lbound==1 .OR. lbound==2) THEN
        ijkN      = IndIJK(lc,l1    ,l2    ,iNOff,ijNOff)
        ijkN1     = IndIJK(lc,l1-1  ,l2    ,iNOff,ijNOff)
      ELSEIF (lbound==3 .OR. lbound==4) THEN
        ijkN      = IndIJK(l2    ,lc,l1    ,iNOff,ijNOff)
        ijkN1     = IndIJK(l2    ,lc,l1-1  ,iNOff,ijNOff)
      ELSEIF (lbound==5 .OR. lbound==6) THEN
        ijkN      = IndIJK(l1    ,l2    ,lc,iNOff,ijNOff)
        ijkN1     = IndIJK(l1-1  ,l2    ,lc,iNOff,ijNOff)
      ENDIF
      arclen1(1) = arclen1(1) + &
            SQRT((xyzRef(XCOORD,ijkN)-xyzRef(XCOORD,ijkN1))**2 + &
                 (xyzRef(YCOORD,ijkN)-xyzRef(YCOORD,ijkN1))**2 + &
                 (xyzRef(ZCOORD,ijkN)-xyzRef(ZCOORD,ijkN1))**2)
    ENDDO ! l1
  ENDDO   ! l2

  DO l1=l1b,l1b
    DO l2=l2b+1,l2e
      IF (lbound==1 .OR. lbound==2) THEN
        ijkN      = IndIJK(lc,l1    ,l2    ,iNOff,ijNOff)
        ijkN2     = IndIJK(lc,l1    ,l2-1  ,iNOff,ijNOff)
      ELSEIF (lbound==3 .OR. lbound==4) THEN
        ijkN      = IndIJK(l2    ,lc,l1    ,iNOff,ijNOff)
        ijkN2     = IndIJK(l2-1  ,lc,l1    ,iNOff,ijNOff)
      ELSEIF (lbound==5 .OR. lbound==6) THEN
        ijkN      = IndIJK(l1    ,l2    ,lc,iNOff,ijNOff)
        ijkN2     = IndIJK(l1    ,l2-1  ,lc,iNOff,ijNOff)
      ENDIF
      arclen2(1) = arclen2(1) + &
            SQRT((xyzRef(XCOORD,ijkN)-xyzRef(XCOORD,ijkN2))**2 + &
                 (xyzRef(YCOORD,ijkN)-xyzRef(YCOORD,ijkN2))**2 + &
                 (xyzRef(ZCOORD,ijkN)-xyzRef(ZCOORD,ijkN2))**2)
    ENDDO ! l2
  ENDDO   ! l1

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ArcLengthPatch


!******************************************************************************
!
! Purpose: calculate face vectors, volumes, cell centroids (optionally),
!          and cell-to-face averaging coefficients
!
! Description: none.
!
! Input: regions%grid = dimensions, grid coordinates.
!
! Output: regions%grid = face vectors (si,sj,sk), volumes (vol),
!                        cell centroids (cofg), avg coeffs. (avgCoI,J,K).
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_CalcGridMetrics( regions )

  USE ModInterfaces, ONLY : RFLO_CalcFaceVectors, RFLO_CalcControlVolumes, &
             RFLO_CalcCellCentroids, RFLO_CalcFaceCentroids, &
             RFLO_InitAvgCoeffs, RFLO_C2fAvgCoeffs, RFLO_C2eAvgCoeffs, &
             RFLO_CheckMetrics
  IMPLICIT NONE

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  REAL(RFREAL) :: skewMin
  TYPE (t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global
  CALL RegisterFunction( global,'RFLO_CalcGridMetrics',&
       'RFLO_ModGridMetrics.F90' )

! initialization some parameters ----------------------------------------------

  global%skewness = 1._RFREAL

! loop over all regions -------------------------------------------------------

  DO iReg=1,regions(1)%global%nRegions
    IF (regions(iReg)%procid==regions(iReg)%global%myProcid & ! reg active and
        .AND. regions(iReg)%active==ACTIVE) THEN              ! on my processor

      CALL RFLO_CalcFaceVectors( regions(iReg) )

      CALL RFLO_CalcControlVolumes( regions(iReg) )

      CALL RFLO_CalcCellCentroids( regions(iReg) )

      CALL RFLO_CalcFaceCentroids( regions(iReg) )

      CALL RFLO_InitAvgCoeffs( regions(iReg) )

      IF (regions(iReg)%mixtInput%faceEdgeAvg==FE_AVG_LINEAR) &
        CALL RFLO_C2fAvgCoeffs( regions(iReg) )

      CALL RFLO_C2eAvgCoeffs( regions(iReg) )

! --- check metrics

      CALL RFLO_CheckMetrics( iReg,regions(iReg) )

    ENDIF     ! region on this processor and active
  ENDDO       ! iReg

! global grid quality measure -------------------------------------------------

  CALL RFLO_GridQualityGlobal( regions )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_CalcGridMetrics


!******************************************************************************
!
! Purpose: global reduction of grid quality measures
!
! Description: none.
!
! Input: regions%grid, regions%global = grid quality data
!
! Output: global%skewness, etc = global skewness, etc
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_GridQualityGlobal( regions )

  IMPLICIT NONE

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  REAL(RFREAL) :: skewMin, volMin
  TYPE (t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global
  CALL RegisterFunction( global,'RFLO_GridQualityGlobal',&
       'RFLO_ModGridMetrics.F90' )

! global skewness and minVol --------------------------------------------------

  skewMin = global%skewness
  volMin  = global%minVol

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &     ! region active and
        regions(iReg)%active==ACTIVE) THEN                ! on my processor
      skewMin = MIN( skewMin,regions(iReg)%levels(1)%grid%skewness )
      volMin  = MIN( volMin, regions(iReg)%levels(1)%grid%minVol   )
    ENDIF   ! active
  ENDDO     ! iReg

  global%skewness = skewMin
  global%minVol   = volMin

#ifdef MPI
  CALL MPI_Allreduce( global%skewness,skewMin,1,MPI_RFREAL,MPI_MIN, &
                      global%mpiComm,global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
       __LINE__ )
  global%skewness = skewMin

  CALL MPI_Allreduce( global%minVol,volMin,1,MPI_RFREAL,MPI_MIN, &
                      global%mpiComm,global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
       __LINE__ )
  global%minVol = volMin
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_GridQualityGlobal


! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModGridMetrics


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModGridMetrics.F90,v $
! Revision 1.7  2008/12/06 08:44:16  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/03/24 00:53:51  wasistho
! fixed loop indexing in patch%arclen1,2
!
! Revision 1.4  2006/03/18 11:01:23  wasistho
! moved some routines to ModGridRegionShape
!
! Revision 1.3  2006/03/18 08:17:01  wasistho
! added arcLengthPatch
!
! Revision 1.2  2006/03/15 06:37:46  wasistho
! added region and global skewness
!
! Revision 1.1  2006/03/04 04:36:41  wasistho
! initial import RFLO_ModGridMetrics
!
!
!
! ******************************************************************************









