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
! Purpose: update face centroids and face normals for RFLO.
!
! Description: none.
!
! Input: regions = data of all regions
!
! Output: regions%levels%plag%fc       = face centroids
!         regions%levels%plag%si,sj,sk = face vectors
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_RFLO_SetMetrics.F90,v 1.3 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLO_SetMetrics( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PLAG_ModInterfaces, ONLY : PLAG_CalcFaceCentroids,    &
                                 PLAG_CopyFaceVectors,      &
                                 PLAG_CECellsFaceCentroids, &
                                 PLAG_CECellsFaceVectors,   &
                                 PLAG_RFLO_FindGridMapping, &
                                 PLAG_RFLO_SendRecvMetrics

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLO_SetMetrics.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PLAG_RFLO_SetMetrics',&
  'PLAG_RFLO_SetMetrics.F90' )

! check if module is active in any region======================================

  IF (.NOT. global%plagUsed) GOTO 999

! Determine the transfomation mapping matrix for corner and edge cells --------

  CALL PLAG_RFLO_FindGridMapping( regions )

! Compute face centroids and normals for corner and edge cells ----------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE ) THEN             ! on my processor

!#ifdef PLAG_DEBUG
      IF ( iReg == 1 ) WRITE(*,*) 'Entering PLAG_CalcFaceCentroids: iReg=',iReg
!#endif
      CALL PLAG_CalcFaceCentroids( regions(iReg) )

!#ifdef PLAG_DEBUG
      IF ( iReg == 1 ) WRITE(*,*) 'Entering PLAG_CopyFaceVectors: iReg=',iReg
!#endif
      CALL PLAG_CopyFaceVectors( regions(iReg) )

    ENDIF ! regions
  ENDDO ! iReg

! Load face centroids for corner and edge cells --------------------------------
! On-processor infrastructure --------------------------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

!#ifdef PLAG_DEBUG
      IF ( iReg == 1 ) WRITE(*,*) 'Entering PLAG_CECellsFaceCentroids: iReg=',iReg
!#endif
      CALL PLAG_CECellsFaceCentroids(   regions, iReg )

!#ifdef PLAG_DEBUG
      IF ( iReg == 1 ) WRITE(*,*) 'Entering PLAG_CECellsFaceVectors: iReg=',iReg
!#endif
      CALL PLAG_CECellsFaceVectors(   regions, iReg )

    ENDIF ! regions
  ENDDO ! iReg

! Communicate buffer data for off-processor regions ---------------------------

  CALL PLAG_RFLO_SendRecvMetrics( regions )

! Reevaluate face centroids and normals to deal with degeneracy case ----------
! Note: Degeneracy case arises when odd number of regions intersect -----------
!       and the node has three surrounding cells as opposed to four -----------
!       It also applies when two regions have communicating faces -------------
!       yet edge-corner search algorithm sets them as edge-corner regions -----

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE ) THEN             ! on my processor
      CALL PLAG_CalcFaceCentroids( regions(iReg) )
      CALL PLAG_CopyFaceVectors(   regions(iReg) )
    ENDIF ! regions
  ENDDO ! iReg

! finalize ====================================================================

999 CONTINUE
  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLO_SetMetrics

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLO_SetMetrics.F90,v $
! Revision 1.3  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:58:15  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.4  2004/02/11 23:18:23  fnajjar
! Added a second call to routines for degeneracy case
!
! Revision 1.3  2004/02/10 21:24:27  fnajjar
! Added capability to determine index mapping between corner-edge regions
!
! Revision 1.2  2004/01/15 21:11:55  fnajjar
! Added MPI-based kernel for corner-edge cell metrics
!
! Revision 1.1  2003/11/12 21:35:58  fnajjar
! Initial import of Corner-Edge cells Infrastructure
!
!******************************************************************************







