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
! Purpose: driver to communicate buffer data of RFLO metrics, face centroids
!          and face normals, for adjacent regions on different processors.
!
! Description: none.
!
! Input: regions = data of all regions.
!
! Output: buffer data for metrics.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_RFLO_SendRecvMetrics.F90,v 1.4 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLO_SendRecvMetrics( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region, t_level, t_dCellTransf
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag
  USE PLAG_ModInterfaces, ONLY: PLAG_RFLO_RecvMetrics,      &
                                PLAG_RFLO_SendMetrics,      &  
                                PLAG_RFLO_ClearSendRequests
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iCorner, iEdge, ijk, iLev, ir, iReg 

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  
  INTEGER :: errorFlag
  INTEGER :: nCorners, nEdges  
  INTEGER :: nFaces, nFaceCentroidSize, nFaceNormalSize, nMetricsBuff

  TYPE(t_level),       POINTER :: pLevel
  TYPE(t_dCellTransf), POINTER :: pSendEcCell, pRecvEcCell
  TYPE(t_plag),        POINTER :: pPlag
  TYPE(t_global),      POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLO_SendRecvMetrics.F90,v $ $Revision: 1.4 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PLAG_RFLO_SendRecvMetrics',&
  'PLAG_RFLO_SendRecvMetrics.F90' )

! Set dimensions --------------------------------------------------------------

  nCorners = 8
  nEdges   = 12
  nFaces   = 6

  nFaceCentroidSize = ZCOORD*KCOORD
  nFaceNormalSize   =      3*KCOORD
  nMetricsBuff      = nFaces*(nFaceCentroidSize + nFaceNormalSize)

! Allocate communicaton buffers ===============================================

  DO iReg = 1, global%nRegions  
    IF ( regions(iReg)%procid==global%myProcid .AND. &   ! region active and
         regions(iReg)%active==ACTIVE          .AND. &   ! on my processor
         global%nProcAlloc > 1                 ) THEN    ! only if multiple processors

      DO iLev=1,regions(iReg)%nGridLevels
        pPlag => regions(iReg)%levels(iLev)%plag
        pPlag%nRequestsMetrics = 0
        
        DO ir=1,global%nRegions
          pSendEcCell => regions(iReg)%levels(iLev)%sendEcCells(ir)
          pRecvEcCell => regions(iReg)%levels(iLev)%recvEcCells(ir)
          
          IF ( pSendEcCell%nCells > 0 ) THEN
            pPlag%nRequestsMetrics = pPlag%nRequestsMetrics+1
            pSendEcCell%iRequestMetrics = pPlag%nRequestsMetrics
            ALLOCATE( pSendEcCell%buffMetrics(pSendEcCell%nCells*nMetricsBuff), &
                      stat=errorFlag )
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) &
              CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
          ENDIF ! pSendEcCell%nCells
          
          IF ( pRecvEcCell%nCells > 0 ) THEN
            pRecvEcCell%iRequestMetrics = -999999
            ALLOCATE( pRecvEcCell%buffMetrics(pRecvEcCell%nCells*nMetricsBuff), &
                      stat=errorFlag )
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) &
            CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
          ENDIF ! pRecvEcCell%nCell 
        ENDDO   ! ir

        ALLOCATE(pPlag%requestsMetrics(pPlag%nRequestsMetrics), stat=errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) &
          CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      ENDDO     ! iLev

    ENDIF       ! regions
  ENDDO         ! iReg

! synchronize through an MPI barrier ------------------------------------------

#ifdef MPI
  CALL MPI_Barrier( global%mpiComm,global%mpierr )
  IF ( global%mpierr /= ERR_NONE ) &
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! communicate buffer ----------------------------------------------------------
                           
! send buffer------------------------------------------------------------------  

  DO iReg = 1, global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor        
      CALL PLAG_RFLO_SendMetrics(regions, iReg)    
    ENDIF       ! regions    
  ENDDO         ! iReg

! receive buffer --------------------------------------------------------------

  DO iReg = 1, global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor        
      CALL PLAG_RFLO_RecvMetrics(regions, iReg)
    ENDIF ! regions
  ENDDO ! iReg 

! wait for data being received by other processors ----------------------------

  DO iReg = 1, global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
      CALL PLAG_RFLO_ClearSendRequests(regions, iReg)
    ENDIF ! regions  
  ENDDO ! iReg

! synchronize through an MPI barrier ------------------------------------------

#ifdef MPI
  CALL MPI_Barrier( global%mpiComm,global%mpierr )
  IF ( global%mpierr /= ERR_NONE ) &
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! deallocate communicaton buffers ---------------------------------------------

  DO iReg = 1, global%nRegions  
    IF ( regions(iReg)%procid==global%myProcid .AND. &   ! region active and
         regions(iReg)%active==ACTIVE          .AND. &   ! on my processor
         global%nProcAlloc > 1                 ) THEN    ! only if multiple processors

      DO iLev=1,regions(iReg)%nGridLevels
        pPlag => regions(iReg)%levels(iLev)%plag

        DO ir=1,global%nRegions
          pSendEcCell => regions(iReg)%levels(iLev)%sendEcCells(ir)
          pRecvEcCell => regions(iReg)%levels(iLev)%recvEcCells(ir)
          
          IF ( pSendEcCell%nCells > 0 ) THEN
            DEALLOCATE( pSendEcCell%buffMetrics, stat=errorFlag )
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) &
              CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
          ENDIF ! pSendEcCell%nCells
          
          IF ( pRecvEcCell%nCells > 0 ) THEN
            DEALLOCATE( pRecvEcCell%buffMetrics, stat=errorFlag )
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) &
            CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )
          ENDIF ! pRecvEcCell%nCell 
        ENDDO   ! ir

        DEALLOCATE(pPlag%requestsMetrics, stat=errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) &
          CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

      ENDDO     ! iLev

    ENDIF       ! regions
  ENDDO         ! iReg

! synchronize through an MPI barrier ------------------------------------------

#ifdef MPI
  CALL MPI_Barrier( global%mpiComm,global%mpierr )
  IF ( global%mpierr /= ERR_NONE ) &
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif
  
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLO_SendRecvMetrics

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLO_SendRecvMetrics.F90,v $
! Revision 1.4  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 20:58:14  fnajjar
! Initial revision after changing case
!
! Revision 1.1  2004/01/15 21:16:48  fnajjar
! Initial import for corner-edge cell metrics
!
!******************************************************************************







