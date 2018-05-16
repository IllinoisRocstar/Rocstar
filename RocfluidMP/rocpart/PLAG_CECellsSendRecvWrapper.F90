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
! Purpose: driver to communicate PLAG buffer data in corner and edge cells
!          for adjacent regions on different processors.
!
! Description: none.
!
! Input: 
!   regions = data of all regions.
!
! Output: buffer data for metrics.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_CECellsSendRecvWrapper.F90,v 1.5 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CECellsSendRecvWrapper( regions )

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModDataStruct, ONLY : t_region, t_level, t_dCellTransf
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag
  USE ModMPI

  USE PLAG_ModParameters
 
  USE PLAG_ModInterfaces, ONLY: PLAG_CECellsRecvData,          &
                                PLAG_CECellsRecvSize,          &
                                PLAG_CECellsSendData,          &
                                PLAG_CECellsSendSize,          &
                                PLAG_CECellsClearRequestsData, &  
                                PLAG_CECellsClearRequestsSize

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: errorFlag,ijk,iLev,ir,iReg,nAiv,nArv,nBuffRecvI,nBuffRecvR,&
             nBuffSendI,nBuffSendR,nCv,nDimBuffI,nDimBuffR,nDv,nTv
  INTEGER :: iBuffI,iBuffR
             
  TYPE(t_level),       POINTER :: pLevel
  TYPE(t_dCellTransf), POINTER :: pSendEcCell, pRecvEcCell
  TYPE(t_plag),        POINTER :: pPlag
  TYPE(t_global),      POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CECellsSendRecvWrapper.F90,v $ $Revision: 1.5 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PLAG_CECellsSendRecv',&
  'PLAG_CECellsSendRecvWrapper.F90' )

! ******************************************************************************
! Set infrastructure for MPI-based communication
! ******************************************************************************

  DO iReg = 1, global%nRegions  
    IF ( regions(iReg)%procid==global%myProcid .AND. &   ! region active and
         regions(iReg)%active==ACTIVE          .AND. &   ! on my processor
         global%plagUsed    .EQV. .TRUE.       .AND. &   ! plag active
         global%nProcAlloc > 1                  ) THEN   ! only if multiple processors

      DO iLev=1,regions(iReg)%nGridLevels
        pPlag => regions(iReg)%levels(iLev)%plag

! =============================================================================  
!       Set dimensions
! =============================================================================  
         
        pPlag%nRequestsCECells  = 0

        nAiv = pPlag%nAiv
        nArv = pPlag%nArv
   
        nCv  = pPlag%nCv
        nDv  = pPlag%nDv
        nTv  = pPlag%nTv  

        nDimBuffI = nAiv  
        nDimBuffR = 2*nArv +4*nCv

! =============================================================================  
!       Loop over all regions 
! =============================================================================         

        DO ir=1,global%nRegions
          pSendEcCell => regions(iReg)%levels(iLev)%sendEcCells(ir)
          pRecvEcCell => regions(iReg)%levels(iLev)%recvEcCells(ir)

! ----------------------------------------------------------------------------- 
!         Initialize buffer sizes
! ----------------------------------------------------------------------------- 
          
          IF ( pSendEcCell%nCells > 0 ) pSendEcCell%nBuffSizePlag = 0 

          IF ( pRecvEcCell%nCells > 0 ) pRecvEcCell%nBuffSizePlag = 0
 
! ----------------------------------------------------------------------------- 
!         Set request sizes
! ----------------------------------------------------------------------------- 
          
          IF ( pSendEcCell%nCells > 0 ) &
            pPlag%nRequestsCECells = pPlag%nRequestsCECells +1 

          pSendEcCell%iRequestPlag = pPlag%nRequestsCECells

          IF ( pRecvEcCell%nCells > 0 ) &
            pRecvEcCell%iRequestPlag = -999999
 
        ENDDO   ! ir

! -----------------------------------------------------------------------------
!       Allocate requests 
! -----------------------------------------------------------------------------

        ALLOCATE( pPlag%requestsCECells(pPlag%nRequestsCECells), &
                  STAT=errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) &
          CALL ErrorStop( global,ERR_ALLOCATE,__LINE__, &
                          'pPlag%requestsCECells' )

        ALLOCATE( pPlag%requestsCECellsI(pPlag%nRequestsCECells), &
                  STAT=errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) &
          CALL ErrorStop( global,ERR_ALLOCATE,__LINE__, &
                          'pPlag%requestsCECellsI' )

        ALLOCATE( pPlag%requestsCECellsR(pPlag%nRequestsCECells), &
                  STAT=errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) &
          CALL ErrorStop( global,ERR_ALLOCATE,__LINE__, &
                          'pPlag%requestsCECellsR' )

      ENDDO     ! iLev
    ENDIF       ! regions
  ENDDO         ! iReg

! ******************************************************************************
! Communicate buffer size 
! ****************************************************************************** 

! =============================================================================                           
! Send buffer size 
! ============================================================================= 

  DO iReg = 1, global%nRegions
    IF ( regions(iReg)%procid==global%myProcid .AND. &  ! region active and
         regions(iReg)%active==ACTIVE          .AND. &  ! on my processor
         global%plagUsed    .EQV. .TRUE.       .AND. &   ! plag active
         global%nProcAlloc > 1                  ) THEN  ! only if multiple processors        

#ifdef PLAG_CECELLS_MPI_DEBUG
  IF (iReg==1) THEN
    WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    WRITE(*,*) ' Entering PLAG_CECellsSendSize: pid, iReg = ',global%myProcId, iReg
    WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  ENDIF ! iReg
#endif

      CALL PLAG_CECellsSendSize(regions, iReg)
    ENDIF       ! regions    
  ENDDO         ! iReg

! =============================================================================  
! Receive buffer size 
! =============================================================================  

  DO iReg = 1, global%nRegions
    IF ( regions(iReg)%procid==global%myProcid .AND. &  ! region active and
         regions(iReg)%active==ACTIVE          .AND. &  ! on my processor
         global%plagUsed    .EQV. .TRUE.       .AND. &   ! plag active
         global%nProcAlloc > 1                  ) THEN  ! only if multiple processors        

#ifdef PLAG_CECELLS_MPI_DEBUG
  IF (iReg==1) THEN
    WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    WRITE(*,*) ' Entering PLAG_CECellsRecvSize: pid, iReg = ',global%myProcId, iReg
    WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  ENDIF ! iReg
#endif

      CALL PLAG_CECellsRecvSize(regions, iReg)
    ENDIF ! regions
  ENDDO ! iReg 

! =============================================================================
! Wait for messages of buffer size being received by other processors 
! =============================================================================
 
  DO iReg = 1, global%nRegions
    IF ( regions(iReg)%procid==global%myProcid .AND. &  ! region active and
         regions(iReg)%active==ACTIVE          .AND. &  ! on my processor
         global%plagUsed    .EQV. .TRUE.       .AND. &   ! plag active
         global%nProcAlloc > 1                  ) THEN  ! only if multiple processors

#ifdef PLAG_CECELLS_MPI_DEBUG
  IF (iReg==1) THEN
    WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    WRITE(*,*) ' Entering PLAG_CECellsClearRequestsSize: pid, iReg = ',global%myProcId, iReg
    WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  ENDIF ! iReg
#endif

      CALL PLAG_CECellsClearRequestsSize(regions, iReg)
    ENDIF ! regions  
  ENDDO ! iReg
  
! ******************************************************************************
! Allocate buffers for MPI-based communication 
! ******************************************************************************        

  DO iReg = 1, global%nRegions  
    IF (  regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE          .AND. &   ! on my processor
          global%plagUsed    .EQV. .TRUE.       .AND. &   ! plag active
          global%nProcAlloc > 1                  ) THEN   ! only if multiple processors

      DO iLev=1,regions(iReg)%nGridLevels
        pPlag => regions(iReg)%levels(iLev)%plag
        DO ir=1,global%nRegions
          pSendEcCell => regions(iReg)%levels(iLev)%sendEcCells(ir)
          pRecvEcCell => regions(iReg)%levels(iLev)%recvEcCells(ir)

! ----------------------------------------------------------------------------- 
!         Allocate send buffers 
! ----------------------------------------------------------------------------- 

          IF ( pSendEcCell%nCells > 0 ) THEN
            nBuffSendI = nDimBuffI*pSendEcCell%nBuffSizePlag
            nBuffSendR = nDimBuffR*pSendEcCell%nBuffSizePlag

            IF ( pSendEcCell%nBuffSizePlag > 0 ) THEN

#ifdef PLAG_CECELLS_MPI_DEBUG
  IF(iReg==1 .OR. iReg==10)THEN
  WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  WRITE(*,*) ' Allocating Send Data Buffers: pid, iReg, ir, nBuffSizePlag, nBuffSendI, nBuffSendR = ', &
        global%myProcId, iReg, ir, pSendEcCell%nBuffSizePlag, nBuffSendI, nBuffSendR
  ENDIF
#endif

              ALLOCATE( pSendEcCell%buffPlagI(nBuffSendI),STAT=errorFlag )
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) &
                CALL ErrorStop( global,ERR_ALLOCATE,__LINE__, &
                                'pSendEcCell%buffPlagI' )

              ALLOCATE( pSendEcCell%buffPlagR(nBuffSendR),STAT=errorFlag )
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) &
                CALL ErrorStop( global,ERR_ALLOCATE,__LINE__, &
                                'pSendEcCell%buffPlagR' )

! ----------------------------------------------------------------------------- 
!             Initialize send buffers 
! ----------------------------------------------------------------------------- 

              DO iBuffI = 1, nBuffSendI
                pSendEcCell%buffPlagI(iBuffI) = 0
              ENDDO ! iBuffI

              DO iBuffR= 1, nBuffSendR
                pSendEcCell%buffPlagR(iBuffR) = 0.0_RFREAL
              ENDDO ! iBuffR
             ENDIF !pSendEcCell%nBuffSizePlag

          ENDIF ! pSendEcCell%nCells

! ----------------------------------------------------------------------------- 
!         Allocate receive buffers 
! ----------------------------------------------------------------------------- 

          IF ( pRecvEcCell%nCells > 0 ) THEN
            nBuffRecvI = nDimBuffI*pRecvEcCell%nBuffSizePlag
            nBuffRecvR = nDimBuffR*pRecvEcCell%nBuffSizePlag

            IF ( pRecvEcCell%nBuffSizePlag > 0 ) THEN 

#ifdef PLAG_CECELLS_MPI_DEBUG
  IF (iReg==1 .OR. iReg==10) THEN
    WRITE(*,*) ' Allocating Receive Data Buffers: pid, iReg, ir, nBuffSizePlag, nBuffRecvI, nBuffRecvR = ', &
      global%myProcId, iReg, ir, pRecvEcCell%nBuffSizePlag, nBuffRecvI, nBuffRecvR
    WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  ENDIF ! iReg
#endif          

              ALLOCATE( pRecvEcCell%buffPlagI(nBuffRecvI),STAT=errorFlag )
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) &
              CALL ErrorStop( global,ERR_ALLOCATE,__LINE__, &
                              'pRecvEcCell%buffPlagI' )
 
              ALLOCATE( pRecvEcCell%buffPlagR(nBuffRecvR),STAT=errorFlag )
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) &
              CALL ErrorStop( global,ERR_ALLOCATE,__LINE__, &
                              'pRecvEcCell%buffPlagR' )

! ----------------------------------------------------------------------------- 
!             Initialize receive buffers 
! ----------------------------------------------------------------------------- 

              DO iBuffI = 1, nBuffRecvI
                pRecvEcCell%buffPlagI(iBuffI) = 0
              ENDDO ! iBuffI

              DO iBuffR= 1, nBuffRecvR
                pRecvEcCell%buffPlagR(iBuffR) = 0.0_RFREAL
              ENDDO ! iBuffR
             ENDIF ! pRecvEcCell%nBuffSizePlag 

          ENDIF ! pRecvEcCell%nCell 
        ENDDO   ! ir
      ENDDO     ! iLev
    ENDIF       ! regions
  ENDDO         ! iReg

! ******************************************************************************
! Synchronize through an MPI barrier 
! ******************************************************************************

#ifdef MPI
  CALL MPI_Barrier( global%mpiComm,global%mpierr )
  IF ( global%mpierr /= ERR_NONE ) &
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! ******************************************************************************
! Communicate buffer data 
! ******************************************************************************                           

! ============================================================================= 
! Send buffer data 
! ============================================================================= 

  DO iReg = 1, global%nRegions
    IF ( regions(iReg)%procid==global%myProcid .AND. &   ! region active and
         regions(iReg)%active==ACTIVE          .AND. &   ! on my processor
         global%plagUsed    .EQV. .TRUE.       .AND. &   ! plag active
         global%nProcAlloc > 1                  ) THEN   ! only if multiple processors        

#ifdef PLAG_CECELLS_MPI_DEBUG
  IF (iReg==1) THEN
    WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    WRITE(*,*) ' Entering PLAG_CECellsSendData: pid, iReg = ',global%myProcId, iReg
    WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  ENDIF ! iReg
#endif

      CALL PLAG_CECellsSendData(regions, iReg)    
    ENDIF       ! regions    
  ENDDO         ! iReg

! =============================================================================
! receive buffer size ---------------------------------------------------------
! =============================================================================

  DO iReg = 1, global%nRegions
    IF ( regions(iReg)%procid==global%myProcid .AND. &   ! region active and
         regions(iReg)%active==ACTIVE          .AND. &   ! on my processor
         global%plagUsed    .EQV. .TRUE.       .AND. &   ! plag active
         global%nProcAlloc > 1                  ) THEN   ! only if multiple processors    

#ifdef PLAG_CECELLS_MPI_DEBUG
  IF (iReg==1) THEN
    WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    WRITE(*,*) ' Entering PLAG_CECellsRecvData: pid, iReg = ',global%myProcId, iReg
    WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  ENDIF ! iReg
#endif

      CALL PLAG_CECellsRecvData(regions, iReg)

    ENDIF ! regions
  ENDDO ! iReg 

! =============================================================================
! Wait for message of buffer data being received by other processors
! =============================================================================

  DO iReg = 1, global%nRegions
    IF ( regions(iReg)%procid==global%myProcid .AND. &   ! region active and
         regions(iReg)%active==ACTIVE          .AND. &   ! on my processor
         global%plagUsed    .EQV. .TRUE.       .AND. &   ! plag active
         global%nProcAlloc > 1                  ) THEN   ! only if multiple processors

#ifdef PLAG_CECELLS_MPI_DEBUG
  IF (iReg==1) THEN
    WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    WRITE(*,*) ' Entering PLAG_CECellsClearRequestsData: pid, iReg = ',global%myProcId, iReg
    WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  ENDIF ! iReg
#endif

      CALL PLAG_CECellsClearRequestsData(regions, iReg)
    ENDIF ! regions  
  ENDDO ! iReg

! ******************************************************************************
! synchronize through an MPI barrier 
! ******************************************************************************

#ifdef MPI
  CALL MPI_Barrier( global%mpiComm,global%mpierr )
  IF ( global%mpierr /= ERR_NONE ) &
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! ******************************************************************************
! Deallocate communicaton buffers 
! ******************************************************************************

  DO iReg = 1, global%nRegions  
    IF ( regions(iReg)%procid==global%myProcid .AND. &   ! region active and
         regions(iReg)%active==ACTIVE          .AND. &   ! on my processor
         global%plagUsed    .EQV. .TRUE.       .AND. &   ! plag active
         global%nProcAlloc > 1                 ) THEN    ! only if multiple processors

      DO iLev=1,regions(iReg)%nGridLevels
        pPlag => regions(iReg)%levels(iLev)%plag

        DO ir=1,global%nRegions
          pSendEcCell => regions(iReg)%levels(iLev)%sendEcCells(ir)
          pRecvEcCell => regions(iReg)%levels(iLev)%recvEcCells(ir)

! ----------------------------------------------------------------------------- 
!         Deallocate send buffers 
! ----------------------------------------------------------------------------- 
          
          IF ( pSendEcCell%nCells > 0 ) THEN
            IF ( pSendEcCell%nBuffSizePlag > 0 ) THEN
              DEALLOCATE( pSendEcCell%buffPlagI, STAT=errorFlag )
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) &
                CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__, &
                                'pSendEcCell%buffCECellsPlagI'  )

              DEALLOCATE( pSendEcCell%buffPlagR, STAT=errorFlag )
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) &
                CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__, &
                                'pSendEcCell%buffCECellsPlagR'  )
             ENDIF ! pSendEcCell%nBuffSizePlag
          ENDIF ! pSendEcCell%nCells

! ----------------------------------------------------------------------------- 
!         Deallocate receive buffers 
! ----------------------------------------------------------------------------- 

          IF ( pRecvEcCell%nCells > 0 ) THEN
            IF ( pRecvEcCell%nBuffSizePlag > 0 ) THEN
              DEALLOCATE( pRecvEcCell%buffPlagI, STAT=errorFlag )
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) &
              CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__, &
                              'pRecvEcCell%buffPlagI'  )

              DEALLOCATE( pRecvEcCell%buffPlagR, STAT=errorFlag )
              global%error = errorFlag
              IF ( global%error /= ERR_NONE ) &
              CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__, &
                              'pRecvEcCell%buffPlagR'  )
             ENDIF ! pRecvEcCell%nBuffSizePlag
          ENDIF ! pRecvEcCell%nCell 
        ENDDO   ! ir

! -----------------------------------------------------------------------------
!       Deallocate requests 
! -----------------------------------------------------------------------------

        DEALLOCATE(pPlag%requestsCECells, STAT=errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) &
          CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

        DEALLOCATE(pPlag%requestsCECellsI, STAT=errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) &
          CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

        DEALLOCATE(pPlag%requestsCECellsR, STAT=errorFlag )
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) &
          CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

      ENDDO     ! iLev

    ENDIF       ! regions
  ENDDO         ! iReg

! ******************************************************************************
! Synchronize through an MPI barrier 
! ******************************************************************************

#ifdef MPI
  CALL MPI_Barrier( global%mpiComm,global%mpierr )
  IF ( global%mpierr /= ERR_NONE ) &
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! ******************************************************************************  
! Finalize
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CECellsSendRecvWrapper

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CECellsSendRecvWrapper.F90,v $
! Revision 1.5  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.2  2005/09/09 17:35:29  fnajjar
! Bug fix for initialization of recv buffers
!
! Revision 1.1  2004/12/01 20:57:19  fnajjar
! Initial revision after changing case
!
! Revision 1.3  2004/03/20 22:31:25  fnajjar
! Initialized buffer sizes after every RK stage
!
! Revision 1.2  2004/03/18 21:40:58  fnajjar
! Added calls to routines for sending-receiving MPI-based buffer data
!
! Revision 1.1  2004/03/10 23:15:33  fnajjar
! Initial import of main MPI-based wrapper
!
!******************************************************************************







