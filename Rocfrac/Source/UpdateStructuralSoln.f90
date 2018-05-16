!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************

MODULE UpdateStructuralSoln

  USE ROCSTAR_RocFrac
  USE ROCSTAR_RocFracComm

CONTAINS

  SUBROUTINE UpdateStructural(glb,NumProcs,Rnet)
    
    IMPLICIT NONE
    
    INCLUDE 'mpif.h'
    
    TYPE(ROCFRAC_GLOBAL) :: glb
    
    INTEGER :: NumNP ! Number of Node Points
    INTEGER :: NumElVol ! Number of Volumetric Elements
    INTEGER :: NumMatVol ! Number of Volumetreic Materials
    INTEGER :: NumMatCoh ! Number of Cohesive Materials
    INTEGER :: NumProcs ! Number of processors
    INTEGER :: iEltype ! Order of element (4:4node, 5:4nodeEnhanced, 10:10node)
    REAL*8,  DIMENSION(1:3*glb%NumNP) :: Rnet
    REAL*8,  POINTER, DIMENSION(:,:) :: MeshCoor
    INTEGER, POINTER, DIMENSION(:) :: MatIdVol
    INTEGER, POINTER, DIMENSION(:,:) :: ElConnVol
    REAL*8,  POINTER, DIMENSION(:) :: E, xnu, rho
    REAL*8,  POINTER, DIMENSION(:) :: Disp ! Nodal Displacement
    REAL*8,  POINTER, DIMENSION(:) :: deltan, deltat
    REAL*8,  POINTER, DIMENSION(:) :: VeloHalf ! velocity
    REAL*8,  POINTER, DIMENSION(:) :: VeloBar, AccelBar
    REAL*8 :: KappaDamp
    LOGICAL :: ALEenabled, DampEnabled
    INTEGER,  POINTER, DIMENSION(:) :: iSolnType
!    LOGICAL :: DefConfig
    INTEGER :: TotNumNeighProcs
    INTEGER :: MPI_COMM_ROCFRAC
!--   Non-block receive, Non-block send request arrays
    INTEGER,  POINTER, DIMENSION(:) :: ReqRcv, ReqSnd
    INTEGER,  POINTER, DIMENSION(:,:) :: StatRcv, StatSnd
    INTEGER :: TotNumNdComm
    
    INTEGER :: NumElPartBndry
    INTEGER, DIMENSION(:), POINTER  :: NumElPartBndryMat
    INTEGER, DIMENSION(:), POINTER  :: NeighProcList
    INTEGER, DIMENSION(:), POINTER  :: NumNdComm
    INTEGER, DIMENSION(:), POINTER  :: NumElVolMat
    
    TYPE(send_buf), DIMENSION(:), POINTER :: NdCommList
    
    TYPE(rcv_buf), DIMENSION(:), POINTER :: RecvDataFrm
    
    TYPE(send_buf), POINTER :: pNdComm
    TYPE(rcv_buf), POINTER :: pRecvDF

    REAL*8, ALLOCATABLE, DIMENSION(:) :: buf

    INTEGER :: j, j1, k , k1, k2
    INTEGER :: ierr,iAnalysisType
    INTEGER :: ElemStart, ElemEnd
    REAL*8,  POINTER, DIMENSION(:) :: SigmaMax, TauMax, Sinit
    REAL*8,  POINTER, DIMENSION(:,:) :: Sthresh1 ,Sthresh2

!------ CALCULATE THE SUBMESH'S BOUNDARY INTERNAL FORCE VECTOR


    NumNP      = glb%NumNP
    NumElVol   = glb%NumElVol
    NumMatVol  = glb%NumMatVol
    NumMatCoh  = glb%NumMatCoh
    iEltype    = glb%iEltype

!    Rnet => glb%Rnet
    MeshCoor   => glb%MeshCoor
    MatIdVol   => glb%MatIdVol
    ElConnVol  => glb%ElConnVol
    Disp       => glb%Disp
    deltan     => glb%deltan
    deltat     => glb%deltat
    SigmaMax   => glb%SigmaMax
    TauMax     => glb%TauMax
    Sthresh1    => glb%Sthresh1
    Sthresh2    => glb%Sthresh2
    Sinit      => glb%Sinit
    VeloHalf   => glb%VeloHalf
    VeloBar    => glb%VeloBar
    AccelBar   => glb%AccelBar
    KappaDamp   = glb%KappaDamp
    ALEenabled  = glb%ALEenabled
    DampEnabled = glb%DampEnabled
    iSolnType   => glb%iSolnType
!    DefConfig = glb%DefConfig
    TotNumNeighProcs = glb%TotNumNeighProcs
    MPI_COMM_ROCFRAC = glb%MPI_COMM_ROCFRAC
    ReqRcv  => glb%ReqRcv
    ReqSnd  => glb%ReqSnd 
    StatRcv => glb%StatRcv
    StatSnd => glb%StatSnd
    TotNumNdComm     = glb%TotNumNdComm
    NumElPartBndry   = glb%NumElPartBndry
    NumElPartBndryMat => glb%NumElPartBndryMat
    NeighProcList     => glb%NeighProcList
    NumNdComm         => glb%NumNdComm
    NumElVolMat       => glb%NumElVolMat
    NdCommList        => glb%NdCommList
    RecvDataFrm       => glb%RecvDataFrm
    E          => glb%E
    xnu        => glb%xnu
    rho        => glb%rho
    
    ElemStart = 1
    
    DO j = 1, NumMatVol
       
       ElemEnd = NumElPartBndryMat(j) + ElemStart - 1
       
       iAnalysisType = iSolnType(j)
       
       IF(iEltype.EQ.4)THEN
          CALL InternalForce_v3d4( glb, Rnet, ElemStart, ElemEnd, iAnalysisType)
       ELSE IF(iEltype.EQ.10)THEN
          CALL InternalForce_v3d10( glb, Rnet, ElemStart, ElemEnd, iAnalysisType)
       ELSE IF(iEltype.EQ.8)THEN
          CALL InternalForce_v3d8( glb, Rnet, ElemStart, ElemEnd, iAnalysisType)
       ENDIF
       
       ElemStart = NumElVolMat(j) + ElemStart
    ENDDO

    DO j = 1, NumMatCoh

       ! Transfer the cohesive solution to the 'Top' Surface

       CALL c3d6nm( glb%nsubn1, glb%nsubf1, &
            glb%nsubn2, glb%nsubf2, &
            glb%sd_coor1, glb%sd_subfaces1, &
            glb%sd_subface_parents1, &
            glb%sd_subface_parents2, &
            glb%sd_subface_counterparts1, &
            glb%sd_subface_nat_coors1, & 
            glb%sd_subface_nat_coors2, & 
            glb%FaceOfVolEl1, & 
            glb%FaceOfVolEl2, &
            glb%NumNp, glb%NumElVol, &
            glb%ElConnVol, &
            glb%nf1,glb%nf2, &
            glb%MapFaceEl2Vol1, &
            glb%MapFaceEl2Vol2, &
            deltan, deltat, SigmaMax, TauMax, Sinit, &
            Rnet, Disp, Sthresh1, glb%NumMatCoh,j,-1.d0)

       ! Transfer the cohesive solution to the 'Bottom' Surface

!!$       CALL c3d6nm( glb%nsubn2, glb%nsubf2, &
!!$            glb%sd_coor2, glb%sd_subfaces2, &
!!$            glb%sd_subface_parents2, &
!!$            glb%sd_subface_counterparts2, &
!!$            glb%sd_subface_nat_coors2, & 
!!$            glb%FaceOfVolEl2, &
!!$            glb%NumNp, glb%NumElVol, &
!!$            glb%ElConnVol, &
!!$            glb%nf2, &
!!$            glb%MapFaceEl2Vol2, &
!!$            deltan, deltat, SigmaMax, TauMax, Sinit, &
!!$            Rnet, Disp, Sthresh1, glb%NumMatCoh,j,1.d0)


    ENDDO

!     
!----- FORM THE BUFFER CONTAINING COMMUNICATED NODAL VALUES
!
    ALLOCATE(buf(1:TotNumNdComm))
    k1 = 1
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)
       pNdComm => NdCommList(j1)
       DO j = 1, NumNdComm(j1)
          k2 = 3*pNdComm%NdId(j)
          buf(k1)   = Rnet( k2 - 2 )
          buf(k1+1) = Rnet( k2 - 1 )
          buf(k1+2) = Rnet( k2 )
          k1 = k1 + 3
       ENDDO
    ENDDO

!     
!-MPI- RECEIVE THE INTERNAL FORCE VECTOR FROM THE NEIGHBORS
!
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)+1
       CALL MPI_IRECV(RecvDataFrm(k)%rcvbuf(1),NumNdComm(j1)*3, &
            MPI_DOUBLE_PRECISION,k-1,10,MPI_COMM_ROCFRAC, &
            ReqRcv(j1),ierr)
    ENDDO
!     
!-MPI- SEND THE INTERNAL FORCE VECTOR TO THE NEIGHBORS
!     
    k2 = 1
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)
       CALL MPI_ISEND(buf(k2), NumNdComm(j1)*3,MPI_DOUBLE_PRECISION, &
            k,10,MPI_COMM_ROCFRAC,ReqSnd(j1),ierr)
       k2 = k2 + NumNdComm(j1)*3
    ENDDO

!     
!-MPI- WAIT FOR INTERNAL FORCE VECTOR COMMUNICATION TO COMPLETE
!

    IF(ALEenabled)THEN
       IF(iEltype.EQ.4)THEN
          CALL v3d4_ale(VeloBar,AccelBar,Disp,VeloHalf,Rnet, &
               E,xnu,rho,NumNP,NumMatVol, &
               NumElVol,MatIdVol,ElConnVol,MeshCoor, &
               NumElPartBndry+1,NumElVol)
       ELSE
          CALL V3D10_ALE(VeloBar,AccelBar,Disp,VeloHalf,Rnet, &
               E,xnu,rho,NumNP,NumMatVol, &
               NumElVol,MatIdVol,ElConnVol,MeshCoor, &
               NumElPartBndry+1,NumElVol)
       ENDIF
    ENDIF

!-- (11) calculate R_in, R_damp

!------ CALCULATE THE SUBMESH'S BOUNDARY INTERNAL FORCE VECTOR

    ElemEnd = 0
    DO j = 1, NumMatVol
       
       ElemStart =  ElemEnd + NumElPartBndryMat(j) + 1
       
       ElemEnd = NumElVolMat(j) + ElemEnd
       
       IF(iEltype.EQ.4)THEN
          CALL InternalForce_v3d4( glb, Rnet, ElemStart, ElemEnd, iAnalysisType)
       ELSE IF(iEltype.EQ.10)THEN
          CALL InternalForce_v3d10( glb, rnet, ElemStart, ElemEnd, iAnalysisType)
       ELSE IF(iEltype.EQ.8)THEN
          CALL InternalForce_v3d8( glb, rnet, ElemStart, ElemEnd, iAnalysisType)
       ENDIF
    ENDDO

    IF(TotNumNeighProcs.GT.0)THEN
       CALL MPI_WAITALL(TotNumNeighProcs,ReqRcv,StatRcv,ierr)
       CALL MPI_WAITALL(TotNumNeighProcs,ReqSnd,StatSnd,ierr)
    ENDIF
    
    DEALLOCATE(buf)

!
!----- ADD NEIGHBOR'S CONTRIBUTION TO THE INTERNAL FORCE VECTOR
!
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)+1
       k1 = 1
       pNdComm => NdCommList(j1)
       pRecvDF => RecvDataFrm(k)
       DO j = 1, NumNdComm(j1)
          k2 = ( pNdComm%NdId(j) )*3
          Rnet(k2-2)= Rnet(k2-2) + pRecvDF%rcvbuf(k1)
          Rnet(k2-1)= Rnet(k2-1) + pRecvDF%rcvbuf(k1+1)
          Rnet(k2)  = Rnet(k2)   + pRecvDF%rcvbuf(k1+2)
          k1 = k1 + 3
       ENDDO
    ENDDO
    
    RETURN
  END SUBROUTINE UpdateStructural

  SUBROUTINE UpdateStructuralHT(glb,NumProcs,Rnet,RnetHT)
    
    IMPLICIT NONE
    
    INCLUDE 'mpif.h'
    
    TYPE(ROCFRAC_GLOBAL) :: glb
    
    INTEGER :: NumNP ! Number of Node Points
    INTEGER :: NumElVol ! Number of Volumetric Elements
    INTEGER :: NumMatVol ! Number of Volumetreic Materials
    INTEGER :: NumProcs ! Number of processors
    INTEGER :: iEltype ! Order of element (4:4node, 5:4nodeEnhanced, 10:10node)
    REAL*8,  DIMENSION(1:3*glb%NumNP) :: Rnet
    REAL*8,  DIMENSION(1:glb%NumNP) :: RnetHT
    REAL*8,  POINTER, DIMENSION(:,:) :: MeshCoor
    INTEGER, POINTER, DIMENSION(:) :: MatIdVol
    INTEGER, POINTER, DIMENSION(:,:) :: ElConnVol
    REAL*8,  POINTER, DIMENSION(:) :: E, xnu, rho, Cp, KappaHT
    REAL*8,  POINTER, DIMENSION(:) :: Disp ! Nodal Displacement
    REAL*8,  POINTER, DIMENSION(:) :: VeloHalf ! velocity
    REAL*8,  POINTER, DIMENSION(:) :: VeloBar, AccelBar
    REAL*8 :: KappaDamp
    LOGICAL :: ALEenabled, DampEnabled
    INTEGER,  POINTER, DIMENSION(:) :: iSolnType
!    LOGICAL :: DefConfig
    INTEGER :: TotNumNeighProcs
    INTEGER :: MPI_COMM_ROCFRAC
!--   Non-block receive, Non-block send request arrays
    INTEGER,  POINTER, DIMENSION(:) :: ReqRcv, ReqSnd
    INTEGER,  POINTER, DIMENSION(:,:) :: StatRcv, StatSnd
    INTEGER :: TotNumNdComm
    
    INTEGER :: NumElPartBndry
    INTEGER, DIMENSION(:), POINTER  :: NumElPartBndryMat
    INTEGER, DIMENSION(:), POINTER  :: NeighProcList
    INTEGER, DIMENSION(:), POINTER  :: NumNdComm
    INTEGER, DIMENSION(:), POINTER  :: NumElVolMat

    integer :: k4
    
    TYPE(send_buf), DIMENSION(:), POINTER :: NdCommList
    
    TYPE(rcv_buf), DIMENSION(:), POINTER :: RecvDataFrm
    
    TYPE(send_buf), POINTER :: pNdComm
    TYPE(rcv_buf), POINTER :: pRecvDF

    REAL*8, ALLOCATABLE, DIMENSION(:) :: buf

    INTEGER :: i, j, j1, k , k1, k2
    INTEGER :: ierr,iAnalysisType
    INTEGER :: ElemStart, ElemEnd
    logical :: HeatTransSoln

    integer :: NumNdsBCHT
    real*8, DIMENSION(:), POINTER :: Temperature
    integer, DIMENSION(:,:), POINTER :: BCFlagHT
    real*8, DIMENSION(:,:), POINTER  :: BCvalueHT

!------ CALCULATE THE SUBMESH'S BOUNDARY INTERNAL FORCE VECTOR


    NumNP      = glb%NumNP
    NumElVol   = glb%NumElVol
    NumMatVol  = glb%NumMatVol
    iEltype    = glb%iEltype

!    Rnet => glb%Rnet
    MeshCoor   => glb%MeshCoor
    MatIdVol   => glb%MatIdVol
    ElConnVol  => glb%ElConnVol
    Disp       => glb%Disp
    VeloHalf   => glb%VeloHalf
    VeloBar    => glb%VeloBar
    AccelBar   => glb%AccelBar
    KappaDamp   = glb%KappaDamp
    ALEenabled  = glb%ALEenabled
    DampEnabled = glb%DampEnabled
    iSolnType   => glb%iSolnType
!    DefConfig = glb%DefConfig
    TotNumNeighProcs = glb%TotNumNeighProcs
    MPI_COMM_ROCFRAC = glb%MPI_COMM_ROCFRAC
    ReqRcv  => glb%ReqRcv
    ReqSnd  => glb%ReqSnd 
    StatRcv => glb%StatRcv
    StatSnd => glb%StatSnd
    TotNumNdComm     = glb%TotNumNdComm
    NumElPartBndry   = glb%NumElPartBndry
    NumElPartBndryMat => glb%NumElPartBndryMat
    NeighProcList     => glb%NeighProcList
    NumNdComm         => glb%NumNdComm
    NumElVolMat       => glb%NumElVolMat
    NdCommList        => glb%NdCommList
    RecvDataFrm       => glb%RecvDataFrm
    E          => glb%E
    xnu        => glb%xnu
    rho        => glb%rho
    Cp        => glb%Cp
    KappaHT        => glb%KappaHT
    HeatTransSoln = glb%HeatTransSoln
    

    NumNdsBCHT = glb%NumNdsBCHT
    BCFlagHT => glb%BCFlagHT
    Temperature => glb%Temperature
    BCvalueHT => glb%BCvalueHT


    ElemStart = 1

! Solve for the Heat Transfer Solution first

    DO j = 1, NumMatVol
       
       ElemEnd = NumElPartBndryMat(j) + ElemStart - 1
       
       iAnalysisType = iSolnType(j)
       
       DO i = 1,glb%NumNdsBCHT
          k4 = glb%BCFlagHT(1,i) ! node
          
          IF (glb%BCFlagHT(2,i).EQ.0) THEN ! impose temperature
             glb%Temperature(k4) = glb%BCvalueHT(1,i) ! *glb%prop
             
          ENDIF
       ENDDO

       IF(iEltype.EQ.4)THEN

          CALL v3d4_thermal(NumElVol, NumNP, ElConnVol, MeshCoor, KappaHT, &
               RnetHT, Temperature, Rho, Cp, MatIdVol, NumMatVol,VeloBar,ElemStart, ElemEnd)

       ELSE
          
          Call v3d10_thermal(NumElVol, NumNP, ElConnVol, MeshCoor, KappaHT, &
               RnetHT, Temperature, Rho, Cp, MatIdVol, NumMatVol, VeloBar,ElemStart, ElemEnd)

       ENDIF
       
       ElemStart = NumElVolMat(j) + ElemStart
    ENDDO
!     
!----- FORM THE BUFFER CONTAINING COMMUNICATED NODAL VALUES
!
    ALLOCATE(buf(1:TotNumNdComm/3))
    k1 = 1
    DO j1 = 1, TotNumNeighProcs
       pNdComm => NdCommList(j1)
       DO j = 1, NumNdComm(j1)
          buf(k1)   = RnetHT(pNdComm%NdId(j))
          k1 = k1 + 1
       ENDDO
    ENDDO

!     
!-MPI- RECEIVE THE INTERNAL FORCE VECTOR FROM THE NEIGHBORS
!
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)+1
       CALL MPI_IRECV(RecvDataFrm(k)%rcvbuf(1),NumNdComm(j1), &
            MPI_DOUBLE_PRECISION,k-1,10,MPI_COMM_ROCFRAC, &
            ReqRcv(j1),ierr)
    ENDDO
!     
!-MPI- SEND THE INTERNAL FORCE VECTOR TO THE NEIGHBORS
!     
    k2 = 1
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)
       CALL MPI_ISEND(buf(k2), NumNdComm(j1),MPI_DOUBLE_PRECISION, &
            k,10,MPI_COMM_ROCFRAC,ReqSnd(j1),ierr)
       k2 = k2 + NumNdComm(j1)
    ENDDO


!------ CALCULATE THE SUBMESH'S BOUNDARY INTERNAL FORCE VECTOR

    ElemEnd = 0
    DO j = 1, NumMatVol
       
       ElemStart =  ElemEnd + NumElPartBndryMat(j) + 1
       
       ElemEnd = NumElVolMat(j) + ElemEnd

       IF(iEltype.EQ.4)THEN

          CALL v3d4_thermal(NumElVol, NumNP, ElConnVol, MeshCoor, KappaHT, &
               RnetHT, Temperature, Rho, Cp, MatIdVol, NumMatVol,VeloBar,ElemStart, ElemEnd)
       ELSE
          
          Call v3d10_thermal(NumElVol, NumNP, ElConnVol, MeshCoor, KappaHT, &
               RnetHT, Temperature, Rho, Cp, MatIdVol, NumMatVol, VeloBar,ElemStart, ElemEnd)

       ENDIF
    ENDDO
    IF(TotNumNeighProcs.GT.0)THEN
       CALL MPI_WAITALL(TotNumNeighProcs,ReqRcv,StatRcv,ierr)
       CALL MPI_WAITALL(TotNumNeighProcs,ReqSnd,StatSnd,ierr)
    ENDIF
    
    DEALLOCATE(buf)

!
!----- ADD NEIGHBOR'S CONTRIBUTION TO THE INTERNAL FORCE VECTOR
!
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)+1
       k1 = 1
       pNdComm => NdCommList(j1)
       pRecvDF => RecvDataFrm(k)
       DO j = 1, NumNdComm(j1)
          k2 = pNdComm%NdId(j)
          RnetHT(k2)  = RnetHT(k2) + pRecvDF%rcvbuf(k1)
          k1 = k1 + 1
       ENDDO
    ENDDO

    Do i = 1, glb%NumNP

       glb%Temperature(i) = glb%Temperature(i) + glb%DT*glb%CapctInv(i)*RnetHT(i)
    ENDDO

    DO i = 1,glb%NumNdsBCHT
       k4 = glb%BCFlagHT(1,i) ! node
       
       IF (glb%BCFlagHT(2,i).EQ.0) THEN ! impose temperature
          glb%Temperature(k4) = glb%BCvalueHT(1,i) !*glb%prop
       ENDIF
    ENDDO


! Solve the Structural solution with temperature expansion


    ElemStart = 1

    DO j = 1, NumMatVol
       
       ElemEnd = NumElPartBndryMat(j) + ElemStart - 1
       
       iAnalysisType = iSolnType(j)
       
       IF(iEltype.EQ.4)THEN

          glb%S11(:,:) = 0.d0
          glb%S22(:,:) = 0.d0
          glb%S33(:,:) = 0.d0

          CALL InternalForce_v3d4HT( glb, rnet, ElemStart, ElemEnd, iAnalysisType, Temperature)

!!$          CALL v3d4_thermalExp(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%ci, &
!!$               glb%S11,glb%S22,glb%S33, &
!!$               glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol, glb%alpha, Temperature, glb%Temperature0)

       ELSE

          glb%S11(:,:) = 0.d0
          glb%S22(:,:) = 0.d0
          glb%S33(:,:) = 0.d0
          

          CALL InternalForce_v3d10( glb, rnet, ElemStart, ElemEnd, iAnalysisType)

          CALL v3d10_thermalExp(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%ci, &
               glb%S11,glb%S22,glb%S33, &
               glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol, glb%alpha, Temperature, glb%Temperature0)

       ENDIF
       
       ElemStart = NumElVolMat(j) + ElemStart
    ENDDO


!     
!----- FORM THE BUFFER CONTAINING COMMUNICATED NODAL VALUES
!
    ALLOCATE(buf(1:TotNumNdComm))
    k1 = 1
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)
       pNdComm => NdCommList(j1)
       DO j = 1, NumNdComm(j1)
          k2 = 3*pNdComm%NdId(j)
          buf(k1)   = Rnet( k2 - 2 )
          buf(k1+1) = Rnet( k2 - 1 )
          buf(k1+2) = Rnet( k2 )
          k1 = k1 + 3
       ENDDO
    ENDDO


!     
!-MPI- RECEIVE THE INTERNAL FORCE VECTOR FROM THE NEIGHBORS
!
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)+1
       CALL MPI_IRECV(RecvDataFrm(k)%rcvbuf(1),NumNdComm(j1)*3, &
            MPI_DOUBLE_PRECISION,k-1,10,MPI_COMM_ROCFRAC, &
            ReqRcv(j1),ierr)
    ENDDO

!     
!-MPI- SEND THE INTERNAL FORCE VECTOR TO THE NEIGHBORS
!     
    k2 = 1
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)
       CALL MPI_ISEND(buf(k2), NumNdComm(j1)*3,MPI_DOUBLE_PRECISION, &
            k,10,MPI_COMM_ROCFRAC,ReqSnd(j1),ierr)
       k2 = k2 + NumNdComm(j1)*3
    ENDDO

!     
!-MPI- WAIT FOR INTERNAL FORCE VECTOR COMMUNICATION TO COMPLETE
!

    IF(ALEenabled)THEN
       IF(iEltype.EQ.4)THEN
          CALL v3d4_ale(VeloBar,AccelBar,Disp,VeloHalf,Rnet, &
               E,xnu,rho,NumNP,NumMatVol, &
               NumElVol,MatIdVol,ElConnVol,MeshCoor, &
               NumElPartBndry+1,NumElVol)
       ELSE
          CALL V3D10_ALE(VeloBar,AccelBar,Disp,VeloHalf,Rnet, &
               E,xnu,rho,NumNP,NumMatVol, &
               NumElVol,MatIdVol,ElConnVol,MeshCoor, &
               NumElPartBndry+1,NumElVol)
       ENDIF
    ENDIF

!-- (11) calculate R_in, R_damp

!------ CALCULATE THE SUBMESH'S BOUNDARY INTERNAL FORCE VECTOR

    ElemEnd = 0
    DO j = 1, NumMatVol
       
       ElemStart =  ElemEnd + NumElPartBndryMat(j) + 1
       
       ElemEnd = NumElVolMat(j) + ElemEnd


       IF(iEltype.EQ.4)THEN

          CALL InternalForce_v3d4HT( glb, Rnet, ElemStart, ElemEnd, iAnalysisType, Temperature)

!!$          CALL v3d4_thermalExp(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%ci, &
!!$               glb%S11,glb%S22,glb%S33, &
!!$               glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol, glb%alpha, Temperature, glb%Temperature0)

       ELSE
          CALL InternalForce_v3d10( glb, rnet, ElemStart, ElemEnd, iAnalysisType)
          
          CALL v3d10_thermalExp(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%ci, &
               glb%S11,glb%S22,glb%S33, &
               glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol, glb%alpha, Temperature, glb%Temperature0)
       ENDIF
    ENDDO

    

    IF(TotNumNeighProcs.GT.0)THEN
       CALL MPI_WAITALL(TotNumNeighProcs,ReqRcv,StatRcv,ierr)
       CALL MPI_WAITALL(TotNumNeighProcs,ReqSnd,StatSnd,ierr)
    ENDIF
    
    DEALLOCATE(buf)

!
!----- ADD NEIGHBOR'S CONTRIBUTION TO THE INTERNAL FORCE VECTOR
!
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)+1
       k1 = 1
       pNdComm => NdCommList(j1)
       pRecvDF => RecvDataFrm(k)
       DO j = 1, NumNdComm(j1)
          k2 = ( pNdComm%NdId(j) )*3
          Rnet(k2-2)= Rnet(k2-2) + pRecvDF%rcvbuf(k1)
          Rnet(k2-1)= Rnet(k2-1) + pRecvDF%rcvbuf(k1+1)
          Rnet(k2)  = Rnet(k2)   + pRecvDF%rcvbuf(k1+2)
          k1 = k1 + 3
       ENDDO
    ENDDO


    RETURN
  END SUBROUTINE UpdateStructuralHT

  SUBROUTINE InternalForce_v3d4( glb, Rnet, ElemStart, ElemEnd, iAnalysisType)
    
    IMPLICIT NONE
    
    TYPE(ROCFRAC_GLOBAL) :: glb
    REAL*8,  DIMENSION(1:3*glb%NumNP) :: Rnet
    
    INTEGER, INTENT(IN) :: ElemStart, ElemEnd, iAnalysisType
    
    IF (iAnalysisType.EQ.0) THEN
       CALL V3D4_NL_ARRUDA_BOYCE(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp,&
            glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
            glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol, &
            glb%xmu,glb%xkappa)
    ELSE IF(iAnalysisType.EQ.1)THEN
       CALL V3D4_NL(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp,glb%ci, &
            glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
            glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol,glb%xmu,glb%xlambda)
    ELSE IF(iAnalysisType.EQ.-1)THEN
       CALL V3D4_NeoHookeanInCompress(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp,&
            glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
            glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol, &
            glb%xmu,glb%xkappa)
    ELSE
!------- Linear elastic
       CALL V3D4(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp,glb%ci, &
            glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
            glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol)
    ENDIF

! Damping Rdamp

    IF(glb%DampEnabled)THEN
       CALL v3d4_damping(glb%VeloHalf,Rnet,glb%NumNP, &
            glb%NumElVol,glb%ElConnVol,glb%MeshCoor, &
            ElemStart,ElemEnd, glb%KappaDamp)
    END IF

       
  END SUBROUTINE INTERNALFORCE_V3D4

  SUBROUTINE InternalForce_v3d4HT( glb, Rnet, ElemStart, ElemEnd, iAnalysisType, Temperature)
    
    IMPLICIT NONE
    
    TYPE(ROCFRAC_GLOBAL) :: glb
    REAL*8,  DIMENSION(1:3*glb%NumNP) :: Rnet
    REAL*8,  DIMENSION(1:glb%NumNP) :: Temperature
    
    INTEGER, INTENT(IN) :: ElemStart, ElemEnd, iAnalysisType
    
    IF (iAnalysisType.EQ.0) THEN
       CALL V3D4_NL_ARRUDA_BOYCE(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp,&
            glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
            glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol, &
            glb%xmu,glb%xkappa)
    ELSE IF(iAnalysisType.EQ.1)THEN
       CALL V3D4_NL(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp,glb%ci, &
            glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
            glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol,glb%xmu,glb%xlambda)
    ELSE IF(iAnalysisType.EQ.-1)THEN
       CALL V3D4_NeoHookeanInCompress(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp,&
            glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
            glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol, &
            glb%xmu,glb%xkappa)
    ELSE
!------- Linear elastic
       CALL V3D4_thermalExp2(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp,glb%ci, &
            glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
            glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol, &
            glb%alpha,Temperature,glb%Temperature0)
    ENDIF

! Damping Rdamp

    IF(glb%DampEnabled)THEN
       CALL v3d4_damping(glb%VeloHalf,Rnet,glb%NumNP, &
            glb%NumElVol,glb%ElConnVol,glb%MeshCoor, &
            ElemStart,ElemEnd, glb%KappaDamp)
    END IF

       
  END SUBROUTINE INTERNALFORCE_V3D4HT


  SUBROUTINE InternalForce_v3d10( glb, Rnet, ElemStart, ElemEnd, iAnalysisType)
    
    IMPLICIT NONE

    TYPE(ROCFRAC_GLOBAL) :: glb
    REAL*8,  DIMENSION(1:3*glb%NumNP) :: Rnet
    
    INTEGER, INTENT(IN) :: ElemStart, ElemEnd, iAnalysisType 

    IF(glb%DebondPart)THEN
       
       CALL V3D10_NL_HUANG(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp, &
            glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
            glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol, &
            glb%STATEV_Part1,glb%STATEV_Part2,glb%NSTATEV,glb%MATRIX,glb%NMATRIX, &
            glb%PARTICLE,glb%NPARTICLE,glb%NPARTICLETYPE,glb%INTERFAC,glb%NINTERFAC,glb%StrainTrace)

    ELSE IF(glb%DebondPart_Matous)THEN


       IF(glb%Debug_State) PRINT*,'Starting v3d10_nl_matous', ElemStart,ElemEnd

       CALL V3D10_NL_Matous(glb%MeshCoor,glb%ElConnVol,Rnet,glb%Disp, &
            glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
            glb%NumNP,ElemStart,ElemEnd,glb%NumElVol, &
            glb%p1, glb%p2, glb%Yin, glb%SoftParam, glb%BulkMod(2), glb%BulkMod(1), &
            glb%ShrMod(2),  glb%ShrMod(1),  glb%PoisRat(2),  glb%PoisRat(1), &
            glb%a_eta, glb%a_zeta,&
            glb%cm, glb%c2, glb%cd, glb%Lo, glb%L_tensor(:,:,1), &
            glb%L_tensor(:,:,2), glb%M_tensor(:,:,1), glb%M_tensor(:,:,2), glb%L_bar, &
            glb%StrainOld)

       IF(glb%Debug_State) PRINT*,'Finished v3d10_nl_matous', ElemStart,ElemEnd


    ELSE


       IF (iAnalysisType.EQ.0) THEN
          
          IF ( glb%ArtificialDamping)THEN
             
             CALL V3D10_NL_ARRUDA_BOYCE_DAMPING(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp, &
                  glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
                  glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol, &
                  glb%xmu,glb%xkappa, glb%rho, glb%cd_fastest, glb%DetF_old, glb%VeloHalf, glb%Dt)
             
          ELSE
             
             CALL V3D10_NL_ARRUDA_BOYCE(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp, &
                  glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
               glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol, &
               glb%xmu,glb%xkappa)
          ENDIF
          
       ELSE IF (iAnalysisType.EQ.1)THEN
       
       IF(glb%iElIntgratn.EQ.0)THEN

       IF ( glb%ArtificialDamping)THEN


          CALL V3D10_NL_DAMPING(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp,glb%ci, &
               glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
               glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol, &
               glb%rho, glb%cd_fastest, glb%DetF_old, glb%VeloHalf, glb%Dt)

       ELSE
          
          CALL V3D10_NL(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp,glb%ci, &
               glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
               glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol)
       ENDIF
       ELSE 
          
          CALL V3D10R_NL(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp,glb%ci,glb%cj, &
               glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
               glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol)
       ENDIF

    ELSE IF (iAnalysisType.EQ.2)THEN


       IF(glb%iElIntgratn.EQ.0)THEN
!------- Linear elastic

          CALL V3D10(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp,glb%ci, &
               glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
               glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol)

       ELSE

          CALL V3D10_B_BAR(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp,glb%ci, &
               glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
               glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol)

       ENDIF
    ENDIF
 endif
    
  END SUBROUTINE InternalForce_v3d10

  SUBROUTINE InternalForce_v3d8( glb, Rnet, ElemStart, ElemEnd, iAnalysisType)
    
    IMPLICIT NONE
    
    TYPE(ROCFRAC_GLOBAL) :: glb
    REAL*8,  DIMENSION(1:3*glb%NumNP) :: Rnet
    
    INTEGER, INTENT(IN) :: ElemStart, ElemEnd, iAnalysisType
    
!------- Linear elastic
    CALL v3d8_me(glb%MeshCoor,glb%MatIdVol,glb%ElConnVol,Rnet,glb%Disp,glb%dmat, &
         glb%S11,glb%S22,glb%S33,glb%S12,glb%S23,glb%S13, &
         glb%NumNP,ElemStart,ElemEnd,glb%NumElVol,glb%NumMatVol,glb%Aenh,glb%enhanced_map,glb%mixed_map)
       
  END SUBROUTINE INTERNALFORCE_V3D8
  
END MODULE UpdateStructuralSoln

