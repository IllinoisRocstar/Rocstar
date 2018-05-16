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
  SUBROUTINE v3d4n(coor,nodes, NumElNeigh,ElConn,Rnet,alpha,disp,ci,&
       numnp,numel,Ahat,numat_vol,&
       nprocs,TotNumNdComm,TotNumNeighProcs,NeighProcList,NumNdComm,neigh_lst)

    USE ROCSTAR_RocFracComm 

    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER :: k1, k2, j1, k

!    INTEGER :: MPI_COMM_WORLD  ! Communicator of Rocfrac
!    INTEGER :: MPI_DOUBLE_PRECISION
    INTEGER :: TotNumNeighProcs, Nprocs,TotNumNdComm
    INTEGER, DIMENSION(1:TotNumNeighProcs) :: NeighProcList
    INTEGER, DIMENSION(1:TotNumNeighProcs) ::NumNdComm
    REAL*8, allocatable, dimension(:) :: buf

    TYPE(rcv_buf), ALLOCATABLE, DIMENSION(:) :: RecvDataFrm
    TYPE(send_buf), DIMENSION(1:TotNumNeighProcs) :: neigh_lst

    INTEGER :: ierr
!--   Non-block receive, Non-block send request arrays
    INTEGER, ALLOCATABLE, DIMENSION(:) :: req_rcv_lst, req_snd_lst
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: stat_rcv_lst,stat_snd_lst

!!$
    INTEGER :: numnp, numel,i
    REAL*8, DIMENSION(1:numnp) :: Ahat
    REAL*8, DIMENSION(1:numnp*3) :: disp, Rnet
    REAL*8, DIMENSION(1:3,1:numnp) :: coor
    INTEGER, DIMENSION(1:numnp) :: NumElNeigh
    INTEGER, DIMENSION(1:numnp,1:40) :: ElConn  ! fix 40 should not be hard coded
    INTEGER,DIMENSION(1:4,1:numel) :: nodes
    REAL*8, DIMENSION(1:4,1:numel) :: alpha
    INTEGER :: numat_vol      ! number of volumetric materials
!--   elastic stiffness consts
    REAL*8, DIMENSION(1:9,1:numat_vol) :: ci

    INTEGER :: j
    REAL*8 :: aaa
!--   coordinate holding variable
    REAL*8 :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
    INTEGER :: k1n1,k1n2,k1n3,k1n4,k2n1,k2n2,k2n3,k2n4
    INTEGER :: k3n1,k3n2,k3n3,k3n4

    REAL*8 :: AhatInv
    INTEGER :: iElNum
    REAL*8 :: SixInv
!--   x, y, and z displacements of nodes
    REAL*8 :: u1,u2,u3,u4,v1,v2,v3,v4,w1,w2,w3,w4

!--  Coordinate subtractions
    REAL*8 :: x14, x24, x34, y14, y24, y34, z14, z24, z34
    REAL*8 :: c11, c21, c31
!   --   6*volume and the volume      
    REAL*8 :: Vx6
!--   spacial derivatives
    REAL*8 :: B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12
!--   strains
    REAL*8 :: E11,E22,E33,E12,E23,E13
    REAL*8 :: B13, B14, B15
    REAL*8 :: SumS11n,SumS22n,SumS33n,SumS12n,SumS23n,SumS13n
    INTEGER :: k3i,k2i,k1i
    INTEGER :: ix
    INTEGER :: n1, n2, n3, n4
    REAL*8, DIMENSION(1:numnp) :: S11n, S22n, S33n, S12n, S23n, S13n

!     CALL MPI_INIT(ierr)
!     CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
!     CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

    INTEGER :: myid 
    IF(TotNumNeighProcs.NE.0)THEN
       ALLOCATE(req_rcv_lst  (1:TotNumNeighProcs)  )
       ALLOCATE(req_snd_lst  (1:TotNumNeighProcs)  )

       ALLOCATE(stat_snd_lst  (1:MPI_STATUS_SIZE,1:TotNumNeighProcs)  )
       ALLOCATE(stat_rcv_lst  (1:MPI_STATUS_SIZE,1:TotNumNeighProcs)  )
    ENDIF

   SixInv = 1.d0/6.d0

    DO i = 1, numnp
       AhatInv = 1.d0/Ahat(i)

       SumS11n = 0.d0
       SumS22n = 0.d0
       SumS33n = 0.d0
       SumS12n = 0.d0
       SumS23n = 0.d0
       SumS13n = 0.d0
       
       DO j = 1, NumElNeigh(i)
            
          iElNum = ElConn(i,j)
            
          n1 = nodes(1,iElNum)
          n2 = nodes(2,iElNum)
          n3 = nodes(3,iElNum)
          n4 = nodes(4,iElNum)

          IF(i.EQ.n1) THEN
             ix = 1
          ELSE IF(i.EQ.n2)THEN
             ix = 2
          ELSE IF(i.EQ.n3)THEN
             ix = 3
          ELSE IF(i.EQ.n4)THEN
             ix = 4
          ENDIF
            
          k3n1 = 3*n1
          k3n2 = 3*n2
          k3n3 = 3*n3
          k3n4 = 3*n4
          
          k2n1 = k3n1 - 1
          k2n2 = k3n2 - 1
          k2n3 = k3n3 - 1
          k2n4 = k3n4 - 1
          
          k1n1 = k3n1 - 2
          k1n2 = k3n2 - 2
          k1n3 = k3n3 - 2
          k1n4 = k3n4 - 2         
          ! k#n# dummy variables replaces:
          u1 = disp(k1n1)        ! (3*n1-2)
          u2 = disp(k1n2)        ! (3*n2-2)
          u3 = disp(k1n3)        ! (3*n3-2)
          u4 = disp(k1n4)        ! (3*n4-2)
          v1 = disp(k2n1)        ! (3*n1-1)
          v2 = disp(k2n2)        ! (3*n2-1)
          v3 = disp(k2n3)        ! (3*n3-1)
          v4 = disp(k2n4)        ! (3*n4-1)
          w1 = disp(k3n1)        ! (3*n1)
          w2 = disp(k3n2)        ! (3*n2)
          w3 = disp(k3n3)        ! (3*n3)
          w4 = disp(k3n4)        ! (3*n4)
           
          x1 = coor(1,n1)
          x2 = coor(1,n2)
          x3 = coor(1,n3)
          x4 = coor(1,n4)
          y1 = coor(2,n1)
          y2 = coor(2,n2)
          y3 = coor(2,n3)
          y4 = coor(2,n4)
          z1 = coor(3,n1)
          z2 = coor(3,n2)
          z3 = coor(3,n3)
          z4 = coor(3,n4)

! See the maple worksheet 'V3D4.mws' for the derivation of [B] 

          B1  = -(y3*z4 - y4*z3 - y2*z4 + y2*z3 + z2*y4 - z2*y3) * SixInv
          B2  =  (x3*z4 - x4*z3 - x2*z4 + x2*z3 + z2*x4 - z2*x3) * SixInv
          B3  = -(x3*y4 - x4*y3 - x2*y4 + x2*y3 + y2*x4 - y2*x3) * SixInv
          B4  =  (y3*z4 - y4*z3 - y1*z4 + y1*z3 + z1*y4 - z1*y3) * SixInv
          B5  = -(x3*z4 - x4*z3 - x1*z4 + x1*z3 + z1*x4 - z1*x3) * SixInv
          B6  =  (x3*y4 - x4*y3 - x1*y4 + x1*y3 + y1*x4 - y1*x3) * SixInv
          B7  = -(y2*z4 - z2*y4 - y1*z4 + y1*z2 + z1*y4 - z1*y2) * SixInv
          B8  =  (x2*z4 - z2*x4 - x1*z4 + x1*z2 + z1*x4 - z1*x2) * SixInv
          B9  = -(x2*y4 - y2*x4 - x1*y4 + x1*y2 + y1*x4 - y1*x2) * SixInv
          B10 =  (y2*z3 - z2*y3 - y1*z3 + y1*z2 + z1*y3 - z1*y2) * SixInv
          B11 = -(x2*z3 - z2*x3 - x1*z3 + x1*z2 + z1*x3 - z1*x2) * SixInv
          B12 =  (x2*y3 - y2*x3 - x1*y3 + x1*y2 + y1*x3 - y1*x2) * SixInv

! calculate the strain            
!                        [E] = [B]{d}
!
          E11 = alpha(ix,iElNum)*(B1*u1 + B4*u2 + B7*u3 + B10*u4)
          E22 = alpha(ix,iElNum)*(B2*v1 + B5*v2 + B8*v3 + B11*v4)
          E33 = alpha(ix,iElNum)*(B3*w1 + B6*w2 + B9*w3 + B12*w4)
          E12 = alpha(ix,iElNum)*(B2*u1 + B1*v1 + B5*u2 + B4*v2 + B8*u3 + B7*v3 + B11*u4 + B10*v4)
          E23 = alpha(ix,iElNum)*(B3*v1 + B2*w1 + B6*v2 + B5*w2 + B9*v3 + B8*w3 + B12*v4 + B11*w4) 
          E13 = alpha(ix,iElNum)*(B3*u1 + B1*w1 + B6*u2 + B4*w2 + B9*u3 + B7*w3 + B12*u4 + B10*w4)
! calculate the stress            -1
!                        [S] = [C]  {E}
!       
          SumS11n = SumS11n + E11*ci(1,1) + E22*ci(2,1) + E33*ci(4,1)
          SumS22n = SumS22n + E11*ci(2,1) + E22*ci(3,1) + E33*ci(5,1)
          SumS33n = SumS33n + E11*ci(4,1) + E22*ci(5,1) + E33*ci(6,1)
          SumS12n = SumS12n + E12*ci(7,1)
          SumS23n = SumS23n + E23*ci(8,1)
          SumS13n = SumS13n + E13*ci(9,1)

       ENDDO
       S11n(i) = AhatInv*SumS11n
       S22n(i) = AhatInv*SumS22n
       S33n(i) = AhatInv*SumS33n
       S12n(i) = AhatInv*SumS12n
       S23n(i) = AhatInv*SumS23n
       S13n(i) = AhatInv*SumS13n ! stress at node

!       print*,myid,i,AhatInv*SumS23n

    ENDDO

    ALLOCATE(RecvDataFrm(0:nprocs-1))

!
!----- FORM THE BUFFER CONTAINING COMMUNICATED STRESS MATRIX NODAL VALUES
    
    ALLOCATE(buf(1:TotNumNdComm*2))
    k1 = 1
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)
       ALLOCATE(RecvDataFrm(k)%rcvbuf(1:NumNdComm(j1)*6))
       DO j = 1, NumNdComm(j1)
          k2 = neigh_lst(j1)%NdID(j) !NdCommList(k)%NdId(j)
!          print*,myid,k2,S23n(k2)
          buf(k1  ) = S11n(k2)
          buf(k1+1) = S22n(k2)
          buf(k1+2) = S33n(k2)
          buf(k1+3) = S12n(k2)
          buf(k1+4) = S23n(k2)
          buf(k1+5) = S13n(k2)
          k1 = k1 + 6
       ENDDO
    ENDDO
!
!-MPI- RECEIVE THE RECIPRICAL MASS MATRIX DIAGONAL FROM THE NEIGHBORS
!
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)
       CALL MPI_IRECV(RecvDataFrm(k)%rcvbuf(1),NumNdComm(j1)*6, &
            MPI_DOUBLE_PRECISION, k, 10, rocstar_communicator,req_rcv_lst(j1),ierr)
    ENDDO
!
!-MPI- SEND THE RECIPRICAL MASS MATRIX DIAGONAL TO THE NEIGHBORS
!
    k2 = 1
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)
       
!       IF(myid.EQ.0) print*,'sending',buf(1:) 
       CALL MPI_ISEND(buf(k2),NumNdComm(j1)*6,&
            MPI_DOUBLE_PRECISION,k,10,rocstar_communicator,req_snd_lst(j1),ierr)
       k2 = k2 + NumNdComm(j1)
    ENDDO
!
!-MPI- WAIT FOR INTERNAL FORCE VECTOR COMMUNICATION TO COMPLETE
!
    IF(TotNumNeighProcs.GT.0)THEN
       CALL MPI_WAITALL(TotNumNeighProcs,req_rcv_lst,stat_rcv_lst,ierr)
       CALL MPI_WAITALL(TotNumNeighProcs,req_snd_lst,stat_snd_lst,ierr)
    ENDIF
    DEALLOCATE(buf)

!
!----- ADD NEIGHBOR'S CONTRIBUTION TO THE INTERNAL FORCE VECTOR
!
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)
       k1 = 1
!       IF(myid.EQ.1) print*,'***',RecvDataFrm(0)%rcvbuf(:)
       DO j = 1, NumNdComm(j1)
          k2 = neigh_lst(j1)%NdID(j)
          S11n(k2) = S11n(k2) + RecvDataFrm(k)%rcvbuf(k1)
          S22n(k2) = S22n(k2) + RecvDataFrm(k)%rcvbuf(k1+1)
          S33n(k2) = S33n(k2) + RecvDataFrm(k)%rcvbuf(k1+2)
          S12n(k2) = S12n(k2) + RecvDataFrm(k)%rcvbuf(k1+3)
          S23n(k2) = S23n(k2) + RecvDataFrm(k)%rcvbuf(k1+4)
          S13n(k2) = S13n(k2) + RecvDataFrm(k)%rcvbuf(k1+5)
          k1 = k1 + 6
       ENDDO
    ENDDO

    DEALLOCATE(RecvDataFrm)

!    IF(myid.eq.0) print*,'1procs',myid,S11n(11),S22n(11),S33n(11),S12n(11),S23n(11),S13n(11)
!    IF(myid.EQ.1) print*,myid,S11n(11),S22n(11),S33n(11),S12n(11),S23n(11),S13n(11)
!    IF(myid.EQ.0) print*,myid,S11n(3),S22n(3),S33n(3),S12n(3),S23n(3),S13n(3)

    DO i = 1, numnp

       k3i = 3*i
       k2i = 3*i-1
       k1i = 3*i-2

       DO j = 1, NumElNeigh(i)
           
            iElNum = ElConn(i,j)
            
            n1 = nodes(1,iElNum)
            n2 = nodes(2,iElNum)
            n3 = nodes(3,iElNum)
            n4 = nodes(4,iElNum)

            k3n1 = 3*n1
            k3n2 = 3*n2
            k3n3 = 3*n3
            k3n4 = 3*n4
            
            k2n1 = k3n1 - 1
            k2n2 = k3n2 - 1
            k2n3 = k3n3 - 1
            k2n4 = k3n4 - 1

            k1n1 = k3n1 - 2
            k1n2 = k3n2 - 2
            k1n3 = k3n3 - 2
            k1n4 = k3n4 - 2 

            x1 = coor(1,n1)
            x2 = coor(1,n2)
            x3 = coor(1,n3)
            x4 = coor(1,n4)
            y1 = coor(2,n1)
            y2 = coor(2,n2)
            y3 = coor(2,n3)
            y4 = coor(2,n4)
            z1 = coor(3,n1)
            z2 = coor(3,n2)
            z3 = coor(3,n3)
            z4 = coor(3,n4)

! See the maple worksheet 'V3D4.mws' for the derivation of [B] 

!!$            B1  = -(y3*z4 - y4*z3 - y2*z4 + y2*z3 + z2*y4 - z2*y3) * SixInv
!!$            B2  =  (x3*z4 - x4*z3 - x2*z4 + x2*z3 + z2*x4 - z2*x3) * SixInv
!!$            B3  = -(x3*y4 - x4*y3 - x2*y4 + x2*y3 + y2*x4 - y2*x3) * SixInv
!!$            B4  =  (y3*z4 - y4*z3 - y1*z4 + y1*z3 + z1*y4 - z1*y3) * SixInv
!!$            B5  = -(x3*z4 - x4*z3 - x1*z4 + x1*z3 + z1*x4 - z1*x3) * SixInv
!!$            B6  =  (x3*y4 - x4*y3 - x1*y4 + x1*y3 + y1*x4 - y1*x3) * SixInv
!!$            B7  = -(y2*z4 - z2*y4 - y1*z4 + y1*z2 + z1*y4 - z1*y2) * SixInv
!!$            B8  =  (x2*z4 - z2*x4 - x1*z4 + x1*z2 + z1*x4 - z1*x2) * SixInv
!!$            B9  = -(x2*y4 - y2*x4 - x1*y4 + x1*y2 + y1*x4 - y1*x2) * SixInv
!!$            B10 =  (y2*z3 - z2*y3 - y1*z3 + y1*z2 + z1*y3 - z1*y2) * SixInv
!!$            B11 = -(x2*z3 - z2*x3 - x1*z3 + x1*z2 + z1*x3 - z1*x2) * SixInv
!!$            B12 =  (x2*y3 - y2*x3 - x1*y3 + x1*y2 + y1*x3 - y1*x2) * SixInv

            B1  = ( (y3-y4)*(z2-z4) - (y2-y4)*(z3-z4) ) * SixInv
            B2  = ( (z3-z4)*(x2-x4) - (z2-z4)*(x3-x4) ) * SixInv
            B3  = ( (x3-x4)*(y2-y4) - (x2-x4)*(y3-y4) ) * SixInv
            B4  = ( (y1-y3)*(z1-z4) - (y1-y4)*(z1-z3) ) * SixInv
            B5  = ( (z1-z3)*(x1-x4) - (z1-z4)*(x1-x3) ) * SixInv
            B6  = ( (x1-x3)*(y1-y4) - (x1-x4)*(y1-y3) ) * SixInv
            B7  = ( (y1-y4)*(z1-z2) - (y1-y2)*(z1-z4) ) * SixInv
            B8  = ( (z1-z4)*(x1-x2) - (z1-z2)*(x1-x4) ) * SixInv
            B9  = ( (x1-x4)*(y1-y2) - (x1-x2)*(y1-y4) ) * SixInv
            B10 = ( (y1-y2)*(z1-z3) - (y1-y3)*(z1-z2) ) * SixInv
            B11 = ( (z1-z2)*(x1-x3) - (z1-z3)*(x1-x2) ) * SixInv
            B12 = ( (x1-x2)*(y1-y3) - (x1-x3)*(y1-y2) ) * SixInv

            IF(n1.EQ.i)THEN
               B13 = B1
               B14 = B2
               B15 = B3
            ELSE IF (n2.EQ.i)THEN
               B13 = B4
               B14 = B5
               B15 = B6
            ELSE IF(n3.EQ.i)THEN
               B13 = B7
               B14 = B8
               b15 = B9
            ELSE IF(n4.EQ.i)THEN
               B13 = B10
               B14 = B11
               B15 = B12
            ENDIF

! ASSEMBLE THE INTERNAL FORCE VECTOR

            Rnet(k1i) = Rnet(k1i) - &
                 (S11n(n1)*B13 + S12n(n1)*B14 + S13n(n1)*B15)*alpha(1,iElNum) - & ! local node 1
                 (S11n(n2)*B13 + S12n(n2)*B14 + S13n(n2)*B15)*alpha(2,iElNum) - & ! local node 2
                 (S11n(n3)*B13 + S12n(n3)*B14 + S13n(n3)*B15)*alpha(3,iElNum) - & ! local node 3
                 (S11n(n4)*B13 + S12n(n4)*B14 + S13n(n4)*B15)*alpha(4,iElNum)     ! local node 4

            Rnet(k2i) = Rnet(k2i) - &
                 (S22n(n1)*B14 + S12n(n1)*B13 + S23n(n1)*B15)*alpha(1,iElNum) - &
                 (S22n(n2)*B14 + S12n(n2)*B13 + S23n(n2)*B15)*alpha(2,iElNum) - &
                 (S22n(n3)*B14 + S12n(n3)*B13 + S23n(n3)*B15)*alpha(3,iElNum) - &
                 (S22n(n4)*B14 + S12n(n4)*B13 + S23n(n4)*B15)*alpha(4,iElNum)

            Rnet(k3i) = Rnet(k3i) - &
                 (S33n(n1)*B15 + S23n(n1)*B14 + S13n(n1)*B13)*alpha(1,iElNum) - &
                 (S33n(n2)*B15 + S23n(n2)*B14 + S13n(n2)*B13)*alpha(2,iElNum) - &
                 (S33n(n3)*B15 + S23n(n3)*B14 + S13n(n3)*B13)*alpha(3,iElNum) - &
                 (S33n(n4)*B15 + S23n(n4)*B14 + S13n(n4)*B13)*alpha(4,iElNum)

!!$            IF(i.EQ.11.and.Myid.EQ.0)print*,'1proc',(S22n(n1)*B14 + S12n(n1)*B13 + S23n(n1)*B15)*alpha(1,iElNum) - &
!!$                 (S22n(n2)*B14 + S12n(n2)*B13 + S23n(n2)*B15)*alpha(2,iElNum) - &
!!$                 (S22n(n3)*B14 + S12n(n3)*B13 + S23n(n3)*B15)*alpha(3,iElNum) - &
!!$                 (S22n(n4)*B14 + S12n(n4)*B13 + S23n(n4)*B15)*alpha(4,iElNum)
!!$
!!$            IF(i.EQ.3.AND.myid.EQ.0)print*,(S22n(n1)*B14 + S12n(n1)*B13 + S23n(n1)*B15)*alpha(1,iElNum) - &
!!$                 (S22n(n2)*B14 + S12n(n2)*B13 + S23n(n2)*B15)*alpha(2,iElNum) - &
!!$                 (S22n(n3)*B14 + S12n(n3)*B13 + S23n(n3)*B15)*alpha(3,iElNum) - &
!!$                 (S22n(n4)*B14 + S12n(n4)*B13 + S23n(n4)*B15)*alpha(4,iElNum)
!!$
!!$            IF(i.EQ.11.AND.myid.EQ.1)print*,(S22n(n1)*B14 + S12n(n1)*B13 + S23n(n1)*B15)*alpha(1,iElNum) - &
!!$                 (S22n(n2)*B14 + S12n(n2)*B13 + S23n(n2)*B15)*alpha(2,iElNum) - &
!!$                 (S22n(n3)*B14 + S12n(n3)*B13 + S23n(n3)*B15)*alpha(3,iElNum) - &
!!$                 (S22n(n4)*B14 + S12n(n4)*B13 + S23n(n4)*B15)*alpha(4,iElNum)
!!$

!!$            IF(myid.EQ.0.and.k3i.eq.9) print*,- &
!!$                 (S33n(n1)*B15 + S23n(n1)*B14 + S13n(n1)*B13)*alpha(1,iElNum) - &
!!$                 (S33n(n2)*B15 + S23n(n2)*B14 + S13n(n2)*B13)*alpha(2,iElNum) - &
!!$                 (S33n(n3)*B15 + S23n(n3)*B14 + S13n(n3)*B13)*alpha(3,iElNum) - &
!!$                 (S33n(n4)*B15 + S23n(n4)*B14 + S13n(n4)*B13)*alpha(4,iElNum)
!!$            IF(myid.EQ.0.and.k3i.eq.33) print*,'1procs',- &
!!$                 (S33n(n1)*B15 + S23n(n1)*B14 + S13n(n1)*B13)*alpha(1,iElNum) - &
!!$                 (S33n(n2)*B15 + S23n(n2)*B14 + S13n(n2)*B13)*alpha(2,iElNum) - &
!!$                 (S33n(n3)*B15 + S23n(n3)*B14 + S13n(n3)*B13)*alpha(3,iElNum) - &
!!$                 (S33n(n4)*B15 + S23n(n4)*B14 + S13n(n4)*B13)*alpha(4,iElNum)
!!$            IF(myid.EQ.1.and.k3i.eq.33) print*,'proc2',- &
!!$                 (S33n(n1)*B15 + S23n(n1)*B14 + S13n(n1)*B13)*alpha(1,iElNum) - &
!!$                 (S33n(n2)*B15 + S23n(n2)*B14 + S13n(n2)*B13)*alpha(2,iElNum) - &
!!$                 (S33n(n3)*B15 + S23n(n3)*B14 + S13n(n3)*B13)*alpha(3,iElNum) - &
!!$                 (S33n(n4)*B15 + S23n(n4)*B14 + S13n(n4)*B13)*alpha(4,iElNum)

         ENDDO
      ENDDO

!      IF(myid.eq.0) print*,'***1procs',myid,Rnet(3*11)
!
!      IF(myid.EQ.0) print*,'***ll',myid,Rnet(3*3)
!      IF(myid.EQ.1) print*,'***ll',myid,Rnet(3*11)
      


!!$
!!$    SixInv = 1.d0/6.d0
!!$
!!$    DO i = 1, numnp
!!$       AhatInv = 1.d0/Ahat(i)
!!$
!!$       SumS11n = 0.d0
!!$       SumS22n = 0.d0
!!$       SumS33n = 0.d0
!!$       SumS12n = 0.d0
!!$       SumS23n = 0.d0
!!$       SumS13n = 0.d0
!!$       
!!$       DO j = 1, NumElNeigh(i)
!!$            
!!$          iElNum = ElConn(i,j)
!!$            
!!$          n1 = nodes(1,iElNum)
!!$          n2 = nodes(2,iElNum)
!!$          n3 = nodes(3,iElNum)
!!$          n4 = nodes(4,iElNum)
!!$
!!$          IF(i.EQ.n1) THEN
!!$             ix = 1
!!$          ELSE IF(i.EQ.n2)THEN
!!$             ix = 2
!!$          ELSE IF(i.EQ.n3)THEN
!!$             ix = 3
!!$          ELSE IF(i.EQ.n4)THEN
!!$             ix = 4
!!$          ENDIF
!!$            
!!$          k3n1 = 3*n1
!!$          k3n2 = 3*n2
!!$          k3n3 = 3*n3
!!$          k3n4 = 3*n4
!!$          
!!$          k2n1 = k3n1 - 1
!!$          k2n2 = k3n2 - 1
!!$          k2n3 = k3n3 - 1
!!$          k2n4 = k3n4 - 1
!!$          
!!$          k1n1 = k3n1 - 2
!!$          k1n2 = k3n2 - 2
!!$          k1n3 = k3n3 - 2
!!$          k1n4 = k3n4 - 2         
!!$          ! k#n# dummy variables replaces:
!!$          u1 = disp(k1n1)        ! (3*n1-2)
!!$          u2 = disp(k1n2)        ! (3*n2-2)
!!$          u3 = disp(k1n3)        ! (3*n3-2)
!!$          u4 = disp(k1n4)        ! (3*n4-2)
!!$          v1 = disp(k2n1)        ! (3*n1-1)
!!$          v2 = disp(k2n2)        ! (3*n2-1)
!!$          v3 = disp(k2n3)        ! (3*n3-1)
!!$          v4 = disp(k2n4)        ! (3*n4-1)
!!$          w1 = disp(k3n1)        ! (3*n1)
!!$          w2 = disp(k3n2)        ! (3*n2)
!!$          w3 = disp(k3n3)        ! (3*n3)
!!$          w4 = disp(k3n4)        ! (3*n4)
!!$           
!!$          x1 = coor(1,n1)
!!$          x2 = coor(1,n2)
!!$          x3 = coor(1,n3)
!!$          x4 = coor(1,n4)
!!$          y1 = coor(2,n1)
!!$          y2 = coor(2,n2)
!!$          y3 = coor(2,n3)
!!$          y4 = coor(2,n4)
!!$          z1 = coor(3,n1)
!!$          z2 = coor(3,n2)
!!$          z3 = coor(3,n3)
!!$          z4 = coor(3,n4)
!!$
!!$! See the maple worksheet 'V3D4.mws' for the derivation of [B] 
!!$
!!$          B1  = -(y3*z4 - y4*z3 - y2*z4 + y2*z3 + z2*y4 - z2*y3) * SixInv
!!$          B2  =  (x3*z4 - x4*z3 - x2*z4 + x2*z3 + z2*x4 - z2*x3) * SixInv
!!$          B3  = -(x3*y4 - x4*y3 - x2*y4 + x2*y3 + y2*x4 - y2*x3) * SixInv
!!$          B4  =  (y3*z4 - y4*z3 - y1*z4 + y1*z3 + z1*y4 - z1*y3) * SixInv
!!$          B5  = -(x3*z4 - x4*z3 - x1*z4 + x1*z3 + z1*x4 - z1*x3) * SixInv
!!$          B6  =  (x3*y4 - x4*y3 - x1*y4 + x1*y3 + y1*x4 - y1*x3) * SixInv
!!$          B7  = -(y2*z4 - z2*y4 - y1*z4 + y1*z2 + z1*y4 - z1*y2) * SixInv
!!$          B8  =  (x2*z4 - z2*x4 - x1*z4 + x1*z2 + z1*x4 - z1*x2) * SixInv
!!$          B9  = -(x2*y4 - y2*x4 - x1*y4 + x1*y2 + y1*x4 - y1*x2) * SixInv
!!$          B10 =  (y2*z3 - z2*y3 - y1*z3 + y1*z2 + z1*y3 - z1*y2) * SixInv
!!$          B11 = -(x2*z3 - z2*x3 - x1*z3 + x1*z2 + z1*x3 - z1*x2) * SixInv
!!$          B12 =  (x2*y3 - y2*x3 - x1*y3 + x1*y2 + y1*x3 - y1*x2) * SixInv
!!$
!!$! calculate the strain            
!!$!                        [E] = [B]{d}
!!$!
!!$          E11 = alpha(ix,iElNum)*(B1*u1 + B4*u2 + B7*u3 + B10*u4)
!!$          E22 = alpha(ix,iElNum)*(B2*v1 + B5*v2 + B8*v3 + B11*v4)
!!$          E33 = alpha(ix,iElNum)*(B3*w1 + B6*w2 + B9*w3 + B12*w4)
!!$          E12 = alpha(ix,iElNum)*(B2*u1 + B1*v1 + B5*u2 + B4*v2 + B8*u3 + B7*v3 + B11*u4 + B10*v4)
!!$          E23 = alpha(ix,iElNum)*(B3*v1 + B2*w1 + B6*v2 + B5*w2 + B9*v3 + B8*w3 + B12*v4 + B11*w4) 
!!$          E13 = alpha(ix,iElNum)*(B3*u1 + B1*w1 + B6*u2 + B4*w2 + B9*u3 + B7*w3 + B12*u4 + B10*w4)
!!$! calculate the stress            -1
!!$!                        [S] = [C]  {E}
!!$!       
!!$          SumS11n = SumS11n + E11*ci(1,1) + E22*ci(2,1) + E33*ci(4,1)
!!$          SumS22n = SumS22n + E11*ci(2,1) + E22*ci(3,1) + E33*ci(5,1)
!!$          SumS33n = SumS33n + E11*ci(4,1) + E22*ci(5,1) + E33*ci(6,1)
!!$          SumS12n = SumS12n + E12*ci(7,1)
!!$          SumS23n = SumS23n + E23*ci(8,1)
!!$          SumS13n = SumS13n + E13*ci(9,1)
!!$
!!$       ENDDO
!!$       S11n(i) = AhatInv*SumS11n
!!$       S22n(i) = AhatInv*SumS22n
!!$       S33n(i) = AhatInv*SumS33n
!!$       S12n(i) = AhatInv*SumS12n
!!$       S23n(i) = AhatInv*SumS23n
!!$       S13n(i) = AhatInv*SumS13n ! stress at node
!!$    ENDDO
!!$
!!$    DO j = 1, NumEl
!!$       
!!$       iElNum = j
!!$       
!!$       n1 = nodes(1,iElNum)
!!$       n2 = nodes(2,iElNum)
!!$       n3 = nodes(3,iElNum)
!!$       n4 = nodes(4,iElNum)
!!$
!!$       k3n1 = 3*n1
!!$       k3n2 = 3*n2
!!$       k3n3 = 3*n3
!!$       k3n4 = 3*n4
!!$       
!!$       k2n1 = k3n1 - 1
!!$       k2n2 = k3n2 - 1
!!$       k2n3 = k3n3 - 1
!!$       k2n4 = k3n4 - 1
!!$
!!$       k1n1 = k3n1 - 2
!!$       k1n2 = k3n2 - 2
!!$       k1n3 = k3n3 - 2
!!$       k1n4 = k3n4 - 2 
!!$
!!$       x1 = coor(1,n1)
!!$       x2 = coor(1,n2)
!!$       x3 = coor(1,n3)
!!$       x4 = coor(1,n4)
!!$       y1 = coor(2,n1)
!!$       y2 = coor(2,n2)
!!$       y3 = coor(2,n3)
!!$       y4 = coor(2,n4)
!!$       z1 = coor(3,n1)
!!$       z2 = coor(3,n2)
!!$       z3 = coor(3,n3)
!!$       z4 = coor(3,n4)
!!$
!!$! See the maple worksheet 'V3D4.mws' for the derivation of [B] 
!!$
!!$       B1  = -(y3*z4 - y4*z3 - y2*z4 + y2*z3 + z2*y4 - z2*y3) * SixInv
!!$       B2  =  (x3*z4 - x4*z3 - x2*z4 + x2*z3 + z2*x4 - z2*x3) * SixInv
!!$       B3  = -(x3*y4 - x4*y3 - x2*y4 + x2*y3 + y2*x4 - y2*x3) * SixInv
!!$       B4  =  (y3*z4 - y4*z3 - y1*z4 + y1*z3 + z1*y4 - z1*y3) * SixInv
!!$       B5  = -(x3*z4 - x4*z3 - x1*z4 + x1*z3 + z1*x4 - z1*x3) * SixInv
!!$       B6  =  (x3*y4 - x4*y3 - x1*y4 + x1*y3 + y1*x4 - y1*x3) * SixInv
!!$       B7  = -(y2*z4 - z2*y4 - y1*z4 + y1*z2 + z1*y4 - z1*y2) * SixInv
!!$       B8  =  (x2*z4 - z2*x4 - x1*z4 + x1*z2 + z1*x4 - z1*x2) * SixInv
!!$       B9  = -(x2*y4 - y2*x4 - x1*y4 + x1*y2 + y1*x4 - y1*x2) * SixInv
!!$       B10 =  (y2*z3 - z2*y3 - y1*z3 + y1*z2 + z1*y3 - z1*y2) * SixInv
!!$       B11 = -(x2*z3 - z2*x3 - x1*z3 + x1*z2 + z1*x3 - z1*x2) * SixInv
!!$       B12 =  (x2*y3 - y2*x3 - x1*y3 + x1*y2 + y1*x3 - y1*x2) * SixInv
!!$
!!$! local node 1
!!$       Rnet(k1n1) = Rnet(k1n1) - (S11n(n1)*B1 + S12n(n1)*B2 + S13n(n1)*B3)*alpha(1,iElNum)
!!$       Rnet(k2n1) = Rnet(k2n1) - (S22n(n1)*B2 + S12n(n1)*B1 + S23n(n1)*B3)*alpha(1,iElNum)
!!$       Rnet(k3n1) = Rnet(k3n1) - (S33n(n1)*B3 + S23n(n1)*B2 + S13n(n1)*B1)*alpha(1,iElNum)
!!$! local node 2 
!!$       Rnet(k1n2) = Rnet(k1n2) - (S11n(n2)*B4 + S12n(n2)*B5 + S13n(n2)*B6)*alpha(2,iElNum)
!!$       Rnet(k2n2) = Rnet(k2n2) - (S22n(n2)*B5 + S12n(n2)*B4 + S23n(n2)*B6)*alpha(2,iElNum)
!!$       Rnet(k3n2) = Rnet(k3n2) - (S33n(n2)*B6 + S23n(n2)*B5 + S13n(n2)*B4)*alpha(2,iElNum)
!!$! local node 3 
!!$       Rnet(k1n3) = Rnet(k1n3) - (S11n(n3)*B7 + S12n(n3)*B8 + S13n(n3)*B9 )*alpha(3,iElNum)
!!$       Rnet(k2n3) = Rnet(k2n3) - (S22n(n3)*B8 + S12n(n3)*B7 + S23n(n3)*B9 )*alpha(3,iElNum)
!!$       Rnet(k3n3) = Rnet(k3n3) - (S33n(n3)*B9 + S23n(n3)*B8 + S13n(n3)*B7 )*alpha(3,iElNum)
!!$! local node 4
!!$       Rnet(k1n4) = Rnet(k1n4) - (S11n(n4)*B10 + S12n(n4)*B11 + S13n(n4)*B12 )*alpha(4,iElNum)
!!$       Rnet(k2n4) = Rnet(k2n4) - (S22n(n4)*B11 + S12n(n4)*B10 + S23n(n4)*B12 )*alpha(4,iElNum)
!!$       Rnet(k3n4) = Rnet(k3n4) - (S33n(n4)*B12 + S23n(n4)*B11 + S13n(n4)*B10 )*alpha(4,iElNum)
!!$
            
!!$! ASSEMBLE THE INTERNAL FORCE VECTOR
!!$
!!$       Rnet(k1i) = Rnet(k1i) - &
!!$            (S11n(n1)*B13 + S12n(n1)*B14 + S13n(n1)*B15)*alpha(1,iElNum) - & ! local node 1
!!$            (S11n(n2)*B13 + S12n(n2)*B14 + S13n(n2)*B15)*alpha(2,iElNum) - & ! local node 2
!!$            (S11n(n3)*B13 + S12n(n3)*B14 + S13n(n3)*B15)*alpha(3,iElNum) - & ! local node 3
!!$            (S11n(n4)*B13 + S12n(n4)*B14 + S13n(n4)*B15)*alpha(4,iElNum) ! local node 4
!!$       
!!$       Rnet(k2i) = Rnet(k2i) - &
!!$            (S22n(n1)*B14 + S12n(n1)*B13 + S23n(n1)*B15)*alpha(1,iElNum) - &
!!$            (S22n(n2)*B14 + S12n(n2)*B13 + S23n(n2)*B15)*alpha(2,iElNum) - &
!!$            (S22n(n3)*B14 + S12n(n3)*B13 + S23n(n3)*B15)*alpha(3,iElNum) - &
!!$            (S22n(n4)*B14 + S12n(n4)*B13 + S23n(n4)*B15)*alpha(4,iElNum)
!!$
!!$       Rnet(k3i) = Rnet(k3i) - &
!!$            (S33n(n1)*B15 + S23n(n1)*B14 + S13n(n1)*B13)*alpha(1,iElNum) - &
!!$            (S33n(n2)*B15 + S23n(n2)*B14 + S13n(n2)*B13)*alpha(2,iElNum) - &
!!$            (S33n(n3)*B15 + S23n(n3)*B14 + S13n(n3)*B13)*alpha(3,iElNum) - &
!!$            (S33n(n4)*B15 + S23n(n4)*B14 + S13n(n4)*B13)*alpha(4,iElNum)
!!$       

!!$   ENDDO
!!$ENDDO

  END SUBROUTINE v3d4n

