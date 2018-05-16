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
SUBROUTINE v3d4n_nl_arruda_boyce(numnp, numel, coor, disp, nodes, Rnet, &
     S11ev, S22ev, S33ev, S12ev, S23ev, S13ev, &
     NumElNeigh,ElConn,alpha, Ahat,&
     xmu,xkappa,NumMatVol, &
     nprocs,TotNumNdComm,TotNumNeighProcs,NeighProcList,NumNdComm,neigh_lst,&
     MPI_STATUS_SIZE,MPI_COMM_ROCFRAC,MPI_DOUBLE_PRECISION, &
     ReqRcv,ReqSnd,StatRcv, StatSnd)

    USE ROCSTAR_RocFracComm 

    IMPLICIT NONE

    INTEGER :: k1, j1, k2

    INTEGER :: TotNumNeighProcs, Nprocs,TotNumNdComm
    INTEGER, DIMENSION(1:TotNumNeighProcs) :: NeighProcList
    INTEGER, DIMENSION(1:TotNumNeighProcs) ::NumNdComm

    TYPE(rcv_buf), pointer, DIMENSION(:) :: RecvDataFrm 
    TYPE(send_buf), DIMENSION(1:TotNumNeighProcs) :: neigh_lst  

    integer :: NumMatVol
    INTEGER :: numnp, numel
    REAL*8, DIMENSION(1:numnp) :: Ahat
    REAL*8, DIMENSION(1:3*numnp) :: disp, Rnet
    REAL*8, DIMENSION(1:3,1:numnp) :: coor
    integer, DIMENSION(1:numnp) ::NumElNeigh
    INTEGER, DIMENSION(1:numnp,1:40) :: ElConn  ! fix 40 should not be hard coded
    INTEGER,DIMENSION(1:4,1:numel) :: nodes
    real*8, dimension(1:4,1:numel) :: alpha
    real*8, dimension(1:NumMatVol) :: xmu, xkappa
    REAL*8, DIMENSION(1:numel) :: S11ev, S22ev, S33ev, S12ev, S23ev, S13ev
    
    INTEGER :: j,i,k
    REAL*8 :: aaa
!--   coordinate holding variable
    REAL*8 :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
    REAL*8 :: x1d,x2d,x3d,x4d,y1d,y2d,y3d,y4d,z1d,z2d,z3d,z4d
    INTEGER :: k1n1,k1n2,k1n3,k1n4,k2n1,k2n2,k2n3,k2n4
    INTEGER :: k3n1,k3n2,k3n3,k3n4

    REAL*8 :: AhatInv
    INTEGER :: iElNum
    REAL*8 :: SixInv
!--   x, y, and z displacements of nodes
    REAL*8 :: u1,u2,u3,u4,v1,v2,v3,v4,w1,w2,w3,w4
!--   partial derivatives of the displacement 
    REAL*8 :: dudx,dvdy,dwdz,dudy,dvdx,dvdz,dwdy,dudz,dwdx

!--  Coordinate subtractions
    REAL*8 :: x14, x24, x34, y14, y24, y34, z14, z24, z34
!-- Added these to speed up B calculation
    
    REAL*8 ::  z12, z13,x12,x13,y12,y13
    REAL*8 :: C11, C12, C13, C21, C22, C23, C31, C32, C33
    REAL*8 :: SH1, SH2, SH3, SH4, SH5, SH6, SH7, SH8, SH9, SH10, SH11, SH12
    REAL*8 :: SH13, SH14, SH15
!   --   6*volume and the volume      
    REAL*8 :: Vx6
!--   spacial derivatives
    REAL*8 :: B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12
!--   strains
    REAL*8 :: E11,E22,E33,E12,E23,E13
    REAL*8 :: B13, B14, B15
    REAL*8 :: SumS11n,SumS22n,SumS33n,SumS12n,SumS23n,SumS13n
    REAL*8 :: SumE11n,SumE22n,SumE33n,SumE12n,SumE23n,SumE13n
    REAL*8, DIMENSION(1:numnp) :: S11n, S22n, S33n,S12n, S23n, S13n
    
    REAL*8 :: E11n, E22n, E33n, E12n, E23n, E13n
    INTEGER :: k3i,k2i,k1i
    INTEGER :: ix
    INTEGER :: n1, n2, n3, n4
    REAL*8 :: VolEl, VolEl0,Vx6inv,vol
    REAL*8 :: S11e, S22e, S33e, S12e, S23e, S13e
    REAL*8 :: VaInv

    REAL*8 :: coeff1, coeff2
    REAL*8 :: Jacob
    REAL*8 :: F11, F12, F13, F21, F22, F23, F31, F32, F33
    REAL*8,DIMENSION(1:numnp) :: F11v, F12v, F13v, F21v, F22v, F23v, F31v, F32v, F33v, voldefv
    REAL*8 :: Vx6def, Vx6invdef, voldef


    INTEGER :: myid
    INTEGER :: ierr

    REAL*8, pointer, dimension(:) :: buf

    integer :: MPI_STATUS_SIZE,MPI_COMM_ROCFRAC,MPI_DOUBLE_PRECISION

!--   Non-block receive, Non-block send request arrays
    INTEGER, DIMENSION(1:TotNumNeighProcs) :: ReqRcv, ReqSnd
    INTEGER, DIMENSION(1:MPI_STATUS_SIZE,1:TotNumNeighProcs) :: StatRcv, StatSnd

    SixInv = 1.d0/6.d0

    DO i = 1, numnp ! for each node

! /* Undeformed Volume of node */

       VaInv = 1.d0/Ahat(i)

 ! /* Initialize to 0, Va, Va0, Fa */

       F11 = 0.d0
       F22 = 0.d0
       F33 = 0.d0
       F12 = 0.d0
       F13 = 0.d0
       F21 = 0.d0
       F23 = 0.d0
       F31 = 0.d0
       F32 = 0.d0 

       voldef = 0.d0
              
       DO j = 1, NumElNeigh(i) ! for each element associated with node i
            
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

          x12 = x1 - x2 ! not used in vol. calc
          x13 = x1 - x3 ! not used in vol. calc
          x14 = x1 - x4
          x24 = x2 - x4
          x34 = x3 - x4
          y12 = y1 - y2 ! not used in vol. calc
          y13 = y1 - y3 ! not used in vol. calc
          y14 = y1 - y4
          y24 = y2 - y4
          y34 = y3 - y4
          z12 = z1 - z2 ! not used in vol. calc
          z13 = z1 - z3 ! not used in vol. calc
          z14 = z1 - z4
          z24 = z2 - z4
          z34 = z3 - z4
          
          c11 =    y24*z34 - z24*y34
          c21 = -( x24*z34 - z24*x34 )
          c31 =    x24*y34 - y24*x34
       
          Vx6 = -( x14*c11 + y14*c21 + z14*c31 )

          Vx6inv = 1.d0/Vx6

          VolEl0 = Vx6/6.d0 ! undeformed volume of element (Ve_O) 

! Compute the Shape functions
! NOTE: Factored for a more equivalent/compact form then maple's

          B1  = (y34*z24 - y24*z34) * Vx6inv
          B2  = (z34*x24 - z24*x34) * Vx6inv
          B3  = (x34*y24 - x24*y34) * Vx6inv
          B4  = (y13*z14 - y14*z13) * Vx6inv
          B5  = (z13*x14 - z14*x13) * Vx6inv
          B6  = (x13*y14 - x14*y13) * Vx6inv
          B7  = (y14*z12 - y12*z14) * Vx6inv
          B8  = (z14*x12 - z12*x14) * Vx6inv
          B9  = (x14*y12 - x12*y14) * Vx6inv
          B10 = (y12*z13 - y13*z12) * Vx6inv
          B11 = (z12*x13 - z13*x12) * Vx6inv
          B12 = (x12*y13 - x13*y12) * Vx6inv
! 
! deformation gradients F
!
          F11 = F11 + alpha(ix,IElNum)*VolEl0*(1.d0 + ( B1*u1 + B4*u2 + B7*u3 + B10*u4 )) ! 1 + ( dudx )
          F22 = F22 + alpha(ix,IElNum)*VolEl0*(1.d0 + ( B2*v1 + B5*v2 + B8*v3 + B11*v4 )) ! 1 + ( dvdy )
          F33 = F33 + alpha(ix,IElNum)*VolEl0*(1.d0 + ( B3*w1 + B6*w2 + B9*w3 + B12*w4 )) ! 1 + ( dwdz )
          F12 = F12 + alpha(ix,IElNum)*VolEl0*(B2*u1 + B5*u2 + B8*u3 + B11*u4) ! dudy
          F21 = F21 + alpha(ix,IElNum)*VolEl0*(B1*v1 + B4*v2 + B7*v3 + B10*v4) ! dvdx
          F23 = F23 + alpha(ix,IElNum)*VolEl0*(B3*v1 + B6*v2 + B9*v3 + B12*v4) ! dvdz
          F32 = F32 + alpha(ix,IElNum)*VolEl0*(B2*w1 + B5*w2 + B8*w3 + B11*w4) ! dwdy
          F13 = F13 + alpha(ix,IElNum)*VolEl0*(B3*u1 + B6*u2 + B9*u3 + B12*u4) ! dudz
          F31 = F31 + alpha(ix,IElNum)*VolEl0*(B1*w1 + B4*w2 + B7*w3 + B10*w4) ! dwdx

       ENDDO

       F11v(i) = F11*VaInv
       F22v(i) = F22*VaInv
       F33v(i) = F33*VaInv
       F12v(i) = F12*VaInv
       F13v(i) = F13*VaInv
       F21v(i) = F21*VaInv
       F23v(i) = F23*VaInv
       F31v(i) = F31*VaInv
       F32v(i) = F32*VaInv

       voldefv(i) = voldef

    ENDDO



    ALLOCATE(RecvDataFrm(0:nprocs-1))
!
!----- FORM THE BUFFER CONTAINING COMMUNICATED STRESS MATRIX NODAL VALUES
    
    ALLOCATE(buf(1:TotNumNdComm/3*10)) ! is really (TotNumNdComm*3 * 3 + 1)
    k1 = 1
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)
       ALLOCATE(RecvDataFrm(k)%rcvbuf(1:NumNdComm(j1)*10))
       DO j = 1, NumNdComm(j1)
          k2 = neigh_lst(j1)%NdId(j) !NdCommList(k)%NdId(j)
          buf(k1  ) = F11v(k2)
          buf(k1+1) = F12v(k2)
          buf(k1+2) = F13v(k2)
          buf(k1+3) = F21v(k2)
          buf(k1+4) = F22v(k2)
          buf(k1+5) = F23v(k2)
          buf(k1+6) = F31v(k2)
          buf(k1+7) = F32v(k2)
          buf(k1+8) = F33v(k2)
          buf(k1+9) = voldefv(k2)
          k1 = k1 + 10
       ENDDO
    ENDDO
!
!-MPI- RECEIVE THE RECIPRICAL MASS MATRIX DIAGONAL FROM THE NEIGHBORS
!
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)
       CALL MPI_IRECV(RecvDataFrm(k)%rcvbuf(1),NumNdComm(j1)*10, &
            MPI_DOUBLE_PRECISION, k, 10, MPI_COMM_ROCFRAC,ReqRcv(j1),ierr)
    ENDDO
!
!-MPI- SEND THE RECIPRICAL MASS MATRIX DIAGONAL TO THE NEIGHBORS
!
    k2 = 1
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)
       
       CALL MPI_ISEND(buf(k2),NumNdComm(j1)*10,&
            MPI_DOUBLE_PRECISION,k,10,MPI_COMM_ROCFRAC,ReqSnd(j1),ierr)
       k2 = k2 + NumNdComm(j1)*10
    ENDDO
!
!-MPI- WAIT FOR INTERNAL FORCE VECTOR COMMUNICATION TO COMPLETE
!
    IF(TotNumNeighProcs.GT.0)THEN
       CALL MPI_WAITALL(TotNumNeighProcs,ReqRcv,StatRcv,ierr)
       CALL MPI_WAITALL(TotNumNeighProcs,ReqSnd,StatSnd,ierr)
    ENDIF
    DEALLOCATE(buf)

!
!----- ADD NEIGHBOR'S CONTRIBUTION TO THE INTERNAL FORCE VECTOR
!
    DO j1 = 1, TotNumNeighProcs
       k = NeighProcList(j1)
       k1 = 1
       DO j = 1, NumNdComm(j1)
          k2 = neigh_lst(j1)%NdId(j)
          F11v(k2)  = F11v(k2)  + RecvDataFrm(k)%rcvbuf(k1)
          F12v(k2)  = F12v(k2)  + RecvDataFrm(k)%rcvbuf(k1+1)
          F13v(k2)  = F13v(k2)  + RecvDataFrm(k)%rcvbuf(k1+2)
          F21v(k2)  = F21v(k2)  + RecvDataFrm(k)%rcvbuf(k1+3)
          F22v(k2)  = F22v(k2)  + RecvDataFrm(k)%rcvbuf(k1+4)
          F23v(k2)  = F23v(k2)  + RecvDataFrm(k)%rcvbuf(k1+5)
          F31v(k2)  = F31v(k2)  + RecvDataFrm(k)%rcvbuf(k1+6)
          F32v(k2)  = F32v(k2)  + RecvDataFrm(k)%rcvbuf(k1+7)
          F33v(k2)  = F33v(k2)  + RecvDataFrm(k)%rcvbuf(k1+8)
          voldefv(k2) = voldefv(k2)  + RecvDataFrm(k)%rcvbuf(k1+9)
          k1 = k1 + 10
       ENDDO
    ENDDO

    DEALLOCATE(RecvDataFrm)
    DO i = 1, numnp ! for each node

       F11 = F11v(i)
       F22 = F22v(i)
       F33 = F33v(i)
       F12 = F12v(i)
       F13 = F13v(i)
       F21 = F21v(i)
       F23 = F23v(i)
       F31 = F31v(i)
       F32 = F32v(i)

       voldef = voldefv(i)

! WRONG FOR MORE THEN ONE MATERIAL

       CALL ARRUDA_BOYCE_CAUCHY(F11,F12,F13,F21,F22,F23,F31,F32,F33, &
            S11n(i),S22n(i),S33n(i),S12n(i),S13n(i),S23n(i),i,xmu(1),xkappa(1))
    ENDDO


    DO iElNum = 1, numel  ! For each element

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
       
       ! k#n# dummy variables replaces:
       u1 = disp(k1n1)           ! (3*n1-2)
       u2 = disp(k1n2)           ! (3*n2-2)
       u3 = disp(k1n3)           ! (3*n3-2)
       u4 = disp(k1n4)           ! (3*n4-2)
       v1 = disp(k2n1)           ! (3*n1-1)
       v2 = disp(k2n2)           ! (3*n2-1)
       v3 = disp(k2n3)           ! (3*n3-1)
       v4 = disp(k2n4)           ! (3*n4-1)
       w1 = disp(k3n1)           ! (3*n1)
       w2 = disp(k3n2)           ! (3*n2)
       w3 = disp(k3n3)           ! (3*n3)
       w4 = disp(k3n4)           ! (3*n4)
       
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
       
       x1 = x1 + u1
       x2 = x2 + u2
       x3 = x3 + u3
       x4 = x4 + u4
       y1 = y1 + v1
       y2 = y2 + v2
       y3 = y3 + v3
       y4 = y4 + v4
       z1 = z1 + w1
       z2 = z2 + w2
       z3 = z3 + w3
       z4 = z4 + w4

       x12 = x1 - x2 ! not used in vol. calc
       x13 = x1 - x3 ! not used in vol. calc
       x14 = x1 - x4
       x24 = x2 - x4
       x34 = x3 - x4
       y12 = y1 - y2 ! not used in vol. calc
       y13 = y1 - y3 ! not used in vol. calc
       y14 = y1 - y4
       y24 = y2 - y4
       y34 = y3 - y4
       z12 = z1 - z2 ! not used in vol. calc
       z13 = z1 - z3 ! not used in vol. calc
       z14 = z1 - z4
       z24 = z2 - z4
       z34 = z3 - z4
            
       c11 =    y24*z34 - z24*y34
       c21 = -( x24*z34 - z24*x34 )
       c31 =    x24*y34 - y24*x34
       
       Vx6def = -( x14*c11 + y14*c21 + z14*c31 )
       
       Vx6invdef = 1.d0 / Vx6def 

! calculate the volume
       voldef = Vx6def/6.d0
            
       IF(Vx6def.LE.0.d0) THEN
          WRITE(*,100) i
          STOP
       ENDIF

    ! equation (14)

       S11e = 0.d0
       S22e = 0.d0
       S33e = 0.d0
       S12e = 0.d0
       S23e = 0.d0
       S13e = 0.d0
       
       DO k = 1, 4
          S11e = S11e + alpha(k,IElNum)*S11n(nodes(k,iElNum))
          S22e = S22e + alpha(k,IElNum)*S22n(nodes(k,iElNum))
          S33e = S33e + alpha(k,IElNum)*S33n(nodes(k,iElNum))
          S12e = S12e + alpha(k,IElNum)*S12n(nodes(k,iElNum))
          S23e = S23e + alpha(k,IElNum)*S23n(nodes(k,iElNum))
          S13e = S13e + alpha(k,IElNum)*S13n(nodes(k,iElNum))   
       ENDDO

       S11ev(iElNum) = S11e
       S22ev(iElNum) = S22e
       S33ev(iElNum) = S33e
       S12ev(iElNum) = S12e
       S23ev(iElNum) = S23e
       S13ev(iElNum) = S13e

       SH1  = (y34*z24 - y24*z34) * Vx6invdef
       SH2  = (z34*x24 - z24*x34) * Vx6invdef
       SH3  = (x34*y24 - x24*y34) * Vx6invdef
       SH4  = (y13*z14 - y14*z13) * Vx6invdef
       SH5  = (z13*x14 - z14*x13) * Vx6invdef
       SH6  = (x13*y14 - x14*y13) * Vx6invdef
       SH7  = (y14*z12 - y12*z14) * Vx6invdef
       SH8  = (z14*x12 - z12*x14) * Vx6invdef
       SH9  = (x14*y12 - x12*y14) * Vx6invdef
       SH10 = (y12*z13 - y13*z12) * Vx6invdef
       SH11 = (z12*x13 - z13*x12) * Vx6invdef
       SH12 = (x12*y13 - x13*y12) * Vx6invdef

! ASSEMBLE THE INTERNAL FORCE VECTOR
!
! local node 1

       Rnet(k1n1) = Rnet(k1n1) - voldef* &
            ( S11e*SH1 + S12e*SH2 + S13e*SH3 )
       Rnet(k2n1) = Rnet(k2n1) - voldef* &
            ( S12e*SH1 + S22e*SH2 + S23e*SH3 )
       Rnet(k3n1) = Rnet(k3n1) - voldef* &
            ( S13e*SH1 + S23e*SH2 + S33e*SH3 )
! local node 2 
       Rnet(k1n2) = Rnet(k1n2) - voldef* &
            ( S11e*SH4 + S12e*SH5 + S13e*SH6 )
       Rnet(k2n2) = Rnet(k2n2) - voldef* &
            ( S12e*SH4 + S22e*SH5 + S23e*SH6 )
       Rnet(k3n2) = Rnet(k3n2) - voldef* &
            ( S13e*SH4 + S23e*SH5 + S33e*SH6 )
! local node 3 
       Rnet(k1n3) = Rnet(k1n3) - voldef* &
            ( S11e*SH7 + S12e*SH8 + S13e*SH9 )
       Rnet(k2n3) = Rnet(k2n3) - voldef* &
            ( S12e*SH7 + S22e*SH8 + S23e*SH9 )
       Rnet(k3n3) = Rnet(k3n3) - voldef* &
            ( S13e*SH7 + S23e*SH8 + S33e*SH9 )
! local node 4  
       Rnet(k1n4) = Rnet(k1n4) - voldef* &
            ( S11e*SH10 + S12e*SH11 + S13e*SH12 )
       Rnet(k2n4) = Rnet(k2n4) - voldef* &
            ( S12e*SH10 + S22e*SH11 + S23e*SH12 )
       Rnet(k3n4) = Rnet(k3n4) - voldef* &
            ( S13e*SH10 + S23e*SH11 + S33e*SH12 )
    
    ENDDO

    RETURN
100 FORMAT(' Negative Jacobian for element: ',i10)

  END SUBROUTINE v3d4n_nl_arruda_boyce

