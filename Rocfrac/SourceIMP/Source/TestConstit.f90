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
PROGRAM TestConst
! xlf90 -qsuffix=f=f90 -c v3d10_nl_huang.f90
! xlf90 -qsuffix=f=f90 -c  HUANG_Const_Model.f90
! xlf90 -qsuffix=f=f90 TestConstit.f90  HUANG_Const_Model.o v3d10_nl_huang.o

  INTEGER,parameter :: NumNP = 10
  INTEGER,parameter :: NumEL = 1
  integer,parameter :: numat_vol = 1
  integer, dimension(1:numat_vol) :: matcstet
  REAL*8, dimension(1:3,1:NumNP) :: Coor
  integer, dimension(1:10,NumEL) :: Conn
  REAL*8, dimension(1:3*NumNP) :: R_in, disp
   REAL*8, dimension(1:4,1:NumEL) :: S11,S22,S33,S12,S23,S13

   REAL*8 :: stretch = 0.001d0


  integer,parameter :: NSTATEV = 1
  integer,parameter :: NMATRIX = 3
  integer,parameter :: NPARTICLE = 4
  integer,parameter :: NPARTICLETYPE = 2
  integer,parameter :: NINTERFAC = 3

  real*8, dimension(1:NSTATEV) :: STATEV
  real*8, dimension(1:NMATRIX) :: MATRIX
  real*8, dimension(1:NPARTICLE,1:NPARTICLETYPE) :: PARTICLE
  real*8, dimension(1:NINTERFAC) :: INTERFAC


  real*8 :: time, dt, strain


  statev = 0.d0

  MATRIX(1) = 1.e6 ! : YOUNG'S MODULUS OF THE MATRIX
  MATRIX(2) = .499 ! : POISSON'S RATIO OF THE MATRIX
  MATRIX(3) = .073 ! : VOLUME FRACTION OF THE MATRIX
  
  PARTICLE(1,1) = 14.23e9    ! : YOUNG'S MODULUS OF TYPE-I PARTICLES
  PARTICLE(2,1) = .31   ! : POISSON'S RATIO OF TYPE-I PARTICLES
  PARTICLE(3,1) = .695    ! : VOLUME FRACTION OF TYPE-I PARTICLES
  PARTICLE(4,1) =125.0e-6    ! : RADIUS OF TYPE-I PARTICLES

  PARTICLE(1,2) = 14.23e9    ! : YOUNG'S MODULUS OF TYPE-I PARTICLES
  PARTICLE(2,2) = .31    ! : POISSON'S RATIO OF TYPE-I PARTICLES
  PARTICLE(3,2) =  .232   ! : VOLUME FRACTION OF TYPE-I PARTICLES
  PARTICLE(4,2) = 4.0e-6   ! : RADIUS OF TYPE-I PARTICLES

  INTERFAC(1) = 1.7e6     ! : STRENGTH OF THE INTERFACE
  INTERFAC(2) = 1550000000.0e6    ! : LINEAR MODULUS OF THE INTERFACE
  INTERFAC(3) = 15000.0e6 ! : SOFTENING MODULUS OF THE INTERFACE


  Conn(1,1) = 1
  Conn(2,1) = 3 
  Conn(3,1) = 2
  Conn(4,1) = 4
  Conn(5,1) = 5
  Conn(6,1) = 6
  Conn(7,1) = 7
  Conn(8,1) = 8
  Conn(9,1) = 9
  Conn(10,1) =10 

 
  Coor(1,1) = 0.000000000E+0;  Coor(2,1) = 1.000000000E+0;  Coor(3,1) =0.000000000E+0

  Coor(1,2) = 1.000000000E+0;  Coor(2,2) = 1.000000000E+0;  Coor(3,2) =1.000000000E+0

  Coor(1,3) = 1.000000000E+0;  Coor(2,3) = 0.000000000E+0;  Coor(3,3) =0.000000000E+0

  Coor(1,4) = 0.000000000E+0;  Coor(2,4) = 0.000000000E+0;  Coor(3,4) =1.000000000E+0

  Coor(1,5) = 5.000000000E-1;  Coor(2,5) = 5.000000000E-1;  Coor(3,5) =0.000000000E+0

  Coor(1,6) = 1.000000000E+0;  Coor(2,6) = 5.000000000E-1;  Coor(3,6) =5.000000000E-1

  Coor(1,7) = 5.000000000E-1;  Coor(2,7) = 1.000000000E+0;  Coor(3,7) =5.000000000E-1

  Coor(1,8) = 0.000000000E+0;  Coor(2,8) = 5.000000000E-1;  Coor(3,8) =5.000000000E-1

  Coor(1,9) = 5.000000000E-1;  Coor(2,9) = 0.000000000E+0;  Coor(3,9) =5.000000000E-1

  Coor(1,10) = 5.000000000E-1;  Coor(2,10) = 5.000000000E-1; Coor(3,10) =1.000000000E+0

  dt = .01d0
  time = 0.d0


  disp(:) = 0.d0

  matcstet(1) = 1 

  DO i = 1, 800

     time = time + dt


     disp(1*3-2) = -stretch*REAL(i)
     disp(1*3-1) =  stretch*REAL(i)
     disp(1*3  ) = -stretch*REAL(i)
     
     disp(2*3-2) = stretch*REAL(i)
     disp(2*3-1) = stretch*REAL(i)
     disp(2*3  ) = stretch*REAL(i)

     disp(3*3-2) = stretch*REAL(i)
     disp(3*3-1) = -stretch*REAL(i)
     disp(3*3  ) = -stretch*REAL(i)

     disp(4*3-2) = -stretch*REAL(i)
     disp(4*3-1) = -stretch*REAL(i)
     disp(4*3  ) =  stretch*REAL(i)
  
     disp(5*3-2) = 0.5d0*( disp(1*3-2) + disp(3*3-2))
     disp(5*3-1) = 0.5d0*( disp(1*3-1) + disp(3*3-1))
     disp(5*3  ) = 0.5d0*( disp(1*3  ) + disp(3*3  ))

     disp(6*3-2) = 0.5d0*( disp(2*3-2) + disp(3*3-2))
     disp(6*3-1) = 0.5d0*( disp(2*3-1) + disp(3*3-1))
     disp(6*3  ) = 0.5d0*( disp(2*3  ) + disp(3*3  ))

     disp(7*3-2) = 0.5d0*( disp(2*3-2) + disp(1*3-2))
     disp(7*3-1) = 0.5d0*( disp(2*3-1) + disp(1*3-1))
     disp(7*3  ) = 0.5d0*( disp(2*3  ) + disp(1*3  ))

     disp(8*3-2) = 0.5d0*( disp(4*3-2) + disp(1*3-2))
     disp(8*3-1) = 0.5d0*( disp(4*3-1) + disp(1*3-1))
     disp(8*3  ) = 0.5d0*( disp(4*3  ) + disp(1*3  ))

     disp(9*3-2) = 0.5d0*( disp(4*3-2) + disp(3*3-2))
     disp(9*3-1) = 0.5d0*( disp(4*3-1) + disp(3*3-1))
     disp(9*3  ) = 0.5d0*( disp(4*3  ) + disp(3*3  ))

     disp(10*3-2) = 0.5d0*( disp(4*3-2) + disp(2*3-2))
     disp(10*3-1) = 0.5d0*( disp(4*3-1) + disp(2*3-1))
     disp(10*3  ) = 0.5d0*( disp(4*3  ) + disp(2*3  ))


     R_in = 0.d0
     CALL V3D10_NL_HUANG(coor,matcstet,Conn,R_in,disp,&
          S11,S22,S33,S12,S23,S13,strain, &
          numnp,1,1,NumEl,numat_vol,&
          STATEV,NSTATEV,MATRIX,NMATRIX, &
          PARTICLE,NPARTICLE,NPARTICLETYPE,INTERFAC,NINTERFAC)

!!$     PRINT*,'time=',time
!!$     print*,S11(1:4,1)
!!$     print*,S22(1:4,1)
!!$     print*,S33(1:4,1)
!!$     print*,S12(1:4,1)
!!$     print*,S13(1:4,1)
!!$     print*,S23(1:4,1)

     write(12,*) strain, S11(1,1)
     
     

  ENDDO


end PROGRAM TestConst

