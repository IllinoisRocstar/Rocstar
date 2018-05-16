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
subroutine ArtificialDamping(numcstet, numnp, &
     rho, cd_fastest, DetFold, dt, F, Fdot, Vx6, StrssVisco)

  implicit none

  integer :: numcstet
  integer :: numnp

  REAL*8 :: rho
  real*8 :: cd_fastest
  real*8 :: DetFold ! change to more general
  real*8 :: dt
  real*8 :: Vx6

  
  REAL*8,DIMENSION(1:3,1:3) :: Id = RESHAPE( &
       (/1.0,0.0, 0.0, &
         0.0,1.0, 0.0, &
         0.0,0.0, 1.0 /),(/3,3/) )


  REAL*8 :: F(3,3),DetF,Finv(3,3), Fdot(3,3), T(3,3), symT(3,3), devT(3,3)
  REAL*8 :: StrssCauchy(3,3), StrssVisco(3,3)

  REAL*8 :: eta, c1, cL, h, Du

  c1 = 1.d0
  cL = 0.1d0 

  Finv(:,:) = F(:,:)

! returns the                           -1
! inverse of the deformation gradient, F  
! and determanate of F

  CALL invert3x3(Finv,DetF)

!   .    -1
!   F * F
!

  T = MATMUL(Fdot,Finv)
!
!        .    -1                         T
! symm(  F * F  ) = symm(T) = 0.5*( T + T )

  symT = 0.5d0*( T + TRANSPOSE(T))

!
! dev( symm(T) ) = symmT - 1/3 * ( tr(symm(T)) * I
!

  devT = symT - 1.d0/3.d0*(symT(1,1)+symT(2,2)+symT(3,3))*Id

!
!  h = ( J* d! * |Vol| ) ^ (1/d)
!  
!       where
!             d is the dimension fo space
!            Vol is the volume of the element
!            J is the Jacobian
!
!

  h = ( DetF * Vx6 )**(1./3.)

!                     n+1        n
!   DelU = h * ( log(J    - log J  ) ) / Delt
!

  Du = h * ( LOG(DetF) - LOG(DetFOld) )/dt



!
!   artifical viscosity,  max(0, -3/4 * h * rho*(c1*DelU - cL*a)

!  IF(Du.GE.0.d0)THEN
!     eta = 0.d0
!  ELSE

  eta = MAX(0.d0,-0.75d0*h*rho/DetF*(c1*Du - cL*cd_fastest))

!  ENDIF

! Cauchy Stress

  StrssCauchy = 2.0*eta*devT

! PK2 stress
     
  StrssVisco = detF*Finv*StrssCauchy*TRANSPOSE(Finv)

! 
! store Jacobian for next time step

  DetFOld = DetF

end subroutine ArtificialDamping

