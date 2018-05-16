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
SUBROUTINE Constant_W(Gm, G2, K2, Km, a_eta, a_zeta, vm, v2, W)


  IMPLICIT NONE

  REAL*8, DIMENSION(1:6,1:6) :: W
! Poisson's ratio
  REAL*8 :: vm, v2
! Shear modulus
  REAL*8 :: Gm, G2
! Bulk modulus
  REAL*8 :: K2, Km, a_eta, a_zeta


  REAL*8 :: e, f, h
  REAL*8 :: A2, B2

  REAL*8 :: intparam
  REAL*8 :: upsilon
  REAL*8 :: gamma

  integer :: i


  e = G2/a_eta
  f = G2/a_zeta
  h = Gm/G2


  CALL Constants_A2_B2(e, f, h, vm, v2, A2, B2)

  intparam = K2/(a_eta)

  upsilon = G2*(21.d0*A2 + 5.d0*B2)/(5.d0*Gm)

  gamma = K2*(3.d0*Km+4.d0*Gm)/( 3.d0*Km*(3.d0*K2 + 4.d0*Gm*(1.d0+3.d0*intparam)) )

! W
  W(:,:) = 0.d0

  W(1,1) = gamma + 2.d0/3.d0*upsilon
  W(1,2) = gamma - 1.d0/3.d0*upsilon
  W(1,3) = W(1,2)
  W(2,1) = W(1,2)
  W(2,2) = W(1,1)
  W(2,3) = W(1,2)
  W(3,1) = W(1,3)
  W(3,2) = W(2,3)
  W(3,3) = W(1,1)
  W(4,4) = upsilon
  W(5,5) = W(4,4)
  W(6,6) = W(4,4)

END SUBROUTINE Constant_W
  

SUBROUTINE Constants_A2_B2(e,f,h,vm,v2,A2,B2)


  IMPLICIT NONE

  REAL*8 :: vm,v2
  REAL*8 :: e, f, h

  REAL*8 :: A2,B2
  

  A2 = -120.d0*(-f+f*vm+e-e*vm)/(-196.d0-196.d0*v2*h*e-224.d0*f*v2*e+166.d0*h*f*v2 &
       -192.d0*v2*h**2*e+364.d0*h*e*vm-255.d0*v2*h*vm-190.d0*v2*h**2*vm+560.d0*f*e*vm+140.d0 &
       *v2*h*e*vm-50.d0*h*f*v2*vm+240.d0*v2*h**2*e*vm-336.d0*h*f*v2*e+266.d0*f*h*vm &
       -160.d0*v2*e*vm+80.d0*f*v2*vm-308.d0*e*h+224.d0*v2*e-273.d0*h-392.d0*e-392.d0*f-784.d0*f*e &
       +112.d0*v2-56.d0*h**2-238.d0*h*f+152.d0*v2*h**2+261.d0*v2*h-112.d0*f*v2+140.d0*vm+315.d0*h &
       *vm+160.d0*f*v2*e*vm+280.d0*f*vm+280.d0*e*vm+70.d0*h**2*vm-80.d0*v2*vm+240.d0*h*f*v2 &
       *e*vm)

  B2 = 15.d0*(-24.d0*v2*h*e-7.d0*h+19.d0*v2*h-28.d0+16.d0*v2-56.d0*f+24.d0*v2*h*e*vm+16.d0*f* &
       v2*vm+28.d0*vm+56.d0*f*vm+7.d0*h*vm-16.d0*f*v2-19.d0*v2*h*vm-16.d0*v2*vm)/(-196.d0-196.d0* &
       v2*h*e-224.d0*f*v2*e+166.d0*h*f*v2-192.d0*v2*h**2*e+364.d0*h*e*vm-255.d0*v2*h*vm- &
       190.d0*v2*h**2*vm+560.d0*f*e*vm+140.d0*v2*h*e*vm-50.d0*h*f*v2*vm+240.d0*v2*h**2*e &
       *vm-336.d0*h*f*v2*e+266.d0*f*h*vm-160.d0*v2*e*vm+80.d0*f*v2*vm-308.d0*e*h+224.d0*v2* &
       e-273.d0*h-392.d0*e-392.d0*f-784.d0*f*e+112.d0*v2-56.d0*h**2-238.d0*h*f+152.d0*v2*h**2+261.d0 &
       *v2*h-112.d0*f*v2+140.d0*vm+315.d0*h*vm+160.d0*f*v2*e*vm+280.d0*f*vm+280.d0*e*vm+70.d0* &
       h**2*vm-80.d0*v2*vm+240.d0*h*f*v2*e*vm)


END SUBROUTINE Constants_A2_B2


!!$
!!$      t1 = e*h
!!$      t3 = 56*e
!!$      t4 = v2*e
!!$      t5 = 32*t4
!!$      t6 = f*e
!!$      t7 = 112*t6
!!$      t8 = v2*h
!!$      t9 = t8*e
!!$      t11 = f*v2
!!$      t12 = t11*e
!!$      t13 = 32*t12
!!$      t14 = h**2
!!$      t15 = v2*t14 !
!!$      t16 = t15*e
!!$      t17 = 24*t16
!!$      t18 = h*f !
!!$      t19 = t18*t4
!!$      t20 = 48*t19
!!$      t21 = t18*v2
!!$      t23 = 21*h
!!$      t24 = 14*t18
!!$      t25 = 7*t14
!!$      t26 = 16*v2
!!$      t27 = 56*f
!!$      t28 = 19*t15
!!$      t29 = 3*t8
!!$      t30 = 16*t11
!!$      t31 = -28*t1+t3-t5+t7+28*t9+t13-t17+t20+28-58*t21-t23-t24-t25-t26+
!!$     #t27+t28-t29+t30
!!$      t36 = e*vm !
!!$      t38 = t4*vm !
!!$      t44 = f*vm
!!$      t47 = v2*vm
!!$      t48 = t18*t47
!!$      t50 = t8*t36
!!$      t58 = -196-392*e-273*h+112*v2-392*f+280*t36+240*t18*t38+152*t15-56
!!$     #*t14-308*t1+280*t44-336*t19-50*t48+140*t50+160*t11*t36+240*t15*t36
!!$     #+140*vm-238*t18 ! 
!!$      t64 = t11*vm
!!$      t68 = t18*vm
!!$      t72 = t8*vm
!!$      t74 = t1*vm
!!$      t76 = h*vm
!!$
!!$      t85 = 166*t21-192*t16-224*t12-196*t9-160*t38+80*t64-190*t15*vm+266 &     
!!$           *t68+560*t6*vm-255*t72+364*t74+315*t76-80*t47+70*t14*vm-112*t11-784 &
!!$           *t6+261*t8+224*t4
!!$
!!$      t87 = 1/(t58+t85) !
!!$      t91 = 28+t27-t17-t25-56*t1+t28+t3-t23-t29-t26+t30
!!$      t98 = t24+t20-20*t48+20*t50+28*t74-28*t68+t7-t5+8*t9-38*t21+t13
!!$      t115 = -24*t9-7*h+19*t8-28+t26-t27+24*t50+16*t64+28*vm+56*t44+7*t76 &
!!$           -t30-19*t72-16*t47
!!$
!!$      A2 = -120*(-f+t44+e-t36)*t87
!!$      B2 = 15*t115*t87




