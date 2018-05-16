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
SUBROUTINE principal_stress(s11,s22,s33,&
     s12,s23,s13, &
     istrgss,NumElVol,SVonMises)

!!****f* Rocfrac/Rocfrac/Source/principal_stress.f90
!!
!!  NAME
!!     principal_stress
!!
!!  FUNCTION
!!
!!   Computes Principal Values of Symmetric Second Rank Tensor
!!
!!   INPUTS
!!
!!        S = Symmetric Second-Rank Tensor Stored as a Vector
!!        P = Principal Values
!!
!!     .. The Components of S Must be Stored in the Following Orders
!!
!!        2-D Problems,  S11,S12,S22
!!        3-D Problems,  S11,S12,S13,S22,S23,S33
!!
!!   OUTPUTS
!!       SVonMises -- VonMises Stress
!!
!!****



  IMPLICIT NONE

  INTEGER :: NumElVol
  INTEGER :: i,j, istrgss
  REAL*8, DIMENSION(1:NumElVol) :: SVonMises
  REAL*8, DIMENSION(1:6) :: S
  REAL*8, DIMENSION(1:istrgss,1:NumElVol) :: s11, s22, s33, s12, s23, s13
  
  REAL*8 :: prin1, prin2, prin3
  
  REAL*8 :: r,x,y,z,t,u,a
  
  REAL*8 :: RT2 = 1.414213562373090
  REAL*8 :: PI23 = 2.094395102393210
  
  DO j = 1, NumElVol
     
     prin1 = 0.0
     prin2 = 0.0
     prin3 = 0.0
     
     DO i = 1, istrgss

        s(1) = s11(i,j)
        s(2) = s12(i,j)
        s(3) = s13(i,j)
        s(4) = s22(i,j)
        s(5) = s23(i,j)
        s(6) = s33(i,j)

!.... 3-D Problem
 
        r = 0.0
        x = (S(1)+S(4)+S(6))/3.0
        y = S(1)*(S(4)+S(6))+S(4)*S(6)-S(2)*S(2)-S(3)*S(3)-S(5)*S(5)
        z = S(1)*S(4)*S(6)+2.0*S(2)*S(3)*S(5)-S(1)*S(5)*S(5) &
             -S(4)*S(3)*S(3)-S(6)*S(2)*S(2)
        t = 3.0*x*x-y
        u = 0.0
        IF(t.lt.1.0e-7.AND.t.gt.-1.0e-7) GOTO 20
        u = SQRT(2.0*t/3.0)
        a = (z + (t-x*x)*x)*rt2/u**3
        r = SQRT(ABS(1.d0 - a*a))
        r = DATAN2(r,a)/3.d0
     
20      continue
        prin1 = prin1 + x + u*rt2*COS(r)
        prin2 = prin2 + x + u*rt2*COS(r - pi23)
        prin3 = prin3 + x + u*rt2*COS(r + pi23)
        
     ENDDO
      
  ! average the guass points for the 10 node tetrahedra

     prin1 = prin1/float(istrgss)
     prin2 = prin2/float(istrgss)
     prin3 = prin3/float(istrgss)

! Von Mises' equivalent stress
  
     SVonMises(j) = SQRT((prin1-prin2)**2 + (prin2-prin3)**2 + (prin3-prin1)**2)/SQRT(2.d0)

  ENDDO


  RETURN
END SUBROUTINE principal_stress

