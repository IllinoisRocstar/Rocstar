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
SUBROUTINE ARRUDA_BOYCE_CAUCHY(F11,F12,F13,F21,F22,F23,F31,F32,F33,&
     Cchy11,Cchy22,Cchy33,Cchy12,Cchy13,Cchy23,ielem,mu,kappa)

!!****f* Rocfrac/Rocfrac/Source/arruda_boyce.cauchy.f90
!!
!!  NAME
!!     ARRUDA_BOYCE_CAUCHY
!!
!!  FUNCTION
!!     Arruda-Boyce constitutive model, returns Cauchy Stress
!!
!!  INPUTS
!!     F11,F12,F13,F21,F22,F23,F31,F32,F33 -- componets 
!!          of the deformation gradient [F]
!!     mu, kappa -- material parameters
!!     ielem -- element id number
!!
!!  OUTPUT
!!     Cchy11, Cchy22, Cchy33,Cchy12,Cchy13,Cchy23  -- componets
!!     of the Cauchy stress tensor
!!
!!  USES
!!     rs, Solve_x
!!
!!
!!***

  IMPLICIT NONE

!--   Variables
  INTEGER :: j,k,l,m,ierr,ielem
  INTEGER :: istep
  REAL*8 :: F11,F12,F13,F21,F22,F23,F31,F32,F33
  REAL*8 :: S11,S22,S33,S12,S13,S23
  REAL*8 :: shear_modulus
  REAL*8 :: delta
  REAL*8 :: Fij(3,3), btens(3,3),princ(3,3)
  REAL*8 :: Cij(3,3), Bulk, sqrt_N, CR, N
  REAL*8 :: e_vec(3,3), Sigma_A(3),stret(1:3)
  REAL*8 :: e_val(3), e_chain, xI3, xmu, xxx, xmu_max, xJay
  REAL*8 :: fv1(3), fv2(3), stretch(3)
  REAL*8 :: fact
  REAL*8 :: mu, kappa
  REAL*8 ::  Cchy11,Cchy22,Cchy33,Cchy12,Cchy13,Cchy23
  REAL*8 :: Pa
  REAL*8 :: sum, sigma(3,3)
  INTEGER :: jd,kd,id


!--(0) initial parameters

!     CR   : (initial shear modulus)/3
!     Bulk : bulk modulus, kappa
!     N    : chain locking stretch ,for example, N=8

  CR   = mu/3.d0
  
  Bulk = 1.d0*kappa
  
  N    = 64.d0
  sqrt_N = SQRT(N)

!--   (1) compute Deformation Gradient Fij 

!     Passed into Subroutine

!                                     
!--   (2) obtains the tensor b, left cauchy-green tensor using equation ???
!
  btens(1,1) = F11**2+F12**2+F13**2
  btens(1,2) = F21*F11+F12*F22+F13*F23
  btens(1,3) = F31*F11+F12*F32+F13*F33
  btens(2,1) = btens(1,2)
  btens(2,2) = F21**2+F22**2+F23**2
  btens(2,3) = F21*F31+F32*F22+F23*F33
  btens(3,1) = btens(1,3)
  btens(3,2) = btens(2,3)
  btens(3,3) = F31**2+F32**2+F33**2

!--  (3) compute eigen values (e_val) and eigen vectors (e_vec) of Cij
!        The rs subroutine seems to be more robust then Jacobi
!
! option 1

! CALL jacobi(btens,stretch,princ)

! option 2

  fv1(1:3) = 0.d0
  fv2(1:3) = 0.d0

  CALL rs(3,3,btens,stretch,1,princ,fv1,fv2,ierr)

  IF(ierr .NE. 0) THEN    
     PRINT *,' error occurs  at element ', ielem
  END IF

!--  (4) calculate stretch
! Not Needed if used subroutine 'jacobi'

  stretch(1:3) = SQRT(stretch(1:3))

!--  (5) compute ramda chain 
 
  e_chain=SQRT(stretch(1)**2+stretch(2)**2+stretch(3)**2) 
  e_chain=1.0d0/SQRT(3.0d0)*(e_chain)

!--  (6) compute I3 and Jay
 
  xJay=  stretch(1)*stretch(2)*stretch(3)
  xI3 = (xJay)**2

!--  (7) compute xmu 
 
  xmu = e_chain/sqrt_N

!--  (8) Solve for x  
 
  CALL Solve_x(xxx,xmu,ielem,stretch)

!--  (9) Compute sigma_A 

  DO k=1,3 
     Sigma_A(k)=CR*sqrt_N*(stretch(k)**2-e_chain**2)/ &
          e_chain*xxx  +  Bulk*LOG(SQRT(xI3))
  END DO

!--  (10a) Compute Cauchy-stress

  DO jd = 1,3     
     DO kd = 1,3
        sum = 0.0d0
        DO id = 1,3
           sum = sum + Sigma_A(id)*princ(jd,id)*princ(kd,id)
        ENDDO
        sigma(jd,kd)=sum
     ENDDO
  ENDDO
         
  Cchy11 = sigma(1,1)
  Cchy22 = sigma(2,2)
  Cchy33 = sigma(3,3)
  Cchy12 = sigma(1,2)
  Cchy23 = sigma(2,3)
  Cchy13 = sigma(1,3)
  RETURN
END SUBROUTINE ARRUDA_BOYCE_CAUCHY

