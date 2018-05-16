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

subroutine HomogenizedMat(Density, Em, PRm, Ep, PRp, cm, cd_fastest )

  IMPLICIT NONE


  INTEGER :: myid
  REAL*8 :: Density
  REAL*8 :: Em, Gm, PRm, Ep,PRp, Gp
  REAL*8 :: cm, c2

  REAL*8 :: cd_fastest

  REAL*8 :: E_homog
  REAL*8 :: PR_homog

  INTEGER, parameter :: numatv = 2

  REAL*8, DIMENSION(1:6,1:6,1:numatv) :: L_tensor, M_tensor


  REAL*8, DIMENSION(1:6,1:6) :: L_bar, M_bar, Lo

  REAL*8, parameter :: alpha1=1.0, alpha2= 0.0

!-----Create material constant matrices


  L_tensor = 0.d0
  M_tensor = 0.d0
  c2 = 1.d0 - cm
  Gm = Em/(2.d0*(1.d0+PRm))
  Gp = Ep/(2.d0*(1.d0+PRp))

  M_tensor(1,1,1) = 1.d0 / Em
  M_tensor(1,2,1) = - PRm / Em
  M_tensor(1,3,1) = - PRm / Em
  M_tensor(2,1,1) = M_tensor(1,2,1) 
  M_tensor(2,2,1) = 1.d0 / Em
  M_tensor(2,3,1) = - PRm / Em
  M_tensor(3,1,1) = M_tensor(1,3,1)
  M_tensor(3,2,1) = M_tensor(2,3,1)
  M_tensor(3,3,1) = 1.d0 /Em
  M_tensor(4,4,1) = 1.d0/Gm
  M_tensor(5,5,1) = 1.d0/Gm
  M_tensor(6,6,1) = 1.d0/Gm

  M_tensor(1,1,2) = 1.d0 / Ep
  M_tensor(1,2,2) = - PRp / Ep
  M_tensor(1,3,2) = - PRp / Ep
  M_tensor(2,1,2) = M_tensor(1,2,2) 
  M_tensor(2,2,2) = 1.d0 / Ep
  M_tensor(2,3,2) = - PRp / Ep
  M_tensor(3,1,2) = M_tensor(1,3,2)
  M_tensor(3,2,2) = M_tensor(2,3,2)
  M_tensor(3,3,2) = 1.d0 /Ep
  M_tensor(4,4,2) = 1.d0/Gp
  M_tensor(5,5,2) = 1.d0/Gp
  M_tensor(6,6,2) = 1.d0/Gp

  CALL invert2(M_tensor(:,:,1), L_tensor(:,:,1), 6)
  CALL invert2(M_tensor(:,:,2), L_tensor(:,:,2), 6)


  CALL CompositeStiffnes(L_tensor(1:6,1:6,1), L_tensor(1:6,1:6,2), M_tensor(1:6,1:6,1), M_tensor(1:6,1:6,2), cm, c2, &
       alpha1, alpha2, L_bar, M_bar, Lo)


  E_homog = 1.d0/M_bar(1,1)
  PR_homog = -M_bar(1,2)*1./M_bar(1,1) ! nu = -M(1,2)*E

!  IF(myid.EQ.0)THEN
  PRINT*,'Homogenized Stiffness =', E_homog
  PRINT*,'Homogenized Poisson'//"'"//'s ratio',PR_homog
!  ENDIF

  cd_fastest = MAX( cd_fastest, &
             SQRT(E_homog*(1.d0-PR_homog)/Density/(1.d0+PR_homog)/(1.d0-2.d0*PR_homog )) )

  RETURN
END SUBROUTINE HOMOGENIZEDMAT

SUBROUTINE CompositeStiffnes(Lm, L2, Mm, M2, cm, c2, alpha1, alpha2, L_bar, M_bar, Lo)

  IMPLICIT NONE


  REAL*8 :: alpha1, alpha2, cm, c2

  REAL*8, DIMENSION(1:6,1:6) :: Lm, L2, Mm, M2
  REAL*8, DIMENSION(1:6,1:6) :: L_bar, M_bar

  REAL*8, DIMENSION(1:6,1:6) :: dident = RESHAPE( &
       (/1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
         0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
         0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, &
         0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, &
         0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, &
         0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0 /),(/6,6/) )

  REAL*8, DIMENSION(1:6,1:6) :: Lo, Mo 
  REAL*8, DIMENSION(1:6,1:6) :: P
  REAL*8, DIMENSION(1:6,1:6) :: S, Sinv
  REAL*8, DIMENSION(1:6,1:6) :: Am, A2, A, B, C, Lst, A3
  
  INTEGER :: i,j

 

!  Comparison medium,  Lo = (a1*L1 + a2*L2)

  Lo = alpha1*Lm + alpha2*L2

  CALL invert2(Lo, Mo, 6)

  
!  /* Type of inclusion */

! Returns Eshelby tensor S

  CALL Eshelby(0, Mo, Lo, P, S)

  
!          -1
! Returns S

  CALL invert2 (S, Sinv, 6)

!  nas_AB (Lo,SI,A,6,6,6); nas_AB (A,B,Lst,6,6,6);
! 
! I - S

 B = dident - S

! Lo * ( I - S)

 A = MATMUL(Lo,Sinv)
!
! Hill's constrant tensor
!
!  *    -1
! L  = S * Lo * ( I - S)

 Lst = MATMUL(A, B)

!!$ PRINT*,'L*'
!!$
!!$ DO i = 1, 6
!!$    PRINT*,Lst(i,:)
!!$ ENDDO


! /* Compute L(bar) */

 A = Lst + Lm
 B = Lst + L2

 CALL invert2(A,Am,6)
 CALL invert2(B,A2,6)
 C = cm*Am + c2*A2
 CALL invert2(C,A,6)
 L_bar = A - Lst

 A = Lst + L_bar

! compute Am

 A2 = Lst + Lm 
 CALL invert2(A2,C,6)
 Am = MATMUL(C,A)
 
! compute A2

 A2 = Lst + L2
 CALL invert2(A2,C,6)
 A2 = MATMUL(C,A)

! Check concentration factors

 DO i = 1, 6
    DO j = 1, 6
       c(i,j) = cm*Am(i,j) + c2*A2(i,j)
       IF ( i .EQ. j .AND. (C(i,j) <= 0.9999 .OR. C(i,j) >=  1.0001))THEN
          PRINT*,'a) Incorrect Concentration tensors Am, A2'
          STOP
       ENDIF
       IF ( i.NE.j .AND.( C(i,j) >= 1.e-10 .OR. C(i,j) <= -1.e-10))THEN
          PRINT*,'b) Incorrect Concentration tensors Am, A2'
          STOP
       ENDIF
    ENDDO
 ENDDO

 CALL invert2(L_Bar, M_bar, 6)

END SUBROUTINE CompositeStiffnes

SUBROUTINE Eshelby(Shape, Mo, Lo, P, S)

  INTEGER :: Shape
  REAL*8, DIMENSION(1:6,1:6) :: Mo, Lo
  REAL*8, DIMENSION(1:6,1:6) :: P

  REAL*8, DIMENSION(1:6,1:6) :: S
  

  REAL*8 :: E, G, nu, K, si, th


  S = 0.d0
  P = 0.d0

  ! Type of inclusion

  SELECT CASE (Shape)

  CASE (0) ! /* Spherical inclusion in Isotropic medium */
     E = 1.d0/Mo(1,1)
     G = 1.d0/Mo(4,4)
     nu = -Mo(1,2)*E
     K = E/(3.d0*(1.d0 - 2.d0*nu))
     si = (1.d0 + nu)/(3.d0*(1.d0 - nu)) 
     th = 2.d0*(4.d0 - 5.d0*nu)/(15.d0*(1.d0 - nu))
     
     P(1,1) = 1.d0/3.d0*(si/3.d0/K + th/G)
     P(1,2) = 1.d0/3.d0*(si/3.d0/K - th/2.d0/G)
     P(1,3) = 1.d0/3.d0*(si/3.d0/K - th/2.d0/G)
     P(2,1) = 1.d0/3.d0*(si/3.d0/K - th/2.d0/G)
     P(2,2) = 1.d0/3.d0*(si/3.d0/K + th/G)
     P(2,3) = 1.d0/3.d0*(si/3.d0/K - th/2.d0/G)
     P(3,1) = 1.d0/3.d0*(si/3.d0/K - th/2.d0/G)
     P(3,2) = 1.d0/3.d0*(si/3.d0/K - th/2.d0/G)
     P(3,3) = 1.d0/3.d0*(si/3.d0/K + th/G)
     P(4,4) = th/G
     P(5,5) = th/G
     P(6,6) = th/G
     S = MATMUL(Lo,P)

  CASE DEFAULT
     PRINT*, "Non existing shape of inclusion"
     STOP
  END SELECT

END SUBROUTINE Eshelby

SUBROUTINE invert2( Inv_in, a, nrow)
                                                                                                          
! *--------------------------------------------------------------------*
! |                                                                    |
! |   Inverts a matrix by gauss jordan elimination and computes the    |
! |   the determinant. The matrix may be symmetric or nonsymmetric.    |
! |                                                                    |
! |    <a>    : a square, double precision matrix                      |
! |    <nrow> : dimensioned number of rows for 'a'. also the actual    |
! |             sizes for inversion computations                       |
! |    <det>  : determinant of the matrix ( zero if matrix is singular)|
! |                                                                    |
! |     NOTE  : the matrix is overwritten by the inverse.              |
! |                                                                    |
! *--------------------------------------------------------------------*
  IMPLICIT NONE
                                                                                                          
!  Argument variables
                                                                                                          
  INTEGER nrow
  REAL*8, DIMENSION(1:nrow,1:nrow) :: a, Inv_in
  REAL*8 :: det
                                                                                 
!  Local variables
                                                                                 
  INTEGER ncol, k, j, i
                                                                                                          
  a = Inv_in
                                                                                 
!  Invert matrix
                                                                                 
  det  = 1.d0
  ncol = nrow
                                                                                 
  DO  k = 1, nrow
                                                                                 
     IF( a(k,k) .EQ. 0.d0 ) THEN ! Check for singular matrix
        det = 0.d0
        PRINT*,' Matrix is singular for inversion'
        STOP
     END IF
                                                                                 
     det = det * a(k, k)
                                                                                 
     DO j = 1, ncol
        IF( j .NE. k ) a(k,j) = a(k,j) / a(k,k)
     END DO
                                                                                 
     a(k,k) = 1.d0 / a(k,k)
                                                                                 
     DO  i = 1, nrow
        IF( i .EQ. k) CYCLE
        DO j = 1, ncol
           IF( j .NE. k) a(i,j) = a(i,j) - a(i,k) * a(k,j)
        END DO
        a(i,k) = -a(i,k) * a(k,k)
     END DO
                                                                                 
  END DO
                                                                                                          
                                                                                 
  RETURN
END SUBROUTINE invert2

