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
SUBROUTINE ARRUDA_BOYCE(Cij, &
     S11,S22,S33,S12,S23,S13,ielem,mu,kappa)

!!****f* Rocfrac/Rocfrac/Source/arruda_boyce.f90
!!
!!  NAME
!!     ARRUDA_BOYCE
!!
!!  FUNCTION
!!     Arruda-Boyce constitutive model, returns 
!!     2nd Piola-Kircheoff Stresses
!!
!!  INPUTS
!!     Cij = right Cauchy-Green deformation tensor
!!     mu, kappa -- material parameters
!!     ielem -- element id number
!!
!!  OUTPUT
!!     S11,S22,S33,S12,S23,S13  -- componets of 
!!           the 2nd Piola-Kircheoff Stresses
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
  REAL*8 :: S11,S22,S33,S12,S13,S23
  REAL*8 :: shear_modulus
  REAL*8 :: delta
  REAL*8 :: Cij(3,3), Bulk, sqrt_N, CR, N
  REAL*8 :: e_vec(3,3), Sigma_A(3)
  REAL*8 :: e_val(3), e_chain, xI3, xmu, xxx, xmu_max, xJay
  REAL*8 :: fv1(3), fv2(3), stretch(3)
  REAL*8 :: fact
  REAL*8 :: mu, kappa

!--   (0) initial parameters

!        CR   : (initial shear modulus)/3
!        Bulk : bulk modulus, kappa
!        N    : chain locking stretch ,for example, N=8

  CR   = mu/3.d0
  
  Bulk = 1.d0*kappa
  
  N    = 4.d0
  
  sqrt_N = SQRT(N)

!                                     
!     obtains the tensor b, left cauchy-green tensor using equation ???
!
!$$$         btens(1,1) = F11**2+F12**2+F13**2
!$$$         btens(1,2) = F21*F11+F12*F22+F13*F23
!$$$         btens(1,3) = F31*F11+F12*F32+F13*F33
!$$$         btens(2,1) = btens(1,2)
!$$$         btens(2,2) = F21**2+F22**2+F23**2
!$$$         btens(2,3) = F21*F31+F32*F22+F23*F33
!$$$         btens(3,1) = btens(1,3)
!$$$         btens(3,2) = btens(2,3)
!$$$         btens(3,3) = F31**2+F32**2+F33**2
!$$$
!$$$         CALL jacobi(btens,stret,princ)


         
!--  (3) compute eigen values (e_val) and eigen vectors (e_vec) of Cij

  DO j=1,3
     stretch(j)=0.0d0
     e_val(j)=0.0d0
     fv1(j)=0.0d0
     fv2(j)=0.0d0
     DO k=1,3
        e_vec(j,k)=0.0d0
     END DO
  END DO
!
  CALL rs(3,3,Cij,e_val,1,e_vec,fv1,fv2,ierr)
  IF(ierr .NE. 0) THEN    
     PRINT *,' error occurs  at element ', ielem
  END IF

!--  (4) calculate stretch

  DO j=1,3
     stretch(j)=SQRT(e_val(j))
  END DO

!--  (5) compute ramda chain 
 
  e_chain=SQRT(stretch(1)**2+stretch(2)**2+stretch(3)**2) 
  e_chain=1.0d0/SQRT(3.0d0)*(e_chain)

!$$$         e_vec(:,:) = princ(:,:)
!$$$
!$$$         stretch(1:3) = stret(1:3)


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

!--  (10) Compute 2nd Piola-Kircheoff Stresses
      
  S11 = 0.0d0
  S22 = 0.0d0
  S33 = 0.0d0
  S12 = 0.0d0
  S23 = 0.0d0
  S13 = 0.0d0

  DO k = 1, 3           
     fact = xJay * Sigma_A(k) / (stretch(k)**2)
     
     S11 = S11 + e_vec(1,k) * e_vec(1,k) * fact
     S22 = S22 + e_vec(2,k) * e_vec(2,k) * fact
     S33 = S33 + e_vec(3,k) * e_vec(3,k) * fact
     S12 = S12 + e_vec(1,k) * e_vec(2,k) * fact
     S23 = S23 + e_vec(2,k) * e_vec(3,k) * fact
     S13 = S13 + e_vec(1,k) * e_vec(3,k) * fact
  END DO
  
  RETURN
END SUBROUTINE ARRUDA_BOYCE


SUBROUTINE Solve_x(x,xmu,i,stretch)

!!****f* Rocfrac/Rocfrac/Source/arruda_boyce.f90/Solve_x
!!
!!  NAME
!!     Solve_x
!!
!!  FUNCTION
!!     Solve for the x
!!
!!  INPUTS
!!     xmu -- mu
!!       i -- element id
!! stretch -- principle stretch
!!
!!  OUTPUT
!!       x  -- x
!!
!!***
  
  IMPLICIT NONE
  
  
  REAL*8 x, xmu, F, Fp, dx, x0
  REAL*8 sinh, cosh, coth 
  REAL*8 stretch(3)
  INTEGER iter, i, k
  
  x0 = 0.1d0
  x  = x0
  iter = 1
  
10 sinh = (EXP(x) - EXP(-x))*0.5d0
  cosh = (EXP(x) + EXP(-x))*0.5d0
  coth = cosh/sinh 
  F = coth -1.d0/x - xmu
  IF(ABS(F) .LE. 1.0d-10) go to 20
  Fp = -1.d0/(sinh**2) + 1.d0/(x**2)
  dx = -F/Fp 
  x = x + dx 
  
  IF(iter .GT. 10000) THEN
     PRINT *,' iteration exceeds 10000'
     PRINT *,' program stop ! '
     PRINT *,' at element ',i 
     WRITE(*,*) (stretch(k),k=1,3) 
     STOP
  END IF
  
  iter=iter+1
  go to 10    
  
20 CONTINUE
  RETURN
END SUBROUTINE Solve_x

SUBROUTINE rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
        
  IMPLICIT NONE

  INTEGER n,nm,ierr,matz
  REAL*8  a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)

!!****f* Rocfrac/Rocfrac/Source/arruda_boyce.f90/rs
!!
!!  NAME
!!     rs
!!
!!  FUNCTION
!!     calls the recommended SEQUENCE of
!!     subroutines from the eigensystem SUBROUTINE package (eispack)
!!     to find the eigenvalues and eigenvectors (IF desired)
!!     of a REAL symmetric matrix.
!!
!!
!!  INPUTS
!!
!!     nm  must be set to the row DIMENSION of the two-dimensional
!!      array parameters as declared in the calling PROGRAM
!!      DIMENSION statement.
!!
!!     n  is the order of the matrix  a.
!!
!!     a  CONTAINS the REAL symmetric matrix.
!!
!!  matz  is an INTEGER variable set equal to zero IF
!!        ONLY eigenvalues are desired.  otherwise it is set to
!!        any non-zero INTEGER for both eigenvalues and eigenvectors.
!!
!!  OUTPUT
!!
!!     w  CONTAINS the eigenvalues in ascending order.
!!
!!     z  CONTAINS the eigenvectors IF matz is not zero.
!!
!!     ierr  is an INTEGER output variable set equal to an error
!!           completion code described in the documentation for tqlrat
!!           and tql2.  the normal completion code is zero.
!!
!!     fv1  and  fv2  are temporary storage arrays.
!!
!!  NOTES
!!
!!     questions and comments should be directed to burton s. garbow,
!!     mathematics and computer science div, argonne national laboratory
!!
!!     this version dated august 1983.
!!
!!***

  IF (n .LE. nm) go to 10
  ierr = 10 * n
  go to 50
!
10 IF (matz .NE. 0) go to 20
!     .......... find eigenvalues ONLY ..........
  CALL  tred1(nm,n,a,w,fv1,fv2)
!  tqlrat encounters catastrophic underflow on the Vax
!     CALL  tqlrat(n,w,fv2,ierr)
CALL  tql1(n,w,fv1,ierr)
go to 50
!     .......... find both eigenvalues and eigenvectors ..........
20 CALL  tred2(nm,n,a,w,fv1,z)
CALL  tql2(nm,n,w,fv1,z,ierr)
50 RETURN
END SUBROUTINE rs

REAL*8 FUNCTION pythag(a,b)

!!****f* Rocfrac/Rocfrac/Source/arruda_boyce.f90/pythag
!!
!!  NAME
!!     pythag
!!
!!  FUNCTION
!!     finds SQRT(a**2+b**2) without overflow or destructive underflow
!!*****

  IMPLICIT NONE
  REAL*8 a,b

  REAL*8 p,r,s,t,u
!     DOUBLE PRECISION p,r,s,t,u
  p = dmax1(ABS(a),ABS(b))
  IF (p .EQ. 0.0d0) go to 20
  r = (dmin1(ABS(a),ABS(b))/p)**2
10 CONTINUE
  t = 4.0d0 + r
  IF (t .EQ. 4.0d0) go to 20
  s = r/t
  u = 1.0d0 + 2.0d0*s
  p = u*p
  r = (s/u)**2 * r
  go to 10
20 pythag = p
  RETURN
END FUNCTION pythag

SUBROUTINE tql1(n,d,e,ierr)

!!****f* Rocfrac/Rocfrac/Source/arruda_boyce.f90/tql1
!!
!!  NAME
!!     tql1
!!
!!  FUNCTION
!!     finds the eigenvalues of a symmetric
!!     tridiagonal matrix by the ql method.
!!
!!   
!!  INPUTS
!!
!!        n is the order of the matrix.
!!
!!        d CONTAINS the diagonal elements of the input matrix.
!!
!!        e CONTAINS the subdiagonal elements of the input matrix
!!          in its last n-1 positions.  e(1) is arbitrary.
!!
!!  OUTPUT
!!
!!        d CONTAINS the eigenvalues in ascending order.  IF an
!!          error EXIT is made, the eigenvalues are correct and
!!          ordered for indices 1,2,...ierr-1, but may not be
!!          the smallest eigenvalues.
!!
!!        e has been destroyed.
!!
!!        ierr is set to
!!          zero       for normal RETURN,
!!          j          IF the j-th eigenvalue has not been
!!                     determined after 30 iterations.
!!   USES
!!    
!!     pythag for SQRT(a*a + b*b) .
!!
!!  NOTES
!!
!!     this SUBROUTINE is a translation of the algol PROCEDURE tql1,
!!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!!     wilkinson.
!!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!!
!!
!!     questions and comments should be directed to burton s. garbow,
!!     mathematics and computer science div, argonne national laboratory
!!
!!     this version dated august 1983.
!!
!!***** 
  IMPLICIT NONE

  INTEGER i,j,l,m,n,ii,l1,l2,mml,ierr
  REAL*8 d(n),e(n)
  REAL*8 c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag



!     ------------------------------------------------------------------
!
  ierr = 0
  IF (n .EQ. 1) go to 1001

  DO 100 i = 2, n
100  e(i-1) = e(i)
     
     f = 0.0d0
     tst1 = 0.0d0
     e(n) = 0.0d0

     DO 290 l = 1, n
        j = 0
        h = ABS(d(l)) + ABS(e(l))
        IF (tst1 .LT. h) tst1 = h
!     .......... look for small sub-diagonal element ..........
        DO 110 m = l, n
           tst2 = tst1 + ABS(e(m))
           IF (tst2 .EQ. tst1) go to 120
!     .......... e(n) is always zero, so there is no EXIT
!                through the bottom of the loop ..........
110        CONTINUE
!
120        IF (m .EQ. l) go to 210
130        IF (j .EQ. 30) go to 1000
           j = j + 1
!     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + SIGN(r,p))
         d(l1) = e(l) * (p + SIGN(r,p))
         dl1 = d(l1)
         h = g - d(l)
         IF (l2 .GT. n) go to 145
!
         DO 140 i = l2, n
  140    d(i) = d(i) - h
!
  145    f = f + h
!     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
!     .......... for i=m-1 step -1 until l DO -- ..........
         DO 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
  200    CONTINUE
!
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + ABS(e(l))
         IF (tst2 .GT. tst1) go to 130
  210    p = d(l) + f
!     .......... order eigenvalues ..........
         IF (l .EQ. 1) go to 250
!     .......... for i=l step -1 until 2 DO -- ..........
         DO 230 ii = 2, l
            i = l + 2 - ii
            IF (p .GE. d(i-1)) go to 270
            d(i) = d(i-1)
  230    CONTINUE
!
  250    i = 1
  270    d(i) = p
  290 CONTINUE
!
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 RETURN
      END
      SUBROUTINE tql2(nm,n,d,e,z,ierr)

!!****f* Rocfrac/Rocfrac/Source/arruda_boyce.f90/tql2
!!
!!  NAME
!!     tql2
!!
!!  FUNCTION
!!  
!!     Finds the eigenvalues and eigenvectors
!!     of a symmetric tridiagonal matrix by the ql method.
!!     the eigenvectors of a full symmetric matrix can also
!!     be found IF  tred2  has been used to reduce this
!!     full matrix to tridiagonal form.

!!  NOTES
!!     this SUBROUTINE is a translation of the algol PROCEDURE tql2,
!!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!!     wilkinson.
!!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!!
!!
!!     questions and comments should be directed to burton s. garbow,
!!     mathematics and computer science div, argonne national laboratory
!!
!!     this version dated august 1983.
!!
!!
!!     INPUTS
!!
!!        nm must be set to the row DIMENSION of two-dimensional
!!          array parameters as declared in the calling PROGRAM
!!          DIMENSION statement.
!!
!!        n is the order of the matrix.
!!
!!        d CONTAINS the diagonal elements of the input matrix.
!!
!!        e CONTAINS the subdiagonal elements of the input matrix
!!          in its last n-1 positions.  e(1) is arbitrary.
!!
!!        z CONTAINS the transformation matrix produced in the
!!          reduction by  tred2, IF performed.  IF the eigenvectors
!!          of the tridiagonal matrix are desired, z must contain
!!          the identity matrix.
!!
!!      OUTPUT
!!
!!        d CONTAINS the eigenvalues in ascending order.  IF an
!!          error EXIT is made, the eigenvalues are correct but
!!          unordered for indices 1,2,...,ierr-1.
!!
!!        e has been destroyed.
!!
!!        z CONTAINS orthonormal eigenvectors of the symmetric
!!          tridiagonal (or full) matrix.  IF an error EXIT is made,
!!          z CONTAINS the eigenvectors associated WITH the stored
!!          eigenvalues.
!!
!!        ierr is set to
!!          zero       for normal RETURN,
!!          j          IF the j-th eigenvalue has not been
!!                     determined after 30 iterations.
!!
!!    USES
!!     pythag for  SQRT(a*a + b*b) .
!!
!!***

      IMPLICIT NONE
!
      INTEGER i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      REAL*8 d(n),e(n),z(nm,n)
      REAL*8 c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
!     DOUBLE PRECISION d(n),e(n),z(nm,n)
!     DOUBLE PRECISION c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
!
      ierr = 0
      IF (n .EQ. 1) go to 1001

      DO 100 i = 2, n
  100 e(i-1) = e(i)

      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0

      DO 240 l = 1, n
         j = 0
         h = ABS(d(l)) + ABS(e(l))
         IF (tst1 .LT. h) tst1 = h
!     .......... look for small sub-diagonal element ..........
         DO 110 m = l, n
            tst2 = tst1 + ABS(e(m))
            IF (tst2 .EQ. tst1) go to 120
!     .......... e(n) is always zero, so there is no EXIT
!                through the bottom of the loop ..........
  110    CONTINUE
!
  120    IF (m .EQ. l) go to 220
  130    IF (j .EQ. 30) go to 1000
         j = j + 1
!     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + SIGN(r,p))
         d(l1) = e(l) * (p + SIGN(r,p))
         dl1 = d(l1)
         h = g - d(l)
         IF (l2 .GT. n) go to 145
!
         DO 140 i = l2, n
  140    d(i) = d(i) - h
!
  145    f = f + h
!     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
!     .......... for i=m-1 step -1 until l DO -- ..........
         DO 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
!     .......... form vector ..........
            DO 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       CONTINUE
!
  200    CONTINUE
!
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + ABS(e(l))
         IF (tst2 .GT. tst1) go to 130
  220    d(l) = d(l) + f
  240 CONTINUE
!     .......... order eigenvalues and eigenvectors ..........
      DO 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
!
         DO 260 j = ii, n
            IF (d(j) .GE. p) go to 260
            k = j
            p = d(j)
  260    CONTINUE
!
         IF (k .EQ. i) go to 300
         d(k) = d(i)
         d(i) = p
!
         DO 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    CONTINUE
!
  300 CONTINUE
!
      go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 RETURN
      END
      SUBROUTINE tred1(nm,n,a,d,e,e2)

      IMPLICIT NONE
!
      INTEGER i,j,k,l,n,ii,nm,jp1
      REAL*8 a(nm,n),d(n),e(n),e2(n)
      REAL*8 f,g,h,scale

!!****f* Rocfrac/Rocfrac/Source/arruda_boyce.f90/tred1
!!
!!  NAME
!!     tred1
!!
!!  FUNCTION 
!!     Reduces a REAL symmetric matrix
!!     to a symmetric tridiagonal matrix using
!!     orthogonal similarity transformations.
!!
!!     INPUTS
!!
!!        nm must be set to the row DIMENSION of two-dimensional
!!          array parameters as declared in the calling PROGRAM
!!          DIMENSION statement.
!!
!!        n is the order of the matrix.
!!
!!        a CONTAINS the REAL symmetric input matrix.  ONLY the
!!          lower triangle of the matrix need be supplied.
!!
!!     OUTPUT
!!
!!        a CONTAINS information about the orthogonal trans-
!!          formations used in the reduction in its strict lower
!!          triangle.  the full upper triangle of a is unaltered.
!!
!!        d CONTAINS the diagonal elements of the tridiagonal matrix.
!!
!!        e CONTAINS the subdiagonal elements of the tridiagonal
!!          matrix in its last n-1 positions.  e(1) is set to zero.
!!
!!        e2 CONTAINS the squares of the corresponding elements of e.
!!          e2 may coincide WITH e IF the squares are not needed.
!!
!!  NOTES
!!
!!     this SUBROUTINE is a translation of the algol PROCEDURE tred1,
!!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!!
!!     questions and comments should be directed to burton s. garbow,
!!     mathematics and computer science div, argonne national laboratory
!!
!!     this version dated august 1983.
!!
!!****

      DO 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 CONTINUE
!     .......... for i=n step -1 until 1 DO -- ..........
      DO 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         IF (l .LT. 1) go to 130
!     .......... scale row (algol tol THEN not needed) ..........
         DO 120 k = 1, l
  120    scale = scale + ABS(d(k))

         IF (scale .NE. 0.0d0) go to 140

         DO 125 j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = 0.0d0
  125    CONTINUE

  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 300

  140    DO 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    CONTINUE

         e2(i) = scale * scale * h
         f = d(l)
         g = -SIGN(SQRT(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         IF (l .EQ. 1) go to 285
!     .......... form a*u ..........
         DO 170 j = 1, l
  170    e(j) = 0.0d0

         DO 240 j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            IF (l .LT. jp1) go to 220

            DO 200 k = jp1, l
               g = g + a(k,j) * d(k)
               e(k) = e(k) + a(k,j) * f
  200       CONTINUE

  220       e(j) = g
  240    CONTINUE
!     .......... form p ..........
         f = 0.0d0

         DO 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    CONTINUE

         h = f / (h + h)
!     .......... form q ..........
         DO 250 j = 1, l
  250    e(j) = e(j) - h * d(j)
!     .......... form reduced a ..........
         DO 280 j = 1, l
            f = d(j)
            g = e(j)

            DO 260 k = j, l
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)

  280    CONTINUE

  285    DO 290 j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
  290    CONTINUE

  300 CONTINUE

      RETURN
      END

      SUBROUTINE tred2(nm,n,a,d,e,z)

!!****f* Rocfrac/Rocfrac/Source/arruda_boyce.f90/tred2
!!
!!  NAME
!!     tred2
!!
!!  FUNCTION
!!     reduces a REAL symmetric matrix to a
!!     symmetric tridiagonal matrix using and accumulating
!!     orthogonal similarity transformations.
!!
!!      INPUTS
!!
!!        nm must be set to the row DIMENSION of two-dimensional
!!          array parameters as declared in the calling PROGRAM
!!          DIMENSION statement.
!!
!!        n is the order of the matrix.
!!
!!        a CONTAINS the REAL symmetric input matrix.  ONLY the
!!          lower triangle of the matrix need be supplied.
!!
!!     OUTPUT
!!
!!        d CONTAINS the diagonal elements of the tridiagonal matrix.
!!
!!        e CONTAINS the subdiagonal elements of the tridiagonal
!!          matrix in its last n-1 positions.  e(1) is set to zero.
!!
!!        z CONTAINS the orthogonal transformation matrix
!!          produced in the reduction.
!!
!!        a and z may coincide.  IF distinct, a is unaltered.
!!
!!   NOTES
!!     this SUBROUTINE is a translation of the algol PROCEDURE tred2,
!!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!!     questions and comments should be directed to burton s. garbow,
!!     mathematics and computer science div, argonne national laboratory
!!
!!     this version dated august 1983.
!!****


      IMPLICIT NONE

      INTEGER i,j,k,l,n,ii,nm,jp1
      REAL*8 a(nm,n),d(n),e(n),z(nm,n)
      REAL*8 f,g,h,hh,scale

      DO 100 i = 1, n

         DO 80 j = i, n
   80    z(j,i) = a(j,i)

         d(i) = a(n,i)
  100 CONTINUE

      IF (n .EQ. 1) go to 510
!     .......... for i=n step -1 until 2 DO -- ..........
      DO 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         IF (l .LT. 2) go to 130
!     .......... scale row (algol tol THEN not needed) ..........
         DO 120 k = 1, l
  120    scale = scale + ABS(d(k))

         IF (scale .NE. 0.0d0) go to 140
  130    e(i) = d(l)

         DO 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    CONTINUE

         go to 290

  140    DO 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    CONTINUE

         f = d(l)
         g = -SIGN(SQRT(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
!     .......... form a*u ..........
         DO 170 j = 1, l
  170    e(j) = 0.0d0

         DO 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            IF (l .LT. jp1) go to 220

            DO 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       CONTINUE

  220       e(j) = g
  240    CONTINUE
!     .......... form p ..........
         f = 0.0d0

         DO 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    CONTINUE

         hh = f / (h + h)
!     .......... form q ..........
         DO 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
!     .......... form reduced a ..........
         DO 280 j = 1, l
            f = d(j)
            g = e(j)

            DO 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)

            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    CONTINUE

  290    d(i) = h
  300 CONTINUE
!     .......... accumulation of transformation matrices ..........
      DO 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         IF (h .EQ. 0.0d0) go to 380

         DO 330 k = 1, l
  330    d(k) = z(k,i) / h

         DO 360 j = 1, l
            g = 0.0d0

            DO 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)

            DO 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    CONTINUE

  380    DO 400 k = 1, l
  400    z(k,i) = 0.0d0

  500 CONTINUE

  510 DO 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 CONTINUE

      z(n,n) = 1.0d0
      e(1) = 0.0d0
      RETURN
      END

