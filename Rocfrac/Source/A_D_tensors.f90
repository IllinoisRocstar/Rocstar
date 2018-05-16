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
SUBROUTINE A_D_tensors(Lo, Lm, L2, Mm, M2, cm, cb, cd,&
     Dbm,Dbb, Dbd, Dmm,Dmb,Dmd, Ddm, Ddb, Ddd, A2, Am, L_bar)

  IMPLICIT NONE


  REAL*8, DIMENSION(1:6,1:6) :: Lo, Lm, L2
  REAL*8 :: cm, cb, cd


  REAL*8, DIMENSION(1:6,1:6) :: Mo
  REAL*8, DIMENSION(1:6,1:6) :: Mm, M2
  REAL*8, DIMENSION(1:6,1:6) :: P
  REAL*8, DIMENSION(1:6,1:6) :: S, SInv,SI
  REAL*8, DIMENSION(1:6,1:6) :: A,B, Lst
  REAL*8 , DIMENSION(1:6,1:6) :: L, M

  REAL*8, DIMENSION(1:6,1:6) :: Dbm, Dbb, Dbd, Dmm,Dmb,Dmd, Ddm, Ddb, Ddd
  REAL*8, DIMENSION(1:6,1:6) :: A2, Am , C

  REAL*8, DIMENSION(1:6,1:6) :: L_bar, LrLbarInv

  INTEGER :: i, j

  LOGICAL :: debug
  

  REAL*8, DIMENSION(1:6,1:6) :: dident = RESHAPE( &
       (/1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
         0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
         0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, &
         0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, &
         0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, &
         0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0 /),(/6,6/) )

  debug = .false.

!!$  DO i = 1, 6
!!$     PRINT*,'Lo',Lo(i,:)
!!$  ENDDO

  CALL invert2 (Lo, Mo, 6)
  
!  /* Type of inclusion */

! Returns Eshelby tensor S

  CALL Eshelby(0,Mo, Lo, P, S)
!          -1
! Returns S
 
  CALL invert2 (S, Sinv, 6)
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
  

 A = Lst + L_bar

! compute Am

 A2 = Lst + Lm 
 CALL invert2(A2,C,6)
 Am = MATMUL(C,A)
 
! compute A2

 A2 = Lst + L2
 CALL invert2(A2,C,6)
 A2 = MATMUL(C,A)

!!$ PRINT*,'am'
!!$ DO i = 1, 6
!!$    PRINT*,am(i,:)
!!$ ENDDO
!!$ PRINT*,'a2'
!!$ DO i = 1, 6
!!$    PRINT*,a2(i,:)
!!$ ENDDO


! Check transformation factors 
 IF(DEBUG)THEN  
  DO i = 1, 6
     DO j = 1, 6
        C(i,j) = cm*Am(i,j) + (cb + cd)*A2(i,j)
        IF (i.EQ.j .AND. (C(i,j) <= 0.9999 .OR. C(i,j) >=  1.00001))THEN
           PRINT*,'Incorrect Transformation tensors Am, A2'
           STOP
        ENDIF
        IF (i .NE. j .AND. (C(i,j) >= 1.000e-4 .OR. C(i,j) <= -1.000e-4))THEN
           PRINT*,'Incorrect Transformation tensors Am, A2'
           STOP
        ENDIF
     ENDDO
  ENDDO
  ENDIF

  B = Lm - L_bar
  CALL invert2(B,LrLbarInv,6)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! /* Dmm */

  A = dident - Am
  C = dident - cm*TRANSPOSE(Am)
  
  B = MATMUL(A, LrLbarInv)
  A = MATMUL(B, C)
  
  B = MATMUL(A,Lm)
  
  Dmm = B
  IF(DEBUG)THEN  
  PRINT*,'Dmm'
  DO i = 1, 6
     PRINT*,Dmm(:,i)
  ENDDO
  endif

! /* Dmb */

 A = dident - Am
 C = -1.d0*cb*TRANSPOSE(A2)

 B = MATMUL(A, LrLbarInv)
 A = MATMUL(B, C)

 B = MATMUL(A,L2)

 Dmb = B
IF(DEBUG)THEN  
 PRINT*,'Dmb'
  DO i = 1, 6
    PRINT*,Dmb(:,i)
 ENDDO
endif

!  /* Dmd */

 A = dident - Am
 C = -1.d0*cd*TRANSPOSE(A2)

 B = MATMUL(A, LrLbarInv)
 A = MATMUL(B, C)

 B = MATMUL(A,L2)

 Dmd = B

 IF(DEBUG)THEN  
    PRINT*,'Dmd'
    DO i = 1, 6
       PRINT*,Dmd(:,i)
    ENDDO
 ENDIF
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! /* Dbb */

 B = L2 - L_bar
 CALL invert2(B,LrLbarInv,6)

 A = dident - A2
 B = L2 - L_bar
 C = dident - cb*TRANSPOSE(A2) 

 B = MATMUL(A, LrLbarInv)
 A = MATMUL(B, C)

 B = MATMUL(A,L2)

 Dbb = B
 IF(DEBUG)THEN  
    PRINT*,'Dbb'
    DO i = 1, 6
       PRINT*,Dbb(:,i)
    ENDDO
 ENDIF


! /* Dbm */

 A = dident - A2
 C = -1.d0*cm*TRANSPOSE(Am)

 B = MATMUL(A, LrLbarInv)
 A = MATMUL(B, C)

 B = MATMUL(A,Lm)

 Dbm = B
 IF(DEBUG)THEN  
    PRINT*,'Dbm'
    DO i = 1, 6
       PRINT*,Dbm(:,i)
    ENDDO
 ENDIF

! /* Dbd */


 A = dident - A2
 C = -1.d0*cd*TRANSPOSE(A2)

 B = MATMUL(A, LrLbarInv)
 A = MATMUL(B, C)

 B = MATMUL(A,L2)

 Dbd = B
 IF(DEBUG)THEN  
    PRINT*,'Dbd'
    DO i = 1, 6
       PRINT*,Dbd(:,i)
    ENDDO
 ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! /* Ddd */

 A = dident - A2
 C = dident - cd*TRANSPOSE(A2)

 B = MATMUL(A, LrLbarInv)
 A = MATMUL(B, C)

 B = MATMUL(A,L2)

 Ddd = B

 IF(DEBUG)THEN  
    PRINT*,'Ddd'
    DO i = 1, 6
       PRINT*,Ddd(:,i)
    ENDDO
 ENDIF

! /* Ddm */

 A = dident - A2
 C = -1.d0*cm*TRANSPOSE(Am)

 B = MATMUL(A, LrLbarInv)
 A = MATMUL(B, C)

 B = MATMUL(A,Lm)

 Ddm = B 
 IF(DEBUG)THEN  
    PRINT*,'Ddm'
    DO i = 1, 6
       PRINT*,Ddm(:,i)
    ENDDO
 ENDIF

! /* Ddb */

 A = dident - A2
 C = - 1.d0*cb*TRANSPOSE(A2)

 B = MATMUL(A, LrLbarInv)
 A = MATMUL(B, C)

 B = MATMUL(A,L2)

 Ddb = B
 IF(DEBUG)THEN  
    PRINT*,'Ddb'
    DO i = 1, 6
       PRINT*,Ddb(:,i)
    ENDDO
 ENDIF

! /* Check concentration factors */

! /* SUM_s=1^N Drs = I - Ar */


 DO i = 1, 6
    DO j = 1, 6
       A(i,j) = Dmm(i,j) + Dmb(i,j) + Dmd(i,j)
       B(i,j) = Dident(i,j) - Am(i,j)
       C(i,j) = Dbm(i,j) +Dbb(i,j) + Dbd(i,j)
       SI(i,j) = dident(i,j) - A2(i,j)
       S(i,j) = Ddm(i,j) + Ddb(i,j) + Ddd(i,j)

       IF(SQRT((A(i,j) - B(i,j))*(A(i,j) - B(i,j))) >= 1.e-10   .OR. &
            SQRT((C(i,j) - SI(i,j))*(C(i,j) - SI(i,j))) >= 1.e-10 .OR. &
            SQRT((S(i,j) - SI(i,j))*(S(i,j) - SI(i,j))) >= 1.e-10 )THEN
          PRINT*,'Incorrect Concentration tensors SUM_s=1^N Drs = I - Ar'
          STOP
       ENDIF
    ENDDO
 ENDDO

       
! /* SUM_s=1^N Drs*Ms = 0 */

 DO i = 1, 6
    DO j = 1, 6
       A(i,j) = Dmm(i,j)
       B(i,j) = Dmb(i,j)
       C(i,j) = Dmd(i,j)
    ENDDO
 ENDDO
 S = MATMUL(A,Mm)
 SI = MATMUL(B,M2)
 A = MATMUL(C,M2)

 DO i = 1, 6
    DO j = 1, 6
       B(i,j) = S(i,j) + SI(i,j) + A(i,j)
       IF( B(i,j) >= 1.000e-10 .OR. B(i,j) <= -1.00e-10)THEN
          PRINT*,'A) Incorrect Concentration tensors | SUM_s=1^N Drs*Ms = 0'
          STOP
       ENDIF
    ENDDO
 ENDDO

 DO i = 1, 6
    DO j = 1, 6
       A(i,j) = Dbm(i,j)
       B(i,j) = Dbb(i,j)
       C(i,j) = Dbd(i,j)
    ENDDO
 ENDDO
 S = MATMUL(A,Mm)
 SI = MATMUL(B,M2)
 A = MATMUL(C,M2)

 DO i = 1, 6
    DO j = 1, 6
       A(i,j) = S(i,j) + SI(i,j) + A(i,j)
       IF( A(i,j) >= 1.000e-10 .OR. A(i,j) <= -1.00e-10)THEN
          PRINT*,'B) Incorrect Concentration tensors | SUM_s=1^N Drs*Ms = 0'
          STOP
       ENDIF
    ENDDO
 ENDDO


 DO i = 1, 6
    DO j = 1, 6
       A(i,j) = Ddm(i,j)
       B(i,j) = Ddb(i,j)
       C(i,j) = Ddd(i,j)
    ENDDO
 ENDDO
 S = MATMUL(A,Mm)
 SI = MATMUL(B,M2)
 A = MATMUL(C,M2)

 DO i = 1, 6
    DO j = 1, 6
       A(i,j) = S(i,j) + SI(i,j) + A(i,j)
       IF( A(i,j) >= 1.000e-10 .OR. A(i,j) <= -1.00e-10)THEN
          PRINT*,'C) Incorrect Concentration tensors | SUM_s=1^N Drs*Ms = 0'
          STOP
       ENDIF
    ENDDO
 ENDDO

! /* cr*Drs*Ms = cs*Mr*Dsr^T */

 A = cm*Dmb
 B = cb*Dbm

 C = MATMUL(A,M2)
 SI = MATMUL(B,Mm)
 
 DO i = 1, 6
    DO j = 1, 6
       IF(SQRT((C(i,j) - SI(i,j))*(C(i,j) - SI(i,j) )) >= 1.e-10)THEN
          PRINT*,'a) Incorrect Concentration tensors | cr*Drs*Ms = cs*Mr*Dsr'
          stop
       ENDIF
    ENDDO
 ENDDO

 DO i = 1, 6
    DO j = 1,6
       A(i,j) = cb*Dbm(i,j)
       B(i,j) = cm*Dmb(i,j)
    ENDDO
 ENDDO

 C = MATMUL(A,Mm)
 SI = MATMUL(M2,B)
 
 DO i = 1, 6
    DO j = 1,6
       IF(SQRT((C(i,j) - SI(i,j))*(C(i,j) - SI(i,j) )) >= 1.e-10)THEN
          PRINT*,'b) Incorrect Concentration tensors | cr*Drs*Ms = cs*Mr*Dsr'
          stop
       ENDIF
    ENDDO
 ENDDO



 DO i = 1,6
    DO j = 1,6
       A(i,j) = cm*Dmm(i,j) + cb*Dbm(i,j) + cd*Ddm(i,j)
       B(i,j) = cm*Dmb(i,j) + cb*Dbb(i,j) + cd*Ddb(i,j)
       C(i,j) = cm*Dmd(i,j) + cb*Dbd(i,j) + cd*Ddd(i,j)
       IF(A(i,j) >= 1.000e-10 .OR. A(i,j) <= -1.000e-10 .OR. &
            B(i,j) >= 1.000e-10 .OR. B(i,j) <= -1.000e-10 .OR. &
            C(i,j) >= 1.000e-10 .OR. C(i,j) <= -1.000e-10)THEN
          PRINT*,'c) Incorrect Concentration tensors | SUM_s=1^N cs*Dsr = 0'
          STOP
       ENDIF
    ENDDO
 ENDDO


END SUBROUTINE A_D_tensors

!-----------------------------------------------------------------INVERT2
 
!SUBROUTINE invert2( Inv_in, a, nrow)
 
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
!!$  IMPLICIT NONE
!!$ 
!!$!  Argument variables
!!$ 
!!$  INTEGER nrow
!!$  REAL*8, DIMENSION(1:nrow,1:nrow) :: a, Inv_in
!!$  REAL*8 :: det
!!$                                                                                
!!$!  Local variables
!!$                                                                                
!!$  INTEGER ncol, k, j, i
!!$
!!$
!!$  a = Inv_in
!!$                                                                                
!!$!  Invert matrix
!!$                                                                                
!!$  det  = 1.d0
!!$  ncol = nrow
!!$                                                                                
!!$  DO  k = 1, nrow
!!$                                                                                
!!$     IF( a(k,k) .EQ. 0.d0 ) THEN ! Check for singular matrix
!!$        det = 0.d0
!!$        PRINT*,' Matrix is singular for inversion'
!!$        STOP
!!$     END IF
!!$                                                                                
!!$     det = det * a(k, k)
!!$                                                                                
!!$     DO j = 1, ncol
!!$        IF( j .NE. k ) a(k,j) = a(k,j) / a(k,k)
!!$     END DO
!!$                                                                                
!!$     a(k,k) = 1.d0 / a(k,k)
!!$                                                                                
!!$     DO  i = 1, nrow
!!$        IF( i .EQ. k) CYCLE
!!$        DO j = 1, ncol
!!$           IF( j .NE. k) a(i,j) = a(i,j) - a(i,k) * a(k,j)
!!$        END DO
!!$        a(i,k) = -a(i,k) * a(k,k)
!!$     END DO
!!$                                                                                
!!$  END DO
!!$
!!$                                                                                
!!$  RETURN
!!$END SUBROUTINE invert2


SUBROUTINE VecMatVecMul(a,B,c,ndim, scaler)

  IMPLICIT NONE
  
  INTEGER :: ndim
  
  REAL*8, DIMENSION(1:ndim,1:ndim) :: B
  REAL*8, DIMENSION(1:ndim) :: a,c

  REAL*8 :: scaler
  
  INTEGER i, j, k

  REAL*8 :: sum1
  REAL*8, DIMENSION(1:ndim) :: sum2
  
  scaler = 0.d0

  DO i = 1, ndim
     DO j = 1, ndim
        sum1 = 0.d0
        DO k = 1, ndim
           sum1 = sum1 + B(j,k)*c(k) 
        ENDDO
        sum2(j) = sum1
     ENDDO
     scaler = scaler + sum2(i)*a(i)
  ENDDO
  
END SUBROUTINE VecMatVecMul

