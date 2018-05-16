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
SUBROUTINE OperatorJ (Lm, Ld, Lb, L2, Lbar,Mm, M2, cm, cb, cd, &
     Dbm, Dbb, Dbd, Dmm,Dmb,Dmd, Ddm, Ddb, Ddd, &
     Am, Ab, Ad, A2, Km, K2, Gm, G2, vm, v2, a_eta, a_zeta, &  
     OpJ, N, Ninv, R, Rinv, Q2, Qm, W)

  IMPLICIT NONE

  integer :: i

  REAL*8 :: cm, cb, cd
  REAL*8 :: a_eta, a_zeta

  REAL*8 :: Gm, G2, Km, K2, vm, v2

  REAL*8, DIMENSION(1:6,1:6) :: Am, Ab, Ad, A2
  REAL*8, DIMENSION(1:6,1:6) :: Lm, Ld, Lb, L2, Lbar
  REAL*8, DIMENSION(1:6,1:6) :: Mm, M2
  REAL*8, DIMENSION(1:6,1:6) :: Dbm, Dbb, Dbd, Dmm,Dmb,Dmd, Ddm, Ddb, Ddd
!                                                       *              *
  
  REAL*8, DIMENSION(1:6,1:6)  :: OpJ , test
  
  REAL*8, DIMENSION(1:6,1:6) :: N, Ninv, R, W, Rinv

  REAL*8, DIMENSION(1:6,1:6) :: Atmp, Lr_LbarInv
  REAL*8, DIMENSION(1:6,1:6) :: Bm, Bb, Bd
  REAL*8, DIMENSION(1:6,1:6) :: Qm, Qb, Qd, Q2

  REAL*8, DIMENSION(1:6,1:6) :: Ident = RESHAPE( &
       (/1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
         0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
         0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, &
         0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, &
         0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, &
         0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0 /),(/6,6/) )

! W


  CALL Constant_W(Gm, G2, K2, Km, a_eta, a_zeta, vm, v2, W)

!  -1
! N  

  N = MATMUL(L2,Ddd - Ident) - MATMUL(W,MATMUL(Lm,Dmd))

  CALL invert2(N, Ninv, 6)

! R

  R = MATMUL(W,MATMUL(Lm,Am)) - MATMUL(L2,A2)

  CALL invert2(R,Rinv,6)

! -- Q --

! Qm
  Atmp = Lm - Lbar
  CALL invert2( Atmp, Lr_LbarInv, 6)

  Qm = MATMUL((Ident - Am),Lr_LbarInv)

!Q2

  Atmp = L2 - Lbar
  CALL invert2( Atmp, Lr_LbarInv, 6)
  
  Q2 = MATMUL(Ident - A2,Lr_LbarInv)

! Qb

  Qb = Q2

! Qd

  Qd = Q2

!                   -1  
! Bm = [ ( Dmd - I)N  (L2 Q2-W Lm Qm )-Qm } A2L2NR

!test  = MATMUL( Dmd * Ninv*(L2*Q2 - W*Lm*Qm) - Qm,TRANSPOSE(A2)*L2*Ninv*R)
!!$  Bb = ( Dbd * Ninv*(L2*Q2 - W*Lm*Qm) - Qb)*TRANSPOSE(A2)*L2*Ninv*R
!!$  Bd = ( (Ddd - Ident) * Ninv*(L2*Q2 - W*Lm*Qm) - Qd)*TRANSPOSE(A2)*L2*Ninv*R



!!$  Bm = MATMUL( MATMUL (Dmd , MATMUL(Ninv,( MATMUL(L2,Q2) - MATMUL(W,MATMUL(Lm,Qm)) ) ) ) - Qm, &
!!$       MATMUL(TRANSPOSE(A2),MATMUL(L2,MATMUL(Ninv,R))) )


  Bm = MATMUL(L2,Q2) - MATMUL(W, MATMUL(Lm,Qm) )

  Bm = MATMUL(Ninv,Bm)

  Bm = MATMUL(Dmd,Bm)

  Bm = Bm - Qm

  Bm = MATMUL(Bm,TRANSPOSE(A2))
  
  Bm = MATMUL(Bm,L2)

  Bm = MATMUL(Bm,Ninv)

  Bm = MATMUL(Bm,R)

!!$  DO i = 1,6 
!!$     PRINT*,'now',Bm(i,1:6)
!!$     PRINT*,'org',test(i,1:6)
!!$  ENDDO
  
  Bb = MATMUL(L2,Q2) - MATMUL(W, MATMUL(Lm,Qm) )

  Bb = MATMUL(Ninv,Bb)

  Bb = MATMUL(Dbd,Bb)

  Bb = Bb - Qb

  Bb = MATMUL(Bb,TRANSPOSE(A2))
  
  Bb = MATMUL(Bb,L2)

  Bb = MATMUL(Bb,Ninv)

  Bb = MATMUL(Bb,R)

!!$  Bb = MATMUL( MATMUL( Dbd , MATMUL(Ninv,( MATMUL(L2,Q2) - MATMUL(W,MATMUL(Lm,Qm)) ) ) ) - Qb , &
!!$       MATMUL(TRANSPOSE(A2),MATMUL(L2,MATMUL(Ninv,R))) )


  Bd = MATMUL(L2,Q2) - MATMUL(W, MATMUL(Lm,Qm) )

  Bd = MATMUL(Ninv,Bd)

  Bd = MATMUL(Ddd - Ident,Bd)

  Bd = Bd - Qd

  Bd = MATMUL(Bd,TRANSPOSE(A2))
  
  Bd = MATMUL(Bd,L2)

  Bd = MATMUL(Bd,Ninv)

  Bd = MATMUL(Bd,R)

!!$  Bd = MATMUL( MATMUL((Ddd - Ident), MATMUL(Ninv, (MATMUL(L2,Q2) - MATMUL(W,MATMUL(Lm,Qm))) ) ) - Qd, &
!!$       MATMUL(TRANSPOSE(A2),MATMUL(L2,MATMUL(Ninv,R))) )

! OpJ = cm*Lm*Bm + cb*L2*Bb + Cd*L2*Bd + L2*( (Ddd - Ident) - Dbd)* Ninv * R

  OpJ = MATMUL( ( Ddd - Ident) - Dbd, Ninv)

  OpJ = MATMUL( L2, OpJ)

  OpJ = MATMUL( OpJ, R)

  
  OpJ = cm*MATMUL(Lm,Bm) + cb*MATMUL(L2,Bb) + cd*MATMUL(L2,Bd) + OpJ



!!$  OpJ = cm*MATMUL(Lm,Bm) + cb*MATMUL(L2,Bb) + cd*MATMUL(L2,Bd) + &
!!$       MATMUL(L2, MATMUL( ( Ddd - Ident) - Dbd, MATMUL(Ninv , R)) )

!!$  DO i = 1, 6
!!$     PRINT*,OpJ(i,1:6)
!!$  ENDDO

END SUBROUTINE OperatorJ

