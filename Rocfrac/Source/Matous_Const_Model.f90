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
SUBROUTINE Matous_Const_Model(cd, cb, &
     p1, p2, Yin, SoftParam, Km, K2, Gm, G2, vm, nu2, a_eta, a_zeta,&
     cm, c2, Lo, Lm, L2, Mm, M2, L_bar, strain, DelStrain,&
     s11, s22, s33, s12, s13, s23,ielem)


  USE Precision

  IMPLICIT NONE

! input

  INTEGER :: ielem
  REAL(KIND=wp) :: cd, cb, SoftParam

  REAL(KIND=wp), DIMENSION(1:6,1:6) :: Lo, L_bar
  REAL(KIND=wp) :: Km, K2, Gm, G2, vm, nu2, a_eta, a_zeta, cm

! added

  REAL(KIND=wp) :: p1, p2, Yin, c2, Y

  REAL(KIND=wp) :: Gy, tot

  REAL(KIND=wp), DIMENSION(1:6,1:6) ::  ScriptJ, tmpM
  INTEGER :: kk,ii,jj
  REAL(KIND=wp), DIMENSION(1:6) :: strain, ScriptJstrain, totv

  REAL(KIND=wp), DIMENSION(1:6,1:6) :: Lm, L2, Mm, M2

  REAL(KIND=wp), DIMENSION(1:6,1:6) :: Dbm,Dbb, Dbd, Dmm,Dmb,Dmd, Ddm, Ddb, Ddd

  REAL(KIND=wp), DIMENSION(1:6,1:6) :: A2, Am

  REAL(KIND=wp), DIMENSION(1:6,1:6) ::  LLbar

  REAL(KIND=wp) :: Skk

  REAL(KIND=wp), DIMENSION(1:6,1:6) :: ScriptAm, ScriptAb, ScriptAd 

  REAL(KIND=wp), DIMENSION(1:6,1:6) :: ScriptK

  REAL(KIND=wp), DIMENSION(1:6,1:6) :: Rinv, Ninv, W, QM, Q2, R, N

  REAL(KIND=wp), DIMENSION(1:6) :: Pt 

  REAL(KIND=wp) :: H, eTKe, Dcd

  REAL(KIND=wp), PARAMETER :: tol = 1.e-6

  REAL(KIND=wp), DIMENSION(1:6,1:6) :: CheckStress1, CheckStress2, test

  REAL(KIND=wp), DIMENSION(1:6) :: mu_d, strain_m, strain_d, strain_b
  REAL(KIND=wp), DIMENSION(1:6) :: stress_d, stress_b, stress_m

  REAL(KIND=wp), DIMENSION(1:6,1:6) :: dident = RESHAPE( &
       (/1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
       0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
       0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
       0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &
       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp /),(/6,6/) )


  REAL(KIND=wp) :: DelScriptK

  REAL(KIND=wp), DIMENSION(1:6) :: DelStrain

  REAL(KIND=wp) :: sum1, sum2, scalar1, scalar2, cd_gauss_old, nor, nor1

  integer :: icnt
! local
  REAL(KIND=wp) :: Dcd_Old

  REAL(KIND=wp) :: s11, s22, s33, s12, s13, s23

  REAL(KIND=wp) :: cd_n, gytmp

  INTEGER,PARAMETER :: MonNd = 175


  Dcd_Old = 0.0_wp

  cb = c2 - cd

  cd_n = cd
     
!  PRINT*,'&&&', ielem,cd, Dcd, cd_n
  CALL A_D_tensors(Lo, Lm, L2, Mm, M2, cm, cb, cd, &
       Dbm,Dbb, Dbd, Dmm,Dmb,Dmd, Ddm, Ddb, Ddd, A2, Am, L_bar)
!  PRINT*,'ljlk'

! present context
!     Lb = L2 ;  Ad = A2
!     Ld = L2 ; Ab = A2
                      !Ld, Lb


  CALL OperatorJ (Lm, L2, L2, L2, L_bar, Mm, M2, cm, cb, cd, &
       Dbm, Dbb, Dbd, Dmm,Dmb,Dmd, Ddm, Ddb, Ddd, &
       Am, A2, A2, A2, Km, K2, Gm, G2, vm, nu2, a_eta, a_zeta, &  
       ScriptJ, N, Ninv, R, Rinv, Q2, Qm, W)

     

!     PRINT*,'Am,A2',Am(1,1),A2(1,1)
!     PRINT*,'Lbar', L_bar(1,1)

! ScriptK = 2*ScriptJ*Rinv*(L2*Q2-W*Lm*Qm)*TRANSPOSE(A2) * L2 * Ninv*R
     
  ScriptK = MATMUL(L2,Q2) - MATMUL(W,MATMUL(Lm,Qm))
  
  ScriptK = MATMUL(Rinv,ScriptK)
  
  ScriptK = 2.d0*MATMUL(ScriptJ,ScriptK)
  
  ScriptK = MATMUL(ScriptK,TRANSPOSE(A2) )
  
  ScriptK = MATMUL(ScriptK,L2 )
  
  ScriptK = MATMUL(ScriptK,Ninv )
  
  ScriptK = MATMUL(ScriptK,R ) 
  
  ScriptAm = Am + MATMUL(Dmd,MATMUL(Ninv,R))
  ScriptAb = A2 + MATMUL(Dbd,MATMUL(Ninv,R))
  ScriptAd = A2 + MATMUL(Ddd-dident, MATMUL(Ninv,R))
  
  LLbar = cm*MATMUL(Lm,ScriptAm) + cb*MATMUL(L2,ScriptAb) + cd*MATMUL(L2,ScriptAd)

  S11 = LLbar(1,1)*Delstrain(1) + LLbar(1,2)*Delstrain(2) + LLbar(1,3)*Delstrain(3) 
  S22 = LLbar(2,1)*Delstrain(1) + LLbar(2,2)*Delstrain(2) + LLbar(2,3)*Delstrain(3) 
  S33 = LLbar(3,1)*Delstrain(1) + LLbar(3,2)*Delstrain(2) + LLbar(3,3)*Delstrain(3)

  Skk = ( S11 + S22 + S33 )/3.d0

!              _ T       _
!  Y = - 0.5 * E   * J * E               ... (24)
!

  CALL VecMatVecMul( strain, ScriptJ, strain, 6, Y)
  
  CALL VecMatVecMul( strain, ScriptJ, Delstrain, 6, scalar1)

  Y = -0.5d0 * Y

! The function G(y) characterizes the damage process given by the Weibull distribution
  IF(Y.LT.Yin)THEN
     Gy = 0.d0
  ELSE
     Gy = c2  - c2 * EXP( - ( (Y - Yin)/(p1*Yin) )**p2)
  ENDIF

!  IF(ielem.EQ.10) PRINT*,Y, cd, Gy, SoftParam, scalar1, Skk
  IF((Gy - SoftParam).GT.0.d0.AND.scalar1<0.d0.AND.skk>0.d0.AND.cd.LT.c2)THEN ! damage
!     IF((Gy - SoftParam).GT.0.d0)THEN ! damage

     
!     IF(i.EQ.357)   PRINT*,'DAMAGE************', 
      
     icnt = 0
     converge: DO
    
!!$        IF(Y.LT.Yin) THEN
!!$           H = 0.d0
!!$        ELSE
!!$           H = c2*EXP( - ( (Y - Yin)/(p1*Yin) )**p2)*p2*(Y-Yin)**(p2-1.d0)/(p1*Yin)**p2
!!$        ENDIF
!!$
!!$        CALL VecMatVecMul( strain, ScriptK, strain, 6, scalar1)         
!!$        scalar1 =  H/(1.d0 + 0.5d0*H*scalar1)
!!$
!!$        CALL VecMatVecMul( strain, ScriptJ, Delstrain, 6, scalar2)
!!$        
!!$        Dcd = - scalar1 * scalar2

        Dcd = Gy - SoftParam

!          IF(i.EQ.MonNd) PRINT*,'cd(1,i) , Dcd', cd(1,i) , Dcd,cd(1,i) + Dcd

!           IF(i.EQ.MonNd) PRINT*,'sc1, sc2',scalar1, scalar2

        cd = cd_n + Dcd


        cb = c2 - cd
           
        CALL A_D_tensors(Lo, Lm, L2, Mm, M2, cm, cb, cd, &
             Dbm,Dbb, Dbd, Dmm,Dmb,Dmd, Ddm, Ddb, Ddd, A2, Am, L_bar)
           
        
        CALL OperatorJ (Lm, L2, L2, L2, L_bar, Mm, M2, cm, cb, cd, &
             Dbm, Dbb, Dbd, Dmm,Dmb,Dmd, Ddm, Ddb, Ddd, &
             Am, A2, A2, A2, Km, K2, Gm, G2, vm, nu2, a_eta, a_zeta, &  
             ScriptJ, N, Ninv, R, Rinv, Q2, Qm, W)
              
! ScriptK = 2*ScriptJ*Rinv*(L2*Q2-W*Lm*Qm)*TRANSPOSE(A2) * L2 * Ninv*R
     
        ScriptK = MATMUL(L2,Q2) - MATMUL(W,MATMUL(Lm,Qm))
        ScriptK = MATMUL(Rinv,ScriptK)
        ScriptK = 2.d0*MATMUL(ScriptJ,ScriptK)
        ScriptK = MATMUL(ScriptK,TRANSPOSE(A2) )
        ScriptK = MATMUL(ScriptK,L2 )
        ScriptK = MATMUL(ScriptK,Ninv )
        ScriptK = MATMUL(ScriptK,R ) 
        
        nor = SQRT( (Dcd - Dcd_Old)**2)
        
        IF(icnt.EQ.0)THEN
           nor1 = nor
           nor = nor/nor1
        ELSE
           nor = nor/nor1
        ENDIF

!              _ T       _
!  Y = - 0.5 * E   * J * E               ... (24)
!

        CALL VecMatVecMul( strain, ScriptJ, strain, 6, Y)
        
        Y = -0.5d0 * Y
        

        IF(Y.LT.Yin)THEN
           Gy = 0.d0
        ELSE
           Gy = c2  - c2 * EXP( - ( (Y - Yin)/(p1*Yin) )**p2)
        ENDIF

!!$        IF(ielem.EQ.MonNd)THEN
!!$           PRINT*, 'nor,cd,gy =', nor, cd, Gy
!!$           PRINT*, 'Y,H,Dcd,kappaH',Y, H, Dcd,DelScriptK
!!$        ENDIF

        IF(Y.LT.Yin)THEN
           cd = cd_n
           EXIT
        ENDIF
           
        icnt = icnt + 1
        
        IF(icnt.GT.100)THEN
           PRINT*,'Did not converge'
           PRINT*,'Element',ielem
           STOP
        ENDIF
        IF( nor .LE. tol) EXIT! converged
        
        Dcd_Old = Dcd
        
     ENDDO converge

     IF(cd.GT.0.99999d0*c2) cd = c2

! Store cd at gauss point

     cb = c2 - cd
     
     IF(dcd.LT.0.d0.OR.cd.GT.c2) THEN
        PRINT*,'Dcd < 0 or cd > c2',ielem
        PRINT*,' Dcd, cd', Dcd, cd
        STOP
     ENDIF

     CALL A_D_tensors(Lo, Lm, L2, Mm, M2, cm, cb, cd, &
          Dbm,Dbb, Dbd, Dmm,Dmb,Dmd, Ddm, Ddb, Ddd, A2, Am, L_bar)
     
     
     CALL OperatorJ (Lm, L2, L2, L2, L_bar, Mm, M2, cm, cb, cd, &
          Dbm, Dbb, Dbd, Dmm,Dmb,Dmd, Ddm, Ddb, Ddd, &
          Am, A2, A2, A2, Km, K2, Gm, G2, vm, nu2, a_eta, a_zeta, &  
          ScriptJ, N, Ninv, R, Rinv, Q2, Qm, W)
        

! ScriptK = 2*ScriptJ*Rinv*(L2*Q2-W*Lm*Qm)*TRANSPOSE(A2) * L2 * Ninv*R
     
     ScriptK = MATMUL(L2,Q2) - MATMUL(W,MATMUL(Lm,Qm))
     ScriptK = MATMUL(Rinv,ScriptK)
     ScriptK = 2.d0*MATMUL(ScriptJ,ScriptK)
     ScriptK = MATMUL(ScriptK,TRANSPOSE(A2) )
     ScriptK = MATMUL(ScriptK,L2 )
     ScriptK = MATMUL(ScriptK,Ninv )
     ScriptK = MATMUL(ScriptK,R )  


!!$! COMPUTE:
!!$
!!$!        _ T        _        _ T      -
!!$! DK = - e  * J *  De - .5 * e  * K * e  * Dcd 
!!$!
!!$
!!$! a)
!!$!      _ T        _ 
!!$!      e  * J *  De 
!!$     
!!$
!!$     CALL VecMatVecMul( strain, ScriptJ, Delstrain,6, scalar1)
!!$
!!$! b)
!!$!         _ T      -
!!$!    .5 * e  * K * e  * Dcd 
!!$!
!!$
!!$
!!$     CALL VecMatVecMul( strain, ScriptK, strain,6, scalar2)
!!$
!!$
!!$! DK
!!$
!!$     DelScriptK = - scalar1 - 0.5d0* scalar2 * Dcd 
!!$        
        
     CALL VecMatVecMul( strain, ScriptJ, strain, 6, Y)
     
     Y = -0.5d0 * Y
     
     IF(Y.LT.Yin) THEN
        Gy = 0.d0
     ELSE
        Gy = c2  - c2 * EXP( - ( (Y - Yin)/(p1*Yin) )**p2)
     ENDIF
     

        ! Xi = DK * H
        
        ! Update Softening Parameter
        !    SoftParam = SoftParam + DelScriptK*H 

     IF(Gy.gt.SoftParam) SoftParam = Gy 
        
!        IF(i.EQ.MonNd) PRINT*,'cd after',cd,SoftParam, Gy

     ScriptAm = Am + MATMUL(Dmd,MATMUL(Ninv,R))
     ScriptAb = A2 + MATMUL(Dbd,MATMUL(Ninv,R))
     ScriptAd = A2 + MATMUL((Ddd-dident),MATMUL(Ninv,R))

  ENDIF ! end of damage loop

! Check micro and macro stresses

  tmpM = MATMUL(Ninv,R)
     
  DO ii = 1, 6
     tot = 0.d0
     DO jj = 1, 6
        tot = tot + tmpM(ii,jj)*strain(jj)
     ENDDO
     mu_d(ii) = tot
  ENDDO
  
  DO ii = 1, 6
     tot = 0.d0
     DO jj = 1, 6
        tot = tot + Am(ii,jj)*strain(jj)
     ENDDO
     strain_m(ii) = tot
  ENDDO
  
  DO ii = 1, 6
     tot = 0.d0
     DO jj = 1, 6
        tot = tot + Dmd(ii,jj)*mu_d(jj)
     ENDDO
     strain_m(ii) = strain_m(ii) + tot
  ENDDO
  
  
  DO ii = 1, 6
     tot = 0.d0
     DO jj = 1, 6
        tot = tot + A2(ii,jj)*strain(jj)
     ENDDO
     strain_b(ii) = tot
  ENDDO
  
  DO ii = 1, 6
     tot = 0.d0
     DO jj = 1, 6
        tot = tot + Dbd(ii,jj)*mu_d(jj)
     ENDDO
     strain_b(ii) = strain_b(ii) + tot
  ENDDO
  
  DO ii = 1, 6
     tot = 0.d0
     DO jj = 1, 6
        tot = tot + A2(ii,jj)*strain(jj)
     ENDDO
     strain_d(ii) = tot
  ENDDO
  
  DO ii = 1, 6
     tot = 0.d0
     DO jj = 1, 6
        tot = tot + Ddd(ii,jj)*mu_d(jj)
     ENDDO
     strain_d(ii) = strain_d(ii) + tot
  ENDDO
  
  DO ii = 1, 6
     tot = 0.d0
     DO jj = 1, 6
        tot = tot + Lm(ii,jj)*strain_m(jj) 
     ENDDO
     stress_m(ii) =  tot
  ENDDO
  
  DO ii = 1, 6
     tot = 0.d0
     DO jj = 1, 6
        tot = tot + L2(ii,jj)*strain_b(jj) 
     ENDDO
     stress_b(ii) =  tot
  ENDDO
  
  DO ii = 1, 6
     tot = 0.d0
     DO jj = 1, 6
        tot = tot + L2(ii,jj)*(strain_d(jj)-mu_d(jj))
     ENDDO
     stress_d(ii) =  tot
  ENDDO
  
  
! calculate the stress            -1
!                        [S] = [C]  {E}
! 

!!$  IF(ielem.EQ.MonNd)THEN
!!$
!!$     DO ii = 1, 6
!!$        PRINT*,llbar(ii,:)
!!$     ENDDO
!!$     PRINT*,'**',cb,cd,cm
!!$     
!!$  ENDIF

  LLbar = cm*MATMUL(Lm,ScriptAm) + cb*MATMUL(L2,ScriptAb) + cd*MATMUL(L2,ScriptAd)
  
  S11 = LLbar(1,1)*strain(1) + LLbar(1,2)*strain(2) + LLbar(1,3)*strain(3)
  S22 = LLbar(2,1)*strain(1) + LLbar(2,2)*strain(2) + LLbar(2,3)*strain(3) 
  S33 = LLbar(3,1)*strain(1) + LLbar(3,2)*strain(2) + LLbar(3,3)*strain(3)
  S12 = strain(6)*LLbar(4,4)
  S23 = strain(4)*LLbar(5,5)
  S13 = strain(5)*LLbar(6,6) 
  
  IF(ABS(S33 - (cm*stress_m(3) + cb*stress_b(3) + cd*stress_d(3))).GT.0.0001d0)THEN
     PRINT*,cd, S33, cm*stress_m(3) + cb*stress_b(3) + cd*stress_d(3)
  ENDIF
 

END SUBROUTINE Matous_Const_Model

