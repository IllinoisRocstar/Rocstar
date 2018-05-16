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
SUBROUTINE VOL_ELEM_MAT_MATOUS(glb)
!
!     Forms the material compliance matix, [C] for an isotropic material
!                    -1
!     Note: [E] = [C]    and {epsilon} = [C]{sigma}
!

  USE ROCSTAR_RocFrac

  IMPLICIT NONE
  TYPE(ROCFRAC_GLOBAL) :: glb

  REAL*8 :: E_l, nu_l

  INTEGER :: i,j                 ! loop counter

!-----Create material constant matrices



  ALLOCATE(glb%L_tensor(1:6,1:6,1:glb%NumMatVol_Part))
  ALLOCATE(glb%M_tensor(1:6,1:6,1:glb%NumMatVol_Part))

  glb%L_tensor = 0.d0
  glb%M_tensor = 0.d0

  DO i = 1, glb%NumMatVol_Part


! -- [E]
!!$     ci(1,i) = e(i)*(1.d0 - xnu(i))/ &
!!$          ( (1.d0+xnu(i))*(1.d0-2.d0*xnu(i)) )
!!$     ci(2,i) = e(i)*xnu(i)/( (1.d0+xnu(i))*(1.d0-2.d0*xnu(i)) )
!!$     ci(3,i) = ci(1,i)
!!$     ci(4,i) = ci(2,i)
!!$     ci(5,i) = ci(2,i)
!!$     ci(6,i) = ci(1,i)
!!$     ci(7,i) = e(i)/( 2.d0 * (1.d0 + xnu(i)) )
!!$     ci(8,i) = ci(7,i)
!!$     ci(9,i) = ci(7,i)

! -- [C]
     glb%M_tensor(1,1,i) = 1.d0 / glb%E1(i)
     glb%M_tensor(1,2,i) = - glb%nu12(i) / glb%E2(i)
     glb%M_tensor(1,3,i) = - glb%nu13(i) / glb%E3(i)
     glb%M_tensor(2,1,i) = glb%M_tensor(1,2,i) 
     glb%M_tensor(2,2,i) = 1.d0 / glb%E2(i)
     glb%M_tensor(2,3,i) = - glb%nu23(i) / glb%E3(i)
     glb%M_tensor(3,1,i) = glb%M_tensor(1,3,i)
     glb%M_tensor(3,2,i) = glb%M_tensor(2,3,i)
     glb%M_tensor(3,3,i) = 1.d0 / glb%E3(i)
     glb%M_tensor(4,4,i) = 1.d0/glb%G12(i)
     glb%M_tensor(5,5,i) = 1.d0/glb%G13(i)
     glb%M_tensor(6,6,i) = 1.d0/glb%G23(i)

! -- 

     CALL invert2(glb%M_tensor(:,:,i), glb%L_tensor(:,:,i), 6)


     PRINT*,' Stiffness Matrix of Material',i
     DO j= 1, 6
        PRINT*,glb%L_tensor(j,:,i)
     ENDDO
     
     PRINT*,' Compliance Matrix of Material',i
     DO j= 1, 6
        PRINT*,glb%M_tensor(j,:,i)
     ENDDO



!!$     PRINT*, 'E = 1/M(1,1) = ', E_l
!!$     PRINT*, 'G = 1/M(4,4) = ', 1.d0/glb%M_tensor(4,4,i)
!!$     PRINT*, 'nu = -M(1,2)*E = ', nu_l

     ! not true for othotropic material: fix

     glb%ShrMod(i) = 1.d0/glb%M_tensor(4,4,i)
     glb%PoisRat(i) = -glb%M_tensor(1,2,i)*1.d0/glb%M_tensor(1,1,i)
     glb%BulkMod(i) = glb%E1(i)/(3.*(1.-2.*glb%PoisRat(i)))

  ENDDO

!  PRINT*,'HOMOGENEOUS MEDIUM : matrix[2], fiber[1] | cm =', cm, ': cf =', cb, ' : cd =', cd
  PRINT*,'Comparison medium Lo = a1*Lm + a2*L2 : a1 =',glb%alpha1, ' | a2 =', glb%alpha2

  
  CALL CompositeStiffnes(glb%L_tensor(:,:,1), glb%L_tensor(:,:,2), glb%M_tensor(:,:,1), glb%M_tensor(:,:,2), glb%cm, glb%c2, &
       glb%alpha1, glb%alpha2, glb%L_bar, glb%M_bar, glb%Lo)

  PRINT*,'Stiffness Matrix of Homogeneous Medium'

  DO j= 1, 6
     PRINT*,glb%L_bar(j,:)
  ENDDO

  PRINT*,'Compliance Matrix of Homogeneous Medium'

    DO j= 1, 6
     PRINT*,glb%M_bar(j,:)
  ENDDO

  PRINT*,''

  PRINT*,'E = 1/M(1,1) = ', 1.d0/glb%M_bar(1,1), 'G = 1/M[4][4] = ', &
       1.d0/glb%M_bar(4,4), 'nu = -M(1,2)*E =,', -glb%M_bar(1,2)*1./glb%M_bar(1,1)

  E_l = 1.d0/glb%M_bar(1,1)
  nu_l = -glb%M_bar(1,2)*1.d0/glb%M_bar(1,1)

  PRINT*,E_l*(1.d0-nu_l),glb%rho(1),(1.d0+nu_l)/(1.d0-2.d0*nu_l)
  glb%cd_fastest = SQRT( E_l*(1.d0-nu_l)/glb%rho(1)/(1.d0+nu_l)/(1.d0-2.d0*nu_l) )

!  ALLOCATE( glb%cd(1:4,1:glb%NumElVol) )

  PRINT*,'done'

  RETURN
END SUBROUTINE VOL_ELEM_MAT_MATOUS

