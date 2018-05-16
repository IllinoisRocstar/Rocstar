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
SUBROUTINE VOL_ELEM_MAT(e,xnu,ci,cj,numat_vol,Integration)

!!****f* Rocfrac/Source/vol_elem_mat.f90
!!
!!  NAME
!!    VOL_ELEM_MAT
!!
!!  FUNCTION
!!
!!     Forms the material compliance matix, [C] for an isotropic material
!!                    -1
!!     Note: [E] = [C]    and {epsilon} = [C]{sigma}
!!
!!  INPUTS
!!
!!     e -- Young's modulus
!!    xnu -- Possion's ratio
!!   numat_vol -- number of volumetric elements
!!   Integration -- = 1 pseudo-Reduced integration (split Cijkl)
!!                  = 0 Cijkl
!!
!!  OUTPUT
!!     ci -- elastic stiffness constants
!!     cj -- split stiffness constants
!!
!!****
!!
      IMPLICIT NONE
      INTEGER :: numat_vol         ! number of materials
!      INTEGER :: nplane            ! type of material
!                                   !   0 = isotropic
!                                   !   1 = orthotropic
!--   young's moduli
      REAL*8, DIMENSION(1:numat_vol) :: e
!--   Poisson's ratios
      REAL*8, DIMENSION(1:numat_vol) :: xnu
!--   elastic stiffness constants [E]
      REAL*8, DIMENSION(1:9,1:numat_vol) :: ci
!--   elastic stiffness constants [E]
      REAL*8, DIMENSION(1:9,1:numat_vol) :: cj
! --  
      REAL*8 :: xmu, xlambda
      REAL*8 :: caux
      integer :: Integration

      INTEGER :: i                 ! loop counter

!-----Create material constant matrices
      DO i = 1, numat_vol
!----- (isotropic)
         IF (Integration.EQ.0) THEN

            caux = e(i) / ( (1.d0+xnu(i))*(1.d0-2.d0*xnu(i)) )

! -- [E]
!
!        [ c1 c2 c4 ..    ] 
!        [    c3 c5 ..    ]
!        [       c6 ..    ]
!        [         c7     ]
!        [           c8   ]
!        [             c9 ]

            ci(1,i) = caux*(1.d0 - xnu(i))
            ci(2,i) = caux*xnu(i)
            ci(3,i) = ci(1,i)
            ci(4,i) = ci(2,i)
            ci(5,i) = ci(2,i)
            ci(6,i) = ci(1,i)
            ci(7,i) = e(i)/( 2.d0 * (1.d0 + xnu(i)) )
            ci(8,i) = ci(7,i)
            ci(9,i) = ci(7,i)


            cj(:,i) = 0.d0

         ELSE IF(Integration.EQ.1)THEN

! -- pseudo-Reduced integration (split Cijkl)
            ! mu
            xmu = e(i)/( 2.d0 * (1.d0 + xnu(i)) )
            ! lambda
            xlambda = e(i)*xnu(i)/( (1.d0+xnu(i))*(1.d0-2.d0*xnu(i)) )

! -- [Emu]
            ci(1,i) = 2.d0*xmu
            ci(2,i) = 0.d0
            ci(3,i) = ci(1,i)
            ci(4,i) = ci(2,i)
            ci(5,i) = ci(2,i)
            ci(6,i) = ci(1,i)
            ci(7,i) = xmu
            ci(8,i) = ci(7,i)
            ci(9,i) = ci(7,i)

! -- [Elambda]

            cj(1,i) = xlambda 
            cj(2,i) = xlambda
            cj(3,i) = cj(1,i)
            cj(4,i) = cj(2,i)
            cj(5,i) = cj(2,i)
            cj(6,i) = cj(1,i)
            cj(7,i) = 0.d0
            cj(8,i) = cj(7,i)
            cj(9,i) = cj(7,i)
            
!----- (orthotropic, not yet implemented)
!         ELSEIF (nplane.GE.1) THEN
!            PRINT*,'ERROR: Material properties not yet implemented'
!            PRINT*,'       Your Choices are: '
!            PRINT*,'                         nplane = 0, isotropic'
!            STOP
         ENDIF
      ENDDO

      RETURN
    END

