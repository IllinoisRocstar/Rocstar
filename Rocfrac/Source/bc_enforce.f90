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
SUBROUTINE bc_enforce(numbound,numnp,id,r,slope,prop,&
     vb, ab, v, a, d, delta, Rnet, xm, DampEnabled,CurrTime)

!!****f* Rocfrac/Rocfrac/Source/bc_enforce.f90
!!
!!  NAME
!!     bc_enforce
!!
!!  FUNCTION
!!
!!     Enforces the structural boundary conditions
!!
!!  INPUTS
!!     numbound -- number of nodes with enforced boundary 
!!            conditions
!!        numnp -- number of nodes
!!           id -- element id
!!            r -- imposed type of boundary condition
!!        slope -- loading amplitude slope
!!        prop  -- proportion of amplitude for loading
!!          vb  -- imposed velocity
!!          ab  -- impoesed acceleration
!!       delta  -- time increment
!!        Rnet  -- sum of forces
!!          xm  -- lumped nodal mass matrix
!! DampEnabled  -- flag for damping
!!     CurrTime -- current time
!!
!!  OUTPUT
!!           v  -- nodal velocity with bc
!!           a  -- nodal acceleration with bc
!!           d  -- nodal displacement with bc
!!****

      IMPLICIT NONE
      
      INTEGER :: numbound, numnp
      INTEGER,DIMENSION(1:4,numbound) :: id
      REAL*8, DIMENSION(1:3*numbound) :: vb, ab
      REAL*8, DIMENSION(1:3,numbound) :: r
      REAL*8 :: delta
      REAL*8, DIMENSION(1:3*numnp) :: Rnet
      REAL*8, DIMENSION(1:3*numnp) :: v, a,d
      REAL*8, DIMENSION(1:numnp) :: xm
      REAL*8 :: slope,prop,CurrTime

      INTEGER :: k1, k2, k3, k4
      INTEGER :: i, i1, i2, i3
      REAL*8 :: a1
      LOGICAL :: DampEnabled
      
      IF(.NOT.(DampEnabled))THEN

         DO i = 1,numbound
            k4 = id(1,i)
            k1 = k4*3 - 2
            k2 = k4*3 - 1
            k3 = k4*3
            i1 = i*3 - 2
            i2 = i*3 - 1
            i3 = i*3
            IF (id(2,i).EQ.0) THEN ! velocity imposed
               vb(i1) = r(1,i)*CurrTime*slope
               ab(i1) = r(1,i)*slope
            ELSE ! force imposed
               a1 = (r(1,i) + Rnet(k1) )*xm(k4)
               vb(i1) = vb(i1) + delta*(ab(i1) + a1)*0.5d0
               ab(i1) = a1
            ENDIF
            IF (id(3,i).EQ.0) THEN ! velocity imposed
               vb(i2) = r(2,i)*CurrTime*slope
               ab(i2) = r(2,i)*slope
            ELSE ! force imposed
               a1 = (r(2,i) + Rnet(k2) )*xm(k4)
               vb(i2) = vb(i2) + delta*(ab(i2)+a1)*0.5d0
               ab(i2) = a1
            ENDIF
            IF (id(4,i).EQ.0) THEN ! velocity imposed
               vb(i3) = r(3,i)*CurrTime*slope
               ab(i3) = r(3,i)*slope
            ELSE ! force imposed
               a1 = (r(3,i) + Rnet(k3) )*xm(k4)
               vb(i3) = vb(i3) + delta*(ab(i3)+a1)*0.5d0
               ab(i3) = a1
            ENDIF
            a(k1) = ab(i1)
            a(k2) = ab(i2)
            a(k3) = ab(i3)
            v(k1) = vb(i1)
            v(k2) = vb(i2)
            v(k3) = vb(i3)
         ENDDO
   ELSE
      DO i = 1,numbound
         k4 = id(1,i)
         k1 = k4*3 - 2
         k2 = k4*3 - 1
         k3 = k4*3
         IF (id(2,i).EQ.0) THEN
            d(k1) = r(1,i)
            v(k1) = r(1,i) ! only correct if = 0, otherwise i would have to store the old boundary velocity
         ENDIF
         IF (id(3,i).EQ.0) THEN
            d(k2) = r(2,i)
            v(k2) = r(2,i)
         ENDIF
         IF (id(4,i).EQ.0) THEN
            d(k3) = r(3,i)
            v(k3) = r(3,i)
         ENDIF
      ENDDO

   ENDIF
      
   RETURN
      
 END SUBROUTINE bc_enforce

