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
SUBROUTINE v3d4_r_bar(d_bar,R_bar,numnp,NumElv, &
     & lmcstet,meshcoor,nstart,nend)

!!****f* Rocfrac/Rocfrac/Source/v3d4_r_bar.f90
!!
!!  NAME
!!    v3d4_r_bar
!!
!!  FUNCTION
!!                 _
!!     Caluculates R for the ALE forumulation for
!!     the 4-node tetrahedral
!! 
!!  INPUTS
!!     d_bar -- Displacement from mesh motion and
!!            surface regression
!!     numnp -- Number of Nodes
!!     NumElv -- Number of volumetric elements
!!     lmcstet -- Element connectivity
!!     meshcoor -- mesh coordinates
!!     nstart -- starting element number
!!     nend -- ending element number 
!!
!!  OUTPUT
!!              _
!!      Rbar -- R additional internal force from ALE
!!
!!****

  IMPLICIT NONE

!--  Number of nodes
  INTEGER :: numnp 
!--  Nodal mesh displacements
  REAL*8, DIMENSION(1:3*numnp) :: d_bar
!--  Internal force
  REAL*8, DIMENSION(1:3*numnp) :: R_bar       
!--  Coordinates of position w.r.t. to global basis of nodes
  REAL*8, DIMENSION(1:3,1:numnp) ::  meshcoor 
!--  Number of LSTs
  INTEGER ::  NumElv
!--  Tetrahedra conectivity table
  INTEGER, DIMENSION(1:4,1:NumElv) :: lmcstet
!--  Displacement holding variable
  REAL*8, DIMENSION(1:12) ::  disp

  INTEGER :: nstart,nend

!--  Coordinate holding variables
  REAL*8 :: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
!--  Corrdinate subtractions
  REAL*8 :: x14, x24, x34, y14, y24, y34, z14, z24, z34
!--  
  REAL*8 :: c11, c12, c13, c21, c22, c23, c31, c32, c33

  REAL*8 :: g11, g14, g17
  REAL*8 :: g44, g47
  REAL*8 :: g77
  
  REAL*8 :: den
  
  REAL*8 :: Vx6
  
  INTEGER :: i, j, i3, j3, ielem
  INTEGER :: nd1, nd2, nd3, nd4
  INTEGER :: nx, ny, nz
  INTEGER :: k
  
  REAL*8, DIMENSION(1:10) :: elstiff
  
  INTEGER, DIMENSION(1:4,1:4) :: map = &
       &     RESHAPE((/1,2,3,4,2,5,6,7,3,6,8,9,4,7,9,10/),(/4,4/))

  DO ielem = nstart, nend
     
     nd1 = lmcstet(1,ielem)
     nd2 = lmcstet(2,ielem)
     nd3 = lmcstet(3,ielem)
     nd4 = lmcstet(4,ielem)
     
     x1 = meshcoor(1,nd1)
     y1 = meshcoor(2,nd1)
     z1 = meshcoor(3,nd1)
     x2 = meshcoor(1,nd2)
     y2 = meshcoor(2,nd2)
     z2 = meshcoor(3,nd2)
     x3 = meshcoor(1,nd3)
     y3 = meshcoor(2,nd3)
     z3 = meshcoor(3,nd3)
     x4 = meshcoor(1,nd4)
     y4 = meshcoor(2,nd4)
     z4 = meshcoor(3,nd4)
     
     DO i = 1, 4
        i3 = lmcstet(i,ielem)*3
        nx = i3 - 2
        ny = i3 - 1
        nz = i3
        disp(i*3-2) = d_bar(nx)
        disp(i*3-1) = d_bar(ny)
        disp(i*3  ) = d_bar(nz)
     ENDDO
         
     x14 = x1 - x4
     x24 = x2 - x4
     x34 = x3 - x4
     y14 = y1 - y4
     y24 = y2 - y4
     y34 = y3 - y4
     z14 = z1 - z4
     z24 = z2 - z4
     z34 = z3 - z4

     c11 =    y24*z34 - z24*y34
     c12 = -( y14*z34 - z14*y34 )
     c13 =    y14*z24 - z14*y24
     c21 = -( x24*z34 - z24*x34 )
     c22 =    x14*z34 - z14*x34
     c23 = -( x14*z24 - z14*x24 )
     c31 =    x24*y34 - y24*x34
     c32 = -( x14*y34 - y14*x34 )
     c33 =    x14*y24 - y14*x24

     g11 = c11**2 + c21**2 + c31**2
     g14 = c11*c12 + c21*c22 + c31*c32
     g17 = c11*c13 + c21*c23 + c31*c33
     
     g44 = c12**2 + c22**2 + c32**2
     g47 = c12*c13 + c22*c23 + c32*c33
     
     g77 = c13**2 + c23**2 + c33**2

! Mathematica output was originally det*( )/6. where det is of the Jacobi,
! but the det Jacobi for the 4-node constant strain tetradedral 
! is equal to  6 * volume
                                ! corresponding term in upper triangle
     elstiff(1) = g11       ! 1
     elstiff(2) = g14       ! 4
     elstiff(3) = g17       ! 7
     elstiff(4) = -g11 - g14 - g17 ! 10
     elstiff(5)= g44       ! 34 ! 13
     elstiff(6)= g47       ! 37 ! 14
     elstiff(7)= -g14 - g44 - g47 !40 ! 15
     elstiff(8)= g77       ! 58 ! 8
     elstiff(9)= -g17 - g47 - g77 ! 61 ! 9
     elstiff(10)= g11 + 2.d0*g14 + 2.d0*g17 + g44 + 2.d0*g47 + g77 ! 73 ! 10

! need the minus sign in order to get the same as the original
         
     Vx6 = -( x14*c11 + y14*c21 + z14*c31 )

!        Vx6 = -( x14*y24*z34 - x14*z24*y34 - y14*x24*z34 + y14*z24*x34
!     $       + z14*x24*y34 - z14*y24*x34)

     den = 1.d0/(6.d0*Vx6)
!     elstiff(1:10) = elstiff(1:10)*den
     elstiff(1) = elstiff(1)*den
     elstiff(2) = elstiff(2)*den
     elstiff(3) = elstiff(3)*den
     elstiff(4) = elstiff(4)*den
     elstiff(5)=  elstiff(5)*den
     elstiff(6)=  elstiff(6)*den
     elstiff(7)=  elstiff(7)*den
     elstiff(8)=  elstiff(8)*den
     elstiff(9)=  elstiff(9)*den
     elstiff(10)= elstiff(10)*den

     DO i = 1, 4
        nz = lmcstet(i,ielem)*3
        nx = nz - 2
        ny = nz - 1
        DO j = 1, 4
           j3 = j*3
           k = map(j,i)
           R_bar(nx) = R_bar(nx) + elstiff(k) * disp(j3-2)
           R_bar(ny) = R_bar(ny) + elstiff(k) * disp(j3-1)
           R_bar(nz) = R_bar(nz) + elstiff(k) * disp(j3)
        ENDDO
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE v3d4_r_bar

