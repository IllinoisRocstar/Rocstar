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
SUBROUTINE FluidPressLoad(NumNP,R_ex, &
     InterfaceSFNumElems, InterfaceSFNumNodes, &
     InterfaceSFElemConn, &
     MapNodeSF,LwrBnd,UppBnd, coor, Disp, MapSFElVolEl,&
     ElConnVol,iElType,NumElVol,TractPress )

  IMPLICIT NONE

!!****f* Rocfrac/Rocfrac/Source/FluidPressLoad.f90
!!
!!  NAME
!!     FluidPressLoad
!!
!!  FUNCTION
!!
!!    Transforms the pressure given in the deformed state to
!!    the undeformed configuration. This subroutine is for
!!    the formulation where all quantities are with respect to
!!    the undeformed configuration.  Assumes the 
!!    pressure (tractions) are constant over the surface of
!!    the triangle.
!!
!!    1-3 for 4 node tet (i.e. 3 node triangles)
!!    4-6 for 10 node tet (i.e. 6 node triangles, mid-side nodes get traction)
!!
!!****


  INTEGER :: numnp          ! number of nodes
  INTEGER :: LwrBnd,UppBnd

  INTEGER :: iElType,NumElVol
  INTEGER, DIMENSION(1:iElType,1:NumElVol) :: ElConnVol

  INTEGER :: InterfaceSFNumElems,InterfaceSFNumNodes
  INTEGER, DIMENSION(1:UppBnd,1:InterfaceSFNumElems) :: InterfaceSFElemConn
  INTEGER, DIMENSION(1:InterfaceSFNumNodes) :: MapNodeSF
  REAL*8, DIMENSION(1:3*NumNp) :: Disp
  
  REAL*8, DIMENSION(1:3*numnp) :: R_ex ! external force
  
  INTEGER :: nx, ny, nz
  INTEGER :: i,j,k
  
  INTEGER,DIMENSION(1:3) :: TriConn
  REAL*8,DIMENSION(1:3)  :: UnifTract

! mesh coordinates
  REAL*8, DIMENSION(1:3,1:numnp) :: coor

  REAL*8 :: x0p,y0p,z0p
  REAL*8 :: x1p,y1p,z1p
  REAL*8 :: x2p,y2p,z2p
  REAL*8 :: x3p,y3p,z3p

! Components sides vector

  REAL*8 :: x1x0, y1y0, z1z0
  REAL*8 :: x2x0, y2y0, z2z0
  
  REAL*8 :: MagNorm, Area, at
  REAL*8 :: gauss =  0.333333333333333d0
  REAL*8 :: weight = 1.0d0
  REAL*8 :: xnorm,ynorm,znorm
  REAL*8 :: u1,v1,w1,u2,v2,w2,u3,v3,w3,u4,v4,w4
  REAL*8 :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
  REAL*8 :: J11,J12,J13,J21,J22,J23
  REAL*8 :: InvDenJplus
  REAL*8 :: Jplus11, Jplus12, Jplus21, Jplus22, Jplus31, Jplus32
  REAL*8 :: GradU11, GradU12, GradU13, GradU21, GradU22, GradU23
  REAL*8 :: GradU31, GradU32, GradU33
  REAL*8 :: F11T,F12T,F13T,F21T,F22T,F23T,F31T,F32T,F33T
  REAL*8 :: J1, UndefNorm(1:3)
  INTEGER :: IdVol

  INTEGER :: n1,n2,n3,n4
! -- dummy and counters
  INTEGER :: k1n1,k1n2,k1n3,k1n4,k2n1,k2n2,k2n3,k2n4
  INTEGER :: k3n1,k3n2,k3n3,k3n4
  REAL*8, DIMENSION(1:InterfaceSFNumElems) :: TractPress
!--  Coordinate subtractions
  REAL*8 :: x14, x24, x34, y14, y24, y34, z14, z24, z34
!--  Coordinate subtractions: These are to speed up B calculation
  REAL*8 :: x12, x13, y12, y13, z12, z13!--  

  REAL*8 :: c11, c21, c31
  REAL*8 :: Vx6, Vx6inv
!--   spacial derivatives
  REAL*8 :: B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12
!--   partial derivatives of the displacement 
  REAL*8 :: dudx,dvdy,dwdz,dudy,dvdx,dvdz,dwdy,dudz,dwdx
  INTEGER, DIMENSION(1:NumElVol) :: MapSFElVolEl
  REAL*8 :: Jac3D
  REAL*8 :: defNorm(1:3)
  REAL*8 :: denFinvT,magnormdef,area0
  REAL*8 :: F11,F12,F13,F21,F22,F23,F31,F32,F33
  REAL*8 :: FinvT11, FinvT12, FinvT13, FinvT21, FinvT22, FinvT23, FinvT31, FinvT32, FinvT33


  DO i = 1, InterfaceSFNumElems

! Compute the Gradient of U (disp) Volumetric 

     IdVol = MapSFElVolEl(i) ! Relates Surface id to Volumetric id

     n1 = ElConnVol(1,IdVol)
     n2 = ElConnVol(2,IdVol)
     n3 = ElConnVol(3,IdVol)
     n4 = ElConnVol(4,IdVol)

     k3n1 = 3*n1
     k3n2 = 3*n2
     k3n3 = 3*n3
     k3n4 = 3*n4
     
     k2n1 = k3n1 - 1
     k2n2 = k3n2 - 1
     k2n3 = k3n3 - 1
     k2n4 = k3n4 - 1
     
     k1n1 = k3n1 - 2
     k1n2 = k3n2 - 2
     k1n3 = k3n3 - 2
     k1n4 = k3n4 - 2

     x1 = coor(1,n1)
     x2 = coor(1,n2)
     x3 = coor(1,n3)
     x4 = coor(1,n4)
     y1 = coor(2,n1)
     y2 = coor(2,n2)
     y3 = coor(2,n3)
     y4 = coor(2,n4)
     z1 = coor(3,n1)
     z2 = coor(3,n2)
     z3 = coor(3,n3)
     z4 = coor(3,n4)
    
 ! k#n# dummy variables replaces:
     u1 = Disp(k1n1)           ! (3*n1-2)
     u2 = Disp(k1n2)           ! (3*n2-2)
     u3 = Disp(k1n3)           ! (3*n3-2)
     u4 = Disp(k1n4)           ! (3*n4-2)
     v1 = Disp(k2n1)           ! (3*n1-1)
     v2 = Disp(k2n2)           ! (3*n2-1)
     v3 = Disp(k2n3)           ! (3*n3-1)
     v4 = Disp(k2n4)           ! (3*n4-1)
     w1 = Disp(k3n1)           ! (3*n1)
     w2 = Disp(k3n2)           ! (3*n2)
     w3 = Disp(k3n3)           ! (3*n3)
     w4 = Disp(k3n4)           ! (3*n4)

     x12 = x1 - x2          ! not used in vol. calc
     x13 = x1 - x3          ! not used in vol. calc
     x14 = x1 - x4
     x24 = x2 - x4
     x34 = x3 - x4
     y12 = y1 - y2          ! not used in vol. calc
     y13 = y1 - y3          ! not used in vol. calc
     y14 = y1 - y4
     y24 = y2 - y4
     y34 = y3 - y4
     z12 = z1 - z2          ! not used in vol. calc
     z13 = z1 - z3          ! not used in vol. calc
     z14 = z1 - z4
     z24 = z2 - z4
     z34 = z3 - z4

     c11 =    y24*z34 - z24*y34
     c21 = -( x24*z34 - z24*x34 )
     c31 =    x24*y34 - y24*x34
     
     Vx6 = -( x14*c11 + y14*c21 + z14*c31 )
     
     Vx6inv = 1.d0 / Vx6

! See the maple worksheet 'V3D4.mws' for the derivation of [B]
! NOTE: Factored for a more equivalent/compact form then maple's

     B1  = (y34*z24 - y24*z34) * Vx6inv
     B2  = (z34*x24 - z24*x34) * Vx6inv
     B3  = (x34*y24 - x24*y34) * Vx6inv
     B4  = (y13*z14 - y14*z13) * Vx6inv
     B5  = (z13*x14 - z14*x13) * Vx6inv
     B6  = (x13*y14 - x14*y13) * Vx6inv
     B7  = (y14*z12 - y12*z14) * Vx6inv
     B8  = (z14*x12 - z12*x14) * Vx6inv
     B9  = (x14*y12 - x12*y14) * Vx6inv
     B10 = (y12*z13 - y13*z12) * Vx6inv
     B11 = (z12*x13 - z13*x12) * Vx6inv
     B12 = (x12*y13 - x13*y12) * Vx6inv
!
! Calculate displacement gradient (H)
!
     dudx = B1*u1 + B4*u2 + B7*u3 + B10*u4
     dvdy = B2*v1 + B5*v2 + B8*v3 + B11*v4
     dwdz = B3*w1 + B6*w2 + B9*w3 + B12*w4
     dudy = B2*u1 + B5*u2 + B8*u3 + B11*u4
     dvdx = B1*v1 + B4*v2 + B7*v3 + B10*v4
     dvdz = B3*v1 + B6*v2 + B9*v3 + B12*v4
     dwdy = B2*w1 + B5*w2 + B8*w3 + B11*w4
     dudz = B3*u1 + B6*u2 + B9*u3 + B12*u4
     dwdx = B1*w1 + B4*w2 + B7*w3 + B10*w4
! 
! deformation gradients F
!
     F11 = 1.d0 + dudx
     F12 = dudy
     F13 = dudz
     F21 = dvdx
     F22 = 1.d0 + dvdy
     F23 = dvdz
     F31 = dwdx
     F32 = dwdy
     F33 = 1.d0 + dwdz

!
!  
!

     denFinvT = F11*(F22*F33-F32*F23) - F12*(F21*F33+F31*F23) + F13*(F21*F32-F31*F22)

     FinvT11 =  (F22*F33-F32*F23)/denFinvT
     FinvT12 = -(F21*F33-F31*F23)/denFinvT
     FinvT13 =  (F21*F32-F31*F22)/denFinvT
     FinvT21 = -(F12*F33-F32*F13)/denFinvT
     FinvT22 =  (F11*F33-F31*F13)/denFinvT
     FinvT23 = -(F11*F32-F31*F12)/denFinvT
     FinvT31 =  (F12*F23-F22*F13)/denFinvT
     FinvT32 = -(F11*F23-F21*F13)/denFinvT
     FinvT33 =  (F11*F22-F21*F12)/denFinvT

     Jac3D = denFinvT 
     
!  Nodes of the Triangle ( Vertex Nodes ) 

     TriConn(1:3) = InterfaceSFElemConn(1:3,i)

! 1) Calculate the Normal in the Deformed Configuration

!  Uses the fact that the norm of the cross product vector
!  is the area of the parallelogram they form.  The triangle they
!  form has half that area.

     x1 = coor(1,MapNodeSF(TriConn(1)))
     y1 = coor(2,MapNodeSF(TriConn(1)))
     z1 = coor(3,MapNodeSF(TriConn(1)))
     
     x2 = coor(1,MapNodeSF(TriConn(2)))
     y2 = coor(2,MapNodeSF(TriConn(2)))
     z2 = coor(3,MapNodeSF(TriConn(2)))
     
     x3 = coor(1,MapNodeSF(TriConn(3)))
     y3 = coor(2,MapNodeSF(TriConn(3)))
     z3 = coor(3,MapNodeSF(TriConn(3)))

     u1 = Disp(3*MapNodeSF(TriConn(1))-2)
     v1 = Disp(3*MapNodeSF(TriConn(1))-1)
     w1 = Disp(3*MapNodeSF(TriConn(1))  )
     
     u2 = Disp(3*MapNodeSF(TriConn(2))-2)
     v2 = Disp(3*MapNodeSF(TriConn(2))-1)
     w2 = Disp(3*MapNodeSF(TriConn(2))  )
     
     u3 = Disp(3*MapNodeSF(TriConn(3))-2)
     v3 = Disp(3*MapNodeSF(TriConn(3))-1)
     w3 = Disp(3*MapNodeSF(TriConn(3))  )

     x0p = x1
     y0p = y1
     z0p = z1 
     
     x1p = x2
     y1p = y2 
     z1p = z2
     
     x2p = x3
     y2p = y3
     z2p = z3

! Find vector componets of the sides

     x1x0 = x1p - x0p
     y1y0 = y1p - y0p
     z1z0 = z1p - z0p

     x2x0 = x2p - x0p
     y2y0 = y2p - y0p
     z2z0 = z2p - z0p     

!  Take the cross product

     x3p = y1y0 * z2z0  -  z1z0  *  y2y0
     y3p = z1z0 * x2x0  -  x1x0  *  z2z0
     z3p = x1x0 * y2y0  -  y1y0  *  x2x0
         
!  Computes the Euclidean norm of a vector in 3D Deformed Config.
         
     MagNorm = SQRT( x3p*x3p + y3p*y3p + z3p*z3p )

     xnorm = x3p/MagNorm
     ynorm = y3p/MagNorm
     znorm = z3p/MagNorm

!     area0 = 0.5*MagNorm

! little n (deformed normal)

     defNorm(1) = Jac3D * (FinvT11*Xnorm+FinvT12*Ynorm+FinvT13*Znorm)
     defNorm(2) = Jac3D * (FinvT21*Xnorm+FinvT22*Ynorm+FinvT23*Znorm)
     defNorm(3) = Jac3D * (FinvT31*Xnorm+FinvT32*Ynorm+FinvT33*Znorm)

     MagNormDef = SQRT(  defNorm(1)**2 +  defNorm(2)**2 +  defNorm(3)**2 )

     defNorm(1) = defNorm(1)/MagNormDef
     defNorm(2) = defNorm(2)/MagNormDef
     defNorm(3) = defNorm(3)/MagNormDef

     x0p = x1 + u1
     y0p = y1 + v1
     z0p = z1 + w1
     
     x1p = x2 + u2
     y1p = y2 + v2
     z1p = z2 + w2
     
     x2p = x3 + u3
     y2p = y3 + v3
     z2p = z3 + w3

! Find vector componets of the sides

     x1x0 = x1p - x0p
     y1y0 = y1p - y0p
     z1z0 = z1p - z0p

     x2x0 = x2p - x0p
     y2y0 = y2p - y0p
     z2z0 = z2p - z0p     

!  Take the cross product

     x3p = y1y0 * z2z0  -  z1z0  *  y2y0
     y3p = z1z0 * x2x0  -  x1x0  *  z2z0
     z3p = x1x0 * y2y0  -  y1y0  *  x2x0
         
!  Computes the Euclidean norm of a vector in 3D Deformed Config.

     MagNormDef = SQRT( x3p*x3p + y3p*y3p + z3p*z3p )

     Area = 0.5*MagNormDef

! Pressure * Jac2D * n * Area0 (Jac2D = Area/Area0)

     UnifTract(1) = -TractPress(i)*defNorm(1)*Area/3.d0
     UnifTract(2) = -TractPress(i)*defNorm(2)*Area/3.d0
     UnifTract(3) = -TractPress(i)*defNorm(3)*Area/3.d0

!  Assembly into global external load vector R_ex

     TriConn(1:3) = InterfaceSFElemConn(LwrBnd:UppBnd,i)

     DO j = 1, 3
        nz = MapNodeSF(TriConn(j))*3 ! volume node
        nx = nz - 2
        ny = nz - 1
        R_ex(nx) = R_ex(nx) + UnifTract(1)
        R_ex(ny) = R_ex(ny) + UnifTract(2)
        R_ex(nz) = R_ex(nz) + UnifTract(3)
     ENDDO

  ENDDO

END SUBROUTINE FluidPressLoad

