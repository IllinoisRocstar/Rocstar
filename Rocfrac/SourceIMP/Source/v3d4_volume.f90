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
SUBROUTINE v3d4_volume(coor,lmcstet,matcstet, &
     rho,numnp,numcstet,numat_vol,Disp,nstart,nend,TotalMass,TotalGeomVolp,TotalGeomUndefVolp,&
     NumVertx)

!!****f* Rocfrac/Rocfrac/Source/vol_elem_mat.f90
!!
!!  NAME
!!    VOL_ELEM_MAT
!!
!!  FUNCTION
!!     Caluculates the volume and surface area for tetraderal
!!     elements
!!
!!****

  IMPLICIT NONE
  integer :: NumVertx       ! order of tetrahedral 4 or 10
  INTEGER :: numnp          ! number of nodes
  INTEGER :: numcstet       ! number of CSTets
  INTEGER :: numat_vol      ! number of materials
  REAL*8 :: DT
!--   densities
  REAL*8, DIMENSION(1:numat_vol) :: rho
!--   material number for CSTet element
  INTEGER, DIMENSION(1:numcstet) :: matcstet
!--   connectivity table for CSTet elem
  INTEGER, DIMENSION(1:NumVertx,1:numcstet) :: lmcstet
!--   global coordinates
  REAL*8, DIMENSION(1:3,1:numnp) :: coor
!-- Displacement
  REAL*8, DIMENSION(1:3*numnp) :: Disp

  REAL*8 :: aa                 ! determinant of jacobian (2*area)
  REAL*8 :: x !,x1,x2,x3         ! dummy variable
  INTEGER :: n1,n2,n3,n4,n5,n6 ! nodes, and dummy vars
  INTEGER :: n7,n8,n9,n10
  INTEGER :: i,nstart,nend     ! loop counter
!--  Coordinate holding variable
  REAL*8 :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
!--  Coordinate subtractions
  REAL*8 :: x14, x24, x34, y14, y24, y34, z14, z24, z34
!--  
  REAL*8 :: c11, c21, c31
!--   6*volume and the volume      
  REAL*8 :: Vx6,volume,Vx6def
  REAL*8 :: TotalMass,TotalGeomVolp,TotalGeomUndefVolp
  real*8,SAVE :: cnt

  TotalGeomUndefVolp = 0.
  TotalGeomVolp = 0.
  
  DO i = nstart, nend

     n1 = lmcstet(1,i)
     n2 = lmcstet(2,i)
     n3 = lmcstet(3,i)
     n4 = lmcstet(4,i)

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
     c21 = -( x24*z34 - z24*x34 )
     c31 =    x24*y34 - y24*x34
     
     Vx6 = -( x14*c11 + y14*c21 + z14*c31 )
     
!        Vx6 = x2*y3*z4 - x2*y4*z3 - y2*x3*z4 + y2*x4*z3 + z2*x3*y4 
!     $        - z2*x4*y3 - x1*y3*z4 + x1*y4*z3 + x1*y2*z4 - x1*y2*z3
!     $        - x1*z2*y4 + x1*z2*y3 + y1*x3*z4 - y1*x4*z3 - y1*x2*z4 
!     $        + y1*x2*z3 + y1*z2*x4 - y1*z2*x3 - z1*x3*y4 + z1*x4*y3
!     $        + z1*x2*y4 - z1*x2*y3 - z1*y2*x4 + z1*y2*x3

     TotalGeomUndefVolp =  TotalGeomUndefVolp + Vx6/6.d0

!!$     IF(x.LT.0.d0)THEN
!!$        PRINT*,'ERROR'
!!$        PRINT*,'NEG, Volume... STOPPING'
!!$        PRINT*,'  ELEMENT =',i
!!$        PRINT*,'  NODES=',n1,n2,n3,n4
!!$        PRINT*,'  x-Coordinates:',x1,x2,x3,x4
!!$        PRINT*,'  y-Coordinates:',y1,y2,y3,y4
!!$        PRINT*,'  z-Coordinates:',z1,z2,z3,z4
!!$        STOP
!!$     ENDIF

     x1 = coor(1,n1) + Disp(3*n1-2)
     x2 = coor(1,n2) + Disp(3*n2-2)
     x3 = coor(1,n3) + Disp(3*n3-2)
     x4 = coor(1,n4) + Disp(3*n4-2)
     y1 = coor(2,n1) + Disp(3*n1-1)
     y2 = coor(2,n2) + Disp(3*n2-1)
     y3 = coor(2,n3) + Disp(3*n3-1)
     y4 = coor(2,n4) + Disp(3*n4-1)
     z1 = coor(3,n1) + Disp(3*n1)
     z2 = coor(3,n2) + Disp(3*n2)
     z3 = coor(3,n3) + Disp(3*n3)
     z4 = coor(3,n4) + Disp(3*n4)
     
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
     c21 = -( x24*z34 - z24*x34 )
     c31 =    x24*y34 - y24*x34

     Vx6def = -( x14*c11 + y14*c21 + z14*c31 )

! Note: The total mass is the TOTAL mass OF all components (not just the propellent)

     TotalGeomVolp =  TotalGeomVolp + Vx6def/6.d0

  ENDDO

  TotalMass = rho(1) * TotalGeomUndefVolp ! assuming propellent is property 1


!  WRITE(545,*) cnt*0.00000095367431640625, &
!       (TotalGeomVolp/(4.d0-cnt*0.00000095367431640625d0*0.0625d0)-1.d0)*100.d0

!  print*,cnt*0.00000095367431640625, &
!       (TotalGeomVolp/(4.d0-cnt*0.00000095367431640625d0*0.0625d0)-1.d0)*100.d0

  cnt = cnt + 1
  
  RETURN
END SUBROUTINE v3d4_volume

SUBROUTINE triangle_area_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3, area )

  IMPLICIT NONE
!
  REAL*8 area
  REAL*8 enorm_3d
  REAL*8 norm
  REAL*8 x1
  REAL*8 x2
  REAL*8 x3
  REAL*8 x4
  REAL*8 y1
  REAL*8 y2
  REAL*8 y3
  REAL*8 y4
  REAL*8 z1
  REAL*8 z2
  REAL*8 z3
  REAL*8 z4
!
  CALL cross_3d ( x2 - x1, y2 - y1, z2 - z1, x3 - x1, y3 - y1, z3 - z1, &
                  x4,      y4,      z4 )


  norm = enorm_3d ( x4, y4, z4 )

  area = 0.5d0 * norm

  RETURN
END SUBROUTINE triangle_area_3d

  SUBROUTINE cross_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!*******************************************************************************
!
!! CROSS_3D computes the cross product of two vectors in 3D.
!
!
!  Definition:
!
!    The cross product in 3D can be regarded as the determinant of the
!    symbolic matrix:
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, the coordinates of the vectors.
!
!    Output, real X3, Y3, Z3, the cross product vector.
!
  IMPLICIT NONE
!
  REAL*8 x1
  REAL*8 x2
  REAL*8 x3
  REAL*8 y1
  REAL*8 y2
  REAL*8 y3
  REAL*8 z1
  REAL*8 z2
  REAL*8 z3
!
  x3 = y1 * z2 - z1 * y2
  y3 = z1 * x2 - x1 * z2
  z3 = x1 * y2 - y1 * x2

  RETURN
END SUBROUTINE cross_3d

FUNCTION enorm_3d ( x1, y1, z1 )
!
!*******************************************************************************
!
!! ENORM_3D computes the Euclidean norm of a vector in 3D.
!
!
!  Modified:
!
!    27 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, the coordinates of the vector.
!
!    Output, real ENORM_3D, the Euclidean norm of the vector.
!
  implicit none
!
  real*8 enorm_3d
  real*8 x1
  real*8 y1
  real*8 z1
!
  enorm_3d = sqrt ( x1 * x1 + y1 * y1 + z1 * z1 )

  return
END FUNCTION enorm_3d

