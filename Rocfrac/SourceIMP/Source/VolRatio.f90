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
  SUBROUTINE VolRatio(n1, n2, n3, n4, AlphaR,coor,numnp,NdMassLump)

!!****f* Rocfrac/Rocfrac/Source/VolRatio.f90
!!
!!  NAME
!!    VolRatio
!!
!!  FUNCTION
!!    Calculates the Volume Ratio for the node based elements
!!    using either the circumcenter of centroid.
!!
!!  INPUTS
!!    n1,n2,n3,n4 -- node number
!!    coor -- nodeal coordinates
!!    numnp -- number of nodes
!!    NdMassLump -- = 1 then circumcenter method
!!                  = 0 then centroid method
!!
!!  OUTPUT
!!
!!    AlphaR -- Ratio of the Volume for a node
!!
!!****
!!     
    IMPLICIT NONE

    CHARACTER*4 :: myid_chr          ! output file number (a character)
    INTEGER :: j
    INTEGER :: j1l, j2l, j3l, j4l
    REAL*8 :: VolEl,TolV
    INTEGER, DIMENSION(1:4,1) :: vx
    REAL*8 :: XCN,YCN,ZCN,XPT,YPT,ZPT
    INTEGER :: NCONT
    REAL*8, DIMENSION(1:3) :: CirCentTet,x,y,z,xmd,ymd,zmd
    REAL*8, DIMENSION(1:4,1:3) :: CirCentTri
    INTEGER :: numnp
    REAL*8, DIMENSION(1:3,1:numnp) :: coor

    INTEGER :: n1, n2, n3, n4
    INTEGER :: in
    REAL*8, DIMENSION(1:4) :: AlphaR

    REAL*8 :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4

    REAL*8 :: x12mp,x13mp,x14mp,x23mp,x24mp,x34mp
    REAL*8 :: y12mp,y13mp,y14mp,y23mp,y24mp,y34mp
    REAL*8 :: z12mp,z13mp,z14mp,z23mp,z24mp,z34mp

!--  Coordinate subtractions
    REAL*8 :: x14v,x24v, x34v, y14v, y24v, y34v, z14v, z24v, z34v
    REAL*8 :: v11,v21,v31
    REAL*8 :: vol

    REAL*8, DIMENSION(1:3) :: norm1, norm2, norm3, norm4, normC1, normC2, normC3,CentroidTet
    REAL*8, DIMENSION(1:3) :: vNorm12, vNorm13, vNorm14, vNorm23, vNorm24, vNorm34
    REAL*8 :: offset12, offset13, offset14, offset23, offset24, offset34

    INTEGER :: same1,same2,same3
    REAL*8 :: offset1, offset2, offset3,offset4, offsetC1, offsetC2, offsetC3,CalcOffset,Offset
    CHARACTER*61 :: qhullstat
    INTEGER :: iopt
    REAL*8 :: aj1, aj2, aj3, aj4
    real*8 feasiblePt(0:2)
    REAL*8, DIMENSION(0:3,0:6) :: Planes
    integer :: NdMassLump

!
!     SET TOLERANCE FOR MINIMUM VOLUME
!
    Tolv       = 1.0E-13

    IF(NdMassLump.EQ.1)THEN
  
       x1 = coor(1,n1) ! get vertex coords.
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
! Volume of Tetrahedra

       x14v = x1 - x4
       x24v = x2 - x4
       x34v = x3 - x4
       y14v = y1 - y4
       y24v = y2 - y4
       y34v = y3 - y4
       z14v = z1 - z4
       z24v = z2 - z4
       z34v = z3 - z4

       v11 =    y24v*z34v - z24v*y34v
       v21 = -( x24v*z34v - z24v*x34v )
       v31 =    x24v*y34v - y24v*x34v
       
       VolEl = -( x14v*v11 + y14v*v21 + z14v*v31 )/6.d0
          
! Coordinates of the edge mid-points

       x12mp = 0.5d0*(x1 + x2)
       x13mp = 0.5d0*(x1 + x3)
       x14mp = 0.5d0*(x1 + x4)
       x23mp = 0.5d0*(x2 + x3)
       x24mp = 0.5d0*(x2 + x4)
       x34mp = 0.5d0*(x3 + x4)
       y12mp = 0.5d0*(y1 + y2)
       y13mp = 0.5d0*(y1 + y3)
       y14mp = 0.5d0*(y1 + y4)
       y23mp = 0.5d0*(y2 + y3)
       y24mp = 0.5d0*(y2 + y4)
       y34mp = 0.5d0*(y3 + y4)
       z12mp = 0.5d0*(z1 + z2)
       z13mp = 0.5d0*(z1 + z3)
       z14mp = 0.5d0*(z1 + z4)
       z23mp = 0.5d0*(z2 + z3)
       z24mp = 0.5d0*(z2 + z4)
       z34mp = 0.5d0*(z3 + z4)

!  normals of the faces
       
       CALL PlaneNorm(x3,y3,z3,x2,y2,z2,x1,y1,z1,norm1(1),norm1(2),norm1(3),same1)
       offset1 = CalcOffset(norm1(1),norm1(2),norm1(3),x1,y1,z1)
       CALL PlaneNorm(x1,y1,z1,x2,y2,z2,x4,y4,z4,norm2(1),norm2(2),norm2(3),same1)
       offset2 = CalcOffset(norm2(1),norm2(2),norm2(3),x1,y1,z1)
       CALL PlaneNorm(x2,y2,z2,x3,y3,z3,x4,y4,z4,norm3(1),norm3(2),norm3(3),same1)
       offset3 = CalcOffset(norm3(1),norm3(2),norm3(3),x2,y2,z2)
       CALL PlaneNorm(x1,y1,z1,x4,y4,z4,x3,y3,z3,norm4(1),norm4(2),norm4(3),same1)
       offset4 = CalcOffset(norm4(1),norm4(2),norm4(3),x1,y1,z1)

       CALL PlaneNorm3D( x1, y1, z1, x2, y2 , z2, vNorm12 )
       offset12 = Offset(vNorm12,x12mp,y12mp,z12mp)
       
       CALL PlaneNorm3D( x1, y1, z1, x3, y3 , z3, vNorm13 )
       offset13 = Offset(vNorm13,x13mp,y13mp,z13mp)
       
       CALL PlaneNorm3D( x1, y1, z1, x4, y4 , z4, vNorm14 )
       offset14 = Offset(vNorm14,x14mp,y14mp,z14mp)
       
       CALL PlaneNorm3D( x2, y2, z2, x3, y3 , z3, vNorm23 )
       offset23 = Offset(vNorm23,x23mp,y23mp,z23mp)
       
       CALL PlaneNorm3D( x2, y2, z2, x4, y4 , z4, vNorm24 )
       offset24 = Offset(vNorm24,x24mp,y24mp,z24mp)
       
       CALL PlaneNorm3D( x3, y3, z3, x4, y4 , z4, vNorm34 )
       offset34 = Offset(vNorm34,x34mp,y34mp,z34mp)

       planes(0:2,0)  = -norm1(1:3)
       planes(  3,0)  = offset1
       planes(0:2,1)  = -norm2(1:3)
       planes(  3,1)  = offset2
       planes(0:2,2)  = -norm3(1:3)
       planes(  3,2)  = offset3
       planes(0:2,3)  = -norm4(1:3)
       planes(  3,3)  = offset4

! Node 1 Volume

       planes(0:2,4) = -vNorm12(1:3)
       planes(3,4)    = offset12
       planes(0:2,5) = -vNorm13(1:3)
       planes(3,5)    = offset13
       planes(0:2,6) = -vNorm14(1:3)
       planes(3,6)    = offset14

       CALL volume_lasserre_file( Planes(0,0), aj1 )
       
! Node 2

       planes(0:2,4) = vNorm12(1:3)
       planes(3,4)   = -offset12
       planes(0:2,5) = -vNorm23(1:3)
       planes(3,5)   = offset23
       planes(0:2,6) = -vNorm24(1:3)
       planes(3,6)   = offset24
       
       CALL volume_lasserre_file( Planes(0,0), aj2 )

! Node 3
       
       planes(0:2,4) = -vNorm34(1:3)
       planes(3,4)   = offset34
       planes(0:2,5) = vNorm13(1:3)
       planes(3,5)   = -offset13
       planes(0:2,6) = vNorm23(1:3)
       planes(3,6)   = -offset23

       CALL volume_lasserre_file( Planes(0,0), aj3)

! Node 4

       planes(0:2,4) = vNorm14(1:3)
       planes(3,4)   = -offset14
       planes(0:2,5) = vNorm24(1:3)
       planes(3,5)   = -offset24
       planes(0:2,6) = vNorm34(1:3)
       planes(3,6)   = -offset34

       CALL volume_lasserre_file( Planes(0,0), aj4)

       IF((aj1+aj2+aj3+aj4)/volEl.GT.1.0005.OR.&
            (aj1+aj2+aj3+aj4)/volEl.LT.0.995)THEN
          PRINT*,'ERROR: Node Volume Calculation Incorrect'
          PRINT*,(aj1+aj2+aj3+aj4)/volEl
          STOP
       ENDIF

       alphaR(1) = Aj1/VolEl
       alphaR(2) = Aj2/VolEl
       alphaR(3) = Aj3/VolEl
       alphaR(4) = Aj4/VolEl
! Option 1

    ELSE IF(NdMassLump.EQ.0)THEN

       alphaR(1) = 0.25d0
       alphaR(2) = 0.25d0
       alphaR(3) = 0.25d0
       alphaR(4) = 0.25d0 

     ENDIF

     RETURN
  END SUBROUTINE VolRatio

