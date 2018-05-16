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
SUBROUTINE HeatFluxLoad(Rnet,NumNP, &
     InterfaceNumElems, InterfaceNumNodes, &
     InterfaceElemConn, &
     MapNode,LwrBnd,UppBnd, Coor, MapSFElVolEl,&
     ElConnVol,iElType,NumElVol,BCqflux)

  IMPLICIT NONE

! 1-3 for 4 node tet (i.e. 3 node triangles)
! 4-6 for 10 node tet (i.e. 6 node triangles, mid-side nodes get traction)

! Computes the traction due to a dummy implied traction

! Applied to the deformed Configuration

  INTEGER :: LwrBnd,UppBnd

  INTEGER :: iElType,NumElVol
  INTEGER, DIMENSION(1:iElType,1:NumElVol) :: ElConnVol

  INTEGER :: InterfaceNumElems,InterfaceNumNodes
  INTEGER, DIMENSION(1:UppBnd,1:InterfaceNumElems) :: InterfaceElemConn
  
  INTEGER :: NumNP          ! number of nodes

  INTEGER, DIMENSION(1:InterfaceNumNodes) :: MapNode
  
  REAL*8, DIMENSION(1:numnp) :: Rnet ! external force
  
  INTEGER :: i
  
  INTEGER,DIMENSION(1:3) :: TriConn

! Array of imposed fluxes on fluid/solid interface
  REAL*8,DIMENSION(1:InterfaceNumElems) :: BCqFlux
  
  INTEGER, DIMENSION(1:NumElVol) :: MapSFElVolEl

! mesh coordinates
  REAL*8, DIMENSION(1:3,1:numnp) :: Coor
!  REAL*8, DIMENSION(1:3,1:numnp) :: disp

  REAL*8 :: x0p,y0p,z0p
  REAL*8 :: x1p,y1p,z1p
  REAL*8 :: x2p,y2p,z2p
  REAL*8 :: x3p,y3p,z3p

! Components sides vector

  REAL*8 :: x1x0, y1y0, z1z0
  REAL*8 :: x2x0, y2y0, z2z0
  
  REAL*8 :: xNorm, Area
  

  
  DO i = 1, InterfaceNumElems
     
!  Nodes of the Triangle ( Vertex Nodes ) 

     TriConn(1:3) = InterfaceElemConn(1:3,i)

!  Uses the fact that the norm of the cross product vector
!  is the area of the parallelogram they form.  The triangle they
!  form has half that area.

     x0p = Coor(1,MapNode(TriConn(1))) !+ Disp(3*MapNode(TriConn(1))-2)
     y0p = Coor(2,MapNode(TriConn(1))) !+ Disp(3*MapNode(TriConn(1))-1)
     z0p = Coor(3,MapNode(TriConn(1))) !+ Disp(3*MapNode(TriConn(1)))
     
     x1p = Coor(1,MapNode(TriConn(2))) !+ Disp(3*MapNode(TriConn(2))-2)
     y1p = Coor(2,MapNode(TriConn(2))) !+ Disp(3*MapNode(TriConn(2))-1)
     z1p = Coor(3,MapNode(TriConn(2))) !+ Disp(3*MapNode(TriConn(2)))
     
     x2p = Coor(1,MapNode(TriConn(3))) !+ Disp(3*MapNode(TriConn(3))-2)
     y2p = Coor(2,MapNode(TriConn(3))) !+ Disp(3*MapNode(TriConn(3))-1)
     z2p = Coor(3,MapNode(TriConn(3))) !+ Disp(3*MapNode(TriConn(3)))
	
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
         
!  Computes the Euclidean norm of a vector in 3D
         
     xNorm = SQRT( x3p*x3p + y3p*y3p + z3p*z3p )
!
!     Area of the triangle is half the area of the parallelogram 
!     determined by the cross product of the vectors P12 and P13
       
     Area = 0.5d0*xNorm

     
!  Assembly into global load vector Rnet

     TriConn(1:3) = InterfaceElemConn(LwrBnd:UppBnd,i)

     Rnet(MapNode(TriConn(1))) = Rnet(MapNode(TriConn(1))) + BCqFlux(i)*Area/3.D0
     Rnet(MapNode(TriConn(2))) = Rnet(MapNode(TriConn(2))) + BCqFlux(i)*Area/3.D0
     Rnet(MapNode(TriConn(3))) = Rnet(MapNode(TriConn(3))) + BCqFlux(i)*Area/3.D0

  ENDDO
  
  RETURN
  
END SUBROUTINE HeatFluxLoad

