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
SUBROUTINE TractLoad(R_ex,numnp, &
     InterfaceSFElemTract, &
     InterfaceSFNumElems, InterfaceSFNumNodes, &
     InterfaceSFElemConn, &
     MapNodeSF,LwrBnd,UppBnd,coor)

  IMPLICIT NONE


! 1-3 for 4 node tet (i.e. 3 node triangles)
! 4-6 for 10 node tet (i.e. 6 node triangles, mid-side nodes get traction)
  INTEGER :: LwrBnd,UppBnd

  INTEGER :: InterfaceSFNumElems,InterfaceSFNumNodes
  INTEGER, DIMENSION(1:UppBnd,1:InterfaceSFNumElems) :: InterfaceSFElemConn
  INTEGER, DIMENSION(1:InterfaceSFNumNodes) :: MapNodeSF
  REAL*8,  DIMENSION(1:3, 1:InterfaceSFNumElems) :: InterfaceSFElemTract
  
  INTEGER :: numnp          ! number of nodes
  
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
  

  REAL*8 :: xNorm, Area, at
  
  at = 0.d0

  DO i = 1, InterfaceSFNumElems
     
!  Nodes of the Triangle ( Vertex Nodes ) 

     TriConn(1:3) = InterfaceSFElemConn(1:3,i)

!  Uses the fact that the norm of the cross product vector
!  is the area of the parallelogram they form.  The triangle they
!  form has half that area.

     x0p = coor(1,MapNodeSF(TriConn(1)))
     y0p = coor(2,MapNodeSF(TriConn(1)))
     z0p = coor(3,MapNodeSF(TriConn(1)))
     
     x1p = coor(1,MapNodeSF(TriConn(2)))
     y1p = coor(2,MapNodeSF(TriConn(2)))
     z1p = coor(3,MapNodeSF(TriConn(2)))
     
     x2p = coor(1,MapNodeSF(TriConn(3)))
     y2p = coor(2,MapNodeSF(TriConn(3)))
     z2p = coor(3,MapNodeSF(TriConn(3)))

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
               
     Area = 0.5d0*xnorm

!  Assign the load weighted equally among the 3 nodes 

     UnifTract(1:3) = InterfaceSFElemTract(1:3,i)*Area/3.d0

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
  
  RETURN
END SUBROUTINE TractLoad

