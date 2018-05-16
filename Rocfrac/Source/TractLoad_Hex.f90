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
SUBROUTINE TractLoad_Hex(R_ex,numnp, &
     InterfaceElemTract, &
     InterfaceNumElems, InterfaceNumNodes, &
     InterfaceElemConn, &
     MapNode,LwrBnd,UppBnd,coor)

  IMPLICIT NONE

  INTEGER :: LwrBnd,UppBnd

  INTEGER :: InterfaceNumElems,InterfaceNumNodes
  INTEGER, DIMENSION(1:UppBnd,1:InterfaceNumElems) :: InterfaceElemConn
  INTEGER, DIMENSION(1:InterfaceNumNodes) :: MapNode
  REAL*8,  DIMENSION(1:InterfaceNumElems) :: InterfaceElemTract
  
  INTEGER :: numnp          ! number of nodes
  
  REAL*8, DIMENSION(1:3*numnp) :: R_ex ! external force
  
  INTEGER :: nx, ny, nz
  INTEGER :: i,j,k
  
  INTEGER,DIMENSION(1:4) :: Conn2D
  real*8,DIMENSION(1:3) :: dlvec
  real*8 :: pmag, MagNormDef, Area

! mesh coordinates
  REAL*8, DIMENSION(1:3,1:numnp) :: coor

  

  REAL*8, DIMENSION(1:3) :: veca, vecb, outwd_surface_normal, &
       half_normal_veca, half_normal_vecb
  

  DO i = 1, InterfaceNumElems
     
!  Nodes of the Quad ( Vertex Nodes ) 

     Conn2D(1:4) = InterfaceElemConn(1:4,i)

!  Uses the fact that the norm of the cross product vector
!  is the area of the parallelogram they form.  The triangle they
!  form has half that area.

     veca(1:3) =  coor(1:3,MapNode(Conn2D(2)))- coor(1:3,MapNode(Conn2D(1)))
     vecb(1:3) =  coor(1:3,MapNode(Conn2D(4)))- coor(1:3,MapNode(Conn2D(1)))
        
     half_normal_veca(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
     half_normal_veca(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
     half_normal_veca(3) = veca(1)*vecb(2) - vecb(1)*veca(2)
        
     veca(1:3) = coor(1:3,MapNode(Conn2D(4))) - coor(1:3,MapNode(Conn2D(3)))
     vecb(1:3) = coor(1:3,MapNode(Conn2D(2))) - coor(1:3,MapNode(Conn2D(3)))
        
     half_normal_vecb(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
     half_normal_vecb(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
     half_normal_vecb(3) = veca(1)*vecb(2) - vecb(1)*veca(2)
        
     outwd_surface_normal = half_normal_veca + half_normal_vecb 

!     MagNormDef = SQRT( outwd_surface_normal(1)**2 + outwd_surface_normal(2)**2 + outwd_surface_normal(3)**2 )

!     Area = 0.5*MagNormDef

     pmag = InterfaceElemTract(i)*0.125d0  ! /4 * /2
        
     dlvec(1:3) = - pmag * outwd_surface_normal(1:3)

! exactly the same results
!     print*,'**', - InterfaceElemTract(i) * Area * 0.25 * outwd_surface_normal(1:3)/MagNormDef
!     print*,'***', dlvec(1:3)


!  Assembly into global external load vector R_ex

     Conn2D(1:4) = InterfaceElemConn(LwrBnd:UppBnd,i)

     DO j = 1, 4
        nz = MapNode(Conn2D(j))*3 ! volume node
        nx = nz - 2
        ny = nz - 1
        R_ex(nx) = R_ex(nx) + dlvec(1)
        R_ex(ny) = R_ex(ny) + dlvec(2)
        R_ex(nz) = R_ex(nz) + dlvec(3)
     ENDDO
  ENDDO
  
  RETURN
END SUBROUTINE TractLoad_Hex

