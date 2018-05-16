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
SUBROUTINE cal_shdx(shdx,meshpos,ndim,nnode,nintk,surf_elem, &
     iface)

!!****f* Rocfrac/Rocfrac/Source/cal_shdx.f90
!!
!!  NAME
!!     cal_shdx
!!
!!  FUNCTION
!!     computes shape fn. values and parametric derivatives at Gauss pts.
!!
!!  INPUTS
!!     meshpos -- mesh coordinates
!!     ndim    -- dimension of problem
!!    nnode    -- number of nodes
!!    nintk    -- number of integration points
!!   surf_eleme --  connectivity of surf_elem
!!     iface  -- face number of tetrahedral
!!
!!  OUTPUT
!!     shdx -- shape fn. derivatives
!!
!!****

  IMPLICIT NONE

  
  INTEGER ndim, nnode, nintk
  DOUBLE PRECISION shdx(ndim,nnode,nintk), &
       meshpos(ndim,nnode)
  DOUBLE PRECISION xi, eta
  DOUBLE PRECISION eighth,one,zero
  INTEGER igauss, knode, i,j, iface, surf_elem(3,4)
  
  DOUBLE PRECISION SHAPE(3), dshape(2,3), jacobian(2,3)
  INTEGER n1,n2,n3 
  

  one = 1.d0
  eighth = 1.d0/8.d0
  zero = 0.d0

!     fill in shape fn. values and parametric derivatives at Gauss pts.

  igauss = 1
  xi  = 1.d0/3.d0
  eta = 1.d0/3.d0
  
  SHAPE(1) = xi
  SHAPE(2) = eta
  SHAPE(3) = one - xi -eta
  
  dshape(1,1) =  one
  dshape(2,1) =  zero
  dshape(1,2) =  zero 
  dshape(2,2) =  one
  dshape(1,3) = -one
  dshape(2,3) = -one 

  n1 = surf_elem(1,iface) 
  n2 = surf_elem(2,iface) 
  n3 = surf_elem(3,iface) 
  
 
  jacobian(1,1) = 1.d0/(meshpos(1,n1)-meshpos(1,n3))
  jacobian(1,2) = 1.d0/(meshpos(2,n1)-meshpos(2,n3))
  jacobian(1,3) = 1.d0/(meshpos(3,n1)-meshpos(3,n3))
  jacobian(2,1) = 1.d0/(meshpos(1,n2)-meshpos(1,n3))
  jacobian(2,2) = 1.d0/(meshpos(2,n2)-meshpos(2,n3))
  jacobian(2,3) = 1.d0/(meshpos(3,n2)-meshpos(3,n3))

!     do i=1,2
!       write(*,'(i5,3e15.7)') i,(jacobian(i,j),j=1,3)
!     end do
!     evaluate shape fn. derivatives
  DO knode = 1, nnode
     DO i = 1,ndim
        shdx(i,knode,igauss) = zero
     END DO
  END DO

  DO knode = 1,3
     DO i = 1,ndim
        shdx(i,surf_elem(knode,iface),igauss) = zero
        DO j = 1,2
           shdx(i,surf_elem(knode,iface),igauss) = &
                shdx(i,surf_elem(knode,iface),igauss) + &
                dshape(j,knode)*jacobian(j,i)
        END DO
     END DO
  END DO

  RETURN
END SUBROUTINE cal_shdx

