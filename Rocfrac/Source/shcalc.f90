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
SUBROUTINE shcalc(shfn, shdx, dv, dvsum, shdx_av, meshpos, &
     ndim,nnode,nintk,xi,eta,zeta,shdxi,wsp2)

!!****f* Rocfrac/Rocfrac/Source/shcalc.f90
!!
!!  NAME
!!     shcalc
!!
!!  FUNCTION
!!     EVALUALTE SHAPE FUNCTION AND ITS DERIVATIVES FOR 4-NODE TET.
!!
!!****

  IMPLICIT NONE
  
  INTEGER ndim, nnode, nintk

  DOUBLE PRECISION shfn(nnode,nintk), shdx(ndim,nnode,nintk), &
       dv(nintk), dvsum, shdx_av(ndim,nnode), &
       meshpos(ndim,nnode),xi(nintk),eta(nintk), &
       zeta(nintk),shdxi(ndim,nnode,nintk)
  
  DOUBLE PRECISION wsp1(3,3),wsp2(3,3),eighth,one,zero
  
  INTEGER igauss, knode, i,j


  one = 1.d0
  eighth = 1.d0/8.d0
  zero = 0.d0

!     fill in shape fn. values and parametric derivatives at Gauss pts.

  igauss = 1
  
  shfn(1,igauss) = zeta(igauss)
  shfn(2,igauss) = xi(igauss) 
  shfn(3,igauss) = one-xi(igauss)-eta(igauss)-zeta(igauss)  
  shfn(4,igauss) = eta(igauss)
  
  shdxi(1,1,igauss) = zero
  shdxi(2,1,igauss) = zero 
  shdxi(3,1,igauss) = one
  
  shdxi(1,2,igauss) = one
  shdxi(2,2,igauss) = zero
  shdxi(3,2,igauss) = zero
  
  shdxi(1,3,igauss) =-one
  shdxi(2,3,igauss) =-one
  shdxi(3,3,igauss) =-one
  
  shdxi(1,4,igauss) = zero
  shdxi(2,4,igauss) = one
  shdxi(3,4,igauss) = zero
  

!    evaluate and store common quantities for all Gauss pts.
!    and shape fn. derivatives w.r.t. spatial coordinates.

  DO i = 1,ndim
     DO j = 1,ndim
        wsp1(i,j) = zero
        DO knode = 1,nnode
           wsp1(i,j) = wsp1(i,j) + meshpos(i,knode)* &
                shdxi(j,knode,igauss)
        END DO
     END DO
  END DO

  CALL ainv(wsp1,wsp2,dv(igauss),ndim)

!        evaluate shape fn. derivatives
  
  DO knode = 1,nnode
     DO i = 1,ndim
        shdx(i,knode,igauss) = zero
        DO j = 1,ndim
           shdx(i,knode,igauss) = shdx(i,knode,igauss) + &
                shdxi(j,knode,igauss)*wsp2(j,i)
        END DO
     END DO
  END DO
  
  RETURN
END SUBROUTINE shcalc

