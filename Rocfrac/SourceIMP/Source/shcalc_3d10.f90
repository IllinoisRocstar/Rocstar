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
!!****f* Rocfrac/Rocfrac/Source/shcalc_3d10.f90
!!
!!  NAME
!!     shcalc_3d10
!!
!!  FUNCTION
!!     EVALUALTE SHAPE FUNCTION AND ITS DERIVATIVES FOR 10-NODE TET.
!!
!!  AUTHOR
!!     written by Changyu Hwang  Nov. 20, 2001
!!****

SUBROUTINE shcalc_3d10(shfn, shdx, dv, dvsum, shdx_av, meshpos, &
     ndim,nnode,nintk,xi,eta,zeta,shdxi,wsp2)

  IMPLICIT NONE

  INTEGER ndim, nnode, nintk

  REAL*8 :: shfn(nnode,nintk), shdx(ndim,nnode,nintk), &
       dv(nintk), dvsum, shdx_av(ndim,nnode), &
       meshpos(ndim,nnode),xi(nintk),eta(nintk), &
       zeta(nintk),shdxi(ndim,nnode,nintk)

  REAL*8 ::  wsp1(3,3),wsp2(3,3),eighth,one,zero,four, two

  REAL*8 ::  sh(10), d_sh(3,10) , r, s, t, u, wi
  INTEGER igauss, knode, i,j


  one = 1.d0
  two = 2.d0
  four = 4.d0
  eighth = 1.d0/8.d0
  zero = 0.d0
  dvsum = 0.d0
  wi = 1.d0/4.d0

!     fill in shape fn. values and parametric derivatives at Gauss pts.
  
  DO igauss = 1, nintk
     
     r =  xi(igauss)
     s = eta(igauss)
     t = zeta(igauss)
     u = 1.d0 - r - s - t

     sh(1) = u * ( two * u - one ) 
     sh(2) = r * ( two * r - one )
     sh(3) = s * ( two * s - one )
     sh(4) = t * ( two * t - one )
     
     sh(5 )= four * u * r
     sh(6 )= four * r * s
     sh(7 )= four * s * u 
     sh(8 )= four * u * t
     sh(9 )= four * r * t 
     sh(10)= four * s * t 
     
     d_sh(1,1) = -four * u + one
     d_sh(2,1) = -four * u + one
     d_sh(3,1) = -four * u + one
     
     d_sh(1,2) = four * r - one
     d_sh(2,2) = zero
     d_sh(3,2) = zero
     
     d_sh(1,3) = zero 
     d_sh(2,3) = four * s - one 
     d_sh(3,3) = zero 
     
     d_sh(1,4) = zero 
     d_sh(2,4) = zero 
     d_sh(3,4) = four * t - one 
     
     d_sh(1,5) = - four * r + four * (one-r-s-t)
     d_sh(2,5) = - four * r 
     d_sh(3,5) = - four * r
     
     d_sh(1,6) = four * s 
     d_sh(2,6) = four * r 
     d_sh(3,6) = zero 
     
     d_sh(1,7) = - four * s
     d_sh(2,7) = - four * s + four * (one-r-s-t)
     d_sh(3,7) = - four * s
     
     d_sh(1,8)= - four * t
     d_sh(2,8)= - four * t
     d_sh(3,8)= - four * t + four * (one-r-s-t)
     
     d_sh(1,9) = four * t 
     d_sh(2,9) = zero 
     d_sh(3,9) = four * r 
     
     d_sh(1,10) = zero 
     d_sh(2,10) = four * t 
     d_sh(3,10) = four * s 
     
     DO j=1,nnode
        shfn(j,igauss) = sh(j)
     END DO
     
     DO i=1,ndim
        DO j=1,nnode
           shdxi(i,j,igauss) = d_sh(i,j)
        END DO
     END DO

!     evaluate and store common quantities for all Gauss pts.
!     and shape fn. derivatives w.r.t. spatial coordinates.
     
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

     dvsum = dvsum + dv(igauss) * wi

!        evaluate shape fn. derivatives

     DO knode = 1,nnode
        DO i = 1,ndim
           shdx(i,knode,igauss) = zero
           DO j = 1,ndim
              shdx(i,knode,igauss) = shdx(i,knode,igauss) +  &
                   shdxi(j,knode,igauss)*wsp2(j,i)
           END DO
        END DO
     END DO
     

  ENDDO
  
  RETURN
END SUBROUTINE shcalc_3d10

