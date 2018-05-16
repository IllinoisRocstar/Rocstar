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
!-----------------------------------------------------------------
!     CODE FOR 3-D ALE, ELASTODYNAMIC, SECOND-ORDER SOLID ELEMENT
!
!     written by Changyu Hwang  Nov. 20, 2001
!
!-----------------------------------------------------------------
SUBROUTINE v3d10_r_bar(d_bar,R_bar, &
     numnp,numlstet,lmlstet,meshcoor,nstart,nend) 

  IMPLICIT NONE

!---  for 10-node tetrahedron element !

  INTEGER numnp                 ! number of nodes
  REAL*8 d_bar(3*numnp)         ! nodal mesh displacements
  REAL*8 R_bar(3*numnp)         ! internal force
  REAL*8 meshcoor(3,numnp)
  INTEGER numlstet              ! number of LSTs
  INTEGER lmlstet(10,numlstet)   ! TET conectivity table
  INTEGER ielem , nx, ny, nz


  INTEGER ndof, nnode, ndim,  nmat, nsurf, nload

!--   for 10-node tetrahedron element !
  REAL*8 meshpos(3,10), meshvel(3,10), meshacc(3,10), &
       sload(4,1,3), bforce(3), elmass(30,30), elstiff(30,30), &
       eldamp(30,30), elload(30), elstiff2(30,30)


  INTEGER ndim_v, ndim_s, nintk_v, nintk_s,knode,  &
       knodr, knodc,i,j,rindx,cindx,k,mm,ll, &
       fp(2,6),iface,igauss, Grindx,Grindy,Grindz

!---  new variables for tetrahedron
  PARAMETER(ndim_v = 3, ndim_s = 2, nintk_v = 4,nintk_s = 3) ! 10-noded tetrahedron


!---  new implementation for 10-node tetrahedron

  DOUBLE PRECISION shfn(10,10), shdx(3,10,10), dv(10), dvsum, &
       delta(3,3), mvgss(3), magss(3), mcgss(3,3), two, &
       wsp1,wsparr(3),trmcgss,wsp2,wsp3,wsp4,one,lambda, &
       mu,third,wsp5,vxi(10),veta(10),vzeta(10),wt_lstet, &
       sxi(4,6),seta(4,6),szeta(4,6),shdxi(3,10,10), &
       wsparr1(3),shdx_av(3,10),zero, half 
  
  DOUBLE PRECISION disp(30), velo(30), Rimsi(30), Dimsi(30)
  DOUBLE PRECISION multiplier , ajacin(3,3), surfnormal(3,4)
  DOUBLE PRECISION factor 
  
  INTEGER :: nstart, nend

  PARAMETER(zero = 0.d0, one = 1.d0, two = 2.d0)
  DATA delta/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/
      
!--   input variables
  ndim  = 3
  nnode = 10 
  ndof  = ndim * nnode
  nsurf = 4

!     because the volume of the parent tetrahedron equals 1/6, the proper form
!     for a quadrature rule in this case is 
!     Integral ( f ) dV =  multiplier * Summation (fi * Ji * Wi)
!     where Ji:determinant of jacobian  Wi:weight  at each integration point i

  multiplier = 1.d0/6.d0	
  wt_lstet = 1.d0/4.d0

!
!--   Gauss integration points for tetrahedron
!
!     reference: Finite element procedures in engineering analysis 
!                by K. Bathe,  Figure 5.27 pp 232 
!     gauss point 1   alpha beta  beta  beta   1/4
!     gauss point 2   beta  alpha beta  beta   1/4
!     gauss point 3   beta  beta  alpha beta   1/4
!     gauss point 4   beta  beta  beta  alpha  1/4
!     alpha = 0.58541020d0
!     beta  = 0.13819660d0

!---  gauss point 1
  vxi(  1) = 0.13819660d0
  veta( 1) = 0.13819660d0
  vzeta(1) = 0.13819660d0

!---  gauss point 2
  vxi(  2) = 0.58541020d0
  veta( 2) = 0.13819660d0
  vzeta(2) = 0.13819660d0

!---  gauss point 3
  vxi(  3) = 0.13819660d0
  veta( 3) = 0.58541020d0
  vzeta(3) = 0.13819660d0

!---  gauss point 4
  vxi(  4) = 0.13819660d0
  veta( 4) = 0.13819660d0
  vzeta(4) = 0.58541020d0

!--   loop for elements

  DO ielem = nstart, nend
     
     DO i=1, nnode
        DO j=1, nnode
           elstiff(i,j) = 0.d0
        ENDDO
     ENDDO
     
     DO i=1, nnode
        nx = lmlstet(i,ielem)*3-2
        ny = lmlstet(i,ielem)*3-1
        nz = lmlstet(i,ielem)*3
        meshpos(1,i)=meshcoor(1,lmlstet(i,ielem))
        meshpos(2,i)=meshcoor(2,lmlstet(i,ielem))
        meshpos(3,i)=meshcoor(3,lmlstet(i,ielem))
        disp(i*3-2) =d_bar(nx)
        disp(i*3-1) =d_bar(ny)
        disp(i*3  ) =d_bar(nz)
     END DO

!--      calculate shape function and its derivatives
     CALL shcalc_3d10(shfn, shdx, dv, dvsum, shdx_av, meshpos, &
          ndim,nnode,nintk_v,vxi,veta,vzeta,shdxi,ajacin)

!        fill in element matrices for contributions from
!        volume integration points.

     factor=multiplier*wt_lstet
     DO igauss = 1, nintk_v
        DO knodr = 1,nnode
           DO knodc = 1,nnode
              elstiff(knodr,knodc) = elstiff(knodr,knodc) + &
                   ( shdx(1,knodr,igauss)*shdx(1,knodc,igauss)+ &
                   shdx(2,knodr,igauss)*shdx(2,knodc,igauss)+ &
                   shdx(3,knodr,igauss)*shdx(3,knodc,igauss) &
                   )*dv(igauss)*factor
           END DO
        END DO
     ENDDO ! for igauss = nintk_v

     DO i=1,nnode
        nx = lmlstet(i,ielem)*3-2
        ny = lmlstet(i,ielem)*3-1
        nz = lmlstet(i,ielem)*3
        DO j=1,nnode
           R_bar(nx)=R_bar(nx)+ elstiff(i,j)*disp(j*3-2)
           R_bar(ny)=R_bar(ny)+ elstiff(i,j)*disp(j*3-1)
           R_bar(nz)=R_bar(nz)+ elstiff(i,j)*disp(j*3-0)
        END DO
     END DO
     
  ENDDO	! for ielem=1,numlstet

  RETURN
END SUBROUTINE v3d10_r_bar

