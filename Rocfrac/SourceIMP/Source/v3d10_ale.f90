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
!     ALE formulation for 10-node tetrahedron element
!
!     written by Changyu Hwang  Nov. 20, 2001
!
!     reference: Finite element procedures in engineering analysis 
!                by K. Bathe,  Figure 5.27 pp 232 
!-----------------------------------------------------------------------
!
!
SUBROUTINE v3d10_ale(v_bar,a_bar,d,vhalf,Rnet, &
     E,xnu,rho,numnp,numat_vol, &
     numlstet,matlstet,lmlstet,meshcoor, &
     nstart,nend)

  IMPLICIT NONE

!---  for tetrahedron element !

  INTEGER numnp                 ! number of nodes
  INTEGER numat_vol             ! number of materials
  REAL*8 v_bar(3*numnp)         ! nodal mesh velocities
  REAL*8 a_bar(3*numnp)         ! nodal mesh accelerations
  REAL*8 d(3*numnp)             ! nodal displacements
  REAL*8 vhalf(3*numnp)         ! nodal velocities
  REAL*8 Rnet(3*numnp)          ! internal force
  REAL*8 E(numat_vol)       ! young's moduli
  REAL*8 xnu(numat_vol) ! Poisson's ratios
  REAL*8 rho(numat_vol)                     ! densities
  REAL*8 meshcoor(3,numnp)
  INTEGER numlstet              ! number of LSTs
  INTEGER matlstet(numlstet)    ! list of all materials
  INTEGER lmlstet(10,numlstet)   ! TET conectivity table
  INTEGER ielem , nx, ny, nz

!     integer ndof, nnode, ndim,  nmat, nsurf, nload, ldflag(nload)
  INTEGER ndof, nnode, ndim,  nmat, nsurf, nload

!--   for 10-node tetrahedron element !
  REAL*8 meshpos(3,10), meshvel(3,10), meshacc(3,10), mat(3), &
       sload(4,1,3), bforce(3), elmass(30,30), elstiff(30,30), &
       eldamp(30,30), elload(30), elstiff2(30,30)

!     INPUT:

!     ndof  -- # of element degrees of freedom
!     nnode -- # of nodes on element
!     ndim  -- # of spatial dimensions
!     nmat  -- # of material parameters (rho, E, nu)
!     nsurf -- # of surfaces on element
!     nload -- # of load types for element (2 - body force, surface traction,
!                        both distributions assumed uniform)
!     meshpos -- coordinates of position w.r.t. to global basis of nodes
!     meshvel -- components of mesh velocity w.r.t to global basis at nodes
!     meshacc --      ''    ''  ''  accelern.  ''  ''   ''     ''  ''   ''
!     mat     -- array containing material parameters (e.g. rho, E, nu)
!     ldflag  -- array indicating if there are body forces and surface loads
!               (in that order) to be taken care of in this element
!               (values:1/0)
!     sload   -- array containing surface loads
!              (face#, pos  or neg if load present or absent, load vec comp)
!     bforce  -- body force components

!     OUTPUT:

!     elmass  -- element mass matrix
!     elstiff -- element stiffness matrix
!     eldamp  -- element damping matrix
!     elload  -- element load vector

  INTEGER ndim_v, ndim_s, nintk_v, nintk_s,knode, &
       knodr, knodc,i,j,rindx,cindx,k,mm,ll, &
       fp(2,6),iface,igauss, Grindx,Grindy,Grindz
!
!---  variables for integration points 
!     parameter(ndim_v = 3, ndim_s = 2, nintk_v = 8,nintk_s = 4) ! 8-noded brick
!     parameter(ndim_v = 3, ndim_s = 2, nintk_v = 4,nintk_s = 1) ! 4-noded tetrahedron
  PARAMETER(ndim_v = 3, ndim_s = 2, nintk_v = 4,nintk_s = 3) ! 10-noded tetrahedron

  DOUBLE PRECISION shfn(10,4), shdx(3,10,4), dv(10), dvsum, &
       delta(3,3), mvgss(3), magss(3), mcgss(3,3), two, &
       wsp1,wsparr(3),trmcgss,wsp2,wsp3,wsp4,one,lambda, &
       mu,third,wsp5,vxi(10),veta(10),vzeta(10),wvar, &
       sxi(4,6),seta(4,6),szeta(4,6),shdxi(3,10,4), &
       wsparr1(3),shdx_av(3,10),zero, half 
  DOUBLE PRECISION disp(30), velo(30), Rimsi(30), Dimsi(30)
  DOUBLE PRECISION multiplier , ajacin(3,3), surfnormal(3,4)
  DOUBLE PRECISION wt_lstet, dA, sum,factor, wt_cstet
  DOUBLE PRECISION sto_wsp1(10),sto_wsp2(10)
  DOUBLE PRECISION sto_wsp3(10),sto_wsp4(10)
  INTEGER vect1(2,10), vect2(2,10), surf_elem(3,4), k2
  REAL*8 R_in(3*numnp)          ! internal force
  REAL*8 R_dp(3*numnp)          ! damping force

  INTEGER :: nstart, nend

  PARAMETER(zero = 0.d0, one = 1.d0, two = 2.d0)

  DATA delta/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/


  ndim  = 3
  nnode = 10 
  ndof  = ndim * nnode
  nsurf = 4
  multiplier = 1.d0/6.d0
  wt_lstet = 1.d0/4.d0
  wvar = 1.d0/3.d0
      

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

!--   three points on each surface: outward normal direction
!--   face 1  : node 1 3 2 
!--   face 2  : node 1 2 4
!--   face 3  : node 2 3 4
!--   face 4  : node 1 4 3

  wvar = 1.d0/dsqrt(3.d0)
  third = 1.d0/3.d0
  half = 0.5d0

!--   face 1  : node 1 3 2   (szeta=t=0)
  sxi(  1,1) = half 
  seta( 1,1) = zero 
  szeta(1,1) = zero 
  sxi(  2,1) = zero 
  seta( 2,1) = half 
  szeta(2,1) = zero 
  sxi(  3,1) = half 
  seta( 3,1) = half 
  szeta(3,1) = zero 
  
!--   face 2  : node 1 2 4   (seta=s=0)
  sxi(  1,2) = half 
  seta( 1,2) = zero 
  szeta(1,2) = zero 
  sxi(  2,2) = zero 
  seta( 2,2) = zero 
  szeta(2,2) = half 
  sxi(  3,2) = half 
  seta( 3,2) = zero 
  szeta(3,2) = half 

!--   face 3  : node 3 4 2
  sxi(  1,3) = half 
  seta( 1,3) = half 
  szeta(1,3) = zero 
  sxi(  2,3) = half 
  seta( 2,3) = zero 
  szeta(2,3) = half 
  sxi(  3,3) = zero 
  seta( 3,3) = half 
  szeta(3,3) = half 
  
!--   face 4  : node 1 4 3   (sxi=r=0)
  sxi(  1,4) = zero 
  seta( 1,4) = half 
  szeta(1,4) = zero 
  sxi(  2,4) = zero 
  seta( 2,4) = zero 
  szeta(2,4) = half 
  sxi(  3,4) = zero 
  seta( 3,4) = half 
  szeta(3,4) = half 

!--   surface normal vectors

  surf_elem(1,1)=1
  surf_elem(2,1)=3
  surf_elem(3,1)=2
  
  surf_elem(1,2)=1
  surf_elem(2,2)=2
  surf_elem(3,2)=4
  
  surf_elem(1,3)=3
  surf_elem(2,3)=4
  surf_elem(3,3)=2
  
  surf_elem(1,4)=1
  surf_elem(2,4)=4
  surf_elem(3,4)=3
  
!--   initialize R_in, R_damping

  DO i = 1, numnp*3 
     R_in(i) = 0.d0
     R_dp(i) = 0.d0
  ENDDO

!     loop for elements
  
  DO ielem = nstart, nend
     
     j = matlstet(ielem)
     mat(1)= rho(j)
     mat(2)= E(j)
     mat(3)= xnu(j)

!--      obtain mesh coordinate, mesh velocity, mesh accel.,
!--      displacement, velocity 

     DO i=1, nnode
        nx = lmlstet(i,ielem)*3-2
        ny = lmlstet(i,ielem)*3-1
        nz = lmlstet(i,ielem)*3
        meshpos(1,i)=meshcoor(1,lmlstet(i,ielem))
        meshpos(2,i)=meshcoor(2,lmlstet(i,ielem))
        meshpos(3,i)=meshcoor(3,lmlstet(i,ielem))
        meshvel(1,i)=v_bar(nx)
        meshvel(2,i)=v_bar(ny)
        meshvel(3,i)=v_bar(nz)
        meshacc(1,i)=a_bar(nx)
        meshacc(2,i)=a_bar(ny)
        meshacc(3,i)=a_bar(nz)
        disp(i*3-2) =d(nx)
        disp(i*3-1) =d(ny)
        disp(i*3  ) =d(nz)
        velo(i*3-2) =vhalf(nx)
        velo(i*3-1) =vhalf(ny)
        velo(i*3  ) =vhalf(nz)
     END DO


!--      evaluate shape function and its derivative at integration points

     CALL shcalc_3d10(shfn, shdx, dv, dvsum, shdx_av, meshpos, &
          ndim,nnode,nintk_v,vxi,veta,vzeta,shdxi,ajacin)


     DO igauss = 1, nintk_v


!          fill in element matrices for contributions from
!          volume integration points.
  
!           evaluate mesh-related quantitites at vol. int. pt.
        DO j = 1,ndim
           mvgss(j) = zero
           magss(j) = zero
           DO knode = 1,nnode
              mvgss(j) = mvgss(j) + meshvel(j,knode)* &
                   shfn(knode,igauss)
              magss(j) = magss(j) + meshacc(j,knode)* &
                   shfn(knode,igauss)
           END DO
        END DO
        
        mcgss(1,1) = zero
        mcgss(1,2) = zero
        mcgss(1,3) = zero
        mcgss(2,1) = zero
        mcgss(2,2) = zero
        mcgss(2,3) = zero
        mcgss(3,1) = zero
        mcgss(3,2) = zero
        mcgss(3,3) = zero
        
        DO knode = 1,nnode
           mcgss(1,1) = mcgss(1,1) + meshvel(1,knode)* &
                                           shdx(1,knode,igauss)
           mcgss(2,1) = mcgss(2,1) + meshvel(2,knode)* &
                shdx(1,knode,igauss)
           mcgss(3,1) = mcgss(3,1) + meshvel(3,knode)* &
                shdx(1,knode,igauss)
           mcgss(1,2) = mcgss(1,2) + meshvel(1,knode)* &
                shdx(2,knode,igauss)
           mcgss(2,2) = mcgss(2,2) + meshvel(2,knode)* &
                shdx(2,knode,igauss)
           mcgss(3,2) = mcgss(3,2) + meshvel(3,knode)* &
                shdx(2,knode,igauss)
           mcgss(1,3) = mcgss(1,3) + meshvel(1,knode)* &
                shdx(3,knode,igauss)
           mcgss(2,3) = mcgss(2,3) + meshvel(2,knode)* &
                shdx(3,knode,igauss)
           mcgss(3,3) = mcgss(3,3) + meshvel(3,knode)* &
                shdx(3,knode,igauss)
        END DO

        trmcgss = mcgss(1,1) + mcgss(2,2) + mcgss(3,3)

        wsparr(1) = mcgss(1,1)*mvgss(1) + &
             mcgss(1,2)*mvgss(2) + &
             mcgss(1,3)*mvgss(3) 
        wsparr(2) = mcgss(2,1)*mvgss(1) + &
             mcgss(2,2)*mvgss(2) + &
             mcgss(2,3)*mvgss(3) 
        wsparr(3) = mcgss(3,1)*mvgss(1) + &
             mcgss(3,2)*mvgss(2) + &
             mcgss(3,3)*mvgss(3) 
        

        DO knodc = 1,nnode
           sto_wsp1(knodc) = &
                mvgss(1)*shdx(1,knodc,igauss) + &
                mvgss(2)*shdx(2,knodc,igauss) + &
                mvgss(3)*shdx(3,knodc,igauss)
           sto_wsp2(knodc) =  &
                wsparr(1)*shdx(1,knodc,igauss) + &
                wsparr(2)*shdx(2,knodc,igauss) + &
                wsparr(3)*shdx(3,knodc,igauss)
           sto_wsp3(knodc) =  &
                mvgss(1)*shdx(1,knodc,igauss) + &
                mvgss(2)*shdx(2,knodc,igauss) + &
                mvgss(3)*shdx(3,knodc,igauss)
           sto_wsp4(knodc) =  &
                magss(1)*shdx(1,knodc,igauss) + &
                magss(2)*shdx(2,knodc,igauss) + &
                magss(3)*shdx(3,knodc,igauss)
        END DO

        factor=mat(1)*dv(igauss)*multiplier*wt_lstet

        DO knodr = 1,nnode
           Grindx = (lmlstet(knodr,ielem)-1)*3 + 1
           Grindy = (lmlstet(knodr,ielem)-1)*3 + 2
           Grindz = (lmlstet(knodr,ielem)-1)*3 + 3
           
!                    damping 
!                    ----------------------------------
           
           R_dp(Grindx) = R_dp(Grindx) -  &
                factor*shfn(knodr,igauss)*two* &
                ( sto_wsp1(1)*velo((1-1)*ndim+1) + &
                sto_wsp1(2)*velo((2-1)*ndim+1) + &
                sto_wsp1(3)*velo((3-1)*ndim+1) + &
                sto_wsp1(4)*velo((4-1)*ndim+1) + &
                sto_wsp1(5)*velo((5-1)*ndim+1) + &
                sto_wsp1(6)*velo((6-1)*ndim+1) + &
                sto_wsp1(7)*velo((7-1)*ndim+1) + &
                sto_wsp1(8)*velo((8-1)*ndim+1) + &
                sto_wsp1(9)*velo((9-1)*ndim+1) + &
                sto_wsp1(10)*velo((10-1)*ndim+1) )
           
           R_dp(Grindy) = R_dp(Grindy) - & 
                factor*shfn(knodr,igauss)*two* &
                ( sto_wsp1(1)*velo((1-1)*ndim+2) + &
                sto_wsp1(2)*velo((2-1)*ndim+2) + &
                sto_wsp1(3)*velo((3-1)*ndim+2) + &
                sto_wsp1(4)*velo((4-1)*ndim+2) + &
                sto_wsp1(5)*velo((5-1)*ndim+2) + &
                sto_wsp1(6)*velo((6-1)*ndim+2) + &
                sto_wsp1(7)*velo((7-1)*ndim+2) + &
                sto_wsp1(8)*velo((8-1)*ndim+2) + &
                sto_wsp1(9)*velo((9-1)*ndim+2) + &
                sto_wsp1(10)*velo((10-1)*ndim+2) )
           
           R_dp(Grindz) = R_dp(Grindz) -  &
                factor*shfn(knodr,igauss)*two* &
                ( sto_wsp1(1)*velo((1-1)*ndim+3) + &
                sto_wsp1(2)*velo((2-1)*ndim+3) + &
                sto_wsp1(3)*velo((3-1)*ndim+3) + &
                sto_wsp1(4)*velo((4-1)*ndim+3) + &
                sto_wsp1(5)*velo((5-1)*ndim+3) + &
                sto_wsp1(6)*velo((6-1)*ndim+3) + &
                sto_wsp1(7)*velo((7-1)*ndim+3) + &
                sto_wsp1(8)*velo((8-1)*ndim+3) + &
                sto_wsp1(9)*velo((9-1)*ndim+3) + &
                sto_wsp1(10)*velo((10-1)*ndim+3) )

!                    stiffness (vol. integration terms)
!                    ----------------------------------

           R_in(Grindx) = R_in(Grindx)  &
                + factor*shfn(knodr,igauss) * ( &
                -trmcgss*sto_wsp1(1)*disp((1-1)*ndim+1) &
                -trmcgss*sto_wsp1(2)*disp((2-1)*ndim+1) &
                -trmcgss*sto_wsp1(3)*disp((3-1)*ndim+1) &
                -trmcgss*sto_wsp1(4)*disp((4-1)*ndim+1) &
                -trmcgss*sto_wsp1(5)*disp((5-1)*ndim+1) &
                -trmcgss*sto_wsp1(6)*disp((6-1)*ndim+1) &
                -trmcgss*sto_wsp1(7)*disp((7-1)*ndim+1) &
                -trmcgss*sto_wsp1(8)*disp((8-1)*ndim+1) &
                -trmcgss*sto_wsp1(9)*disp((9-1)*ndim+1) &
                -trmcgss*sto_wsp1(10)*disp((10-1)*ndim+1) &
                +        sto_wsp2(1)*disp((1-1)*ndim+1) &
                +        sto_wsp2(2)*disp((2-1)*ndim+1) &
                +        sto_wsp2(3)*disp((3-1)*ndim+1) &
                +        sto_wsp2(4)*disp((4-1)*ndim+1) &
                +        sto_wsp2(5)*disp((5-1)*ndim+1) &
                +        sto_wsp2(6)*disp((6-1)*ndim+1) &
                +        sto_wsp2(7)*disp((7-1)*ndim+1) &
                +        sto_wsp2(8)*disp((8-1)*ndim+1) &
                +        sto_wsp2(9)*disp((9-1)*ndim+1) &
                +        sto_wsp2(10)*disp((10-1)*ndim+1) &
                -        sto_wsp4(1)*disp((1-1)*ndim+1) &
                -        sto_wsp4(2)*disp((2-1)*ndim+1) &
                -        sto_wsp4(3)*disp((3-1)*ndim+1) &
                -        sto_wsp4(4)*disp((4-1)*ndim+1) &
                -        sto_wsp4(5)*disp((5-1)*ndim+1) &
                -        sto_wsp4(6)*disp((6-1)*ndim+1) &
                -        sto_wsp4(7)*disp((7-1)*ndim+1) &
                -        sto_wsp4(8)*disp((8-1)*ndim+1) &
                -        sto_wsp4(9)*disp((9-1)*ndim+1) &
                -        sto_wsp4(10)*disp((10-1)*ndim+1) ) &
                + factor*sto_wsp3(knodr) * (          &
                -        sto_wsp1(1)*disp((1-1)*ndim+1) &
                -        sto_wsp1(2)*disp((2-1)*ndim+1) &
                -        sto_wsp1(3)*disp((3-1)*ndim+1) &
                -        sto_wsp1(4)*disp((4-1)*ndim+1) &
                -        sto_wsp1(5)*disp((5-1)*ndim+1) &
                -        sto_wsp1(6)*disp((6-1)*ndim+1) &
                -        sto_wsp1(7)*disp((7-1)*ndim+1) &
                -        sto_wsp1(8)*disp((8-1)*ndim+1) &
                -        sto_wsp1(9)*disp((9-1)*ndim+1) &
                -        sto_wsp1(10)*disp((10-1)*ndim+1) )


           R_in(Grindy) = R_in(Grindy)  &
                + factor*shfn(knodr,igauss) * ( &
                -trmcgss*sto_wsp1(1)*disp((1-1)*ndim+2) &
                -trmcgss*sto_wsp1(2)*disp((2-1)*ndim+2) &
                -trmcgss*sto_wsp1(3)*disp((3-1)*ndim+2) &
                -trmcgss*sto_wsp1(4)*disp((4-1)*ndim+2) &
                -trmcgss*sto_wsp1(5)*disp((5-1)*ndim+2) &
                -trmcgss*sto_wsp1(6)*disp((6-1)*ndim+2) &
                -trmcgss*sto_wsp1(7)*disp((7-1)*ndim+2) &
                -trmcgss*sto_wsp1(8)*disp((8-1)*ndim+2) &
                -trmcgss*sto_wsp1(9)*disp((9-1)*ndim+2) &
                -trmcgss*sto_wsp1(10)*disp((10-1)*ndim+2) &
                +        sto_wsp2(1)*disp((1-1)*ndim+2) &
                +        sto_wsp2(2)*disp((2-1)*ndim+2) &
                +        sto_wsp2(3)*disp((3-1)*ndim+2) &
                +        sto_wsp2(4)*disp((4-1)*ndim+2) &
                +        sto_wsp2(5)*disp((5-1)*ndim+2) &
                +        sto_wsp2(6)*disp((6-1)*ndim+2) &
                +        sto_wsp2(7)*disp((7-1)*ndim+2) &
                +        sto_wsp2(8)*disp((8-1)*ndim+2) &
                +        sto_wsp2(9)*disp((9-1)*ndim+2) &
                +        sto_wsp2(10)*disp((10-1)*ndim+2) &
                -        sto_wsp4(1)*disp((1-1)*ndim+2) &
                -        sto_wsp4(2)*disp((2-1)*ndim+2) &
                -        sto_wsp4(3)*disp((3-1)*ndim+2) &
                -        sto_wsp4(4)*disp((4-1)*ndim+2) &
                -        sto_wsp4(5)*disp((5-1)*ndim+2) &
                -        sto_wsp4(6)*disp((6-1)*ndim+2) &
                -        sto_wsp4(7)*disp((7-1)*ndim+2) &
                -        sto_wsp4(8)*disp((8-1)*ndim+2) &
                -        sto_wsp4(9)*disp((9-1)*ndim+2) &
                -        sto_wsp4(10)*disp((10-1)*ndim+2) ) &
                + factor*sto_wsp3(knodr) * (         & 
                -        sto_wsp1(1)*disp((1-1)*ndim+2) &
                -        sto_wsp1(2)*disp((2-1)*ndim+2) &
                -        sto_wsp1(3)*disp((3-1)*ndim+2) &
                -        sto_wsp1(4)*disp((4-1)*ndim+2) &
                -        sto_wsp1(5)*disp((5-1)*ndim+2) &
                -        sto_wsp1(6)*disp((6-1)*ndim+2) &
                -        sto_wsp1(7)*disp((7-1)*ndim+2) &
                -        sto_wsp1(8)*disp((8-1)*ndim+2) &
                -        sto_wsp1(9)*disp((9-1)*ndim+2) &
                -        sto_wsp1(10)*disp((10-1)*ndim+2) )	
                
           R_in(Grindz) = R_in(Grindz)  &
                + factor*shfn(knodr,igauss) * ( &
                -trmcgss*sto_wsp1(1)*disp((1-1)*ndim+3) &
                -trmcgss*sto_wsp1(2)*disp((2-1)*ndim+3) &
                -trmcgss*sto_wsp1(3)*disp((3-1)*ndim+3) &
                -trmcgss*sto_wsp1(4)*disp((4-1)*ndim+3) &
                -trmcgss*sto_wsp1(5)*disp((5-1)*ndim+3) &
                -trmcgss*sto_wsp1(6)*disp((6-1)*ndim+3) &
                -trmcgss*sto_wsp1(7)*disp((7-1)*ndim+3) &
                -trmcgss*sto_wsp1(8)*disp((8-1)*ndim+3) &
                -trmcgss*sto_wsp1(9)*disp((9-1)*ndim+3) &
                -trmcgss*sto_wsp1(10)*disp((10-1)*ndim+3) &
                +        sto_wsp2(1)*disp((1-1)*ndim+3) &
                +        sto_wsp2(2)*disp((2-1)*ndim+3) &
                +        sto_wsp2(3)*disp((3-1)*ndim+3) &
                +        sto_wsp2(4)*disp((4-1)*ndim+3) &
                +        sto_wsp2(5)*disp((5-1)*ndim+3) &
                +        sto_wsp2(6)*disp((6-1)*ndim+3) &
                +        sto_wsp2(7)*disp((7-1)*ndim+3) &
                +        sto_wsp2(8)*disp((8-1)*ndim+3) &
                +        sto_wsp2(9)*disp((9-1)*ndim+3) &
                +        sto_wsp2(10)*disp((10-1)*ndim+3) &
                -        sto_wsp4(1)*disp((1-1)*ndim+3) &
                -        sto_wsp4(2)*disp((2-1)*ndim+3) &
                -        sto_wsp4(3)*disp((3-1)*ndim+3) &
                -        sto_wsp4(4)*disp((4-1)*ndim+3) &
                -        sto_wsp4(5)*disp((5-1)*ndim+3) &
                -        sto_wsp4(6)*disp((6-1)*ndim+3) &
                -        sto_wsp4(7)*disp((7-1)*ndim+3) &
                -        sto_wsp4(8)*disp((8-1)*ndim+3) &
                -        sto_wsp4(9)*disp((9-1)*ndim+3) &
                -        sto_wsp4(10)*disp((10-1)*ndim+3) ) &
                + factor*sto_wsp3(knodr) * (          &
                -        sto_wsp1(1)*disp((1-1)*ndim+3) &
                -        sto_wsp1(2)*disp((2-1)*ndim+3) &
                -        sto_wsp1(3)*disp((3-1)*ndim+3) &
                -        sto_wsp1(4)*disp((4-1)*ndim+3) &
                -        sto_wsp1(5)*disp((5-1)*ndim+3) &
                -        sto_wsp1(6)*disp((6-1)*ndim+3) &
                -        sto_wsp1(7)*disp((7-1)*ndim+3) &
                -        sto_wsp1(8)*disp((8-1)*ndim+3) &
                -        sto_wsp1(9)*disp((9-1)*ndim+3) &
                -        sto_wsp1(10)*disp((10-1)*ndim+3) )

!                    --------------------------------------------
!                    all stiffness terms above are from physical
!                    acceleration. The isotropic linear
!                    elastic stiffness is calculated in CST
!                    --------------------------------------------

        END DO

     ENDDO	! for igauss=1,nintk_v


!     contribution to stiffness matrix from surface terms
!     ---------------------------------------------------

     DO  iface = 1,nsurf

!--      shdx is constant and already calculated !!

        CALL shcalc_3d10(shfn,shdx,dv,dvsum,shdx_av,meshpos,ndim, &
             nnode,nintk_s,sxi(1,iface),seta(1,iface), &
             szeta(1,iface),shdxi,ajacin)

        DO igauss = 1,nintk_s

           DO j = 1,ndim
              mvgss(j) = zero
              DO knode = 1,nnode
                 mvgss(j) = mvgss(j) + meshvel(j,knode)* &
                      shfn(knode,igauss)
              END DO
           END DO

!--         unit out-normal vector on each surface

           DO j = 1, 3
              wsparr(j) =meshpos(j,surf_elem(2,iface)) -  &
                   meshpos(j,surf_elem(1,iface))
              wsparr1(j)=meshpos(j,surf_elem(3,iface)) - &
                   meshpos(j,surf_elem(1,iface))
           END DO

!           do cross product: wsparr = nk * da
           
           wsp1 = wsparr(2)*wsparr1(3)  - wsparr(3)*wsparr1(2)
           wsp2 = wsparr(3)*wsparr1(1)  - wsparr(1)*wsparr1(3)
           wsp3 = wsparr(1)*wsparr1(2)  - wsparr(2)*wsparr1(1)
           wsparr(1) = wsp1 * 0.5d0
           wsparr(2) = wsp2 * 0.5d0
           wsparr(3) = wsp3 * 0.5d0


!           evaluate terms required in stiffness

           wsp1 = mvgss(1)*wsparr(1) &
                + mvgss(2)*wsparr(2) &
                + mvgss(3)*wsparr(3)

!           fill stiffness matrix


           DO knodc = 1,nnode
              sto_wsp2(knodc) =  &
                   mvgss(1)*shdx(1,knodc,igauss) + &
                   mvgss(2)*shdx(2,knodc,igauss) + &
                   mvgss(3)*shdx(3,knodc,igauss)
           END DO

           factor = third      
           DO knodr = 1, nnode
              Grindx = (lmlstet(knodr,ielem)-1)*3 + 1
              Grindy = (lmlstet(knodr,ielem)-1)*3 + 2
              Grindz = (lmlstet(knodr,ielem)-1)*3 + 3

              R_in(Grindx) = R_in(Grindx)  &
                   + mat(1)*shfn(knodr,igauss)*wsp1*(   &       
                   +   sto_wsp2(1)*disp((1-1)*ndim+1) &
                   +   sto_wsp2(2)*disp((2-1)*ndim+1) &
                   +   sto_wsp2(3)*disp((3-1)*ndim+1) &
                   +   sto_wsp2(4)*disp((4-1)*ndim+1) &
                   +   sto_wsp2(5)*disp((5-1)*ndim+1) &
                   +   sto_wsp2(6)*disp((6-1)*ndim+1) &
                   +   sto_wsp2(7)*disp((7-1)*ndim+1) &
                   +   sto_wsp2(8)*disp((8-1)*ndim+1) &
                   +   sto_wsp2(9)*disp((9-1)*ndim+1) &
                   +   sto_wsp2(10)*disp((10-1)*ndim+1) ) &
                   *    factor
              
              R_in(Grindy) = R_in(Grindy)  &
                   + mat(1)*shfn(knodr,igauss)*wsp1*(  &        
                   +   sto_wsp2(1)*disp((1-1)*ndim+2) &
                   +   sto_wsp2(2)*disp((2-1)*ndim+2) &
                   +   sto_wsp2(3)*disp((3-1)*ndim+2) &
                   +   sto_wsp2(4)*disp((4-1)*ndim+2) &
                   +   sto_wsp2(5)*disp((5-1)*ndim+2) &
                   +   sto_wsp2(6)*disp((6-1)*ndim+2) &
                   +   sto_wsp2(7)*disp((7-1)*ndim+2) &
                   +   sto_wsp2(8)*disp((8-1)*ndim+2) &
                   +   sto_wsp2(9)*disp((9-1)*ndim+2) &
                   +   sto_wsp2(10)*disp((10-1)*ndim+2) ) &
                   *   factor
              
              R_in(Grindz) = R_in(Grindz)  &
                   + mat(1)*shfn(knodr,igauss)*wsp1*(  &        
                   +   sto_wsp2(1)*disp((1-1)*ndim+3) &
                   +   sto_wsp2(2)*disp((2-1)*ndim+3) &
                   +   sto_wsp2(3)*disp((3-1)*ndim+3) &
                   +   sto_wsp2(4)*disp((4-1)*ndim+3) &
                   +   sto_wsp2(5)*disp((5-1)*ndim+3) &
                   +   sto_wsp2(6)*disp((6-1)*ndim+3) &
                   +   sto_wsp2(7)*disp((7-1)*ndim+3) &
                   +   sto_wsp2(8)*disp((8-1)*ndim+3) &
                   +   sto_wsp2(9)*disp((9-1)*ndim+3) &
                   +   sto_wsp2(10)*disp((10-1)*ndim+3) ) &
                   *   factor
           END DO
           
        ENDDO  	! igauss = 1, nintk_s 
     ENDDO	! iface=1,nsurf
            
  ENDDO	! for ielem=1,numlstet

  DO i=1,numnp*3
     Rnet(i) = Rnet(i) - R_in(i) - R_dp(i)
  ENDDO
  
111 FORMAT(i5,5e15.7)

  RETURN
END SUBROUTINE v3d10_ale

