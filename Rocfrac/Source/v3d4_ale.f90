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
SUBROUTINE v3d4_ale(v_bar,a_bar,d,vhalf,Rnet, &
     E,xnu,rho,numnp,numat_vol, &
     numcstet,matcstet,lmcstet,meshcoor, &
     nstart,nend)

  IMPLICIT NONE
  
  INTEGER :: numnp          ! number of nodes
  INTEGER :: numat_vol      ! number of materials
  INTEGER :: numcstet       ! number of tetrahedrals
! nodal mesh velocities
  REAL*8, DIMENSION(1:3*numnp) :: v_bar
! nodal mesh accelerations
  REAL*8, DIMENSION(1:3*numnp) :: a_bar
! nodal displacements
  REAL*8, DIMENSION(1:3*numnp) :: d
! nodal velocities
  REAL*8, DIMENSION(1:3*numnp) :: vhalf
! Contibution to the Net Force Vector
  REAL*8, DIMENSION(1:3*numnp) :: Rnet
! young's moduli
  REAL*8, DIMENSION(1:numat_vol) :: E
! Poisson's ratios
  REAL*8, DIMENSION(1:numat_vol) :: xnu
! densities 
  REAL*8, DIMENSION(1:numat_vol) :: rho
! mesh coordinates
  REAL*8, DIMENSION(1:3,1:numnp) :: meshcoor
! Tet connectivity table
  INTEGER, DIMENSION(1:4,1:numcstet) :: lmcstet 
  
  INTEGER matcstet(numcstet)    ! list of all materials
  INTEGER ielem , nx, ny, nz 

!     CODE FOR 3-D ALE, ELASTODYNAMIC, FIRST-ORDER SOLID ELEMENT

!     integer ndof, nnode, ndim,  nmat, nsurf
  INTEGER ndof, nnode, ndim,  nmat, nsurf

!--   for tetrahedron element !
  REAL*8 meshpos(3,4), meshvel(3,4), meshacc(3,4), mat(3)

!     INPUT:

!     ndof -- # of element degrees of freedom
!     nnode -- # of nodes on element
!     ndim -- # of spatial dimensions
!     nmat -- # of material parameters (rho, E, nu)
!     nsurf -- # of surfaces on element
!     meshpos -- coordinates of position w.r.t. to global basis of nodes
!     meshvel -- components of mesh velocity w.r.t to global basis at nodes
!     meshacc --      ''    ''  ''  accelern.  ''  ''   ''     ''  ''   ''
!     mat -- array containing material parameters (e.g. rho, E, nu)

  INTEGER  nintk_v, knode, &
       knodr, knodc,i,j,rindx,cindx,k,mm,ll, &
       iface, igauss, Grindx,Grindy,Grindz
!
!---  new variables for tetrahedron
  PARAMETER( nintk_v = 1) ! 4-noded tetrahedron

  DOUBLE PRECISION shfn(8,8), shdx(3,8,8), dv(8), dvsum, &
       mvgss(3), magss(3), mcgss(3,3), &
       wsp1,wsparr(3),trmcgss,wsp2,wsp3,wsp4, &
       third,wsp5,vxi(8),veta(8),vzeta(8),wvar, &
       sxi(4,6),seta(4,6),szeta(4,6),shdxi(3,8,8), &
       wsparr1(3),shdx_av(3,8)
  DOUBLE PRECISION disp(12), velo(12), Rimsi(12), Dimsi(12)
  DOUBLE PRECISION multiplier , ajacin(3,3), surfnormal(3,4)
  DOUBLE PRECISION wt_cst, dA, sum,factor
  DOUBLE PRECISION sto_wsp1(4),sto_wsp2(4),sto_wsp3(4),sto_wsp4(4)
  INTEGER  surf_elem(3,4), k2

  INTEGER :: nstart, nend

  ndim  = 3
  nnode = 4
  ndof  = ndim * nnode
  nsurf = 4
  multiplier = 1.d0/6.d0
      

!     initialize integration pt.coords
!     one point integration rule for tetrahedron
!     reference: The Finite Element Method by Thomas J.R. Hughes
!     Chap.3 pp 170-174 

  vxi(1)   = 0.25d0 
  veta(1)  = 0.25d0 
  vzeta(1) = 0.25d0 

!--   one point on each surface
!     face 1  : node 1 2 3 
!     face 2  : node 1 2 4
!     face 3  : node 2 3 4
!     face 4  : node 1 3 4 

  wvar = 1.d0/dsqrt(3.d0)
  third = 1.d0/3.d0
  
  sxi(1,1)   = third 
  seta(1,1)  = 0.d0 
  szeta(1,1) = third  
  
  sxi(1,2)   = third  
  seta(1,2)  = third  
  szeta(1,2) = third  
  
  sxi(1,3)   = third 
  seta(1,3)  = third 
  szeta(1,3) = 0.d0 
  
  sxi(1,4)   = 0.d0
  seta(1,4)  = third 
  szeta(1,4) = third 
  
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
  

!     loop for elements
  igauss = 1
  DO ielem = nstart, nend
     j = matcstet(ielem)
     mat(1)= rho(j)
     mat(2)= E(j)
     mat(3)= xnu(j)
     
     DO i=1, nnode
        nx = lmcstet(i,ielem)*3-2
        ny = lmcstet(i,ielem)*3-1
        nz = lmcstet(i,ielem)*3
        meshpos(1,i) = meshcoor(1,lmcstet(i,ielem))
        meshpos(2,i) = meshcoor(2,lmcstet(i,ielem))
        meshpos(3,i) = meshcoor(3,lmcstet(i,ielem))
        meshvel(1,i) = v_bar(nx)
        meshvel(2,i) = v_bar(ny)
        meshvel(3,i) = v_bar(nz)
        meshacc(1,i) = a_bar(nx)
        meshacc(2,i) = a_bar(ny)
        meshacc(3,i) = a_bar(nz)
        disp(i*3-2) = d(nx)
        disp(i*3-1) = d(ny)
        disp(i*3  ) = d(nz)
        velo(i*3-2) = vhalf(nx)
        velo(i*3-1) = vhalf(ny)
        velo(i*3  ) = vhalf(nz)
     END DO


     CALL shcalc(shfn, shdx, dv, dvsum, shdx_av, meshpos, &
             ndim,nnode,nintk_v,vxi,veta,vzeta,shdxi,ajacin)
     

!        fill in element matrices for contributions from
!        volume integration points.

  
!           evaluate mesh-related quantitites at vol. int. pt.

     mvgss(1) = 0.25d0*(meshvel(1,1)+meshvel(1,2)+ &
          meshvel(1,3)+meshvel(1,4))
     mvgss(2) = 0.25d0*(meshvel(2,1)+meshvel(2,2)+ &
          meshvel(2,3)+meshvel(2,4))
     mvgss(3) = 0.25d0*(meshvel(3,1)+meshvel(3,2)+ &
          meshvel(3,3)+meshvel(3,4))
     magss(1) = 0.25d0*(meshacc(1,1)+meshacc(1,2)+ &
          meshacc(1,3)+meshacc(1,4))
     magss(2) = 0.25d0*(meshacc(2,1)+meshacc(2,2)+ &
          meshacc(2,3)+meshacc(2,4))
     magss(3) = 0.25d0*(meshacc(3,1)+meshacc(3,2)+ &
          meshacc(3,3)+meshacc(3,4))

     mcgss(1,1) = 0.d0
     mcgss(1,2) = 0.d0
     mcgss(1,3) = 0.d0
     mcgss(2,1) = 0.d0
     mcgss(2,2) = 0.d0
     mcgss(2,3) = 0.d0
     mcgss(3,1) = 0.d0
     mcgss(3,2) = 0.d0
     mcgss(3,3) = 0.d0
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

     factor=mat(1)*dv(1)*multiplier
     DO knodr = 1,nnode
        Grindx = (lmcstet(knodr,ielem)-1)*3 + 1
        Grindy = (lmcstet(knodr,ielem)-1)*3 + 2
        Grindz = (lmcstet(knodr,ielem)-1)*3 + 3

!                    damping 
!                    ----------------------------------

        Rnet(Grindx) = Rnet(Grindx) +  &
             factor*shfn(knodr,igauss)*2.d0* &
             ( sto_wsp1(1)*velo((1-1)*ndim+1) + &
             sto_wsp1(2)*velo((2-1)*ndim+1) + &
             sto_wsp1(3)*velo((3-1)*ndim+1) + &
             sto_wsp1(4)*velo((4-1)*ndim+1) )
        
        Rnet(Grindy) = Rnet(Grindy) +  &
             factor*shfn(knodr,igauss)*2.d0* &
             ( sto_wsp1(1)*velo((1-1)*ndim+2) + &
             sto_wsp1(2)*velo((2-1)*ndim+2) + &
             sto_wsp1(3)*velo((3-1)*ndim+2) + &
             sto_wsp1(4)*velo((4-1)*ndim+2) )
        
        Rnet(Grindz) = Rnet(Grindz) +  &
             factor*shfn(knodr,igauss)*2.d0* &
             ( sto_wsp1(1)*velo((1-1)*ndim+3) + &
             sto_wsp1(2)*velo((2-1)*ndim+3) + &
             sto_wsp1(3)*velo((3-1)*ndim+3) + &
             sto_wsp1(4)*velo((4-1)*ndim+3) )
        
!                    stiffness (vol. integration terms)
!                    ----------------------------------

        Rnet(Grindx) = Rnet(Grindx)  &
             - factor*shfn(knodr,igauss) * ( &
             -trmcgss*sto_wsp1(1)*disp((1-1)*ndim+1) &
             -trmcgss*sto_wsp1(2)*disp((2-1)*ndim+1) &
             -trmcgss*sto_wsp1(3)*disp((3-1)*ndim+1) &
             -trmcgss*sto_wsp1(4)*disp((4-1)*ndim+1) &
             +        sto_wsp2(1)*disp((1-1)*ndim+1) &
             +        sto_wsp2(2)*disp((2-1)*ndim+1) &
             +        sto_wsp2(3)*disp((3-1)*ndim+1) &
             +        sto_wsp2(4)*disp((4-1)*ndim+1) &
             -        sto_wsp4(1)*disp((1-1)*ndim+1) &
             -        sto_wsp4(2)*disp((2-1)*ndim+1) &
             -        sto_wsp4(3)*disp((3-1)*ndim+1) &
             -        sto_wsp4(4)*disp((4-1)*ndim+1) ) &
             - factor*sto_wsp3(knodr) * (   &       
             -        sto_wsp1(1)*disp((1-1)*ndim+1) &
             -        sto_wsp1(2)*disp((2-1)*ndim+1) &
             -        sto_wsp1(3)*disp((3-1)*ndim+1) &
             -        sto_wsp1(4)*disp((4-1)*ndim+1) )
        

        Rnet(Grindy) = Rnet(Grindy)  &
             - factor*shfn(knodr,igauss) * ( &
             -trmcgss*sto_wsp1(1)*disp((1-1)*ndim+2) &
             -trmcgss*sto_wsp1(2)*disp((2-1)*ndim+2) &
             -trmcgss*sto_wsp1(3)*disp((3-1)*ndim+2) &
             -trmcgss*sto_wsp1(4)*disp((4-1)*ndim+2) &
             +        sto_wsp2(1)*disp((1-1)*ndim+2) &
             +        sto_wsp2(2)*disp((2-1)*ndim+2) &
             +        sto_wsp2(3)*disp((3-1)*ndim+2) &
             +        sto_wsp2(4)*disp((4-1)*ndim+2) &
             -        sto_wsp4(1)*disp((1-1)*ndim+2) &
             -        sto_wsp4(2)*disp((2-1)*ndim+2) &
             -        sto_wsp4(3)*disp((3-1)*ndim+2) &
             -        sto_wsp4(4)*disp((4-1)*ndim+2) ) &
             - factor*sto_wsp3(knodr) * (  &        
             -        sto_wsp1(1)*disp((1-1)*ndim+2) &
             -        sto_wsp1(2)*disp((2-1)*ndim+2) &
             -        sto_wsp1(3)*disp((3-1)*ndim+2) &
             -        sto_wsp1(4)*disp((4-1)*ndim+2) )
        
        Rnet(Grindz) = Rnet(Grindz)  &
             - factor*shfn(knodr,igauss) * ( &
             -trmcgss*sto_wsp1(1)*disp((1-1)*ndim+3) &
             -trmcgss*sto_wsp1(2)*disp((2-1)*ndim+3) &
             -trmcgss*sto_wsp1(3)*disp((3-1)*ndim+3) &
             -trmcgss*sto_wsp1(4)*disp((4-1)*ndim+3) &
             +        sto_wsp2(1)*disp((1-1)*ndim+3) &
             +        sto_wsp2(2)*disp((2-1)*ndim+3) &
             +        sto_wsp2(3)*disp((3-1)*ndim+3) &
             +        sto_wsp2(4)*disp((4-1)*ndim+3) &
             -        sto_wsp4(1)*disp((1-1)*ndim+3) &
             -        sto_wsp4(2)*disp((2-1)*ndim+3) &
             -        sto_wsp4(3)*disp((3-1)*ndim+3) &
             -        sto_wsp4(4)*disp((4-1)*ndim+3) ) &
             - factor*sto_wsp3(knodr) * (          &
             -        sto_wsp1(1)*disp((1-1)*ndim+3) &
             -        sto_wsp1(2)*disp((2-1)*ndim+3) &
             -        sto_wsp1(3)*disp((3-1)*ndim+3) &
             -        sto_wsp1(4)*disp((4-1)*ndim+3) )

!                    --------------------------------------------
!                    all stiffness terms above are from physical
!                    acceleration. The isotropic linear
!                    elastic stiffness is calculated in CST
!                    --------------------------------------------

     END DO


!     contribution to stiffness matrix from surface terms
!     ---------------------------------------------------

     DO iface = 1,nsurf

!--         shdx is constant and already calculated !!

        shfn(1,1)=szeta(1,iface)
        shfn(2,1)=sxi(1,iface)
        shfn(3,1)=1.d0-sxi(1,iface)-seta(1,iface)-szeta(1,iface)
        shfn(4,1)=seta(1,iface)

        DO j = 1,ndim
           mvgss(j) = 0.d0
           DO knode = 1,nnode
              mvgss(j) = mvgss(j) + meshvel(j,knode)* &
                   shfn(knode,igauss)
           END DO
        END DO
        
!--  unit out-normal vector on each surface

        DO j = 1, 3
           wsparr(j) =meshpos(j,surf_elem(2,iface)) - &
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

            
        DO knodr = 1, nnode
           Grindx = (lmcstet(knodr,ielem)-1)*3 + 1
           Grindy = (lmcstet(knodr,ielem)-1)*3 + 2
           Grindz = (lmcstet(knodr,ielem)-1)*3 + 3

           Rnet(Grindx) = Rnet(Grindx)  &
                - mat(1)*shfn(knodr,igauss)*wsp1*(  &        
                +   sto_wsp2(1)*disp((1-1)*ndim+1) &
                +   sto_wsp2(2)*disp((2-1)*ndim+1) &
                +   sto_wsp2(3)*disp((3-1)*ndim+1) &
                +   sto_wsp2(4)*disp((4-1)*ndim+1) )

           Rnet(Grindy) = Rnet(Grindy)  &
                - mat(1)*shfn(knodr,igauss)*wsp1*(  &        
                +   sto_wsp2(1)*disp((1-1)*ndim+2) &
                +   sto_wsp2(2)*disp((2-1)*ndim+2) &
                +   sto_wsp2(3)*disp((3-1)*ndim+2) &
                +   sto_wsp2(4)*disp((4-1)*ndim+2) )

           Rnet(Grindz) = Rnet(Grindz)  &
                - mat(1)*shfn(knodr,igauss)*wsp1*(  &        
                +   sto_wsp2(1)*disp((1-1)*ndim+3) &
                +   sto_wsp2(2)*disp((2-1)*ndim+3) &
                +   sto_wsp2(3)*disp((3-1)*ndim+3) &
                +   sto_wsp2(4)*disp((4-1)*ndim+3) )
        END DO
            
     END DO	! for iface=1,nsurf

1112 CONTINUE

  END DO	! for ielem=1,numcstet

111 FORMAT(i5,5e15.7)

  RETURN
END SUBROUTINE v3d4_ale

