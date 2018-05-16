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
SUBROUTINE v3d8_me(coor,MatType,ElConnVol,Rnet,disp,dmat, &
     S11,S22,S33,S12,S23,S13, &
     numnp,nstart,nend,NumEL,NumMatType, Aenh,enhanced_map,mixed_map)

!!****f* Rocfrac/Source/v3d8_me.f90
!!
!!  NAME
!!    v3d8_me
!!
!!  FUNCTION
!!
!!    Computes the internal force vector for an 8-node hexahedral
!!    mixed enhanced element formulation. 
!!
!!    Assumes small deformation, linear elastic.
!!   
!!    Original author: Ertugrul Taciroglu
!!
!!       Based on the CSAR internal document,
!!        "A Linear Mixed-Enhanced Strain Element: 
!!          Formulation, Algorithms and Verification"
!!
!!    Foundation of Theory:
!!    
!!      "A mixed-enhanced strain method Part I: Geometraically linear problems"
!!
!!
!!  INPUTS
!!
!!   NumNP   -- Number of nodes
!!   NumEL   -- Number of elements
!!   Coor    -- number of coordinates
!!   MatType -- Material id
!!   disp    -- Nodal Displacement
!!   ElConnVol  -- Nodal connectivity
!!   Dmat    -- material D matrix
!!   S11, S22, S33, S12, S23, S13 -- Stress
!!   nstart, nend -- element beginning and end loop counter
!!   NumMatVol -- number of materials
!!   Aenh -- Element enhancement parameter
!!   enhanced_map -- enhancement map
!!   mixed_map    --  mixed map
!!  
!!
!!  OUTPUT
!!
!!    Rnet -- internal force vector
!!
!!****
 
  IMPLICIT NONE
!-----Global variables
  INTEGER :: numnp          ! number of nodes
  INTEGER :: numat_vol      ! number of volumetric materials
  INTEGER :: NumEl       ! number of LSTets
  integer :: NumMatType
!--   coordinate array
  REAL*8, DIMENSION(1:3,1:numnp) :: coor
!--   connectivity table for Brick 
  INTEGER, DIMENSION(1:8,1:NumEl) :: ElConnVol
  INTEGER, DIMENSION(1:NumEl) :: MatType
!--   elastic stiffness consts
  REAL*8, DIMENSION(1:9,NumMatType) :: ci
!--   internal force
  REAL*8, DIMENSION(1:3*numnp) :: R_in
!--   displacement vector
  REAL*8, DIMENSION(1:3*numnp) :: disp
!---- Local variables
!--   node numbers
  INTEGER :: n1,n2,n3,n4,n5,n6,n7,n8,n9,n10

  REAL*8 :: Rnet(1:3*numnp)

  INTEGER :: k1n1,k1n2,k1n3,k1n4,k1n5,k1n6,k1n7,k1n8
  INTEGER :: k2n1,k2n2,k2n3,k2n4,k2n5,k2n6,k2n7,k2n8
  INTEGER :: k3n1,k3n2,k3n3,k3n4,k3n5,k3n6,k3n7,k3n8

  REAL*8, DIMENSION(1:3,1:8) :: coord
  REAL*8, DIMENSION(1:24) :: Udisp
  
  REAL*8 :: element_volume
  INTEGER :: ielem, imat
  integer :: igpt
  REAL*8, DIMENSION(1:24) :: fint
  REAL*8 nNn(8), dn(8,3), jac(3,3), jacinv(3,3), &
       dn20(20,3),nNn20(20),&
       t(3,3), tinv(3,3), &
       dmat(NumMatType,9,9), bmat(9,24), bmat2(1:8,1:9,1:24), &
       grad(9,24), tmp38(3,8), coeff(8),&
       tkront(9,9), tkrontinv(9,9),&
       e_mixed(9,12),e_mixedT(12,9), &
       gsm(12,24), gbg(12,12), bmat_avg(9,24),&
       bu(9,24), ba(9,9), bsm(12,24),&                                          ! print edge
       kuu(24,24), kau(9,24), kaa(9,9),atol, baVec(1:9,1:9,1:8), buVec(1:9,1:24,1:8),&
       tmp249(24,9),tmp91(9,1),kaa_copy(9,9),enh2(9),dnorm,bmatigs(1:8,1:9,1:24),&
       stresstmp(9,1), tmp9(9),dmatinf(NumMatType,9,9),stress(1:8,1:9),sumVect(1:9),tkrontinvT(1:9,1:9),&
       tkrontT(1:9,1:9)
  

! weight for hex element

  REAL*8, DIMENSION(1:8) :: wi = &
     (/1.000000000000000,1.000000000000000,1.000000000000000, &
       1.000000000000000,1.000000000000000,1.000000000000000, &
       1.000000000000000,1.000000000000000/)

! root3 = 1./sqrt(3.)

  REAL*8, DIMENSION(1:3,1:8) :: ri = RESHAPE( &
       (/-0.577350269189626,-0.577350269189626,-0.577350269189626, &
          0.577350269189626,-0.577350269189626,-0.577350269189626, &
          0.577350269189626, 0.577350269189626,-0.577350269189626, &
         -0.577350269189626, 0.577350269189626,-0.577350269189626, &
         -0.577350269189626,-0.577350269189626, 0.577350269189626, &
          0.577350269189626,-0.577350269189626, 0.577350269189626, &
          0.577350269189626, 0.577350269189626, 0.577350269189626, &
         -0.577350269189626, 0.577350269189626, 0.577350269189626/),(/3,8/) )

  REAL*8, DIMENSION(1:9,1:9) :: sop = RESHAPE( &
       (/1.000000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,0.500000000000000,0.000000000000000, &
         0.500000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,0.500000000000000, &
         0.000000000000000,0.000000000000000,0.000000000000000, &
         0.500000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,0.500000000000000,0.000000000000000, &
         0.500000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,1.000000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,0.500000000000000, &
         0.000000000000000,0.500000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,0.500000000000000, &
         0.000000000000000,0.000000000000000,0.000000000000000, &
         0.500000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,0.500000000000000, &
         0.000000000000000,0.500000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,1.000000000000000/),(/9,9/) )

  REAL*8, DIMENSION(1:3,1:3) :: dident = RESHAPE( &
       (/1.000000000000000,0.000000000000000,0.000000000000000, &
         0.000000000000000,1.000000000000000,0.000000000000000, &
         0.000000000000000,0.000000000000000,1.000000000000000/),(/3,3/) )

  REAL*8 :: zero = 0.0, one = 1.d0, detj,det, three = 3.d0
  LOGICAL error, debug
  REAL*8 :: alpha2
  REAL*8, DIMENSION(1:9) :: fenh
  REAL*8, DIMENSION(1:9,1:9) :: e_enhanced
  REAL*8 :: maxDisp, con1, TP, sum
  REAL*8 :: ThreeEighth = 3./8.
  REAL*8 :: Denh
  REAL*8 :: Aenh(1:9,1:NumEl)
  INTEGER :: igauss
  REAL*8, DIMENSION(8) :: xi, eta, zeta
  REAL*8, DIMENSION(8) :: xiE, etaE, zetaE
  REAL*8, DIMENSION(1:8,1:9,1:12) :: mixed_map
  REAL*8, DIMENSION(1:8,1:9,1:9)  :: enhanced_map
  REAL*8, DIMENSION(1:8,1:8) :: ShapeFun,ShapeFunE
  INTEGER :: otdev,mcrd,nnode
  INTEGER,parameter :: ngpts = 8
  REAL*8 :: strainEnh(1:ngpts,1:9,1:NumEl)
  REAL*8, dimension(1:NumEl) :: S11,S22,S33,S12,S23,S13 
  INTEGER :: i,k,j
  INTEGER :: nstart, nend



  DO ielem = nstart, nend ! Loop over elements

        imat = MatType(ielem)

        n1 = ElConnVol(1,ielem)
        n2 = ElConnVol(2,ielem)
        n3 = ElConnVol(3,ielem)
        n4 = ElConnVol(4,ielem)
        n5 = ElConnVol(5,ielem)
        n6 = ElConnVol(6,ielem)
        n7 = ElConnVol(7,ielem)
        n8 = ElConnVol(8,ielem)

        k3n1  = 3*n1
        k3n2  = 3*n2
        k3n3  = 3*n3
        k3n4  = 3*n4
        k3n5  = 3*n5
        k3n6  = 3*n6
        k3n7  = 3*n7
        k3n8  = 3*n8
     
        k2n1  = k3n1  - 1
        k2n2  = k3n2  - 1
        k2n3  = k3n3  - 1
        k2n4  = k3n4  - 1
        k2n5  = k3n5  - 1
        k2n6  = k3n6  - 1
        k2n7  = k3n7  - 1
        k2n8  = k3n8  - 1
     
        k1n1  = k3n1  - 2
        k1n2  = k3n2  - 2
        k1n3  = k3n3  - 2
        k1n4  = k3n4  - 2 
        k1n5  = k3n5  - 2
        k1n6  = k3n6  - 2
        k1n7  = k3n7  - 2
        k1n8  = k3n8  - 2

        coord(1,1) = coor(1,n1)
        coord(2,1) = coor(2,n1)
        coord(3,1) = coor(3,n1)
        
        coord(1,2) = coor(1,n2)
        coord(2,2) = coor(2,n2)
        coord(3,2) = coor(3,n2)
        
        coord(1,3) = coor(1,n3)
        coord(2,3) = coor(2,n3)
        coord(3,3) = coor(3,n3)
        
        coord(1,4) = coor(1,n4)
        coord(2,4) = coor(2,n4)
        coord(3,4) = coor(3,n4)
        
        coord(1,5) = coor(1,n5)
        coord(2,5) = coor(2,n5)
        coord(3,5) = coor(3,n5)
        
        coord(1,6) = coor(1,n6)
        coord(2,6) = coor(2,n6)
        coord(3,6) = coor(3,n6)
        
        coord(1,7) = coor(1,n7)
        coord(2,7) = coor(2,n7)
        coord(3,7) = coor(3,n7)
        
        coord(1,8) = coor(1,n8)
        coord(2,8) = coor(2,n8)
        coord(3,8) = coor(3,n8)
        
        Udisp(1) = disp(k1n1)
        Udisp(2) = disp(k2n1)
        Udisp(3) = disp(k3n1)

        Udisp(4) = disp(k1n2)
        Udisp(5) = disp(k2n2)
        Udisp(6) = disp(k3n2) 

        Udisp(7) = disp(k1n3)
        Udisp(8) = disp(k2n3)
        Udisp(9) = disp(k3n3) 

        Udisp(10) = disp(k1n4)
        Udisp(11) = disp(k2n4)
        Udisp(12) = disp(k3n4) 

        Udisp(13) = disp(k1n5)
        Udisp(14) = disp(k2n5)
        Udisp(15) = disp(k3n5)  

        Udisp(16) = disp(k1n6)
        Udisp(17) = disp(k2n6)
        Udisp(18) = disp(k3n6)   

        Udisp(19) = disp(k1n7)
        Udisp(20) = disp(k2n7)
        Udisp(21) = disp(k3n7)   

        Udisp(22) = disp(k1n8)
        Udisp(23) = disp(k2n8)
        Udisp(24) = disp(k3n8)

        fint = 0.d0

!  Compute average element quantities 

        element_volume = 0.                      ! Initialize
        t(1:3,1:3) = 0.
        bmat_avg(1:9,1:24) = 0.

!!$!     Calculate the volume of each element. First calculate the Jacobian
!!$!     at the "center" of the element.
!!$
!!$           jac(1,1) = eighth*( coord(1,1) - coord(1,2) - coord(1,3) + coord(1,4) +  &
!!$                coord(1,5) - coord(1,6) - coord(1,7) + coord(1,8) )  
!!$        
!!$           jac(1,2) = eighth*( coord(2,1) - coord(2,2) - coord(2,3) + coord(2,4) +  &
!!$                coord(2,5) - coord(2,6) - coord(2,7) + coord(2,8) )  
!!$           
!!$           jac(1,3) = eighth*( coord(3,1) - coord(3,2) - coord(3,3) + coord(3,4) +  & 
!!$                coord(3,5) - coord(3,6) - coord(3,7) + coord(3,8) )  
!!$           
!!$           jac(2,1) = eighth*( coord(1,1) + coord(1,2) - coord(1,3) - coord(1,4) +  &
!!$                coord(1,5) + coord(1,6) - coord(1,7) - coord(1,8) )  
!!$           
!!$           jac(2,2) = eighth*( coord(2,1) + coord(2,2) - coord(2,3) - coord(2,4) +  &
!!$                coord(2,5) + coord(2,6) - coord(2,7) - coord(2,8) )  
!!$           
!!$           jac(2,3) = eighth*( coord(3,1) + coord(3,2) - coord(3,3) - coord(3,4) +  &
!!$                coord(3,5) + coord(3,6) - coord(3,7) - coord(3,8) )
!!$           
!!$           jac(3,1) = eighth*( coord(1,1) + coord(1,2) + coord(1,3) + coord(1,4) -  &
!!$                coord(1,5) - coord(1,6) - coord(1,7) - coord(1,8) )  
!!$           
!!$           jac(3,2) = eighth*( coord(2,1) + coord(2,2) + coord(2,3) + coord(2,4) -  & 
!!$                coord(2,5) - coord(2,6) - coord(2,7) - coord(2,8) )  
!!$           
!!$           jac(3,3) = eighth*( coord(3,1) + coord(3,2) + coord(3,3) + coord(3,4) -  &
!!$                coord(3,5) - coord(3,6) - coord(3,7) - coord(3,8) )  
!!$
!!$        
!!$           temp_det = jac(1,1)*( jac(2,2)*jac(3,3) - jac(2,3)*jac(3,2) ) -  &
!!$                jac(1,2)*( jac(2,1)*jac(3,3) - jac(2,3)*jac(3,1) ) +  &
!!$                jac(1,3)*( jac(2,1)*jac(3,2) - jac(2,2)*jac(3,1) )
!!$        
!!$           volume = eight*temp_det            
          ! 
! 1st LOOP OVER GUAUSS POINTS

           DO igpt = 1, 8 ! LOOP over gauss points

!!$!   Calculate the derivatives of the shape functions at the gauss points
!!$            
!!$             dn(1,1) =  eighth*(one + eta(igauss))*(one + zeta(igauss))
!!$             dn(1,2) =  eighth*(one + xi(igauss) )*(one + zeta(igauss))
!!$             dn(1,3) =  eighth*(one + xi(igauss) )*(one + eta(igauss) )
!!$             
!!$             dn(2,1) = -eighth*(one + eta(igauss))*(one + zeta(igauss))
!!$             dn(2,2) =  eighth*(one - xi(igauss) )*(one + zeta(igauss))
!!$             dn(2,3) =  eighth*(one - xi(igauss) )*(one + eta(igauss) )
!!$             
!!$             dn(3,1) = -eighth*(one - eta(igauss))*(one + zeta(igauss))
!!$             dn(3,2) = -eighth*(one - xi(igauss) )*(one + zeta(igauss))
!!$             dn(3,3) =  eighth*(one - xi(igauss) )*(one - eta(igauss) )
!!$             
!!$             dn(4,1) =  eighth*(one - eta(igauss))*(one + zeta(igauss))
!!$             dn(4,2) = -eighth*(one + xi(igauss) )*(one + zeta(igauss))
!!$             dn(4,3) =  eighth*(one + xi(igauss) )*(one - eta(igauss) )
!!$             
!!$             dn(5,1) =  eighth*(one + eta(igauss))*(one - zeta(igauss))
!!$             dn(5,2) =  eighth*(one + xi(igauss) )*(one - zeta(igauss))
!!$             dn(5,3) = -eighth*(one + xi(igauss) )*(one + eta(igauss) )
!!$             
!!$             dn(6,1) = -eighth*(one + eta(igauss))*(one - zeta(igauss))
!!$             dn(6,2) =  eighth*(one - xi(igauss) )*(one - zeta(igauss))
!!$             dn(6,3) = -eighth*(one - xi(igauss) )*(one + eta(igauss) )
!!$             
!!$             dn(7,1) = -eighth*(one - eta(igauss))*(one - zeta(igauss))
!!$             dn(7,2) = -eighth*(one - xi(igauss) )*(one - zeta(igauss))
!!$             dn(7,3) = -eighth*(one - xi(igauss) )*(one - eta(igauss) )
!!$             
!!$             dn(8,1) =  eighth*(one - eta(igauss))*(one - zeta(igauss))
!!$             dn(8,2) = -eighth*(one + xi(igauss) )*(one - zeta(igauss))
!!$             dn(8,3) = -eighth*(one + xi(igauss) )*(one - eta(igauss) )
!!$
!!$!   Calculate the Jacobian, its determinant and inverse for each gauss point.
!!$            
!!$             jac(1,1) =  dn(1,1)*coord(1,1) + dn(2,1)*coord(1,2)  &
!!$                       + dn(3,1)*coord(1,3) + dn(4,1)*coord(1,4)  &
!!$                       + dn(5,1)*coord(1,5) + dn(6,1)*coord(1,6)  &
!!$                       + dn(7,1)*coord(1,7) + dn(8,1)*coord(1,8)  
!!$          
!!$             jac(2,1) =  dn(1,1)*coord(2,1) + dn(2,1)*coord(2,2)  &
!!$                       + dn(3,1)*coord(2,3) + dn(4,1)*coord(2,4)  & 
!!$                       + dn(5,1)*coord(2,5) + dn(6,1)*coord(2,6)  &
!!$                       + dn(7,1)*coord(2,7) + dn(8,1)*coord(2,8)  
!!$                                     
!!$             jac(3,1) =  dn(1,1)*coord(3,1) + dn(2,1)*coord(3,2)  &
!!$                       + dn(3,1)*coord(3,3) + dn(4,1)*coord(3,4)  &
!!$                       + dn(5,1)*coord(3,5) + dn(6,1)*coord(3,6)  & 
!!$                       + dn(7,1)*coord(3,7) + dn(8,1)*coord(3,8)  
!!$          
!!$             jac(1,2) =  dn(1,2)*coord(1,1) + dn(2,2)*coord(1,2)  &
!!$                       + dn(3,2)*coord(1,3) + dn(4,2)*coord(1,4)  & 
!!$                       + dn(5,2)*coord(1,5) + dn(6,2)*coord(1,6)  &
!!$                       + dn(7,2)*coord(1,7) + dn(8,2)*coord(1,8)  
!!$                                                         
!!$             jac(2,2) =  dn(1,2)*coord(2,1) + dn(2,2)*coord(2,2)  &
!!$                       + dn(3,2)*coord(2,3) + dn(4,2)*coord(2,4)  &
!!$                       + dn(5,2)*coord(2,5) + dn(6,2)*coord(2,6)  &
!!$                       + dn(7,2)*coord(2,7) + dn(8,2)*coord(2,8)
!!$                                                         
!!$             jac(3,2) =  dn(1,2)*coord(3,1) + dn(2,2)*coord(3,2)  &
!!$                       + dn(3,2)*coord(3,3) + dn(4,2)*coord(3,4)  & 
!!$                       + dn(5,2)*coord(3,5) + dn(6,2)*coord(3,6)  &
!!$                       + dn(7,2)*coord(3,7) + dn(8,2)*coord(3,8)
!!$          
!!$             jac(1,3) =  dn(1,3)*coord(1,1) + dn(2,3)*coord(1,2)  &
!!$                       + dn(3,3)*coord(1,3) + dn(4,3)*coord(1,4)  & 
!!$                       + dn(5,3)*coord(1,5) + dn(6,3)*coord(1,6)  &
!!$                       + dn(7,3)*coord(1,7) + dn(8,3)*coord(1,8)
!!$                                                           
!!$             jac(2,3) =  dn(1,3)*coord(2,1) + dn(2,3)*coord(2,2)  &
!!$                       + dn(3,3)*coord(2,3) + dn(4,3)*coord(2,4)  &
!!$                       + dn(5,3)*coord(2,5) + dn(6,3)*coord(2,6)  & 
!!$                       + dn(7,3)*coord(2,7) + dn(8,3)*coord(2,8)
!!$                                                           
!!$             jac(3,3) =  dn(1,3)*coord(3,1) + dn(2,3)*coord(3,2)  &
!!$                       + dn(3,3)*coord(3,3) + dn(4,3)*coord(3,4)  & 
!!$                       + dn(5,3)*coord(3,5) + dn(6,3)*coord(3,6)  &
!!$                       + dn(7,3)*coord(3,7) + dn(8,3)*coord(3,8)
!!$
!!$            det(ielem,igauss) = jac(1,1)*( jac(2,2)*jac(3,3) - jac(2,3)*jac(3,2) ) -  &
!!$                                 jac(1,2)*( jac(2,1)*jac(3,3) - jac(2,3)*jac(3,1) ) +  &
!!$                                 jac(1,3)*( jac(2,1)*jac(3,2) - jac(2,2)*jac(3,1) )
!!$          
!!$             jacinv(1,1) =  ( jac(2,2)*jac(3,3) - jac(2,3)*jac(3,2) ) / det(ielem,igauss)
!!$             jacinv(1,2) = -( jac(1,2)*jac(3,3) - jac(1,3)*jac(3,2) ) / det(ielem,igauss)
!!$             jacinv(1,3) =  ( jac(1,2)*jac(2,3) - jac(1,3)*jac(2,2) ) / det(ielem,igauss)
!!$             jacinv(2,1) = -( jac(2,1)*jac(3,3) - jac(2,3)*jac(3,1) ) / det(ielem,igauss)
!!$             jacinv(2,2) =  ( jac(1,1)*jac(3,3) - jac(1,3)*jac(3,1) ) / det(ielem,igauss)
!!$             jacinv(2,3) = -( jac(1,1)*jac(2,3) - jac(1,3)*jac(2,1) ) / det(ielem,igauss)
!!$             jacinv(3,1) =  ( jac(2,1)*jac(3,2) - jac(2,2)*jac(3,1) ) / det(ielem,igauss)
!!$             jacinv(3,2) = -( jac(1,1)*jac(3,2) - jac(1,2)*jac(3,1) ) / det(ielem,igauss)
!!$             jacinv(3,3) =  ( jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1) ) / det(ielem,igauss)  

! Get shape functions and derivatives
              CALL get_shape(ri,nNn,dn,igpt)  

!     Calculate the volume of each element. First calculate the Jacobian
!     at the "center" of the element.

! Compute jacobien
              CALL get_jacobien(coord,3,8,dn,jac,jacinv,detj,error)
            
              coeff(igpt) = detj ! * wi(igpt) = 1.d0

              element_volume = element_volume + coeff(igpt) 

              t(1:3,1:3) = t(1:3,1:3) + jac * coeff(igpt)

!  Compute conventional strain displacement matrix <bmat> at each gauss point
              
!              CALL TENSORMUL3(jacinv,3,3,'t',dn,8,3,'t',tmp38,3,8,one,zero)

              DO i = 1, 3
                 DO  j = 1, 8
                    sum = 0.d0
                    DO k = 1, 3
                       sum = sum + jacinv(k,i) * dn(j,k)
                    ENDDO
                    tmp38(i,j) = sum
                 end DO
              enddo

              CALL kronecker_product(tmp38,3,8,dident,3,3,grad,9)
              
!             CALL TENSORMUL3(sop,9,9,'n',grad,9,24,'n',bmat,9,24,one,zero)

              DO i = 1, 9
                 DO  j = 1, 24
                    sum = 0.d0
                    DO k = 1, 9
                       sum = sum + sop(i,k) * grad(k,j)
                    ENDDO
                    bmat(i,j) = sum
                 end DO
              enddo

              bmat2(igpt,:,:) = bmat(:,:)

              bmat_avg(:,:) = bmat_avg(:,:) + bmat(:,:) * coeff(igpt)
           
           END DO  ! END 1st LOOP over gauss points
           
           t(1:3,1:3) = t(1:3,1:3) / element_volume       ! Compute average jacobien <t>
           tinv(1:3,1:3) = t(1:3,1:3)                     !   and its inverse
           CALL invert3x3(tinv,det)

! Compute average <bmat>
           bmat_avg = bmat_avg / element_volume
        
! Compute product matrices
           CALL kronecker_product(t,3,3,t,3,3,tkront,9)  
           CALL kronecker_product(tinv,3,3,tinv,3,3,tkrontinv,9)

!  Compute the mixed strain-displacement matrix <bsm>

           gsm(1:12,1:24) = 0.

!! 2nd LOOP OVER GAUSS POINTS

           DO igpt = 1, 8                           ! LOOP over gauss points

              alpha2 = coeff(igpt)

              e_mixed(1:9,1:12) = mixed_map(igpt,1:9,1:12)

              bmat(1:9,1:24) = bmat2(igpt,1:9,1:24) - bmat_avg(1:9,1:24)     
!              CALL TENSORMUL3(tkront,9,9,'t',bmat,9,24,'n',grad,9,24,one,zero)

              tkrontT = transpose(tkront)
              DO i = 1, 9
                 DO  j = 1, 24
                    sum = 0.d0
                    DO k = 1, 9
                       sum = sum + tkrontT(i,k) * bmat(k,j)
                    ENDDO
                    grad(i,j) = sum
                 end DO
              enddo
           

! Accumilate G small
!              CALL TENSORMUL3(e_mixed,9,12,'t',grad,9,24,'n',gsm,12,24,alpha2,one)
              
!              e_mixedT = transpose(e_mixed)
              DO i = 1, 12
                 DO  j = 1, 24
                    sum = 0.d0
                    DO k = 1, 9
                       sum = sum + e_mixed(k,i) * grad(k,j)
                    ENDDO
                    gsm(i,j) = gsm(i,j) + alpha2*sum
                 end DO
              enddo


           END DO                  ! END LOOP over gauss points

! Compute 12x24 matrix b in (Eqn. 2.27)
!
! (Eqn. 2.29)

           bsm(1,1:24)  = gsm(1,1:24)*ThreeEighth
           bsm(2,1:24)  = gsm(2,1:24)*ThreeEighth
           bsm(3,1:24)  = gsm(3,1:24)* three*ThreeEighth
           bsm(4,1:24)  = gsm(4,1:24)*ThreeEighth
           bsm(5,1:24)  = gsm(5,1:24)*ThreeEighth
           bsm(6,1:24)  = gsm(6,1:24)* three*ThreeEighth
           bsm(7,1:24)  = gsm(7,1:24)*ThreeEighth
           bsm(8,1:24)  = gsm(8,1:24)*ThreeEighth
           bsm(9,1:24)  = gsm(9,1:24)* three*ThreeEighth
           bsm(10,1:24) = gsm(10,1:24)*ThreeEighth*0.5
           bsm(11,1:24) = gsm(11,1:24)*ThreeEighth*0.5
           bsm(12,1:24) = gsm(12,1:24)*ThreeEighth*0.5


!  Compute stiffness matrices <kaa>

           kaa(1:9, 1:9)  = 0.d0
           fenh(:) = 0.d0

!  Compute mixed/enhanced strain-displacement matrices <bu> & <ba> at 
!  each gauss point

           DO igpt = 1, 8   ! Start 3rd loop over gauss points 
           
              detj = coeff(igpt)  ! / wi(igpt) but wi(igpt) = 1. for all gauss points   
           
!             CALL get_mixed_map(e_mixed,ri,igpt)
              e_mixed(1:9,1:12) = mixed_map(igpt,1:9,1:12)

!             CALL get_enhanced_map(e_enhanced,ri,igpt)
              e_enhanced(1:9,1:9) = enhanced_map(igpt,1:9,1:9)
           
! Compute <bu> (9 x 24) (Eqn. 2.27)

!              CALL TENSORMUL3(e_mixed,9,12,'n',bsm,12,24,'n',grad,9,24,one,zero)

              DO i = 1, 9
                 DO  j = 1, 24
                    sum = 0.d0
                    DO k = 1, 12
                       sum = sum + e_mixed(i,k) * bsm(k,j)
                    ENDDO
                    grad(i,j) = sum
                 end DO
              enddo


              ! CALL TENSORMUL3(tkrontinv,9,9,'t',grad,9,24,'n',bu,9,24,one,zero)

              tkrontinvT = Transpose(tkrontinv)
              DO i = 1, 9
                 DO  j = 1, 24
                    sum = 0.d0
                    DO k = 1, 9
                       sum = sum + tkrontinvT(i,k) * grad(k,j)
                    ENDDO
                    bu(i,j) = sum
                 end DO
              enddo

              bu(1:9,1:24) = bu(1:9,1:24) / detj + bmat_avg(1:9,1:24)

! Compute <ba> (9 x 9) (Eqn. 2.31)

!              CALL TENSORMUL3(tkrontinv,9,9,'t',e_enhanced,9,9,'n',ba,9,9,one,zero)

              DO i = 1, 9
                 DO  j = 1, 9
                    sum = 0.d0
                    DO k = 1, 9
                       sum = sum + tkrontinvT(i,k) * e_enhanced(k,j)
                    ENDDO
                    ba(i,j) = sum
                 end DO
              enddo

! Compute the enhanced strain field
!               gauss pt, component, element

              DO i = 1, 9
                 sum = 0.
                 DO k = 1, 24
                    sum = sum + bu(i,k)*Udisp(k)
                 ENDDO
                 DO k = 1, 9
                    sum = sum + ba(i,k)*Aenh(k,ielem)/detj ! times ba()/detj
                 ENDDO
                 StrainEnh(igpt,i,ielem) = sum
              ENDDO
              ba(1:9,1:9) = ba(1:9,1:9) /detj ! for later
                 
!            tmp9(1:9) = Aenh(:,ielem)
!            StrainEnh(igpt,:,ielem) = MATMUL(bu,Udisp) + MATMUL(ba,tmp9)

! Compute the enhanced stress field
              
              DO i = 1, 9
                 sum = 0.
                 DO k = 1,9
                    sum = sum + dmat(imat,i,k)*StrainEnh(igpt,k,ielem)
                 ENDDO
                 stress(igpt,i) = sum
              ENDDO

!             stress(igpt,:,ielem) = MATMUL(dmat,StrainEnh(igpt,:,ielem))

! Compute the enhanced field parameter

!         T
! F    = B  S dV      (Accumulate Fenh from each guass point)
!  enh    a      

              DO i = 1, 9
                 sum = 0.d0
                 DO k = 1, 9
                    sum = sum + ba(k,i)*stress(igpt,k)
                 ENDDO
                 fenh(i) = fenh(i) + detj*sum
              ENDDO

!             CALL TENSORMUL(ba,9,9,'t',stresstmp,9,1,'n',fenh,9,detj,one)

! Compute the enhanced field parameter using the enhanced stress
         
!         T
! F    = B  S dV      (Accumulate Fint from each guass point)
!  int    u      

              DO i = 1, 24
                 sum = 0.
                 DO k = 1, 9
                    sum = sum + bu(k,i)*stress(igpt,k)
                 ENDDO
                 fint(i) = fint(i) + detj*sum
              ENDDO

!             CALL TENSORMUL(bu,9,24,'t',stresstmp,9,1,'n',fint,24,detj,one)
!
! Compute <kaa>
!

!              CALL TENSORMUL3(dmat(imat,:,:),9,9,'n',ba,9,9,'n',tkront,9,9,one,zero)

              DO i = 1, 9
                 DO  j = 1, 9
                    sum = 0.d0
                    DO k = 1, 9
                       sum = sum + dmat(imat,i,k) * ba(k,j)
                    ENDDO
                    tkront(i,j) = sum
                 end DO
              enddo


!              CALL TENSORMUL3(ba,9,9,'t',tkront,9,9,'n',kaa,9,9,detj,one)

              DO i = 1, 9
                 DO  j = 1, 9
                    sum = 0.d0
                    DO k = 1, 9
                       sum = sum + ba(k,i) * tkront(k,j)
                    ENDDO
                    kaa(i,j) = kaa(i,j) + detj*sum
                 end DO
              enddo

           END DO


!        -1
! [ K   ]
!    aa 

        CALL invert(kaa,9,det)

! Compute Eqn. (2.44) modified for explicit integration
!
!       -1                 enh
! -[ K  ]   * ( F    ) = DA    = Denh
!     aa         enh

        DO i = 1, 9

           Denh = 0.

! Matrix-Vector multiplication

           DO k = 1, 9
              Denh = Denh - kaa(i,k)*fenh(k)
           ENDDO
!
!  enh     enh   enh
! A     = A  + DA
!  i+1     i     i
!
           Aenh(i,ielem) = Aenh(i,ielem) + Denh

        ENDDO


! ASSEMBLE THE INTERNAL + Enhanced FORCE VECTOR
!
! local node 1
        Rnet(k1n1)  = Rnet(k1n1)  - fint(1)
        Rnet(k2n1)  = Rnet(k2n1)  - fint(2)
        Rnet(k3n1)  = Rnet(k3n1)  - fint(3)
! local node 2
        Rnet(k1n2)  = Rnet(k1n2)  - fint(4)
        Rnet(k2n2)  = Rnet(k2n2)  - fint(5)
        Rnet(k3n2)  = Rnet(k3n2)  - fint(6)
! local node 3
        Rnet(k1n3)  = Rnet(k1n3)  - fint(7)
        Rnet(k2n3)  = Rnet(k2n3)  - fint(8)
        Rnet(k3n3)  = Rnet(k3n3)  - fint(9)
! local node 4
        Rnet(k1n4)  = Rnet(k1n4)  - fint(10)
        Rnet(k2n4)  = Rnet(k2n4)  - fint(11)
        Rnet(k3n4)  = Rnet(k3n4)  - fint(12)
! local node 5
        Rnet(k1n5)  = Rnet(k1n5)  - fint(13)
        Rnet(k2n5)  = Rnet(k2n5)  - fint(14)
        Rnet(k3n5)  = Rnet(k3n5)  - fint(15)
! local node 6
        Rnet(k1n6)  = Rnet(k1n6)  - fint(16)
        Rnet(k2n6)  = Rnet(k2n6)  - fint(17)
        Rnet(k3n6)  = Rnet(k3n6)  - fint(18)
! local node 7
        Rnet(k1n7)  = Rnet(k1n7)  - fint(19)
        Rnet(k2n7)  = Rnet(k2n7)  - fint(20)
        Rnet(k3n7)  = Rnet(k3n7)  - fint(21)
! local node 8
        Rnet(k1n8)  = Rnet(k1n8)  - fint(22)
        Rnet(k2n8)  = Rnet(k2n8)  - fint(23)
        Rnet(k3n8)  = Rnet(k3n8)  - fint(24)

! average stress over an element

        sumVect = 0.d0
        DO k = 1, 9
           DO i = 1,8
              sumVect(k) = sumVect(k) + stress(i,k)*.125d0
           enddo
        enddo

        S11(ielem) = sumVect(1)
        S22(ielem) = sumVect(5)
        S33(ielem) = sumVect(9)
        S12(ielem) = sumVect(4)
        S23(ielem) = sumVect(8)
        S13(ielem) = sumVect(7)

     ENDDO
     
!     PRINT*,'Max Aenh =',MAXVAL(Aenh(:,:))

   END SUBROUTINE v3d8_me



!--------------------------------------------------------------GET_SHAPE

SUBROUTINE get_shape(r,n,dn,igpt)
! *--------------------------------------------------------------------*
! |                                                                    |
! |    Returns the shape function N and its derivative matrix DN       |
! |                                                                    |
! |    <r>    : 3 x  ngpts vector of coordinates of the isoparametric  |
! |             point (usually an integration point)                   |
! |    <igpt> : no. of the particular gauss point                      |
! |    <n>    : 8 x 1 shape function matrix  evaluated at <ri>         |
! |    <dn>   : 8 x 3 shape function derivative matrix                 |
! |             evaluated at <ri>                                      |
! |                                                                    |
! *--------------------------------------------------------------------*
      
  IMPLICIT real*8 (a-h, o-z)
!      include 'ABA_PARAM.INC'
      
!  Argument variables

  DOUBLE PRECISION r(3,*), n(*), dn(8,*)

!  Data statements

  DATA one /1.000000000000000/

!  Initializations

  r1 = r(1,igpt)
  r2 = r(2,igpt)
  r3 = r(3,igpt)

!  Compute shape functions <n>

  n(1) = -(r1 - one) * (r2 - one) * (r3 - one) 
  n(2) =  (r1 + one) * (r2 - one) * (r3 - one) 
  n(3) = -(r1 + one) * (r2 + one) * (r3 - one) 
  n(4) =  (r1 - one) * (r2 + one) * (r3 - one) 
  n(5) =  (r1 - one) * (r2 - one) * (r3 + one) 
  n(6) = -(r1 + one) * (r2 - one) * (r3 + one) 
  n(7) =  (r1 + one) * (r2 + one) * (r3 + one)
  n(8) = -(r1 - one) * (r2 + one) * (r3 + one) 

  n(1:8) = 0.125000000000000 * n(1:8)
  
!  Compute shape function derivative <dn>

  dn(1, 1) =  -(r2 - one) * (r3 - one) 
  dn(2, 1) =   (r2 - one) * (r3 - one) 
  dn(3, 1) =  -(r2 + one) * (r3 - one)
  dn(4, 1) =   (r2 + one) * (r3 - one)
  dn(5, 1) =   (r2 - one) * (r3 + one)
  dn(6, 1) =  -(r2 - one) * (r3 + one)
  dn(7, 1) =   (r2 + one) * (r3 + one)
  dn(8, 1) =  -(r2 + one) * (r3 + one)
  
  dn(1, 2) =  -(r1 - one) * (r3 - one)
  dn(2, 2) =   (r1 + one) * (r3 - one)
  dn(3, 2) =  -(r1 + one) * (r3 - one)
  dn(4, 2) =   (r1 - one) * (r3 - one)
  dn(5, 2) =   (r1 - one) * (r3 + one)
  dn(6, 2) =  -(r1 + one) * (r3 + one)
  dn(7, 2) =   (r1 + one) * (r3 + one)
  dn(8, 2) =  -(r1 - one) * (r3 + one)
  
  dn(1, 3) =  -(r1 - one) * (r2 - one) 
  dn(2, 3) =   (r1 + one) * (r2 - one) 
  dn(3, 3) =  -(r1 + one) * (r2 + one)
  dn(4, 3) =   (r1 - one) * (r2 + one) 
  dn(5, 3) =   (r1 - one) * (r2 - one)
  dn(6, 3) =  -(r1 + one) * (r2 - one) 
  dn(7, 3) =   (r1 + one) * (r2 + one) 
  dn(8, 3) =  -(r1 - one) * (r2 + one)
  
  dn(1:8, 1:3) = 0.125000000000000 * dn(1:8, 1:3)

  RETURN
END SUBROUTINE get_shape

SUBROUTINE get_shape20(r,n,dn,igpt)
! *--------------------------------------------------------------------*
! |                                                                    |
! |    Returns the shape function N and its derivative matrix DN       |
! |                                                                    |
! |    <r>    : 3 x  ngpts vector of coordinates of the isoparametric  |
! |             point (usually an integration point)                   |
! |    <igpt> : no. of the particular gauss point                      |
! |    <n>    : 8 x 1 shape function matrix  evaluated at <ri>         |
! |    <dn>   : 8 x 3 shape function derivative matrix                 |
! |             evaluated at <ri>                                      |
! |                                                                    |
! *--------------------------------------------------------------------*
      
  IMPLICIT real*8 (a-h, o-z)
!      include 'ABA_PARAM.INC'
      
!  Argument variables

  DOUBLE PRECISION r(3,*), n(*), dn(20,*)

  INTEGER :: xii(20), etai(20), zetai(20)

!  Data statements

  DATA one /1.000000000000000/

!  Initializations

  xi = r(1,igpt)
  eta = r(2,igpt)
  zeta = r(3,igpt)

  xii=(/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
  etai=(/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
  zetai=(/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)

  DO l = 1, 20
     xi0=xi*xii(l)
     eta0=eta*etai(l)
     zeta0=zeta*zetai(l)

!  Compute shape functions <n>
     IF(l==4.OR.l==8.OR.l==16.OR.l==20) THEN
        n(l)=.25*(1.-xi*xi)*(1.+eta0)*(1.+zeta0)
     ELSE IF(l>=9.AND.l<=12)THEN
        n(l)=.25*(1.+xi0)*(1.-eta*eta)*(1.+zeta0)
     ELSE IF(l==2.OR.l==6.OR.l==14.OR.l==18) THEN
        n(l)=.25*(1.+xi0)*(1.+eta0)*(1.-zeta*zeta)
     ELSE
        n(l)=.125*(1.+xi0)*(1.+eta0)*(1.+zeta0)*(xi0+eta0+zeta0-2)
     END IF

!  Compute shape function derivative <dn>
     IF(l==4.OR.l==8.OR.l==16.OR.l==20) THEN
        dn(l,1)=-.5*xi*(1.+eta0)*(1.+zeta0)
        dn(l,2)=.25*etai(l)*(1.-xi*xi)*(1.+zeta0)
        dn(l,3)=.25*zetai(l)*(1.-xi*xi)*(1.+eta0)
     ELSE IF(l>=9.AND.l<=12)THEN
        dn(l,1)=.25*xii(l)*(1.-eta*eta)*(1.+zeta0)
        dn(l,2)=-.5*eta*(1.+xi0)*(1.+zeta0)
        dn(l,3)=.25*zetai(l)*(1.+xi0)*(1.-eta*eta)
     ELSE IF(l==2.OR.l==6.OR.l==14.OR.l==18) THEN
        dn(l,1)=.25*xii(l)*(1.+eta0)*(1.-zeta*zeta)
        dn(l,2)=.25*etai(l)*(1.+xi0)*(1.-zeta*zeta)
        dn(l,3)=-.5*zeta*(1.+xi0)*(1.+eta0)
     ELSE
        dn(l,1)=.125*xii(l)*(1.+eta0)*(1.+zeta0)*(2.*xi0+eta0+zeta0-1.)
        dn(l,2)=.125*etai(l)*(1.+xi0)*(1.+zeta0)*(xi0+2.*eta0+zeta0-1.)
        dn(l,3)=.125*zetai(l)*(1.+xi0)*(1.+eta0)*(xi0+eta0+2.*zeta0-1.)
     END IF
  ENDDO

  RETURN
END SUBROUTINE get_shape20



!-----------------------------------------------------------GET_JACOBIEN
SUBROUTINE get_jacobien(coords,mcrd,nnode,dn,&
     jac,jacinv,detj,error)
! *--------------------------------------------------------------------*
! |                                                                    |
! |    Returns the jacobien, its inverse and its determinant           |
! |                                                                    |
! |    <coords> : mcrd x nnode matrix of element nodal coordinates     |
! |    <dn>     : 8 x 3  shape function derivative matrix              |
! |    <jac>    :   x    jacobien matrix                               |
! |    <dn>     : 8 x 3  inverse jacobien matrix                       |
! |    <detj>   : determinant of the jacobien matrix                   |
! |    <error>  : logical varible; set to true if determinant is zero  |
! |                                                                    |
! |     calls   :  TENSORMUL, INVERT3x3                                   |
! |                                                                    |
! *--------------------------------------------------------------------*

      IMPLICIT real*8 (a-h, o-z)
!      include 'ABA_PARAM.INC'
      
!  Argument variables

  DOUBLE PRECISION jac(mcrd,mcrd), jacinv(mcrd,mcrd)
  DIMENSION coords(1:3,1:nnode), dn(nnode,mcrd)
  LOGICAL error

!  Data statements

  DATA zero,one /0.000000000000000,1.000000000000000/

!  Initialize variables

  det = zero
  error = .TRUE.

!  Compute jacobien

  CALL TENSORMUL(coords, mcrd, nnode,'n', &
       dn, nnode, mcrd, 'n', jac, mcrd, one, zero)

!  call printmat(jac,mcrd,mcrd,7,' jac ')

!  Compute inverse jacobien and determinant

  jacinv = jac
  
  CALL invert3x3(jacinv,detj)

!      call printmat(jacinv,mcrd,mcrd,7,' jinv')

!  PRINT*, 'detj =', detj
  IF(detj.NE.zero) error = .FALSE.

  RETURN
END SUBROUTINE get_jacobien

!-----------------------------------------------------------------TENSORMUL
SUBROUTINE TENSORMUL(a,lda,n,achar,b,ldb,m,bchar,c,ldc,alpha,beta)
! *--------------------------------------------------------------------*
! |    Matrix multiplication routine                                   |
! |                                                                    |
! |    computes one of the following :                                 |
! |                C = alpha * A B + beta C                            |
! |                C = alpha * Trans(A) B  + beta C                    |
! |                C = alpha * A Trans(B) + beta C                     |
! |                C = alpha * Trans(A) Trans(B) + beta C              |
! |    depending the character input                                   |
! |                aform : 'n' 'N' 't' 'T'                             |
! |                bform : 'n' 'N' 't' 'T'                             |
! |    where 'n' stands for normal 't' stands for transpose            |
! |                                                                    |
! |    Note that:                                                      |
! |                                                                    |
! |        INPUTs----lda : leading dimension of A                      |
! |                    n : the other dimension of A                    |
! |                  ldb : leading dimension of B                      |
! |                    m : the other dimension of B                    |
! |                  ldc : leading dim. of C (must be large enough)    |
! |                alpha : a scalar multiplier to initialize C         |
! |                                                                    |
! |        OUTPUTs-----C : is the product matrix                       |
! |                                                                    |
! *--------------------------------------------------------------------*


      IMPLICIT REAL*8  (a-h,o-z)

!  Argument variables

      DIMENSION a(lda,*), b(ldb,*), c(ldc,*)
      CHARACTER*1 achar,bchar

!  Local variables

      PARAMETER (max = 30)
      DIMENSION anew(max,max), bnew(max,max)
      INTEGER crow, ccol

!  Check character switches and setup indices

      IF(achar.EQ.'n'.OR.achar.EQ.'N') THEN
         ia = lda
         ja = n
         crow = lda
         DO i=1,ia
            DO j=1,ja
               anew(i,j) = a(i,j)
            END DO
         END DO
      ELSEIF(achar.EQ.'t'.OR.achar.EQ.'T') THEN
         ia = n
         ja = lda
         crow = n
         DO i=1,ia
            DO j=1,ja
               anew(i,j) = a(j,i)
            END DO
         END DO
      ELSE
         WRITE(7,*) ' >> unrecognized character, ',achar
         STOP
      END IF


      IF(bchar.EQ.'n'.OR.bchar.EQ.'N') THEN
         ib = ldb
         jb = m
         ccol = m
         DO i=1,ib
            DO j=1,jb
               bnew(i,j) = b(i,j)
            END DO
         END DO
      ELSEIF(bchar.EQ.'t'.OR.bchar.EQ.'T') THEN
         ib = m
         jb = ldb
         ccol = ldb
         DO i=1,ib
            DO j=1,jb
               bnew(i,j) = b(j,i)
            END DO
         END DO
      ELSE
         WRITE(7,*) ' >> unrecognized character, ',bchar
         STOP
      END IF

      IF (ja.NE.ib) THEN
         WRITE(7,*) ' >> incompatible matrices'
         STOP
      END IF

!  Compute C

      c(1:crow,1:ccol) = beta * c(1:crow,1:ccol)
      DO i=1, crow
         DO j=1,ccol
            DO k=1, ja
               c(i,j) = c(i,j) + alpha * anew(i,k)*bnew(k,j)
            END DO
         END DO
      END DO

      RETURN
    END SUBROUTINE TENSORMUL

!!$!-----------------------------------------------------------------TENSORMUL
!!$    SUBROUTINE TENSORMUL3(a,lda,n,achar,b,ldb,m,bchar,c,ldc,nn,alpha,beta)
!!$! *--------------------------------------------------------------------*
!!$! |    Matrix multiplication routine                                   |
!!$! |                                                                    |
!!$! |    computes one of the following :                                 |
!!$! |                C = alpha * A B + beta C                            |
!!$! |                C = alpha * Trans(A) B  + beta C                    |
!!$! |                C = alpha * A Trans(B) + beta C                     |
!!$! |                C = alpha * Trans(A) Trans(B) + beta C              |
!!$! |    depending the character input                                   |
!!$! |                aform : 'n' 'N' 't' 'T'                             |
!!$! |                bform : 'n' 'N' 't' 'T'                             |
!!$! |    where 'n' stands for normal 't' stands for transpose            |
!!$! |                                                                    |
!!$! |    Note that:                                                      |
!!$! |                                                                    |
!!$! |        INPUTs----lda : leading dimension of A                      |
!!$! |                    n : the other dimension of A                    |
!!$! |                  ldb : leading dimension of B                      |
!!$! |                    m : the other dimension of B                    |
!!$! |                  ldc : leading dim. of C (must be large enough)    |
!!$! |                alpha : a scalar multiplier to initialize C         |
!!$! |                                                                    |
!!$! |        OUTPUTs-----C : is the product matrix                       |
!!$! |                                                                    |
!!$! *--------------------------------------------------------------------*
!!$
!!$
!!$      IMPLICIT REAL*8(a-h,o-z)
!!$
!!$!  Argument variables
!!$
!!$      DIMENSION a(lda,n), b(ldb,m), c(ldc,nn)
!!$      CHARACTER*1 achar,bchar
!!$
!!$!  Local variables
!!$
!!$      PARAMETER (max = 30)
!!$      DIMENSION anew(max,max), bnew(max,max)
!!$      INTEGER crow, ccol
!!$
!!$!  Check character switches and setup indices
!!$
!!$      IF(achar.EQ.'n'.OR.achar.EQ.'N') THEN
!!$         ia = lda
!!$         ja = n
!!$         crow = lda
!!$         DO i=1,ia
!!$            DO j=1,ja
!!$               anew(i,j) = a(i,j)
!!$            END DO
!!$         END DO
!!$      ELSEIF(achar.EQ.'t'.OR.achar.EQ.'T') THEN
!!$         ia = n
!!$         ja = lda
!!$         crow = n
!!$         DO i=1,ia
!!$            DO j=1,ja
!!$               anew(i,j) = a(j,i)
!!$            END DO
!!$         END DO
!!$      ELSE
!!$         WRITE(7,*) ' >> unrecognized character, ',achar
!!$         STOP
!!$      END IF
!!$
!!$
!!$      IF(bchar.EQ.'n'.OR.bchar.EQ.'N') THEN
!!$         ib = ldb
!!$         jb = m
!!$         ccol = m
!!$         DO i=1,ib
!!$            DO j=1,jb
!!$               bnew(i,j) = b(i,j)
!!$            END DO
!!$         END DO
!!$      ELSEIF(bchar.EQ.'t'.OR.bchar.EQ.'T') THEN
!!$         ib = m
!!$         jb = ldb
!!$         ccol = ldb
!!$         DO i=1,ib
!!$            DO j=1,jb
!!$               bnew(i,j) = b(j,i)
!!$            END DO
!!$         END DO
!!$      ELSE
!!$         WRITE(7,*) ' >> unrecognized character, ',bchar
!!$         STOP
!!$      END IF
!!$
!!$      IF (ja.NE.ib) THEN
!!$         WRITE(7,*) ' >> incompatible matrices'
!!$         STOP
!!$      END IF
!!$
!!$!  Compute C
!!$
!!$      c(1:crow,1:ccol) = beta * c(1:crow,1:ccol)
!!$      DO i=1, crow
!!$         DO j=1,ccol
!!$            DO k=1, ja
!!$               c(i,j) = c(i,j) + alpha * anew(i,k)*bnew(k,j)
!!$            END DO
!!$         END DO
!!$      END DO
!!$
!!$      RETURN
!!$    END SUBROUTINE TENSORMUL3
!!$
!!$!-----------------------------------------------------------------TENSORMUL3
!!$SUBROUTINE TENSORMULSQR(a,lda,n,achar,b,ldb,m,bchar,c,ldc,alpha,beta)
!!$! *--------------------------------------------------------------------*
!!$! |    Matrix multiplication routine                                   |
!!$! |                                                                    |
!!$! |    computes one of the following :                                 |
!!$! |                C = alpha * A B + beta C                            |
!!$! |                C = alpha * Trans(A) B  + beta C                    |
!!$! |                C = alpha * A Trans(B) + beta C                     |
!!$! |                C = alpha * Trans(A) Trans(B) + beta C              |
!!$! |    depending the character input                                   |
!!$! |                aform : 'n' 'N' 't' 'T'                             |
!!$! |                bform : 'n' 'N' 't' 'T'                             |
!!$! |    where 'n' stands for normal 't' stands for transpose            |
!!$! |                                                                    |
!!$! |    Note that:                                                      |
!!$! |                                                                    |
!!$! |        INPUTs----lda : leading dimension of A                      |
!!$! |                    n : the other dimension of A                    |
!!$! |                  ldb : leading dimension of B                      |
!!$! |                    m : the other dimension of B                    |
!!$! |                  ldc : leading dim. of C (must be large enough)    |
!!$! |                alpha : a scalar multiplier to initialize C         |
!!$! |                                                                    |
!!$! |        OUTPUTs-----C : is the product matrix                       |
!!$! |                                                                    |
!!$! *--------------------------------------------------------------------*
!!$
!!$
!!$      IMPLICIT REAL*8  (a-h,o-z)
!!$
!!$!  Argument variables
!!$
!!$      DIMENSION a(9,9), b(9,9), c(9,9)
!!$      CHARACTER*1 achar,bchar
!!$
!!$!  Local variables
!!$
!!$      PARAMETER (max = 30)
!!$      DIMENSION anew(max,max), bnew(max,max)
!!$      INTEGER crow, ccol
!!$
!!$!  Check character switches and setup indices
!!$
!!$      IF(achar.EQ.'n'.OR.achar.EQ.'N') THEN
!!$         ia = lda
!!$         ja = n
!!$         crow = lda
!!$         DO i=1,ia
!!$            DO j=1,ja
!!$               anew(i,j) = a(i,j)
!!$            END DO
!!$         END DO
!!$      ELSEIF(achar.EQ.'t'.OR.achar.EQ.'T') THEN
!!$         ia = n
!!$         ja = lda
!!$         crow = n
!!$         DO i=1,ia
!!$            DO j=1,ja
!!$               anew(i,j) = a(j,i)
!!$            END DO
!!$         END DO
!!$      ELSE
!!$         WRITE(7,*) ' >> unrecognized character, ',achar
!!$         STOP
!!$      END IF
!!$
!!$
!!$      IF(bchar.EQ.'n'.OR.bchar.EQ.'N') THEN
!!$         ib = ldb
!!$         jb = m
!!$         ccol = m
!!$         DO i=1,ib
!!$            DO j=1,jb
!!$               bnew(i,j) = b(i,j)
!!$            END DO
!!$         END DO
!!$      ELSEIF(bchar.EQ.'t'.OR.bchar.EQ.'T') THEN
!!$         ib = m
!!$         jb = ldb
!!$         ccol = ldb
!!$         DO i=1,ib
!!$            DO j=1,jb
!!$               bnew(i,j) = b(j,i)
!!$            END DO
!!$         END DO
!!$      ELSE
!!$         WRITE(7,*) ' >> unrecognized character, ',bchar
!!$         STOP
!!$      END IF
!!$
!!$      IF (ja.NE.ib) THEN
!!$         WRITE(7,*) ' >> incompatible matrices'
!!$         STOP
!!$      END IF
!!$
!!$!  Compute C
!!$
!!$      c(1:crow,1:ccol) = beta * c(1:crow,1:ccol)
!!$      DO i=1, crow
!!$         DO j=1,ccol
!!$            DO k=1, ja
!!$               c(i,j) = c(i,j) + alpha * anew(i,k)*bnew(k,j)
!!$            END DO
!!$         END DO
!!$      END DO
!!$
!!$      RETURN
!!$    END SUBROUTINE TENSORMULSQR
!!$

!-----------------------------------------------------------------INVERT

SUBROUTINE invert(a,nrow,det)

! *--------------------------------------------------------------------*
! |                                                                    |
! |   Inverts a matrix by gauss jordan elimination and computes the    |
! |   the determinant. The matrix may be symmetric or nonsymmetric.    |
! |                                                                    |
! |    <a>    : a square, double precision matrix                      |
! |    <nrow> : dimensioned number of rows for 'a'. also the actual    |
! |             sizes for inversion computations                       |
! |    <det>  : determinant of the matrix ( zero if matrix is singular)|
! |                                                                    |
! |     NOTE  : the matrix is overwritten by the inverse.              |
! |                                                                    |
! *--------------------------------------------------------------------*
  IMPLICIT NONE

!  Argument variables

  INTEGER nrow
  REAL*8, DIMENSION(1:nrow,1:nrow) :: a
  REAL*8 :: det

!  Local variables

  INTEGER ncol, k, j, i

!  Invert matrix

  det  = 1.
  ncol = nrow
  
  DO  k = 1, nrow
     
     IF( a(k,k) .EQ. 0. ) THEN ! Check for singular matrix
        det = 0.
        RETURN
     END IF
     
     det = det * a(k, k)
     
     DO j = 1, ncol
        IF( j .NE. k ) a(k,j) = a(k,j) / a(k,k)
     END DO
     
     a(k,k) = 1. / a(k,k)
     
     DO  i = 1, nrow
        IF( i .EQ. k) CYCLE
        DO j = 1, ncol
           IF( j .NE. k) a(i,j) = a(i,j) - a(i,k) * a(k,j)
        END DO
        a(i,k) = -a(i,k) * a(k,k)
     END DO
     
  END DO

  RETURN
END SUBROUTINE invert

!--------------------------------------------------------------INVERT3x3

SUBROUTINE invert3x3(a,det)
! *--------------------------------------------------------------------*
! |                                                                    |
! |   Inverts 3 by 3 matrix and computes the determinant               |
! |   The matrix may be symmetric or nonsymmetric.                     |
! |                                                                    |
! |    <a>    : 3x3 input matrix                                       |
! |    <det>  : determinant of the matrix ( zero if matrix is singular)|
! |                                                                    |
! |     NOTE  : the matrix is overwritten by the inverse.              |
! |                                                                    |
! *--------------------------------------------------------------------*
  IMPLICIT NONE

!  Argument variables

      DOUBLE PRECISION a(3,3), det

!  Local variables

      DOUBLE PRECISION acopy(3,3), zero

!  Data statements

      DATA zero /0.000000000000000/

!  Invert matrix

      acopy = a

      det = -acopy(1,3)*acopy(2,2)*acopy(3,1) &
           +acopy(1,2)*acopy(2,3)*acopy(3,1) &
           +acopy(1,3)*acopy(2,1)*acopy(3,2) &
           -acopy(1,1)*acopy(2,3)*acopy(3,2) &
           -acopy(1,2)*acopy(2,1)*acopy(3,3) &
           +acopy(1,1)*acopy(2,2)*acopy(3,3)

      IF(det.EQ.zero) THEN      ! Check for singular matrix
         a(1:3,1:3) = zero
         RETURN
      END IF

      a(1,1) = -acopy(2,3)*acopy(3,2) + acopy(2,2)*acopy(3,3)
      a(2,1) =  acopy(2,3)*acopy(3,1) - acopy(2,1)*acopy(3,3)
      a(3,1) = -acopy(2,2)*acopy(3,1) + acopy(2,1)*acopy(3,2)  
      a(1,2) =  acopy(1,3)*acopy(3,2) - acopy(1,2)*acopy(3,3)
      a(2,2) = -acopy(1,3)*acopy(3,1) + acopy(1,1)*acopy(3,3)
      a(3,2) =  acopy(1,2)*acopy(3,1) - acopy(1,1)*acopy(3,2)
      a(1,3) = -acopy(1,3)*acopy(2,2) + acopy(1,2)*acopy(2,3)
      a(2,3) =  acopy(1,3)*acopy(2,1) - acopy(1,1)*acopy(2,3)
      a(3,3) = -acopy(1,2)*acopy(2,1) + acopy(1,1)*acopy(2,2) 

      a = a / det

      RETURN
    END SUBROUTINE invert3x3
!----------------------------------------------------------GET_MIXED_MAP
      SUBROUTINE get_mixed_map(e_mixed,r,igpt)
! *--------------------------------------------------------------------*
! |                                                                    |
! |    Returns the linear mapping matrix <e_mixed>                     |
! |    evaluated at a given isoparametric coordinate <r>               |
! |                                                                    |
! |    <r>          : 3 x ? vector of coordinates of the isoparametric |
! |                   point (usually an integration point)             |
! |    <e_mixed>    : 9 x 12 matrix that interpolates the mixed field  |
! |    <igpt>       : no. of the particular gauss point                |
! |                                                                    |
! |     NOTE      : 9  is no. of stress/strain components              |
! |                 12 is no. of mixed field parameters                |
! |                                                                    |
! *--------------------------------------------------------------------*

      IMPLICIT Real*8 (a-h, o-z)
!      include 'ABA_PARAM.INC'
      
!  Argument variables

      DIMENSION e_mixed(9,12), r(3,*)

!  Data statements

      DATA zero /0.000000000000000/

!  Initialize variables

      r1 = r(1,igpt)
      r2 = r(2,igpt)
      r3 = r(3,igpt)

      e_mixed = zero

!  Evaluate <e_mixed>

      e_mixed(1,1)  = r2
      e_mixed(1,2)  = r3
      e_mixed(1,3)  = r2 * r3

      e_mixed(2,10) = r3
      e_mixed(3,12) = r2
      e_mixed(4,10) = r3

      e_mixed(5,4)  = r1        
      e_mixed(5,5)  = r3
      e_mixed(5,6)  = r1 * r3

      e_mixed(6,11) = r1
      e_mixed(7,12) = r2
      e_mixed(8,11) = r1

      e_mixed(9,7)  = r1        
      e_mixed(9,8)  = r2
      e_mixed(9,9)  = r1 * r2

      RETURN
    END SUBROUTINE get_mixed_map

!-------------------------------------------------------GET_ENHANCED_MAP
    SUBROUTINE get_enhanced_map(e_enhanced,r,igpt)
! *--------------------------------------------------------------------*
! |                                                                    |
! |    Returns the linear mapping matrix <e_enhanced>                  |
! |    evaluated at a given isoparametric coordinate <r>               |
! |                                                                    |
! |    <r>          : 3 x ? vector of coordinates of the isoparametric |
! |                   point (usually an integration point)             |
! |    <e_enhanced> : 9 x 9  matrix that interpolates the enhanced     |
! |                   field                                            |
! |    <igpt>       : no. of the particular gauss point                |
! |                                                                    |
! |     NOTE      : 9  is no. of stress/strain components              |
! |                 9  is no. of enhanced field parameters             |
! |                                                                    |
! *--------------------------------------------------------------------*

      IMPLICIT real*8 (a-h, o-z)
!      include 'ABA_PARAM.INC'
      
!  Argument variables

      DIMENSION e_enhanced(9,9), r(3,*)

!  Data statements

      DATA zero /0.000000000000000/

!  Initialize variables

      r1 = r(1,igpt)
      r2 = r(2,igpt)
      r3 = r(3,igpt)

      e_enhanced = zero

!  Evaluate <e_enhanced>

      e_enhanced(1,1)  = r1    
      e_enhanced(1,2)  = r1 * r2
      e_enhanced(1,3)  = r1 * r3

      e_enhanced(5,4)  = r2
      e_enhanced(5,5)  = r2 * r3
      e_enhanced(5,6)  = r2 * r1

      e_enhanced(9,7)  = r3        
      e_enhanced(9,8)  = r3 * r2
      e_enhanced(9,9)  = r3 * r1

      RETURN
    END SUBROUTINE get_enhanced_map
!------------------------------------------------------KRONECKER_PRODUCT
      SUBROUTINE kronecker_product(a,lda,nda,b,ldb,ndb,c,ldc)
! *--------------------------------------------------------------------*
! |  This subroutine takes the Kronecker Product of the two matrices   |
! |                                                                    |
! |              C = A  <kron_prod> B                                  |
! |                                                                    |
! |      INPUTs----lda : leading dimension of A                        |
! |                nda : the other dimension of A                      |
! |                ldb : leading dimension of B                        |
! |                ndb : the other dimension of B                      |
! |                ldc : leading dim. of C ( must be >= lda*ldb )      |
! |                                                                    |
! |      OUTPUTs-----C : is the product matrix                         |
! |                                                                    |
! *--------------------------------------------------------------------*

!     include 'ABA_PARAM.INC'
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

! Argument variables

      DIMENSION a(lda,nda), b(ldb,ndb), c(ldc,*)

!  Local variables

      INTEGER ia, ja, ist, iend, jst, jend

!  Compute Kronecker Product
      
      DO ja=1,nda
         DO ia=1,lda
         
            ist  = (ia - 1) * ldb + 1
            iend = ist + ldb - 1
            jst  = (ja - 1) * ndb + 1
            jend = jst + ndb - 1

            c(ist:iend,jst:jend) = a(ia,ja) * b(1:ldb,1:ndb)

         END DO
      END DO

      RETURN
    END SUBROUTINE kronecker_product

