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

!!****
!!
!!  NAME
!!     implicit_v3d8_me_K
!!
!!  FUNCTION
!!     This subroutine constructs the stiffness matrix associated
!!     with all of the elements that are in this processors .inp file.
!!
!!  INPUTS
!!     coor -- Coordinates of the nodes
!!     ElConnVol -- The connectivity table
!!     dmat -- Material stiffness matrix
!!     numnp -- Number of nodes
!!     nstart -- First element to compute the stiffness matrix for
!!     nend -- Last element to compute the stiffness matrix for
!!     NumEl -- Number of elements
!!     MatType -- Mapping of materials to elements
!!     NumMatType -- Total number of materials
!!     enhanced_map -- The enhanced map
!!     mixed_map -- The mixed map
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE implicit_v3d8_me_k(coor,ElConnVol,ci, & 
     numnp,nstart,nend,NumEL,MatType,NumMatType,enhanced_map,mixed_map,NumMatVol) 
 
 
  ! dmat is the 9x9 material matrix 
 
  USE comp_row_global 
  USE implicit_global 
  USE Precision 
 
  IMPLICIT NONE 
 
!-----Global variables 
  INTEGER :: numnp          ! number of nodes 
  INTEGER :: numat_vol      ! number of volumetric materials 
  INTEGER :: NumEl       ! number of LSTets 
  integer :: NumMatType 
  INTEGER :: NumMatVol
!--   coordinate array 
  REAL*8, DIMENSION(1:3,1:numnp) :: coor 
!--   connectivity table for Brick  
  INTEGER, DIMENSION(1:8,1:NumEl) :: ElConnVol 
  INTEGER, DIMENSION(1:NumEl) :: MatType 
!--   elastic stiffness consts 
  REAL*8, DIMENSION(1:6,1:6,NumMatType) :: ci 
!--   internal force 
  REAL*8, DIMENSION(1:3*numnp) :: R_in 
!--   displacement vector 
  REAL*8, DIMENSION(1:3*numnp) :: disp 
!---- Local variables 
!--   node numbers 
  INTEGER :: n1,n2,n3,n4,n5,n6,n7,n8,n9,n10 
 
  INTEGER :: k1n1,k1n2,k1n3,k1n4,k1n5,k1n6,k1n7,k1n8 
  INTEGER :: k2n1,k2n2,k2n3,k2n4,k2n5,k2n6,k2n7,k2n8 
  INTEGER :: k3n1,k3n2,k3n3,k3n4,k3n5,k3n6,k3n7,k3n8 
 
  REAL*8, DIMENSION(1:3,1:8) :: coord 
  !REAL*8, DIMENSION(1:24) :: Udisp 
   
  REAL*8 :: element_volume 
  INTEGER :: ielem, imat 
  integer :: igpt 
  REAL*8 nNn(8), dn(8,3), jac(3,3), jacinv(3,3), & 
       dn20(20,3),nNn20(20),& 
       t(3,3), tinv(3,3), & 
       dmat(9,9), bmat(9,24), bmat2(1:8,1:9,1:24), & 
       grad(9,24), tmp38(3,8), coeff(8),& 
       tkront(9,9), tkrontinv(9,9),& 
       e_mixed(9,12), & 
       gsm(12,24), gbg(12,12), bmat_avg(9,24),& 
       bu(9,24), ba(9,9), bsm(12,24),&                                          ! print edge 
       kuu(24,24), kau(9,24), kaa(9,9),atol, baVec(1:9,1:9,1:8), buVec(1:9,1:24,1:8),& 
       tmp249(24,9),tmp91(9,1),kaa_copy(9,9),enh2(9),dnorm,bmatigs(1:8,1:9,1:24),& 
       stresstmp(9,1), tmp9(9),dmatinf(5,9,9) 
 
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
 
  REAL*8 :: zero = 0.0, one = 1.0, detj,det, three = 3. 
  LOGICAL error, debug 
  REAL*8 :: alpha2, alpha 
  REAL*8, DIMENSION(1:9) :: fenh 
  REAL*8, DIMENSION(1:9,1:9) :: e_enhanced 
  REAL*8 :: maxDisp, con1, TP, sum 
  REAL*8 :: ThreeEighth = 3./8. 
  REAL*8 :: Denh 
  INTEGER :: igauss 
  REAL*8, DIMENSION(8) :: xi, eta, zeta 
  REAL*8, DIMENSION(8) :: xiE, etaE, zetaE 
  REAL(KIND=8), DIMENSION(1:8,1:9,1:12) :: mixed_map 
  REAL(KIND=8), DIMENSION(1:8,1:9,1:9)  :: enhanced_map 
  REAL(KIND=8), DIMENSION(1:8,1:8) :: ShapeFun,ShapeFunE 
  INTEGER :: otdev,mcrd,nnode 
  INTEGER,parameter :: ngpts = 8 
  REAL(KIND=8) :: strainEnh(1:ngpts,1:9,1:NumEl) 
  !REAL*8 :: stress(1:ngpts,1:9,1:NumEl) 
  INTEGER :: i,k,j 
  INTEGER :: nstart, nend 
 
  REAL(KIND=8), DIMENSION(1:24,1:24) :: amatrx 
 
 
!  INTEGER :: GlNumNp 
!  INTEGER, ALLOCATABLE, DIMENSION(:) :: rp_k, cval_k 
!  REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: aval_k 
  INTEGER :: nnz_k 
  INTEGER :: m 
  INTEGER :: dof1, dof2 
  INTEGER :: ldof1, ldof2 
  INTEGER :: gdof1, gdof2 
  INTEGER :: counter 
  REAL(kind=wp) :: tempval
 
 
  ! Initialize the K matrix 
  CALL numKnnz(NumNp,NumEl,ElConnVol,nnz_k) 
  ALLOCATE(rp_k(1:3*GNumNp+1)) 
  ALLOCATE(cval_k(1:nnz_k)) 
  ALLOCATE(aval_k(1:nnz_k)) 
  CALL initK(NumNp,NumEl,ElConnVol,nnz_k,rp_k,cval_k,aval_k) 
  !print*,myid,' nnz_k = ',nnz_k 
 
  DO ielem = nstart, nend ! Loop over elements 
 
        imat = MatType(ielem)

        CALL map_ci2dmat(ci(1:6,1:6,imat),dmat)
 
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
         
!!$        Udisp(1) = disp(k1n1) 
!!$        Udisp(2) = disp(k2n1) 
!!$        Udisp(3) = disp(k3n1) 
!!$ 
!!$        Udisp(4) = disp(k1n2) 
!!$        Udisp(5) = disp(k2n2) 
!!$        Udisp(6) = disp(k3n2)  
!!$ 
!!$        Udisp(7) = disp(k1n3) 
!!$        Udisp(8) = disp(k2n3) 
!!$        Udisp(9) = disp(k3n3)  
!!$ 
!!$        Udisp(10) = disp(k1n4) 
!!$        Udisp(11) = disp(k2n4) 
!!$        Udisp(12) = disp(k3n4)  
!!$ 
!!$        Udisp(13) = disp(k1n5) 
!!$        Udisp(14) = disp(k2n5) 
!!$        Udisp(15) = disp(k3n5)   
!!$ 
!!$        Udisp(16) = disp(k1n6) 
!!$        Udisp(17) = disp(k2n6) 
!!$        Udisp(18) = disp(k3n6)    
!!$ 
!!$        Udisp(19) = disp(k1n7) 
!!$        Udisp(20) = disp(k2n7) 
!!$        Udisp(21) = disp(k3n7)    
!!$ 
!!$        Udisp(22) = disp(k1n8) 
!!$        Udisp(23) = disp(k2n8) 
!!$        Udisp(24) = disp(k3n8) 
 
 
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
 
              CALL TENSORMUL(jacinv,3,3,'t',dn,8,3,'t',tmp38,3,one,zero) 
              CALL kronecker_product(tmp38,3,8,dident,3,3,grad,9) 
              CALL TENSORMUL(sop,9,9,'n',grad,9,24,'n',bmat,9,one,zero) 
 
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
 
           gbg = 0.   
           gsm(1:12,1:24) = 0. 
 
!! 2nd LOOP OVER GAUSS POINTS 
 
 
           DO igpt = 1, 8                           ! LOOP over gauss points 
 
              alpha2 = coeff(igpt) 
 
              e_mixed(1:9,1:12) = mixed_map(igpt,1:9,1:12) 
              CALL TENSORMUL(e_mixed,9,12,'t',e_mixed,9,12,'n',gbg,12,one,one) 
 
 
              bmat(1:9,1:24) = bmat2(igpt,1:9,1:24) - bmat_avg(1:9,1:24)      
              CALL TENSORMUL (tkront,9,9,'t',bmat,9,24,'n',grad,9,one,zero) 
            
 
! Accumilate G small 
              CALL TENSORMUL(e_mixed,9,12,'t',grad,9,24,'n',gsm,12,alpha2,one) 
 
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
           kuu(1:24,1:24) = 0.                       ! initialize 
           kau(1:9, 1:24) = 0. 
           kaa(1:9, 1:9)  = 0. 
 
!  Compute mixed/enhanced strain-displacement matrices <bu> & <ba> at  
!  each gauss point 
 
 
           DO igpt = 1, 8   ! Start 3rd loop over gauss points  
            
              detj = coeff(igpt)  ! / wi(igpt) but wi(igpt) = 1. for all gauss points    
            
!             CALL get_mixed_map(e_mixed,ri,igpt) 
              e_mixed(1:9,1:12) = mixed_map(igpt,1:9,1:12) 
 
!             CALL get_enhanced_map(e_enhanced,ri,igpt) 
              e_enhanced(1:9,1:9) = enhanced_map(igpt,1:9,1:9) 
            
! Compute <bu> (9 x 24) (Eqn. 2.27) 
 
              CALL TENSORMUL(e_mixed,9,12,'n',bsm,12,24,'n',grad,9,one,zero) 
              CALL TENSORMUL(tkrontinv,9,9,'t',grad,9,24,'n',bu,9,one,zero) 
              bu(1:9,1:24) = bu(1:9,1:24) / detj + bmat_avg(1:9,1:24) 
 
! Compute <ba> (9 x 9) (Eqn. 2.31) 
 
              CALL TENSORMUL(tkrontinv,9,9,'t',e_enhanced,9,9,'n',ba,9,one,zero) 
 
! Compute the enhanced strain field 
!               gauss pt, component, element 
 
              ba(1:9,1:9) = ba(1:9,1:9) /detj ! for later 
 
! 
! Compute <kaa> 
! 
              alpha = detj ! coeff(igpt) 
 
              CALL TENSORMUL(dmat(:,:),9,9,'n',bu,9,24,'n',grad,9,one,zero) 
! Compute <kuu> 
              CALL TENSORMUL(bu,9,24,'t',grad,9,24,'n',kuu,24,detj,one) 
! Compute <kau> 
              CALL TENSORMUL(ba,9,9, 't',grad,9,24,'n',kau,9,detj,one)   
 
              CALL TENSORMUL(dmat(:,:),9,9,'n',ba,9,9,'n',tkront,9,one,zero) 
 
              CALL TENSORMUL(ba,9,9,'t',tkront,9,9,'n',kaa,9,detj,one) 
 
           END DO 
 
!  Compute <amatrx> 
 
      amatrx = kuu ! Static condensation 
!        -1 
! [ K   ] 
!    aa  
 
      call invert(kaa,9,det) 
 
      CALL TENSORMUL(kaa,9,9,'n',kau,9,24,'n',grad,9,one,zero) 
      call TENSORMUL(kau,9,24,'t',grad,9,24,'n',amatrx,24,-one,one) 
 
      DO i=1,24 
         DO j=1,24 
            atol = ABS(amatrx(i,j)-amatrx(j,i))  
            IF(atol.GE.1.0E-03) THEN 
               WRITE(7,*) '>>> amatrx is UNSYMMETRIC, fatal error', & 
                    ' i =',i,' j =',j, ' diff =', atol 
            END IF 
         END DO 
      END DO 
 
!************************************************************** 
! Assemble local K (i.e. amatrix) into Global stiffness matrix 
!************************************************************** 

      DO dof1 = 1, 3  ! dof1 loop
         DO dof2 = 1, 3  ! dof2 loop
            DO i = 1, 8  ! node1 loop
               gdof1 = 3 * Local2Global(ElConnVol(i,ielem)) - 3 + dof1
               ldof1 = 3 * i - 3 + dof1
               DO j = 1, 8  ! node2 loop
                  ldof2 = 3 * j - 3 + dof2
                  gdof2 = 3 * Local2Global(ElConnVol(j,ielem)) - 3 + dof2

                  ! Make sure the matrix is exactly symmetric by using only values in the upper triangle of amatrx
                  IF ( ldof1 <= ldof2 ) THEN
                     tempval = amatrx(ldof1,ldof2)
                  ELSE
                     tempval = amatrx(ldof2,ldof1)
                  ENDIF

                  ! Place the nonzeros into the K matrix
                  IF (tempval /= 0.0) THEN  
                     counter = 0
                     DO m = rp_k(gdof1)+1, rp_k(gdof1+1)
                        IF (cval_k(m) == gdof2-1) THEN
                           aval_k(m) = aval_k(m) + tempval
                           counter = 1
                        ENDIF
                     ENDDO
                     if(counter==0) print*,'WARNING:  Unable to add value to K matrix at (',gdof1,',',gdof2, &
                          ') on processor ',myid
                  ELSE
                     print*,myid,'ZERO DETECTED IN LOCAL STIFFNESS MATRIX!',gdof1,gdof2
                  END IF

               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
 
!  Compute enhanced field parameters <enh> 
 
!     enh = zero 
!      
!     do i=1, maxiters 
!      
!     call tensormul(kau,9,24,'n',u,ndofel,1,'n',fenh,9,one,zero) 
!     call tensormul(kaa_copy,9,9,'n',enh,9,1,'n',fenh,9,one,one) 
!      
!     call tensormul(kaa,9,9,'n',fenh,9,1,'n',denh,9,-one,zero) 
!      
!     enh = enh + denh 
!      
!     call tensormul(denh,9,1,'t',denh,9,1,'n',dnorm,1,one,zero) 
!      
!     write(7,*) ' >>> iterno = ',i,' dnorm =',dnorm 
!      
!     if(dnorm.le.0.0000001) then  
!     write(7,*) '>>> CONVERGED' 
!     call printmat(enh,9,1,otdev,' enh ') 
!     goto 1 
!     end if 
!      
!     end do 
!      
!     1    continue                  ! Convergence achieved or max. iters exhausted 
 
   ENDDO 
 
   ALLOCATE(rp_temp(1:3*GNumNp+1),cval_temp(1:nnz_k),aval_temp(1:nnz_k)) 
   nnz_temp = nnz_k 
   rp_temp = rp_k 
   cval_temp = cval_k 
   aval_temp = aval_k 
   DEALLOCATE(rp_k,cval_k,aval_k) 
 
 
 END SUBROUTINE implicit_v3d8_me_k


 SUBROUTINE map_ci2dmat(ci,dmat)


   REAL*8, DIMENSION(1:6,1:6) :: ci
   REAL*8, DIMENSION(1:9,1:9) :: dmat


   dmat(1,1) = ci(1,1)
   dmat(1,2) = ci(1,4)
   dmat(1,3) = ci(1,5)
   dmat(1,4) = ci(1,4)
   dmat(1,5) = ci(1,2)
   dmat(1,6) = ci(1,6)
   dmat(1,7) = ci(1,5)
   dmat(1,8) = ci(1,6)
   dmat(1,9) = ci(1,3)
   
   dmat(2,1) = ci(1,4)
   dmat(2,2) = ci(4,4)
   dmat(2,3) = ci(4,5)
   dmat(2,4) = ci(4,4)
   dmat(2,5) = ci(2,4)
   dmat(2,6) = ci(4,6)
   dmat(2,7) = ci(4,5)
   dmat(2,8) = ci(4,6)
   dmat(2,9) = ci(3,4)
   
   dmat(3,1) = ci(1,5)
   dmat(3,2) = ci(4,5)
   dmat(3,3) = ci(5,5)
   dmat(3,4) = ci(4,5)
   dmat(3,5) = ci(2,5)
   dmat(3,6) = ci(5,6)
   dmat(3,7) = ci(5,5)
   dmat(3,8) = ci(5,6)
   dmat(3,9) = ci(3,5)
   
   dmat(4,1) = ci(1,4)
   dmat(4,2) = ci(4,4)
   dmat(4,3) = ci(4,5)
   dmat(4,4) = ci(4,4)
   dmat(4,5) = ci(2,4)
   dmat(4,6) = ci(4,6)
   dmat(4,7) = ci(4,5)
   dmat(4,8) = ci(4,6)
   dmat(4,9) = ci(3,4)
   
   dmat(5,1) = ci(1,2)
   dmat(5,2) = ci(2,4)
   dmat(5,3) = ci(2,5)
   dmat(5,4) = ci(2,4)
   dmat(5,5) = ci(2,2)
   dmat(5,6) = ci(2,6)
   dmat(5,7) = ci(2,5)
   dmat(5,8) = ci(2,6)
   dmat(5,9) = ci(2,3)
   
   dmat(6,1) = ci(1,6)
   dmat(6,2) = ci(4,6)
   dmat(6,3) = ci(5,6)
   dmat(6,4) = ci(4,6)
   dmat(6,5) = ci(2,6)
   dmat(6,6) = ci(6,6)
   dmat(6,7) = ci(5,6)
   dmat(6,8) = ci(6,6)
   dmat(6,9) = ci(3,6)
   
   dmat(7,1) = ci(1,5)
   dmat(7,2) = ci(4,5)
   dmat(7,3) = ci(5,5)
   dmat(7,4) = ci(4,5)
   dmat(7,5) = ci(2,5)
   dmat(7,6) = ci(5,6)
   dmat(7,7) = ci(5,5)
   dmat(7,8) = ci(5,6)
   dmat(7,9) = ci(3,5)
   
   
   dmat(8,1) = ci(1,6)
   dmat(8,2) = ci(4,6)
   dmat(8,3) = ci(5,6)
   dmat(8,4) = ci(4,6)
   dmat(8,5) = ci(2,6)
   dmat(8,6) = ci(6,6)
   dmat(8,7) = ci(5,6)
   dmat(8,8) = ci(6,6)
   dmat(8,9) = ci(3,6)
   
   
   dmat(9,1) = ci(1,3)
   dmat(9,2) = ci(3,4)
   dmat(9,3) = ci(3,5)
   dmat(9,4) = ci(3,4)
   dmat(9,5) = ci(2,3)
   dmat(9,6) = ci(3,6)
   dmat(9,7) = ci(3,5)
   dmat(9,8) = ci(3,6)
   dmat(9,9) = ci(3,3)

 END SUBROUTINE map_ci2dmat

