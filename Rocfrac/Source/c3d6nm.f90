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

SUBROUTINE c3d6nm( nsubn, nsubf, &
     nsubn2, nsubf2, &
     CoorOverlay, ElConnOverlay, &
     sd_subface_parents, &
     sd_subface_parents2, &
     sd_subface_mate, &
     sd_subface_nat_coors, &
     sd_subface_nat_coors_mate, &
     FaceOfVolEl, &
     FaceOfVolEl2, &
     NumNp, NumEl, &
     ElConn, &
     nf,nf2, &
     MapFaceEl2Vol, &
     MapFaceEl2Vol2, &
     deltan, deltat, sigmax, taumax, Sinit, &
     Rnet, Disp, Sthresh, NumMatCoh,MatID,SignFlag)

  IMPLICIT NONE

!-----global variables
  INTEGER :: nsubn ! number of Overlay nodes
  INTEGER :: nsubf ! number of Overlay elements
  INTEGER :: nsubn2, nsubf2
  REAL*8  :: SignFlag ! -1 for top material, +1 for bottom material

  INTEGER :: nf,nf2

  REAL*8, DIMENSION(1:3,1:nsubn)  :: CoorOverlay         ! Overlay Coordinates
  INTEGER, DIMENSION(1:3,1:nsubf) :: ElConnOverlay ! Overlay Connectivity
  INTEGER, DIMENSION(1:nsubf)     :: sd_subface_parents
  INTEGER, DIMENSION(1:nsubf)     :: sd_subface_mate
  INTEGER, DIMENSION(1:nsubf2)    :: sd_subface_parents2
  REAL*4, DIMENSION(1:6,1:nsubf)  :: sd_subface_nat_coors
  REAL*4, DIMENSION(1:6,1:nsubf2)  :: sd_subface_nat_coors_mate 
  INTEGER, DIMENSION(1:nf)     :: MapFaceEl2Vol
  INTEGER, DIMENSION(1:nf)     :: FaceOfVolEL 
  INTEGER, DIMENSION(1:nf2)     :: MapFaceEl2Vol2
  INTEGER, DIMENSION(1:nf2)     :: FaceOfVolEL2

  REAL*8, DIMENSION(1:3,1:nsubf) :: Sthresh ! the threshold, then damage of face

  INTEGER :: NumNp      ! number of nodal points in the overlay mesh
  INTEGER :: NumEl      ! number of cohesive overlay elements
  INTEGER :: NumMatCoh  ! number of cohesive materials
  INTEGER, DIMENSION(1:4,1:NumEL) :: ElConn        ! connectivity tables for cohesive el


  INTEGER, DIMENSION(1:NumEL) :: Map2VolNd

  REAL*8 :: Rnet(3*NumNp)         ! Net force at each node(currently=Rcoh)
  REAL*8 :: disp(3*NumNp)         ! nodal displacements
  REAL*8 :: deltan(NumMatCoh)   ! normal characteristic length
  REAL*8 :: deltat(NumMatCoh)   ! tangential characteristic length
  REAL*8 :: sigmax(NumMatCoh)   ! maximum normal stress
  REAL*8 :: taumax(NumMatCoh)   ! maximum normal stress
  REAL*8 :: Sinit(NumMatCoh)    ! initial Sthresh
!-----local variables
  REAL*8 :: NN1,NN2,NN3         ! Numbers to input into NN
  REAL*8 :: deltn, deltt        ! the char length of current element
  REAL*8 :: sig                 ! max normal stress of element
  REAL*8 :: tau                 ! max shearing stress of element
  INTEGER :: MatID                  ! material number of current element
  INTEGER :: n1,n2,n3,n4,n5,n6           ! twice the node numbers of element
  INTEGER :: i,j,k              ! loop counters
  REAL*8 :: x                   ! dummy variables
  REAL*8 :: g(3),weight(3)      ! integration point locations and weights
  REAL*8, DIMENSION(1:3,1:18) :: NN           ! Shape function matrix
  REAL*8 :: dloc(18)            ! local displacements (u1x u1y u2x u2y ...)
  REAL*8 :: delta1(2)           ! normal and tangential displacement at gauss point
  REAL*8 :: Tn,Tt               ! normal and tangential tractions
  REAL*8 :: T(3)                ! contains Tx Ty Tz
  REAL*8, DIMENSION(1:18) :: rcohloc         ! local cohesive reaction
  REAL*8 :: v12(3),v13(3)       ! vectors from nodes 1 to 2 and 1 to 3
  REAL*8 :: normal(3)           ! normal vector to element-v12 cross v13
  REAL*8 :: tangential(3)       ! tangential vector
  REAL*8 :: magnitude           ! magnitude of a vector
  REAL*8 :: delta(3)            ! dx dy dz
  REAL*8 :: aa(18)              ! dummy variable
  REAL*8 :: area, norm
  INTEGER :: face1, face2
  INTEGER :: FaceEL

  INTEGER :: nd1_Overlay,nd2_Overlay,nd3_Overlay
  INTEGER ::  VolEL
  INTEGER :: side
  INTEGER :: ii,FaceEL_counterparts


  weight(1)=0.33333333333333d0         !
  weight(2)=0.33333333333333d0         ! Weights of each gauss point
  weight(3)=0.33333333333333d0         !

!-----The cohesive law


  DO i = 1, nsubf            ! Loop over all the cohesive elements


     deltn = 1.d0/deltan(MatID)  ! deltn and deltt are the same as
     deltt = 1.d0/deltat(MatID)  ! 1/dnc and 1/dtc

     sig = sigmax(MatID)         ! maximum sigma and tau before failure
     tau = taumax(MatID) 

     nd1_Overlay = ElConnOverlay(1,i)
     nd2_Overlay = ElConnOverlay(2,i)
     nd3_Overlay = ElConnOverlay(3,i)


! Build vectors from nodes 1 to 2 and 1 to 3

     DO j=1,3
        v12(j) = CoorOverlay(j,nd2_Overlay)-CoorOverlay(j,nd1_Overlay)
        v13(j) = CoorOverlay(j,nd3_Overlay)-CoorOverlay(j,nd1_Overlay)
     END DO

     magnitude = SQRT((v12(2)*v13(3)-v12(3)*v13(2))**2+ &
                &(v12(1)*v13(3)-v12(3)*v13(1))**2+ &
                &(v12(1)*v13(2)-v12(2)*v13(1))**2)

! Find vector normal to element


     normal(1) =  (v12(2)*v13(3)-v12(3)*v13(2))
     normal(2) = -(v12(1)*v13(3)-v12(3)*v13(1))
     normal(3) =  (v12(1)*v13(2)-v12(2)*v13(1))


!    Use the fact that the norm of the cross product vector
!    is the area of the parallelogram they form.  The triangle they
!    form has half that area.
!
!  Reference:
!
!    Adrian Bowyer and John Woodwark,
!    A Programmer's Geometry,
!    Butterworths, 1983.

     norm = SQRT( normal(1)*normal(1) + normal(2)*normal(2) + normal(3)*normal(3) )

     area = 0.5d0 * norm

     normal(1) = normal(1)/magnitude
     normal(2) = normal(2)/magnitude
     normal(3) = normal(3)/magnitude


     FaceEL = sd_subface_parents(i)

     VolEL = MapFaceEl2Vol(FaceEL)

     side = FaceOfVolEl(FaceEl)

     IF(side.EQ.1)THEN
        n1 = 1; n2 = 2; n3 = 3
     ELSE IF(side.eq.2)THEN
        n1 = 1; n2 = 2; n3 = 4
     ELSE IF(side.eq.3)THEN
        n1 = 2; n2 = 3; n3 = 4
     ELSE IF(side.eq.4)THEN
        n1 = 4; n2 = 3; n3 = 1
     ENDIF

     n1 = ElConn(n1,VolEL)          ! n1 = node number 1 = n1m
     n2 = ElConn(n2,VolEL)          ! n2 = node number 2 = n2m
     n3 = ElConn(n3,VolEL)          ! n3 = node number 3 = n3m

     
     FaceEL_counterparts = sd_subface_parents2(sd_subface_mate(i))

     VolEL = MapFaceEl2Vol2(FaceEL_counterparts)

     side = FaceOfVolEl2(FaceEl_counterparts)


     IF(side.EQ.1)THEN
        n4 = 1; n5 = 2; n6 = 3
     ELSE IF(side.eq.2)THEN
        n4 = 1; n5 = 2; n6 = 4
     ELSE IF(side.eq.3)THEN
        n4 = 2; n5 = 3; n6 = 4
     ELSE IF(side.eq.4)THEN
        n4 = 4; n5 = 3; n6 = 1
     ENDIF


     n4 = ElConn(n4,VolEL)          ! n4 = node number 4 = n1p
     n5 = ElConn(n5,VolEL)          ! n5 = node number 4 = n2p
     n6 = ElConn(n6,VolEL)          ! n6 = node number 4 = n3p

     dloc(1) = disp (n1*3-2)    !
     dloc(2) = disp (n1*3-1)    !
     dloc(3) = disp (n1*3)      !
     dloc(4) = disp (n2*3-2)    !
     dloc(5) = disp (n2*3-1)    !
     dloc(6) = disp (n2*3)      !
     dloc(7) = disp (n3*3-2)    !
     dloc(8) = disp (n3*3-1)    ! Get a local displacement vector from
     dloc(9) = disp (n3*3)      ! the global displacement vector
     dloc(10) = disp (n4*3-2)   !
     dloc(11) = disp (n4*3-1)   !
     dloc(12) = disp (n4*3)     !
     dloc(13) = disp (n5*3-2)   !
     dloc(14) = disp (n5*3-1)   !
     dloc(15) = disp (n5*3)     !
     dloc(16) = disp (n6*3-2)   !
     dloc(17) = disp (n6*3-1)   !
     dloc(18) = disp (n6*3)     !


! loop over the 3 gauss points

     rcohloc(1:18) = 0.d0

     ii = 1
     DO k= 1,3


! Define shape function values
        NN1 = 1.d0 - DBLE(sd_subface_nat_coors(ii,i))-DBLE(sd_subface_nat_coors(ii+1,i))  !
        NN2 = DBLE(sd_subface_nat_coors(ii,i))        !
        NN3 = DBLE(sd_subface_nat_coors(ii+1,i))

! Build NN

        NN(1:3,1:9) = 0.d0

        NN(1,1) = -NN1   !
        NN(2,2) = -NN1   !
        NN(3,3) = -NN1   !
        NN(1,4) = -NN2   !
        NN(2,5) = -NN2   !
        NN(3,6) = -NN2   ! Assign values to the shape function
        NN(1,7) = -NN3   ! matrix
        NN(2,8) = -NN3   !
        NN(3,9) = -NN3   !

! Define shape function values
        NN1 = 1.d0 - DBLE(sd_subface_nat_coors_mate(ii,FaceEL_counterparts))- &
             DBLE(sd_subface_nat_coors_mate(ii+1,FaceEL_counterparts))  !
        NN2 = DBLE(sd_subface_nat_coors_mate(ii,FaceEL_counterparts))        !
        NN3 = DBLE(sd_subface_nat_coors_mate(ii+1,FaceEL_counterparts))
        ii = ii + 2

! Build NN

        NN(1:3,10:18) = 0.d0

        NN(1,10) = NN1   !
        NN(2,11) = NN1   !
        NN(3,12) = NN1   !
        NN(1,13) = NN2   !
        NN(2,14) = NN2   !
        NN(3,15) = NN2   ! Assign values to the shape function
        NN(1,16) = NN3   ! matrix
        NN(2,17) = NN3   !
        NN(3,18) = NN3   !

        delta = MATMUL(NN,dloc)     ! dx dy dz at each gauss point

        ii = ii + 2

!        IF(SignFlag.EQ.-1.d0)PRINT*,delta


! Find delta1 <dn dt> at each gauss point

        delta1(1) = DOT_PRODUCT(delta,normal)

        DO j=1,3
           tangential(j) = delta(j)-delta1(1)*normal(j)
        END DO

        delta1(2) = SQRT(tangential(1)**2+tangential(2)**2+tangential(3)**2)

! Caluate S value at gauss point

        x = 1.d0-SQRT((delta1(1)*deltn)**2+(delta1(2)*deltt)**2)
        Sthresh(k,i)=MAX(0.d0,MIN(x,Sthresh(k,i)))

!        PRINT*,Sthresh(k,i)

                             ! deltn = 1/dnc, deltt = 1/dtc

! TRACTIONS

        IF (delta1(1)>=0.d0) THEN
           Tn=Sthresh(k,i)/(1.d0-Sthresh(k,i))*sig*delta1(1)*deltn/Sinit(MatID)
        ELSE
           Tn=sig*delta1(1)*deltn/(1.d0-Sinit(MatID))
        END IF

        Tt=Sthresh(k,i)/(1.d0-Sthresh(k,i))*tau*delta1(2)*deltt/Sinit(MatID)


! Build a traction vector <Tx Ty Tz>

        DO j=1,3
           T(j)=Tn*normal(j)+Tt*tangential(j)
        END DO


        aa = MATMUL(TRANSPOSE(NN),T)  !  to make following calc easier

! Build rcohloc
!! Check this calculation.  Should weight be in there?  should it be * 1/3?
        DO j=1,18
           rcohloc(j) = rcohloc(j) - area*aa(j)*weight(k)
        END DO

     END DO      ! End of looping over integration points
     

! Build Rnet


     Rnet(n1*3-2) = Rnet(n1*3-2) + rcohloc(1)
     Rnet(n1*3-1) = Rnet(n1*3-1) + rcohloc(2)
     Rnet(n1*3)   = Rnet(n1*3)   + rcohloc(3)
     Rnet(n2*3-2) = Rnet(n2*3-2) + rcohloc(4)
     Rnet(n2*3-1) = Rnet(n2*3-1) + rcohloc(5)
     Rnet(n2*3)   = Rnet(n2*3)   + rcohloc(6)
     Rnet(n3*3-2) = Rnet(n3*3-2) + rcohloc(7)
     Rnet(n3*3-1) = Rnet(n3*3-1) + rcohloc(8)
     Rnet(n3*3)   = Rnet(n3*3)   + rcohloc(9)


     Rnet(n4*3-2) = Rnet(n4*3-2) + rcohloc(10)
     Rnet(n4*3-1) = Rnet(n4*3-1) + rcohloc(11)
     Rnet(n4*3)   = Rnet(n4*3)   + rcohloc(12)
     Rnet(n5*3-2) = Rnet(n5*3-2) + rcohloc(13)
     Rnet(n5*3-1) = Rnet(n5*3-1) + rcohloc(14)
     Rnet(n5*3)   = Rnet(n5*3)   + rcohloc(15)
     Rnet(n6*3-2) = Rnet(n6*3-2) + rcohloc(16)
     Rnet(n6*3-1) = Rnet(n6*3-1) + rcohloc(17)
     Rnet(n6*3)   = Rnet(n6*3)   + rcohloc(18)


  END DO         ! End of looping over elements

END SUBROUTINE C3D6NM


SUBROUTINE crossprod_3d ( x1, y1, z1, x2, y2, z2, x3, y3, z3 )
!
!*******************************************************************************
!
!! CROSS_3D computes the cross product of two vectors in 3D.
!
!
!  Definition:
!
!    The cross product in 3D can be regarded as the determinant of the
!    symbolic matrix:
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
!
!  Modified:
!
!    19 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X1, Y1, Z1, X2, Y2, Z2, the coordinates of the vectors.
!
!    Output, real X3, Y3, Z3, the cross product vector.
!
  IMPLICIT NONE
!
  REAL*8 :: x1
  REAL*8 :: x2
  REAL*8 :: x3
  REAL*8 :: y1
  REAL*8 :: y2
  REAL*8 :: y3
  REAL*8 :: z1
  REAL*8 :: z2
  REAL*8 :: z3
!
  x3 = y1 * z2 - z1 * y2
  y3 = z1 * x2 - x1 * z2
  z3 = x1 * y2 - y1 * x2
  
  RETURN
END SUBROUTINE crossprod_3d


