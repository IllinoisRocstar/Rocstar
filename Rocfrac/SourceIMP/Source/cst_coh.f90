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
      SUBROUTINE CST_COH(coor,lmc,matc,R_co,d,deltan,deltat, &
          sigmax,taumax,Sinit,Sthresh,nstep,numnp,numco,numat_coh, &
          delta, numploadelem, idpressload, pressload, &
          numboundmesh, idmesh, rmesh, nboundtype,num_fail_coh,  &
          index_fail_coh, npress, nmesh, ielem_coh, iface_coh)

      IMPLICIT NONE

!-----global variables
      INTEGER :: nstep                 ! time step
      INTEGER :: idstrt                ! time step to start output
      INTEGER :: idint                 ! time interval for output
      INTEGER :: numnp                 ! number of nodal points
      INTEGER :: numco                 ! number of CST cohesive elements
      INTEGER :: numat_coh             ! number of cohesive materials
      INTEGER :: start_ccst            ! number of cohesive elements 
                                       ! on border
      INTEGER :: iflag

!-----S/F interaction variables (ALE explicit code)
      INTEGER :: npress
      INTEGER :: nmesh
      INTEGER :: numploadelem 
      INTEGER :: numboundmesh
      INTEGER :: num_fail_coh
      INTEGER, DIMENSION(1:5,1:npress) :: idpressload
      REAL*8 , DIMENSION(1:3,1:npress) :: pressload
      INTEGER, DIMENSION(1:4,1:nmesh)  :: idmesh
      REAL*8 , DIMENSION(1:3,1:nmesh)  :: rmesh
      INTEGER, DIMENSION(1:numnp) :: nboundtype
      INTEGER, DIMENSION(1:numco) :: index_fail_coh
      INTEGER, DIMENSION(1:2,1:numco) :: ielem_coh
      INTEGER, DIMENSION(1:2,1:numco) :: iface_coh

! --  global coordinates
      REAL*8, DIMENSION(1:3,1:numnp) :: coor
! --  connectivity table for cohesive el
      INTEGER, DIMENSION(1:6,1:numco) :: lmc
! --  mat number for cohesive element
      INTEGER, DIMENSION(1:numco) :: matc
! --  cohesive traction loads
      real*8,DIMENSION(1:3*numnp) :: R_co
! --  nodal displacements
      REAL*8, DIMENSION(1:3*numnp) :: d
! --  the threshold, then damage of edge
      REAL*8, DIMENSION(1:3,1:numco) :: Sthresh
! --  normal characteristic length
      REAL*8, DIMENSION(1:numat_coh) :: deltan
! --  tangential characteristic length
      REAL*8, DIMENSION(1:numat_coh) :: deltat
! --  maximum normal stress
      REAL*8, DIMENSION(1:numat_coh) :: sigmax
! --  maximum normal stress
      REAL*8, DIMENSION(1:numat_coh) :: taumax
! --  initial Sthresh
      REAL*8, DIMENSION(1:numat_coh) :: Sinit
      REAL*8 :: delta              ! time step
!-----local variables
      REAL*8 :: deltn, deltt    ! the char length of current element
      REAL*8 :: sig             ! max normal stress of element
      REAL*8 :: tau             ! max shearing stress of element
      REAL*8 :: area            ! area of face
      REAL*8 :: Dx1,Dy1,Dz1     ! nd 1 openings in the x,y,z direction
      REAL*8 :: Dx2,Dy2,Dz2     ! nd 2 openings in the x,y,z direction
      REAL*8 :: Dx3,Dy3,Dz3     ! nd 3 openings in the x,y,z direction
      REAL*8 :: Dn1,Dn2,Dn3     ! openings in the normal direction
      REAL*8 :: Dt1,Dt2,Dt3     ! openings in the tan direction
      REAL*8 :: Dn,Dt           ! opening at the sample point
      REAL*8 :: Tn,Tt           ! cohesive traction at sample point
      REAL*8 :: Tx1,Tx2,Tx3     ! cohesive traction at sample point
      REAL*8 :: Ty1,Ty2,Ty3     ! cohesive traction at sample point
      REAL*8 :: Tz1,Tz2,Tz3     ! cohesive traction at sample point
      REAL*8 :: Rx1,Rx2,Rx3     ! cohesive forces at nodes 1,2,3  
      REAL*8 :: Ry1,Ry2,Ry3     ! cohesive forces at nodes 1,2,3
      REAL*8 :: Rz1,Rz2,Rz3     ! cohesive forces at nodes 1,2,3
! --  openings in tagent directions
      REAL*8 :: Dt_x, Dt_y, Dt_z
! --  Components of the tangent vector
      REAL*8 :: xtangent_1,ytangent_1, ztangent_1
      REAL*8 :: xtangent_2,ytangent_2, ztangent_2
      REAL*8 :: xtangent_3,ytangent_3, ztangent_3
! --  Components of two vectors parallel to plane
      REAL*8 :: vec1_x,vec1_y,vec1_z,vec2_x,vec2_y,vec2_z
! --  Componets of the normal vector of the face
      REAL*8 :: xnormal,ynormal,znormal
! --  Norm of the normal
      REAL*8 :: xnorm
! --  
 
      REAL*8 :: g1,g2           ! gauss quadrature points
      REAL*8 :: w1              ! weights
      INTEGER :: m              ! material number of current element
      INTEGER :: n1,n2,n3,n4,n5,n6 ! 3 times the node number of element
      INTEGER :: i              ! loop counter
      INTEGER :: j              ! loop counter
      INTEGER :: iaux
      REAL*8 :: x               ! dummy
      parameter (g1 = 0.666666666666667, &
          g2 = 0.166666666666667, &
          w1 = 0.333333333333333)
!-----The cohesive law
      DO i = 1,numco
      IF (index_fail_coh(i) .eq. 0) then

         m = matc(i)
         deltn = 1.d0/deltan(m)
         deltt = 1.d0/deltat(m)
         sig = sigmax(m)
         tau = taumax(m)
         n1 = lmc(1,i)*3        ! n1 = node number 1 (*3)
         n2 = lmc(2,i)*3        ! n2 = node number 2 (*3)
         n3 = lmc(3,i)*3        ! n3 = node number 3 (*3)
         n4 = lmc(4,i)*3        ! n4 = node number 4 (*3)
         n5 = lmc(5,i)*3        ! n5 = node number 5 (*3)
         n6 = lmc(6,i)*3        ! n6 = node number 6 (*3)

! --- Calculate the normal
!
!     Components of two vectors parallel to plane
         vec1_x = coor(1,lmc(3,i)) - coor(1,lmc(1,i))
         vec1_y = coor(2,lmc(3,i)) - coor(2,lmc(1,i))
         vec1_z = coor(3,lmc(3,i)) - coor(3,lmc(1,i))
         vec2_x = coor(1,lmc(2,i)) - coor(1,lmc(1,i))
         vec2_y = coor(2,lmc(2,i)) - coor(2,lmc(1,i))
         vec2_z = coor(3,lmc(2,i)) - coor(3,lmc(1,i))
!     Take the cross product
         xnormal = vec1_y*vec2_z - vec2_y*vec1_z
         ynormal = - (vec1_x*vec2_z - vec2_x*vec1_z)
         znormal = vec1_x*vec2_y - vec2_x*vec1_y
!     Length of the normal
         xnorm = SQRT(xnormal**2 + ynormal**2 + znormal**2)

!     Normalizing the normal
!             _      _
!    _       Vec1 x Vec2
!    n = ------------------
!         || Vec1 x Vec2 ||
!
         xnormal = xnormal/xnorm
         ynormal = ynormal/xnorm
         znormal = znormal/xnorm
         
         Dx1 = d(n5-2) - d(n1-2)
         Dy1 = d(n5-1) - d(n1-1)
         Dz1 = d(n5)   - d(n1)
!
!                  _      _
!     COD    =    COD dot n
!        n
!
         Dn1 = Dx1*xnormal + Dy1*ynormal + Dz1*znormal
!
!      _           _           _
!     COD    =    COD  -  COD  n
!        t                   n
!
         Dt_x = Dx1 - Dn1*xnormal
         Dt_y = Dy1 - Dn1*ynormal
         Dt_z = Dz1 - Dn1*znormal
!
!                    _
!     COD    =   || COD || 
!        t             t
!
         Dt1 = SQRT( Dt_x**2 + Dt_y**2 + Dt_z**2)
!                    _ 
!                   COD
!                      t
!     Tangent = t = ----
!                   COD
!                      t
!
         IF(Dt1.EQ.0.d0)THEN
            xtangent_1 = 0.d0
            ytangent_1 = 0.d0
            ztangent_1 = 0.d0
         ELSE
            xtangent_1 = Dt_x / Dt1
            ytangent_1 = Dt_y / Dt1
            ztangent_1 = Dt_z / Dt1
         ENDIF
         
         Dt1 =  Dt1*deltt       ! 
         Dn1 =  Dn1*deltn       ! normalized openings
         
         Dx2 = d(n6-2) - d(n2-2)
         Dy2 = d(n6-1) - d(n2-1)
         Dz2 = d(n6)   - d(n2)
         
         Dn2 = Dx2*xnormal + Dy2*ynormal + Dz2*znormal
         
         Dt_x = Dx2 - Dn2*xnormal
         Dt_y = Dy2 - Dn2*ynormal
         Dt_z = Dz2 - Dn2*znormal
         Dt2 = SQRT( Dt_x**2 + Dt_y**2 + Dt_z**2)
         
         IF(Dt2.EQ.0.d0)THEN
            xtangent_2 = 0.d0
            ytangent_2 = 0.d0
            ztangent_2 = 0.d0
         ELSE
            xtangent_2 = Dt_x / Dt2
            ytangent_2 = Dt_y / Dt2
            ztangent_2 = Dt_z / Dt2
         ENDIF 
         
         Dt2 = Dt2*deltt
         Dn2 = Dn2*deltn
         
         Dx3 = d(n4-2) - d(n3-2)
         Dy3 = d(n4-1) - d(n3-1)
         Dz3 = d(n4)   - d(n3)
         Dn3 = Dx3*xnormal + Dy3*ynormal + Dz3*znormal
         
         Dt_x = Dx3 - Dn3*xnormal
         Dt_y = Dy3 - Dn3*ynormal
         Dt_z = Dz3 - Dn3*znormal
         Dt3 = SQRT( Dt_x**2 + Dt_y**2 + Dt_z**2)
         
         IF(Dt3.EQ.0.d0)THEN
            xtangent_3 = 0.d0
            ytangent_3 = 0.d0
            ztangent_3 = 0.d0
         ELSE
            xtangent_3 = Dt_x / Dt3
            ytangent_3 = Dt_y / Dt3
            ztangent_3 = Dt_z / Dt3
         ENDIF
            
         Dt3 = Dt3*deltt  
         Dn3 = Dn3*deltn 

!----gauss point 1 
         Dt = g1*Dt1 + g2*Dt2 + g2*Dt3
         Dn = g1*Dn1 + g2*Dn2 + g2*Dn3
         x = 1.d0 - SQRT( Dn*Dn + Dt*Dt )
         IF (x.LE.0.d0) THEN
            Sthresh(1,i) = 0.d0
         ELSEIF (x.LE.Sthresh(1,i)) THEN
            Sthresh(1,i) = x
         ENDIF
         IF (Dn.GT.0.d0) THEN
            Tn = Sthresh(1,i)/(1.d0-Sthresh(1,i))*sig*Dn
         ELSE
            Tn = Sinit(m)/(1.d0-Sinit(m))*sig*Dn
         ENDIF
         Tt = Sthresh(1,i)/(1.d0-Sthresh(1,i))*tau*Dt
!     _      _          _
!     T = T  n   +   T  t
!          n          t
!
         Tx1 = Tn*xnormal + Tt*xtangent_1
         Ty1 = Tn*ynormal + Tt*ytangent_1
         Tz1 = Tn*znormal + Tt*ztangent_1
!----gauss point 2
         Dt = g2*Dt1 + g1*Dt2 + g2*Dt3
         Dn = g2*Dn1 + g1*Dn2 + g2*Dn3
         x = 1.d0 - SQRT( Dn*Dn + Dt*Dt )
         IF (x.LE.0.d0) THEN
            Sthresh(2,i) = 0.d0
         ELSEIF (x.LE.Sthresh(2,i)) THEN
            Sthresh(2,i) = x
         ENDIF
         IF (Dn.GT.0.d0) THEN
            Tn = Sthresh(2,i)/(1.d0-Sthresh(2,i))*sig*Dn
         ELSE
            Tn = Sinit(m)/(1.d0-Sinit(m))*sig*Dn
         ENDIF
         Tt = Sthresh(2,i)/(1.d0-Sthresh(2,i))*tau*Dt
!     _      _          _
!     T = T  n   +   T  t
!          n          t
!         
         Tx2 = Tn*xnormal + Tt*xtangent_2
         Ty2 = Tn*ynormal + Tt*ytangent_2
         Tz2 = Tn*znormal + Tt*ztangent_2
!----gauss point 3
         Dt = g2*Dt1 + g2*Dt2 + g1*Dt3
         Dn = g2*Dn1 + g2*Dn2 + g1*Dn3
         x = 1.d0 - SQRT( Dn*Dn + Dt*Dt )
         IF (x.LE.0.d0) THEN
            Sthresh(3,i) = 0.d0
         ELSEIF (x.LE.Sthresh(3,i)) THEN
            Sthresh(3,i) = x
         ENDIF
         IF (Dn.GT.0.d0) THEN
            Tn = Sthresh(3,i)/(1.d0-Sthresh(3,i))*sig*Dn
         ELSE
            Tn = Sinit(m)/(1.d0-Sinit(m))*sig*Dn
         ENDIF
         Tt = Sthresh(3,i)/(1.d0-Sthresh(3,i))*tau*Dt

!     _      _          _
!     T = T  n   +   T  t
!          n          t
!         
         Tx3 = Tn*xnormal + Tt*xtangent_3
         Ty3 = Tn*ynormal + Tt*ytangent_3
         Tz3 = Tn*znormal + Tt*ztangent_3
!
!     Area of the triangle is half the area of the parallelogram 
!     determined by the cross product of the vectors P12 and P13
         
         area = 0.5d0*xnorm
         
         Rx1 = area*w1*(g1*Tx1+g2*Tx2+g2*Tx3)
         Ry1 = area*w1*(g1*Ty1+g2*Ty2+g2*Ty3)
         Rz1 = area*w1*(g1*Tz1+g2*Tz2+g2*Tz3)
         Rx2 = area*w1*(g2*Tx1+g1*Tx2+g2*Tx3)
         Ry2 = area*w1*(g2*Ty1+g1*Ty2+g2*Ty3)
         Rz2 = area*w1*(g2*Tz1+g1*Tz2+g2*Tz3)
         Rx3 = area*w1*(g2*Tx1+g2*Tx2+g1*Tx3)
         Ry3 = area*w1*(g2*Ty1+g2*Ty2+g1*Ty3)
         Rz3 = area*w1*(g2*Tz1+g2*Tz2+g1*Tz3)
         
         R_co(n1-2) = R_co(n1-2) + Rx1
         R_co(n1-1) = R_co(n1-1) + Ry1
         R_co(n1)   = R_co(n1)   + Rz1
         R_co(n5-2) = R_co(n5-2) - Rx1
         R_co(n5-1) = R_co(n5-1) - Ry1
         R_co(n5)   = R_co(n5)   - Rz1
         R_co(n2-2) = R_co(n2-2) + Rx2
         R_co(n2-1) = R_co(n2-1) + Ry2
         R_co(n2)   = R_co(n2)   + Rz2
         R_co(n6-2) = R_co(n6-2) - Rx2
         R_co(n6-1) = R_co(n6-1) - Ry2
         R_co(n6)   = R_co(n6)   - Rz2
         R_co(n3-2) = R_co(n3-2) + Rx3
         R_co(n3-1) = R_co(n3-1) + Ry3
         R_co(n3)   = R_co(n3)   + Rz3
         R_co(n4-2) = R_co(n4-2) - Rx3
         R_co(n4-1) = R_co(n4-1) - Ry3
         R_co(n4)   = R_co(n4)   - Rz3
!
!        find out free surface as a result of crack propagation 
!
!------  Update mesh BC conditions

         IF ((Sthresh(1,i) .le. 1.0E-3) .AND.&
            (Sthresh(2,i) .le. 1.0E-3) .AND.&
            (Sthresh(3,i) .le. 1.0E-3) ) Then

            IF (index_fail_coh(i) .eq. 0) then
                num_fail_coh = num_fail_coh + 1
                index_fail_coh(i) = 1

                print *,' cohesive element ',i,' failed'
                print *,' total number of failed coh elems ',&
                         num_fail_coh

!------------ Solid Propellant part
              numploadelem  = numploadelem + 1
              idpressload(1,numploadelem) = ielem_coh(1,i)
              idpressload(2,numploadelem) = iface_coh(1,i)
              idpressload(3,numploadelem) = lmc(1,i)
              idpressload(4,numploadelem) = lmc(2,i)
              idpressload(5,numploadelem) = lmc(3,i)
              pressload(1,numploadelem) = 1.0d0
              pressload(2,numploadelem) = 1.0d0
              pressload(3,numploadelem) = 1.0d0

!------------ Case part 
              numploadelem  = numploadelem + 1
!-- be careful of node numbering !
              idpressload(1,numploadelem) = ielem_coh(2,i)
              idpressload(2,numploadelem) = iface_coh(2,i)
              idpressload(3,numploadelem) = lmc(5,i)
              idpressload(4,numploadelem) = lmc(4,i)
              idpressload(5,numploadelem) = lmc(6,i)
              pressload(1,numploadelem) = 1.0d0
              pressload(2,numploadelem) = 1.0d0
              pressload(3,numploadelem) = 1.0d0

!------------ store newly-generated crack face subjected to burning

              do j = 1, 3	! only for SP 
                 numboundmesh = numboundmesh + 1
                 idmesh(1,numboundmesh) =  lmc(j,i)
                 if(nboundtype(lmc(j,i)) .eq. 8) then ! corner node
                    idmesh(2,numboundmesh) =  0
                    idmesh(3,numboundmesh) =  0
                    idmesh(4,numboundmesh) =  1
                    rmesh( 1,numboundmesh) =  1.0d0
                    rmesh( 2,numboundmesh) = -1.0d0 
                    rmesh( 3,numboundmesh) =  0.0d0 
                 else 
                    idmesh(2,numboundmesh) =  1
                    idmesh(3,numboundmesh) =  0
                    idmesh(4,numboundmesh) =  1
                    rmesh( 1,numboundmesh) =  0.0d0
                    rmesh( 2,numboundmesh) = -1.0d0 
                    rmesh( 3,numboundmesh) =  0.0d0 
                 endif
              enddo

            ENDIF 
         ENDIF 

      ENDIF
      ENDDO

      print *,' total number of failed coh elems ',num_fail_coh

      RETURN
      END

