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
      SUBROUTINE v3d4n_mass( &
           coor, &               ! Mesh Coordinates
           lmcstet, &            ! Connectivity array
           matcstet, &           ! Element Material ID
           rho, &                ! Density 
           xm, &                 ! mass matrix
           numnp, &
           numcstet, &
           numat_vol, &
           nstart, &
           nend,TotalMass,NumElNeigh,ElConn,Alpha,Vhat)

!-----Performs displacement based computations for 2-d, triangular
!-----linear elastic elements with linear interpolation functions.
!-----(constant strain triangles)
!-----The mass matrixes surface tractions and the strain and stress 
!-----calculations are all done here for the CSTs.
      IMPLICIT NONE
      INTEGER :: numnp          ! number of nodes
      INTEGER :: numcstet       ! number of CSTets
      INTEGER :: numat_vol      ! number of materials
!--   densities
      REAL*8, DIMENSION(1:numat_vol) :: rho
!--   reciprical of mass matrix diagonal
      REAL*8, DIMENSION(1:numnp) :: xm
!--   material number for CSTet element
      INTEGER, DIMENSION(1:numcstet) :: matcstet
!--   connectivity table for CSTet elem
      INTEGER, DIMENSION(1:4,1:numcstet) :: lmcstet
!--   global coordinates
      REAL*8, DIMENSION(1:3,1:numnp) :: coor
      integer, DIMENSION(1:numnp) ::NumElNeigh
      INTEGER, DIMENSION(1:numnp,1:40) :: ElConn ! fix 40 should not be hard coded
      REAL*8 :: aa                 ! determinant of jacobian (2*area)
      REAL*8 :: x !,x1,x2,x3         ! dummy variable
      INTEGER :: m                 ! current element's material number
      INTEGER :: n1,n2,n3,n4,n5,n6 ! nodes, and dummy vars
      INTEGER :: n7,n8,n9,n10
      INTEGER :: i,j,nstart,nend     ! loop counter
!--  Coordinate holding variable
      REAL*8 :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
!--  Coordinate subtractions
      REAL*8 :: x14, x24, x34, y14, y24, y34, z14, z24, z34
!--  
      REAL*8 :: c11, c21, c31
!--   6*volume and the volume      
      REAL*8 :: Vx6,volume
      REAL*8 :: TotalMass,VolEl, aux
      real*8, dimension(1:4,1:numcstet) :: alpha
      integer :: iElNum
      real*8 :: sumx
      REAL*8, DIMENSION(1:numnp) :: Vhat


      DO i = 1, numnp
         sumx = 0.d0
         DO j = 1, NumElNeigh(i)
            iElNum = ElConn(i,j)
            m = matcstet(iElNum)

            n1 = lmcstet(1,iElNum)
            n2 = lmcstet(2,iElNum)
            n3 = lmcstet(3,iElNum)
            n4 = lmcstet(4,iElNum) 

            x1 = coor(1,n1)
            x2 = coor(1,n2)
            x3 = coor(1,n3)
            x4 = coor(1,n4)
            y1 = coor(2,n1)
            y2 = coor(2,n2)
            y3 = coor(2,n3)
            y4 = coor(2,n4)
            z1 = coor(3,n1)
            z2 = coor(3,n2)
            z3 = coor(3,n3)
            z4 = coor(3,n4)
         
            x14 = x1 - x4
            x24 = x2 - x4
            x34 = x3 - x4
            y14 = y1 - y4
            y24 = y2 - y4
            y34 = y3 - y4
            z14 = z1 - z4
            z24 = z2 - z4
            z34 = z3 - z4

            c11 =    y24*z34 - z24*y34
            c21 = -( x24*z34 - z24*x34 )
            c31 =    x24*y34 - y24*x34

            Vx6 = -( x14*c11 + y14*c21 + z14*c31 )
          
            VolEl = Vx6/6.d0

            IF( VolEl.LT.0.d0)THEN
               PRINT*,'ERROR'
               PRINT*,'NEG, Volume... STOPPING'
               PRINT*,'  ELEMENT =',i
               PRINT*,'  NODES=',n1,n2,n3,n4
               PRINT*,'  x-Coordinates:',x1,x2,x3,x4
               PRINT*,'  y-Coordinates:',y1,y2,y3,y4
               PRINT*,'  z-Coordinates:',z1,z2,z3,z4
               STOP
            ENDIF
          
            IF(n1.EQ.i)THEN
               aux =  alpha(1,iElNum)*VolEl
               sumx = sumx + aux
!               print*,aux,rho(m)
               xm(n1) = xm(n1) + aux*rho(m)
               TotalMass = TotalMass + aux*rho(m)
            ELSE IF(n2.EQ.i)THEN
               aux = alpha(2,iElNum)*VolEl
               sumx = sumx + aux 
!               print*,aux,rho(m)
               xm(n2) = xm(n2) + aux*rho(m)
               TotalMass = TotalMass + aux*rho(m)
            ELSE IF(n3.EQ.i)THEN
               aux = alpha(3,iElNum)*VolEl
               sumx = sumx +  aux 
!               print*,aux,rho(m)
               xm(n3) = xm(n3) + aux*rho(m)
               TotalMass = TotalMass + aux*rho(m)
            ELSE IF(n4.EQ.i)THEN
               aux = alpha(4,iElNum)*VolEl
               sumx = sumx + aux
!               print*,aux,rho(m)
               xm(n4) = xm(n4) + aux*rho(m) 
               TotalMass = TotalMass + aux*rho(m)
            ENDIF
         ENDDO 
         Vhat(i) = sumx
      ENDDO

! fix need to change if have more then one material, I.e. don't want to figure in the case

      RETURN
      END

