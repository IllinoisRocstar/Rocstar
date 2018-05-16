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
SUBROUTINE jacobi(btens,stret,princ)

!!****f* Rocfrac/Rocfrac/Source/jacobi.f90
!!
!!  NAME
!!    jacobi
!!
!!  FUNCTION
!!
!!     Evaluates the stretches and principal directions given the b matrix
!!     using the Jacobi iteration. Adapted from numerical recpies
!!
!!  INPUTS
!!
!!     btens  -->  left Cauchy-Green tensor
!!
!!  OUTPUT
!!     stret  -->  vector containing the stretches
!!     princ  -->  matrix containing the three principal column vectors
!!
!!****

  IMPLICIT DOUBLE PRECISION (a-h,o-z)
  DIMENSION btens(3,3),stret(3),princ(3,3)
!
!
!     Initialise princ to the identity
!     
  DO i=1,3
     DO j=1,3
        princ(i,j)=0.d0
     ENDDO
     princ(i,i)=1.d0
     stret(i)=btens(i,i)
  ENDDO
!
!     Starts sweeping. 
!
      DO is=1,50
         sum=0.d0
         DO ip=1,2
            DO iq=ip+1,3
               sum=sum+ABS(btens(ip,iq))
            ENDDO
         ENDDO
!
!     IF the sum of off-diagonal terms is zero evaluates the 
!     stretches and returns
!
         IF(sum.EQ.0.d0) THEN
            DO i=1,3
               stret(i)=SQRT(stret(i))
            ENDDO
            RETURN
         ENDIF
!
!     Performs the sweep in three rotations. One per off diagonal term     
!
         DO ip=1,2
            DO iq=ip+1,3
               od=100.d0*ABS(btens(ip,iq))
               IF((od+ABS(stret(ip)).NE.ABS(stret(ip))).AND. &
                    (od+ABS(stret(iq)).NE.ABS(stret(iq)))) THEN
                  hd=stret(iq)-stret(ip)
!
!    Evaluates the rotation angle 
!
                  IF(ABS(hd)+od.EQ.ABS(hd)) THEN
                    t=btens(ip,iq)/hd
                  ELSE
                    theta=0.5d0*hd/btens(ip,iq)
                     t=1.d0/(ABS(theta)+SQRT(1.d0+theta**2))
                     IF(theta.LT.0.d0) t=-t
                  ENDIF
!
!     Re-evaluates the diagonal terms
!
                  c=1.d0/SQRT(1.d0+t**2)
                  s=t*c
                  tau=s/(1.d0+c)
                  h=t*btens(ip,iq)
                  stret(ip)=stret(ip)-h
                  stret(iq)=stret(iq)+h
!
!     Re-evaluates the remaining off-diagonal terms         
!
                  ir=6-ip-iq
                  g=btens(MIN(ir,ip),MAX(ir,ip))
                  h=btens(MIN(ir,iq),MAX(ir,iq))
                  btens(MIN(ir,ip),MAX(ir,ip))=g-s*(h+g*tau)
                  btens(MIN(ir,iq),MAX(ir,iq))=h+s*(g-h*tau)
!
!     Rotates the eigenvectors
!
                  DO ir=1,3
                     g=princ(ir,ip)
                     h=princ(ir,iq)
                     princ(ir,ip)=g-s*(h+g*tau)
                     princ(ir,iq)=h+s*(g-h*tau)
                  ENDDO
               ENDIF
               btens(ip,iq)=0.d0
            ENDDO
         ENDDO
      ENDDO
!
!     IF convergence is not achieved stops
!
      WRITE(6,100)
      STOP
 100  FORMAT(' Jacobi iteration unable to converge')
      END

