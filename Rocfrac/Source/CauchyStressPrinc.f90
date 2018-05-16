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
SUBROUTINE CauchyStressPrinc(ndime,xmu,xlamb,detf,stret,princ,sigma,sprin)

!!****f* Rocfrac/Rocfrac/Source/CauchyStressPrinc.f90
!!
!!  NAME
!!     CauchyStressPrinc
!!
!!  FUNCTION
!!
!!     Determines the Cauchy stress tensor for materials defined in
!!     principal directions
!!
!!  INPUTS
!!
!!     ndime  -->  number of dimensions
!!     xmu    -->  mu coefficient
!!     xlamb  -->  lambda coefficient
!!     detf   -->  determinant of F, i.e. J
!!     stret  -->  vector containing the stretches
!!     btens  -->  matrix containing the the principal directions
!!
!!  OUTPUT
!!     sigma  -->  Cauchy stress tensor
!!     sprin  -->  principal stresses
!!
!!****

  IMPLICIT DOUBLE PRECISION (a-h,o-z)
  DIMENSION sigma(3,3),princ(3,3),stret(3),sprin(3)
      
!
!     Obtains the principal stress using equation (5.103) (5.90)
!
  detfinv = 1.d0/detf

  DO id=1,ndime
     sprin(id) = detfinv*(xlamb*LOG(detf) + 2.d0*xmu*LOG(stret(id)))
  ENDDO
!
!     Obtains the Cartesian Cauchy stress tensor using equation (5.47)
!
  DO jd=1,ndime      
     DO kd=1,ndime
        sum=0.0d0
        DO id=1,ndime
           sum=sum+sprin(id)*princ(jd,id)*princ(kd,id)
        ENDDO
        sigma(jd,kd)=sum
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE CauchyStressPrinc

