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
!------------------------------------------------------GET_MAT_STIFFNESS
SUBROUTINE get_mat_stiffness(e,dnu,dmat)
! *--------------------------------------------------------------------*
! |                                                                    |
! |    Returns the material stiffness matrix, given the material       |
! |    properties (FOR LINEAR ELASTICITY ONLY).                        |
! |                                                                    |
! |    <dmat>      : 9 x 9 matrix of the material tensor               |
! |    <props>     : vector that contains material properties          |
! |                  props(1) : Young's Modulus                        |
! |                  props(2) : Poisson's Ratio                        |
! |                                                                    |
! *--------------------------------------------------------------------*

  IMPLICIT NONE

  REAL*8 , DIMENSION(1:9,1:9) :: dmat
  REAL*8  :: e, dnu
  REAL*8  :: d1111,d1122,d1212
  
!  Initialize coefficients

  d1111 = e * (1.d0-dnu) / ((1.d0+dnu)*(1.d0-2.d0*dnu) )
  d1122 = e * dnu / ((1.d0+dnu)*(1.d0-2.d0*dnu) )
  d1212 = e / ((1.d0+dnu)*2.d0)

      
!  Compute <mat_c>
      
  dmat(1:9,1:9) = 0.d0
  
  dmat(1,1) = d1111
  dmat(1,5) = d1122
  dmat(1,9) = d1122
  
  dmat(2,2) = d1212
  dmat(2,4) = d1212
  
  dmat(3,3) = d1212
  dmat(3,7) = d1212
  
  dmat(4,2) = d1212
  dmat(4,4) = d1212
  
  dmat(5,1) = d1122
  dmat(5,5) = d1111
  dmat(5,9) = d1122
  
  dmat(6,6) = d1212
  dmat(6,8) = d1212
  
  dmat(7,3) = d1212
  dmat(7,7) = d1212
  
  dmat(8,6) = d1212
  dmat(8,8) = d1212
  
  dmat(9,1) = d1122
  dmat(9,5) = d1122
  dmat(9,9) = d1111
      
  RETURN
END SUBROUTINE get_mat_stiffness

