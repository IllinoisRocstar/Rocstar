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
subroutine enhanced_elem_maps_hex(mixed_map,enhanced_map)

  implicit none

  integer :: igpt

  REAL*8, DIMENSION(1:8,1:9,1:12) :: mixed_map
  REAL*8, DIMENSION(1:8,1:9,1:9)  :: enhanced_map

  REAL*8, DIMENSION(1:8) :: xi, eta, zeta
  REAL*8, parameter :: one = 1.d0, three = 3.d0
  


  xi   = (/ one/SQRT(three), -one/SQRT(three), -one/SQRT(three), &
       one/SQRT(three),  one/SQRT(three), -one/SQRT(three), &
       -one/SQRT(three),  one/SQRT(three)     /)
  
  eta  = (/ one/SQRT(three),  one/SQRT(three), -one/SQRT(three), &
       -one/SQRT(three),  one/SQRT(three),  one/SQRT(three), &
       -one/SQRT(three), -one/SQRT(three)     /)
  
  zeta = (/ one/SQRT(three),  one/SQRT(three),  one/SQRT(three), &
       one/SQRT(three), -one/SQRT(three), -one/SQRT(three), &  
       -one/SQRT(three), -one/SQRT(three)     /)

!  Compute static things for gauss points

  DO igpt =1, 8

!  Compute mixed mapping matrix <mixed_map> at each gauss point

     mixed_map(igpt,1,1)  = eta(igpt)
     mixed_map(igpt,1,2)  = zeta(igpt)
     mixed_map(igpt,1,3)  = eta(igpt) * zeta(igpt)
     mixed_map(igpt,2,10) = zeta(igpt)
     mixed_map(igpt,3,12) = eta(igpt)
     mixed_map(igpt,4,10) = zeta(igpt)
     mixed_map(igpt,5,4)  = xi(igpt)        
     mixed_map(igpt,5,5)  = zeta(igpt)
     mixed_map(igpt,5,6)  = xi(igpt) * zeta(igpt)
     mixed_map(igpt,6,11) = xi(igpt)
     mixed_map(igpt,7,12) = eta(igpt)
     mixed_map(igpt,8,11) = xi(igpt)
     mixed_map(igpt,9,7)  = xi(igpt)        
     mixed_map(igpt,9,8)  = eta(igpt)
     mixed_map(igpt,9,9)  = xi(igpt) * eta(igpt)

!  Compute enhanced mapping matrix <enhanced_map> at each gauss point

     enhanced_map(igpt,1,1)  = xi(igpt)    
     enhanced_map(igpt,1,2)  = xi(igpt) * eta(igpt)
     enhanced_map(igpt,1,3)  = xi(igpt) * zeta(igpt)
     enhanced_map(igpt,5,4)  = eta(igpt)
     enhanced_map(igpt,5,5)  = eta(igpt) * zeta(igpt)
     enhanced_map(igpt,5,6)  = eta(igpt) * xi(igpt)
     enhanced_map(igpt,9,7)  = zeta(igpt)        
     enhanced_map(igpt,9,8)  = zeta(igpt) * eta(igpt)
     enhanced_map(igpt,9,9)  = zeta(igpt) * xi(igpt)

  END DO


end subroutine enhanced_elem_maps_hex

