! *********************************************************************
! * Rocstar Simulation Suite                                          *
! * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
! *                                                                   *
! * Illinois Rocstar LLC                                              *
! * Champaign, IL                                                     *
! * www.illinoisrocstar.com                                           *
! * sales@illinoisrocstar.com                                         *
! *                                                                   *
! * License: See LICENSE file in top level of distribution package or *
! * http://opensource.org/licenses/NCSA                               *
! *********************************************************************
! *********************************************************************
! * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
! * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
! * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
! * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
! * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
! * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
! * Arising FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
! * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
! *********************************************************************
!******************************************************************************
!
! Purpose: grid generation based on the linear transfinite interpolation (TFI).
!
! Description: file contains the following routines:
!              - RFLO_Tfint1d = 1-D TFI
!              - RFLO_Tfint2d = 2-D TFI
!
! Input for Tfint1d: s   = normalized arclength
!                    p1  = coordinates of first edge point
!                    p2  = coordinates of last edge point.
!
! Input for Tfint2d: s1-4 = normalized arclength along edges 1-4
!                    e1-4 = coordinates along edges 1-4
!                    p1-4 = coordinates at corners 1-4.
!
! Output: xyz = interpolated point.
!
! Notes: the edges and corners for Tfint2d are enumbered as follows:
!
!          4    4 (j=je)    3
!           ----------------
!           |              |
!  1 (i=is) |              | 2 (i=ie)
!           |              |
!           ----------------
!          1    3 (j=js)    2
!
!******************************************************************************
!
! $Id: RFLO_Tfint.F90,v 1.3 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_Tfint1d( s,p1,p2,xyz )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: s, p1(3), p2(3), xyz(3)

!******************************************************************************

  xyz(:) = (1._RFREAL-s)*p1(:) + s*p2(:)

END SUBROUTINE RFLO_Tfint1d

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_Tfint2d( s1,s2,s3,s4,e1,e2,e3,e4,p1,p2,p3,p4,xyz )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: s1, s2, s3, s4
  REAL(RFREAL) :: e1(3), e2(3), e3(3), e4(3)
  REAL(RFREAL) :: p1(3), p2(3), p3(3), p4(3)
  REAL(RFREAL) :: xyz(3)

! ... local variables
  REAL(RFREAL) :: a1, a2, phij1, phij2, phii1, phii2

!******************************************************************************

  a1     = s1 - s2
  a2     = s3 - s4
  phij1  = (s1-a1*s3)/(1._RFREAL-a1*a2)
  phij2  = 1._RFREAL - phij1
  phii1  = (s3-a2*s1)/(1._RFREAL-a1*a2)
  phii2  = 1._RFREAL - phii1
  xyz(:) =       phij2*e3(:) +       phij1*e4(:) + &
                 phii2*e1(:) +       phii1*e2(:) - &
           phij2*phii2*p1(:) - phij1*phii1*p3(:) - &
           phij2*phii1*p2(:) - phij1*phii2*p4(:)

END SUBROUTINE RFLO_Tfint2d

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_Tfint.F90,v $
! Revision 1.3  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:21  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:17  wasistho
! lower to upper case
!
! Revision 1.1  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
!******************************************************************************






