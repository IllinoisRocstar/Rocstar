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
! Purpose: average cell centred values around a grid node.
!
! Description: there are the following functions:
!              - Aver    = averages one variable
!              - AverDiv = divides variable 1 by variable 2 and
!                          conducts the averaging.
!
! Input: 
!
! Output: Aver    = average value
!         AverDiv = average value divided by second variable.
!
! Notes: simple arithmetical averaging is used.
!
!******************************************************************************
!
! $Id: Average.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

DOUBLE PRECISION FUNCTION Aver( cell,iEq,var )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  INTEGER :: cell(8), iEq

  REAL(RFREAL), POINTER :: var(:,:)

!******************************************************************************

  Aver = 0.125_RFREAL * (var(iEq,cell(1))+var(iEq,cell(2))+ &
                         var(iEq,cell(3))+var(iEq,cell(4))+ &
                         var(iEq,cell(5))+var(iEq,cell(6))+ &
                         var(iEq,cell(7))+var(iEq,cell(8)))

END FUNCTION Aver

! #############################################################################

DOUBLE PRECISION FUNCTION AverDiv( cell,iEq1,var1,iEq2,var2 )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  INTEGER :: cell(8), iEq1, iEq2

  REAL(RFREAL), POINTER :: var1(:,:), var2(:,:)

!******************************************************************************

  AverDiv = 0.125_RFREAL * &
            (var1(iEq1,cell(1))/var2(iEq2,cell(1))+ &
             var1(iEq1,cell(2))/var2(iEq2,cell(2))+ &
             var1(iEq1,cell(3))/var2(iEq2,cell(3))+ &
             var1(iEq1,cell(4))/var2(iEq2,cell(4))+ &
             var1(iEq1,cell(5))/var2(iEq2,cell(5))+ &
             var1(iEq1,cell(6))/var2(iEq2,cell(6))+ &
             var1(iEq1,cell(7))/var2(iEq2,cell(7))+ &
             var1(iEq1,cell(8))/var2(iEq2,cell(8)))

END FUNCTION AverDiv

!******************************************************************************
!
! RCS Revision history:
!
! $Log: Average.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:29:09  fnajjar
! Initial revision after changing case
!
! Revision 1.1  2003/09/25 15:40:22  jferry
! Implented Rocsmoke post-processing
!
!
!******************************************************************************






