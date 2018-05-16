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
! Purpose: two point averaging of vector variable
!
! Description: vector variable fvar is averaged from two reference values
!              at location ijkN1 and ijkN2, and the result is stored as
!              fvar at location ijkD
!
! Input: ijkD  = location where averaged vector is stored
!        ijkN1 = location of first referenced vector
!        ijkN2 = location of second referenced vector
!        iFBeg = begin index of vector variable component
!        iFEnd = end index of vector variable component
!        fvar  = vector variable array to be averaged from ijkN1 and ijkN2
!
! Output: fvar = averaged vector variable stored at ijkD
!
! Notes: this routine is meant as service routine for RFLO_CopyEdgeFaceNorm
!        and RFLO_CopyEdgeFaceParal in rocflo/, but can be employed by
!        other routines if needed 
!
!******************************************************************************
!
! $Id: AverageVecVar.F90,v 1.3 2008/12/06 08:44:08 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE AverageVecVar( ijkD,ijkN1,ijkN2,iFBeg,iFEnd,fvar )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  INTEGER               :: ijkD, ijkN1, ijkN2, iFBeg, iFEnd
  REAL(RFREAL), POINTER :: fvar(:,:)

!******************************************************************************

  fvar(iFBeg:iFEnd,ijkD) = 0.5_RFREAL*(fvar(iFBeg:iFEnd,ijkN1)+ &
                                       fvar(iFBeg:iFEnd,ijkN2))

END SUBROUTINE AverageVecVar

!******************************************************************************
!
! RCS Revision history:
!
! $Log: AverageVecVar.F90,v $
! Revision 1.3  2008/12/06 08:44:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:47:51  haselbac
! Initial revision after changing case
!
! Revision 1.4  2002/08/30 01:47:58  jblazek
! Added support for moving grids.
!
! Revision 1.3  2002/08/22 01:45:53  jblazek
! And also removed ModError ...
!
! Revision 1.2  2002/08/22 01:43:04  jblazek
! Removed RegisterFunction and DeregisterFunction to speed things up.
!
! Revision 1.1  2002/07/30 02:41:34  wasistho
! promoted from contained routine to own file
!
!******************************************************************************






