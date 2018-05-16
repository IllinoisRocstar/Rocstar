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
! ******************************************************************************
!
! Purpose: Collect relations for gas specific heat at constant volume.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: MixtPerf_Cv.F90,v 1.3 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

FUNCTION MixtPerf_Cv_CpR(Cp,R)

  USE ModDataTypes
  
  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: Cp,R
  REAL(RFREAL) :: MixtPerf_Cv_CpR
  
  MixtPerf_Cv_CpR = Cp - R  

END FUNCTION MixtPerf_Cv_CpR

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtPerf_Cv.F90,v $
! Revision 1.3  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2006/03/26 20:20:55  haselbac
! Initial revision
!
!
! ******************************************************************************






