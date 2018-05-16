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
! Purpose: Collect relations for molar mass.
!
! Description: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: MixtPerf_M.F90,v 1.3 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

FUNCTION MixtPerf_M_R(R)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: R
  REAL(RFREAL) :: MixtPerf_M_R
   
  MixtPerf_M_R = 8314.3_RFREAL/R

END FUNCTION MixtPerf_M_R

!******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtPerf_M.F90,v $
! Revision 1.3  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:49:11  haselbac
! Initial revision after changing case
!
! Revision 1.1  2002/06/05 18:30:20  haselbac
! Initial revision
!
!******************************************************************************






