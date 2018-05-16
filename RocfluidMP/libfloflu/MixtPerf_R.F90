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
! Purpose: Collect relations for gas constant.
!
! Description: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: MixtPerf_R.F90,v 1.3 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

FUNCTION MixtPerf_R_CpG(Cp,G)

  USE ModDataTypes
  
  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: Cp,G
  REAL(RFREAL) :: MixtPerf_R_CpG
  
  MixtPerf_R_CpG = Cp - Cp/G  

END FUNCTION MixtPerf_R_CpG

!------------------------------------------------------------------------------

FUNCTION MixtPerf_R_M(M)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: M
  REAL(RFREAL) :: MixtPerf_R_M
   
  MixtPerf_R_M = 8314.3_RFREAL/M

END FUNCTION MixtPerf_R_M

!******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtPerf_R.F90,v $
! Revision 1.3  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:49:25  haselbac
! Initial revision after changing case
!
! Revision 1.3  2003/03/15 16:20:58  haselbac
! Simplified computation, eliminates diffs for || computations...
!
! Revision 1.2  2002/06/05 18:30:50  haselbac
! Added function mixtPerf_R_CpG
!
! Revision 1.1  2002/05/04 16:16:52  haselbac
! Initial revision
!
!******************************************************************************






