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
! Purpose: Collect relations for static and total temperature for perfect 
!   gases.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: MixtPerf_T.F90,v 1.6 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

FUNCTION MixtPerf_T_CGR(C,G,R)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: C,G,R
  REAL(RFREAL) :: MixtPerf_T_CGR
   
  MixtPerf_T_CGR = C*C/(G*R)

END FUNCTION MixtPerf_T_CGR

!------------------------------------------------------------------------------

FUNCTION MixtPerf_T_CpHoVm2(Cp,Ho,Vm2)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: Cp,Ho,Vm2
  REAL(RFREAL) :: MixtPerf_T_CpHoVm2
   
  MixtPerf_T_CpHoVm2 = (Ho-0.5_RFREAL*Vm2)/Cp

END FUNCTION MixtPerf_T_CpHoVm2

!------------------------------------------------------------------------------

FUNCTION MixtPerf_T_CvEoVm2(Cv,Eo,Vm2)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: Cv,Eo,Vm2
  REAL(RFREAL) :: MixtPerf_T_CvEoVm2
   
  MixtPerf_T_CvEoVm2 = (Eo-0.5_RFREAL*Vm2)/Cv

END FUNCTION MixtPerf_T_CvEoVm2

! ------------------------------------------------------------------------------

FUNCTION MixtPerf_T_DPR(D,P,R)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: D,P,R
  REAL(RFREAL) :: MixtPerf_T_DPR
   
  MixtPerf_T_DPR = P/(D*R)

END FUNCTION MixtPerf_T_DPR

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_T_GMaTo(G,Ma,To)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: G,Ma,To
  REAL(RFREAL) :: MixtPerf_T_GMaTo
   
  MixtPerf_T_GMaTo = To/(1.0_RFREAL + 0.5_RFREAL*(G - 1.0_RFREAL)*Ma*Ma)

END FUNCTION MixtPerf_T_GMaTo

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_To_CpTUVW(Cp,T,U,V,W)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: Cp,T,U,V,W
  REAL(RFREAL) :: MixtPerf_To_CpTUVW
   
  MixtPerf_To_CpTUVW = T + 0.5_RFREAL*(U*U + V*V + W*W)/Cp

END FUNCTION MixtPerf_To_CpTUVW

!******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtPerf_T.F90,v $
! Revision 1.6  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/05/01 20:58:15  haselbac
! Added MixtPerf_T_CpHoVm2 function
!
! Revision 1.3  2006/03/26 20:21:17  haselbac
! Added fuction
!
! Revision 1.2  2005/07/14 21:39:56  haselbac
! Added MixtPerf_To_CpTUVW
!
! Revision 1.1  2004/12/01 16:49:30  haselbac
! Initial revision after changing case
!
! Revision 1.3  2002/06/10 21:16:27  haselbac
! Added MixtPerf_T_GMaTo function
!
! Revision 1.2  2002/06/05 18:31:14  haselbac
! Added function mixtPerf_T_DPR
!
! Revision 1.1  2002/05/04 16:16:52  haselbac
! Initial revision
!
! ******************************************************************************






