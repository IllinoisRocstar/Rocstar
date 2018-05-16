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
! Purpose: Collect relations for static and total internal energy per unit
!   mass.
!
! Description: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: MixtPerf_E.F90,v 1.3 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
!******************************************************************************

FUNCTION MixtPerf_Eo_DGPUVW(D,G,P,U,V,W)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: D,G,P,U,V,W
  REAL(RFREAL) :: MixtPerf_Eo_DGPUVW
   
  MixtPerf_Eo_DGPUVW = P/(D*(G - 1.0_RFREAL)) + 0.5_RFREAL*(U*U + V*V + W*W)

END FUNCTION MixtPerf_Eo_DGPUVW

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_Eo_DGPVm(D,G,P,Vm)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: D,G,P,Vm
  REAL(RFREAL) :: MixtPerf_Eo_DGPVm
   
  MixtPerf_Eo_DGPVm = P/(D*(G - 1.0_RFREAL)) + 0.5_RFREAL*Vm*Vm

END FUNCTION MixtPerf_Eo_DGPVm

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_Eo_GRTUVW(G,R,T,U,V,W)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: G,R,T,U,V,W
  REAL(RFREAL) :: MixtPerf_Eo_GRTUVW
   
  MixtPerf_Eo_GRTUVW = R*T/(G - 1.0_RFREAL) + 0.5_RFREAL*(U*U + V*V + W*W)

END FUNCTION MixtPerf_Eo_GRTUVW

!******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtPerf_E.F90,v $
! Revision 1.3  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:48:54  haselbac
! Initial revision after changing case
!
! Revision 1.2  2004/04/01 21:26:20  haselbac
! Added MixtPerf_E_GRTUVW
!
! Revision 1.1  2002/05/04 16:16:52  haselbac
! Initial revision
!
!******************************************************************************






