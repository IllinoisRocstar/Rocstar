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
! Purpose: Collect relations for static and total speed of sound for perfect
!   gases.
!
! Description: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: MixtPerf_C.F90,v 1.3 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

FUNCTION MixtPerf_C_Co2GUVW(Co2,G,U,V,W)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: Co2,G,U,V,W
  REAL(RFREAL) :: MixtPerf_C_Co2GUVW
   
  MixtPerf_C_Co2GUVW = & 
    SQRT(Co2 - 0.5_RFREAL*(G - 1.0_RFREAL)*(U*U + V*V + W*W))

END FUNCTION MixtPerf_C_Co2GUVW

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_C_DGP(D,G,P)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: D,G,P
  REAL(RFREAL) :: MixtPerf_C_DGP
   
  MixtPerf_C_DGP = SQRT(G*P/D)

END FUNCTION MixtPerf_C_DGP

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_C_GHoVm2(G,Ho,Vm2)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: G,Ho,Vm2
  REAL(RFREAL) :: MixtPerf_C_GHoVm2
   
  MixtPerf_C_GHoVm2 = SQRT((G - 1.0_RFREAL)*(Ho - 0.5_RFREAL*Vm2))

END FUNCTION MixtPerf_C_GHoVm2

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_C_GRT(G,R,T)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: G,R,T
  REAL(RFREAL) :: MixtPerf_C_GRT
   
  MixtPerf_C_GRT = SQRT(G*R*T)

END FUNCTION MixtPerf_C_GRT

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_C2_GRT(G,R,T)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: G,R,T
  REAL(RFREAL) :: MixtPerf_C2_GRT
   
  MixtPerf_C2_GRT = G*R*T

END FUNCTION MixtPerf_C2_GRT

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_Co2_CGUVW(C,G,U,V,W)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: C,G,U,V,W
  REAL(RFREAL) :: MixtPerf_Co2_CGUVW
   
  MixtPerf_Co2_CGUVW = C*C + 0.5_RFREAL*(G - 1.0_RFREAL)*(U*U + V*V + W*W)

END FUNCTION MixtPerf_Co2_CGUVW

!******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtPerf_C.F90,v $
! Revision 1.3  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:48:47  haselbac
! Initial revision after changing case
!
! Revision 1.2  2002/05/28 13:44:44  haselbac
! Added new functions
!
! Revision 1.1  2002/05/04 16:16:52  haselbac
! Initial revision
!
!******************************************************************************






