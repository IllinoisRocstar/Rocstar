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
! Purpose: Collect relations for mixture speed of sound for gas-liquid model.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: MixtGasLiq_C.F90,v 1.3 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

FUNCTION MixtGasLiq_C(Cvm,D,P,Dl,Dv,Dg,VFl,VFv,VFg,Cl2,Cv2,Cg2,Bl2,Bv2,Bg2)

   USE ModDataTypes

   IMPLICIT NONE

   REAL(RFREAL), INTENT(IN) :: Bg2,Bl2,Bv2,Cl2,Cv2,Cg2,Cvm,D,Dg,Dl,Dv,P,VFg, & 
                               VFl,VFv
   REAL(RFREAL) :: MixtGasLiq_C
   
   REAL(RFREAL) ::  denom,numer,term1,term2,term3

   term1 = Bl2*VFl*Dv*Cv2*Dg*Cg2
   term2 = Bv2*VFv*Dl*Cl2*Dg*Cg2
   term3 = Bg2*VFg*Dl*Cl2*Dv*Cv2

   numer  = D*Cvm*Dl*Cl2*Dv*Cv2*Dg*Cg2 + P*(term1 + term2 + term3)
   denom  = D*D*Cvm*(VFl*Dv*Cv2*Dg*Cg2 + VFv*Dl*Cl2*Dg*Cg2 + VFg*Dl*Cl2*Dv*Cv2)

   MixtGasLiq_C = SQRT(numer/denom)

END FUNCTION MixtGasLiq_C

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtGasLiq_C.F90,v $
! Revision 1.3  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2006/03/26 20:20:46  haselbac
! Initial revision
!
! ******************************************************************************






