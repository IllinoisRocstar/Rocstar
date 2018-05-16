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
! Purpose: Collect relations for mixture pressure for gas-liquid model.
!
! Description: None.
!
! Notes: 
!   1. Use Dz instead of Do because of conflict with DO keyword in Fortran.
!
! ******************************************************************************
!
! $Id: MixtGasLiq_P.F90,v 1.3 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

FUNCTION MixtGasLiq_P(DYl,DYv,DYg,Cl2,Cv2,Cg2,D,Dz,Po,To,Bp,Bt,T)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: Bp,Bt,Cg2,Cl2,Cv2,D,Dz,DYg,DYl,DYv,Po,T,To
  REAL(RFREAL) :: MixtGasLiq_P

  REAL(RFREAL) :: term1,term2
     
  term1 = Cl2*(Dz + Bt*(T-To) - Bp*Po - DYl) - (DYv*Cv2 + DYg*Cg2)
  term2 = 4.0_RFREAL*Cl2*(Dz + Bt*(T-To) - Bp*Po)*(DYv*Cv2 + DYg*Cg2)

  MixtGasLiq_P = 0.5_RFREAL*(-term1 + SQRT(term1*term1 + term2))

END FUNCTION MixtGasLiq_P

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtGasLiq_P.F90,v $
! Revision 1.3  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2006/03/26 20:20:50  haselbac
! Initial revision
!
! ******************************************************************************







