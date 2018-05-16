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
! Purpose: Collect relations for static and total density for perfect gases.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: MixtPerf_D.F90,v 1.4 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************


FUNCTION MixtPerf_D_CGP(C,G,P)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: C,G,P
  REAL(RFREAL) :: MixtPerf_D_CGP
   
  MixtPerf_D_CGP = G*P/(C*C)

END FUNCTION MixtPerf_D_CGP

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_D_DoGMa(Do,G,Ma)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: Do,G,Ma
  REAL(RFREAL) :: MixtPerf_D_DoGMa
   
  MixtPerf_D_DoGMa = Do/ &
    (1.0_RFREAL + 0.5_RFREAL*(G-1.0_RFREAL)*Ma*Ma)**(1.0_RFREAL/(G-1.0_RFREAL))

END FUNCTION MixtPerf_D_DoGMa

! -----------------------------------------------------------------------------

FUNCTION MixtPerf_D_PRT(P,R,T)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: P,R,T 
  REAL(RFREAL) :: MixtPerf_D_PRT
   
  MixtPerf_D_PRT = P/(R*T)

END FUNCTION MixtPerf_D_PRT

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtPerf_D.F90,v $
! Revision 1.4  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/03/15 20:43:17  haselbac
! Added MixtPerf_D_CGP
!
! Revision 1.1  2004/12/01 16:48:50  haselbac
! Initial revision after changing case
!
! Revision 1.2  2003/09/16 01:23:45  haselbac
! Added function MixtPerf_D_DoGMa
!
! Revision 1.1  2002/05/04 16:16:52  haselbac
! Initial revision
!
! ******************************************************************************






