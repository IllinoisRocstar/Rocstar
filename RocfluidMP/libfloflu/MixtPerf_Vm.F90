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
! Purpose: Collect relations for velocity magnitude for perfect gases.
!
! Description: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: MixtPerf_Vm.F90,v 1.3 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

FUNCTION MixtPerf_Vm_C2Co2G(C2,Co2,G)

  USE ModDataTypes

  IMPLICIT NONE
  
  REAL(RFREAL), INTENT(IN) :: C2,Co2,G
  REAL(RFREAL) :: MixtPerf_Vm_C2Co2G
   
  IF ( Co2 > C2 ) THEN  
    MixtPerf_Vm_C2Co2G = SQRT(2.0_RFREAL/(G - 1.0_RFREAL)*(Co2 - C2))
  ELSE 
    MixtPerf_Vm_C2Co2G = 0.0_RFREAL
  END IF ! Co2

END FUNCTION MixtPerf_Vm_C2Co2G

!******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtPerf_Vm.F90,v $
! Revision 1.3  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:49:36  haselbac
! Initial revision after changing case
!
! Revision 1.3  2003/12/04 03:23:00  haselbac
! Added check against taking SQRT of negative number
!
! Revision 1.2  2002/05/28 13:46:10  haselbac
! Corrected name of MixtPerf_Vm_C2Co2
!
! Revision 1.1  2002/05/04 16:16:52  haselbac
! Initial revision
!
!******************************************************************************






