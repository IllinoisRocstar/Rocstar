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
! Purpose: compute gas variables for a thermally and calorically
!          perfect gas.
!
! Description: none.
!
! Input: inBeg    = first value to update
!        inEnd    = last value to update
!        indCp    = indicates if cp varies over cells (=1) or is constant (=0)
!        indMol   = indicates if the mol mass varies over cells (=1) or is
!                   constant (=0)
!        refCp    = reference value of specific heat at const. pressure
!        refGamma = reference value of the specific heat ratio
!        cv       = conservative variables
!
! Output: gv = gas variables (cp, molecular mass)
!
! Notes: it is assumed that the values are constant and given by the
!        reference values.
!
!******************************************************************************
!
! $Id: PerfgasGasVars.F90,v 1.3 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PerfgasGasVars( inBeg,inEnd,indCp,indMol,refCp,refGamma,cv,gv )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModParameters
  USE ModInterfaces, ONLY: MixtPerf_M_R, MixtPerf_R_CpG
  IMPLICIT NONE

! ... parameters
  INTEGER :: inBeg, inEnd, indCp, indMol

  TYPE(t_global), POINTER :: global

  REAL(RFREAL)          :: refCp, refGamma
  REAL(RFREAL), POINTER :: cv(:,:), gv(:,:)

! ... loop variables
  INTEGER :: ic

! ... local variables
  REAL(RFREAL) :: rgas

!******************************************************************************

  rgas = MixtPerf_R_CpG( refCp,refGamma )

  DO ic=inBeg,inEnd
    gv(GV_MIXT_MOL,ic*indMol) = MixtPerf_M_R( rgas )
    gv(GV_MIXT_CP ,ic*indCp ) = refCp
  ENDDO

END SUBROUTINE PerfgasGasVars

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PerfgasGasVars.F90,v $
! Revision 1.3  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:50:00  haselbac
! Initial revision after changing case
!
! Revision 1.8  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/07/05 23:20:46  jblazek
! Corrected bug in perfgasDependentVars.F90; did some cosmetics.
!
! Revision 1.3  2002/06/05 18:33:39  haselbac
! Converted to use of mixtPerf routines
!
! Revision 1.2  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.1  2002/01/10 00:02:06  jblazek
! Added calculation of mixture properties.
!
!******************************************************************************






