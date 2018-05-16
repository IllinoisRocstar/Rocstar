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
! Purpose: Add source terms to the residual.
!
! Description: None.
!
! Input: 
!   region      data of current region.
!
! Output: None.
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: SourceTerms.F90,v 1.6 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SourceTerms(region)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  
  IMPLICIT NONE

#ifdef RFLO
#include "Indexing.h"
#endif

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region) :: region

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: ibc,iec,ic
#ifdef RFLO
  INTEGER :: idcbeg,idcend,jdcbeg,jdcend,kdcbeg,kdcend
  INTEGER :: iLev,iCOff,ijCOff
#endif
  REAL(RFREAL) :: fac,fvol,rhoVol
  REAL(RFREAL), POINTER :: cv(:,:),cvOld(:,:),sDual(:,:),rhs(:,:),vol(:)
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'SourceTerms',&
  'SourceTerms.F90' )

! ******************************************************************************
! Get dimensions and pointers
! ******************************************************************************

#ifdef RFLO
  iLev = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)
  
  cv    => region%levels(iLev)%mixt%cv
  cvOld => region%levels(iLev)%mixt%cvOld
  sDual => region%levels(iLev)%mixt%sDual
  rhs   => region%levels(iLev)%mixt%rhs
  vol   => region%levels(iLev)%grid%vol
#endif

#ifdef RFLU
  ibc = 1
  iec = region%grid%nCells

  cv  => region%mixt%cv
  rhs => region%mixt%rhs
  vol => region%grid%vol 
#endif

#ifdef RFLO
! ******************************************************************************
! Source term due to dual time-stepping
! ******************************************************************************

  IF (global%solverType==SOLV_IMPLICIT .AND. global%dualTstSource) THEN
    fac = 1.5_RFREAL/global%dtMin
    DO ic=ibc,iec
      fvol = fac*vol(ic)
      rhs(CV_MIXT_DENS,ic) = rhs(CV_MIXT_DENS,ic) - sDual(CV_MIXT_DENS,ic) + &
                             fvol*cv(CV_MIXT_DENS,ic)
      rhs(CV_MIXT_XMOM,ic) = rhs(CV_MIXT_XMOM,ic) - sDual(CV_MIXT_XMOM,ic) + &
                             fvol*cv(CV_MIXT_XMOM,ic)
      rhs(CV_MIXT_YMOM,ic) = rhs(CV_MIXT_YMOM,ic) - sDual(CV_MIXT_YMOM,ic) + &
                             fvol*cv(CV_MIXT_YMOM,ic)
      rhs(CV_MIXT_ZMOM,ic) = rhs(CV_MIXT_ZMOM,ic) - sDual(CV_MIXT_ZMOM,ic) + &
                             fvol*cv(CV_MIXT_ZMOM,ic)
      rhs(CV_MIXT_ENER,ic) = rhs(CV_MIXT_ENER,ic) - sDual(CV_MIXT_ENER,ic) + &
                             fvol*cv(CV_MIXT_ENER,ic)
    ENDDO
  ENDIF
#endif

! ******************************************************************************
! Source term due to acceleration
! ******************************************************************************

  IF ( global%accelOn .EQV. .TRUE. ) THEN
    DO ic = ibc,iec
      rhoVol = vol(ic)*cv(CV_MIXT_DENS,ic)
      rhs(CV_MIXT_XMOM,ic) = rhs(CV_MIXT_XMOM,ic) - global%accelX*rhoVol
      rhs(CV_MIXT_YMOM,ic) = rhs(CV_MIXT_YMOM,ic) - global%accelY*rhoVol
      rhs(CV_MIXT_ZMOM,ic) = rhs(CV_MIXT_ZMOM,ic) - global%accelZ*rhoVol
      rhs(CV_MIXT_ENER,ic) = rhs(CV_MIXT_ENER,ic) - &
                             vol(ic)*(global%accelX*cv(CV_MIXT_XMOM,ic)+ &
                                      global%accelY*cv(CV_MIXT_YMOM,ic)+ &
                                      global%accelZ*cv(CV_MIXT_ZMOM,ic))
    END DO ! ic
  END IF ! global%accelOn

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SourceTerms

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SourceTerms.F90,v $
! Revision 1.6  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.3  2005/11/10 01:58:37  haselbac
! Added Rocflu support for accel terms, clean-up
!
! Revision 1.2  2004/12/04 06:13:31  wasistho
! cvOld to cv in updating rhs by DualTST source terms
!
! Revision 1.1  2004/12/01 16:51:24  haselbac
! Initial revision after changing case
!
! Revision 1.9  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/08/28 20:05:38  jblazek
! Added acceleration terms.
!
! Revision 1.5  2003/07/03 21:48:44  jblazek
! Implemented dual-time stepping.
!
! Revision 1.4  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.3  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
!
! *****************************************************************************







