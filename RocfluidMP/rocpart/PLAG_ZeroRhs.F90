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
! Purpose: zero out residuals.
!
! Description: none.
!
! Input: region = current region.
!
! Output: region%levels%plag%rhs
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_ZeroRhs.F90,v 1.3 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_zeroRhs( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModMixture, ONLY    : t_mixt_input
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! ... loop variables
  INTEGER :: iCont, iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

#ifdef RFLO
  INTEGER  :: iLev
#endif
  INTEGER  :: nCont, nPcls

  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pRhs

  TYPE(t_plag),   POINTER :: pPlag
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ZeroRhs.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global, 'PLAG_ZeroRhs',&
  'PLAG_ZeroRhs.F90' )

! Get dimensions --------------------------------------------------------------

#ifdef RFLO
  iLev  = region%currLevel
  nPcls = region%levels(iLev)%plag%nPcls
#endif
#ifdef RFLU
  nPcls = region%plag%nPcls
#endif
  nCont = region%plagInput%nCont

! Set pointers ----------------------------------------------------------------

#ifdef RFLO
  pPlag => region%levels(iLev)%plag
#endif
#ifdef RFLU
  pPlag => region%plag
#endif

  pRhs  => pPlag%rhs

! Zero out residuals ----------------------------------------------------------

  DO iPcls = 1, nPcls
    pRhs(:,iPcls) = 0.0_RFREAL
  END DO  ! iPcls

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_ZeroRhs

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ZeroRhs.F90,v $
! Revision 1.3  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:58:23  fnajjar
! Initial revision after changing case
!
! Revision 1.4  2004/03/25 21:16:43  jferry
! fixed Vapor Energy bug
!
! Revision 1.3  2004/02/06 21:18:20  fnajjar
! Initial Integration of Rocpart with Rocflu
!
! Revision 1.2  2003/01/16 20:27:58  f-najjar
! Removed iRegionGlobal
!
! Revision 1.1  2002/10/25 14:20:32  f-najjar
! Initial Import of Rocpart
!
!******************************************************************************







