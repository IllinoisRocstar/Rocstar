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
! Purpose: set dummycells to zero.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: regions%levels%radi%rhs = residual in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_FlimZeroDummyCells.F90,v 1.4 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FlimZeroDummyCells( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_ZeroDummyCells
#endif
#ifdef RFLU
  USE ModInterfaces, ONLY: RFLU_ZeroVirtualCellVars
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), TARGET :: region

! ... local variables
#ifdef RFLO
  INTEGER :: iLev
  REAL(RFREAL), POINTER   :: rhs(:,:)
#endif
  TYPE(t_global), POINTER :: global
#ifdef RFLU
  TYPE(t_region), POINTER :: pRegion
#endif

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RADI_FlimZeroDummyCells',&
  'RADI_FlimZeroDummyCells.F90' )

! zero out residuals in dummy cells -------------------------------------------

#ifdef RFLO
  iLev = region%currLevel
  IF (ASSOCIATED( region%levels(iLev)%radi%rhs )) THEN
    rhs => region%levels(iLev)%radi%rhs
    CALL RFLO_ZeroDummyCells( region,rhs )
  ENDIF
#endif

#ifdef RFLU
  pRegion => region

  IF (ASSOCIATED(pRegion%radi%rhs)) THEN
    CALL RFLU_ZeroVirtualCellVars(pRegion,pRegion%radi%rhs)
  ENDIF
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_FlimZeroDummyCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_FlimZeroDummyCells.F90,v $
! Revision 1.4  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/05/16 20:46:00  haselbac
! Changed RFLU_ZeroDummyCells to RFLU_ZeroVirtualCellVars
!
! Revision 1.1  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
!
!******************************************************************************







