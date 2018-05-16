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
! Output: regions%levels%turb%rhs = residual in dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_RansZeroDummyCells.F90,v 1.5 2009/08/26 12:28:52 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RansZeroDummyCells( region ) ! PUBLIC

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
  REAL(RFREAL), POINTER   :: rhs(:,:), dsterm(:,:)
#endif
  TYPE(t_global), POINTER :: global
#ifdef RFLU
  TYPE(t_region), POINTER :: pRegion
#endif

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'TURB_RansZeroDummyCells',&
  'TURB_RansZeroDummyCells.F90' )

! zero out residuals in dummy cells -------------------------------------------

#ifdef RFLO
  iLev = region%currLevel
  IF (ASSOCIATED( region%levels(iLev)%turb%rhs ) .eqv. .true.) THEN
    rhs => region%levels(iLev)%turb%rhs
    CALL RFLO_ZeroDummyCells( region,rhs )
  ENDIF
  IF (ASSOCIATED( region%levels(iLev)%turb%dsterm ) .eqv. .true.) THEN
    dsterm => region%levels(iLev)%turb%dsterm
    CALL RFLO_ZeroDummyCells( region,dsterm )
  ENDIF
#endif

#ifdef RFLU
  pRegion => region

  IF (ASSOCIATED(pRegion%turb%rhs) .eqv. .true.) THEN
    CALL RFLU_ZeroVirtualCellVars(pRegion,pRegion%turb%rhs)
  ENDIF
  IF (ASSOCIATED(pRegion%turb%dsterm) .eqv. .true.) THEN
    CALL RFLU_ZeroVirtualCellVars(pRegion,pRegion%turb%dsterm)
  ENDIF
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_RansZeroDummyCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_RansZeroDummyCells.F90,v $
! Revision 1.5  2009/08/26 12:28:52  mtcampbe
! Ported to Hera.   Fixed logical expression syntax errors.  Replaced all
! IF (logical_variable)  with IF (logical_variable .eqv. .true.) as
! consistent with the specification.  Also changed: IF( ASSOCIATED(expr) )
! to IF ( ASSOCIATED(expr) .eqv. .true. ).   Intel compilers produce code
! which silently fails for some mal-formed expressions, so these changes
! are a net which should ensure that they are evaluated as intended.
!
! Revision 1.4  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/05/16 20:46:00  haselbac
! Changed RFLU_ZeroDummyCells to RFLU_ZeroVirtualCellVars
!
! Revision 1.1  2004/03/19 02:54:58  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/11 03:26:33  wasistho
! changed rocturb nomenclature
!
! Revision 1.1  2004/03/08 23:35:46  wasistho
! changed turb nomenclature
!
! Revision 1.2  2003/10/21 20:31:36  wasistho
! added dt relaxation in steady flow due to RANS source term
!
! Revision 1.1  2003/10/16 20:20:03  wasistho
! installed RaNS in steady state flow (Exp.Mult.Stg)
!
!
!******************************************************************************







