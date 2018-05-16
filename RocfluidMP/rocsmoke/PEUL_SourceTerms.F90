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
! Purpose: compute source terms for smoke and add them to the residual
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%peul%rhs = source terms.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_SourceTerms.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_SourceTerms( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters

  USE PEUL_ModInterfaces, ONLY : PEUL_SourceEqEul

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: ipt

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: nPtypes

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_SourceTerms.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_SourceTerms',&
  'PEUL_SourceTerms.F90' )

! begin -----------------------------------------------------------------------

  nPtypes = region%peulInput%nPtypes

  DO ipt = 1,nPtypes

    SELECT CASE (region%peulInput%ptypes(ipt)%methodV)

    CASE (PEUL_METHV_EQEUL)
      CALL PEUL_SourceEqEul(region,ipt)

    END SELECT ! methodV

  END DO ! ipt

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_SourceTerms

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_SourceTerms.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:10:00  haselbac
! Initial revision after changing case
!
! Revision 1.1  2004/05/04 20:27:37  jferry
! Implemented equilibrium Eulerian method
!
!******************************************************************************







