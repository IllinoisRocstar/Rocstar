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
! Purpose: compute numerical dissipation for RaNS equations model.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%turb%diss = turb dissipation.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_RansNumericalDissipation.F90,v 1.3 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RansNumericalDissipation( region ) ! PUBLIC

  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE TURB_ModParameters
#ifdef RFLO
  USE TURB_ModInterfaces, ONLY : TURB_FloRansCentralDissipation
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... local variables
  INTEGER :: discrType, discrOrder

!******************************************************************************

  CALL RegisterFunction( region%global,'TURB_RansNumericalDissipation', &
                         'TURB_RansNumericalDissipation.F90' )

  IF (region%turbInput%modelClass /= MODEL_RANS) GOTO 999

! get parameters --------------------------------------------------------------

  discrType  = region%turbInput%spaceDiscr
  discrOrder = region%turbInput%spaceOrder

! dissipation schemes ---------------------------------------------------------

#ifdef RFLO
  IF ( discrType ==RANS_DISCR_CEN .AND. &
      (discrOrder==RANS_DISCR_ORD1 .OR. discrOrder==RANS_DISCR_ORD2)) THEN
    CALL TURB_FloRansCentralDissipation( region )
  ENDIF
#endif

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( region%global )

END SUBROUTINE TURB_RansNumericalDissipation

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_RansNumericalDissipation.F90,v $
! Revision 1.3  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/03/20 00:29:17  wasistho
! set turb_rflo_ransNumericalDiss to turb_ransNumerical..
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/11 03:26:33  wasistho
! changed rocturb nomenclature
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.2  2003/12/23 21:16:46  wasistho
! split call registerfunction in 2 lines
!
! Revision 1.1  2003/10/27 04:54:38  wasistho
! added RaNS upwind schemes
!
!
!******************************************************************************







