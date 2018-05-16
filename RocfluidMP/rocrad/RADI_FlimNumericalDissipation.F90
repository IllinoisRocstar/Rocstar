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
! Purpose: compute numerical dissipation for FLD transport equation.
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%radi%diss = radi dissipation.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_FlimNumericalDissipation.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FlimNumericalDissipation( region ) ! PUBLIC

  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE RADI_ModParameters
#ifdef RFLO
  USE RADI_ModInterfaces, ONLY : RADI_FloFlimCentralDissipation
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... local variables
  INTEGER :: discrType, discrOrder

!******************************************************************************

  CALL RegisterFunction( region%global,'RADI_FlimNumericalDissipation', &
                         'RADI_FlimNumericalDissipation.F90' )

  IF (region%radiInput%radiModel /= RADI_MODEL_FLDTRAN) GOTO 999

! get parameters --------------------------------------------------------------

  discrType  = region%radiInput%spaceDiscr
  discrOrder = region%radiInput%spaceOrder

! dissipation schemes ---------------------------------------------------------

#ifdef RFLO
  IF ( discrType ==FLD_DISCR_CEN .AND. &
      (discrOrder==FLD_DISCR_ORD1 .OR. discrOrder==FLD_DISCR_ORD2)) THEN
    CALL RADI_FloFlimCentralDissipation( region )
  ENDIF
#endif

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( region%global )

END SUBROUTINE RADI_FlimNumericalDissipation

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_FlimNumericalDissipation.F90,v $
! Revision 1.3  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
!
!
!******************************************************************************







