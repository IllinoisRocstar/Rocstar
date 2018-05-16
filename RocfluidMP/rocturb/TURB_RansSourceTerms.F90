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
! Purpose: compute source terms, if applicable, turbulence transport equations
!          depending on model selected (only for RaNS and DES) and added to 
!          turb RHS
!
! Description: this routine calls source term routine specific to chosen model
!
! Input: regions = data of all regions
!
! Output: region%levels%turb%rhs updated by source term routine being called
!
! Notes: this routine updates RHS of turbulence transport equations, not
!        of main fluid equations (EULER or NS), hence the name 
!        TURB_RansSourceTerms instead of TURB_SourceTerms.
!
!******************************************************************************
!
! $Id: TURB_RansSourceTerms.F90,v 1.6 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RansSourceTerms( region ) ! PUBLIC

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE TURB_ModInterfaces, ONLY : TURB_RansSASourceTerms
#ifdef RFLU
  USE TURB_ModInterfaces, ONLY : TURB_FluCv2Cons, TURB_FluCv2Prim  
#endif
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... local variables
  TYPE(t_global), POINTER :: global

#ifdef RFLU
  INTEGER :: prevCvState
#endif

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'TURB_RansSourceTerms',&
  'TURB_RansSourceTerms.F90' )

  IF (region%turbInput%modelClass /= MODEL_RANS) GOTO 999

#ifdef RFLU
! Specific Rocflu ------------------------------------------------------------ 
! check the state of cv first and convert to conservative if not yet

  IF (region%mixt%cvState /= CV_MIXT_STATE_CONS) THEN
    prevCvState = region%mixt%cvState
    CALL TURB_FluCv2Cons(region,CV_MIXT_STATE_CONS)
  ENDIF
#endif

! add source terms of any selected turbulence transport to turb%rhs

  IF ((region%mixtInput%turbModel == TURB_MODEL_SA) .OR. &
      (region%mixtInput%turbModel == TURB_MODEL_DESSA) .OR. &
      (region%mixtInput%turbModel == TURB_MODEL_HDESSA)) THEN
    CALL TURB_RansSASourceTerms( region )
  ENDIF

#ifdef RFLU
! convert cv back to the previous state before entering this routine

  IF (region%mixt%cvState /= prevCvState) &
      CALL TURB_FluCv2Prim( region,prevCvState )  
#endif

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_RansSourceTerms

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_RansSourceTerms.F90,v $
! Revision 1.6  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/03/09 06:35:55  wasistho
! incorporated HDESSA
!
! Revision 1.3  2004/03/29 21:10:01  wasistho
! add flu routines
!
! Revision 1.2  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







