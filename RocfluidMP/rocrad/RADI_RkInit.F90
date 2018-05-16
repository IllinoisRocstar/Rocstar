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
! Purpose: set initial FLD radiation solution field or parameters at 1st RKstg. 
!
! Description: none.
!
! Input: region = data of current region,
!        iStage = current RK stage.
!
! Output: region%levels%radi or global% = new solution or parameters at
!         beginning of each time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_RkInit.F90,v 1.3 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_RkInit( region, istage ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE RADI_ModParameters
  USE RADI_ModInterfaces, ONLY : RADI_ExtinctionCoef, RADI_CalcEffTemp, &
                                 RADI_FluxLimiter, RADI_FlimRkInit
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  
  INTEGER, INTENT(IN) :: istage
  
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RADI_RkInit',&
  'RADI_RkInit.F90' )

! global values (currently none) ----------------------------------------------

! extinction (emission, absorbsion, scattering) coefficients ------------------

  CALL RADI_ExtinctionCoef( region )

! radiation effective temperature

  CALL RADI_CalcEffTemp( region )

! which model -----------------------------------------------------------------

  IF ((region%radiInput%radiModel == RADI_MODEL_ROSS) .OR. &
      (region%radiInput%radiModel == RADI_MODEL_FLDSRC) .OR. &
      (region%radiInput%radiModel == RADI_MODEL_FLDTRAN)) THEN
    CALL RADI_FluxLimiter( region )
  ENDIF ! radiModel

  IF (region%radiInput%radiModel == RADI_MODEL_FLDTRAN) THEN
    CALL RADI_FlimRkInit( region, iStage )
  ENDIF ! radiModel

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_RkInit

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_RkInit.F90,v $
! Revision 1.3  2008/12/06 08:44:38  mtcampbe
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







