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
! Purpose: update solution/quantities pertinent to radiation if desired
!
! Description: Depending user input parameters, different quantities are 
!              updated, such as vorticity components.
!
! Input: region = data of current region
!
! Output: specified RADI variables updated
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_SolutionUpdate.F90,v 1.3 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_SolutionUpdate( region ) ! PUBLIC

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE RADI_ModInterfaces, ONLY : RADI_FlimRkUpdate, RADI_FlimEmsUpdate
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... local variables
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RADI_SolutionUpdate',&
  'RADI_SolutionUpdate.F90' )

! update FLD radiation solution -----------------------------------------------

  IF (region%radiInput%radiModel == RADI_MODEL_FLDTRAN) THEN
    IF (global%flowType == FLOW_UNSTEADY .AND. &
        global%solverType == SOLV_EXPLICIT) THEN
      CALL RADI_FlimRkUpdate( region )
    ELSE
      CALL RADI_FlimEmsUpdate( region )
    ENDIF
  ENDIF

! finalize -------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_SolutionUpdate

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_SolutionUpdate.F90,v $
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







