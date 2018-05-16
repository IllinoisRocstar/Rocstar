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
! Purpose: set initial RaNS or LES solution field or parameters at 1st EMSstg. 
!
! Description: none.
!
! Input: region = data of current region,
!        iStage = current RK stage.
!
! Output: region%levels%turb or global% = new solution or parameters at
!         beginning of each time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_EmsInit.F90,v 1.4 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_EmsInit( region, istage ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE TURB_ModParameters
  USE TURB_ModInterfaces, ONLY : TURB_RansEmsInit

  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  
  INTEGER, INTENT(IN) :: istage
  
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'TURB_EmsInit',&
  'TURB_EmsInit.F90' )
  
! which class of model --------------------------------------------------------

  IF (region%turbInput%modelClass == MODEL_RANS) THEN
    CALL TURB_RansEmsInit( region, iStage )
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_EmsInit

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_EmsInit.F90,v $
! Revision 1.4  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/09 22:17:28  wasistho
! removed LES rkInit
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.1  2004/02/26 21:31:17  wasistho
! install TURB_emsInit
!
!
!******************************************************************************







