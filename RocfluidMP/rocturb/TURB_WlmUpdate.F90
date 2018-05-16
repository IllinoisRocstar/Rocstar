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
! Purpose: Call routines that update wall stresses and wall heat flux, 
!          depending on wall layer model selected. 
!
! Description: If log-layer model is selected then call to TURB_UpdateLogLay
!              followed by TURB_ReyAnalogy for the heat transfer, while if
!              boundary-layer model is chosen also call to TURB_UpdateBndLay
!              in between to add the pressure gradient and usteady terms.
!
! Input: region  = data of current region.
!        patch   = current patch.
!
! Output: wall stresses and wall heat flux in body fitted coordinate updated.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_WlmUpdate.F90,v 1.5 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_WlmUpdate( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE TURB_ModInterfaces, ONLY : TURB_WlmUpdateBndlay, TURB_WlmReyAnalogy
#ifdef RFLO
  USE TURB_ModInterfaces, ONLY : TURB_FloWlmUpdateLoglay
#endif
  USE ModError
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_WlmUpdate.F90,v $ $Revision: 1.5 $'

  global => region%global
  CALL RegisterFunction( global,'TURB_WlmUpdate',&
  'TURB_WlmUpdate.F90' )

#ifdef RFLO
  CALL TURB_FloWlmUpdateLoglay( region,patch )
#endif
#ifdef RFLU
!  CALL TURB_FluWlmUpdateLoglay( region,patch )
#endif

  IF (patch%valBola%switches(WLM_INPUT_MODEL) == WLM_MODEL_BNDLAY) THEN
    CALL TURB_WlmUpdateBndlay( region,patch )
  ENDIF

  CALL TURB_WlmReyAnalogy( region,patch )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_WlmUpdate

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_WlmUpdate.F90,v $
! Revision 1.5  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/03/24 03:37:03  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/05 04:37:01  wasistho
! changed nomenclature
!
! Revision 1.3  2004/03/02 03:51:18  wasistho
! forgot colon after Id and Log
!
!
!******************************************************************************







