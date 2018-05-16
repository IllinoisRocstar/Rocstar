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
! Purpose: Obtain viscous fluxes based the fixed Smagorinsky model.
!
! Description: Get eddy viscosity from the fixed Smagorinsky model then
!              add the viscous fluxes to the dissipation residual, mixt%diss.
!              The eddy viscosity is obtained by calling LesCalcEddyVis then
!              the viscous fluxes by calling VisFluxEddy.
!
! Input: region  = data of current region
!        ibn,ien = begin and end node index 
!
! Output: eddy viscosity (mueT) and viscous fluxes.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_LesFluxFixSmag.F90,v 1.6 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_LesFluxFixSmag( region,ibn,ien )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE TURB_ModInterfaces, ONLY : TURB_LesCalcEddyVis
                                 
  USE ModError
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
#ifdef RFLO
  TYPE(t_region)          :: region
#endif
#ifdef RFLU
  TYPE(t_region), POINTER :: region
#endif
  INTEGER                 :: ibn, ien

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_LesFluxFixSmag.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_LesFluxFixSmag',&
  'TURB_LesFluxFixSmag.F90' )

! obtain eddy viscosity in i, j and k faces ----------------------------------
  CALL TURB_LesCalcEddyVis( region,ibn,ien,DIRI )
#ifdef RFLO
  CALL TURB_LesCalcEddyVis( region,ibn,ien,DIRJ )
  CALL TURB_LesCalcEddyVis( region,ibn,ien,DIRK )
#endif

#ifdef RFLU
!  apply boundary conditions for LES variables
!  CALL TURB_FluLesBndConditions( targetFlag )
#endif

! get viscous fluxes
  CALL TURB_VisFluxEddy( region )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_LesFluxFixSmag

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_LesFluxFixSmag.F90,v $
! Revision 1.6  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2004/05/28 02:03:58  wasistho
! update unstructured grid LES
!
! Revision 1.3  2004/03/25 04:40:41  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/19 02:49:56  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.2  2003/10/09 23:07:34  wasistho
! renamed CalcEddyVis to LesCalcEddyVis
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************







