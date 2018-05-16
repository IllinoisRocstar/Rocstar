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
! Purpose: set explicit interfaces to subroutines and functions.
!
! Description: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModInterfacesInteract.F90,v 1.8 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfacesInteract

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! Multiphase Interactions
! =============================================================================

  SUBROUTINE INRT_BuildVersionString( versionString )
    CHARACTER(*) :: versionString
  END SUBROUTINE INRT_BuildVersionString

  SUBROUTINE INRT_BurnStatusUpdate( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
  END SUBROUTINE INRT_BurnStatusUpdate

  SUBROUTINE INRT_PrintMaterialInput( global )
    USE ModGlobal,    ONLY : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE INRT_PrintMaterialInput

  SUBROUTINE INRT_PrintUserInput( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(IN) :: region
  END SUBROUTINE INRT_PrintUserInput

  SUBROUTINE INRT_ReadMaterialInput( global )
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE INRT_ReadMaterialInput

  SUBROUTINE INRT_SetMaterial(global,material,name)
    USE ModGlobal,    ONLY : t_global
    USE ModMaterials, ONLY : t_material
    TYPE(t_global),   POINTER    :: global
    TYPE(t_material), POINTER    :: material
    CHARACTER(*),     INTENT(in) :: name
  END SUBROUTINE INRT_SetMaterial

  SUBROUTINE INRT_SetParticleTemp( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE INRT_SetParticleTemp

  SUBROUTINE INRT_SourceTerms( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE INRT_SourceTerms

  SUBROUTINE INRT_TwoDimAverage( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE INRT_TwoDimAverage

  SUBROUTINE INRT_UserInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE INRT_UserInput

  SUBROUTINE INRT_VaporEnergyConversion( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT), TARGET :: region
  END SUBROUTINE INRT_VaporEnergyConversion

  END INTERFACE

END MODULE ModInterfacesInteract

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfacesInteract.F90,v $
! Revision 1.8  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2004/12/01 00:09:25  wasistho
! added BuildVersionString
!
! Revision 1.5  2004/07/27 21:27:14  jferry
! removed rocinteract allocation routines (moved to rocpart)
!
! Revision 1.4  2004/03/02 21:47:28  jferry
! Added After Update interactions
!
! Revision 1.3  2003/03/24 23:25:48  jferry
! moved some routines from libfloflu to rocinteract
!
! Revision 1.2  2003/03/11 15:59:16  jferry
! Added routines to rocinteract
!
! Revision 1.1  2003/03/04 22:12:34  jferry
! Initial import of Rocinteract
!
!******************************************************************************






