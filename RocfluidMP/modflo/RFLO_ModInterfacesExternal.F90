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
! $Id: RFLO_ModInterfacesExternal.F90,v 1.9 2008/12/06 08:44:16 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************
  
MODULE RFLO_ModInterfacesExternal

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! Interfaces to external code
! =============================================================================

  SUBROUTINE RFLO_GetBoundaryValues( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetBoundaryValues

  SUBROUTINE RFLO_GetDeformation( region,boundMoved,dNode )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    LOGICAL :: boundMoved(6)
    REAL(RFREAL), POINTER :: dNode(:,:)
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetDeformation

  SUBROUTINE RFLO_SendBoundaryValues( region,initialize )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    LOGICAL :: initialize
  END SUBROUTINE RFLO_SendBoundaryValues

  SUBROUTINE RFLO_SendBoundaryValuesAlpha( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    LOGICAL :: initialize
  END SUBROUTINE RFLO_SendBoundaryValuesAlpha

#ifdef GENX
  SUBROUTINE Fluid_finalize( globalGenx )
    USE ModRocstar, ONLY : t_globalGenx
    TYPE(t_globalGenx), POINTER :: globalGenx
  END SUBROUTINE Fluid_finalize

  SUBROUTINE Fluid_preHdfOutput( globalGenx )
    USE ModRocstar, ONLY : t_globalGenx
    TYPE(t_globalGenx), POINTER :: globalGenx
  END SUBROUTINE Fluid_preHdfOutput

  SUBROUTINE Fluid_postHdfOutput( globalGenx )
    USE ModRocstar, ONLY : t_globalGenx
    TYPE(t_globalGenx), POINTER :: globalGenx
  END SUBROUTINE Fluid_postHdfOutput

  SUBROUTINE RFLO_FlowSolverDummy( globalGenx,timeSystem,dTimeSystem, &
                                   genxHandleBc,genxHandleGm )
    USE ModDataTypes
    USE ModRocstar, ONLY : t_globalGenx
    INTEGER, INTENT(in) :: genxHandleBc, genxHandleGm
    DOUBLE PRECISION, INTENT(in) :: timeSystem, dTimeSystem
    TYPE(t_globalGenx), POINTER  :: globalGenx
  END SUBROUTINE RFLO_FlowSolverDummy

  SUBROUTINE RFLO_InitGenxInterface( regions,handle,solver,inSurf,inVol, &
                                     obtain_attribute )
    USE ModDataStruct, ONLY  : t_region
    CHARACTER(*) :: inSurf, inVol
    INTEGER :: handle, solver, obtain_attribute
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_InitGenxInterface

  SUBROUTINE RFLO_UpdateInbuffGm( globalGenx,dAlpha )
    USE ModDataTypes
    USE ModRocstar, ONLY : t_globalGenx
    DOUBLE PRECISION, INTENT(in) :: dAlpha
    TYPE(t_globalGenx), POINTER  :: globalGenx
  END SUBROUTINE RFLO_UpdateInbuffGm

#ifdef PEUL
  SUBROUTINE PEUL_InitGenxInterface( regions,wins,winv )
    USE ModDataTypes
    USE ModDataStruct, ONLY  : t_region
    TYPE(t_region), POINTER :: regions(:)
    CHARACTER(CHRLEN) :: wins, winv
  END SUBROUTINE PEUL_InitGenxInterface
#endif

#ifdef PLAG
  SUBROUTINE PLAG_InitGenxInterface( regions, wins, inPlag, obtain_attribute )
    USE ModDataTypes
    USE ModDataStruct, ONLY  : t_region
    TYPE(t_region), POINTER :: regions(:)
    CHARACTER(CHRLEN) :: wins, inPlag
    INTEGER :: obtain_attribute
  END SUBROUTINE PLAG_InitGenxInterface

  SUBROUTINE PLAG_SetSizeGenx( region )
    USE ModDataTypes
    USE ModDataStruct, ONLY  : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PLAG_SetSizeGenx
#endif

#ifdef RADI
  SUBROUTINE RADI_InitGenxInterface( regions,wins,winv )
    USE ModDataTypes
    USE ModDataStruct, ONLY  : t_region
    TYPE(t_region), POINTER :: regions(:)
    CHARACTER(CHRLEN) :: wins, winv
  END SUBROUTINE RADI_InitGenxInterface
#endif
#ifdef TURB
  SUBROUTINE TURB_InitGenxInterface( regions,wins,winv )
    USE ModDataTypes
    USE ModDataStruct, ONLY  : t_region
    TYPE(t_region), POINTER :: regions(:)
    CHARACTER(CHRLEN) :: wins, winv
  END SUBROUTINE TURB_InitGenxInterface
#endif

  SUBROUTINE randInitGenxInterface( regions,wins,winv )
    USE ModDataTypes
    USE ModDataStruct, ONLY  : t_region
    TYPE(t_region), POINTER :: regions(:)
    CHARACTER(CHRLEN) :: wins, winv
  END SUBROUTINE randInitGenxInterface
#endif

  END INTERFACE

END MODULE RFLO_ModInterfacesExternal

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModInterfacesExternal.F90,v $
! Revision 1.9  2008/12/06 08:44:16  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2005/12/08 19:57:37  wasistho
! added postHdfOutput
!
! Revision 1.6  2005/12/08 00:20:59  wasistho
! added Fluid_preHdfOutput
!
! Revision 1.5  2004/07/02 22:04:22  fnajjar
! Updated and added PLAG interface calls
!
! Revision 1.4  2004/06/29 23:58:29  wasistho
! migrated to Roccom-3
!
! Revision 1.3  2003/11/21 22:37:16  fnajjar
! Added PLAG and PEUL GenX interfaces
!
! Revision 1.2  2003/08/09 02:02:46  wasistho
! added TURB and RADI_initGenxInterface
!
! Revision 1.1  2002/12/27 22:07:14  jblazek
! Splitted up RFLO_ModInterfaces and ModInterfaces.
!
!******************************************************************************






