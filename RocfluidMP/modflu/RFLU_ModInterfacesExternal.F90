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
! ******************************************************************************
!
! Purpose: Set explicit interfaces to external subroutines and functions for 
!   Rocflu.
!
! Description: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModInterfacesExternal.F90,v 1.8 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************
  
MODULE RFLU_ModInterfacesExternal

  IMPLICIT NONE

  INTERFACE

#ifdef GENX
  SUBROUTINE RFLU_CheckCouplingInput(regions)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), DIMENSION(:) :: regions   
  END SUBROUTINE RFLU_CheckCouplingInput

  SUBROUTINE Fluid_compute_integrals(globalGenx,integ)
    USE ModRocstar
    DOUBLE PRECISION, DIMENSION(MAN_INTEG_SIZE) :: integ
    TYPE(t_globalGenx), POINTER :: globalGenx
  END SUBROUTINE Fluid_compute_integrals
 
  SUBROUTINE Fluid_finalize(globalGenx)
    USE ModRocstar, ONLY: t_globalGenx
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
 
  SUBROUTINE RFLU_FlowSolverDummy(globalGenx,timeSystem,dTimeSystem, &
                                  genxHandleBc,genxHandleGm)
    USE ModDataTypes
    USE ModRocstar, ONLY: t_globalGenx
    INTEGER, INTENT(in) :: genxHandleBc, genxHandleGm
    DOUBLE PRECISION, INTENT(IN) :: dTimeSystem,timeSystem
    TYPE(t_globalGenx), POINTER :: globalGenx
  END SUBROUTINE RFLU_FlowSolverDummy

  SUBROUTINE RFLU_GetBoundaryValues(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLU_GetBoundaryValues

  SUBROUTINE RFLU_GetDeformation(region)
    USE ModDataTypes
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLU_GetDeformation

  SUBROUTINE RFLU_PutBoundaryValues(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLU_PutBoundaryValues

  SUBROUTINE RFLU_PutBoundaryValuesAlpha(region)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLU_PutBoundaryValuesAlpha

  SUBROUTINE RFLU_UpdateInbuffGm(globalGenx,dAlpha)
    USE ModDataTypes
    USE ModRocstar, ONLY: t_globalGenx
    DOUBLE PRECISION, INTENT(in) :: dAlpha
    TYPE(t_globalGenx), POINTER :: globalGenx
  END SUBROUTINE RFLU_UpdateInbuffGm
#endif

  END INTERFACE

END MODULE RFLU_ModInterfacesExternal

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModInterfacesExternal.F90,v $
! Revision 1.8  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2006/01/07 10:18:16  wasistho
! added Fluid_pre/postHdfOutput
!
! Revision 1.5  2004/10/19 19:28:06  haselbac
! Removed interfaces of routines moved into GENX modules, cosmetics
!
! Revision 1.4  2003/10/03 20:45:08  haselbac
! Removed interface for RFLU_UpdateInbuffWrapper
!
! Revision 1.3  2003/05/01 14:10:18  haselbac
! Added RFLU_CheckCouplingInput
!
! Revision 1.2  2003/04/24 15:38:10  haselbac
! Deleted input argument in RFLU_PutBoundaryValues
!
! Revision 1.1  2003/04/10 14:37:10  haselbac
! Initial revision
!
! ******************************************************************************






