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
! Description: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: SPEC_ModInterfaces.F90,v 1.5 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2001-2003 by the University of Illinois
!
!******************************************************************************

MODULE SPEC_ModInterfaces

  IMPLICIT NONE

  INTERFACE

    SUBROUTINE SPEC_CheckUserInput(regions)
      USE ModDataStruct, ONLY: t_region
      TYPE(t_region), DIMENSION(:), POINTER :: regions
    END SUBROUTINE SPEC_CheckUserInput

    SUBROUTINE SPEC_DerivedInputValues(regions)
      USE ModDataStruct, ONLY: t_region
      TYPE(t_region), DIMENSION(:), POINTER :: regions
    END SUBROUTINE SPEC_DerivedInputValues

    SUBROUTINE SPEC_EqEulCorrPatch(pRegion,pPatch,iSpec)
      USE ModDataStruct, ONLY : t_region
      USE ModBndPatch, ONLY: t_patch
      TYPE(t_region), POINTER :: pRegion
      TYPE(t_patch), POINTER :: pPatch
      INTEGER, INTENT(IN) :: iSpec
    END SUBROUTINE SPEC_EqEulCorrPatch

    SUBROUTINE SPEC_InitInputValues(regions)
      USE ModDataStruct, ONLY: t_region
      TYPE(t_region), DIMENSION(:), POINTER :: regions
    END SUBROUTINE SPEC_InitInputValues

    SUBROUTINE SPEC_InitInputValuesSpecType(region,iSpecType)
      USE ModDataStruct, ONLY: t_region
      INTEGER, INTENT(IN) :: iSpecType
      TYPE(t_region) :: region
    END SUBROUTINE SPEC_InitInputValuesSpecType

    SUBROUTINE SPEC_ReadInputFile(regions)
      USE ModDataStruct, ONLY: t_region
      TYPE(t_region), DIMENSION(:), POINTER :: regions
    END SUBROUTINE SPEC_ReadInputFile

    SUBROUTINE SPEC_ReadSpecSection(regions)
      USE ModDataStruct, ONLY: t_region
      TYPE(t_region), DIMENSION(:), POINTER :: regions
    END SUBROUTINE SPEC_ReadSpecSection

    SUBROUTINE SPEC_ReadSpecTypeSection(regions,iSpecType)
      USE ModDataStruct, ONLY: t_region
      INTEGER :: iSpecType
      TYPE(t_region), DIMENSION(:), POINTER :: regions
    END SUBROUTINE SPEC_ReadSpecTypeSection

  END INTERFACE

END MODULE SPEC_ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_ModInterfaces.F90,v $
! Revision 1.5  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/07/30 22:47:37  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.2  2003/11/25 21:08:33  haselbac
! Added interfaces
!
! Revision 1.1.1.1  2001/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************






