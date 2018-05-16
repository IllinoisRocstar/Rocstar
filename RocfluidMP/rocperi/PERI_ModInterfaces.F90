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
! $Id $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************
  
MODULE PERI_ModInterfaces

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! Interfaces to external code
! =============================================================================

  SUBROUTINE PERI_CheckParamInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PERI_CheckParamInput

  SUBROUTINE PERI_CnlForceTerm( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PERI_CnlForceTerm

  SUBROUTINE PERI_DerivedInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PERI_DerivedInputValues

  SUBROUTINE PERI_InitInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PERI_InitInputValues

  SUBROUTINE PERI_ReadInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PERI_ReadInputFile

  SUBROUTINE PERI_ReadPeriSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PERI_ReadPeriSection

  SUBROUTINE PERI_CoCnlInitSolution( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PERI_CoCnlInitSolution

  SUBROUTINE PERI_CoCprInitSolution( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PERI_CoCprInitSolution

  SUBROUTINE PERI_CoCprSlowTerms( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PERI_CoCprSlowTerms

  SUBROUTINE PERI_CoPgradUpdate( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PERI_CoPgradUpdate

  END INTERFACE

END MODULE PERI_ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_ModInterfaces.F90,v $
! Revision 1.6  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2004/06/08 23:52:31  wasistho
! changed nomenclature
!
!
!******************************************************************************






