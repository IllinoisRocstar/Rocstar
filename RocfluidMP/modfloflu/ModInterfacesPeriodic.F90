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
! $Id: ModInterfacesPeriodic.F90,v 1.7 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfacesPeriodic

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! periodic flows
! =============================================================================

  SUBROUTINE PERI_AllocateMemory( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PERI_AllocateMemory

  SUBROUTINE PERI_BuildVersionString( versionString )
    CHARACTER(*) :: versionString
  END SUBROUTINE PERI_BuildVersionString

  SUBROUTINE PERI_SourceTerms( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PERI_SourceTerms

  SUBROUTINE PERI_InitSolution( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER         :: iReg
  END SUBROUTINE PERI_InitSolution

  SUBROUTINE PERI_PrintUserInput( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PERI_PrintUserInput

  SUBROUTINE PERI_SolutionUpdate( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE PERI_SolutionUpdate

  SUBROUTINE PERI_UserInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PERI_UserInput

#ifdef RFLU
! =============================================================================
! Rocflu-specific routines
! =============================================================================

  SUBROUTINE PERI_RFLU_DeallocateMemory( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: region
  END SUBROUTINE PERI_RFLU_DeallocateMemory
#endif

  END INTERFACE

END MODULE ModInterfacesPeriodic

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfacesPeriodic.F90,v $
! Revision 1.7  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2004/12/01 00:08:52  wasistho
! added BuildVersionString
!
! Revision 1.4  2004/06/18 15:47:14  wasistho
! correct the location of endif RFLU
!
! Revision 1.3  2004/06/17 23:07:04  wasistho
! added PERI_RFLU_DeallocateMemory
!
! Revision 1.2  2003/04/05 02:03:15  wasistho
! regions to region in PERI_solutionUpdate
!
! Revision 1.1  2003/03/29 03:29:13  wasistho
! install ROCPERI
!
!
!******************************************************************************






