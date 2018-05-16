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
! $Id: ModInterfacesRadiation.F90,v 1.8 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfacesRadiation

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! radiation
! =============================================================================

  SUBROUTINE RADI_AllocateMemory( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_AllocateMemory

  SUBROUTINE RADI_BuildVersionString( versionString )
    CHARACTER(*) :: versionString
  END SUBROUTINE RADI_BuildVersionString

  SUBROUTINE RADI_EmsInit( region,iStage )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER :: iStage
  END SUBROUTINE RADI_EmsInit

  SUBROUTINE RADI_FlimConvectiveFluxes( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_FlimConvectiveFluxes

  SUBROUTINE RADI_FlimNumericalDissipation( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_FlimNumericalDissipation

  SUBROUTINE RADI_FlimSourceTerms( region ) 
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_FlimSourceTerms

  SUBROUTINE RADI_FlimZeroDummyCells( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_FlimZeroDummyCells

  SUBROUTINE RADI_InitSolution( iReg,region )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iReg
    TYPE(t_region) :: region
  END SUBROUTINE RADI_InitSolution

  SUBROUTINE RADI_PrintUserInput( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_PrintUserInput

  SUBROUTINE RADI_ReadSolution( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RADI_ReadSolution

  SUBROUTINE RADI_RkInit( region, iStage )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER :: iStage
  END SUBROUTINE RADI_RkInit

  SUBROUTINE RADI_SolutionUpdate( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_SolutionUpdate

  SUBROUTINE RADI_SourceTerms( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_SourceTerms

  SUBROUTINE RADI_UserInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RADI_UserInput

#ifdef RFLO
! =============================================================================
! Rocflo-specific routines
! =============================================================================

  SUBROUTINE RADI_RFLO_FlimAllocDataBuffers( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE RADI_RFLO_FlimAllocDataBuffers

  SUBROUTINE RADI_RFLO_FlimBndConditionsRecv( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE RADI_RFLO_FlimBndConditionsRecv

  SUBROUTINE RADI_RFLO_FlimBndConditionsSend( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE RADI_RFLO_FlimBndConditionsSend

  SUBROUTINE RADI_RFLO_FlimBndConditionsSet( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE RADI_RFLO_FlimBndConditionsSet

  SUBROUTINE RADI_RFLO_FlimClearSendRequests( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE RADI_RFLO_FlimClearSendRequests

  SUBROUTINE RADI_RFLO_FlimResSmoothing( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_RFLO_FlimResSmoothing

  SUBROUTINE RADI_RFLO_FlimResSmoothingCoeff( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_RFLO_FlimResSmoothingCoeff

  SUBROUTINE RADI_RFLO_FlimSpectralRadii( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RADI_RFLO_FlimSpectralRadii

  SUBROUTINE RADI_RFLO_ReadSolution( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RADI_RFLO_ReadSolution

  SUBROUTINE RADI_RFLO_WriteSolution( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RADI_RFLO_WriteSolution
#endif

  END INTERFACE

END MODULE ModInterfacesRadiation

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfacesRadiation.F90,v $
! Revision 1.8  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2004/12/01 00:08:46  wasistho
! added BuildVersionString
!
! Revision 1.5  2004/09/30 17:07:07  wasistho
! prepared for full FLD radiation model
!
! Revision 1.4  2004/09/23 03:49:05  wasistho
! changed RADI_WriteSol.. to RADI_RFLO_WriteSol..
!
! Revision 1.3  2003/07/22 02:54:49  wasistho
! alphabetical ordering
!
! Revision 1.2  2003/07/17 01:01:01  wasistho
! initial activation rocrad
!
! Revision 1.1  2002/12/27 22:07:14  jblazek
! Splitted up RFLO_ModInterfaces and ModInterfaces.
!
!******************************************************************************






