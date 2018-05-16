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
! $Id: ModInterfacesEulerian.F90,v 1.10 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfacesEulerian

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! Eulerian particles
! =============================================================================

  SUBROUTINE PEUL_AllocateDataBuffers( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER    :: regions(:)
    INTEGER,        INTENT(IN) :: iReg
  END SUBROUTINE PEUL_AllocateDataBuffers
  
  SUBROUTINE PEUL_AllocateMemory( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE PEUL_AllocateMemory

  SUBROUTINE PEUL_BoundaryConditionsRecv( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER    :: regions(:)
    INTEGER,        INTENT(IN) :: iReg
  END SUBROUTINE PEUL_BoundaryConditionsRecv

  SUBROUTINE PEUL_BoundaryConditionsSend( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER    :: regions(:)
    INTEGER,        INTENT(IN) :: iReg
  END SUBROUTINE PEUL_BoundaryConditionsSend
  
  SUBROUTINE PEUL_BoundaryConditionsSet( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER    :: regions(:)
    INTEGER,        INTENT(IN) :: iReg
  END SUBROUTINE PEUL_BoundaryConditionsSet

  SUBROUTINE PEUL_BuildVersionString( versionString )
    CHARACTER(*) :: versionString
  END SUBROUTINE PEUL_BuildVersionString

  SUBROUTINE PEUL_CentralDissipation( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE PEUL_CentralDissipation

  SUBROUTINE PEUL_ClearSendRequests( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER    :: regions(:)
    INTEGER,        INTENT(IN) :: iReg
  END SUBROUTINE PEUL_ClearSendRequests

  SUBROUTINE PEUL_ConvectiveFluxes( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE PEUL_ConvectiveFluxes

  SUBROUTINE PEUL_EnforcePositivity( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE PEUL_EnforcePositivity

  SUBROUTINE PEUL_InitSolution( iReg,region )
    USE ModDataStruct, ONLY : t_region
    INTEGER,        INTENT(IN)    :: iReg
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE PEUL_InitSolution

  SUBROUTINE PEUL_PrintUserInput( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(IN) :: region
  END SUBROUTINE PEUL_PrintUserInput

  SUBROUTINE PEUL_ReadSolution( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PEUL_ReadSolution

  SUBROUTINE PEUL_ReadSolutionRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER,        INTENT(IN) :: iReg
    TYPE(t_region), POINTER    :: regions(:)
  END SUBROUTINE PEUL_ReadSolutionRegion

  SUBROUTINE PEUL_ResidualSmoothing( region ) 
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE PEUL_ResidualSmoothing

  SUBROUTINE PEUL_ResidualSmoothingCoeffs( region ) 
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE PEUL_ResidualSmoothingCoeffs

  SUBROUTINE PEUL_SourceTerms( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE PEUL_SourceTerms

  SUBROUTINE PEUL_SpectralRadii( region ) 
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), INTENT(INOUT) :: region
  END SUBROUTINE PEUL_SpectralRadii

  SUBROUTINE PEUL_StatMapping( global )
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE PEUL_StatMapping

  SUBROUTINE PEUL_UserInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PEUL_UserInput

  SUBROUTINE PEUL_WriteSolution( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE PEUL_WriteSolution

  END INTERFACE

END MODULE ModInterfacesEulerian

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfacesEulerian.F90,v $
! Revision 1.10  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2004/12/29 23:27:06  wasistho
! prepared statistics for PLAG and PEUL
!
! Revision 1.7  2004/12/01 00:09:14  wasistho
! added BuildVersionString
!
! Revision 1.6  2004/05/03 15:09:41  jferry
! added equilibrium Eulerian capability for smoke
!
! Revision 1.5  2004/03/02 21:44:51  jferry
! Added clipping options
!
! Revision 1.4  2003/05/07 15:09:23  jferry
! Added interface for PEUL_ReadSolutionRegion
!
! Revision 1.3  2003/04/09 14:12:44  fnajjar
! Added Interfaces of MPI communication routines
!
! Revision 1.2  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
! Revision 1.1  2002/12/27 22:07:14  jblazek
! Splitted up RFLO_ModInterfaces and ModInterfaces.
!
!******************************************************************************






