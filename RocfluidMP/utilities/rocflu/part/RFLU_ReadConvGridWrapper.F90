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
! Purpose: Read grid file and convert if necessary.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ReadConvGridWrapper.F90,v 1.3 2008/12/06 08:45:03 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ReadConvGridWrapper(pRegion)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region

  USE RFLU_ModCENTAUR
  USE RFLU_ModCOBALT
  USE RFLU_ModGAMBIT
  USE RFLU_ModMESH3D
  USE RFLU_ModTETMESH
  USE RFLU_ModVGRIDns

#ifdef GENX
  USE RFLU_ModRocstarIO, ONLY: RFLU_GENX_DecideReadFile, & 
                            RFLU_GENX_GetGrid
#endif 

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================
   
  TYPE(t_region), POINTER :: pRegion 
   
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ReadConvGridWrapper.F90,v $ $Revision: 1.3 $'

  global => pRegion%global
  
  CALL RegisterFunction(global,'RFLU_ReadConvGridWrapper',&
  'RFLU_ReadConvGridWrapper.F90')

! ******************************************************************************
! Read grid file and convert if necessary 
! ******************************************************************************

  IF ( global%gridSource == GRID_SRC_CENTAUR_ASCII ) THEN    
    CALL RFLU_ReadGridCENTAURASCII(pRegion)
    CALL RFLU_ConvCENTAUR2ROCFLU(pRegion)
  ELSE IF ( global%gridSource == GRID_SRC_CENTAUR_BINARY ) THEN    
    CALL RFLU_ReadGridCENTAURBinary(pRegion)
    CALL RFLU_ConvCENTAUR2ROCFLU(pRegion)
  ELSE IF ( global%gridSource == GRID_SRC_VGRIDNS ) THEN 
    CALL RFLU_ReadGridVGRIDns(pRegion)
    CALL RFLU_ConvVGRIDns2ROCFLU(pRegion)
  ELSE IF ( global%gridSource == GRID_SRC_MESH3D ) THEN
    CALL RFLU_ReadGridMESH3D(pRegion)
    CALL RFLU_ConvMESH3D2ROCFLU(pRegion) 
  ELSE IF ( global%gridSource == GRID_SRC_TETMESH ) THEN 
    CALL RFLU_ReadGridTETMESH(pRegion)
    CALL RFLU_ConvTETMESH2ROCFLU(pRegion)
  ELSE IF ( global%gridSource == GRID_SRC_COBALT ) THEN 
    CALL RFLU_ReadGridCOBALT(pRegion)
    CALL RFLU_ConvCOBALT2ROCFLU(pRegion)
  ELSE IF ( global%gridSource == GRID_SRC_GAMBIT_NEUTRAL ) THEN 
    CALL RFLU_ReadGridGAMBITNeutral(pRegion)
    CALL RFLU_ConvGAMBIT2ROCFLU(pRegion)    
  ELSE 
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)  
  END IF ! global%gridSource

! ******************************************************************************
! End
! ******************************************************************************
 
  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ReadConvGridWrapper


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ReadConvGridWrapper.F90,v $
! Revision 1.3  2008/12/06 08:45:03  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2005/04/15 15:09:18  haselbac
! Initial revision
!
! Revision 1.2  2004/11/03 15:06:30  haselbac
! Added GAMBIT grid conversion option
!
! Revision 1.1  2004/10/19 19:30:36  haselbac
! Initial revision
!
! ******************************************************************************







