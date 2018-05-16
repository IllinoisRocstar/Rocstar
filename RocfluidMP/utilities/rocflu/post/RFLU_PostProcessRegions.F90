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
! Purpose: Wrapper for postprocessing results without merging the regions.
!
! Description: None.
!
! Input: 
!   levels      Level data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_PostProcessRegions.F90,v 1.26 2008/12/06 08:45:06 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PostProcessRegions(levels)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global 
  USE ModError
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModParameters

  USE RFLU_ModPlottingVars

#ifndef NO_TECPLOT  
  USE RFLU_ModTECPLOT
#endif

#ifdef PLAG
  USE PLAG_ModSurfStats
#endif

  USE ModInterfaces, ONLY: RFLU_ClosePostInfo, &
                           RFLU_DestroyGrid, &
                           RFLU_OpenPostInfo, &
                           RFLU_PostProcessRegionsCommon1, & 
                           RFLU_PostProcessRegionsCommon2, &                           
                           RFLU_ReadPostInfo, & 
                           RFLU_SetPatchPlotFlags
    
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_level), POINTER :: levels(:)
  
! ==============================================================================
! Local variables
! ==============================================================================

  LOGICAL :: postInfoFileExists
  CHARACTER(CHRLEN) :: casename,RCSIdentString
  INTEGER :: errorFlag,iPatch,iReg
  TYPE(t_region), POINTER :: pRegion
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PostProcessRegions.F90,v $ $Revision: 1.26 $'

  global => levels(1)%regions(1)%global

  CALL RegisterFunction(global,'RFLU_PostProcessRegions',&
  'RFLU_PostProcessRegions.F90')

! ******************************************************************************  
! Open files
! ******************************************************************************

  CALL RFLU_OpenPostInfo(global,FILE_STATUS_OLD,postInfoFileExists) 

! ******************************************************************************
! Open files for plotting. NOTE for TECPLOT, initialize interface and need to 
! count number of plotting variables first.
! ******************************************************************************

  SELECT CASE ( global%postOutputFormat ) 
#ifndef NO_TECPLOT
    CASE ( POST_OUTPUT_FORMAT_TECPLOT )  
      CALL RFLU_TEC_Init(global)

      IF ( global%nRegions == 1 ) THEN ! single region
        pRegion => levels(1)%regions(0)
      ELSE 
        pRegion => levels(1)%regions(1)
      END IF ! global%nRegions

      CALL RFLU_CountPlottingVars(pRegion)
      CALL RFLU_CreatePlottingVarMaps(pRegion)
      CALL RFLU_BuildPlottingVarMaps(pRegion)
      CALL RFLU_PrintPlottingVarsInfo(pRegion)

      CALL RFLU_TEC_OpenFileField(pRegion)

      CALL RFLU_DestroyPlottingVarMaps(pRegion)

      IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN
        CALL RFLU_TEC_OpenFilePatch(pRegion)
      END IF ! global%patchCoeffFlag  
    
#ifdef PLAG
      IF ( global%plagUsed .EQV. .TRUE. ) THEN
        CALL RFLU_TEC_OpenFilePnt(pRegion)
        CALL RFLU_TEC_OpenFilePatchStats(pRegion)
      END IF ! global%plagUsed
#endif
#endif
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%postOutputFormat  

! ******************************************************************************  
! Loop over regions
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    IF ( global%nRegions == 1 ) THEN ! single region
      pRegion => levels(1)%regions(0)
    ELSE 
      pRegion => levels(1)%regions(iReg)
    END IF ! global%nRegions

    IF ( postInfoFileExists .EQV. .TRUE. ) THEN 
      CALL RFLU_ReadPostInfo(pRegion,INFOFILE_READMODE_FLAG)    
    END IF ! postInfoFileExists  

! ==============================================================================  
!   Check for activation status and write files
! ==============================================================================  

    IF ( pRegion%postActiveFlag .EQV. .TRUE. ) THEN 
      CALL RFLU_PostProcessRegionsCommon1(pRegion,postInfoFileExists)
   
      SELECT CASE ( global%postOutputFormat )   
#ifndef NO_TECPLOT   
        CASE ( POST_OUTPUT_FORMAT_TECPLOT )  
          IF ( global%postPlotVolFlag .EQV. .TRUE. ) THEN      
            CALL RFLU_TEC_BuildDataFieldVol(pRegion)
            CALL RFLU_TEC_WriteFileFieldVol(pRegion)
            CALL RFLU_TEC_DestroyDataFieldVol(pRegion)                       
          END IF ! global%postPlotVolFlag

          CALL RFLU_SetPatchPlotFlags(pRegion)

          DO iPatch = 1,pRegion%grid%nPatches
            pPatch => pRegion%patches(iPatch)

            IF ( pPatch%plotFlag .EQV. .TRUE. ) THEN 
              CALL RFLU_TEC_BuildDataFieldSurf(pRegion,pPatch)               
              CALL RFLU_TEC_WriteFileFieldSurf(pRegion,pPatch)
              CALL RFLU_TEC_DestroyDataFieldSurf(pRegion,pPatch)
            END IF ! pPatch%plotFlag
          END DO ! iPatch

          IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN 
            DO iPatch = 1,pRegion%grid%nPatches
              pPatch => pRegion%patches(iPatch)

              CALL RFLU_TEC_BuildDataPatch(pRegion,pPatch)  
              CALL RFLU_TEC_WriteFilePatch(pRegion,pPatch)
              CALL RFLU_TEC_DestroyDataPatch(pRegion,pPatch) 
            END DO ! iPatch 
          END IF ! global%patchCoeffFlag                                    

#ifdef PLAG
          IF ( global%plagUsed .EQV. .TRUE. ) THEN
            CALL RFLU_TEC_WriteFilePnt(pRegion)

            DO iPatch = 1,pRegion%grid%nPatches
              pPatch => pRegion%patches(iPatch)

              IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
                CALL RFLU_TEC_BuildDataPatchStats(pRegion,pPatch)  
                CALL RFLU_TEC_WriteFilePatch(pRegion,pPatch)
                CALL RFLU_TEC_DestroyDataPatch(pRegion,pPatch) 
              END IF ! pPatch%plotStatsFlag
            END DO ! iPatch 
          END IF ! global%plagUsed
#endif
#endif                     
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! global%postOutputFormat

      CALL RFLU_PostProcessRegionsCommon2(pRegion)
    END IF ! pRegion%postActiveFlag                       
  END DO ! iReg
  
! ******************************************************************************
! Close files
! ******************************************************************************  
 
  SELECT CASE ( global%postOutputFormat )   
#ifndef NO_TECPLOT   
    CASE ( POST_OUTPUT_FORMAT_TECPLOT )  
      CALL RFLU_TEC_CloseFileField(global)

      IF ( global%patchCoeffFlag .EQV. .TRUE. ) THEN
        CALL RFLU_TEC_CloseFilePatch(global)
      END IF ! global%patchCoeffFlag 

#ifdef PLAG
      IF ( global%plagUsed .EQV. .TRUE. ) THEN
        CALL RFLU_TEC_CloseFilePnt(global)  
        CALL RFLU_TEC_CloseFilePatchStats(global)
      END IF ! global%plagUsed
#endif
#endif      
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%postOutputFormat
  
  IF ( postInfoFileExists .EQV. .TRUE. ) THEN 
    CALL RFLU_ClosePostInfo(global) 
  END IF ! postInfoFileExists  
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PostProcessRegions

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PostProcessRegions.F90,v $
! Revision 1.26  2008/12/06 08:45:06  mtcampbe
! Updated license.
!
! Revision 1.25  2008/11/19 22:18:17  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.24  2007/03/19 21:44:41  haselbac
! Adapted to changes related to plotting variables
!
! Revision 1.23  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.22  2005/12/01 23:38:55  haselbac
! Bug fix: Added IFs on plagUsed
!
! Revision 1.21  2005/10/05 20:53:00  haselbac
! Fixed missing module bugs
!
! Revision 1.20  2005/10/05 20:22:40  haselbac
! Partial rewrite as part of adding ENSIGHT
!
! Revision 1.19  2005/10/05 16:20:11  haselbac
! Bug fix: Added missing USE RFLU_ModStencilCells
!
! Revision 1.18  2005/10/05 14:32:18  haselbac
! Adapted to changes in stencil modules, added use of vertex list module
!
! Revision 1.17  2005/09/23 19:02:51  haselbac
! Added writing of patch stats files
!
! Revision 1.16  2005/08/18 18:48:16  haselbac
! Added postVortCoreFlag to IF statements
!
! Revision 1.15  2005/08/10 00:36:13  haselbac
! Adapted to changes in RFLU_ModPlottingVars
!
! Revision 1.14  2005/08/09 01:10:49  haselbac
! Added IFs for patch coeffs, writing serial files, plotting patches
!
! Revision 1.13  2005/07/25 12:23:30  haselbac
! Added vorticity plotting variables
!
! Revision 1.12  2005/07/01 21:29:01  haselbac
! Bug fix: Added counting of plotting vars before creating them
!
! Revision 1.11  2005/05/31 16:44:58  haselbac
! Enclosed point TEC routines within ifdef PLAG
!
! Revision 1.10  2005/05/18 22:26:41  fnajjar
! ACH: Adapted point files to multiple regions
!
! Revision 1.9  2005/05/09 20:35:22  haselbac
! Bug fix: Added calls for PLAG patch data back in
!
! Revision 1.8  2005/05/03 20:39:34  haselbac
! Only compute plotting vars if postDiscFlag is TRUE
!
! Revision 1.7  2005/05/01 14:22:50  haselbac
! Added postprocessing of plotting vars
!
! Revision 1.6  2005/04/29 12:53:15  haselbac
! Added calls to compute errors at probe locations
!
! Revision 1.5  2005/04/15 15:08:35  haselbac
! Changed call to RFLU_SetVars
!
! Revision 1.4  2005/01/30 22:04:55  haselbac
! Bug fix for no interpolation of data
!
! Revision 1.3  2005/01/20 14:53:27  haselbac
! Cosmetics only
!
! Revision 1.2  2005/01/03 16:05:49  haselbac
! Adapted to changes in RFLU_ModStencils
!
! Revision 1.1  2004/12/29 20:58:54  haselbac
! Initial revision
!
! ******************************************************************************







