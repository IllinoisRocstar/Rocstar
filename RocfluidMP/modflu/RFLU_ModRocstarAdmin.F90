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
! Purpose: Collection of routines related to GENX interaction.
!
! Description: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModGENXAdmin.F90,v 1.26 2009/05/12 20:20:56 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModRocstarAdmin

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid 
  USE ModBorder, ONLY: t_border
  USE ModMPI 

  USE RFLU_ModRocstarUtils
  USE RFLU_ModRocstarTools, ONLY: RFLU_GENX_InitBFLAG

  IMPLICIT NONE

  INCLUDE 'comf90.h'

  PRIVATE
  PUBLIC :: RFLU_GENX_BuildGridSurf, &
            RFLU_GENX_BuildPConn, & 
            RFLU_GENX_CloseRocinCtrlFiles, &
            RFLU_GENX_CreateAttrDisp, & 
            RFLU_GENX_CreateAttrFlow, & 
            RFLU_GENX_CreateAttrGridSurf, &
            RFLU_GENX_CreateAttrGSpeeds, &
            RFLU_GENX_CreateAttrInterf, & 
            RFLU_GENX_CreateAttrStats, & 
            RFLU_GENX_CreateAttrTurb, & 
            RFLU_GENX_CreateAttrWrapper, & 
            RFLU_GENX_CreateDataInterf, &
            RFLU_GENX_CreateGridSurf, & 
            RFLU_GENX_CreatePConn, & 
            RFLU_GENX_CreateWindows, &
            RFLU_GENX_CreateWindowsDone, & 
            RFLU_GENX_DestroyDataInterf, &            
            RFLU_GENX_DestroyGridSurf, &
            RFLU_GENX_DestroyPConn, & 
            RFLU_GENX_HardCodeWindowName, &
            RFLU_GENX_InitRocman, & 
            RFLU_GENX_ReadCtrlFile, &
            RFLU_GENX_RegisterDataFlow, & 
            RFLU_GENX_RegisterDataGSpeeds, & 
            RFLU_GENX_RegisterDataInterf, & 
            RFLU_GENX_RegisterDataStats, & 
            RFLU_GENX_RegisterDataTurb, & 
            RFLU_GENX_RegisterDataWrapper, & 
            RFLU_GENX_RegisterDisp, & 
            RFLU_GENX_RegisterGrid, &
            RFLU_GENX_RegisterGridSurf, & 
            RFLU_GENX_RegisterGridVol, & 
            RFLU_GENX_SetConnSize, & 
            RFLU_GENX_StoreCommunicator, &
            RFLU_GENX_StoreNamesHandles

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================  
! Private
! ==============================================================================  

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModGENXAdmin.F90,v $ $Revision: 1.26 $'       
  
! ==============================================================================  
! Public
! ==============================================================================  
  

! ******************************************************************************
! Contained routines
! ******************************************************************************

  CONTAINS 
  
  
  


! ******************************************************************************
!
! Purpose: Build surface grid.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Roccom requires locally-numbered boundary-face lists.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_BuildGridSurf(pRegion)
  
    USE RFLU_ModBoundLists, ONLY: RFLU_BuildBFaceLocLists
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    INTEGER :: errorFlag,iPatch,ivg,ivl
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global
    
    pGrid => pRegion%grid
    
    CALL RegisterFunction(global,'RFLU_GENX_BuildGridSurf',&
  'RFLU_ModRocstarAdmin.F90')
    
! ******************************************************************************
!   Create locally-numbered boundary face lists
! ******************************************************************************

    CALL RFLU_BuildBFaceLocLists(pRegion)
  
! ******************************************************************************
!   Copy coordinates. NOTE no longer setting coordinates of virtual vertices to 
!   CRAZY_VALUE_INT because this makes testing of service modules difficult. 
! ******************************************************************************  

    DO iPatch=1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
                
      DO ivl = 1,pPatch%nBVertTot
        ivg = pPatch%bv(ivl)
      
        pPatch%xyz(XCOORD,ivl) = pGrid%xyz(XCOORD,ivg)
        pPatch%xyz(YCOORD,ivl) = pGrid%xyz(YCOORD,ivg)
        pPatch%xyz(ZCOORD,ivl) = pGrid%xyz(ZCOORD,ivg)                
      END DO ! ivl
    END DO ! iPatch
  
! ******************************************************************************
!   End
! ******************************************************************************
  
    CALL DeregisterFunction(global)  
  
  END SUBROUTINE RFLU_GENX_BuildGridSurf

  
  
  



! ******************************************************************************
!
! Purpose: Build pconn dataitem. 
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_BuildPConn(pRegion)
    
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    INTEGER :: errorFlag,iBorder,icl,ipc,ivl,paneId
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_border), POINTER :: pBorder
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global
    
    pGrid => pRegion%grid

    CALL RegisterFunction(global,'RFLU_GENX_BuildPConn',&
  'RFLU_ModRocstarAdmin.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building pconn...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Fill pconn array
! ******************************************************************************  

    ipc = 0 

! ==============================================================================
!   Block 1: Shared vertices
! ==============================================================================

    ipc = ipc + 1
    pGrid%pconn(ipc) = pGrid%nBorders

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ------------------------------------------------------------------------------
!     Pane id
! ------------------------------------------------------------------------------

      CALL RFLU_GENX_BuildPaneId(pBorder%iRegionGlobal,0,paneId)

      ipc = ipc + 1
      pGrid%pconn(ipc) = paneId

! ------------------------------------------------------------------------------
!     Dimension
! ------------------------------------------------------------------------------

      ipc = ipc + 1
      pGrid%pconn(ipc) = pBorder%nVertShared

! ------------------------------------------------------------------------------
!     Data
! ------------------------------------------------------------------------------

      DO ivl = 1,pBorder%nVertShared
        ipc = ipc + 1
        pGrid%pconn(ipc) = pBorder%ivgShared(ivl)
      END DO ! ivl
    END DO ! iBorder      

! ==============================================================================
!   Block 2: Vertices to send
! ==============================================================================

    ipc = ipc + 1
    pGrid%pconn(ipc) = pGrid%nBorders

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ------------------------------------------------------------------------------
!     Pane id
! ------------------------------------------------------------------------------

      CALL RFLU_GENX_BuildPaneId(pBorder%iRegionGlobal,0,paneId)

      ipc = ipc + 1
      pGrid%pconn(ipc) = paneId

! ------------------------------------------------------------------------------
!     Dimension
! ------------------------------------------------------------------------------

      ipc = ipc + 1
      pGrid%pconn(ipc) = pBorder%nVertSend

! ------------------------------------------------------------------------------
!     Data
! ------------------------------------------------------------------------------

      DO ivl = 1,pBorder%nVertSend
        ipc = ipc + 1
        pGrid%pconn(ipc) = pBorder%ivgSend(ivl)
      END DO ! ivl
    END DO ! iBorder      

! ==============================================================================
!   Block 3: Vertices to recv
! ==============================================================================

    ipc = ipc + 1
    pGrid%pconn(ipc) = pGrid%nBorders

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ------------------------------------------------------------------------------
!     Pane id
! ------------------------------------------------------------------------------

      CALL RFLU_GENX_BuildPaneId(pBorder%iRegionGlobal,0,paneId)

      ipc = ipc + 1
      pGrid%pconn(ipc) = paneId

! ------------------------------------------------------------------------------
!     Dimension
! ------------------------------------------------------------------------------

      ipc = ipc + 1
      pGrid%pconn(ipc) = pBorder%nVertRecv

! ------------------------------------------------------------------------------
!     Data
! ------------------------------------------------------------------------------

      DO ivl = 1,pBorder%nVertRecv
        ipc = ipc + 1
        pGrid%pconn(ipc) = pBorder%ivgRecv(ivl)
      END DO ! ivl
    END DO ! iBorder      

! ==============================================================================
!   Block 4: Cells to send
! ==============================================================================

    ipc = ipc + 1
    pGrid%pconn(ipc) = pGrid%nBorders

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ------------------------------------------------------------------------------
!     Pane id
! ------------------------------------------------------------------------------


      CALL RFLU_GENX_BuildPaneId(pBorder%iRegionGlobal,0,paneId)

      ipc = ipc + 1
      pGrid%pconn(ipc) = paneId

! ------------------------------------------------------------------------------
!     Dimension
! ------------------------------------------------------------------------------

      ipc = ipc + 1
      pGrid%pconn(ipc) = pBorder%nCellsSend

! ------------------------------------------------------------------------------
!     Data
! ------------------------------------------------------------------------------

      DO icl = 1,pBorder%nCellsSend
        ipc = ipc + 1
        pGrid%pconn(ipc) = pBorder%icgSend(icl)
      END DO ! ivl
    END DO ! iBorder      

! ==============================================================================
!   Block 5: Cells to recv
! ==============================================================================

    ipc = ipc + 1
    pGrid%pconn(ipc) = pGrid%nBorders

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ------------------------------------------------------------------------------
!     Pane id
! ------------------------------------------------------------------------------

      CALL RFLU_GENX_BuildPaneId(pBorder%iRegionGlobal,0,paneId)

      ipc = ipc + 1
      pGrid%pconn(ipc) = paneId

! ------------------------------------------------------------------------------
!     Dimension
! ------------------------------------------------------------------------------

      ipc = ipc + 1
      pGrid%pconn(ipc) = pBorder%nCellsRecv

! ------------------------------------------------------------------------------
!     Data
! ------------------------------------------------------------------------------

      DO icl = 1,pBorder%nCellsRecv
        ipc = ipc + 1
        pGrid%pconn(ipc) = pBorder%icgRecv(icl)
      END DO ! ivl
    END DO ! iBorder   

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building pconn done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_GENX_BuildPConn



  
  

! ******************************************************************************
!
! Purpose: Close control files for Rocin.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_CloseRocinCtrlFiles(global)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_global), POINTER :: global

! ==============================================================================  
!   Locals 
! ==============================================================================  

    CHARACTER(CHRLEN) :: iFileName
    CHARACTER(GENX_TIME_STRING_LEN) :: timeString
    INTEGER :: errorFlag    
        
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_GENX_CloseRocinCtrlFiles', & 
                          'RFLU_ModRocstarAdmin.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing Rocin control files...'
    END IF ! global%verbLevel
     
! ******************************************************************************
!   Close files
! ******************************************************************************

! ==============================================================================  
!   Volume file 
! ==============================================================================  
 
    CALL RFLU_GENX_BuildTimeString(global%currentTime,timeString)
    iFileName = './Rocflu/Rocin/fluid_in_'//TRIM(timeString)//'.txt'

    CLOSE(IF_CTRL_VOL,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName) 
    END IF ! global%error

! ==============================================================================  
!   Surface file 
! ==============================================================================  

    CALL RFLU_GENX_BuildTimeString(global%currentTime,timeString)
    iFileName = './Rocflu/Rocin/ifluid_in_'//TRIM(timeString)//'.txt'

    CLOSE(IF_CTRL_SURF,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName) 
    END IF ! global%error
 
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Closing Rocin control files done.'
    END IF ! global%verbLevel
      
    CALL DeregisterFunction(global)  
  
  END SUBROUTINE RFLU_GENX_CloseRocinCtrlFiles
  
  
  
  
  

! ******************************************************************************
!
! Purpose: Create new dataitems for displacements.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_CreateAttrDisp(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winName
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
! ******************************************************************************
!   Create new dataitems
! ******************************************************************************
    
    winName = global%volWinName
    
    CALL COM_new_dataitem(TRIM(winName)//'.disp','n',COM_DOUBLE_PRECISION,3,'')

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_CreateAttrDisp
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Create new dataitems for mixture solution.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_CreateAttrFlow(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winName
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global
    
! ******************************************************************************
!   Create new dataitems
! ******************************************************************************

! ==============================================================================  
!   Volume
! ==============================================================================  

    winName = global%volWinName
    
    CALL COM_new_dataitem(TRIM(winName)//'.rhof' ,'e',COM_DOUBLE_PRECISION, &
                           1,'kg/(m^3)')
    CALL COM_new_dataitem(TRIM(winName)//'.rhovf','e',COM_DOUBLE_PRECISION, &
                           3,'kg/(m^2 s)')
    CALL COM_new_dataitem(TRIM(winName)//'.rhoEf','e',COM_DOUBLE_PRECISION, &
                           1,'(J/(m^3))')  
    CALL COM_new_dataitem(TRIM(winName)//'.pf','e',COM_DOUBLE_PRECISION,1, & 
                           'N/(m^2)')  
    CALL COM_new_dataitem(TRIM(winName)//'.Tf','e',COM_DOUBLE_PRECISION,1, &
                           'K')
    CALL COM_new_dataitem(TRIM(winName)//'.af','e',COM_DOUBLE_PRECISION,1, &
                           'm/s')
  
! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_CreateAttrFlow  
  
  
  
  
  


! ******************************************************************************
!
! Purpose: Create new dataitems for surface grid.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Must not access dimensions of grid, because not necessarily known.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_CreateAttrGridSurf(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winName
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
! ******************************************************************************
!   Create new dataitems
! ******************************************************************************

    winName = global%surfWinName

    CALL COM_new_dataitem(TRIM(winName)//'.bcflag'     ,'p',COM_INTEGER,1,'')
    CALL COM_new_dataitem(TRIM(winName)//'.patchNo'    ,'p',COM_INTEGER,1,'')    
    CALL COM_new_dataitem(TRIM(winName)//'.cnstr_type' ,'p',COM_INTEGER,1,'')
  
    CALL COM_new_dataitem(TRIM(winName)//'.t3g:real'   ,'p',COM_INTEGER,3,'')
    CALL COM_new_dataitem(TRIM(winName)//'.t3g:virtual','p',COM_INTEGER,3,'')    

    CALL COM_new_dataitem(TRIM(winName)//'.q4g:real'   ,'p',COM_INTEGER,4,'')
    CALL COM_new_dataitem(TRIM(winName)//'.q4g:virtual','p',COM_INTEGER,4,'')    

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_CreateAttrGridSurf
  
  

 




! ******************************************************************************
!
! Purpose: Create new dataitems for grid speeds.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_CreateAttrGSpeeds(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winName
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global
    
! ******************************************************************************
!   Create new dataitems
! ******************************************************************************

! ==============================================================================  
!   Volume
! ==============================================================================  

    winName = global%volWinName
    
    CALL COM_new_dataitem(TRIM(winName)//'.gs','p',COM_DOUBLE_PRECISION, & 
                           1,'m/s') 

! ==============================================================================  
!   Surface
! ==============================================================================  

    winName = global%surfWinName
  
    CALL COM_new_dataitem(TRIM(winName)//'.gs','e',COM_DOUBLE_PRECISION, & 
                           1,'m/s')
  
! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_CreateAttrGSpeeds







! ******************************************************************************
!
! Purpose: Create new dataitems for interface data.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_CreateAttrInterf(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winName
    INTEGER :: comd,comi
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global

    winName = global%surfWinName
    
    comd = COM_DOUBLE_PRECISION 
    comi = COM_INTEGER 
    
! ******************************************************************************
!   Create new dataitems
! ******************************************************************************

! ==============================================================================  
!   Incoming data
! ==============================================================================  
    
    CALL COM_new_dataitem(TRIM(winName)//'.du_alp'    ,'n',comd,3,'m'        )
    CALL COM_new_dataitem(TRIM(winName)//'.mdot_alp'  ,'e',comd,1,'kg/(m^2s)')
    CALL COM_new_dataitem(TRIM(winName)//'.rhofvf_alp','e',comd,3,'kg/(m^2s)')
    CALL COM_new_dataitem(TRIM(winName)//'.Tflm_alp'  ,'e',comd,1,'K'        )
    CALL COM_new_dataitem( TRIM(winName)//'.zoomFact' ,'w',comd,1,'none' )
    CALL COM_new_dataitem(TRIM(winName)//'.Tb_alp'    ,'e',comd,1,'K'        )

! ==============================================================================  
!   Outgoing data
! ==============================================================================  

    CALL COM_new_dataitem(TRIM(winName)//'.nf_alp'  ,'e',comd,3,''          )
    CALL COM_new_dataitem(TRIM(winName)//'.pf'      ,'e',comd,1,'Pa'        )
    CALL COM_new_dataitem(TRIM(winName)//'.tf'      ,'e',comd,3,'Pa'        )  
    CALL COM_new_dataitem(TRIM(winName)//'.qc'      ,'e',comd,1,'kgK/(m^2s)')
    CALL COM_new_dataitem(TRIM(winName)//'.qr'      ,'e',comd,1,'kgK/(m^2s)')
    CALL COM_new_dataitem(TRIM(winName)//'.rhof_alp','e',comd,1,'kg/m^3'    )
    CALL COM_new_dataitem(TRIM(winName)//'.Tf'      ,'e',comd,1,'K'         )
!    CALL COM_new_dataitem(TRIM(winName)//'.Tv'      ,'e',comd,1,'K'         )
!    CALL COM_new_dataitem(TRIM(winName)//'.dn'      ,'e',comd,1,'m'         )

! ==============================================================================  
!   Interface data (utility data used in the interface)
! ==============================================================================  

    CALL COM_new_dataitem(TRIM(winName)//'.bflag'   ,'e',comi,1,''          )
      
! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_CreateAttrInterf






  
! ******************************************************************************
!
! Purpose: Create new dataitems for statistics data.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_CreateAttrStats(pRegion)
#ifdef STATS  
    USE STAT_RFLU_ModRocstarAdmin, ONLY: STAT_RFLU_GenxCreateAttr
#endif
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion
    
! ******************************************************************************
!   Create new dataitems
! ******************************************************************************
#ifdef STATS
    CALL STAT_RFLU_GenxCreateAttr(pRegion)
#endif  
! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_CreateAttrStats






  
! ******************************************************************************
!
! Purpose: Create new dataitems for turbulence.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_CreateAttrTurb(pRegion)
#ifdef TURB  
    USE TURB_RFLU_ModRocstarAdmin, ONLY: TURB_RFLU_GenxCreateAttr
#endif
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion
    
! ******************************************************************************
!   Create new dataitems
! ******************************************************************************
#ifdef TURB
    CALL TURB_RFLU_GenxCreateAttr(pRegion)
#endif  
! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_CreateAttrTurb






  
! ******************************************************************************
!
! Purpose: Wrapper for new dataitems routines.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_CreateAttrWrapper(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion
    
! ******************************************************************************
!   Create new dataitems
! ******************************************************************************

! ==============================================================================  
!   Mixture
! ==============================================================================  

    CALL RFLU_GENX_CreateAttrFlow(pRegion)
    CALL RFLU_GENX_CreateAttrGridSurf(pRegion)
    CALL RFLU_GENX_CreateAttrGSpeeds(pRegion)     
    CALL RFLU_GENX_CreateAttrInterf(pRegion)
    CALL RFLU_GENX_CreateAttrDisp(pRegion)

! ==============================================================================  
!   Statistics
! ==============================================================================  

    CALL RFLU_GENX_CreateAttrStats(pRegion)

! ==============================================================================  
!   Turbulence
! ==============================================================================  

    CALL RFLU_GENX_CreateAttrTurb(pRegion)
  
! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_CreateAttrWrapper







! ******************************************************************************
!
! Purpose: Create interface data arrays.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_CreateDataInterf(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    INTEGER :: errorFlag,ifl,iPatch,ivl    
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
      
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

    pGrid => pRegion%grid 

    CALL RegisterFunction(global,'RFLU_GENX_CreateDataInterf',&
  'RFLU_ModRocstarAdmin.F90')

! ******************************************************************************
!   Allocate memory 
! ******************************************************************************  
  
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
    
! ==============================================================================
!     Surface grid deformation
! ==============================================================================

      ALLOCATE(pPatch%dXyz(XCOORD:ZCOORD,pPatch%nBVertTot),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%dXyz')
      END IF ! global%error      

      DO ivl = 1,pPatch%nBVertTot   
        pPatch%dXyz(XCOORD,ivl) = 0.0_RFREAL
        pPatch%dXyz(YCOORD,ivl) = 0.0_RFREAL
        pPatch%dXyz(ZCOORD,ivl) = 0.0_RFREAL            
      END DO ! ivl

! ==============================================================================
!     Input data
! ==============================================================================      
      
! ------------------------------------------------------------------------------
!     All patches      
! ------------------------------------------------------------------------------
      
      ALLOCATE(pPatch%duAlp(XCOORD:ZCOORD,pPatch%nBVertTot),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%duAlp')
      END IF ! global

      DO ivl = 1,pPatch%nBVertTot 
        pPatch%duAlp(XCOORD,ivl) = 0.0_RFREAL
        pPatch%duAlp(YCOORD,ivl) = 0.0_RFREAL
        pPatch%duAlp(ZCOORD,ivl) = 0.0_RFREAL      
      END DO ! ivl

! ------------------------------------------------------------------------------
!     Interacting patches (burning or non-burning)       
! ------------------------------------------------------------------------------

      IF ( pPatch%bcCoupled /= BC_NOT_COUPLED ) THEN 
        ALLOCATE(pPatch%rhofvfAlp(XCOORD:ZCOORD,pPatch%nBFacesTot), & 
                 STAT=errorFlag)
        global%error = errorFlag         
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%rhofvfAlp')
        END IF ! global

        DO ifl = 1,pPatch%nBFacesTot 
          pPatch%rhofvfAlp(XCOORD,ifl) = 0.0_RFREAL
          pPatch%rhofvfAlp(YCOORD,ifl) = 0.0_RFREAL
          pPatch%rhofvfAlp(ZCOORD,ifl) = 0.0_RFREAL        
        END DO ! ifl
        
        ALLOCATE(pPatch%tbAlp(pPatch%nBFacesTot),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%tbAlp')
        END IF ! global

        DO ifl = 1,pPatch%nBFacesTot 
          pPatch%tbAlp(ifl) = 0.0_RFREAL
        END DO ! ifl           
      END IF ! pPatch

! ------------------------------------------------------------------------------
!     Burning patches 
! ------------------------------------------------------------------------------

      IF ( pPatch%bcCoupled == BC_BURNING ) THEN 
        ALLOCATE(pPatch%mdotAlp(pPatch%nBFacesTot),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mdotAlp')
        END IF ! global

        DO ifl = 1,pPatch%nBFacesTot 
          pPatch%mdotAlp(ifl) = 0.0_RFREAL
        END DO ! ifl

        ALLOCATE(pPatch%tflmAlp(pPatch%nBFacesTot),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%tflmAlp')
        END IF ! global

        DO ifl = 1,pPatch%nBFacesTot 
          pPatch%tflmAlp(ifl) = 0.0_RFREAL
        END DO ! ifl      
      END IF ! pPatch  
      
! ==============================================================================
!     Output data
! ==============================================================================      

! ------------------------------------------------------------------------------
!     Interacting patches (burning or non-burning)       
! ------------------------------------------------------------------------------

      IF ( pPatch%bcCoupled /= BC_NOT_COUPLED ) THEN 
        ALLOCATE(pPatch%nfAlp(XCOORD:ZCOORD,pPatch%nBFacesTot),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%nfAlp')
        END IF ! global             

        DO ifl = 1,pPatch%nBFacesTot 
          pPatch%nfAlp(XCOORD,ifl) = 0.0_RFREAL
          pPatch%nfAlp(YCOORD,ifl) = 0.0_RFREAL
          pPatch%nfAlp(ZCOORD,ifl) = 0.0_RFREAL                
        END DO ! ifl

        ALLOCATE(pPatch%rhofAlp(pPatch%nBFacesTot),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%rhofAlp')
        END IF ! global                    

        DO ifl = 1,pPatch%nBFacesTot 
          pPatch%rhofAlp(ifl) = 0.0_RFREAL
        END DO ! ifl

        ALLOCATE(pPatch%pf(pPatch%nBFacesTot),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%pf')
        END IF ! global            

        DO ifl = 1,pPatch%nBFacesTot 
          pPatch%pf(ifl) = 0.0_RFREAL
        END DO ! ifl

        ALLOCATE(pPatch%tracf(XCOORD:ZCOORD,pPatch%nBFacesTot),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%tracf')
        END IF ! global

        DO ifl = 1,pPatch%nBFacesTot 
          pPatch%tracf(XCOORD,ifl) = 0.0_RFREAL
          pPatch%tracf(YCOORD,ifl) = 0.0_RFREAL
          pPatch%tracf(ZCOORD,ifl) = 0.0_RFREAL                             
        END DO ! ifl      
        
        ALLOCATE(pPatch%qc(pPatch%nBFacesTot),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%qc')
        END IF ! global             

        DO ifl = 1,pPatch%nBFacesTot 
          pPatch%qc(ifl) = 0.0_RFREAL
        END DO ! ifl

        ALLOCATE(pPatch%qr(pPatch%nBFacesTot),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%qr')
        END IF ! global             

        DO ifl = 1,pPatch%nBFacesTot 
          pPatch%qr(ifl) = 0.0_RFREAL
        END DO ! ifl        
      END IF ! pPatch%bcCoupled

! ------------------------------------------------------------------------------
!     Burning patches
! ------------------------------------------------------------------------------

      IF ( pPatch%bcCoupled == BC_BURNING ) THEN     
        ALLOCATE(pPatch%tempf(pPatch%nBFacesTot),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%tempf')
        END IF ! global  

        DO ifl = 1,pPatch%nBFacesTot 
          pPatch%tempf(ifl) = 0.0_RFREAL
        END DO ! ifl

! TEMPORARY
!       bFlag will need to be removed from Rocflu eventually. For the 
!       moment, leave in, but always set bFlag to 1. This means that igniting 
!       computations cannot be run.
    
        ALLOCATE(pPatch%bFlag(pPatch%nBFacesTot),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bFlag')
        END IF ! global  

        DO ifl = 1,pPatch%nBFaces 
          pPatch%bFlag(ifl) = 1
        END DO ! ifl   

        DO ifl = pPatch%nBFaces+1,pPatch%nBFacesTot 
          pPatch%bFlag(ifl) = CRAZY_VALUE_INT
        END DO ! ifl        
! END TEMPORARY 
      END IF ! pPatch%bcType
    END DO ! iPatch  
  
! ******************************************************************************
!   End
! ******************************************************************************

    CALL RFLU_GENX_InitBFLAG(pRegion)

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_GENX_CreateDataInterf










! ******************************************************************************
!
! Purpose: Create surface grid.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Roccom requires locally-numbered boundary-face lists.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_CreateGridSurf(pRegion)
  
    USE RFLU_ModBoundLists, ONLY: RFLU_CreateBFaceLocLists
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    INTEGER :: errorFlag,iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global
    
    pGrid => pRegion%grid

    CALL RegisterFunction(global,'RFLU_GENX_CreateGridSurf',&
  'RFLU_ModRocstarAdmin.F90')
    
! ******************************************************************************
!   Create locally-numbered boundary face lists
! ******************************************************************************

    CALL RFLU_CreateBFaceLocLists(pRegion)
  
! ******************************************************************************
!   Allocate memory and initialize 
! ******************************************************************************  

    DO iPatch=1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      ALLOCATE(pPatch%bcFlag(1),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bcFlag')
      END IF ! global            
                
      pPatch%bcFlag(1) = pPatch%bcCoupled
                
      ALLOCATE(pPatch%patchNo(1),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%patchNo')
      END IF ! global                 

      pPatch%patchNo(1) = pPatch%iPatchGlobal          
                
      ALLOCATE(pPatch%cnstrType(1),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%cnstrType')
      END IF ! global                 

      pPatch%cnstrType(1) = RFLU_GENX_SetCnstrType(pPatch%movePatchDir)         
                
      ALLOCATE(pPatch%xyz(XCOORD:ZCOORD,pPatch%nBVertTot),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%xyz')
      END IF ! global%error  
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_GENX_CreateGridSurf







! ******************************************************************************
!
! Purpose: Create pconn dataitem.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_CreatePConn(pRegion)
    
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    INTEGER :: errorFlag,iBorder
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_border), POINTER :: pBorder
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global
    
    pGrid => pRegion%grid

    CALL RegisterFunction(global,'RFLU_GENX_CreatePConn',&
  'RFLU_ModRocstarAdmin.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating pconn...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Determine size of pconn dataitem
! ******************************************************************************  

    pGrid%pconnSizeTot   = 0
    pGrid%pconnSizeGhost = 0 

! ==============================================================================
!   Block 1: Shared vertices
! ==============================================================================

    pGrid%pconnSizeTot = pGrid%pconnSizeTot + 1

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      pGrid%pconnSizeTot = pGrid%pconnSizeTot + pBorder%nVertShared + 2  
    END DO ! iBorder      

! ==============================================================================
!   Block 2: Vertices to send
! ==============================================================================

    pGrid%pconnSizeTot   = pGrid%pconnSizeTot   + 1
    pGrid%pconnSizeGhost = pGrid%pconnSizeGhost + 1

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      pGrid%pconnSizeTot   = pGrid%pconnSizeTot   + pBorder%nVertSend + 2 
      pGrid%pconnSizeGhost = pGrid%pconnSizeGhost + pBorder%nVertSend + 2 
    END DO ! iBorder   

! ==============================================================================
!   Block 3: Vertices to recv
! ==============================================================================

    pGrid%pconnSizeTot   = pGrid%pconnSizeTot   + 1
    pGrid%pconnSizeGhost = pGrid%pconnSizeGhost + 1

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      pGrid%pconnSizeTot   = pGrid%pconnSizeTot   + pBorder%nVertRecv  + 2 
      pGrid%pconnSizeGhost = pGrid%pconnSizeGhost + pBorder%nVertRecv  + 2 
    END DO ! iBorder 

! ==============================================================================
!   Block 4: Cells to send
! ==============================================================================

    pGrid%pconnSizeTot   = pGrid%pconnSizeTot   + 1
    pGrid%pconnSizeGhost = pGrid%pconnSizeGhost + 1

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      pGrid%pconnSizeTot   = pGrid%pconnSizeTot   + pBorder%nCellsSend + 2 
      pGrid%pconnSizeGhost = pGrid%pconnSizeGhost + pBorder%nCellsSend + 2 
    END DO ! iBorder   

! ==============================================================================
!   Block 5: Cells to recv
! ==============================================================================

    pGrid%pconnSizeTot   = pGrid%pconnSizeTot   + 1
    pGrid%pconnSizeGhost = pGrid%pconnSizeGhost + 1

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

      pGrid%pconnSizeTot   = pGrid%pconnSizeTot   + pBorder%nCellsRecv + 2 
      pGrid%pconnSizeGhost = pGrid%pconnSizeGhost + pBorder%nCellsRecv + 2 
    END DO ! iBorder 

! ******************************************************************************
!   Allocate memory 
! ******************************************************************************  

    ALLOCATE(pGrid%pconn(pGrid%pconnSizeTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%pconn')
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating pconn done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_GENX_CreatePConn








! ******************************************************************************
!
! Purpose: Create windows.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   communicator        Communicator
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_CreateWindows(pRegion,communicator)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    INTEGER, INTENT(IN), OPTIONAL :: communicator
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global
    
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating windows...'
    END IF ! global%myProcid
    
! ******************************************************************************
!   Set names and create windows
! ******************************************************************************

    global%surfWinName = TRIM(global%winName)//'_surf'
    global%volWinName = TRIM(global%winName)//'_vol'

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'Surface window name:', & 
                                    TRIM(global%surfWinName)
      WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'Volume window name:', & 
                                    TRIM(global%volWinName)                                    
    END IF ! global%myProcid
    
    IF ( PRESENT(communicator) ) THEN 
      CALL COM_new_window(TRIM(global%surfWinName),communicator)
      CALL COM_new_window(TRIM(global%volWinName),communicator)
    ELSE 
      CALL COM_new_window(TRIM(global%surfWinName))
      CALL COM_new_window(TRIM(global%volWinName))      
    END IF ! PRESENT

! ******************************************************************************
!   End
! ******************************************************************************
  
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating windows done.'
    END IF ! global%myProcid  
  
  END SUBROUTINE RFLU_GENX_CreateWindows  







! ******************************************************************************
!
! Purpose: Inform Roccom that window creation is done.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_CreateWindowsDone(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global
    
! ******************************************************************************
!   Inform Roccom that window creation done
! ******************************************************************************

    global%surfWinName = TRIM(global%winName)//'_surf'
    CALL COM_window_init_done(TRIM(global%surfWinName))

    global%volWinName = TRIM(global%winName)//'_vol'
    CALL COM_window_init_done(TRIM(global%volWinName))

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_CreateWindowsDone








! ******************************************************************************
!
! Purpose: Destroy interface data arrays.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_DestroyDataInterf(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    INTEGER :: errorFlag,iPatch   
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
      
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

    pGrid => pRegion%grid 

    CALL RegisterFunction(global,'RFLU_GENX_DestroyDataInterf',&
  'RFLU_ModRocstarAdmin.F90')

! ******************************************************************************
!   Allocate memory 
! ******************************************************************************  
  
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
    
! ==============================================================================
!     Surface grid deformation
! ==============================================================================

      DEALLOCATE(pPatch%dXyz,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%dXyz')
      END IF ! global%error      

! ==============================================================================
!     Input data
! ==============================================================================      
      
! ------------------------------------------------------------------------------
!     All patches       
! ------------------------------------------------------------------------------
      
      DEALLOCATE(pPatch%duAlp,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%duAlp')
      END IF ! global

! ------------------------------------------------------------------------------
!     Interacting patches (burning or non-burning)       
! ------------------------------------------------------------------------------

      IF ( pPatch%bcCoupled /= BC_NOT_COUPLED ) THEN 
        DEALLOCATE(pPatch%rhofvfAlp,STAT=errorFlag)
        global%error = errorFlag         
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%rhofvfAlp')
        END IF ! global
        
        DEALLOCATE(pPatch%tbAlp,STAT=errorFlag)
        global%error = errorFlag         
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%tbfAlp')
        END IF ! global        
      END IF ! pPatch

! ------------------------------------------------------------------------------
!     Burning patches
! ------------------------------------------------------------------------------

      IF ( pPatch%bcCoupled == BC_BURNING ) THEN 
        DEALLOCATE(pPatch%mdotAlp,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%mdotAlp')
        END IF ! global

        DEALLOCATE(pPatch%tflmAlp,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%tflmAlp')
        END IF ! global
      END IF ! pPatch  
      
! ==============================================================================
!     Output data
! ==============================================================================      

! ------------------------------------------------------------------------------
!     Interacting patches (burning or non-burning)       
! ------------------------------------------------------------------------------

      IF ( pPatch%bcCoupled /= BC_NOT_COUPLED ) THEN 
        DEALLOCATE(pPatch%nfAlp,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%nfAlp')
        END IF ! global             

        DEALLOCATE(pPatch%rhofAlp,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%rhofAlp')
        END IF ! global                    

        DEALLOCATE(pPatch%pf,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%pf')
        END IF ! global            

        DEALLOCATE(pPatch%tracf,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%tracf')
        END IF ! global
        
        DEALLOCATE(pPatch%qc,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%qc')
        END IF ! global             

        DEALLOCATE(pPatch%qr,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%qr')
        END IF ! global             
      END IF ! pPatch%bcCoupled

! ------------------------------------------------------------------------------
!     Burning patches 
! ------------------------------------------------------------------------------

      IF ( pPatch%bcCoupled == BC_BURNING ) THEN     
        DEALLOCATE(pPatch%tempf,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%tempf')
        END IF ! global
        
! TO BE REMOVED - Needs discussion    
        DEALLOCATE(pPatch%bFlag,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bFlag')
        END IF ! global          
! END TO BE REMOVED
      END IF ! pPatch%bcType
    END DO ! iPatch  
  
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_GENX_DestroyDataInterf








! ******************************************************************************
!
! Purpose: Destroy surface grid.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_DestroyGridSurf(pRegion)
  
    USE RFLU_ModBoundLists, ONLY: RFLU_DestroyBFaceLocLists  
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    INTEGER :: errorFlag,iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global
    
    pGrid => pRegion%grid

    CALL RegisterFunction(global,'RFLU_GENX_DestroyGridSurf',&
  'RFLU_ModRocstarAdmin.F90')
    
! ******************************************************************************
!   Create locally-numbered boundary face lists
! ******************************************************************************

    CALL RFLU_DestroyBFaceLocLists(pRegion)
  
! ******************************************************************************
!   Allocate memory 
! ******************************************************************************  

    DO iPatch=1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
     
      DEALLOCATE(pPatch%bcFlag,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bcFlag')
      END IF ! global
      
      DEALLOCATE(pPatch%patchNo,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%patchNo')
      END IF ! global          
             
      DEALLOCATE(pPatch%cnstrType,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%cnstrType')
      END IF ! global        
             
      DEALLOCATE(pPatch%xyz,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%xyz')
      END IF ! global%error  
    END DO ! iPatch
  
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_GENX_DestroyGridSurf










! ******************************************************************************
!
! Purpose: Destroy pconn dataitem.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_DestroyPConn(pRegion)
    
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    INTEGER :: errorFlag
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global
    
    pGrid => pRegion%grid

    CALL RegisterFunction(global,'RFLU_GENX_DestroyPConn',&
  'RFLU_ModRocstarAdmin.F90')     

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying pconn...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Deallocate memory 
! ******************************************************************************  

    DEALLOCATE(pGrid%pconn,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%pconn')
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying pconn done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_GENX_DestroyPConn











! ******************************************************************************
!
! Purpose: Hard-code window name.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!
! Output: None.
!
! Notes: 
!   1. Routine is needed because some utilities are not called through Roccom.
!      This means that rocflu_load_module is not called in which the window name
!      is set. This window name must be set before the windows are created.
!   2. This routine MUST ONLY BE CALLED IF the parent code is not called through
!      Roccom. 
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_HardCodeWindowName(global)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_global), POINTER :: global
        
! ******************************************************************************
!   Set names and register windows
! ******************************************************************************

    global%winName = 'ROCFLU'

! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_HardCodeWindowName








! ******************************************************************************
!
! Purpose: Initialize Rocman.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   handle              Handle
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_InitRocman(pRegion,handle)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    INTEGER, INTENT(IN) :: handle
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    INTEGER :: solverType
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing Rocman...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Initialize Rocman
! ******************************************************************************

    IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN 
      solverType = 1
    ELSE  
      solverType = 0
    END IF ! pRegion       

    CALL COM_call_function(handle,3,TRIM(global%surfWinName), & 
                           TRIM(global%volWinName),solverType)
   
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing Rocman done.'
    END IF ! global%verbLevel 
  
  END SUBROUTINE RFLU_GENX_InitRocman








! ******************************************************************************
!
! Purpose: Read Rocflu control file for GENX runs.
!
! Description: None.
!
! Input:
!   global      global type
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_ReadCtrlFile(global)

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: errorFlag,verbLevelCOM

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_GENX_ReadCtrlFile',&
  'RFLU_ModRocstarAdmin.F90')

! ==============================================================================
!   Open file
! ==============================================================================

    iFileName = './Rocflu/RocfluControl.txt'

    OPEN(IF_CONTROL,FILE=iFileName,FORM='FORMATTED',STATUS='OLD', & 
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error

! ==============================================================================
!   Read file
! ==============================================================================

! ------------------------------------------------------------------------------
!   Case name
! ------------------------------------------------------------------------------

    READ(IF_CONTROL,'(A)',IOSTAT=errorFlag) global%casename
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop( global,ERR_FILE_READ,__LINE__,'File: '//TRIM(iFileName) )
    END IF ! global%error

! ------------------------------------------------------------------------------
!   Input and output directories
! ------------------------------------------------------------------------------

    READ(IF_CONTROL,'(A)',IOSTAT=errorFlag) global%inDir
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      global%inDir = './'
    ELSE
      IF ( global%inDir(LEN_TRIM(global%inDir): & 
                        LEN_TRIM(global%inDir)) /= '/' ) THEN
        global%inDir = TRIM(global%inDir)//'/'
      END IF ! global%inDir
    END IF ! global%error

    READ(IF_CONTROL,'(A)',IOSTAT=errorFlag) global%outDir
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      global%outDir = './'
    ELSE
      IF ( global%outDir(LEN_TRIM(global%outDir): & 
                         LEN_TRIM(global%outDir)) /=  '/' ) THEN
        global%outDir = TRIM(global%outDir)//'/'
      END IF ! global%outDir
    END IF ! global%error

    READ(IF_CONTROL,'(A)',IOSTAT=errorFlag) global%outDirHDF
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      global%outDirHDF = './'
    ELSE
      IF ( global%outDirHDF(LEN_TRIM(global%outDirHDF): & 
                            LEN_TRIM(global%outDirHDF)) /=  '/' ) THEN
        global%outDirHDF = TRIM(global%outDirHDF)//'/'
      END IF ! global%outDirHDF
    END IF ! global%error

! ------------------------------------------------------------------------------
!   Verbosity and checking levels. NOTE do not report error if Roccom verbosity
!   level does not exist for backward compatibility. NOTE cannot write warning
!   if Roccom verbosity level does not exit because global%myProcid is not set 
!   yet.
! ------------------------------------------------------------------------------

    READ(IF_CONTROL,*,IOSTAT=errorFlag) global%verbLevel
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      global%verbLevel = VERBOSE_HIGH
    ELSE
      global%verbLevel = MAX(global%verbLevel,VERBOSE_NONE)
      global%verbLevel = MIN(global%verbLevel,VERBOSE_HIGH)
    END IF ! global%error

    READ(IF_CONTROL,*,IOSTAT=errorFlag) global%checkLevel
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      global%checkLevel = CHECK_HIGH
    ELSE
      global%checkLevel = MAX(global%checkLevel,CHECK_NONE)
      global%checkLevel = MIN(global%checkLevel,CHECK_HIGH)
    END IF ! global%error

    READ(IF_CONTROL,*,IOSTAT=errorFlag) verbLevelCOM
    IF ( errorFlag == ERR_NONE ) THEN
      global%verbLevelCOM = verbLevelCOM
    END IF ! global%error

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(IF_CONTROL,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_GENX_ReadCtrlFile
 






! ******************************************************************************
!
! Purpose: Register flow data.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_RegisterDataFlow(pRegion)
  
    USE ModMixture, ONLY: t_mixt_input
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winName
    INTEGER :: paneId,sz,ng
    REAL(RFREAL), POINTER :: pReal    
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_mixt_input), POINTER :: pMixtInput
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

    pGrid      => pRegion%grid
    pMixtInput => pRegion%mixtInput    

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering flow data...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Register data
! ******************************************************************************

! ==============================================================================  
!   Volume
! ==============================================================================  

    winName = global%volWinName       

    CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,0,paneId) 

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'Window name:',TRIM(winName)
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Pane id:', paneId           
    END IF ! global%verbLevel
   
! ------------------------------------------------------------------------------
!   Conserved variables
! ------------------------------------------------------------------------------   
   
    pReal => pRegion%mixt%cv(CV_MIXT_DENS,1)
  
    CALL COM_set_array(TRIM(winName)//'.rhof',paneId,pReal,pMixtInput%nCv, &
                       pGrid%nCellsTot)

 
    pReal => pRegion%mixt%cv(CV_MIXT_XMOM,1)
   
    CALL COM_set_array(TRIM(winName)//'.rhovf',paneId,pReal,pMixtInput%nCv, &
                       pGrid%nCellsTot)
 
    pReal => pRegion%mixt%cv(CV_MIXT_ENER,1)
   
    CALL COM_set_array(TRIM(winName)//'.rhoEf',paneId,pReal,pMixtInput%nCv, &
                       pGrid%nCellsTot)                       

! ------------------------------------------------------------------------------
!   Dependent variables
! ------------------------------------------------------------------------------   
   
    pReal => pRegion%mixt%dv(DV_MIXT_PRES,1)
   
    CALL COM_set_array(TRIM(winName)//'.pf',paneId,pReal,pMixtInput%nDv, &
                       pGrid%nCellsTot)
 
    pReal => pRegion%mixt%dv(DV_MIXT_TEMP,1)
   
    CALL COM_set_array(TRIM(winName)//'.Tf',paneId,pReal,pMixtInput%nDv, &
                       pGrid%nCellsTot)

    pReal => pRegion%mixt%dv(DV_MIXT_SOUN,1)
   
    CALL COM_set_array(TRIM(winName)//'.af',paneId,pReal,pMixtInput%nDv, &
                       pGrid%nCellsTot)
   
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering flow data done.'
    END IF ! global%verbLevel 
    
    CALL COM_get_size(TRIM(winName)//'.rhof',paneId,sz,ng)
  
  END SUBROUTINE RFLU_GENX_RegisterDataFlow







! ******************************************************************************
!
! Purpose: Register grid speeds.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_RegisterDataGSpeeds(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winName
    INTEGER :: iPatch,paneId
    REAL(RFREAL), POINTER :: pReal    
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

    pGrid => pRegion%grid

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering grid speeds...'
    END IF ! global%verbLevel
        
! ******************************************************************************
!   Register data
! ******************************************************************************

! ==============================================================================  
!   Volume
! ==============================================================================  

    winName = global%volWinName       

    CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,0,paneId) 

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Volume data...'
      WRITE(STDOUT,'(A,5X,A,1X,A)') SOLVER_NAME,'Window name:',TRIM(winName)
      WRITE(STDOUT,'(A,5X,A,1X,I5.5)') SOLVER_NAME,'Pane id:', paneId           
    END IF ! global%verbLevel
   
    pReal => pGrid%gs(1)

    CALL COM_set_size(TRIM(winName)//'.gs',paneId,pGrid%nFacesTot,0)
    CALL COM_set_array(TRIM(winName)//'.gs',paneId,pReal,1,pGrid%nFacesTot)

! ==============================================================================  
!   Surface
! ==============================================================================  

    winName = global%surfWinName

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Surface data...'
      WRITE(STDOUT,'(A,5X,A,1X,A)') SOLVER_NAME,'Window name:',TRIM(winName)   
    END IF ! global%verbLevel
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,iPatch,paneId) 

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        WRITE(STDOUT,'(A,5X,A,1X,I2)')   SOLVER_NAME,'Patch:',iPatch
        WRITE(STDOUT,'(A,7X,A,1X,I5.5)') SOLVER_NAME,'Pane id:', paneId      
      END IF ! global%verbLevel
      
      pReal => pPatch%gs(1)
            
      CALL COM_set_array(TRIM(winName)//'.gs',paneId,pReal,1,pPatch%nBFacesTot)
    END DO ! iPatch    
   
! ******************************************************************************
!   End
! ******************************************************************************
 
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering grid speeds done.'
    END IF ! global%verbLevel 
  
  END SUBROUTINE RFLU_GENX_RegisterDataGSpeeds







! ******************************************************************************
!
! Purpose: Register interface data.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_RegisterDataInterf(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winName
    INTEGER :: iPatch,paneId
    REAL(RFREAL), POINTER :: pReal    
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

    pGrid => pRegion%grid

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering interface data...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Register interface data
! ******************************************************************************

    winName = global%surfWinName

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'Window name:',TRIM(winName)   
    END IF ! global%verbLevel
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,iPatch,paneId) 

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I2)')   SOLVER_NAME,'Patch:',iPatch
        WRITE(STDOUT,'(A,5X,A,1X,I5.5)') SOLVER_NAME,'Pane id:',paneId      
      END IF ! global%verbLevel
      
! ==============================================================================  
!     Input data
! ==============================================================================  

! ------------------------------------------------------------------------------
!     All patches
! ------------------------------------------------------------------------------
      
      CALL COM_set_array(TRIM(winName)//'.du_alp',paneId,pPatch%duAlp)

! ------------------------------------------------------------------------------
!     Interacting patches (burning or non-burning)      
! ------------------------------------------------------------------------------

      IF ( pPatch%bcCoupled /= BC_NOT_COUPLED ) THEN 
        CALL COM_set_array(TRIM(winName)//'.rhofvf_alp',paneId,pPatch%rhofvfAlp)
        CALL COM_set_array(TRIM(winName)//'.Tb_alp',paneId,pPatch%tbAlp)
      END IF ! pPatch%bcCoupled

! ------------------------------------------------------------------------------
!     Burning patches 
! ------------------------------------------------------------------------------

      IF ( pPatch%bcCoupled == BC_BURNING ) THEN  
        CALL COM_set_array(TRIM(winName)//'.mdot_alp',paneId,pPatch%mdotAlp)        
        CALL COM_set_array(TRIM(winName)//'.Tflm_alp',paneId,pPatch%tflmAlp)
      END IF ! pPatch%bcCoupled

! ==============================================================================  
!     Output data
! ==============================================================================  

! ------------------------------------------------------------------------------
!     Interacting patches (burning or non-burning)      
! ------------------------------------------------------------------------------

      IF ( pPatch%bcCoupled /= BC_NOT_COUPLED ) THEN         
        CALL COM_set_array(TRIM(winName)//'.nf_alp',paneId,pPatch%nfAlp)        
        CALL COM_set_array(TRIM(winName)//'.rhof_alp',paneId,pPatch%rhofAlp)        
        CALL COM_set_array(TRIM(winName)//'.pf',paneId,pPatch%pf)
        CALL COM_set_array(TRIM(winName)//'.tf',paneId,pPatch%tracf)
        CALL COM_set_array(TRIM(winName)//'.qc',paneId,pPatch%qc)
        CALL COM_set_array(TRIM(winName)//'.qr',paneId,pPatch%qr)        
      END IF ! pPatch%bcCoupled

! ------------------------------------------------------------------------------
!     Burning patches 
! ------------------------------------------------------------------------------

      IF ( pPatch%bcCoupled == BC_BURNING ) THEN         
        CALL COM_set_array(TRIM(winName)//'.Tf',paneId,pPatch%tempf)    
        CALL COM_set_array(TRIM(winName)//'.bflag',paneId,pPatch%bFlag)        
      END IF ! pPatch%bcCoupled                  
    END DO ! iPatch    
   
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering interface data done.'
    END IF ! global%verbLevel
  
  END SUBROUTINE RFLU_GENX_RegisterDataInterf






  
! ******************************************************************************
!
! Purpose: Statistics data registration.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_RegisterDataStats(pRegion)
#ifdef STATS  
    USE STAT_RFLU_ModRocstarAdmin, ONLY: STAT_RFLU_GenxRegisterData
#endif
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion
    
! ******************************************************************************
!   Register data
! ******************************************************************************
#ifdef STATS
    CALL STAT_RFLU_GenxRegisterData(pRegion)
#endif  
! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_RegisterDataStats






  
! ******************************************************************************
!
! Purpose: Turbulence data registration.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_RegisterDataTurb(pRegion)
#ifdef TURB  
    USE TURB_RFLU_ModRocstarAdmin, ONLY: TURB_RFLU_GenxRegisterData
#endif
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion
    
! ******************************************************************************
!   Register data
! ******************************************************************************
#ifdef TURB
    IF (pRegion%global%turbActive) CALL TURB_RFLU_GenxRegisterData(pRegion)
#endif
! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_RegisterDataTurb






  
! ******************************************************************************
!
! Purpose: Wrapper for data registration routines.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_RegisterDataWrapper(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion
    
! ******************************************************************************
!   Register data
! ******************************************************************************

! ==============================================================================  
!   Mixture
! ==============================================================================  

    CALL RFLU_GENX_RegisterDataFlow(pRegion)
    CALL RFLU_GENX_RegisterDataGSpeeds(pRegion) 
    CALL RFLU_GENX_RegisterDataInterf(pRegion)
    CALL RFLU_GENX_RegisterDisp(pRegion)  

! ==============================================================================  
!   Statistics
! ==============================================================================  

    CALL RFLU_GENX_RegisterDataStats(pRegion)

! ==============================================================================  
!   Turbulence
! ==============================================================================  

    CALL RFLU_GENX_RegisterDataTurb(pRegion)
  
! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_RegisterDataWrapper






! ******************************************************************************
!
! Purpose: Register volume grid.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!    1. Must be called after new dataitems were created.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_RegisterDisp(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winName
    INTEGER :: paneId  
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

    pGrid => pRegion%grid

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering displacements...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Register displacements
! ******************************************************************************

    winName = global%volWinName
    
    CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,0,paneId) 

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN         
      WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'Window name:',TRIM(winName)
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Pane id:',paneId   
    END IF ! global%verbLevel

    CALL COM_set_array(TRIM(winName)//'.disp',paneId,pGrid%disp)
   
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering displacements done.'
    END IF ! global%verbLevel
  
  END SUBROUTINE RFLU_GENX_RegisterDisp








! ******************************************************************************
!
! Purpose: Register grid.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!    1. Must be called after new dataitems were created.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_RegisterGrid(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
     
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

    pGrid => pRegion%grid

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering grid...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Register grid
! ******************************************************************************

! ==============================================================================  
!   Volume
! ==============================================================================  

    CALL RFLU_GENX_RegisterGridVol(pRegion)
   
! ==============================================================================  
!   Surface 
! ==============================================================================  
    
    IF ( pGrid%nPatches > 0 ) THEN 
      CALL RFLU_GENX_RegisterGridSurf(pRegion)
    END IF ! pGrid%nPatches

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering grid done.'
    END IF ! global%verbLevel
  
  END SUBROUTINE RFLU_GENX_RegisterGrid  








! ******************************************************************************
!
! Purpose: Register surface grid.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!    1. Must be called after new dataitems were created.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_RegisterGridSurf(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winName
    INTEGER :: iPatch,paneId
    INTEGER, POINTER :: pInt    
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

    pGrid => pRegion%grid

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering surface grid...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Register surface grid
! ******************************************************************************
       
    winName = global%surfWinName

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'Window name:',TRIM(winName)   
    END IF ! global%verbLevel
   
! ==============================================================================  
!   Loop over patches 
! ==============================================================================  
 
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      
      CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,iPatch,paneId) 

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I2)')   SOLVER_NAME,'Patch:',iPatch
        WRITE(STDOUT,'(A,5X,A,1X,I2)')   SOLVER_NAME,'Global patch id:', & 
                                         pPatch%iPatchGlobal
        WRITE(STDOUT,'(A,5X,A,1X,I5.5)') SOLVER_NAME,'Pane id:', paneId      
        WRITE(STDOUT,'(A,5X,A,1X,2X,I3)') SOLVER_NAME,'bcFlag:   ', & 
                                          pPatch%bcFlag(1)
        WRITE(STDOUT,'(A,5X,A,1X,2X,I4)') SOLVER_NAME,'cnstrType:', & 
                                          pPatch%cnstrType(1)   
      END IF ! global%verbLevel

! ------------------------------------------------------------------------------
!     Patch quantities
! ------------------------------------------------------------------------------  

      CALL COM_set_size(TRIM(winName)//'.bcflag',paneId,1)
      CALL COM_set_array(TRIM(winName)//'.bcflag',paneId,pPatch%bcFlag)

      CALL COM_set_size(TRIM(winName)//'.patchNo',paneId,1)
      CALL COM_set_array(TRIM(winName)//'.patchNo',paneId,pPatch%patchNo)

      CALL COM_set_size(TRIM(winName)//'.cnstr_type',paneId,1)
      CALL COM_set_array(TRIM(winName)//'.cnstr_type',paneId,pPatch%cnstrType)

! ------------------------------------------------------------------------------
!     Coordinates
! ------------------------------------------------------------------------------  

      CALL COM_set_size(TRIM(winName)//'.nc',paneId,pPatch%nBVertTot, &
                        pPatch%nBVertTot-pPatch%nBVert)
      CALL COM_set_array(TRIM(winName)//'.nc',paneId,pPatch%xyz)
      
! ------------------------------------------------------------------------------
!     Connectivity. IMPORTANT: Note that locally-numbered boundary-face lists 
!     are registered as connectivity. Globally-numbered boundary-face lists are
!     registered as additional data because they are needed for the basic 
!     grid description. 
! ------------------------------------------------------------------------------  

! --- Actual faces -------------------------------------------------------------

      IF ( pPatch%nBTris > 0 ) THEN 
        pInt => pPatch%bTri2vLoc(1,1)

        CALL COM_set_size(TRIM(winName)//'.:t3:real',paneId,pPatch%nBTris,0)
        CALL COM_set_array(TRIM(winName)//'.:t3:real',paneId,pInt) 
        
        pInt => pPatch%bTri2v(1,1)

        CALL COM_set_size(TRIM(winName)//'.t3g:real',paneId,pPatch%nBTris,0)
        CALL COM_set_array(TRIM(winName)//'.t3g:real',paneId,pInt)                                    
      END IF ! pPatch%nBTris

      IF ( pPatch%nBQuads > 0 ) THEN 
        pInt => pPatch%bQuad2vLoc(1,1)

        CALL COM_set_size(TRIM(winName)//'.:q4:real',paneId,pPatch%nBQuads,0)
        CALL COM_set_array(TRIM(winName)//'.:q4:real',paneId,pInt)                             

        pInt => pPatch%bQuad2v(1,1)

        CALL COM_set_size(TRIM(winName)//'.q4g:real',paneId,pPatch%nBQuads,0)
        CALL COM_set_array(TRIM(winName)//'.q4g:real',paneId,pInt)                             
      END IF ! pPatch%nBTris      
      
! --- Virtual faces ------------------------------------------------------------

      IF ( pPatch%nBTrisTot > pPatch%nBTris ) THEN 
        pInt => pPatch%bTri2vLoc(1,pPatch%nBTris+1)

        CALL COM_set_size(TRIM(winName)//'.:t3:virtual',paneId, & 
                          pPatch%nBTrisTot-pPatch%nBTris, &
                          pPatch%nBTrisTot-pPatch%nBTris)
        CALL COM_set_array(TRIM(winName)//'.:t3:virtual',paneId,pInt)
        
        pInt => pPatch%bTri2v(1,pPatch%nBTris+1)

        CALL COM_set_size(TRIM(winName)//'.t3g:virtual',paneId, & 
                          pPatch%nBTrisTot-pPatch%nBTris, &
                          pPatch%nBTrisTot-pPatch%nBTris)
        CALL COM_set_array(TRIM(winName)//'.t3g:virtual',paneId,pInt)                                     
      END IF ! pPatch%nBTris

      IF ( pPatch%nBQuadsTot > pPatch%nBQuads ) THEN 
        pInt => pPatch%bQuad2vLoc(1,pPatch%nBQuads+1)

        CALL COM_set_size(TRIM(winName)//'.:q4:virtual',paneId, &
                          pPatch%nBQuadsTot-pPatch%nBQuads, &
                          pPatch%nBQuadsTot-pPatch%nBQuads)
        CALL COM_set_array(TRIM(winName)//'.:q4:virtual',paneId,pInt) 
        
        pInt => pPatch%bQuad2v(1,pPatch%nBQuads+1)

        CALL COM_set_size(TRIM(winName)//'.q4g:virtual',paneId, &
                          pPatch%nBQuadsTot-pPatch%nBQuads, &
                          pPatch%nBQuadsTot-pPatch%nBQuads)
        CALL COM_set_array(TRIM(winName)//'.q4g:virtual',paneId,pInt)                                    
      END IF ! pPatch%nBQuads      
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering surface grid done.'
    END IF ! global%verbLevel
  
  END SUBROUTINE RFLU_GENX_RegisterGridSurf
  
  
  
  



! ******************************************************************************
!
! Purpose: Register volume grid.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!    1. Must be called after new dataitems were created.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_RegisterGridVol(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winName
    INTEGER :: paneId
    INTEGER, POINTER :: pInt    
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

    pGrid => pRegion%grid

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering volume grid...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Register grid
! ******************************************************************************

    winName = global%volWinName
    
    CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,0,paneId) 

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN         
      WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'Window name:',TRIM(winName)
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Pane id:',paneId   
    END IF ! global%verbLevel

! ==============================================================================  
!   Coordinates
! ==============================================================================  

    CALL COM_set_size(TRIM(winName)//'.nc',paneId,pGrid%nVertTot, &
                      pGrid%nVertTot-pGrid%nVert)
    CALL COM_set_array(TRIM(winName)//'.nc',paneId,pGrid%xyz)

! ==============================================================================  
!   Connectivity
! ==============================================================================  

! ------------------------------------------------------------------------------
!   pconn dataitem
! ------------------------------------------------------------------------------

    CALL COM_set_size(TRIM(winName)//'.pconn',paneId,pGrid%pconnSizeTot, &
                      pGrid%pconnSizeGhost)
    CALL COM_set_array(TRIM(winName)//'.pconn',paneId,pGrid%pconn)    

! ------------------------------------------------------------------------------
!   Actual cells
! ------------------------------------------------------------------------------

    IF ( pGrid%nTets > 0 ) THEN 
      pInt => pGrid%tet2v(1,1)

      CALL COM_set_size(TRIM(winName)//'.:T4:real',paneId,pGrid%nTets,0)
      CALL COM_set_array(TRIM(winName)//'.:T4:real',paneId,pInt)  
    END IF ! pGrid%nTets 

    IF ( pGrid%nHexs > 0 ) THEN 
      pInt => pGrid%hex2v(1,1)

      CALL COM_set_size(TRIM(winName)//'.:H8:real',paneId,pGrid%nHexs,0)
      CALL COM_set_array(TRIM(winName)//'.:H8:real',paneId,pInt)  
    END IF ! pGrid%nHexs 

    IF ( pGrid%nPris > 0 ) THEN 
      pInt => pGrid%pri2v(1,1)

      CALL COM_set_size(TRIM(winName)//'.:W6:real',paneId,pGrid%nPris,0)
      CALL COM_set_array(TRIM(winName)//'.:W6:real',paneId,pInt)  
    END IF ! pGrid%nPris 

    IF ( pGrid%nPyrs > 0 ) THEN 
      pInt => pGrid%pyr2v(1,1)

      CALL COM_set_size(TRIM(winName)//'.:P5:real',paneId,pGrid%nPyrs,0)
      CALL COM_set_array(TRIM(winName)//'.:P5:real',paneId,pInt)  
    END IF ! pGrid%nPyrs     

! ------------------------------------------------------------------------------
!   Virtual cells
! ------------------------------------------------------------------------------

    IF ( pGrid%nTetsTot > pGrid%nTets ) THEN 
      pInt => pGrid%tet2v(1,pGrid%nTets+1)

      CALL COM_set_size(TRIM(winName)//'.:T4:virtual',paneId, &
                        pGrid%nTetsTot-pGrid%nTets,pGrid%nTetsTot-pGrid%nTets)
      CALL COM_set_array(TRIM(winName)//'.:T4:virtual',paneId,pInt) 
    END IF ! pGrid%nTetsTot

    IF ( pGrid%nHexsTot > pGrid%nHexs ) THEN 
      pInt => pGrid%hex2v(1,pGrid%nHexs+1)

      CALL COM_set_size(TRIM(winName)//'.:H8:virtual',paneId, &
                        pGrid%nHexsTot-pGrid%nHexs,pGrid%nHexsTot-pGrid%nHexs)
      CALL COM_set_array(TRIM(winName)//'.:H8:virtual',paneId,pInt)  
    END IF ! pGrid%nHexsTot

    IF ( pGrid%nPrisTot > pGrid%nPris ) THEN 
      pInt => pGrid%pri2v(1,pGrid%nPris+1)

      CALL COM_set_size(TRIM(winName)//'.:W6:virtual',paneId, &
                        pGrid%nPrisTot-pGrid%nPris,pGrid%nPrisTot-pGrid%nPris)
      CALL COM_set_array(TRIM(winName)//'.:W6:virtual',paneId,pInt)  
    END IF ! pGrid%nPrisTot

    IF ( pGrid%nPyrsTot > pGrid%nPyrs ) THEN 
      pInt => pGrid%pyr2v(1,pGrid%nPyrs+1)

      CALL COM_set_size(TRIM(winName)//'.:P5:virtual',paneId, &
                        pGrid%nPyrsTot-pGrid%nPyrs,pGrid%nPyrsTot-pGrid%nPyrs)
      CALL COM_set_array(TRIM(winName)//'.:P5:virtual',paneId,pInt)  
    END IF ! pGrid%nPrisTot

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering volume grid done.'
    END IF ! global%verbLevel
  
  END SUBROUTINE RFLU_GENX_RegisterGridVol

  








! ******************************************************************************
!
! Purpose: Register size of connectivity tables.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!    1. Must be called in initializor if have not registered grid so Roccom
!       knows how many cells exist in grid.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_SetConnSize(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winName
    INTEGER :: paneId   
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

    pGrid => pRegion%grid

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting conn size...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Register conn size
! ******************************************************************************

    winName = global%volWinName
    
    CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,0,paneId) 

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN         
      WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'Window name:',TRIM(winName)
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Pane id:',paneId   
    END IF ! global%verbLevel

    CALL COM_set_size(TRIM(winName)//'.conn',paneId,pGrid%nCellsTot, &
                      pGrid%nCellsTot-pGrid%nCells)   

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Setting conn size done.'
    END IF ! global%verbLevel
  
  END SUBROUTINE RFLU_GENX_SetConnSize













! ******************************************************************************
!
! Purpose: Store communicator.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   communicator        Communicator
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_StoreCommunicator(global,communicator)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    INTEGER, INTENT(IN) :: communicator
    TYPE(t_global), POINTER :: global
        
! ******************************************************************************
!   Store communicator
! ******************************************************************************

    global%communicator = communicator
 
! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE RFLU_GENX_StoreCommunicator










! ******************************************************************************
!
! Purpose: Set input window names.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Input volume window name is assumed to be subdivided into separate 
!      strings for the mixture and multi-physics modules by empty spaces. 
!
! ******************************************************************************

  SUBROUTINE RFLU_GENX_StoreNamesHandles(global,surfWinNameInput, & 
                                         volWinNameInput,handleObtain)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    CHARACTER(*), INTENT(IN) :: surfWinNameInput,volWinNameInput
    INTEGER, INTENT(IN) :: handleObtain
    TYPE(t_global), POINTER :: global

! ==============================================================================  
!   Locals 
! ==============================================================================  
        
    CHARACTER(CHRLEN) :: tempString1,tempString2
    INTEGER :: errorFlag
        
! ******************************************************************************
!   Set names 
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Storing names and handles...'
    END IF ! global%myProcid

! ==============================================================================  
!   Surface
! ==============================================================================  

    global%surfWinNameInput = surfWinNameInput

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'Input surface window name:', & 
                                    TRIM(global%surfWinNameInput)
    END IF ! global%myProcid
    
! ==============================================================================  
!   Volume. NOTE separate input volume window into substrings.
! ==============================================================================  
  
    READ(volWinNameInput,*,IOSTAT=errorFlag) tempString1,tempString2
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_STRING_READ,__LINE__)
    END IF ! global%error

! ------------------------------------------------------------------------------
!   Mixture
! ------------------------------------------------------------------------------
  
    global%volWinNameInput = TRIM(tempString1)

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'Input volume window name:', & 
                                    TRIM(global%volWinNameInput)
    END IF ! global%myProcid
    
! TO DO 
! When PLAG is being integrated, a new window name will need to be defined and
! it will need to be set to tempString2
! END TO DO 
      
! ******************************************************************************
!   Set handles
! ******************************************************************************

    global%handleObtain = handleObtain   

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Storing names and handles done.'
    END IF ! global%myProcid
      
  END SUBROUTINE RFLU_GENX_StoreNamesHandles









END MODULE RFLU_ModRocstarAdmin

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModGENXAdmin.F90,v $
! Revision 1.26  2009/05/12 20:20:56  mtcampbe
! Rocon integration changes
!
! Revision 1.25  2008/12/06 08:44:21  mtcampbe
! Updated license.
!
! Revision 1.24  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.23  2007/07/08 21:45:03  gzheng
! changed the PRESENT is used for PGI compiler
!
! Revision 1.22  2007/04/20 16:07:48  mtcampbe
! Updating for burnout support function RFLU_GENX_InitBFLAG
!
! Revision 1.21  2007/04/14 14:12:53  mtcampbe
! Mods for TZ
!
! Revision 1.20  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.19  2006/01/04 20:11:18  wasistho
! added turbActive condition
!
! Revision 1.18  2006/01/03 06:29:53  wasistho
! added CreateAttrWrapper and RegisterDataWrapper
!
! Revision 1.17  2005/10/14 14:06:53  haselbac
! Added alloc/dealloc and setting of tbAlp array
!
! Revision 1.16  2005/10/14 13:07:36  haselbac
! Removed tbAlp - added erroneously when checking in other changes
!
! Revision 1.15  2005/10/13 17:27:32  haselbac
! Fixed bug in format statement, lead to *** output
!
! Revision 1.14  2005/07/01 15:14:36  haselbac
! Added reading of verbLevelCOM
!
! Revision 1.13  2005/06/09 20:19:52  haselbac
! Added cnstr_type to various routines
!
! Revision 1.12  2005/05/04 03:34:10  haselbac
! Added writing of more info when registering surface grid
!
! Revision 1.11  2005/04/15 15:06:55  haselbac
! Added routines to create, build, and destroy pconn dataitem
!
! Revision 1.10  2005/03/09 15:07:41  haselbac
! Erroneous check-in, removed debug code
!
! Revision 1.9  2005/03/09 15:06:03  haselbac
! Added 2d option
!
! Revision 1.8  2004/11/03 17:02:16  haselbac
! Removal of vertex and cell flag arrays, removed RFLU_GENX_CreateAttrGridVol (no longer necessary)
!
! Revision 1.7  2004/11/02 14:03:34  haselbac
! Bug fix: Added missing declaration of t_mixt_input
!
! Revision 1.6  2004/11/02 02:31:33  haselbac
! Replaced CV_MIXT_NEQS and DV_MIXT_NVAR
!
! Revision 1.5  2004/10/27 12:26:52  jiao
! Fixed definition of bcflag, which used to be defined three times.
!
! Revision 1.4  2004/10/27 05:46:45  jiao
! Added COM_set_size back for .gs for volume window.
!
! Revision 1.3  2004/10/27 05:24:55  jiao
! Removed COM_set_size on du_alp and gs.
!
! Revision 1.2  2004/10/22 14:00:57  haselbac
! Bug fix: False type used when creating attr for disp
!
! Revision 1.1  2004/10/19 19:27:25  haselbac
! Initial revision
!
! ******************************************************************************
















