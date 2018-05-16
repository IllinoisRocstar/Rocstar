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
! $Id: TURB_RFLU_ModRocstarAdmin.F90,v 1.7 2008/12/06 08:44:01 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE TURB_RFLU_ModRocstarAdmin

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid 
  USE ModMPI 

  USE RFLU_ModRocstarUtils

  IMPLICIT NONE

  INCLUDE 'comf90.h'

  PRIVATE

  PUBLIC :: TURB_RFLU_GenxCreateAttr, &
            TURB_RFLU_GenxRegisterData, &
            TURB_RFLU_GenxGetData

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================  
! Private
! ==============================================================================  

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: TURB_RFLU_ModRocstarAdmin.F90,v $ $Revision: 1.7 $'       
  
! ==============================================================================  
! Public
! ==============================================================================  
  

! ******************************************************************************
! Contained routines
! ******************************************************************************

  CONTAINS 
  
  
  
! ******************************************************************************
!
! Purpose: Create new dataitems for turbulence solution.
!
! Description: None.
!
! Input:
!   pRegion		Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE TURB_RFLU_GenxCreateAttr(pRegion)
  
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
    
    CHARACTER(CHRLEN) :: winv
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global
  
! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE TURB_RFLU_GenxCreateAttr
  
  
  
  
  


! ******************************************************************************
!
! Purpose: Register turbulence data.
!
! Description: None.
!
! Input:
!   pRegion		Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE TURB_RFLU_GenxRegisterData(pRegion)

    USE TURB_ModParameters

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
    
    CHARACTER(CHRLEN) :: winv
    INTEGER :: paneId, ilb, ndim, nTv, errorFlag
    REAL(RFREAL), POINTER :: pReal
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

    pGrid  => pRegion%grid

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering turbulence data...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Create new dataitems
! ******************************************************************************

! ==============================================================================  
!   Volume
! ==============================================================================  

    winv = global%volWinName

! ------------------------------------------------------------------------------
!   Global
! ------------------------------------------------------------------------------
    CALL COM_new_dataitem(TRIM(winv)//'.esg1Sum','w',COM_DOUBLE,1,'J/(m^3s)')
    CALL COM_new_dataitem(TRIM(winv)//'.esg4Sum','w',COM_DOUBLE,1,'J/(m^3s)')
    CALL COM_set_array( TRIM(winv)//'.esg1Sum' ,0, global%esg1Sum )
    CALL COM_set_array( TRIM(winv)//'.esg4Sum' ,0, global%esg4Sum )

! ------------------------------------------------------------------------------
!   Model and instantaneous variables
! ------------------------------------------------------------------------------

    CALL COM_new_dataitem( TRIM(winv)//'.lens','e',COM_DOUBLE,1,'m' )
    CALL COM_new_dataitem( TRIM(winv)//'.mut' ,'e',COM_DOUBLE,1,'kg/(ms)' )
    CALL COM_new_dataitem( TRIM(winv)//'.vort','e',COM_DOUBLE,3,'1/s' )

    
! ******************************************************************************
!   Register data
! ******************************************************************************

! ==============================================================================  
!   Volume
! ==============================================================================  

    winv = global%volWinName       

    CALL RFLU_GENX_BuildPaneId(pRegion%iRegionGlobal,0,paneId) 

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'Window name:',TRIM(winv)
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Pane id:', paneId           
    END IF ! global%verbLevel
   
! ------------------------------------------------------------------------------
!   Initialisation
! ------------------------------------------------------------------------------   
    IF (pRegion%turbInput%modelClass /= MODEL_RANS) THEN
      ALLOCATE( pRegion%turb%lens(pGrid%nCellsTot),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    ENDIF
    IF (pRegion%mixtInput%turbModel == TURB_MODEL_NONE) THEN
      ALLOCATE( pRegion%turb%vort(3,pGrid%nCellsTot),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    ENDIF

    pRegion%mixt%tv   = 0._RFREAL
    pRegion%turb%lens = 0._RFREAL
    pRegion%turb%vort = 0._RFREAL

! ------------------------------------------------------------------------------
!   Model and instantaneous variables
! ------------------------------------------------------------------------------
    ilb  = 1
    ndim = pGrid%nCellsTot
    nTv  = pRegion%mixtInput%nTv 

    pReal => pRegion%turb%lens(1)
    CALL COM_set_array( TRIM(winv)//'.lens' ,paneId, pReal )

    pReal => pRegion%mixt%tv(TV_MIXT_MUET,ilb)
    CALL COM_set_array( TRIM(winv)//'.mut',paneId, pReal,nTv )

    pReal => pRegion%turb%vort(1,ilb)
    CALL COM_set_array( TRIM(winv)//'.vort',paneId,pReal )
   
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering turbulence data done.'
    END IF ! global%verbLevel 
  
  END SUBROUTINE TURB_RFLU_GenxRegisterData








! ******************************************************************************
!
! Purpose: Get turbulence data through Roccom.
!
! Description: None.
!
! Input:
!   pRegion		Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE TURB_RFLU_GenxGetData(pRegion)
  
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Loop
! ==============================================================================  
    INTEGER :: iStat

! ==============================================================================  
!   Locals 
! ==============================================================================  
    
    CHARACTER(CHRLEN) :: winName,winNameIn
    CHARACTER(CHRLEN), POINTER :: statNm(:,:,:)
    INTEGER :: handleIn,handleObtain,handleOut
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global

    winNameIn = global%volWinNameInput
    winName   = global%volWinName  
    
! ******************************************************************************
!   First check if it is new or restart turbulence. If new, skip this routine.
! ******************************************************************************

    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.mut')
    IF (handleIn <= 0) THEN
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Starting new turbulence...'
      END IF ! global%verbLevel
      GOTO 999
    ENDIF

! ******************************************************************************
!   Get data
! ******************************************************************************

    handleObtain = global%handleObtain

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting turbulence data...'
    END IF ! global%verbLevel 
        
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid
    
! ==============================================================================
!   Global
! ==============================================================================   
    
    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.esg1Sum')
    handleOut = COM_get_dataitem_handle(TRIM(winName)//'.esg1Sum')
    CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    
    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.esg4Sum')
    handleOut = COM_get_dataitem_handle(TRIM(winName)//'.esg4Sum')
    CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    
! ==============================================================================
!   Turbulence variables
! ==============================================================================   
    
    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.lens')
    handleOut = COM_get_dataitem_handle(TRIM(winName)//'.lens')
    CALL COM_call_function(handleObtain,2,handleIn,handleOut)

    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.mut')
    handleOut = COM_get_dataitem_handle(TRIM(winName)//'.mut')
    CALL COM_call_function(handleObtain,2,handleIn,handleOut)

    handleIn  = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.vort')    
    handleOut = COM_get_dataitem_handle(TRIM(winName)//'.vort')
    CALL COM_call_function(handleObtain,2,handleIn,handleOut)


! ******************************************************************************
!   End
! ******************************************************************************
 
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting turbulence data done.'
    END IF ! global%verbLevel

999 CONTINUE
  
  END SUBROUTINE TURB_RFLU_GenxGetData







! ******************************************************************************
! End Module
! ******************************************************************************




END MODULE TURB_RFLU_ModRocstarAdmin

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_RFLU_ModRocstarAdmin.F90,v $
! Revision 1.7  2008/12/06 08:44:01  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/02/11 06:45:06  wasistho
! identify new or restart turbulence in genxGetData
!
! Revision 1.4  2006/01/10 05:00:05  wasistho
! added GenxGetData
!
! Revision 1.3  2006/01/04 20:06:02  wasistho
! modified data registration
!
! Revision 1.2  2006/01/03 09:51:35  wasistho
! get rid of ifdef rflu
!
! Revision 1.1  2006/01/03 06:34:41  wasistho
! initial import
!
!
! ******************************************************************************






