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
! $Id: STAT_RFLU_ModRocstarAdmin.F90,v 1.8 2008/12/06 08:44:01 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE STAT_RFLU_ModRocstarAdmin

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

  PUBLIC :: STAT_RFLU_GenxCreateAttr, &
            STAT_RFLU_GenxRegisterData, &
            STAT_RFLU_GenxGetData

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================  
! Private
! ==============================================================================  

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: STAT_RFLU_ModRocstarAdmin.F90,v $ $Revision: 1.8 $'       
  
! ==============================================================================  
! Public
! ==============================================================================  
  

! ******************************************************************************
! Contained routines
! ******************************************************************************

  CONTAINS 
  
  
  
! ******************************************************************************
!
! Purpose: Create new dataitems for time-averaged statistics.
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

  SUBROUTINE STAT_RFLU_GenxCreateAttr(pRegion)
  
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
    
    CHARACTER(CHRLEN) :: winv
    CHARACTER(CHRLEN), POINTER :: statNm(:,:,:)
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global
  
! ******************************************************************************
!   End
! ******************************************************************************
  
  END SUBROUTINE STAT_RFLU_GenxCreateAttr
  
  
  
  
  


! ******************************************************************************
!
! Purpose: Register statistics data.
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

  SUBROUTINE STAT_RFLU_GenxRegisterData(pRegion)
  
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
    
    CHARACTER(CHRLEN) :: winv
    CHARACTER(CHRLEN), POINTER :: statNm(:,:,:)
    INTEGER :: paneId, ilb
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start, set pointers
! ******************************************************************************
    
    global => pRegion%global

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering statistics data...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Create new dataitems
! ******************************************************************************

! ==============================================================================  
!   Volume
! ==============================================================================  

    winv = global%volWinName
    IF ((global%flowType == FLOW_UNSTEADY).AND.(global%doStat == ACTIVE)) THEN

! ------------------------------------------------------------------------------
!     Global
! ------------------------------------------------------------------------------
      CALL COM_new_dataitem( TRIM(winv)//'.tStat','w',COM_DOUBLE,1,'s' )
      CALL COM_set_array( TRIM(winv)//'.tStat',0, global%integrTime )

! ------------------------------------------------------------------------------
!     Mixture
! ------------------------------------------------------------------------------
      IF (global%mixtNStat > 0) THEN
        statNm => global%mixtStatNm
        DO iStat=1,global%mixtNStat
          CALL COM_new_dataitem( TRIM(winv)//'.'//TRIM(statNm(1,1,iStat)),'e',&
                                  COM_DOUBLE,1,TRIM(statNm(1,2,iStat)) )
        ENDDO
      ENDIF

#ifdef TURB
      IF ((global%turbActive .EQV. .true.) .AND. (global%turbNStat > 0)) THEN
        statNm => global%turbStatNm
        DO iStat=1,global%turbNStat
          CALL COM_new_dataitem( TRIM(winv)//'.'//TRIM(statNm(1,1,iStat)),'e',&
                                  COM_DOUBLE,1,TRIM(statNm(1,2,iStat)) )
        ENDDO
      ENDIF  ! turbNStat
#endif
    ENDIF   ! unsteady and dostat

    
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
   
    IF ((global%flowType==FLOW_UNSTEADY) .AND. (global%doStat==ACTIVE)) THEN
      ilb = 1
! ------------------------------------------------------------------------------
!     Mixture statistics
! ------------------------------------------------------------------------------   
      IF (global%mixtNStat > 0) THEN
        statNm => global%mixtStatNm
        DO iStat=1,global%mixtNStat
          CALL COM_set_array( TRIM(winv)//'.'//TRIM(statNm(1,1,iStat)),paneId,&
             pRegion%mixt%tav(iStat,ilb), global%mixtNStat)
        ENDDO
      ENDIF  ! mixtNStat

! ------------------------------------------------------------------------------
!     Turbulence statistics
! ------------------------------------------------------------------------------   
#ifdef TURB
      IF ((global%turbActive .EQV. .true.) .AND. (global%turbNStat > 0)) THEN
        statNm => global%turbStatNm
        DO iStat=1,global%turbNStat
          CALL COM_set_array( TRIM(winv)//'.'//TRIM(statNm(1,1,iStat)),paneId,&
             pRegion%turb%tav(iStat,ilb), global%turbNStat)
        ENDDO
      ENDIF  ! turbNStat
#endif   
    ENDIF    ! flowType and doStat

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Registering statistics data done.'
    END IF ! global%verbLevel 
  
  END SUBROUTINE STAT_RFLU_GenxRegisterData








! ******************************************************************************
!
! Purpose: Get statistics data through Roccom.
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

  SUBROUTINE STAT_RFLU_GenxGetData(pRegion)
  
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

    handleObtain = global%handleObtain

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting statistics data...'
    END IF ! global%verbLevel 
        
    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
    END IF ! global%myProcid

! ******************************************************************************
!   Get data
! ******************************************************************************

    winNameIn = global%volWinNameInput
    winName   = global%volWinName  
    
! ==============================================================================
!   Global
! ==============================================================================   

    handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.tStat')
    handleOut = COM_get_dataitem_handle(TRIM(winName)//'.tStat')
    CALL COM_call_function(handleObtain,2,handleIn,handleOut)
    
! ==============================================================================
!   Statistics variables
! ==============================================================================   
    
! ------------------------------------------------------------------------------
!   Mixture statistics
! ------------------------------------------------------------------------------

    IF (global%mixtNStat > 0) THEN
      statNm => global%mixtStatNm
      DO iStat=1,global%mixtNStat

        handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.'// &
                                                  TRIM(statNm(1,1,iStat)))
        handleOut = COM_get_dataitem_handle(TRIM(winName)//'.'// &
                                             TRIM(statNm(1,1,iStat)))
        CALL COM_call_function(handleObtain,2,handleIn,handleOut)

      ENDDO
    ENDIF  ! mixtNStat

    
! ------------------------------------------------------------------------------
!   Turbulence statistics
! ------------------------------------------------------------------------------

#ifdef TURB
    IF (global%turbNStat > 0) THEN
      statNm => global%turbStatNm
      DO iStat=1,global%turbNStat

        handleIn = COM_get_dataitem_handle_const(TRIM(winNameIn)//'.'// &
                                                  TRIM(statNm(1,1,iStat)))
        handleOut = COM_get_dataitem_handle(TRIM(winName)//'.'// &
                                             TRIM(statNm(1,1,iStat)))
        CALL COM_call_function(handleObtain,2,handleIn,handleOut)

      ENDDO
    ENDIF  ! turbNStat
#endif   
    
! ==============================================================================
!   Set variables back to accumulated mode
! ==============================================================================   
    
! ------------------------------------------------------------------------------
!   Mixture statistics
! ------------------------------------------------------------------------------

    IF (global%mixtNStat > 0) THEN
      DO iStat=1,global%mixtNStat
        pRegion%mixt%tav(iStat,:) = &
        pRegion%mixt%tav(iStat,:)*global%integrTime
      ENDDO
    ENDIF ! mixtNstat
    
! ------------------------------------------------------------------------------
!   Turbulence statistics
! ------------------------------------------------------------------------------

#ifdef TURB
    IF ((global%turbActive .EQV. .true.) .AND. (global%turbNStat > 0)) THEN
      DO iStat=1,global%turbNStat
        pRegion%turb%tav(iStat,:) = &
        pRegion%turb%tav(iStat,:)*global%integrTime
      ENDDO
    ENDIF ! turbNstat
#endif

! ******************************************************************************
!   End
! ******************************************************************************
 
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting statistics data done.'
    END IF ! global%verbLevel   
  
  END SUBROUTINE STAT_RFLU_GenxGetData






! ******************************************************************************
! End Module
! ******************************************************************************




END MODULE STAT_RFLU_ModRocstarAdmin

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: STAT_RFLU_ModRocstarAdmin.F90,v $
! Revision 1.8  2008/12/06 08:44:01  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2006/02/05 03:46:00  wasistho
! added global%turbActive
!
! Revision 1.5  2006/01/10 06:32:39  wasistho
! mixt%tav to turb%tav within turb stats
!
! Revision 1.4  2006/01/10 05:00:15  wasistho
! added GenxGetData
!
! Revision 1.3  2006/01/04 20:06:13  wasistho
! modified data registration
!
! Revision 1.2  2006/01/03 09:51:43  wasistho
! get rid of ifdef rflu
!
! Revision 1.1  2006/01/03 06:34:41  wasistho
! initial import
!
!
! ******************************************************************************






