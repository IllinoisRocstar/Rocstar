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
! Purpose: Driver routine for postprocessor.
!
! Description: None.
!
! Input:
!   caseString  String with casename
!   stampString String with iteration or time stamp
!   verbLevel   Verbosity level
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: rflupost.F90,v 1.10 2008/12/06 08:45:06 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE rflupost(caseString,stampString,verbLevel)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global 
  USE ModError
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModParameters
  USE ModMPI
  
  USE RFLU_ModRegionMapping
  
  USE ModInterfaces, ONLY: RFLU_BuildDataStruct, &
                           RFLU_GetUserInput, &
                           RFLU_InitGlobal, &
                           RFLU_MergePostProcessRegions, & 
                           RFLU_PostProcessRegions, &
                           RFLU_PostProcessRegions_ENS, & 
                           RFLU_PrintHeader, &
                           RFLU_PrintWarnInfo, &
                           RFLU_ReadInputFile, &
                           RFLU_SetModuleType, & 
                           RFLU_WriteVersionString

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments      
! ==============================================================================

  CHARACTER(*) :: caseString,stampString
  INTEGER, INTENT(IN) :: verbLevel

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: casename,nRegions,RCSIdentString,stamp
  INTEGER :: errorFlag,iReg
  TYPE(t_region), POINTER :: pRegion
  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER :: levels(:)

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: rflupost.F90,v $ $Revision: 1.10 $'

! ******************************************************************************
! Initialize global data
! ******************************************************************************
  
  ALLOCATE(global,STAT=errorFlag)
  IF ( errorFlag /= ERR_NONE ) THEN 
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR - Pointer allocation failed.'
    STOP
  END IF ! errorFlag   
  
  casename = caseString(1:LEN(caseString))
  stamp    = stampString(1:LEN(stampString))

  CALL RFLU_InitGlobal(casename,verbLevel,CRAZY_VALUE_INT,global)

  CALL RegisterFunction(global,'rflupost', &
                        'rflupost.F90')

  CALL RFLU_SetModuleType(global,MODULE_TYPE_POSTPROC)

! ******************************************************************************
! Print header and write version string
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC ) THEN
    CALL RFLU_WriteVersionString(global)     
    IF ( global%verbLevel /= VERBOSE_NONE ) THEN
      CALL RFLU_PrintHeader(global)
    END IF ! global%verbLevel
  END IF ! global%myProcid

! ******************************************************************************
! Read mapping file, impose serial mapping, and build basic data structure
! ******************************************************************************

  CALL RFLU_ReadRegionMappingFile(global,MAPFILE_READMODE_PEEK,global%myProcid)
  CALL RFLU_SetRegionMappingSerial(global)
  CALL RFLU_CreateRegionMapping(global,MAPTYPE_REG)
  CALL RFLU_ImposeRegionMappingSerial(global)
  
  CALL RFLU_BuildDataStruct(global,levels)
  CALL RFLU_ApplyRegionMapping(global,levels)
  CALL RFLU_DestroyRegionMapping(global,MAPTYPE_REG)
  
! ******************************************************************************  
! Read input file
! ******************************************************************************

  CALL RFLU_GetUserInput(levels(1)%regions)  

  IF ( global%flowType == FLOW_STEADY ) THEN
    READ(stamp,*) global%currentIter
  ELSE
    READ(stamp,*) global%currentTime
  END IF ! global%flowType  

! ******************************************************************************
! Check for TECPLOT library. NOTE if not available, continue anyway if have 
! solution initialized from hardcode, because this is likely to mean that errors
! will be computed in RFLU_PostProcess* routines.
! ******************************************************************************

#ifdef NO_TECPLOT
  IF ( global%postOutputFormat == POST_OUTPUT_FORMAT_TECPLOT ) THEN
    IF ( global%initFlowFlag /= INITFLOW_FROMHARDCODE ) THEN 
      CALL ErrorStop(global,ERR_POST_OUTPUT_FORMAT_INVALID,__LINE__)
    ELSE 
      global%warnCounter = global%warnCounter + 1
      WRITE(STDOUT,'(A,1X,A,1X,A)') SOLVER_NAME,'*** WARNING ***', &
            'TECPLOT library not available. Continuing anyway.'
    END IF ! global%initFlowFlag
  END IF ! global%postOutputFormat
#endif

! ******************************************************************************
! Call appropriate routine to write output
! ******************************************************************************

  IF ( global%postMergeFlag .EQV. .TRUE.) THEN
    IF ( global%nRegionsLocal > 1 ) THEN  
      CALL RFLU_MergePostProcessRegions(levels)
    ELSE 
      SELECT CASE ( global%postOutputFormat ) 
        CASE ( POST_OUTPUT_FORMAT_TECPLOT )
          CALL RFLU_PostProcessRegions(levels)
        CASE ( POST_OUTPUT_FORMAT_ENSIGHT )
          CALL RFLU_PostProcessRegions_ENS(levels)
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)        
      END SELECT ! global%postOutputFormat
    END IF ! global%nRegionsLocal
  ELSE 
    SELECT CASE ( global%postOutputFormat ) 
      CASE ( POST_OUTPUT_FORMAT_TECPLOT )
        CALL RFLU_PostProcessRegions(levels)
      CASE ( POST_OUTPUT_FORMAT_ENSIGHT )
        CALL RFLU_PostProcessRegions_ENS(levels)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)        
    END SELECT ! global%postOutputFormat
  END IF ! global%postMergeFlag
  
! ******************************************************************************
! Print info about warnings
! ******************************************************************************

  CALL RFLU_PrintWarnInfo(global)
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE rflupost

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: rflupost.F90,v $
! Revision 1.10  2008/12/06 08:45:06  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:18:17  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.7  2006/02/06 23:55:55  haselbac
! Added comm argument to RFLU_InitGlobal
!
! Revision 1.6  2005/11/27 02:07:38  haselbac
! Added setting of module type (for reading of EEv data)
!
! Revision 1.5  2005/10/05 20:15:12  haselbac
! Added support for ENSIGHT
!
! Revision 1.4  2005/05/03 03:12:21  haselbac
! Converted to C++ reading of command-line
!
! Revision 1.3  2005/04/15 15:08:38  haselbac
! Removed Charm/FEM stuff
!
! Revision 1.2  2005/01/14 21:02:02  haselbac
! Changed call for peeking at reg map
!
! Revision 1.1  2004/12/29 20:58:54  haselbac
! Initial revision
!
! ******************************************************************************







