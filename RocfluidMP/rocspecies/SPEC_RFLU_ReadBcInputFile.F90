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
! Purpose: Read in user input related to boundary conditions for species 
!   (done on all processors).
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
! $Id: SPEC_RFLU_ReadBcInputFile.F90,v 1.7 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_ReadBcInputFile(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModError
  USE ModParameters
  USE ModMPI
  
  USE ModBuildFileNames, ONLY: BuildFileNamePlain
  
  USE ModInterfacesSpecies, ONLY: SPEC_RFLU_ReadBcFarfSection, & 
                                  SPEC_RFLU_ReadBcInflowSection, & 
                                  SPEC_RFLU_ReadBcInjectSection, & 
                                  SPEC_RFLU_ReadBcSectionDummy
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: iFileName,RCSIdentString
  CHARACTER(256) :: line
  INTEGER :: errorFlag,iPatch,loopCounter
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_grid) :: grid
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_RFLU_ReadBcInputFile.F90,v $ $Revision: 1.7 $'

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_ReadBcInputFile',&
  'SPEC_RFLU_ReadBcInputFile.F90')

! ******************************************************************************
! Open file
! ******************************************************************************

  CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.bc',iFileName)
  
  OPEN(IF_INPUT,FILE=iFileName,FORM='FORMATTED',STATUS='OLD',IOSTAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! Read file looking for keywords
! ******************************************************************************

  loopCounter = 0

  KeyWordLoop: DO
    READ(IF_INPUT,'(A256)',IOSTAT=errorFlag) line
    
    IF ( errorFlag > 0 ) THEN ! Error occurred
      CALL ErrorStop(global,ERR_FILE_READ,__LINE__,'File: '//TRIM(iFileName))
    ELSE IF ( errorFlag < 0 ) THEN ! Encountered end of file
      EXIT KeyWordLoop
    END IF ! errorFlag    
    
    SELECT CASE( TRIM(line) )
      CASE ('# BC_SLIPW')
        CALL SPEC_RFLU_ReadBcSectionDummy(pRegion)
      CASE ('# BC_NOSLIP')
        CALL SPEC_RFLU_ReadBcSectionDummy(pRegion)      
! TEMPORARY - Keep this for backward compatibility
      CASE ('# BC_INFLOW')
        CALL SPEC_RFLU_ReadBcInflowSection(pRegion)
! END TEMPORARY    
      CASE ('# BC_INFLOW_TOTANG')
        CALL SPEC_RFLU_ReadBcInflowSection(pRegion)
      CASE ('# BC_INFLOW_VELTEMP')
        CALL SPEC_RFLU_ReadBcInflowSection(pRegion)             
      CASE ('# BC_OUTFLOW')
        CALL SPEC_RFLU_ReadBcSectionDummy(pRegion)      
      CASE ('# BC_FARF')
        CALL SPEC_RFLU_ReadBcFarfSection(pRegion)
      CASE ('# BC_INJECT')
        CALL SPEC_RFLU_ReadBcInjectSection(pRegion)
      CASE ('# END')
        EXIT
    END SELECT ! TRIM(line)

    loopCounter = loopCounter + 1 ! Prevent infinite loop
    IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
      CALL ErrorStop(global,ERR_INFINITE_LOOP ,__LINE__)
    END IF ! loopCounter
  END DO KeyWordLoop

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_INPUT,IOSTAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_RFLU_ReadBcInputFile

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_ReadBcInputFile.F90,v $
! Revision 1.7  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/04/07 15:19:25  haselbac
! Removed tabs
!
! Revision 1.4  2005/04/27 02:13:37  haselbac
! Added routine to read INFLOW_VELTEMP
!
! Revision 1.3  2005/03/31 17:20:03  haselbac
! Fixed problem with crashes by reading with dummy routine
!
! Revision 1.2  2004/06/16 20:01:23  haselbac
! Added use of ModBuildFileNames, cosmetics
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
! ******************************************************************************







