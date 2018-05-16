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
! Purpose: Pick regions.
!
! Description: None.
!
! Input: 
!   regions		Region data
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_PickRegions.F90,v 1.5 2008/12/06 08:45:04 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_PickRegions(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI
  USE ModParameters
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER :: infoType
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iReg
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PickRegions.F90,v $ $Revision: 1.5 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_PickRegions',&
  'RFLU_PickRegions.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Picking regions...'
  END IF ! global%verbLevel

! ******************************************************************************
! Enter information
! ******************************************************************************

  WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Enter information on regions:'
  WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'a - Pick all regions'   
  WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'s - Pick some regions'
  WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Enter information type:'
  READ(STDIN,'(A)') infoType
   
  SELECT CASE ( infoType ) 
  
! ==============================================================================
!   All regions
! ==============================================================================
    
    CASE ( 'a' ) 
      DO iReg = 1,global%nRegionsLocal
        regions(iReg)%activeFlag = .TRUE.
      END DO ! iReg    

! ==============================================================================
!   Some regions
! ==============================================================================

    CASE ( 's' )   
      DO iReg = 1,global%nRegionsLocal ! Must initialize first
        regions(iReg)%activeFlag = .FALSE.
      END DO ! iReg      
      
      DO     
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
                                 'Enter global region index (< 1 to exit):'       
        READ(STDIN,*) iReg
        
        IF ( iReg < 1 ) THEN 
          EXIT
        ELSE IF ( iReg > global%nRegionsLocal ) THEN
          global%warnCounter = global%warnCounter + 1         
         
          WRITE(STDOUT,'(A,5X,A,1X,A,I7,1X,A)') SOLVER_NAME, & 
                       '*** WARNING *** Invalid input.', & 
                       'There are only ',global%nRegionsLocal,'regions.'
        END IF ! iReg
        
        regions(iReg)%activeFlag = .TRUE.
      END DO ! <empty>

! ==============================================================================
!   Default
! ==============================================================================
      
    CASE DEFAULT 
      global%warnCounter = global%warnCounter + 1     
    
      WRITE(STDOUT,'(A,5X,A,1X,A)') SOLVER_NAME, &
                                    '*** WARNING *** Invalid input.', & 
                                    'No regions selected.' 
                                    
      DO iReg = 1,global%nRegionsLocal
        regions(iReg)%activeFlag = .FALSE.
      END DO ! iReg             
  END SELECT ! infoType
  
! ******************************************************************************
! End
! ******************************************************************************

  IF ( (global%myProcid == MASTERPROC) .AND. & 
       (global%verbLevel > VERBOSE_NONE) ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Picking regions done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PickRegions

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PickRegions.F90,v $
! Revision 1.5  2008/12/06 08:45:04  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:18:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/06/14 17:48:47  haselbac
! Bug fix: Missing SOLVER_NAME
!
! Revision 1.2  2003/07/22 02:08:33  haselbac
! Added global%warnCounter
!
! Revision 1.1.1.1  2003/06/04 22:31:20  haselbac
! Initial revision
!
!******************************************************************************







