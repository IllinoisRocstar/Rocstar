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
! Purpose: Read file with with post-processor information.
!
! Description: None.
!
! Input: 
!   pRegion              Pointer to region data
!   readMode             Reading mode
!
! Output: None.
!
! Notes: 
!   1. Need to have two reading modes because for the activation flag to be
!      useful, need to do this before any quantities are read or allocated, 
!      but this means that nCellsSpecial will be overwritten when grid is 
!      created. So read twice, first for activation flag only, and then after 
!      grid is created, read nCellsSpecial and the actual indices of the 
!      special cells.
!
! ******************************************************************************
!
! $Id: RFLU_ReadPostInfo.F90,v 1.7 2008/12/06 08:45:06 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ReadPostInfo(pRegion,readMode)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModMPI
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: readMode
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: dummyLogical
  CHARACTER(CHRLEN) :: dummyString,dummyString2,iRegionString,RCSIdentString
  INTEGER :: dummyInteger,dummyInteger2,ics,iFile,ifs,indx,iReg,iRegionGlobal, &
             nCellsSpecial
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ReadPostInfo.F90,v $ $Revision: 1.7 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ReadPostInfo',&
  'RFLU_ReadPostInfo.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading post-processor info...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal 
  END IF ! global%verbLevel

! ******************************************************************************
! Set grid pointer and variables
! ******************************************************************************

  pGrid => pRegion%grid

  iFile = IF_POSTINFO

! ******************************************************************************
! Rewind file - must be done because stays open while read several times
! ******************************************************************************

  REWIND(iFile)

! ******************************************************************************
! Read from file
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    READ(iFile,'(A)') dummyString

    WRITE(iRegionString,'(I5.5)') pRegion%iRegionGlobal
    indx = INDEX(dummyString,TRIM(iRegionString))
    
    IF ( indx /= 0 ) THEN
      dummyString2 = dummyString(indx:indx+LEN_TRIM(iRegionString)-1)
      READ(dummyString2,*) iRegionGlobal
    ELSE 
      iRegionGlobal = CRAZY_VALUE_INT ! Anything but pRegion%iRegionGlobal
    END IF ! indx
        
! ==============================================================================
!   Found entry for my region
! ==============================================================================

    IF ( iRegionGlobal == pRegion%iRegionGlobal ) THEN
    
! ------------------------------------------------------------------------------
!     Read activation flag only
! ------------------------------------------------------------------------------
    
      IF ( readMode == INFOFILE_READMODE_FLAG ) THEN
        READ(iFile,'(L1)') pRegion%postActiveFlag
        
        IF ( pRegion%postActiveFlag .EQV. .TRUE. ) THEN
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Region is active.'
          END IF ! global%verbLevel          
                         
          READ(iFile,*) dummyInteger ! Special cells

          DO ics = 1,dummyInteger
            READ(iFile,*) dummyInteger2
          END DO ! ics
          
          READ(iFile,*) dummyInteger ! Special faces

          DO ifs = 1,dummyInteger
            READ(iFile,*) dummyInteger2
          END DO ! ifs          
        ELSE 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Region is not active.'
          END IF ! global%verbLevel         
        END IF ! pRegion%postActiveFlag   

! ------------------------------------------------------------------------------
!     Read special cells and faces only
! ------------------------------------------------------------------------------

      ELSE IF ( readMode == INFOFILE_READMODE_DATA ) THEN
        READ(iFile,'(L1)') dummyLogical
            
        IF ( dummyLogical .EQV. .TRUE. ) THEN 
          READ(iFile,*) pGrid%nCellsSpecial

          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
                                           'Number of special cells:', &
                                           pGrid%nCellsSpecial
          END IF ! global%verbLevel

          DO ics = 1,pGrid%nCellsSpecial
            READ(iFile,*) pGrid%cellsSpecial(ics)
          END DO ! ics 
          
          READ(iFile,*) pGrid%nFacesSpecial

          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
                                           'Number of special faces:', &
                                           pGrid%nFacesSpecial
          END IF ! global%verbLevel

          DO ifs = 1,pGrid%nFacesSpecial
            READ(iFile,*) pGrid%facesSpecial(1,ifs), & 
                          pGrid%facesSpecial(2,ifs)
          END DO ! ifs                                        
        ELSE 
          pGrid%nCellsSpecial = 0
          pGrid%nFacesSpecial = 0          
        END IF ! dummyLogical
      ELSE 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! readMode
      
! ==============================================================================
!   Found entry for other region(s) - skip by reading dummy variables
! ==============================================================================

    ELSE 
      READ(iFile,'(L1)') dummyLogical
    
      IF ( dummyLogical .EQV. .TRUE. ) THEN 
        READ(iFile,*) dummyInteger ! Special cells

        DO ics = 1,dummyInteger
          READ(iFile,*) dummyInteger2
        END DO ! ics 
        
        READ(iFile,*) dummyInteger ! Special faces

        DO ifs = 1,dummyInteger
          READ(iFile,*) dummyInteger2
        END DO ! ifs          
      END IF ! dummyLogical    
    END IF ! iRegionGlobal
  END DO ! iReg

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading post-processor info done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ReadPostInfo

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ReadPostInfo.F90,v $
! Revision 1.7  2008/12/06 08:45:06  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:18:17  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.4  2004/09/27 01:43:50  haselbac
! Modified to read info about special faces
!
! Revision 1.3  2004/03/23 03:17:34  haselbac
! Changed format statements
!
! Revision 1.2  2003/08/07 15:37:16  haselbac
! Changed var names
!
! Revision 1.1  2003/06/04 22:44:14  haselbac
! Initial revision
!
! ******************************************************************************







