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
! Purpose: Check for posivity of scalar variables
!
! Description: None.
!
! Input: 
!   pRegion     Region data
!   moduleType  Type of module
!   nVarScal    Number of scalar variables
!   cvScal      Conserved scalar variables
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_ScalarCheckPositivity.F90,v 1.4 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ScalarCheckPositivity(pRegion,moduleType,nVarScal,cvScal)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid  
  USE ModParameters
  USE ModMPI
    
  USE ModInterfaces, ONLY: RFLU_PrintLocInfo  
    
  IMPLICIT NONE

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: moduleType,nVarScal
  REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cvScal
  TYPE(t_region), POINTER :: pRegion

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER, PARAMETER :: MAX_NEGATIVE_LOCS = 10
  INTEGER :: icg,iVarScal,nLocs
  INTEGER :: loc(MAX_NEGATIVE_LOCS,MIN_VAL:MAX_VAL)
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ScalarCheckPositivity.F90,v $ $Revision: 1.4 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ScalarCheckPositivity',&
  'RFLU_ScalarCheckPositivity.F90')

! *****************************************************************************
! Set pointers and variables
! *****************************************************************************

  pGrid => pRegion%grid

  nLocs = 0

! *****************************************************************************
! Loop over cells and check for positivity
! *****************************************************************************

  cellLoop: DO icg = 1,pGrid%nCells        
    varLoop: DO iVarScal = 1,nVarScal
      IF ( cvScal(iVarScal,icg) < 0.0_RFREAL ) THEN
        nLocs = nLocs + 1

        IF ( nLocs == 1 ) THEN 
          WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                'Negative positive-definite variables detected!'
                
          SELECT CASE ( moduleType ) 
            CASE ( FTYPE_SPEC ) 
              WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Module: Species.'
            CASE DEFAULT 
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! moduleType                    

          IF ( global%flowType == FLOW_UNSTEADY ) THEN 
            WRITE(STDOUT,'(A,3X,A,1X,1PE12.5)') SOLVER_NAME,'Current time:', &
                                                global%currentTime              
          ELSE 
            WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME, &
                                           'Current iteration number:', &
                                           global%currentIter           
          END IF ! global%flowType                 

          WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                           pRegion%iRegionGlobal 
!          WRITE(STDOUT,'(A,6X,A,6(1X,A))') SOLVER_NAME,'#', &
!                                           '   Density   ', &
!                                           '  x-velocity ', &
!                                           '  y-velocity ', &
!                                           '  z-velocity ', &
!                                           '   Pressure  ', &
!                                           ' Temperature '       
        END IF ! nLocs

        IF ( nLocs <= MAX_NEGATIVE_LOCS ) THEN 
          WRITE(STDOUT,'(A,4X,I3,6(1X,E13.6))') SOLVER_NAME,nLocs, & 
                                                cvScal(1:nVarScal,icg)

          loc(nLocs,MIN_VAL:MAX_VAL) = icg
          
          EXIT varLoop ! NOTE this EXIT statement                                   
        END IF ! nLocs
      END IF ! cvScal
    END DO varLoop    
  END DO cellLoop

! *****************************************************************************
! Write out message and call error handling routine
! *****************************************************************************

  IF ( nLocs > 0 ) THEN 
    IF ( nLocs > MAX_NEGATIVE_LOCS ) THEN 
       WRITE(STDOUT,'(A,3X,A,1X,I3,1X,A,1X,I9,1X,A)') SOLVER_NAME, & 
             'Only wrote the first',MAX_NEGATIVE_LOCS,'of',nLocs, & 
             'cells with negative positive-definite variables.'    
      CALL RFLU_PrintLocInfo(pRegion,loc,MAX_NEGATIVE_LOCS, & 
                             LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
    ELSE 
      CALL RFLU_PrintLocInfo(pRegion,loc(1:nLocs,MIN_VAL:MAX_VAL),nLocs, & 
                             LOCINFO_MODE_SILENT,OUTPUT_MODE_ANYBODY)
    END IF ! nLocs
    
    CALL ErrorStop(global,ERR_NEGATIVE_POSDEF,__LINE__)   
  END IF ! nLocs

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ScalarCheckPositivity

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ScalarCheckPositivity.F90,v $
! Revision 1.4  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.1  2004/01/29 22:56:05  haselbac
! Initial revision
!
!******************************************************************************







