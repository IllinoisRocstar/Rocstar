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
! Purpose: Collection of routines for particle statistics on Eulerian grid.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: PLAG_ModStats.F90,v 1.4 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_ModStats

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY: t_plag 
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI
  USE PLAG_ModParameters

#ifdef RFLO
#include "Indexing.h"  
  USE ModInterfaces, ONLY: RFLO_GetCellOffSet, &
                           RFLO_GetDimensDummy
#endif
  
  USE PLAG_ModEulerian, ONLY: PLAG_CalcEulerianField
 
  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: PLAG_ModStats.F90,v $ $Revision: 1.4 $'        

! ==============================================================================
! Public functions
! ==============================================================================
  
  PUBLIC :: PLAG_CreateStat,  & 
            PLAG_DestroyStat, & 
            PLAG_InitStat

! ==============================================================================
! Private functions
! ==============================================================================

!  PRIVATE :: 

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
  
  
  
  
  
  
  

  





! *******************************************************************************
!
! Purpose: Create statistics.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPlag               Pointer to plag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_CreateStat(pRegion,pPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
#ifdef RFLO
    TYPE(t_region) :: pRegion 
#endif

#ifdef RFLU
    TYPE(t_region), POINTER :: pRegion
#endif

    TYPE(t_plag), POINTER :: pPlag

! ==============================================================================
!   Locals
! ==============================================================================

#ifdef RFLO
    INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend  
    INTEGER :: iLev,iCellOffset, ijCellOffset   
#endif

    INTEGER :: errorFlag,ibc,iec,nTav
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_CreateStat',&
  'PLAG_ModStats.F90')

    nTav =global%plagNStat     

#ifdef RFLO
    iLev  = pRegion%currLevel
    pGrid => pRegion%levels(iLev)%grid

    CALL RFLO_GetDimensDummy( pRegion,iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( pRegion,iLev,iCellOffset,ijCellOffset )
    ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCellOffset,ijCellOffset)
    iec = IndIJK(idcend,jdcend,kdcend,iCellOffset,ijCellOffset)
#endif

#ifdef RFLU
    pGrid => pRegion%grid
    ibc = 1
    iec = pGrid%nCells
#endif
  
! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    ALLOCATE(pPlag%tav(nTav,ibc:iec),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%tav')
    END IF ! global%error 

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_CreateStat 






! *******************************************************************************
!
! Purpose: Destroy statistics.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPlag               Pointer to plag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_DestroyStat(pRegion,pPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
#ifdef RFLO
    TYPE(t_region) :: pRegion 
#endif

#ifdef RFLU
    TYPE(t_region), POINTER :: pRegion
#endif

    TYPE(t_plag), POINTER :: pPlag
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag
    TYPE(t_global), POINTER :: global    

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_DestroyStat',&
  'PLAG_ModStats.F90')

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    DEALLOCATE(pPlag%tav,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPlag%tav')
    END IF ! global%error
  
! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    NULLIFY(pPlag%tav)
 
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_DestroyStat  





! *******************************************************************************
!
! Purpose: Initialize statistics.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPlag               Pointer to plag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_InitStat(pRegion,pPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
     
#ifdef RFLO
    TYPE(t_region) :: pRegion
#endif

#ifdef RFLU
    TYPE(t_region), POINTER :: pRegion
#endif

    TYPE(t_plag), POINTER :: pPlag
  
! ==============================================================================
!   Locals
! ==============================================================================

#ifdef RFLO
    INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend  
    INTEGER :: iLev,iCellOffset, ijCellOffset   
#endif

    INTEGER :: errorFlag,ibc,iec,iCell,iVar,nVars
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'PLAG_InitStat',&
  'PLAG_ModStats.F90')

    nVars = global%plagNStat      

#ifdef RFLO
    iLev  = pRegion%currLevel
    pGrid => pRegion%levels(iLev)%grid

    CALL RFLO_GetDimensDummy( pRegion,iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( pRegion,iLev,iCellOffset,ijCellOffset )
    ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCellOffset,ijCellOffset)
    iec = IndIJK(idcend,jdcend,kdcend,iCellOffset,ijCellOffset)
#endif

#ifdef RFLU
    pGrid => pRegion%grid
    ibc = 1
    iec = pGrid%nCells
#endif

! ******************************************************************************
!   Initialize memory
! ******************************************************************************
        
    DO iCell = ibc,iec
      DO iVar = 1, nVars
        pPlag%tav(iVar,iCell) = 0.0_RFREAL
      END DO ! iVar                  
    END DO ! iCell                      

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_InitStat








END MODULE PLAG_ModStats

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModStats.F90,v $
! Revision 1.4  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.1  2005/01/08 20:44:32  fnajjar
! Initial import for PLAG statistics
!
! ******************************************************************************









