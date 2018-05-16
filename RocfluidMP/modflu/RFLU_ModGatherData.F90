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
! Purpose: Gather and scatter solution data.
!
! Description: None
!
! Notes:
!   1. This module collects routines to gather and scatter solution data
!      when partitioning or merging solution data for parallel multi-physics
!      calculations. The way the partitioning and merging is done means that
!      a single solution vector must be used. Since the solution data for
!      multi-physics calculations in RocfluMP is stored separately, it needs
!      to be gathered before partitioning, and scattered after merging.
!
! ******************************************************************************
!
! $Id: RFLU_ModGatherData.F90,v 1.8 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModGatherData

  USE ModDataTypes
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModError

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_CountGatheredData, &
            RFLU_CreateGatheredData, &
            RFLU_DestroyGatheredData, &             
            RFLU_GatherData, &
            RFLU_ScatterGatheredData

  INTEGER, PARAMETER, PUBLIC :: GATHER_MODE_ACTUAL_ONLY    = 1, &
                                GATHER_MODE_ACTUAL_VIRTUAL = 2

  CONTAINS


! ******************************************************************************
!
! Purpose: Count gathered data - determine number of variables
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: 
!   nVarsOut    Number of variables to be gathered
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_CountGatheredData(pRegion,nVarsOut)

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(OUT) :: nVarsOut
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: nVars
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CountGatheredData',&
  'RFLU_ModGatherData.F90')

! ******************************************************************************
!   Determine number of variables
! ******************************************************************************

! ==============================================================================
!   Mixture
! ==============================================================================

    nVars = pRegion%mixtInput%nCv

! ==============================================================================
!   Physical modules
! ==============================================================================

#ifdef SPEC
    IF ( global%specUsed .EQV. .TRUE. ) THEN
      nVars = nVars + pRegion%specInput%nSpecies
    END IF ! global%specUsed
#endif

! ******************************************************************************
!   Set nVarsOut
! ******************************************************************************

    nVarsOut = nVars

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_CountGatheredData





! ******************************************************************************
!
! Purpose: Create gathered data
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   icType      Cell type
!   gatherMode  Gather mode
!
! Output: 
!   nVarsOut    Number of variables to be gathered
!   cv          Array to hold gathered data
!
! Notes:
!   1. Gather mode indicates whether actual cells or virtual cells are to 
!      be gathered.
!
! ******************************************************************************

  SUBROUTINE RFLU_CreateGatheredData(pRegion,icType,gatherMode,nVarsOut,cv)

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: gatherMode,icType
    INTEGER, INTENT(OUT) :: nVarsOut
    REAL(RFREAL), DIMENSION(:,:), POINTER :: cv
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,nCells,nVars
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateGatheredData',&
  'RFLU_ModGatherData.F90')

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Determine number of variables
! ******************************************************************************

    CALL RFLU_CountGatheredData(pRegion,nVars)

! ******************************************************************************
!   Determine number of cells
! ******************************************************************************

    SELECT CASE ( gatherMode )
      CASE ( GATHER_MODE_ACTUAL_ONLY )
        SELECT CASE ( icType )
          CASE ( CELL_TYPE_TET )
            nCells = pGrid%nTets
          CASE ( CELL_TYPE_HEX )
            nCells = pGrid%nHexs
          CASE ( CELL_TYPE_PRI )
            nCells = pGrid%nPris
          CASE ( CELL_TYPE_PYR )
            nCells = pGrid%nPyrs
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! icType
      CASE ( GATHER_MODE_ACTUAL_VIRTUAL )
        SELECT CASE ( icType )
          CASE ( CELL_TYPE_TET )
            nCells = pGrid%nTetsTot
          CASE ( CELL_TYPE_HEX )
            nCells = pGrid%nHexsTot
          CASE ( CELL_TYPE_PRI )
            nCells = pGrid%nPrisTot
          CASE ( CELL_TYPE_PYR )
            nCells = pGrid%nPyrsTot
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! icType
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! gatherMode

! ******************************************************************************
!   Set nVarsOut
! ******************************************************************************

    nVarsOut = nVars

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    ALLOCATE(cv(nVars,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cv')
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_CreateGatheredData







! ******************************************************************************
!
! Purpose: Destroy gathered data
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: 
!   cv          Array to be destroyed
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_DestroyGatheredData(pRegion,cv)

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), DIMENSION(:,:), POINTER :: cv
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyGatheredData',&
  'RFLU_ModGatherData.F90')

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    DEALLOCATE(cv,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cv')
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_DestroyGatheredData






! ******************************************************************************
!
! Purpose: Gather data
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   icType      Cell type
!   gatherMode  Gather mode
!
! Output: 
!   cv          Gathered data
!
! Notes:
!   1. Gather mode indicates whether actual cells or virtual cells are to 
!      be gathered.
!
! ******************************************************************************

  SUBROUTINE RFLU_GatherData(pRegion,icType,gatherMode,cv)

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: gatherMode,icType
    REAL(RFREAL), DIMENSION(:,:), POINTER :: cv
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg,icl,iVar,nCells,nVars
    INTEGER, DIMENSION(:), POINTER :: pXyz2CellGlob
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pVar
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_GatherData',&
  'RFLU_ModGatherData.F90')

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Determine number of cells and set renumbering pointer
! ******************************************************************************

    SELECT CASE ( gatherMode )
      CASE ( GATHER_MODE_ACTUAL_ONLY )
        SELECT CASE ( icType )
          CASE ( CELL_TYPE_TET )
            nCells = pGrid%nTets
            pXyz2CellGlob => pGrid%tet2CellGlob
          CASE ( CELL_TYPE_HEX )
            nCells = pGrid%nHexs
            pXyz2CellGlob => pGrid%hex2CellGlob
          CASE ( CELL_TYPE_PRI )
            nCells = pGrid%nPris
            pXyz2CellGlob => pGrid%pri2CellGlob
          CASE ( CELL_TYPE_PYR )
            nCells = pGrid%nPyrs
            pXyz2CellGlob => pGrid%pyr2CellGlob
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! icType
      CASE ( GATHER_MODE_ACTUAL_VIRTUAL )
        SELECT CASE ( icType )
          CASE ( CELL_TYPE_TET )
            nCells = pGrid%nTetsTot
            pXyz2CellGlob => pGrid%tet2CellGlob
          CASE ( CELL_TYPE_HEX )
            nCells = pGrid%nHexsTot
            pXyz2CellGlob => pGrid%hex2CellGlob
          CASE ( CELL_TYPE_PRI )
            nCells = pGrid%nPrisTot
            pXyz2CellGlob => pGrid%pri2CellGlob
          CASE ( CELL_TYPE_PYR )
            nCells = pGrid%nPyrsTot
            pXyz2CellGlob => pGrid%pyr2CellGlob
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! icType
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! gatherMode

! ******************************************************************************
!   Gather data
! ******************************************************************************

    DO icl = 1,nCells
      icg = pXyz2CellGlob(icl)

      nVars = 0

! ==============================================================================
!     Mixture
! ==============================================================================

      pVar => pRegion%mixt%cv

      DO iVar = CV_MIXT_DENS,CV_MIXT_ENER
        nVars = nVars + 1
        cv(nVars,icl) = pVar(iVar,icg)
      END DO ! iVar

! ==============================================================================
!     Physical modules
! ==============================================================================

#ifdef SPEC
      IF ( global%specUsed .EQV. .TRUE. ) THEN
        pVar => pRegion%spec%cv

        DO iVar = 1,pRegion%specInput%nSpecies
          nVars = nVars + 1
          cv(nVars,icl) = pVar(iVar,icg)
        END DO ! iVar
      END IF ! global%specUsed
#endif
    END DO ! icl

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_GatherData






! ******************************************************************************
!
! Purpose: Scatter gathered data
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   icType      Cell type
!   gatherMode  Gather mode
!
! Output: 
!   cv          Array holding gathered data
!
! Notes:
!   1. Gather mode indicates whether actual cells or virtual cells are to 
!      be scattered.
!
! ******************************************************************************

  SUBROUTINE RFLU_ScatterGatheredData(pRegion,icType,gatherMode,cv)

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: gatherMode,icType
    REAL(RFREAL), DIMENSION(:,:), POINTER :: cv
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg,icl,iVar,nCells,nVars
    INTEGER, DIMENSION(:), POINTER :: pXyz2CellGlob
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pVar
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ScatterGatheredData',&
  'RFLU_ModGatherData.F90')

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Determine number of cells and set renumbering pointer
! ******************************************************************************

    SELECT CASE ( gatherMode )
      CASE ( GATHER_MODE_ACTUAL_ONLY )
        SELECT CASE ( icType )
          CASE ( CELL_TYPE_TET )
            nCells = pGrid%nTets
            pXyz2CellGlob => pGrid%tet2CellGlob
          CASE ( CELL_TYPE_HEX )
            nCells = pGrid%nHexs
            pXyz2CellGlob => pGrid%hex2CellGlob
          CASE ( CELL_TYPE_PRI )
            nCells = pGrid%nPris
            pXyz2CellGlob => pGrid%pri2CellGlob
          CASE ( CELL_TYPE_PYR )
            nCells = pGrid%nPyrs
            pXyz2CellGlob => pGrid%pyr2CellGlob
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! icType
      CASE ( GATHER_MODE_ACTUAL_VIRTUAL )
        SELECT CASE ( icType )
          CASE ( CELL_TYPE_TET )
            nCells = pGrid%nTetsTot
            pXyz2CellGlob => pGrid%tet2CellGlob
          CASE ( CELL_TYPE_HEX )
            nCells = pGrid%nHexsTot
            pXyz2CellGlob => pGrid%hex2CellGlob
          CASE ( CELL_TYPE_PRI )
            nCells = pGrid%nPrisTot
            pXyz2CellGlob => pGrid%pri2CellGlob
          CASE ( CELL_TYPE_PYR )
            nCells = pGrid%nPyrsTot
            pXyz2CellGlob => pGrid%pyr2CellGlob
          CASE DEFAULT
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END SELECT ! icType
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! gatherMode

! ******************************************************************************
!   Scatter data for actual cells
! ******************************************************************************

    DO icl = 1,nCells
      icg = pXyz2CellGlob(icl)

      nVars = 0

! ==============================================================================
!     Mixture
! ==============================================================================

      pVar => pRegion%mixt%cv

      DO iVar = CV_MIXT_DENS,CV_MIXT_ENER
        nVars = nVars + 1
        pVar(iVar,icg) = cv(nVars,icl)
      END DO ! iVar

! ==============================================================================
!     Physical modules
! ==============================================================================

#ifdef SPEC
      IF ( global%specUsed .EQV. .TRUE. ) THEN
        pVar => pRegion%spec%cv

        DO iVar = 1,pRegion%specInput%nSpecies
          nVars = nVars + 1
          pVar(iVar,icg) = cv(nVars,icl)
        END DO ! iVar
      END IF ! global%specUsed
#endif
    END DO ! icl

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ScatterGatheredData




END MODULE RFLU_ModGatherData

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModGatherData.F90,v $
! Revision 1.8  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.5  2004/11/02 02:31:56  haselbac
! Replaced CV_MIXT_NEQS
!
! Revision 1.4  2004/10/19 19:27:56  haselbac
! Clean-up
!
! Revision 1.3  2004/07/28 15:29:20  jferry
! created global variable for spec use
!
! Revision 1.2  2004/01/22 16:03:59  haselbac
! Made contents of modules PRIVATE, only procs PUBLIC, to avoid errors on ALC 
! and titan
!
! Revision 1.1  2003/11/25 21:03:30  haselbac
! Initial revision
!
! ******************************************************************************











