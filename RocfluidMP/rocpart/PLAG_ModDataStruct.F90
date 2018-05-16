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
! Purpose: Suite of routines for manipulation of particle data structure.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: PLAG_ModDataStruct.F90,v 1.5 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2007 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_ModDataStruct
  
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY: t_plag
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI

  USE PLAG_ModParameters

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: PLAG_DSTR_BuildCell2PclList, &
            PLAG_DSTR_CopyParticleWrapper, & 
            PLAG_DSTR_CreatePclListCSR, &
            PLAG_DSTR_DestroyCell2PclList, &
            PLAG_DSTR_DestroyPclListCSR, &
            PLAG_DSTR_MergeParticleWrapper
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: PLAG_ModDataStruct.F90,v $ $Revision: 1.5 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS







! ******************************************************************************
!
! Purpose: Build cell-to-particle list in CSR format.
!
! Description: None.
!
! Input: 
!   pRegion	Pointer to region
!
! Output: None.
!
! Notes: 
!   1. This routine builds the list of cells which have non-zero particles 
!      and the list of particles in CSR format as well as an access array for
!      the CSR list. 
!   2. In this routine, only memory for the list of cells and the access array
!      is allocated, NOT for the CSR list itself (see next note).
!   3. The actual CSR list is expected to have been created before this 
!      routine is called.
!
! ******************************************************************************

SUBROUTINE PLAG_DSTR_BuildCell2PclList(pRegion)

  USE ModSortSearch

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion
  
! ==============================================================================  
! Locals
! ==============================================================================  

  INTEGER :: errorFlag,icg,iCSR,iLoc,iPcl
  INTEGER, DIMENSION(:), POINTER :: nPclsPerCell
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_DSTR_BuildCell2PclList',&
  'PLAG_ModDataStruct.F90')

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building cell-to-particle list...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                     pRegion%iRegionGlobal
  END IF ! global%verbLevel

! ******************************************************************************
! Set pointers and values
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

  pPlag%nCellsNzPcl    = 0 
  pPlag%nCellsNzPclMax = MIN(1000,pGrid%nCells) ! Initial guess
  
! ******************************************************************************
! Build list of cells with non-zero particles and CSR info array which will be
! used to access actual CSR data structure
! ******************************************************************************

  CALL PLAG_DSTR_CreateCell2PclList(pRegion)
  
  DO iPcl = 1,pPlag%nPcls
    icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)
  
    IF ( pPlag%nCellsNzPcl > 0 ) THEN 
      CALL BinarySearchInteger(pPlag%icgNzPcl(1:pPlag%nCellsNzPcl), & 
                               pPlag%nCellsNzPcl,icg,iLoc)

      IF ( iLoc == ELEMENT_NOT_FOUND ) THEN
        IF ( pPlag%nCellsNzPcl == pPlag%nCellsNzPclMax ) THEN 
          CALL PLAG_DSTR_RecreateCell2PclList(pRegion)
        END IF ! pPlag%nCellsNzPcl
       
        pPlag%nCellsNzPcl = pPlag%nCellsNzPcl + 1
        pPlag%icgNzPcl(pPlag%nCellsNzPcl) = icg
        pPlag%iPclPerCellCSRInfo(pPlag%nCellsNzPcl) = 1
        
        CALL QuickSortIntegerInteger(pPlag%icgNzPcl(1:pPlag%nCellsNzPcl), &
             pPlag%iPclPerCellCSRInfo(1:pPlag%nCellsNzPcl),pPlag%nCellsNzPcl)
      ELSE 
        pPlag%iPclPerCellCSRInfo(iLoc) = pPlag%iPclPerCellCSRInfo(iLoc) + 1
      END IF ! iLoc                         
    ELSE 
      pPlag%nCellsNzPcl = 1
      pPlag%icgNzPcl(1) = icg
      pPlag%iPclPerCellCSRInfo(1) = 1 
    END IF ! pGrid%nCellsNzPcl
  END DO ! iPcl

  IF ( global%verbLevel > VERBOSE_LOW ) THEN 
    WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME,'Number of cells with '// &
                                   'non-zero particles:',pPlag%nCellsNzPcl
  END IF ! global%verbLevel

! ******************************************************************************
! Build actual CSR data structure
! ******************************************************************************

! ==============================================================================
! Finalize CSR info array and check for consistency 
! ==============================================================================

  DO icg = 2,pPlag%nCellsNzPcl
    pPlag%iPclPerCellCSRInfo(icg) = pPlag%iPclPerCellCSRInfo(icg) & 
                                  + pPlag%iPclPerCellCSRInfo(icg-1)
  END DO ! icg
  
  IF ( pPlag%iPclPerCellCSRInfo(pPlag%nCellsNzPcl) /= pPlag%nPcls ) THEN 
    CALL ErrorStop(global,ERR_PLAG_DSTR_INVALID,__LINE__)
  END IF ! pPlag%iPclPerCellCSRInfo

  ALLOCATE(nPclsPerCell(pGrid%nCells),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'nPclsPerCell')
  END IF ! global%error

  DO icg = 1,pGrid%nCells
    nPclsPerCell(icg) = 0
  END DO ! icg

  DO iPcl = 1,pPlag%nPcls
    icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)

    CALL BinarySearchInteger(pPlag%icgNzPcl(1:pPlag%nCellsNzPcl), &
                             pPlag%nCellsNzPcl,icg,iLoc)
    
    IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
      nPclsPerCell(icg) = nPclsPerCell(icg) + 1
      
      IF ( iLoc > 1 ) THEN 
        iCSR = pPlag%iPclPerCellCSRInfo(iLoc-1) + 1
      ELSE 
        iCSR = 1
      END IF ! iLoc

      iCSR = iCSR + nPclsPerCell(icg) - 1
      
      IF ( iCSR > pPlag%nPcls ) THEN ! Defensive coding
        CALL ErrorStop(global,ERR_PLAG_DSTR_INVALID,__LINE__)
      END IF ! iCSR
      
      pPlag%iPclPerCellCSR(iCSR) = iPcl
    ELSE 
      CALL ErrorStop(global,ERR_PLAG_DSTR_INVALID,__LINE__)
    END IF ! iLoc
  END DO ! iPcl

  DEALLOCATE(nPclsPerCell,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'nPclsPerCell')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************
    
  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building cell-to-particle list done.'
  END IF ! global%verbLevel
    
  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_DSTR_BuildCell2PclList







! ******************************************************************************
!
! Purpose: Add particle to data structure.
!
! Description: None.
!
! Input: 
!  global	Pointer to global data
!  pPlag	Pointer to particle data structure (origin)
!  pPlag2	Pointer to particle data structure (target)
!  iPcl 	Index of particle in data structure (origin)
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE PLAG_DSTR_CopyParticle(global,pPlag,pPlag2,iPcl)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  INTEGER, INTENT(IN) ::iPcl
  TYPE(t_global), POINTER :: global
  TYPE(t_plag), POINTER :: pPlag,pPlag2
  
! ==============================================================================  
! Locals
! ==============================================================================  

  INTEGER :: iVar

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'PLAG_DSTR_CopyParticle',&
  'PLAG_ModDataStruct.F90')

! ******************************************************************************
! Start, increment counter, and check dimensions
! ******************************************************************************

  pPlag2%nPcls = pPlag2%nPcls + 1
  
  IF ( pPlag2%nPcls > pPlag2%nPclsMax ) THEN 
    CALL ErrorStop(global,ERR_PLAG_MEMOVERFLOW,__LINE__)
  END IF ! pPlag2%nPcls 

  IF ( (LBOUND(pPlag2%cv,1) /= LBOUND(pPlag%cv,1)) .OR. & 
       (UBOUND(pPlag2%cv,1) /= UBOUND(pPlag%cv,1)) ) THEN 
    CALL ErrorStop(global,ERR_LUBOUND_MISMATCH,__LINE__)
  END IF ! LBOUND 
  
  IF ( (LBOUND(pPlag2%arv,1) /= LBOUND(pPlag%arv,1)) .OR. & 
       (UBOUND(pPlag2%arv,1) /= UBOUND(pPlag%arv,1)) ) THEN 
    CALL ErrorStop(global,ERR_LUBOUND_MISMATCH,__LINE__)
  END IF ! LBOUND
  
  IF ( (LBOUND(pPlag2%aiv,1) /= LBOUND(pPlag%aiv,1)) .OR. & 
       (UBOUND(pPlag2%aiv,1) /= UBOUND(pPlag%aiv,1)) ) THEN 
    CALL ErrorStop(global,ERR_LUBOUND_MISMATCH,__LINE__)
  END IF ! LBOUND    

! ******************************************************************************
! Copy data
! ******************************************************************************

  DO iVar = 1,pPlag2%nCv
    pPlag2%cv(iVar,pPlag2%nPcls) = pPlag%cv(iVar,iPcl)
  END DO ! iVar

  DO iVar = 1,pPlag2%nArv
    pPlag2%arv(iVar,pPlag2%nPcls) = pPlag%arv(iVar,iPcl)
  END DO ! iVar

  DO iVar = 1,pPlag2%nAiv
    pPlag2%aiv(iVar,pPlag2%nPcls) = pPlag%aiv(iVar,iPcl)
  END DO ! iVar
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_DSTR_CopyParticle








! ******************************************************************************
!
! Purpose: Wrapper for adding particles to data structure.
!
! Description: None.
!
! Input: 
!  global	Pointer to global data
!  pPlag	Pointer to particle data structure (origin)
!  pPlag2	Pointer to particle data structure (target)
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE PLAG_DSTR_CopyParticleWrapper(global,pPlag,pPlag2)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_global), POINTER :: global
  TYPE(t_plag), POINTER :: pPlag,pPlag2
  
! ==============================================================================  
! Locals
! ==============================================================================  

  INTEGER :: iPcl

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'PLAG_DSTR_CopyParticleWrapper',&
  'PLAG_ModDataStruct.F90')

! ******************************************************************************
! Start
! ******************************************************************************

  DO iPcl = 1,pPlag%nPcls
    CALL PLAG_DSTR_CopyParticle(global,pPlag,pPlag2,iPcl)
  END DO ! iPcl
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_DSTR_CopyParticleWrapper







! ******************************************************************************
!
! Purpose: Allocate memory for cell-to-particle list.
!
! Description: None.
!
! Input: 
!   pRegion	Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE PLAG_DSTR_CreateCell2PclList(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion
  
! ==============================================================================  
! Locals
! ==============================================================================  

  INTEGER :: errorFlag,icg
  TYPE(t_global), POINTER :: global
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_DSTR_CreateCell2PclList',&
  'PLAG_ModDataStruct.F90')

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  pPlag => pRegion%plag

  ALLOCATE(pPlag%icgNzPcl(pPlag%nCellsNzPclMax),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%icgNzPcl')
  END IF ! global%error
  
  ALLOCATE(pPlag%iPclPerCellCSRInfo(pPlag%nCellsNzPclMax),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%iPclPerCellCSRInfo')
  END IF ! global%error  
  
  DO icg = 1,pPlag%nCellsNzPclMax
    pPlag%iPclPerCellCSRInfo(icg) = 0
  END DO ! icg
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_DSTR_CreateCell2PclList







! ******************************************************************************
!
! Purpose: Allocate memory for cell-to-particle list in CSR format.
!
! Description: None.
!
! Input: 
!   pRegion	Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE PLAG_DSTR_CreatePclListCSR(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion
  
! ==============================================================================  
! Locals
! ==============================================================================  

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_DSTR_CreatePclListCSR',&
  'PLAG_ModDataStruct.F90')

! ******************************************************************************
! Allocate memory
! ******************************************************************************
  
  pPlag => pRegion%plag
  
  IF ( pPlag%nPcls > 0 ) THEN 
    ALLOCATE(pPlag%iPclPerCellCSR(pPlag%nPcls),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%iPclPerCellCSR')
    END IF ! global%error  
  END IF ! pPlag%nPcls

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_DSTR_CreatePclListCSR






! ******************************************************************************
!
! Purpose: Deallocate memory for cell-to-particle list.
!
! Description: None.
!
! Input: 
!   pRegion	Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE PLAG_DSTR_DestroyCell2PclList(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion
  
! ==============================================================================  
! Locals
! ==============================================================================  

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_DSTR_DestroyCell2PclList',&
  'PLAG_ModDataStruct.F90')

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  pPlag => pRegion%plag

  DEALLOCATE(pPlag%icgNzPcl,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPlag%icgNzPcl')
  END IF ! global%error
  
  DEALLOCATE(pPlag%iPclPerCellCSRInfo,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPlag%iPclPerCellCSRInfo')
  END IF ! global%error  
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_DSTR_DestroyCell2PclList






! ******************************************************************************
!
! Purpose: Deallocate memory for cell-to-particle list in CSR format.
!
! Description: None.
!
! Input: 
!   pRegion	Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE PLAG_DSTR_DestroyPclListCSR(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion
  
! ==============================================================================  
! Locals
! ==============================================================================  

  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_DSTR_DestroyPclListCSR',&
  'PLAG_ModDataStruct.F90')

! ******************************************************************************
! Allocate memory
! ******************************************************************************
  
  pPlag => pRegion%plag
  
  IF ( ASSOCIATED(pPlag%iPclPerCellCSR) .EQV. .TRUE. ) THEN 
    DEALLOCATE(pPlag%iPclPerCellCSR,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPlag%iPclPerCellCSR')
    END IF ! global%error  
  END IF ! ASSOCIATED

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_DSTR_DestroyPclListCSR






! ******************************************************************************
!
! Purpose: Wrapper for merging particle data structures.
!
! Description: None.
!
! Input: 
!  global	Pointer to global data
!  pGrid	Pointer to grid data structure
!  pPlag	Pointer to particle data structure (origin)
!  pPlag2	Pointer to particle data structure (target)
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE PLAG_DSTR_MergeParticleWrapper(pRegion,pPlag,pPlag2)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_plag), POINTER :: pPlag,pPlag2
  TYPE(t_region), POINTER :: pRegion
  
! ==============================================================================  
! Locals
! ==============================================================================  

  INTEGER :: icg,icg2,iPcl,nPclsNew,nPclsOld
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_DSTR_MergeParticleWrapper',&
  'PLAG_ModDataStruct.F90')

! ******************************************************************************
! Set pointers
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Copy particles and fix cell and region indices
! ******************************************************************************

  nPclsOld = pPlag2%nPcls

  CALL PLAG_DSTR_CopyParticleWrapper(global,pPlag,pPlag2)

  nPclsNew = pPlag2%nPcls

  DO iPcl = nPclsOld+1,nPclsNew
    icg2 = pPlag2%aiv(AIV_PLAG_ICELLS,iPcl)
    icg  = pGrid%pc2sc(icg2) 
  
    pPlag2%aiv(AIV_PLAG_ICELLS,iPcl) = icg
    pPlag2%aiv(AIV_PLAG_REGINI,iPcl) = pRegion%iRegionGlobal
  END DO ! iPcl
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_DSTR_MergeParticleWrapper







! ******************************************************************************
!
! Purpose: Re-allocate memory for cell-to-particle list.
!
! Description: None.
!
! Input: 
!   pRegion	Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE PLAG_DSTR_RecreateCell2PclList(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion
  
! ==============================================================================  
! Locals
! ==============================================================================  

  INTEGER :: errorFlag,icl
  INTEGER, DIMENSION(:), POINTER :: tempArray1,tempArray2
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_DSTR_RecreateCell2PclList',&
  'PLAG_ModDataStruct.F90')

! ******************************************************************************
! Allocate temporary array and copy into it
! ******************************************************************************

  pGrid => pRegion%grid
  pPlag => pRegion%plag

  ALLOCATE(tempArray1(pPlag%nCellsNzPcl),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'tempArray1')
  END IF ! global%error
  
  ALLOCATE(tempArray2(pPlag%nCellsNzPcl),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'tempArray2')
  END IF ! global%error
  
  DO icl = 1,pPlag%nCellsNzPcl
    tempArray1(icl) = pPlag%icgNzPcl(icl)
    tempArray2(icl) = pPlag%iPclPerCellCSRInfo(icl)
  END DO ! icl

! ******************************************************************************
! Reallocate array, expand size
! ******************************************************************************
  
  CALL PLAG_DSTR_DestroyCell2PclList(pRegion)
  
  pPlag%nCellsNzPclMax = MIN(2*pPlag%nCellsNzPclMax,pGrid%nCells)
  
  CALL PLAG_DSTR_CreateCell2PclList(pRegion)

! ******************************************************************************
! Copy into expanded array and deallocate temporary array
! ******************************************************************************
  
  DO icl = 1,pPlag%nCellsNzPcl
    pPlag%icgNzPcl(icl)           = tempArray1(icl) 
    pPlag%iPclPerCellCSRInfo(icl) = tempArray2(icl) 
  END DO ! icl  
  
  DEALLOCATE(tempArray1,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'tempArray1')
  END IF ! global%error  
  
  DEALLOCATE(tempArray2,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'tempArray2')
  END IF ! global%error  
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_DSTR_RecreateCell2PclList






! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE PLAG_ModDataStruct

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModDataStruct.F90,v $
! Revision 1.5  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/03/27 00:41:47  haselbac
! Bug fix: Incorrect loop limit in copying particle data
!
! Revision 1.2  2007/03/27 00:20:49  haselbac
! Substantial additions to allow faster initialization
!
! Revision 1.1  2007/03/12 23:33:29  haselbac
! Initial revision
!
! ******************************************************************************















