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
! Purpose: Suite of routines for coloring.
!
! Description: None.
!
! Notes: 
!   1. Use of overloading.
!
! ******************************************************************************
!
! $Id: RFLU_ModColoring.F90,v 1.7 2008/12/06 08:44:20 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModColoring

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid  
  USE ModBndPatch, ONLY: t_patch
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_COL_BuildColoring, & 
            RFLU_COL_CreateColoring, &
            RFLU_COL_DestroyColoring, & 
            RFLU_COL_ReadColoring, & 
            RFLU_COL_WriteColoring

! ******************************************************************************
! Interface declaration for overloaded procedures
! ******************************************************************************  
        
  INTERFACE RFLU_COL_BuildColoring
    MODULE PROCEDURE RFLU_COL_BuildColoringP,RFLU_COL_BuildColoringS
  END INTERFACE 
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModColoring.F90,v $ $Revision: 1.7 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  







! ******************************************************************************
!
! Purpose: Build coloring.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pRegionSerial       Pointer to serial region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_COL_BuildColoringP(pRegion,pRegionSerial)

  USE RFLU_ModCopyData, ONLY: RFLU_COPY_CellDataS2P_I1D

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion,pRegionSerial

! ==============================================================================
! Locals
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_COL_BuildColoringP', &
                        'RFLU_ModColoring.F90')

  pGrid => pRegion%grid

! ******************************************************************************
! Copy coloring
! ******************************************************************************

  CALL RFLU_COPY_CellDataS2P_I1D(global,pGrid,pGrid%col,pRegionSerial%grid%col)
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_COL_BuildColoringP








! ******************************************************************************
!
! Purpose: Build coloring.
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

SUBROUTINE RFLU_COL_BuildColoringS(pRegion)

  USE RFLU_ModResidual, ONLY: RFLU_GetResidualSupport1, & 
                              RFLU_GetResidualSupport2

  USE ModSortSearch

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  LOGICAL :: needNewColor
  INTEGER :: cntr,errorFlag,icg,icg2,icl,iLoc,iSoc,nCellMembsMin,nSocMax, & 
             offs,rsSize,rsSizeMax
  INTEGER, DIMENSION(:), ALLOCATABLE :: cellMembsTemp,rs
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid


! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_COL_BuildColoringS', &
                        'RFLU_ModColoring.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building coloring...' 
  END IF ! global%myProcid

  pGrid => pRegion%grid

! ******************************************************************************
! Initialize
! ******************************************************************************

  rsSizeMax     = 10000 ! Maximum allowed size of residual support
  nSocMax       =  2048 ! Maximum allowed number of struct orthog columns
  nCellMembsMin =   100 ! Minimum number of cell members

  pGrid%nSoc    = 0
  pGrid%nSocMax = nSocMax
 
! ******************************************************************************
! Allocate temporary memory
! ******************************************************************************

  ALLOCATE(rs(rsSizeMax),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'rs')
  END IF ! global%errorFlag  

  ALLOCATE(pGrid%soc(pGrid%nSocMax),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%soc')
  END IF ! global%error

  DO iSoc = 1,nSocMax
    ALLOCATE(pGrid%soc(iSoc)%cellMembs(nCellMembsMin),STAT=errorFlag)
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%soc%cellMembs')
    END IF ! global%error
  END DO ! iSoc

! ******************************************************************************
! Loop over cells
! ******************************************************************************

  DO icg = 1,pGrid%nCells

! ==============================================================================
!   Get stencil members and sort
! ==============================================================================

    SELECT CASE ( pRegion%mixtInput%spaceOrder )
      CASE ( 1 ) 
        CALL RFLU_GetResidualSupport1(pRegion,icg,rs,rsSizeMax,rsSize)
      CASE ( 2 ) 
! TEMPORARY
!        CALL RFLU_GetResidualSupport2(pRegion,icg,rs,rsSizeMax,rsSize)
        CALL RFLU_GetResidualSupport1(pRegion,icg,rs,rsSizeMax,rsSize)
! END TEMPORARY
      CASE DEFAULT 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%spaceOrder
                 
! ==============================================================================
!   Loop over colors
! ==============================================================================

    needNewColor = .TRUE.

    colorLoop: DO iSoc = 1,pGrid%nSoc
      
! ------------------------------------------------------------------------------
!     Determine whether stencil members already in this color
! ------------------------------------------------------------------------------      
      
      cntr = 0
      
      needNewColor = .FALSE.
      
      rsLoop: DO icl = 1,rsSize
        icg2 = rs(icl)

        IF ( pGrid%soc(iSoc)%nCellMembs > 0 ) THEN
          CALL BinarySearchInteger(pGrid%soc(iSoc)%cellMembs(1:pGrid%soc(iSoc)%nCellMembs), & 
                                   pGrid%soc(iSoc)%nCellMembs,icg2,iLoc)
        ELSE 
          iLoc = ELEMENT_NOT_FOUND
        END IF ! pGrid%soc%nCellMembs

        IF ( iLoc == ELEMENT_NOT_FOUND ) THEN 
          cntr = cntr + 1
        ELSE
          IF ( iSoc /= pGrid%nSoc ) THEN
            EXIT rsLoop
          ELSE 
            needNewColor = .TRUE.
            
            EXIT colorLoop
          END IF ! iSoc           
        END IF ! iLoc                         
      END DO rsLoop
     
! ------------------------------------------------------------------------------
!     Add to this color if none of the stencil members already present
! ------------------------------------------------------------------------------      
           
      IF ( (needNewColor .EQV. .FALSE.) .AND. (cntr == rsSize) ) THEN
        offs = pGrid%soc(iSoc)%nCellMembs
       
! ----- Array too small, so reallocate -----------------------------------------       
       
        IF ( (pGrid%soc(iSoc)%nCellMembs + rsSize) > SIZE(pGrid%soc(iSoc)%cellMembs,1) ) THEN  
          ALLOCATE(cellMembsTemp(pGrid%soc(iSoc)%nCellMembs),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cellMembsTemp')
          END IF !global%error
       
          DO icl = 1,pGrid%soc(iSoc)%nCellMembs
            cellMembsTemp(icl) = pGrid%soc(iSoc)%cellMembs(icl)
          END DO ! icl
       
          DEALLOCATE(pGrid%soc(iSoc)%cellMembs,STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%soc%cellMembs')
          END IF !global%error          
       
          nCellMembsMin = 2*(pGrid%soc(iSoc)%nCellMembs + rsSize)
       
          ALLOCATE(pGrid%soc(iSoc)%cellMembs(nCellMembsMin),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%soc%cellMembs')
          END IF !global%error        
       
          DO icl = 1,pGrid%soc(iSoc)%nCellMembs
            pGrid%soc(iSoc)%cellMembs(icl) = cellMembsTemp(icl)
          END DO ! icl       
       
          DEALLOCATE(cellMembsTemp,STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cellMembsTemp')
          END IF !global%error       
        END IF ! pGrid%soc%nCellMembs
                
! ----- Add stencil members and sort -------------------------------------------                
        
        pGrid%col(icg) = iSoc
        
        DO icl = 1,rsSize
          pGrid%soc(iSoc)%cellMembs(offs+icl) = rs(icl)
        END DO ! icl
        
        pGrid%soc(iSoc)%nCellMembs = pGrid%soc(iSoc)%nCellMembs + rsSize
        
        CALL QuickSortInteger(pGrid%soc(iSoc)%cellMembs(1:pGrid%soc(iSoc)%nCellMembs), & 
                              pGrid%soc(iSoc)%nCellMembs)     
                              
        EXIT colorLoop                        
      END IF ! needNewColor      
    END DO colorLoop

! ==============================================================================
!   Add new color
! ==============================================================================

    IF ( needNewColor .EQV. .TRUE. ) THEN
      IF ( pGrid%nSoc < pGrid%nSocMax ) THEN  
        pGrid%nSoc = pGrid%nSoc + 1
      ELSE 
! TEMPORARY
        WRITE(*,*) 'ERROR! About to exceed dimensions of soc!'
        STOP
! END TEMPORARY      
      END IF ! nSoc

      pGrid%col(icg) = pGrid%nSoc
      
      IF ( rsSize > SIZE(pGrid%soc(pGrid%nSoc)%cellMembs,1) ) THEN
        DEALLOCATE(pGrid%soc(pGrid%nSoc)%cellMembs,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%soc%cellMembs')
        END IF ! global%error

        ALLOCATE(pGrid%soc(pGrid%nSoc)%cellMembs(2*rsSize),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%soc%cellMembs')
        END IF ! global%error
      END IF ! rsSize

      DO icl = 1,rsSize
        pGrid%soc(pGrid%nSoc)%cellMembs(icl) = rs(icl)
      END DO ! icl   
      
      pGrid%soc(pGrid%nSoc)%nCellMembs = rsSize
      
      CALL QuickSortInteger(pGrid%soc(pGrid%nSoc)%cellMembs(1:pGrid%soc(pGrid%nSoc)%nCellMembs), & 
                            pGrid%soc(pGrid%nSoc)%nCellMembs)         
    END IF ! needNewColor
  END DO ! icg

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME,'Number of colors:',pGrid%nSoc
  END IF ! global%myProcid

! ******************************************************************************
! Deallocate temporary memory
! ******************************************************************************
  
  DEALLOCATE(rs,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'rs')
  END IF ! global%errorFlag
  
  DO iSoc = 1,pGrid%nSocMax
    DEALLOCATE(pGrid%soc(iSoc)%cellMembs,STAT=errorFlag)
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%soc%cellMembs')
    END IF ! global%error
  END DO ! iSoc

  DEALLOCATE(pGrid%soc,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%soc')
  END IF ! global%error

  pGrid%nSoc    = 0
  pGrid%nSocMax = 0  
  
! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building coloring done.' 
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_COL_BuildColoringS






! ******************************************************************************
!
! Purpose: Create coloring.
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

SUBROUTINE RFLU_COL_CreateColoring(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
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
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_COL_CreateColoring', &
                        'RFLU_ModColoring.F90')

  pGrid => pRegion%grid

! ******************************************************************************
! Allocate memory
! ******************************************************************************

  ALLOCATE(pGrid%col(pGrid%nCellsTot),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%col')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_COL_CreateColoring








! ******************************************************************************
!
! Purpose: Destroy coloring.
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

SUBROUTINE RFLU_COL_DestroyColoring(pRegion)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
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
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_COL_DestroyColoring', &
                        'RFLU_ModColoring.F90')

  pGrid => pRegion%grid

! ******************************************************************************
! Destroy memory
! ******************************************************************************

  DEALLOCATE(pGrid%col,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%col')
  END IF ! global%error

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_COL_DestroyColoring








! ******************************************************************************
!
! Purpose: Read coloring.
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

  SUBROUTINE RFLU_COL_ReadColoring(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion  
   
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,iFile,loopCounter,nCellsTot
    CHARACTER(CHRLEN) :: iFileName,sectionString
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_COL_ReadColoring',&
  'RFLU_ModColoring.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading coloring...' 
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN     
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid

! ==============================================================================
!   Open file
! ==============================================================================

    iFile = IF_RNMB

    CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.col', & 
                            pRegion%iRegionGlobal,iFileName) 

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)   
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%myProcid

    READ(iFile,'(A)') sectionString 
    IF ( TRIM(sectionString) /= '# ROCFLU coloring file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM

! ==============================================================================
!   Dimensions
! ==============================================================================
  
    pGrid => pRegion%grid  

    READ(iFile,'(A)') sectionString 
    IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM

    READ(iFile,'(I8)') nCellsTot

! ==============================================================================
!   Check dimensions (against those read from dimensions file)
! ==============================================================================

    IF ( nCellsTot /= pGrid%nCellsTot ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nCellsTot    
    
! ==============================================================================
!   Rest of file
! ==============================================================================

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1
    
      READ(iFile,'(A)') sectionString
    
      SELECT CASE ( TRIM(sectionString) ) 

! ------------------------------------------------------------------------------
!       Vertex renumbering
! ------------------------------------------------------------------------------

        CASE ( '# Coloring' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Coloring...'
          END IF ! global%myProcid 
              
          READ(iFile,'(10(I8))') (pGrid%col(icg),icg=1,pGrid%nCellsTot)
    
! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------ 
      
        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%myProcid  

          EXIT
      
! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------ 
      
        CASE DEFAULT
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN
            WRITE(STDOUT,'(3X,A)') sectionString
          END IF ! global%myProcid      

          CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
      END SELECT ! TRIM
  
! ==============================================================================
!     Guard against infinite loop - might be unnecessary because of read errors?
! ==============================================================================  

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter  
    END DO ! <empty>
   
! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading coloring done.' 
    END IF ! global%myProcid

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_COL_ReadColoring








! ******************************************************************************
!
! Purpose: Recreate cell list.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   nVertPerCell        Number of vertices per cell
!   nCellsMax           Maximum number of cells
!   x2v                 Connectivity array
!   x2cg                Cell mapping array
!
! Output: 
!   nCellsMax           Increased maximum number of cells
!   x2v                 Enlarged connectivity array
!   x2cg                Enlarged cell mapping array
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_COL_RecreateCellList(global,nVertPerCell,nCellsMax,x2v,x2cg)

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: nVertPerCell
    INTEGER, INTENT(INOUT) :: nCellsMax
    INTEGER, DIMENSION(:), POINTER :: x2cg
    INTEGER, DIMENSION(:,:), POINTER :: x2v
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icl,ivl,nCellsMaxOld    
    INTEGER, DIMENSION(:), ALLOCATABLE:: x2cgTemp
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: x2vTemp      

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_COL_RecreateCellList',&
  'RFLU_ModColoring.F90')

! ******************************************************************************
!   Increase maximum number of cells
! ******************************************************************************

    nCellsMaxOld =   nCellsMax 
    nCellsMax    = 2*nCellsMax

! ******************************************************************************
!   Copy existing arrays into larger arrays
! ******************************************************************************

! ==============================================================================
!   Connectivity array
! ==============================================================================

    ALLOCATE(x2vTemp(nVertPerCell,nCellsMaxOld),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'x2vTemp')
    END IF ! global%error
          
    DO icl = 1,nCellsMaxOld
      DO ivl = 1,nVertPerCell
        x2vTemp(ivl,icl) = x2v(ivl,icl)
      END DO ! ivl
    END DO ! icl

    DEALLOCATE(x2v,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'x2v')
    END IF ! global%error

    ALLOCATE(x2v(nVertPerCell,nCellsMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'x2v')
    END IF ! global%error

    DO icl = 1,nCellsMaxOld
      DO ivl = 1,nVertPerCell
        x2v(ivl,icl) = x2vTemp(ivl,icl)
      END DO ! ivl
    END DO ! icl    

    DEALLOCATE(x2vTemp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'x2vTemp')
    END IF ! global%error

! ==============================================================================
!   Cell mapping array
! ==============================================================================

    ALLOCATE(x2cgTemp(nCellsMaxOld),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'x2cgTemp')
    END IF ! global%error
          
    DO icl = 1,nCellsMaxOld
      x2cgTemp(icl) = x2cg(icl)
    END DO ! icl

    DEALLOCATE(x2cg,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'x2cg')
    END IF ! global%error

    ALLOCATE(x2cg(nCellsMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'x2cg')
    END IF ! global%error

    DO icl = 1,nCellsMaxOld
      x2cg(icl) = x2cgTemp(icl)
    END DO ! icl    

    DEALLOCATE(x2cgTemp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'x2cgTemp')
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_COL_RecreateCellList








! ******************************************************************************
!
! Purpose: Write coloring.
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

  SUBROUTINE RFLU_COL_WriteColoring(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion  
   
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,iFile
    CHARACTER(CHRLEN) :: iFileName,sectionString
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_COL_WriteColoring',&
  'RFLU_ModColoring.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing coloring...'
    END IF ! global%myProcid

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal          
    END IF ! global%myProcid 

! ==============================================================================
!   Open file
! ==============================================================================

    iFile = IF_COLOR

    CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.col', & 
                            pRegion%iRegionGlobal,iFileName) 

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", & 
         IOSTAT=errorFlag)   
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%myProcid

    sectionString = '# ROCFLU coloring file'  
    WRITE(iFile,'(A)') TRIM(sectionString)  

! ==============================================================================
!   Dimensions
! ==============================================================================
  
    pGrid => pRegion%grid  

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Dimensions...'
    END IF ! global%myProcid

    sectionString = '# Dimensions'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(I8)') pGrid%nCellsTot

! ==============================================================================
!   Coloring
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Cells...'
    END IF ! global%myProcid

    sectionString = '# Coloring'  
    WRITE(iFile,'(A)') TRIM(sectionString)      
    WRITE(iFile,'(10(I8))') (pGrid%col(icg),icg=1,pGrid%nCellsTot)    
   
! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%myProcid

    sectionString = '# End'
    WRITE(iFile,'(A)') TRIM(sectionString) 

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing coloring done.'
    END IF ! global%myProcid

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_COL_WriteColoring








! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModColoring

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModColoring.F90,v $
! Revision 1.7  2008/12/06 08:44:20  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:31  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/04/07 15:19:18  haselbac
! Removed tabs
!
! Revision 1.4  2005/09/22 17:10:31  hdewey2
! Modified coloring so the Jacobian is always colored based on the 1st order cell stencil.
!
! Revision 1.3  2005/08/24 01:36:15  haselbac
! Fixed bug and extended to second-order
!
! Revision 1.2  2005/08/19 02:35:26  haselbac
! Changed for cColors to col, modified existing and added I/O routines
!
! Revision 1.1  2005/08/17 20:04:52  hdewey2
! Initial revision
!
! ******************************************************************************













