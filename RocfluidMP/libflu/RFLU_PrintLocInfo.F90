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
! Purpose: Print out locations of cells (in terms of coordinates and whether
!   they are interior or boundary cells)
!
! Description: None.
!
! Input: 
!   pRegion             Pointer to region
!   locUnsorted         Unsorted location array
!   nLocUnsorted        First dimension of locUnsorted array
!   locInfoMode         Indicating verbosity of output 
!   outputMode          Indicating who may write to screen
! 
! Output: N/A.
!
! Notes: 
!   1. Yes yes yes, GOTO is not very nice, but it does the job here...
!   2. Introduced additional argument so that do not get a lot of information 
!      on cell connectivity and vertex coordinates every time this routine 
!      is called. The additional information is mainly interesting if running 
!      with grid motion.
!   3. The additional information on cell connectivity and vertex coordinates
!      can be useful because a cell may not a boundary cell, but still contain
!      vertices which are on boundary...
!
! ******************************************************************************
!
! $Id: RFLU_PrintLocInfo.F90,v 1.18 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2000-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PrintLocInfo(pRegion,locUnsorted,nLocUnsorted,locInfoMode, & 
                             outputMode)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModMPI  

  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModSortSearch

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
! Parameters
! ==============================================================================   
  
  INTEGER, INTENT(IN) :: locInfoMode,nLocUnsorted,outputMode
  INTEGER, INTENT(INOUT) :: locUnsorted(1:nLocUnsorted,MIN_VAL:MAX_VAL)
  TYPE(t_region), POINTER :: pRegion
  
! ==============================================================================
! Local variables
! ==============================================================================

  LOGICAL :: outputFlag
  CHARACTER(CHRLEN) :: cellTypeString,locString,RCSIdentString,tempString
  INTEGER :: errorFlag,il,ib,ic,icCntr,icFlag,icl,ict,id,iPatch,iv,ivCntr, &
             ivFlag,ivg,jl,nLocSorted,nLocSortedEst,nVertSimplified, &
             nVertUnsorted,vLen
  INTEGER :: v(8)
  INTEGER, DIMENSION(:), ALLOCATABLE :: bf2cSorted,locSorted,locBound, & 
                                        vertSorted,vertUnsorted
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

  RCSIdentString = '$RCSfile: RFLU_PrintLocInfo.F90,v $ $Revision: 1.18 $'

! ******************************************************************************
! Start, set output flag
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_PrintLocInfo',&
  'RFLU_PrintLocInfo.F90')

  IF ( outputMode == OUTPUT_MODE_MASTER_ONLY ) THEN 
    outputFlag = (global%myProcid == MASTERPROC)
  ELSE IF ( outputMode == OUTPUT_MODE_ANYBODY ) THEN 
    outputFlag = .TRUE.
  ELSE 
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END IF ! outputMode

  IF ( (outputFlag .EQV. .TRUE.) .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN        
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing location information...'
    WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                     pRegion%iRegionGlobal 
  END IF ! global

! ==============================================================================
! Set grid pointer
! ==============================================================================

  pGrid => pRegion%grid

! ==============================================================================
! Check whether cell centroids allocated and computed
! ==============================================================================

  IF ( ASSOCIATED(pGrid%cofg) .EQV. .FALSE. ) THEN
    global%warnCounter = global%warnCounter + 1  
  
    IF ( (outputFlag .EQV. .TRUE.) .AND. &
         global%verbLevel >= VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,A,1X,A)') SOLVER_NAME, & 
                   '*** WARNING *** Array cofg not allocated.', &
                   'Returning to calling procedure.'
    END IF ! global

    GOTO 1
  END IF ! ASSOCIATED
  
  IF ( (MINVAL(pGrid%cofg) == 0.0_RFREAL) .AND. & 
       (MAXVAL(pGrid%cofg) == 0.0_RFREAL) ) THEN
    global%warnCounter = global%warnCounter + 1       
        
    IF ( (outputFlag .EQV. .TRUE.) .AND. &
         global%verbLevel >= VERBOSE_NONE ) THEN        
      WRITE(STDOUT,'(A,3X,A,1X,A,1X,A)') SOLVER_NAME,'*** WARNING ***', & 
                   'Geometry apparently not computed yet.', & 
                   'Returning to calling procedure.'
    END IF ! global 
    
    GOTO 1         
  END IF ! MINVAL

! ==============================================================================
! Check whether boundary face data structure allocated and meaningful
! ==============================================================================
  
  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)      
    
    IF ( ASSOCIATED(pPatch%bf2c) .EQV. .FALSE. ) THEN
      global%warnCounter = global%warnCounter + 1    
    
      IF ( (outputFlag .EQV. .TRUE.) .AND. &
           global%verbLevel >= VERBOSE_NONE ) THEN        
        WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME, &
                     '*** WARNING *** Array bf2c not allocated.', & 
                     'Returning to calling procedure.'
      END IF ! global
      
      GOTO 1      
    END IF ! ASSOCIATED
    
    IF ( (MINVAL(pPatch%bf2c) == 0) .AND. (MAXVAL(pPatch%bf2c) == 0) ) THEN
      global%warnCounter = global%warnCounter + 1    
    
      IF ( (outputFlag .EQV. .TRUE.) .AND. &
           global%verbLevel >= VERBOSE_NONE ) THEN        
        WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME, &
                     '*** WARNING *** Array bf2c apparently not yet filled.', & 
                     'Returning to calling procedure.'
      END IF ! global
      
      GOTO 1   
    END IF ! ASSOCIATED    
  END DO ! iPatch

! ******************************************************************************
! Sort and merge array of extrema locations
! ******************************************************************************

  CALL QuickSortInteger(locUnsorted(1:nLocUnsorted,MIN_VAL),nLocUnsorted)
  CALL QuickSortInteger(locUnsorted(1:nLocUnsorted,MAX_VAL),nLocUnsorted)

  nLocSortedEst = 2*nLocUnsorted

  ALLOCATE(locSorted(nLocSortedEst),STAT=errorFlag)
  global%error = errorFlag    
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'locSorted')
  END IF ! global%error
  
  locSorted(:) = 0
     
  CALL MergeSortedIntegers(global,nLocUnsorted,nLocUnsorted, & 
                           locUnsorted(1:nLocUnsorted,MIN_VAL), & 
                           locUnsorted(1:nLocUnsorted,MAX_VAL), & 
                           nLocSortedEst,nLocSorted,locSorted)
                                                 
  IF ( (outputFlag .EQV. .TRUE.) .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Cell location information:'                                 
    WRITE(STDOUT,'(A,6X,A,6X,A,3X,A,2(4X,A),3X,A)') SOLVER_NAME,'#','Cell', & 
                 'x-coordinate','y-coordinate','z-coordinate','Location'
  END IF ! global
  
! ******************************************************************************
! Loop over locations, determine more information on cells
! ******************************************************************************

  ALLOCATE(locBound(pGrid%nPatches),STAT=errorFlag)
  global%error = errorFlag    
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'locBound')
  END IF ! global%error

  DO il = 1,nLocSorted
    ic = locSorted(il)
    
    icCntr = 0
    
! ==============================================================================
!   Loop over boundary face-to-cell lists and search for cell. NOTE bf2c list 
!   is not sorted, so need to sort first for binary search to work.
! ==============================================================================    
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
            
      IF ( pPatch%nBFaces > 0 ) THEN 
        ALLOCATE(bf2cSorted(pPatch%nBFaces),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'bf2cSorted')
        END IF ! global%error

        bf2cSorted(1:pPatch%nBFaces) = pPatch%bf2c(1:pPatch%nBFaces)

        CALL QuickSortInteger(bf2cSorted,pPatch%nBFaces)      
        CALL BinarySearchInteger(bf2cSorted,pPatch%nBFaces,ic,icFlag)

        IF ( icFlag /= ELEMENT_NOT_FOUND ) THEN ! vertex found
          icCntr = icCntr + 1 
          locBound(icCntr) = pPatch%iPatchGlobal        
        END IF ! ivFlag
        
        DEALLOCATE(bf2cSorted,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'bf2cSorted')
        END IF ! global%error        
      END IF ! pPatch%nBFaces
    END DO ! iPatch  
          
! ==============================================================================
!   Cell found, write out information
! ==============================================================================          
          
    IF ( icCntr == 0 ) THEN ! interior cell
      WRITE(locString,'(A)') 'Interior'
    ELSE ! cell adjacent to patch
      IF ( icCntr == 1 ) THEN 
        WRITE(locString,'(A,1X,I3)') 'Global patch:',locBound(icCntr)
      ELSE
        WRITE(locString,'(A)') 'Global patches: '
        DO jl = 1,icCntr
          WRITE(tempString,'(1X,I3)') locBound(jl)
          locString = TRIM(locString)//TRIM(tempString)
        END DO ! jl
      END IF ! icCntr
    END IF ! icCntr       
          
    IF ( (outputFlag .EQV. .TRUE.) .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN       
      WRITE(STDOUT,'(A,4X,I3,1X,I9,3(1X,E15.8),2X,A)') SOLVER_NAME,il,ic, & 
                   pGrid%cofg(XCOORD:ZCOORD,ic),TRIM(locString)
    END IF ! global
  END DO ! il

! ******************************************************************************
! If in verbose mode, print more information
! ******************************************************************************

  IF ( locInfoMode == LOCINFO_MODE_VERBOSE ) THEN 

! ==============================================================================
!   Determine more information on vertices which make up cells
! ==============================================================================

    IF ( (outputFlag .EQV. .TRUE.) .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A,2(1X,A))') SOLVER_NAME,'Cell connectivity', &
                                       'information:'
      WRITE(STDOUT,'(A,6X,A,6X,A,2X,A,24X,A)') SOLVER_NAME,'#','Cell', & 
                                               'Type','Vertices'                                 
    END IF ! global

    ALLOCATE(vertUnsorted(8*nLocSorted),STAT=errorFlag)
    global%error = errorFlag    
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vertUnsorted')
    END IF ! global%error

    nVertUnsorted = 0

    DO il = 1,nLocSorted
      ic = locSorted(il)

      ict = pGrid%cellGlob2Loc(1,ic)
      icl = pGrid%cellGlob2Loc(2,ic) 

      SELECT CASE ( ict ) 
        CASE ( CELL_TYPE_TET ) 
          cellTypeString = 'Tetrahedron'
          vLen = 4
          v(1:vLen) = pGrid%tet2v(1:vLen,icl)

          IF ( (outputFlag .EQV. .TRUE.) .AND. &
               global%verbLevel >= VERBOSE_HIGH ) THEN       
            WRITE(STDOUT,'(A,4X,I3,1X,I9,2X,A11,4(1X,I9))') SOLVER_NAME,il, & 
                  ic,cellTypeString,v(1:vLen)
          END IF ! global  

          vertUnsorted(nVertUnsorted+1:nVertUnsorted + vLen) = v(1:vLen)
          nVertUnsorted = nVertUnsorted + vLen                                
        CASE ( CELL_TYPE_HEX )
          cellTypeString = 'Hexahedron' 
          vLen = 8 
          v(1:vLen) = pGrid%hex2v(1:vLen,icl) 

          IF ( (outputFlag .EQV. .TRUE.) .AND. &
               global%verbLevel >= VERBOSE_HIGH ) THEN       
            WRITE(STDOUT,'(A,4X,I3,1X,I9,2X,A11,8(1X,I9))') SOLVER_NAME,il, & 
                  ic,cellTypeString,v(1:vLen)
          END IF ! global    

          vertUnsorted(nVertUnsorted+1:nVertUnsorted + vLen) = v(1:vLen)
          nVertUnsorted = nVertUnsorted + vLen                          
        CASE ( CELL_TYPE_PRI ) 
          cellTypeString = 'Prism'

          vLen = 6 
          v(1:vLen) = pGrid%pri2v(1:vLen,icl) 

          IF ( (outputFlag .EQV. .TRUE.) .AND. &
               global%verbLevel >= VERBOSE_HIGH ) THEN       
            WRITE(STDOUT,'(A,4X,I3,1X,I9,2X,A11,6(1X,I9))') SOLVER_NAME,il, & 
                  ic,cellTypeString,v(1:vLen)
          END IF ! global        

          vertUnsorted(nVertUnsorted+1:nVertUnsorted + vLen) = v(1:vLen)
          nVertUnsorted = nVertUnsorted + vLen                              
        CASE ( CELL_TYPE_PYR ) 
          cellTypeString = 'Pyramid'          

          vLen = 5
          v(1:vLen) = pGrid%pyr2v(1:vLen,icl) 

          IF ( (outputFlag .EQV. .TRUE.) .AND. &
               global%verbLevel >= VERBOSE_HIGH ) THEN       
            WRITE(STDOUT,'(A,4X,I3,1X,I9,2X,A11,5(1X,I9))') SOLVER_NAME,il, & 
                  ic,cellTypeString,v(1:vLen)
          END IF ! global                

          vertUnsorted(nVertUnsorted+1:nVertUnsorted + vLen) = v(1:vLen)
          nVertUnsorted = nVertUnsorted + vLen          
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! ict
    END DO ! il

! ==============================================================================
!   Printing vertex information
! ==============================================================================

    IF ( (outputFlag .EQV. .TRUE.) .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A,2(1X,A))') SOLVER_NAME,'Vertex location', & 
                                       'information:'       
      WRITE(STDOUT,'(A,6X,A,4X,A,3X,A,2(4X,A),3X,A)') SOLVER_NAME,'#', & 
            'Vertex','x-coordinate','y-coordinate','z-coordinate','Location'
    END IF ! global

    CALL QuickSortInteger(vertUnsorted,nVertUnsorted)
    CALL SimplifySortedIntegers(vertUnsorted,nVertUnsorted,nVertSimplified)

    DO iv = 1,nVertSimplified
      ivg = vertUnsorted(iv)

! ------------------------------------------------------------------------------
!     Loop over boundary vertex lists and search for vertex
! ------------------------------------------------------------------------------    
     
      ivCntr = 0 

      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        IF ( pPatch%nBVert > 0 ) THEN
          CALL BinarySearchInteger(pPatch%bv,pPatch%nBVert,ivg,ivFlag)
        ELSE 
          ivFlag = ELEMENT_NOT_FOUND
        END IF ! pPatch%nBVert

        IF ( ivFlag /= ELEMENT_NOT_FOUND ) THEN ! vertex found
          ivCntr = ivCntr + 1 
          locBound(ivCntr) = iPatch
        END IF ! ivFlag
      END DO ! iPatch  

! ------------------------------------------------------------------------------    
!     Vertex found, write out information
! ------------------------------------------------------------------------------    
          
      IF ( ivCntr == 0 ) THEN ! interior vertex
        WRITE(locString,'(A)') 'Interior'
      ELSE ! boundary vertex
        IF ( ivCntr == 1 ) THEN 
          WRITE(locString,'(A,1X,I3)') 'Boundary:',locBound(ivCntr)
        ELSE
          WRITE(locString,'(A)') 'Boundaries: '
          DO jl = 1,ivCntr
            WRITE(tempString,'(1X,I3)') locBound(jl)
            locString = TRIM(locString)//TRIM(tempString)
          END DO ! jl
        END IF ! ivCntr
      END IF ! ivCntr       

      IF ( (outputFlag .EQV. .TRUE.) .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN       
        WRITE(STDOUT,'(A,4X,I3,1X,I9,3(1X,E15.8),2X,A)') SOLVER_NAME,iv, & 
              ivg,pGrid%xyz(XCOORD:ZCOORD,ivg),TRIM(locString)
      END IF ! global
    END DO ! iv  

    DEALLOCATE(vertUnsorted,STAT=errorFlag)
    global%error = errorFlag    
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vertUnsorted')
    END IF ! global%error
  END IF ! locInfoMode

  DEALLOCATE(locBound,STAT=errorFlag)
  global%error = errorFlag    
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'locBound')
  END IF ! global%error

  DEALLOCATE(locSorted,STAT=errorFlag)
  global%error = errorFlag  
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'locSorted')
  END IF ! global%error  
  
1 CONTINUE    
  
! ******************************************************************************
! End
! ******************************************************************************
  
  IF ( (outputFlag .EQV. .TRUE.) .AND. &
       global%verbLevel >= VERBOSE_HIGH ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Printing location information done.'  
  END IF ! global
  
  CALL DeregisterFunction(global)    
  
END SUBROUTINE RFLU_PrintLocInfo

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PrintLocInfo.F90,v $
! Revision 1.18  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.17  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.16  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.15  2005/03/23 17:08:02  haselbac
! Bug fix: Only call search when have nonzero vertices
!
! Revision 1.14  2004/10/19 19:25:12  haselbac
! Cosmetics only
!                                     
! Revision 1.13  2003/12/04 03:23:57  haselbac                        
! Changed formatting of output                                        
!
! Revision 1.12  2003/07/22 01:57:21  haselbac                        
! Added global%warnCounter and changed ip to iPatch                   
!
! Revision 1.11  2003/06/04 22:42:08  haselbac                        
! Added argument to indicate who can write to screen                  
!
! Revision 1.10  2003/04/28 22:40:26  haselbac                        
! Changed region to pRegion, was actually a bug                       
!
! Revision 1.9  2003/03/15 17:05:58  haselbac                         
! Added argument, 2 bug fixes, increased output                       
!
! Revision 1.8  2003/01/28 15:40:00  haselbac                         
! Cosmetics only                                                      
!
! Revision 1.7  2002/10/27 18:53:41  haselbac                         
! Changed declaration of locUnsorted                                  
!
! Revision 1.6  2002/10/17 19:57:57  haselbac                         
! Cosmetic changes to output                                          
!
! Revision 1.5  2002/10/08 15:48:56  haselbac                         
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem  
!
! Revision 1.4  2002/10/05 18:48:20  haselbac                         
! Cosmetic changes to output formatting                               
!
! Revision 1.3  2002/09/09 14:15:01  haselbac                         
! global now under regions, changed interface to MergeSortedIntegers  
!
! Revision 1.2  2002/06/17 13:31:22  haselbac                         
! Prefixed SOLVER_NAME to all screen output                           
!
! Revision 1.1  2002/05/04 16:09:00  haselbac                         
! Initial revision                                                    
!
! ******************************************************************************







