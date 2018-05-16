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
! Purpose: Suite of topology routines.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModTopologyUtils.F90,v 1.9 2008/12/06 08:44:24 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModTopologyUtils

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  USE ModSortSearch

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: RFLU_BuildCellVertList, & 
            RFLU_BuildConnVertList, & 
            RFLU_BuildFaceVertList, & 
            RFLU_BuildVertCellNghbList

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModTopologyUtils.F90,v $ $Revision: 1.9 $'


! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS





! ******************************************************************************
!
! Purpose: Build list of vertices given list of cells.
!
! Description: None.
!
! Input: 
!   global              Pointer to global data
!   pGrid               Pointer to grid 
!   cList               List of cells
!   cListDim            Number of cells
!   vListDimMax         Maximum number of vertices
!
! Output: 
!   vList               List of vertices
!   vListDim            Actual number of vertices
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildCellVertList(global,pGrid,cList,cListDim,vList, & 
                                    vListDimMax,vListDim)

    USE RFLU_ModHashTable

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: cListDim,vListDimMax
    INTEGER, INTENT(OUT) :: vListDim
    INTEGER, INTENT(IN) :: cList(cListDim)
    INTEGER, INTENT(OUT) :: vList(vListDimMax)
    TYPE(t_grid), POINTER :: pGrid 
    TYPE(t_global), POINTER :: global     
      
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: i,icg,icl,ict,ivl,key,errorFlag
      
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_BuildCellVertList',&
  'RFLU_ModTopologyUtils.F90')

! ******************************************************************************
!   Initialize
! ******************************************************************************
           
    vListDim = 0   
    
    DO i = 1,vListDimMax 
      vList(i) = 0
    END DO ! i     

! ******************************************************************************
!   Create hash table
! ******************************************************************************

    CALL RFLU_CreateHashTable(global,vListDimMax) 
    
! ******************************************************************************
!   Build vertex list
! ******************************************************************************
                
    DO i = 1,cListDim
      icg = cList(i)
      
      ict = pGrid%cellGlob2Loc(1,icg)
      icl = pGrid%cellGlob2Loc(2,icg)
      
      SELECT CASE ( ict ) 
        CASE ( CELL_TYPE_TET ) 
          DO ivl = 1,4
            CALL RFLU_HashBuildKey(pGrid%tet2v(ivl:ivl,icl),1,key)          
            CALL RFLU_HashVertex(global,key,pGrid%tet2v(ivl,icl),vListDim, & 
                                 vList)             
          END DO ! ivl
        CASE ( CELL_TYPE_HEX ) 
          DO ivl = 1,8
            CALL RFLU_HashBuildKey(pGrid%hex2v(ivl:ivl,icl),1,key)          
            CALL RFLU_HashVertex(global,key,pGrid%hex2v(ivl,icl),vListDim, & 
                                 vList)           
          END DO ! ivl
        CASE ( CELL_TYPE_PRI ) 
          DO ivl = 1,6
            CALL RFLU_HashBuildKey(pGrid%pri2v(ivl:ivl,icl),1,key)          
            CALL RFLU_HashVertex(global,key,pGrid%pri2v(ivl,icl),vListDim, & 
                                 vList)      
          END DO ! ivl
        CASE ( CELL_TYPE_PYR ) 
          DO ivl = 1,5
            CALL RFLU_HashBuildKey(pGrid%pyr2v(ivl:ivl,icl),1,key)          
            CALL RFLU_HashVertex(global,key,pGrid%pyr2v(ivl,icl),vListDim, & 
                                 vList)       
          END DO ! ivl
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! ict    
    END DO ! i
                                
! ******************************************************************************
!   Destroy hash table
! ******************************************************************************
       
    CALL RFLU_DestroyHashTable(global)         

! ******************************************************************************
!   Sort list of vertices
! ******************************************************************************

    CALL QuickSortInteger(vList(1:vListDim),vListDim)
                                
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_BuildCellVertList








! ******************************************************************************
!
! Purpose: Build list of vertices given list of cells.
!
! Description: None.
!
! Input: 
!   global              Pointer to global data
!   connList            List of entities in connectivity list
!   connListDim1        Number of entries per entity
!   connListDim2        Number of entities in connectivity 
!   vListDimMax         Maximum number of vertices
!
! Output: 
!   vList               List of vertices
!   vListDim            Actual number of vertices
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildConnVertList(global,connList,connListDim1,connListDim2, & 
                                    vList,vListDimMax,vListDim)

    USE RFLU_ModHashTable

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: connListDim1,connListDim2,vListDimMax
    INTEGER, INTENT(OUT) :: vListDim
    INTEGER, INTENT(IN) :: connList(connListDim1,connListDim2)
    INTEGER, INTENT(OUT) :: vList(vListDimMax)
    TYPE(t_global), POINTER :: global     
      
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: i,i1,i2,ivl,key
      
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_BuildConnVertList',&
  'RFLU_ModTopologyUtils.F90')

! ******************************************************************************
!   Initialize
! ******************************************************************************
           
    vListDim = 0   
    
    DO i = 1,vListDimMax 
      vList(i) = 0
    END DO ! i     

! ******************************************************************************
!   Create hash table
! ******************************************************************************

    CALL RFLU_CreateHashTable(global,vListDimMax) 
    
! ******************************************************************************
!   Build vertex list
! ******************************************************************************
                
    DO i2 = 1,connListDim2
      DO i1 = 1,connListDim1
        CALL RFLU_HashBuildKey(connList(i1:i1,i2),1,key)          
        CALL RFLU_HashVertex(global,key,connList(i1,i2),vListDim,vList)             
      END DO ! i1
    END DO ! i2
                                
! ******************************************************************************
!   Destroy hash table
! ******************************************************************************
       
    CALL RFLU_DestroyHashTable(global)         

! ******************************************************************************
!   Sort list of vertices
! ******************************************************************************

    CALL QuickSortInteger(vList(1:vListDim),vListDim)
                                
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_BuildConnVertList









! ******************************************************************************
!
! Purpose: Build list of vertices given list of faces.
!
! Description: None.
!
! Input: 
!   global              Pointer to global data
!   pGrid               Pointer to grid 
!   fList               List of faces
!   fListDim            Number of faces
!   vListDimMax         Maximum number of vertices
!
! Output: 
!   vList               List of vertices
!   vListDim            Actual number of vertices
!   errorFlag           Error flag 
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildFaceVertList(global,pGrid,fList,fListDim,vList, & 
                                    vListDimMax,vListDim,errorFlag)

    USE RFLU_ModHashTable

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: fListDim,vListDimMax
    INTEGER, INTENT(OUT) :: errorFlag,vListDim
    INTEGER, INTENT(IN) :: fList(fListDim)
    INTEGER, INTENT(OUT) :: vList(vListDimMax)
    TYPE(t_grid), POINTER :: pGrid 
    TYPE(t_global), POINTER :: global     
      
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: i,ifg,key
      
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_BuildFaceVertList',&
  'RFLU_ModTopologyUtils.F90')

! ******************************************************************************
!   Initialize
! ******************************************************************************
           
    vListDim = 0   
    
    DO i = 1,vListDimMax 
      vList(i) = 0
    END DO ! i     

! ******************************************************************************
!   Create hash table
! ******************************************************************************

    CALL RFLU_CreateHashTable(global,vListDimMax) 
    
! ******************************************************************************
!   Build vertex list
! ******************************************************************************
                
    iLoop: DO i = 1,fListDim
      ifg = fList(i)
      
      CALL RFLU_HashBuildKey(pGrid%f2v(1,ifg:ifg),1,key)          
      CALL RFLU_HashVertex(global,key,pGrid%f2v(1,ifg),vListDim,vList,errorFlag)       

      IF ( errorFlag /= ERR_NONE ) THEN 
        EXIT iLoop
      END IF ! errorFlag

      CALL RFLU_HashBuildKey(pGrid%f2v(2,ifg:ifg),1,key)          
      CALL RFLU_HashVertex(global,key,pGrid%f2v(2,ifg),vListDim,vList,errorFlag)
      
      IF ( errorFlag /= ERR_NONE ) THEN 
        EXIT iLoop
      END IF ! errorFlag

      CALL RFLU_HashBuildKey(pGrid%f2v(3,ifg:ifg),1,key)          
      CALL RFLU_HashVertex(global,key,pGrid%f2v(3,ifg),vListDim,vList,errorFlag)
            
      IF ( errorFlag /= ERR_NONE ) THEN 
        EXIT iLoop
      END IF ! errorFlag

      IF ( pGrid%f2v(4,ifg) /= VERT_NONE ) THEN 
        CALL RFLU_HashBuildKey(pGrid%f2v(4,ifg:ifg),1,key)          
        CALL RFLU_HashVertex(global,key,pGrid%f2v(4,ifg),vListDim,vList, & 
                             errorFlag)                

        IF ( errorFlag /= ERR_NONE ) THEN 
          EXIT iLoop
        END IF ! errorFlag
      END IF ! pGrid%f2v
    END DO iLoop
                                
! ******************************************************************************
!   Destroy hash table
! ******************************************************************************
       
    CALL RFLU_DestroyHashTable(global)         

! ******************************************************************************
!   Sort list of vertices
! ******************************************************************************

    IF ( errorFlag == ERR_NONE ) THEN
      CALL QuickSortInteger(vList(1:vListDim),vListDim)
    END IF ! errorFlag
                                
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_BuildFaceVertList









! ******************************************************************************
!
! Purpose: Build list of cell neighbors given list of vertices.
!
! Description: None.
!
! Input: 
!   global              Pointer to global data  
!   pGrid               Pointer to grid
!   vListOrig           List of vertices
!   vListOrigDim        Number of vertices in list
!   nLayers             Number of cell layers to be added
!   iReg                Target region (see notes)
!   cList               List of cells (empty)
!   cListDimMax         Dimension of list of cells
!
! Output: 
!   cList               List of cells
!   cListDim            Actual number of cells in list
!
! Notes: 
!    1. Cells are added in a different way depending on the value of iReg:
!       - If iReg = 0 add cell in any case.
!       - If iReg > 0 add cell only if region index of cell after partitioning 
!         is EQUAL to iReg.
!       - If iReg < 0 add cell only if region index of cell after partitioning
!         is NOT EQUAL to iReg.  
!    2. Routine must work even if sc2r mapping not existent, so have special
!       treatment for this case, which can arise if this routine is called 
!       in connection with building sype virtual cells.    
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildVertCellNghbList(global,pGrid,vListOrig,vListOrigDim, &
                                        nLayers,iReg,cList,cListDimMax,cListDim)

    IMPLICIT NONE
      
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: cListDimMax,iReg,nLayers,vListOrigDim
    INTEGER, INTENT(OUT) :: cListDim
    INTEGER, INTENT(IN) :: vListOrig(vListOrigDim)
    INTEGER, INTENT(OUT) :: cList(cListDimMax)
    TYPE(t_grid), POINTER :: pGrid 
    TYPE(t_global), POINTER :: global     
      
! ==============================================================================
!   Locals
! ==============================================================================

    LOGICAL :: check1,check2,check3
    INTEGER :: cListDimNew,cListTempDim,cListTempDimMax,errorFlag,icg,icl, & 
               iLayer,iLoc,ivg,ivl,iv2c,vListDim,vListTempDim,vListTempDim2, & 
               vListTempDimMax
    INTEGER, DIMENSION(:), ALLOCATABLE :: cListTemp,vList,vListTemp
      
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_BuildVertCellNghbList',&
  'RFLU_ModTopologyUtils.F90')

! ******************************************************************************
!   Allocate temporary memory
! ******************************************************************************
                
    cListTempDimMax = cListDimMax/nLayers
    
    ALLOCATE(cListTemp(cListTempDimMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cListTemp')
    END IF ! global%error        

    vListDim = vListOrigDim

    ALLOCATE(vList(vListDim),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vList')
    END IF ! global%error    

! ******************************************************************************
!   Initialize
! ******************************************************************************             
                
    cListDim = 0
    
    DO icl = 1,cListDimMax
      cList(icl) = 0 
    END DO ! icl         

    DO ivl = 1,vListDim
      vList(ivl) = vListOrig(ivl)
    END DO ! ivl

! ******************************************************************************
!   Loop over layers of cells to be added
! ******************************************************************************             
                                 
    DO iLayer = 1,nLayers    

! ==============================================================================
!     Initialize temporary cell list
! ==============================================================================              

      cListTempDim = 0
      
      DO icl = 1,cListTempDimMax
        cListTemp(icl) = 0
      END DO ! icl  

! ==============================================================================
!     Loop through vertex list and build list of new cells
! ==============================================================================              
        
      DO ivl = 1,vListDim
        ivg = vList(ivl)

! ------------------------------------------------------------------------------
!       For given vertex, loop through adjacent cells
! ------------------------------------------------------------------------------

        DO iv2c = pGrid%v2cInfo(V2C_BEG,ivg),pGrid%v2cInfo(V2C_END,ivg)
          icg = pGrid%v2c(iv2c)

          IF ( iReg == 0 ) THEN 
            check1 = .TRUE.
          ELSE 
            check1 = .FALSE.
          END IF ! iReg
                    
          IF ( ASSOCIATED(pGrid%sc2r) .EQV. .TRUE. ) THEN 
            IF ( (iReg > 0) .AND. (pGrid%sc2r(icg) == iReg) ) THEN
              check2 = .TRUE.
            ELSE 
              check2 = .FALSE.
            END IF ! iReg
            
            IF ( (iReg < 0) .AND. (pGrid%sc2r(icg) /= ABS(iReg)) ) THEN
              check3 = .TRUE.
            ELSE 
              check3 = .FALSE.
            END IF ! iReg            
          ELSE
            check2 = .FALSE.
            check3 = .FALSE.
          END IF ! ASSOCIATED(pGrid%sc2r)

! ------- Check whether cell is in target region -------------------------------

          IF ( (check1 .EQV. .TRUE.) .OR. & 
               (check2 .EQV. .TRUE.) .OR. &
               (check3 .EQV. .TRUE.) ) THEN 

! --------- Add cell to temporary cell list if not already in lists 
      
            IF ( cListDim > 0 ) THEN                             
              CALL BinarySearchInteger(cList(1:cListDim),cListDim,icg,iLoc)
            ELSE 
              iLoc = ELEMENT_NOT_FOUND
            END IF ! cListDim

            IF ( iLoc == ELEMENT_NOT_FOUND ) THEN
              IF ( cListTempDim > 0 ) THEN 
                CALL BinarySearchInteger(cListTemp(1:cListTempDim), & 
                                         cListTempDim,icg,iLoc)
              ELSE 
                iLoc = ELEMENT_NOT_FOUND
              END IF ! cListTempDim

              IF ( iLoc == ELEMENT_NOT_FOUND ) THEN 
                IF ( cListTempDim < cListTempDimMax ) THEN                                                                                      
                  cListTempDim = cListTempDim + 1

                  cListTemp(cListTempDim) = icg

                  IF ( cListTempDim > 1 ) THEN 
                    CALL QuickSortInteger(cListTemp(1:cListTempDim), & 
                                          cListTempDim)
                  END IF ! cListTempDim
                ELSE 
                  CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'cListTemp')
                END IF ! cListTempDim
              END IF ! iLoc            
            END IF ! iLoc
                                  
          END IF ! iReg          
        END DO ! iv2c
      END DO ! ivl    
 
! ==============================================================================
!     Merge cell lists
! ==============================================================================              
      
      IF ( cListDim + cListTempDim <= cListDimMax ) THEN
        DO icl = cListDim+1,cListDim+cListTempDim
          cList(icl) = cListTemp(icl-cListDim)
        END DO ! icl
        
        cListDim = cListDim + cListTempDim
        
        IF ( cListDim /= cListTempDim ) THEN ! Added cells, so need to sort
          CALL QuickSortInteger(cList(1:cListDim),cListDim)
        END IF ! cListDim
      ELSE 
        CALL ErrorStop(global,ERR_EXCEED_DIMENS,__LINE__,'cList')              
      END IF ! cListDim 
      
! ==============================================================================
!     Build new list of vertices from which cells can be added 
! ==============================================================================              
      
      IF( iLayer < nLayers ) THEN 

! ------------------------------------------------------------------------------
!       Allocate temporary memory
! ------------------------------------------------------------------------------

        vListTempDimMax = 8*cListTempDim ! TEMPORARY
        
        ALLOCATE(vListTemp(vListTempDimMax),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vListTemp')
        END IF ! global%error

! ------------------------------------------------------------------------------
!       Build list of vertices from last layer of cells
! ------------------------------------------------------------------------------

        CALL RFLU_BuildCellVertList(global,pGrid,cListTemp(1:cListTempDim), & 
                                    cListTempDim,vListTemp,vListTempDimMax, & 
                                    vListTempDim)

! ------------------------------------------------------------------------------
!       Remove vertices which are shared with current list of vertices
! ------------------------------------------------------------------------------

        CALL RemoveCommonSortedIntegers(vList,vListDim, & 
                                        vListTemp(1:vListTempDim), & 
                                        vListTempDim,vListTempDim2)

! ------------------------------------------------------------------------------
!       Build new list of vertices
! ------------------------------------------------------------------------------

        DEALLOCATE(vList,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vList')
        END IF ! global%error     

        vListDim = vListTempDim2

        ALLOCATE(vList(vListDim),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vList')
        END IF ! global%error   

        DO ivl = 1,vListDim
          vList(ivl) = vListTemp(ivl)
        END DO ! ivl

! ------------------------------------------------------------------------------
!       Deallocate temporary memory
! ------------------------------------------------------------------------------

        DEALLOCATE(vListTemp,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vListTemp')
        END IF ! global%error 
      END IF ! iLayer
    END DO ! iLayer

! ******************************************************************************
!   Deallocate temporary memory
! ******************************************************************************

    DEALLOCATE(vList,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vList')
    END IF ! global%error
        
    DEALLOCATE(cListTemp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cListTemp')
    END IF ! global%error                    
          
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_BuildVertCellNghbList





! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModTopologyUtils


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModTopologyUtils.F90,v $
! Revision 1.9  2008/12/06 08:44:24  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:35  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.6  2006/03/25 21:59:46  haselbac
! Modified building of vert cell nghb list bcos of sype changes
!
! Revision 1.5  2005/10/05 20:09:32  haselbac
! Added fn for building vert list from conn table
!
! Revision 1.4  2005/06/15 20:49:23  haselbac
! Bug fix: Duplicate declaration of errorFlag in RFLU_BuildFaceVertList
!
! Revision 1.3  2005/06/14 17:47:24  haselbac
! Adapted RFLU_BuildFaceVertList for adaptive memory allocation
!
! Revision 1.2  2005/01/17 19:55:20  haselbac
! Added routine to build list of vertices from list of faces, cosmetics
!
! Revision 1.1  2004/12/29 20:58:17  haselbac
! Initial revision
!
! ******************************************************************************











