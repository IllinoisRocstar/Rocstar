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
! Purpose: Suite of routines to construct cell stencils.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModStencilsCells.F90,v 1.18 2008/12/06 08:44:24 mtcampbe Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModStencilsCells

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY: t_grid
  USE ModStencil, ONLY: t_stencil
  USE ModSortSearch
  USE ModMPI

  USE RFLU_ModCellFaceEdgeInfo
  USE RFLU_ModStencilsUtils
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_BuildC2CStencilWrapper, &  
            RFLU_BuildListCC2CStencil, &    
            RFLU_CreateC2CStencilWrapper, &  
            RFLU_DestroyC2CStencilWrapper, & 
            RFLU_DestroyListCC2CStencil, &  
            RFLU_SetInfoC2CStencilWrapper
            
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModStencilsCells.F90,v $ $Revision: 1.18 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS 
 
 
 
 
 

! *******************************************************************************
!
! Purpose: Build 1D cell-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   fnDir       Direction of stencil
!   icgBeg      Beginning global cell index
!   icgEnd      Ending global cell index
!
! Output: None.
!
! Notes: 
!   1. Restricted to hexahedra.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildC2CStencil_1D(pRegion,fnDir,icgBeg,icgEnd)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: fnDir,icgBeg,icgEnd
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: c1,c2,c2cs1DBeg,c2cs1DEnd,degr,errorFlag,icg,icg2, &
               icl,ict,ifg,ifl,iLayer,iloc,iPatch,isl, &
               nCellMembsInfoMax,nCellMembsInfoMaxLoc,nCellMembsInfoMin, &
               nCellMembsInfoMinLoc,nFaces,nLayersInfoMax,nLayersInfoMaxLoc, &
               nLayersInfoMin,nLayersInfoMinLoc,nLayersMax,stencilSizeMax, &
               stencilSizeMin
    INTEGER, DIMENSION(:), ALLOCATABLE :: c2cs1D
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: layerInfo
    REAL(RFREAL) :: fn
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildC2CStencil_1D',&
  'RFLU_ModStencilsCells.F90')

    IF ( (global%myProcid == MASTERPROC) .AND. &
         (global%verbLevel > VERBOSE_NONE) .AND. & 
         (icgEnd > icgBeg) ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building 1D cell-to-cell stencil...' 
      WRITE(STDOUT,'(A,3X,A,1X,I2)') SOLVER_NAME,'Direction:',fnDir           
    END IF ! global%myProcid
    
! ******************************************************************************
!   Set grid pointer and check for required arrays
! ******************************************************************************

    pGrid => pRegion%grid

    IF ( ASSOCIATED(pGrid%hex2f) .EQV. .FALSE. ) THEN 
      CALL ErrorStop(global,ERR_ASSOCIATED,__LINE__,'pGrid%hex2f')
    END IF ! ASSOCIATED

! ******************************************************************************
!   Set variables
! ******************************************************************************

    nLayersMax     = pGrid%c2csInfo%nLayersMax
    stencilSizeMax = pGrid%c2csInfo%nCellMembsMax
    stencilSizeMin = pGrid%c2csInfo%nCellMembsMin

    nCellMembsInfoMax = 0
    nCellMembsInfoMin = HUGE(1)

    nLayersInfoMax = 0
    nLayersInfoMin = HUGE(1)

! ******************************************************************************
!   Allocate temporary memory
! ******************************************************************************

    ALLOCATE(c2cs1D(stencilSizeMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'c2cs1D')
    END IF ! global%error          

    ALLOCATE(layerInfo(X2CS_LAYER_BEG:X2CS_LAYER_END,nLayersMax), &
             STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'layerInfo')
    END IF ! global%error             

! ******************************************************************************
!   Loop over cells
! ******************************************************************************

    DO icg = icgBeg,icgEnd
      ict = pGrid%cellGlob2Loc(1,icg) ! local cell index
      
      IF ( ict /= CELL_TYPE_HEX ) THEN 
        CALL ErrorStop(global,ERR_STENCILMEMBER_INVALID,__LINE__)     
      END IF ! ict
      
      icl = pGrid%cellGlob2Loc(2,icg) ! local cell index

! ==============================================================================
!     Initialize variables
! ==============================================================================

      degr = 0

      DO isl = 1,stencilSizeMax
        c2cs1D(isl) = 0
      END DO ! isl

      DO iLayer = 1,nLayersMax
        layerInfo(X2CS_LAYER_BEG,iLayer) = 0
        layerInfo(X2CS_LAYER_END,iLayer) = 0          
      END DO ! iLayer 

      pGrid%c2cs1D(fnDir,icg)%nLayers = 1

! ******************************************************************************
!     Select cell type and set pointer to cell-to-face connectivity array
! ******************************************************************************

      nFaces = SIZE(pGrid%hex2f,2)
   
! ******************************************************************************
!     Loop over faces of cell
! ******************************************************************************

      DO ifl = 1,nFaces
        iPatch = pGrid%hex2f(1,ifl,icl)
        ifg    = pGrid%hex2f(2,ifl,icl)

        IF ( iPatch == 0 ) THEN ! Interior face
          c1 = pGrid%f2c(1,ifg)
          c2 = pGrid%f2c(2,ifg)

          IF ( c1 == icg ) THEN 
            fn =  pGrid%fn(fnDir,ifg)
          ELSE IF ( c2 == icg ) THEN
            fn = -pGrid%fn(fnDir,ifg)
          ELSE ! defensive programming
            CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END IF ! c1

          IF ( ABS(fn) >= 0.999_RFREAL ) THEN 
            IF ( c1 == icg ) THEN 
              icg2 = c2
            ELSE IF ( c2 == icg ) THEN 
              icg2 = c1
            ELSE ! defensive programming
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END IF ! c1

            IF ( degr > 0 ) THEN ! Search whether already member
              CALL BinarySearchInteger(c2cs1D(1:degr),degr,icg2,iloc)

              IF ( iloc == ELEMENT_NOT_FOUND ) THEN
                IF ( degr < stencilSizeMax ) THEN  
                  degr = degr + 1
                  c2cs1D(degr) = icg2
                  CALL QuickSortInteger(c2cs1D(1:degr),degr)
                END IF ! degr 
              END IF ! iloc
            ELSE ! First member                
              degr = degr + 1
              c2cs1D(degr) = icg2                   
            END IF ! degr
          END IF ! ABS(fn)                  
        ELSE IF ( iPatch > 0 ) THEN ! Boundary face  
! TO DO 
!       Currently do not add boundary faces for constrained reconstruction
! END TO DO 
        ELSE ! Defensive programming
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! iPatch
      END DO ! ifl

      layerInfo(X2CS_LAYER_BEG,1) = 1
      layerInfo(X2CS_LAYER_END,1) = degr

! ==============================================================================             
!     Extend basic stencil. NOTE for 1D stencil do not have to check weight 
!     singularity
! ==============================================================================             

      DO iLayer = 2,nLayersMax
        IF ( degr < stencilSizeMin ) THEN 
          c2cs1DBeg = layerInfo(X2CS_LAYER_BEG,iLayer-1)
          c2cs1DEnd = layerInfo(X2CS_LAYER_END,iLayer-1)              

          CALL RFLU_AddCellLayer_1D(global,pGrid,stencilSizeMax,icg,degr, &
                                    c2cs1DBeg,c2cs1DEnd,c2cs1D,fnDir)

          pGrid%c2cs1D(fnDir,icg)%nLayers = pGrid%c2cs1D(fnDir,icg)%nLayers + 1

          layerInfo(X2CS_LAYER_BEG,iLayer) = &
            layerInfo(X2CS_LAYER_END,iLayer-1) + 1 
          layerInfo(X2CS_LAYER_END,iLayer) = degr
        ELSE 
          EXIT
        END IF ! degr       
      END DO ! iLayer 

! ==============================================================================        
!     Store stencil 
! ==============================================================================

      pGrid%c2cs1D(fnDir,icg)%nCellMembs = degr

      ALLOCATE(pGrid%c2cs1D(fnDir,icg)%cellMembs( &
               pGrid%c2cs1D(fnDir,icg)%nCellMembs),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%c2cs1D%cellMembs')
      END IF ! global%error

      DO isl = 1,pGrid%c2cs1D(fnDir,icg)%nCellMembs
        pGrid%c2cs1D(fnDir,icg)%cellMembs(isl) = c2cs1D(isl)
      END DO ! isl

      ALLOCATE(pGrid%c2cs1D(fnDir,icg)%layerInfo(X2CS_LAYER_BEG: &
               X2CS_LAYER_END,pGrid%c2cs1D(fnDir,icg)%nLayers),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%c2cs1D%layerInfo')
      END IF ! global%error        

      DO iLayer = 1,pGrid%c2cs1D(fnDir,icg)%nLayers
        pGrid%c2cs1D(fnDir,icg)%layerInfo(X2CS_LAYER_BEG,iLayer) = &
          layerInfo(X2CS_LAYER_BEG,iLayer)
        pGrid%c2cs1D(fnDir,icg)%layerInfo(X2CS_LAYER_END,iLayer) = &
          layerInfo(X2CS_LAYER_END,iLayer)            
      END DO ! iLayer 

! ==============================================================================
!     Extract information for later printing 
! ==============================================================================      

      IF ( pGrid%c2cs1D(fnDir,icg)%nLayers < nLayersInfoMin ) THEN 
        nLayersInfoMin    = pGrid%c2cs1D(fnDir,icg)%nLayers
        nLayersInfoMinLoc = icg
      END IF ! pGrid%c2cs1D(fnDir,icg)%nLayers                   

      IF ( pGrid%c2cs1D(fnDir,icg)%nLayers > nLayersInfoMax ) THEN 
        nLayersInfoMax    = pGrid%c2cs1D(fnDir,icg)%nLayers
        nLayersInfoMaxLoc = icg
      END IF ! pGrid%c2cs1D(fnDir,icg)%nLayers 
                              
      IF ( pGrid%c2cs1D(fnDir,icg)%nCellMembs < nCellMembsInfoMin ) THEN 
        nCellMembsInfoMin    = pGrid%c2cs1D(fnDir,icg)%nCellMembs
        nCellMembsInfoMinLoc = icg
      END IF ! pGrid%c2cs1D(fnDir,icg)%nCellMembs                   

      IF ( pGrid%c2cs1D(fnDir,icg)%nCellMembs > nCellMembsInfoMax ) THEN 
        nCellMembsInfoMax = pGrid%c2cs1D(fnDir,icg)%nCellMembs
        nCellMembsInfoMaxLoc = icg        
      END IF ! pGrid%c2cs1D(fnDir,icg)%nCellMembs                        
    END DO ! icg      

! ******************************************************************************
!   Deallocate temporary memory
! ******************************************************************************

    DEALLOCATE(c2cs1D,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'c2cs1D')
    END IF ! global%error 

    DEALLOCATE(layerInfo,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'layerInfo')
    END IF ! global%error                

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( (global%myProcid == MASTERPROC) .AND. &
         (global%verbLevel > VERBOSE_NONE) .AND. & 
         (icgEnd > icgBeg) ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Building 1D cell-to-cell stencil done.'            
    END IF ! global%myProcid   

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildC2CStencil_1D

 
 
  
 
 
 


! *******************************************************************************
!
! Purpose: Build cell-to-cell stencil based on geometry.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   icgBeg              Beginning global cell index
!   icgEnd              Ending global cell index
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildC2CStencil_1D_G(pRegion,dir,icgBeg,icgEnd)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: dir,icgBeg,icgEnd
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: c2cs1DBeg,c2cs1DEnd,degr,errorFlag,icg,icg2,icl,ict, &
               iLayer,iloc,isl,ivg,ivl,iv2c,nBFaceMembs,nBFaceMembsMax, & 
               nBFaceMembsMaxTemp,nCellMembsInfoMax,nCellMembsInfoMaxLoc, &
               nCellMembsInfoMin,nCellMembsInfoMinLoc,nLayersInfoMax, &
               nLayersInfoMaxLoc,nLayersInfoMin,nLayersInfoMinLoc, &
               nLayersMax,nRows,order,orderNominal,stencilSizeMax, &
               stencilSizeMin
    INTEGER, DIMENSION(:), ALLOCATABLE :: c2cs1D
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: bFaceMembs,layerInfo
    REAL(RFREAL) :: rc(XCOORD:ZCOORD)      
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildC2CStencil_1D_G',&
  'RFLU_ModStencilsCells.F90')

    IF ( (global%myProcid == MASTERPROC) .AND. (icgEnd > icgBeg) ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building cell-to-cell stencil...'            

        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
          WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                           pRegion%iRegionGlobal
        END IF ! global%verbLevel
      END IF ! global%verbLevel
    END IF ! global%myProcid        

! ******************************************************************************
!   Set grid pointer and set constants
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Set variables
! ******************************************************************************

    orderNominal   = pGrid%c2csInfo%orderNominal 
    nLayersMax     = pGrid%c2csInfo%nLayersMax     
    stencilSizeMax = pGrid%c2csInfo%nCellMembsMax           
    stencilSizeMin = pGrid%c2csInfo%nCellMembsMin                 

    nCellMembsInfoMax = 0
    nCellMembsInfoMin = HUGE(1)

    nLayersInfoMax = 0
    nLayersInfoMin = HUGE(1)

! ******************************************************************************
!   Allocate temporary memory
! ******************************************************************************

    ALLOCATE(c2cs1D(stencilSizeMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'c2cs1D')
    END IF ! global%error          

    ALLOCATE(layerInfo(X2CS_LAYER_BEG:X2CS_LAYER_END,nLayersMax), &
             STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'layerInfo')
    END IF ! global%error             

! ******************************************************************************
!   Loop over cells and build stencil
! ******************************************************************************

    DO icg = icgBeg,icgEnd
      ict = RFLU_GetGlobalCellType(global,pGrid,icg)
      icl = pGrid%cellGlob2Loc(2,icg)

      rc(XCOORD) = pGrid%cofg(XCOORD,icg)
      rc(YCOORD) = pGrid%cofg(YCOORD,icg)
      rc(ZCOORD) = pGrid%cofg(ZCOORD,icg)      

! ==============================================================================             
!     Initialize
! ==============================================================================             

      degr = 0

      DO isl = 1,stencilSizeMax
        c2cs1D(isl) = 0
      END DO ! isl

      DO iLayer = 1,nLayersMax
        layerInfo(X2CS_LAYER_BEG,iLayer) = 0
        layerInfo(X2CS_LAYER_END,iLayer) = 0          
      END DO ! iLayer 

      pGrid%c2cs1D(dir,icg)%nLayers = 1

! ==============================================================================             
!     Build basic stencil
! ==============================================================================             

      CALL RFLU_BuildC2CStencilBasic_1D(pRegion,icg,stencilSizeMax,dir,degr, &
                                        c2cs1D)

      layerInfo(X2CS_LAYER_BEG,1) = 1
      layerInfo(X2CS_LAYER_END,1) = degr

! ==============================================================================             
!     Extend basic stencil
! ==============================================================================             

      DO iLayer = 2,nLayersMax
        order  = orderNominal

        IF ( degr < stencilSizeMin ) THEN        
          c2cs1DBeg = layerInfo(X2CS_LAYER_BEG,iLayer-1)
          c2cs1DEnd = layerInfo(X2CS_LAYER_END,iLayer-1)              

          CALL RFLU_AddCellLayer_1D_G(global,pGrid,stencilSizeMax,icg,degr, &
                                      c2cs1DBeg,c2cs1DEnd,c2cs1D,rc,dir)

          pGrid%c2cs1D(dir,icg)%nLayers = pGrid%c2cs1D(dir,icg)%nLayers + 1

          layerInfo(X2CS_LAYER_BEG,iLayer) = & 
            layerInfo(X2CS_LAYER_END,iLayer-1) + 1 
          layerInfo(X2CS_LAYER_END,iLayer) = degr
        ELSE 
          EXIT
        END IF ! sCount       
      END DO ! iLayer 

! ==============================================================================        
!     Store stencil 
! ==============================================================================

      pGrid%c2cs1D(dir,icg)%nCellMembs = degr

      ALLOCATE(pGrid%c2cs1D(dir,icg)%cellMembs( &
               pGrid%c2cs1D(dir,icg)%nCellMembs),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%c2cs1D%cellMembs')
      END IF ! global%error

      DO isl = 1,pGrid%c2cs1D(dir,icg)%nCellMembs
        pGrid%c2cs1D(dir,icg)%cellMembs(isl) = c2cs1D(isl)
      END DO ! isl

      ALLOCATE(pGrid%c2cs1D(dir,icg)%layerInfo(X2CS_LAYER_BEG:X2CS_LAYER_END, &
               pGrid%c2cs1D(dir,icg)%nLayers),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%c2cs1D%layerInfo')
      END IF ! global%error        

      DO iLayer = 1,pGrid%c2cs1D(dir,icg)%nLayers
        pGrid%c2cs1D(dir,icg)%layerInfo(X2CS_LAYER_BEG,iLayer) = &
          layerInfo(X2CS_LAYER_BEG,iLayer)
        pGrid%c2cs1D(dir,icg)%layerInfo(X2CS_LAYER_END,iLayer) = &
          layerInfo(X2CS_LAYER_END,iLayer)            
      END DO ! iLayer 
      
! ==============================================================================
!     Extract information for later printing 
! ==============================================================================      

      IF ( pGrid%c2cs1D(dir,icg)%nLayers < nLayersInfoMin ) THEN 
        nLayersInfoMin    = pGrid%c2cs1D(dir,icg)%nLayers
        nLayersInfoMinLoc = icg
      END IF ! pGrid%c2cs1D(dir,icg)%nLayers                   

      IF ( pGrid%c2cs1D(dir,icg)%nLayers > nLayersInfoMax ) THEN 
        nLayersInfoMax    = pGrid%c2cs1D(dir,icg)%nLayers
        nLayersInfoMaxLoc = icg
      END IF ! pGrid%c2cs1D(dir,icg)%nLayers 
                              
      IF ( pGrid%c2cs1D(dir,icg)%nCellMembs < nCellMembsInfoMin ) THEN 
        nCellMembsInfoMin    = pGrid%c2cs1D(dir,icg)%nCellMembs
        nCellMembsInfoMinLoc = icg
      END IF ! pGrid%c2cs1D(dir,icg)%nCellMembs                   

      IF ( pGrid%c2cs1D(dir,icg)%nCellMembs > nCellMembsInfoMax ) THEN 
        nCellMembsInfoMax = pGrid%c2cs1D(dir,icg)%nCellMembs
        nCellMembsInfoMaxLoc = icg        
      END IF ! pGrid%c2cs1D(dir,icg)%nCellMembs    
    END DO ! icg      

! ******************************************************************************
!   Deallocate temporary memory
! ******************************************************************************

    DEALLOCATE(c2cs1D,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'c2cs')
    END IF ! global%error 

    DEALLOCATE(layerInfo,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'layerInfo')
    END IF ! global%error                

! ******************************************************************************
!   Write out information on stencils
! ******************************************************************************

    IF ( (global%myProcid == MASTERPROC) .AND. & 
         (global%verbLevel > VERBOSE_LOW) .AND. & 
         (icgEnd > icgBeg) ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Statistics:'         
      WRITE(STDOUT,'(A,5X,A,2(1X,I3),2(1X,I9))') SOLVER_NAME, &
            'Minimum/maximum number of cell layers: ',nLayersInfoMin, &
            nLayersInfoMax,nLayersInfoMinLoc,nLayersInfoMaxLoc 
      WRITE(STDOUT,'(A,5X,A,2(1X,I3),2(1X,I9))') SOLVER_NAME, &
            'Minimum/maximum number of cell members:',nCellMembsInfoMin, &
            nCellMembsInfoMax,nCellMembsInfoMinLoc,nCellMembsInfoMaxLoc         
    END IF ! global%myProcid

#ifdef CHECK_DATASTRUCT
! ******************************************************************************
!   Data structure output for checking
! ******************************************************************************

    WRITE(STDOUT,'(A)') SOLVER_NAME
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Cell-to-cell stencils'
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Maximum number of layers:', & 
                                   pGrid%c2csInfo%nLayersMax       
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Minimum stencil size:', & 
                                   pGrid%c2csInfo%nCellMembsMin 

    DO icg = 1,pGrid%nCellsTot        
      WRITE(STDOUT,'(A,1X,I6,2(1X,I3),3X,20(1X,I6))') SOLVER_NAME,icg, & 
            pGrid%c2cs1D(dir,icg)%nLayers,pGrid%c2cs1D(dir,icg)%nCellMembs, &
            pGrid%c2cs1D(dir,icg)%cellMembs(1:pGrid%c2cs1D(dir,icg)%nCellMembs)
    END DO ! icg

    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'    
    WRITE(STDOUT,'(A)') SOLVER_NAME
#endif          

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( (global%myProcid == MASTERPROC) .AND. &
         (global%verbLevel > VERBOSE_NONE) .AND. & 
         (icgEnd > icgBeg) ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building cell-to-cell stencil done.'            
    END IF ! global%myProcid      

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildC2CStencil_1D_G
 
 
 
 
 
 
  
  
! *******************************************************************************
!
! Purpose: Build cell-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   icgBeg              Beginning global cell index
!   icgEnd              Ending global cell index
!   addBFaces           Flag indicating whether should add boundary faces
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildC2CStencil(pRegion,icgBeg,icgEnd,addBFaces)

    USE ModInterfaces, ONLY: RFLU_InvertMatrixSVD

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    LOGICAL, INTENT(IN) :: addBFaces
    INTEGER, INTENT(IN) :: icgBeg,icgEnd
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: c2csBeg,c2csEnd,degr,errorFlag,icg,icg2,icl,ict, &
               iLayer,iloc,isl,ivg,ivl,iv2c,nBFaceMembs,nBFaceMembsMax, & 
               nBFaceMembsMaxTemp,nCellMembsInfoMax,nCellMembsInfoMaxLoc, &
               nCellMembsInfoMin,nCellMembsInfoMinLoc,nLayersInfoMax, &
               nLayersInfoMaxLoc,nLayersInfoMin,nLayersInfoMinLoc, &
               nLayersMax,nRows,order,orderNominal,sCount,stencilSizeMax, &
               stencilSizeMin,nCols,iRow,iCol
    INTEGER, DIMENSION(:), ALLOCATABLE :: c2cs,c2csTemp
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: bFaceMembs,layerInfo        
    REAL(RFREAL) :: dx,dy,dz,term
    REAL(RFREAL) :: colMax(4)     
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: a,aInv       
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildC2CStencil',&
  'RFLU_ModStencilsCells.F90')

    IF ( (global%myProcid == MASTERPROC) .AND. (icgEnd > icgBeg) ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building cell-to-cell stencil...'            

        IF ( global%verbLevel > VERBOSE_LOW ) THEN 
          WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                           pRegion%iRegionGlobal
        END IF ! global%verbLevel
      END IF ! global%verbLevel
    END IF ! global%myProcid        

! ******************************************************************************
!   Set grid pointer and set constants
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Set variables
! ******************************************************************************

    orderNominal   = pGrid%c2csInfo%orderNominal 
    nLayersMax     = pGrid%c2csInfo%nLayersMax     
    nBFaceMembsMax = pGrid%c2csInfo%nBFaceMembsMax  
    stencilSizeMax = pGrid%c2csInfo%nCellMembsMax           
    stencilSizeMin = pGrid%c2csInfo%nCellMembsMin                 

    nCellMembsInfoMax = 0
    nCellMembsInfoMin = HUGE(1)

    nLayersInfoMax = 0
    nLayersInfoMin = HUGE(1)
        
    nBFaceMembsMaxTemp = 2*nBFaceMembsMax    
                                              
! ******************************************************************************
!   Allocate temporary memory
! ******************************************************************************

    ALLOCATE(c2cs(stencilSizeMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'c2cs')
    END IF ! global%error          

    ALLOCATE(bFaceMembs(2,nBFaceMembsMaxTemp),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'bFaceMembs')
    END IF ! global%error         

    ALLOCATE(layerInfo(X2CS_LAYER_BEG:X2CS_LAYER_END,nLayersMax), &
             STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'layerInfo')
    END IF ! global%error             

! ******************************************************************************
!   Loop over cells and build stencil
! ******************************************************************************

    DO icg = icgBeg,icgEnd
      ict = RFLU_GetGlobalCellType(global,pGrid,icg)
      icl = pGrid%cellGlob2Loc(2,icg)

! ==============================================================================             
!     Initialize
! ==============================================================================             

      degr = 0

      DO isl = 1,stencilSizeMax
        c2cs(isl) = 0
      END DO ! isl

      DO iLayer = 1,nLayersMax
        layerInfo(X2CS_LAYER_BEG,iLayer) = 0
        layerInfo(X2CS_LAYER_END,iLayer) = 0          
      END DO ! iLayer 

      pGrid%c2cs(icg)%nLayers = 1

! ==============================================================================             
!     Build basic stencil
! ==============================================================================             

      CALL RFLU_BuildC2CStencilBasic(pRegion,icg,stencilSizeMax,degr,c2cs)

      layerInfo(X2CS_LAYER_BEG,1) = 1
      layerInfo(X2CS_LAYER_END,1) = degr

! ==============================================================================             
!     Extend basic stencil
! ==============================================================================             

      DO iLayer = 2,nLayersMax
        order  = orderNominal
        sCount = 0 

! ------------------------------------------------------------------------------
!       Check whether stencil weights are singular
! ------------------------------------------------------------------------------

        IF ( degr >= stencilSizeMin ) THEN 
          nRows = degr
          nCols = pRegion%mixtInput%dimens + 1

          ALLOCATE(a(nRows,nCols),STAT=errorFlag)
          global%error = errorFlag   
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'a')
          END IF ! global%error 

          ALLOCATE(aInv(nRows,nCols),STAT=errorFlag)
          global%error = errorFlag   
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'aInv')
          END IF ! global%error          

          SELECT CASE ( pRegion%mixtInput%dimens ) 
            CASE ( 2 )
              DO isl = 1,degr
                icg2 = c2cs(isl)

                dx = pGrid%cofg(XCOORD,icg2) - pGrid%cofg(XCOORD,icg)
                dy = pGrid%cofg(YCOORD,icg2) - pGrid%cofg(YCOORD,icg)

                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)

                a(isl,1) = term
                a(isl,2) = term*dx
                a(isl,3) = term*dy                                                
              END DO ! isl             
            CASE ( 3 )
              DO isl = 1,degr
                icg2 = c2cs(isl)

                dx = pGrid%cofg(XCOORD,icg2) - pGrid%cofg(XCOORD,icg)
                dy = pGrid%cofg(YCOORD,icg2) - pGrid%cofg(YCOORD,icg)
                dz = pGrid%cofg(ZCOORD,icg2) - pGrid%cofg(ZCOORD,icg)

                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)

                a(isl,1) = term
                a(isl,2) = term*dx
                a(isl,3) = term*dy
                a(isl,4) = term*dz                                
              END DO ! isl            
            CASE DEFAULT
              CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
          END SELECT ! pRegion%mixtInput%dimens

          DO iCol = 1,nCols           
            colMax(iCol) = -HUGE(1.0_RFREAL)

            DO iRow = 1,nRows
              colMax(iCol) = MAX(colMax(iCol),ABS(a(iRow,iCol)))
            END DO ! iRow

            DO iRow = 1,nRows
              a(iRow,iCol) = a(iRow,iCol)/colMax(iCol) 
            END DO ! iRow                     
          END DO ! iCol                     

          CALL RFLU_InvertMatrixSVD(global,nRows,nCols,a,aInv,sCount)

          DEALLOCATE(a,STAT=errorFlag)
          global%error = errorFlag   
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'a')
          END IF ! global%error 

          DEALLOCATE(aInv,STAT=errorFlag)
          global%error = errorFlag   
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'aInv')
          END IF ! global%error                    
        END IF ! degr

! ------------------------------------------------------------------------------
!       Check whether to reject or accept stencil. If singular or too small, 
!       add layer of cells. 
! ------------------------------------------------------------------------------

        IF ( sCount /= 0 .OR. degr < stencilSizeMin ) THEN        
          c2csBeg = layerInfo(X2CS_LAYER_BEG,iLayer-1)
          c2csEnd = layerInfo(X2CS_LAYER_END,iLayer-1)              

          CALL RFLU_AddCellLayer(global,pGrid,stencilSizeMax,icg,degr, &
                                 c2csBeg,c2csEnd,c2cs)
          pGrid%c2cs(icg)%nLayers = pGrid%c2cs(icg)%nLayers + 1

          layerInfo(X2CS_LAYER_BEG,iLayer) = & 
            layerInfo(X2CS_LAYER_END,iLayer-1) + 1 
          layerInfo(X2CS_LAYER_END,iLayer) = degr
        ELSE 
          EXIT
        END IF ! sCount       
      END DO ! iLayer 

! ==============================================================================        
!     Store stencil 
! ==============================================================================

      pGrid%c2cs(icg)%nCellMembs = degr

      ALLOCATE(pGrid%c2cs(icg)%cellMembs(pGrid%c2cs(icg)%nCellMembs), & 
               STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%c2cs%cellMembs')
      END IF ! global%error

      DO isl = 1,pGrid%c2cs(icg)%nCellMembs
        pGrid%c2cs(icg)%cellMembs(isl) = c2cs(isl)
      END DO ! isl

      ALLOCATE(pGrid%c2cs(icg)%layerInfo(X2CS_LAYER_BEG:X2CS_LAYER_END, &
               pGrid%c2cs(icg)%nLayers),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%c2cs%layerInfo')
      END IF ! global%error        

      DO iLayer = 1,pGrid%c2cs(icg)%nLayers
        pGrid%c2cs(icg)%layerInfo(X2CS_LAYER_BEG,iLayer) = &
          layerInfo(X2CS_LAYER_BEG,iLayer)
        pGrid%c2cs(icg)%layerInfo(X2CS_LAYER_END,iLayer) = &
          layerInfo(X2CS_LAYER_END,iLayer)            
      END DO ! iLayer 

! ==============================================================================        
!     Add boundary faces to stencil. If the stencil contains boundary faces,
!     sort them by distance and pick <nBFaceMembsMax> closest ones. NOTE build
!     temporary c2csTemp array which includes cell icg itself in addition to 
!     cells already in stencil. This is done so that boundary face(s) attached
!     to cell icg will also show up in boundary-face stencil.
! ==============================================================================        

      nBFaceMembs = 0

      IF ( addBFaces .EQV. .TRUE. ) THEN 
        ALLOCATE(c2csTemp(pGrid%c2cs(icg)%nCellMembs+1),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'c2csTemp')
        END IF ! global%error

        c2csTemp(1) = icg

        DO isl = 1,pGrid%c2cs(icg)%nCellMembs
          c2csTemp(isl+1) = pGrid%c2cs(icg)%cellMembs(isl)
        END DO ! isl

        CALL RFLU_AddBFaces(pRegion,nBFaceMembsMaxTemp, & 
                            pGrid%c2cs(icg)%nCellMembs+1, &
                            c2csTemp,nBFaceMembs,bFaceMembs)

        DEALLOCATE(c2csTemp,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'c2csTemp')
        END IF ! global%error
      END IF ! addBFaces

      IF ( nBFaceMembs > 0 ) THEN 
        CALL RFLU_SortBFaces(pRegion,pGrid%cofg(XCOORD:ZCOORD,icg), &
                             nBFaceMembs,bFaceMembs(1:2,1:nBFaceMembs))        

        pGrid%c2cs(icg)%nBFaceMembs = MIN(nBFaceMembs,nBFaceMembsMax)

        ALLOCATE(pGrid%c2cs(icg)%bFaceMembs(2,pGrid%c2cs(icg)%nBFaceMembs), & 
                 STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%c2cs%bFaceMembs')
        END IF ! global%error

        DO isl = 1,pGrid%c2cs(icg)%nBFaceMembs
          pGrid%c2cs(icg)%bFaceMembs(1,isl) = bFaceMembs(1,isl)
          pGrid%c2cs(icg)%bFaceMembs(2,isl) = bFaceMembs(2,isl)
        END DO ! isl
      ELSE 
        pGrid%c2cs(icg)%nBFaceMembs = 0
        
        NULLIFY(pGrid%c2cs(icg)%bFaceMembs)
      END IF ! nBFaceMembs 
      
! ==============================================================================
!     Extract information for later printing 
! ==============================================================================      

      IF ( pGrid%c2cs(icg)%nLayers < nLayersInfoMin ) THEN 
        nLayersInfoMin    = pGrid%c2cs(icg)%nLayers
        nLayersInfoMinLoc = icg
      END IF ! pGrid%c2cs(icg)%nLayers                   

      IF ( pGrid%c2cs(icg)%nLayers > nLayersInfoMax ) THEN 
        nLayersInfoMax    = pGrid%c2cs(icg)%nLayers
        nLayersInfoMaxLoc = icg
      END IF ! pGrid%c2cs(icg)%nLayers 
                              
      IF ( pGrid%c2cs(icg)%nCellMembs < nCellMembsInfoMin ) THEN 
        nCellMembsInfoMin    = pGrid%c2cs(icg)%nCellMembs
        nCellMembsInfoMinLoc = icg
      END IF ! pGrid%c2cs(icg)%nCellMembs                   

      IF ( pGrid%c2cs(icg)%nCellMembs > nCellMembsInfoMax ) THEN 
        nCellMembsInfoMax = pGrid%c2cs(icg)%nCellMembs
        nCellMembsInfoMaxLoc = icg        
      END IF ! pGrid%c2cs(icg)%nCellMembs    
    END DO ! icg      

! ******************************************************************************
!   Deallocate temporary memory
! ******************************************************************************

    DEALLOCATE(c2cs,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'c2cs')
    END IF ! global%error 

    DEALLOCATE(bFaceMembs,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'bFaceMembs')
    END IF ! global%error       

    DEALLOCATE(layerInfo,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'layerInfo')
    END IF ! global%error                

! ******************************************************************************
!   Write out information on stencils
! ******************************************************************************

    IF ( (global%myProcid == MASTERPROC) .AND. & 
         (global%verbLevel > VERBOSE_LOW) .AND. & 
         (icgEnd > icgBeg) ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Statistics:'         
      WRITE(STDOUT,'(A,5X,A,2(1X,I3),2(1X,I9))') SOLVER_NAME, &
            'Minimum/maximum number of cell layers: ',nLayersInfoMin, &
            nLayersInfoMax,nLayersInfoMinLoc,nLayersInfoMaxLoc 
      WRITE(STDOUT,'(A,5X,A,2(1X,I3),2(1X,I9))') SOLVER_NAME, &
            'Minimum/maximum number of cell members:',nCellMembsInfoMin, &
            nCellMembsInfoMax,nCellMembsInfoMinLoc,nCellMembsInfoMaxLoc         
    END IF ! global%myProcid

#ifdef CHECK_DATASTRUCT
! ******************************************************************************
!   Data structure output for checking
! ******************************************************************************

    WRITE(STDOUT,'(A)') SOLVER_NAME
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Cell-to-cell stencils'
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Maximum number of layers:', & 
                                   pGrid%c2csInfo%nLayersMax       
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Minimum stencil size:', & 
                                   pGrid%c2csInfo%nCellMembsMin 

    DO icg = 1,pGrid%nCellsTot        
      WRITE(STDOUT,'(A,1X,I6,2(1X,I3),3X,20(1X,I6))') SOLVER_NAME,icg, & 
            pGrid%c2cs(icg)%nLayers,pGrid%c2cs(icg)%nCellMembs, &
            pGrid%c2cs(icg)%cellMembs(1:pGrid%c2cs(icg)%nCellMembs)
    END DO ! icg

    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'    
    WRITE(STDOUT,'(A)') SOLVER_NAME
#endif          

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( (global%myProcid == MASTERPROC) .AND. &
         (global%verbLevel > VERBOSE_NONE) .AND. & 
         (icgEnd > icgBeg) ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building cell-to-cell stencil done.'            
    END IF ! global%myProcid      

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildC2CStencil
  
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Build basic cell-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   icgBeg              Beginning global cell index
!   icgEnd              Ending global cell index
!   addBFaces           Flag indicating whether should add boundary faces
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildC2CStencilBasic(pRegion,icg,stencilSizeMax,degr,c2cs)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: icg,stencilSizeMax
    INTEGER, INTENT(OUT) :: degr
    INTEGER, INTENT(OUT) :: c2cs(stencilSizeMax)
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER, PARAMETER :: NVERTMAX = 8
    INTEGER :: errorFlag,icg2,icl,iclTemp,icl2,ict,ict2,iLoc,ivg,ivl,ivl2, & 
               iv2c,nCells,nCells2,nVert,nVert2,nVert3               
    INTEGER, DIMENSION(NVERTMAX) :: vTemp,vTemp2,vTemp3
    INTEGER, DIMENSION(:), ALLOCATABLE :: icgTemp
    INTEGER, DIMENSION(:,:), POINTER :: x2v  
    REAL(RFREAL) :: iNVert    
    REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: nVertSharedFrac    
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildC2CStencilBasic',&
  'RFLU_ModStencilsCells.F90')

! ******************************************************************************
!   Set grid pointer 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Build basic stencil
! ******************************************************************************

    ict = RFLU_GetGlobalCellType(global,pGrid,icg)
    icl = pGrid%cellGlob2Loc(2,icg)

! ==============================================================================             
!   Build temporary list of cells which are vertex-neighbors of current cell 
! ==============================================================================             

    SELECT CASE ( ict ) 
      CASE ( CELL_TYPE_TET ) 
        x2v => pGrid%tet2v
      CASE ( CELL_TYPE_HEX ) 
        x2v => pGrid%hex2v
      CASE ( CELL_TYPE_PRI )
        x2v => pGrid%pri2v
      CASE ( CELL_TYPE_PYR ) 
        x2v => pGrid%pyr2v
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! ict
    
    nVert  = SIZE(x2v,1)
    iNVert = 1.0_RFREAL/REAL(nVert,KIND=RFREAL) 
    
    nCells = 0

    DO ivl = 1,nVert
      ivg = x2v(ivl,icl)

      nCells = nCells + pGrid%v2cInfo(V2C_END,ivg) & 
                      - pGrid%v2cInfo(V2C_BEG,ivg) & 
                      + 1
    END DO ! ivl

    ALLOCATE(icgTemp(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'icgTemp')
    END IF ! global%error

    ALLOCATE(nVertSharedFrac(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'nVertSharedFrac')
    END IF ! global%error          

    nCells = 0

    DO ivl = 1,nVert
      ivg = x2v(ivl,icl)

      DO iv2c = pGrid%v2cInfo(V2C_BEG,ivg),pGrid%v2cInfo(V2C_END,ivg)
        nCells = nCells + 1

        icgTemp(nCells) = pGrid%v2c(iv2c)
      END DO ! iv2c
    END DO ! ivl

! ==============================================================================             
!   Sort temporary list and remove duplicates as well as current cell
! ==============================================================================             

    CALL QuickSortInteger(icgTemp(1:nCells),nCells)
    CALL SimplifySortedIntegers(icgTemp(1:nCells),nCells,nCells2)                    
    CALL BinarySearchInteger(icgTemp(1:nCells2),nCells2,icg,iLoc)

    IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
      CALL RemoveInteger(icgTemp(1:nCells2),nCells2,iLoc)
    ELSE 
! TEMPORARY
      WRITE(*,*) 'ERROR! Cell icg not found in own stencil!'
      STOP
! END TEMPORARY    
    END IF ! iLoc                    

! ==============================================================================             
!   Determine number of shared vertices between current cell and neighboring
!   cells and store fraction of shared vertices
! ==============================================================================             

    vTemp(1:nVert) = x2v(1:nVert,icl)
    CALL QuickSortInteger(vTemp(1:nVert),nVert)
                    
    DO iclTemp = 1,nCells2
      icg2 = icgTemp(iclTemp)

      ict2 = RFLU_GetGlobalCellType(global,pGrid,icg2)
      icl2 = pGrid%cellGlob2Loc(2,icg2)            

      SELECT CASE ( ict2 ) 
        CASE ( CELL_TYPE_TET )
          nVert2 = 4
          vTemp2(1:4) = pGrid%tet2v(1:4,icl2)
        CASE ( CELL_TYPE_HEX )
          nVert2 = 8
          vTemp2(1:8) = pGrid%hex2v(1:8,icl2)
        CASE ( CELL_TYPE_PRI ) 
          nVert2 = 6
          vTemp2(1:6) = pGrid%pri2v(1:6,icl2)
        CASE ( CELL_TYPE_PYR )
          nVert2 = 5
          vTemp2(1:5) = pGrid%pyr2v(1:5,icl2)
        CASE DEFAULT
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! ict2

      CALL QuickSortInteger(vTemp2(1:nVert2),nVert2) 
      CALL FindCommonSortedIntegers(vTemp(1:nVert),nVert,vTemp2(1:nVert2), &
                                    nVert2,vTemp3,NVERTMAX,nVert3,errorFlag)
! TEMPORARY
      IF ( errorFlag /= ERR_NONE ) THEN 
        WRITE(*,*) 'ERROR!',errorFlag
      END IF ! errorFlag

      IF ( nVert3 == 0 ) THEN 
        WRITE(*,*) 'ERROR! nVert3 = 0!!!'
      END IF ! nVert3                    
! END TEMPORARY

      nVertSharedFrac(iclTemp) = nVert3*iNVert                          
    END DO ! iclTemp     

! ==============================================================================             
!   Sort fraction of shared vertices and store as basic stencil
! ==============================================================================             

    CALL QuickSortRFREALInteger(nVertSharedFrac(1:nCells2),icgTemp(1:nCells2), &
                                nCells2)
        
    degr = 0            
        
    DO icl2 = nCells2,MAX(nCells2-stencilSizeMax+1,1),-1 
      degr = degr + 1      
      c2cs(degr) = icgTemp(icl2)    
    END DO ! icl2    

    CALL QuickSortInteger(c2cs(1:degr),degr)

! ==============================================================================             
!   Deallocate temporary lists
! ==============================================================================             

    DEALLOCATE(icgTemp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'icgTemp')
    END IF ! global%error

    DEALLOCATE(nVertSharedFrac,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'nVertSharedFrac')
    END IF ! global%error                    
      
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildC2CStencilBasic    
  
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Build basic cell-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   icg                 Global cell index
!   stencilSizeMax	Maximum allowed stencil size
!   dir			Direction in which stencil is to be built
!
! Output: 
!   degr		Stencil size
!   c2cs1D		Stencil
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildC2CStencilBasic_1D(pRegion,icg,stencilSizeMax,dir, &
                                          degr,c2cs1D)

    USE RFLU_ModGeometryTools, ONLY: RFLU_TestVectorCartAxisAligned

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: dir,icg,stencilSizeMax
    INTEGER, INTENT(OUT) :: degr
    INTEGER, INTENT(OUT) :: c2cs1D(stencilSizeMax)
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg2,icl,ict,iLoc,ivg,ivl,iv2c,nCells,nCells2
    INTEGER, DIMENSION(:), ALLOCATABLE :: icgTemp
    REAL(RFREAL) :: dr(XCOORD:ZCOORD),rc(XCOORD:ZCOORD)
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildC2CStencilBasic_1D',&
  'RFLU_ModStencilsCells.F90')

! ******************************************************************************
!   Set grid pointer 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Build basic stencil
! ******************************************************************************

    ict = RFLU_GetGlobalCellType(global,pGrid,icg)
    icl = pGrid%cellGlob2Loc(2,icg)

    IF ( ict /= CELL_TYPE_HEX ) THEN 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! ict
    
    IF ( icl /= icg ) THEN 
! TEMPORARY
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
! END TEMPORARY    
    END IF ! icl

    rc(XCOORD) = pGrid%cofg(XCOORD,icg)
    rc(YCOORD) = pGrid%cofg(YCOORD,icg)
    rc(ZCOORD) = pGrid%cofg(ZCOORD,icg)

! ==============================================================================             
!   Build temporary list of cells which are vertex-neighbors of current cell and
!   aligned with selected direction
! ============================================================================== 
    
    nCells = 0

    DO ivl = 1,8
      ivg = pGrid%hex2v(ivl,icg)

      nCells = nCells + pGrid%v2cInfo(V2C_END,ivg) & 
                      - pGrid%v2cInfo(V2C_BEG,ivg) & 
                      + 1
    END DO ! ivl

    ALLOCATE(icgTemp(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'icgTemp')
    END IF ! global%error

    nCells = 0

    DO ivl = 1,8
      ivg = pGrid%hex2v(ivl,icl)

      DO iv2c = pGrid%v2cInfo(V2C_BEG,ivg),pGrid%v2cInfo(V2C_END,ivg)
        icg2 = pGrid%v2c(iv2c)
        
        dr(XCOORD) = pGrid%cofg(XCOORD,icg2) - rc(XCOORD)
        dr(YCOORD) = pGrid%cofg(YCOORD,icg2) - rc(YCOORD)
        dr(ZCOORD) = pGrid%cofg(ZCOORD,icg2) - rc(ZCOORD)
                
        IF ( RFLU_TestVectorCartAxisAligned(global,dr,dir) .EQV. .TRUE. ) THEN  
          nCells = nCells + 1

          icgTemp(nCells) = icg2
        END IF ! RFLU_TestVectorCartAxisAligned
      END DO ! iv2c
    END DO ! ivl

! ==============================================================================             
!   Sort temporary list and remove duplicates as well as current cell
! ==============================================================================             

    CALL QuickSortInteger(icgTemp(1:nCells),nCells)
    CALL SimplifySortedIntegers(icgTemp(1:nCells),nCells,nCells2)                    
    CALL BinarySearchInteger(icgTemp(1:nCells2),nCells2,icg,iLoc)

    IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
      CALL RemoveInteger(icgTemp(1:nCells2),nCells2,iLoc)
    ELSE 
! TEMPORARY
      WRITE(*,*) 'ERROR! Cell icg not found in own stencil!'
      STOP
! END TEMPORARY    
    END IF ! iLoc                    

! ==============================================================================             
!   Store as basic stencil
! ==============================================================================             
        
    degr = 0            
        
    DO icl = nCells2,MAX(nCells2-stencilSizeMax+1,1),-1 
      degr = degr + 1      
      c2cs1D(degr) = icgTemp(icl)    
    END DO ! icl    

    CALL QuickSortInteger(c2cs1D(1:degr),degr)

! ==============================================================================             
!   Deallocate temporary lists
! ==============================================================================             

    DEALLOCATE(icgTemp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'icgTemp')
    END IF ! global%error
      
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildC2CStencilBasic_1D
  
  
  
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Wrapper routine for building cell-to-cell stencils.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   icgInput            Global cell index
!   constrInput         Flag indicating whether have constrained reconstruction
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildC2CStencilWrapper(pRegion,icgInput,constrInput)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN), OPTIONAL :: icgInput,constrInput
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    LOGICAL :: addBFaces
    INTEGER :: icgBeg,icgEnd
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildC2CStencilWrapper',&
  'RFLU_ModStencilsCells.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

    IF ( .NOT. PRESENT(icgInput) ) THEN 
      icgBeg = 1
      icgEnd = pGrid%nCellsTot
    ELSE 
      icgBeg = icgInput
      icgEnd = icgInput
    END IF ! PRESENT   

    IF ( .NOT. PRESENT(constrInput) ) THEN 
      addBFaces = .TRUE.
    ELSE 
      IF ( constrInput == CONSTR_NONE ) THEN 
        addBFaces = .FALSE.
      ELSE 
        addBFaces = .TRUE.
      END IF ! constrInput 
    END IF ! PRESENT       

! ******************************************************************************
!   Call routines to build stencils
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%stencilDimensCells )
      CASE ( 1 ) 
        CALL RFLU_BuildC2CStencil_1D_G(pRegion,XCOORD,icgBeg,icgEnd)
        
        IF ( pRegion%mixtInput%dimens > 1 ) THEN
          CALL RFLU_BuildC2CStencil_1D_G(pRegion,YCOORD,icgBeg,icgEnd)

          IF ( pRegion%mixtInput%dimens > 2 ) THEN
            CALL RFLU_BuildC2CStencil_1D_G(pRegion,ZCOORD,icgBeg,icgEnd)
          END IF ! pRegion%mixtInput%dimens
        END IF ! pRegion%mixtInput%dimens
      CASE ( 2,3 ) 
        CALL RFLU_BuildC2CStencil(pRegion,icgBeg,icgEnd,addBFaces)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%stencilDimensCells

! ******************************************************************************
!   End
! ******************************************************************************
    
    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildC2CStencilWrapper
    
  
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Build list of cell-to-cell stencils which are constrained.
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

  SUBROUTINE RFLU_BuildListCC2CStencil(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildListCC2CStencil',&
  'RFLU_ModStencilsCells.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Building list of constrained ', & 
                                 'cell-to-cell stencil...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Count and build list of constrained cell-to-cell stencils
! ******************************************************************************

    pGrid%nCellsConstr = 0

    IF ( pRegion%mixtInput%cReconstCells > CONSTR_NONE ) THEN 
      DO icg = 1,pGrid%nCellsTot
        IF ( pGrid%c2cs(icg)%nBFaceMembs > 0 ) THEN
          pGrid%nCellsConstr = pGrid%nCellsConstr + 1
        END IF ! pGrid%c2cs(icg)%nBFaceMembs
      END DO ! icg
      
      IF ( pGrid%nCellsConstr > 0 ) THEN
        ALLOCATE(pGrid%icgConstr(pGrid%nCellsConstr),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%icgConstr')
        END IF ! global%error
        
        pGrid%nCellsConstr = 0
        
        DO icg = 1,pGrid%nCellsTot
          IF ( pGrid%c2cs(icg)%nBFaceMembs > 0 ) THEN
            pGrid%nCellsConstr = pGrid%nCellsConstr + 1

            pGrid%icgConstr(pGrid%nCellsConstr) = icg
          END IF ! pGrid%c2cs(icg)%nBFaceMembs
        END DO ! icg
      ELSE 
        NULLIFY(pGrid%icgConstr)        
      END IF ! pGrid%nCellsConstr
    ELSE 
      NULLIFY(pGrid%icgConstr)
    END IF ! pRegion%mixtInput%cReconstCells
        
! ******************************************************************************
!   Print info
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A,A,1X,I5)') SOLVER_NAME,'Number of constrained ', & 
                                       'cell-to-cell stencils:', &
                                       pGrid%nCellsConstr           
    END IF ! global%verbLevel     
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Building list of constrained ', & 
                                 'cell-to-cell stencil done.'   
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildListCC2CStencil   
   
  
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Create 1D cell-to-cell stencil.
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

  SUBROUTINE RFLU_CreateC2CStencil_1D(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,fnDir,fnDirEnd,icg
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateC2CStencil_1D',&
  'RFLU_ModStencilsCells.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating 1D cell-to-cell stencil...'
    END IF ! global%myProcid 

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyC2CStencil_1D(pRegion)

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory and initialize
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%dimens ) 
      CASE ( 1 ) 
        fnDirEnd = 1
      CASE ( 2 ) 
        fnDirEnd = 2
      CASE ( 3 ) 
        fnDirEnd = 3
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%dimens

    ALLOCATE(pGrid%c2cs1D(XCOORD:fnDirEnd,pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%c2cs1D')
    END IF ! global%error               

    DO fnDir = XCOORD,fnDirEnd
      DO icg = 1,pGrid%nCellsTot
        pGrid%c2cs1D(fnDir,icg)%nCellMembs  = 0
        pGrid%c2cs1D(fnDir,icg)%nBFaceMembs = 0
      END DO ! icg
    END DO ! fnDir

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Creating 1D cell-to-cell stencil done.'            
    END IF ! global%myProcid 

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_CreateC2CStencil_1D
    
  
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Create cell-to-cell stencil.
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

  SUBROUTINE RFLU_CreateC2CStencil(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateC2CStencil',&
  'RFLU_ModStencilsCells.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating cell-to-cell stencil...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyC2CStencil(pRegion)

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory and initialize
! ******************************************************************************

    ALLOCATE(pGrid%c2cs(pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%c2cs')
    END IF ! global%error               

    DO icg = 1,pGrid%nCellsTot
      pGrid%c2cs(icg)%nCellMembs  = 0
      pGrid%c2cs(icg)%nBFaceMembs = 0
    END DO ! icg

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Creating cell-to-cell stencil done.'            
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_CreateC2CStencil  
  
   
  
  
  
  
! *******************************************************************************
!
! Purpose: Wrapper routine for creating cell-to-cell stencils.
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

  SUBROUTINE RFLU_CreateC2CStencilWrapper(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
      
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateC2CStencilWrapper',&
  'RFLU_ModStencilsCells.F90')

! ******************************************************************************
!   Call routines to create stencils
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%stencilDimensCells )
      CASE ( 1 ) 
        CALL RFLU_CreateC2CStencil_1D(pRegion)
      CASE ( 2,3 ) 
        CALL RFLU_CreateC2CStencil(pRegion)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%stencilDimensCells

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_CreateC2CStencilWrapper  
  
  
  
   
  
  
  
  
! *******************************************************************************
!
! Purpose: Destroy 1D cell-to-cell stencil.
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

  SUBROUTINE RFLU_DestroyC2CStencil_1D(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,fnDir,fnDirEnd,icg
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyC2CStencil_1D',&
  'RFLU_ModStencilsCells.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Destroying 1D cell-to-cell stencil...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%dimens ) 
      CASE ( 1 ) 
        fnDirEnd = 1
      CASE ( 2 ) 
        fnDirEnd = 2
      CASE ( 3 ) 
        fnDirEnd = 3
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%dimens

    DO fnDir = XCOORD,fnDirEnd
      DO icg = 1,pGrid%nCellsTot
        IF ( pGrid%c2cs1D(fnDir,icg)%nCellMembs > 0 ) THEN       
          DEALLOCATE(pGrid%c2cs1D(fnDir,icg)%cellMembs,STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__, &
                          'pGrid%c2cs1D%cellMembs')
          END IF ! global%error
        
          pGrid%c2cs1D(fnDir,icg)%nCellMembs = 0
        END IF ! pGrid%c2cs1D%nCellMembs

        IF ( pGrid%c2cs1D(fnDir,icg)%nBFaceMembs > 0 ) THEN 
          DEALLOCATE(pGrid%c2cs1D(fnDir,icg)%bFaceMembs,STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__, &
                          'pGrid%c2cs1D%bFaceMembs')
          END IF ! global%error
        
          pGrid%c2cs1D(fnDir,icg)%nBFaceMembs = 0
        END IF ! pGrid%c2cs1D%nBFaceMembs                  
      END DO ! icg
    END DO ! fnDir

    DEALLOCATE(pGrid%c2cs1D,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%c2cs1D')
    END IF ! global%error               

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyC2CStencil_1D(pRegion)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Destroying 1D cell-to-cell stencil done.'
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyC2CStencil_1D 


  
  
  



! *******************************************************************************
!
! Purpose: Destroy cell-to-cell stencil.
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

  SUBROUTINE RFLU_DestroyC2CStencil(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyC2CStencil',&
  'RFLU_ModStencilsCells.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Destroying cell-to-cell stencil...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    DO icg = 1,pGrid%nCellsTot
      IF ( pGrid%c2cs(icg)%nCellMembs > 0 ) THEN       
        DEALLOCATE(pGrid%c2cs(icg)%cellMembs,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%c2cs%cellMembs')
        END IF ! global%error
        
        pGrid%c2cs(icg)%nCellMembs = 0
      END IF ! pGrid%c2cs%nCellMembs

      IF ( pGrid%c2cs(icg)%nBFaceMembs > 0 ) THEN 
        DEALLOCATE(pGrid%c2cs(icg)%bFaceMembs,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%c2cs%bFaceMembs')
        END IF ! global%error
        
        pGrid%c2cs(icg)%nBFaceMembs = 0
      END IF ! pGrid%c2cs%nBFaceMembs                  
    END DO ! icg

    DEALLOCATE(pGrid%c2cs,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%c2cs')
    END IF ! global%error               

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyC2CStencil(pRegion)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Destroying cell-to-cell stencil done.'
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyC2CStencil 

  
    
  
  


! *******************************************************************************
!
! Purpose: Wrapper routine for destroying cell-to-cell stencils.
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

  SUBROUTINE RFLU_DestroyC2CStencilWrapper(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
      
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyC2CStencilWrapper',&
  'RFLU_ModStencilsCells.F90')

! ******************************************************************************
!   Call routines to destroy stencils
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%stencilDimensCells )
      CASE ( 1 ) 
        CALL RFLU_DestroyC2CStencil_1D(pRegion)
      CASE ( 2,3 ) 
        CALL RFLU_DestroyC2CStencil(pRegion)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%stencilDimensCells

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyC2CStencilWrapper  
  
  
  
  
  
  

! *******************************************************************************
!
! Purpose: Destroy list of cell-to-cell stencils which are constrained.
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

  SUBROUTINE RFLU_DestroyListCC2CStencil(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyListCC2CStencil',&
  'RFLU_ModStencilsCells.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Destroying list of ', & 
                                 'constrained cell-to-cell stencil...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Destroy list of constrained cell-to-cell stencils
! ******************************************************************************

    IF ( pRegion%mixtInput%cReconstCells > CONSTR_NONE ) THEN 
      IF ( pGrid%nCellsConstr > 0 ) THEN
        DEALLOCATE(pGrid%icgConstr,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%icgConstr')
        END IF ! global%error
        
        pGrid%nCellsConstr = 0    
      END IF ! pGrid%nCellsConstr
    END IF ! pRegion%mixtInput%cReconstCells
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Destroying list of ', & 
                                 'constrained cell-to-cell stencil done.'   
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyListCC2CStencil   
    
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Nullify 1D cell-to-cell stencil.
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
  
  SUBROUTINE RFLU_NullifyC2CStencil_1D(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global 

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_NullifyC2CStencil_1D',&
  'RFLU_ModStencilsCells.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Nullifying 1D cell-to-cell stencil...'            
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    NULLIFY(pGrid%c2cs1D)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Nullifying 1D cell-to-cell stencil done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_NullifyC2CStencil_1D  
    
  
  
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Nullify cell-to-cell stencil.
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
  
  SUBROUTINE RFLU_NullifyC2CStencil(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global 

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_NullifyC2CStencil',&
  'RFLU_ModStencilsCells.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Nullifying cell-to-cell stencil...'            
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    NULLIFY(pGrid%c2cs)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Nullifying cell-to-cell stencil done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_NullifyC2CStencil 
  
  
    
  



! *******************************************************************************
!
! Purpose: Set 1D cell-to-cell stencil information.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   orderNominal        Nominal order of accuracy
!
! Output: None.
!
! Notes: 
!   1. VERY IMPORTANT: orderNominal is the polynomial order, i.e., the order
!      of the polynomial which is represented exactly. So orderNominal=1 really
!      means a second-order accurate representation. 
!
! ******************************************************************************

  SUBROUTINE RFLU_SetInfoC2CStencil_1D(pRegion,orderNominal)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: orderNominal
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: nBFaceMembsMax,nLayersMax,stencilSizeMax,stencilSizeMin
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_SetInfoC2CStencil_1D',&
  'RFLU_ModStencilsCells.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Setting 1D cell-to-cell stencil information...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Set stencil information
! ******************************************************************************

    nLayersMax     = orderNominal+1
    nBFaceMembsMax = 0              ! TEMPORARY
    stencilSizeMin = orderNominal+1 ! No difference between min and max value
    stencilSizeMax = orderNominal+1

    pGrid%c2csInfo%orderNominal   = orderNominal
    pGrid%c2csInfo%nLayersMax     = nLayersMax      
    pGrid%c2csInfo%nBFaceMembsMax = nBFaceMembsMax 
    pGrid%c2csInfo%nCellMembsMax  = stencilSizeMax    
    pGrid%c2csInfo%nCellMembsMin  = stencilSizeMin

! ******************************************************************************
!   Print stencil information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
            'Maximum allowed number of cell layers in 1D stencil:  ', &
             nLayersMax
      WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
            'Minimum required number of cell members in 1D stencil:', &
             stencilSizeMin                          
      WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
            'Maximum allowed number of cell members in 1D stencil: ', &
             stencilSizeMax   
    END IF ! global%myProcid
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                             'Setting 1D cell-to-cell stencil information done.'
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_SetInfoC2CStencil_1D







! *******************************************************************************
!
! Purpose: Set cell-to-cell stencil information.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   orderNominal        Nominal order of accuracy
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetInfoC2CStencil(pRegion,orderNominal)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: orderNominal
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: nBFaceMembsMax,nLayersMax,stencilSizeMax,stencilSizeMin
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_SetInfoC2CStencil',&
  'RFLU_ModStencilsCells.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Setting cell-to-cell stencil information...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Set stencil information. NOTE nBFaceMembsMax must be one less than the 
!   number of unknowns (columns), otherwise LAPACK routine used for constrained  
!   least-squares problem always gives trivial solution.
! ******************************************************************************

    nLayersMax     = 6
! TEMPORARY
!    nBFaceMembsMax = RFLU_ComputeStencilSize(global,pRegion%mixtInput%dimens, &
!                                             1,orderNominal) - 1  
    nBFaceMembsMax = 9
! END TEMPORARY
    stencilSizeMin = RFLU_ComputeStencilSize(global,pRegion%mixtInput%dimens, &
                                             1,orderNominal)
    stencilSizeMax = 10*stencilSizeMin       

    pGrid%c2csInfo%orderNominal   = orderNominal
    pGrid%c2csInfo%nLayersMax     = nLayersMax      
    pGrid%c2csInfo%nBFaceMembsMax = nBFaceMembsMax 
    pGrid%c2csInfo%nCellMembsMax  = stencilSizeMax    
    pGrid%c2csInfo%nCellMembsMin  = stencilSizeMin

! ******************************************************************************
!   Print stencil information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
            'Maximum allowed number of cell layers:           ',nLayersMax
      WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
            'Minimum required number of cell members:         ',stencilSizeMin                          
      WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
            'Maximum allowed number of cell members:          ',stencilSizeMax 
      WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
            'Maximum allowed number of boundary face members: ',nBFaceMembsMax    
    END IF ! global%myProcid
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Setting cell-to-cell stencil information done.'
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_SetInfoC2CStencil







! *******************************************************************************
!
! Purpose: Wrapper routine for setting info for cell-to-cell stencils.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   orderNominal        Nominal order of accuracy
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetInfoC2CStencilWrapper(pRegion,orderNominal)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: orderNominal
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
      
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_SetInfoC2CStencilWrapper',&
  'RFLU_ModStencilsCells.F90')

! ******************************************************************************
!   Call routines to set info for stencils
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%stencilDimensCells )
      CASE ( 1 ) 
        CALL RFLU_SetInfoC2CStencil_1D(pRegion,orderNominal)
      CASE ( 2,3 ) 
        CALL RFLU_SetInfoC2CStencil(pRegion,orderNominal)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%stencilDimensCells

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_SetInfoC2CStencilWrapper 






! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModStencilsCells


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModStencilsCells.F90,v $
! Revision 1.18  2008/12/06 08:44:24  mtcampbe
! Updated license.
!
! Revision 1.17  2008/11/19 22:17:35  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.16  2007/07/08 21:45:03  gzheng
! changed the PRESENT is used for PGI compiler
!
! Revision 1.15  2007/03/28 18:17:49  haselbac
! Bug fix: Incorrect format statement, too many close parentheses
!
! Revision 1.14  2007/03/06 18:07:56  haselbac
! Added capability to build 1D stencils based on geometry (not on hex2f)
!
! Revision 1.13  2007/02/27 13:07:06  haselbac
! Enabled 1d computations
!
! Revision 1.12  2006/12/15 13:26:36  haselbac
! Fixed bug in format statement, found by ifort
!
! Revision 1.11  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.10  2006/04/07 14:51:21  haselbac
! Adapted to new stencilDimens param
!
! Revision 1.9  2006/04/01 16:41:05  haselbac
! Cosmetics only
!
! Revision 1.8  2006/04/01 15:56:52  haselbac
! Cosmetics only
!
! Revision 1.7  2006/03/20 13:55:37  haselbac
! Completely redone building of basic stencil (see notes)
!
! Revision 1.6  2006/03/09 15:04:40  haselbac
! Bug fix and put routines in right order
!
! Revision 1.5  2006/03/09 14:08:23  haselbac
! Cosmetics, removed CC2C list routines to prevent huge output from rflupart
!
! Revision 1.4  2006/01/06 22:13:43  haselbac
! Renamed routines with shorter names bcos of new 1D routines
!
! Revision 1.3  2005/12/25 15:33:22  haselbac
! Added cell-specific constraint flag
!
! Revision 1.2  2005/10/27 19:19:36  haselbac
! Changed names, clean-up
!
! Revision 1.1  2005/10/05 14:33:44  haselbac
! Initial revision
!
! ******************************************************************************

























