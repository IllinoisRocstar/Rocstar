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
! Purpose: Suite of routines to construct boundary-face stencils.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModStencilsBFaces.F90,v 1.10 2008/12/06 08:44:24 mtcampbe Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModStencilsBFaces

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
  PUBLIC :: RFLU_BuildBF2CStencilWrapper, &
            RFLU_CreateBF2CStencilWrapper, &  
            RFLU_DestroyBF2CStencilWrapper, &
            RFLU_SetInfoBF2CStencilWrapper
            
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModStencilsBFaces.F90,v $ $Revision: 1.10 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
 
 
 

! *******************************************************************************
!
! Purpose: Build 1D boundary-face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPatch              Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildBF2CStencil_1D(pRegion,pPatch)

    USE RFLU_ModPatchUtils, ONLY: RFLU_GetPatchNormalDirection

    USE ModTools, ONLY: FloatEqual

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_patch), POINTER :: pPatch     
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    LOGICAL :: fnDirFlag
    INTEGER :: fnDir,f2cs1DBeg,f2cs1DEnd,degr,errorFlag,ifg,iLayer,iloc, &
               iPatch,isl,nBFaceMembs,nBFaceMembsMax, &
               nCellMembsInfoMax,nCellMembsInfoMaxLoc,nCellMembsInfoMin, &
               nCellMembsInfoMinLoc,nFaces,nLayersInfoMax,nLayersInfoMaxLoc, &
               nLayersInfoMin,nLayersInfoMinLoc,nLayersMax,stencilSizeMax, &
               stencilSizeMin
    INTEGER, DIMENSION(:), ALLOCATABLE :: f2cs1D
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: layerInfo      
    REAL(RFREAL) :: nx,ny,nz,nm  
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid            

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildBF2CStencil_1D',&
  'RFLU_ModStencilsBFaces.F90')

    IF ( global%myProcid == MASTERPROC ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                                 'Building 1D boundary-face-to-cell stencil...'
                               
        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A,1X,I4)') SOLVER_NAME,'Patch:', &
                                         pPatch%iPatchLocal 
        END IF ! global%verbLevel
      END IF ! global%verbLevel                                                                
    END IF ! global%myProcid

! ******************************************************************************
!   Set grid pointer and check for required arrays
! ******************************************************************************

    pGrid => pRegion%grid

    IF ( ASSOCIATED(pGrid%hex2f) .EQV. .FALSE. ) THEN 
      CALL ErrorStop(global,ERR_ASSOCIATED,__LINE__,'pGrid%hex2f')
    END IF ! ASSOCIATED
    
! ******************************************************************************
!   For non-virtual patches, build stencil
! ******************************************************************************

    IF ( pPatch%bcType /= BC_VIRTUAL ) THEN 
      nLayersMax     = pPatch%bf2cs1DInfo%nLayersMax     
      nBFaceMembsMax = pPatch%bf2cs1DInfo%nBFaceMembsMax  
      stencilSizeMax = pPatch%bf2cs1DInfo%nCellMembsMax           
      stencilSizeMin = pPatch%bf2cs1DInfo%nCellMembsMin

      nCellMembsInfoMax = 0
      nCellMembsInfoMin = HUGE(1)

      nLayersInfoMax = 0
      nLayersInfoMin = HUGE(1) 

      IF ( pPatch%flatFlag .EQV. .FALSE. ) THEN 
        CALL ErrorStop(global,ERR_PATCH_NOT_FLAT,__LINE__)
      ELSE 
        CALL RFLU_GetPatchNormalDirection(global,pPatch,fnDir,fnDirFlag)

        IF ( fnDirFlag .EQV. .FALSE. ) THEN  
          CALL ErrorStop(global,ERR_PATCH_NOT_ALIGNED,__LINE__)  
        END IF ! FloatEqual                          
      END IF ! pPatch%flatFlag

      IF ( pPatch%renumFlag .EQV. .FALSE. ) THEN 
        CALL ErrorStop(global,ERR_PATCH_RENUMFLAG,__LINE__)
      END IF ! pPatch%renumFlag 

! ==============================================================================             
!     Allocate temporary memory
! ==============================================================================             

      ALLOCATE(f2cs1D(stencilSizeMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'f2cs')
      END IF ! global%error 

      ALLOCATE(layerInfo(X2CS_LAYER_BEG:X2CS_LAYER_END,nLayersMax), &
               STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'layerInfo')
      END IF ! global%error            

! ==============================================================================             
!     Loop over faces
! ==============================================================================             

      DO ifg = 1,pPatch%nBFaces

! ------------------------------------------------------------------------------
!       Initialize variables
! ------------------------------------------------------------------------------

        degr = 0

        DO isl = 1,stencilSizeMax
          f2cs1D(isl) = 0
        END DO ! isl

        DO iLayer = 1,nLayersMax
          layerInfo(X2CS_LAYER_BEG,iLayer) = 0
          layerInfo(X2CS_LAYER_END,iLayer) = 0          
        END DO ! iLayer               

! ------------------------------------------------------------------------------
!       Build basic stencil consisting of cells abutting face
! ------------------------------------------------------------------------------

        degr = 1            
      
        f2cs1D(1) = pPatch%bf2c(ifg)

        pPatch%bf2cs1D(ifg)%nLayers = 1        

        layerInfo(X2CS_LAYER_BEG,1) = 1
        layerInfo(X2CS_LAYER_END,1) = degr

! ------------------------------------------------------------------------------
!       Extend basic stencil. NOTE for 1D stencil do not have to check weight 
!       singularity
! ------------------------------------------------------------------------------

        DO iLayer = 2,nLayersMax
          IF ( degr < stencilSizeMin ) THEN 
            f2cs1DBeg = layerInfo(X2CS_LAYER_BEG,iLayer-1)
            f2cs1DEnd = layerInfo(X2CS_LAYER_END,iLayer-1)              

            CALL RFLU_AddCellLayer_1D(global,pGrid,stencilSizeMax,0,degr, &
                                      f2cs1DBeg,f2cs1DEnd,f2cs1D,fnDir)

            pPatch%bf2cs1D(ifg)%nLayers = pPatch%bf2cs1D(ifg)%nLayers + 1

            layerInfo(X2CS_LAYER_BEG,iLayer) = &
              layerInfo(X2CS_LAYER_END,iLayer-1) + 1 
            layerInfo(X2CS_LAYER_END,iLayer) = degr
          ELSE 
            EXIT
          END IF ! degr       
        END DO ! iLayer 

! ------------------------------------------------------------------------------
!       Store stencil and layer info
! ------------------------------------------------------------------------------

        pPatch%bf2cs1D(ifg)%nCellMembs = degr

        ALLOCATE(pPatch%bf2cs1D(ifg)%cellMembs(pPatch%bf2cs1D(ifg)%nCellMembs), &
                 STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, & 
                         'pPatch%bf2cs1D%cellMembs')
        END IF ! global%error

        DO isl = 1,pPatch%bf2cs1D(ifg)%nCellMembs
          pPatch%bf2cs1D(ifg)%cellMembs(isl) = f2cs1D(isl)
        END DO ! isl

        ALLOCATE(pPatch%bf2cs1D(ifg)%layerInfo(X2CS_LAYER_BEG:X2CS_LAYER_END, &
                 pPatch%bf2cs1D(ifg)%nLayers),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, &
                         'pPatch%bf2cs1D%layerInfo')
        END IF ! global%error        

        DO iLayer = 1,pPatch%bf2cs1D(ifg)%nLayers
          pPatch%bf2cs1D(ifg)%layerInfo(X2CS_LAYER_BEG,iLayer) = &
            layerInfo(X2CS_LAYER_BEG,iLayer)
          pPatch%bf2cs1D(ifg)%layerInfo(X2CS_LAYER_END,iLayer) = &
            layerInfo(X2CS_LAYER_END,iLayer)            
        END DO ! iLayer 

! ------------------------------------------------------------------------------
!       Add boundary faces to stencil
! ------------------------------------------------------------------------------

        nBFaceMembs = 0      

! ------------------------------------------------------------------------------
!       Extract information for later printing 
! ------------------------------------------------------------------------------

        IF ( pPatch%bf2cs1D(ifg)%nLayers < nLayersInfoMin ) THEN 
          nLayersInfoMin    = pPatch%bf2cs1D(ifg)%nLayers
          nLayersInfoMinLoc = ifg
        END IF ! pPatch%bf2cs(ifg)%nLayers                   

        IF ( pPatch%bf2cs1D(ifg)%nLayers > nLayersInfoMax ) THEN 
          nLayersInfoMax    = pPatch%bf2cs1D(ifg)%nLayers
          nLayersInfoMaxLoc = ifg
        END IF ! pPatch%bf2cs(ifg)%nLayers

        IF ( pPatch%bf2cs1D(ifg)%nCellMembs < nCellMembsInfoMin ) THEN 
          nCellMembsInfoMin    = pPatch%bf2cs1D(ifg)%nCellMembs
          nCellMembsInfoMinLoc = ifg
        END IF ! pPatch%bf2cs(ifg)%nCellMembs                   

        IF ( pPatch%bf2cs1D(ifg)%nCellMembs > nCellMembsInfoMax ) THEN 
          nCellMembsInfoMax    = pPatch%bf2cs1D(ifg)%nCellMembs
          nCellMembsInfoMaxLoc = ifg        
        END IF ! pPatch%bf2cs1D(ifg)%nCellMembs                                                                  
      END DO ! ifg  

! ==============================================================================             
!     Write out information on stencils
! ==============================================================================             

      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel > VERBOSE_LOW ) THEN 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Statistics:'         
        WRITE(STDOUT,'(A,7X,A,2(1X,I3),2(1X,I9))') SOLVER_NAME, &
              'Minimum/maximum number of cell layers: ',nLayersInfoMin, &
              nLayersInfoMax,nLayersInfoMinLoc,nLayersInfoMaxLoc 
        WRITE(STDOUT,'(A,7X,A,2(1X,I3),2(1X,I9))') SOLVER_NAME, &
              'Minimum/maximum number of cell members:',nCellMembsInfoMin, &
              nCellMembsInfoMax,nCellMembsInfoMinLoc,nCellMembsInfoMaxLoc         
      END IF ! global%myProcid

! ==============================================================================             
!     Deallocate temporary memory
! ==============================================================================             

      DEALLOCATE(f2cs1D,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'f2cs1D')
      END IF ! global%error

      DEALLOCATE(layerInfo,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'layerInfo')
      END IF ! global%error

#ifdef CHECK_DATASTRUCT
! ==============================================================================             
!     Data structure output for checking
! ==============================================================================             

      WRITE(STDOUT,'(A)') SOLVER_NAME
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Boundary face-to-cell stencils'
      WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Maximum number of layers:', & 
                                     pPatch%bf2cs1DInfo%nLayersMax      
      WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Minimum stencil size:', & 
                                     pPatch%bf2cs1DInfo%nCellMembsMin

      DO ifg = 1,pPatch%nBFaces
        WRITE(STDOUT,'(A,1X,I6,2(1X,I3),3X,20(1X,I6))') SOLVER_NAME,ifg, & 
              pPatch%bf2cs1D(ifg)%nLayers,pPatch%bf2cs1D(ifg)%nCellMembs, & 
              pPatch%bf2cs1D(ifg)%cellMembs(1:pPatch%bf2cs1D(ifg)%nCellMembs)
      END DO ! ifg

      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'    
      WRITE(STDOUT,'(A)') SOLVER_NAME
#endif         
    END IF ! pPatch%bcType

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Building 1D boundary-face-to-cell stencil done.'            
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildBF2CStencil_1D

 
 
 
 
 
 
! *******************************************************************************
!
! Purpose: Build boundary-face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPatch              Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildBF2CStencil(pRegion,pPatch)

    USE ModInterfaces, ONLY: RFLU_InvertMatrixSVD

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_patch), POINTER :: pPatch     
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: degr,errorFlag,f2csBeg,f2csEnd,icg,ifg,ifg2,iLayer,iloc, &
               isl,ivl,iv2c,nBFaceMembs,nBFaceMembsMax,nBFaceMembsMaxTemp, & 
               nCellMembsInfoMax,nCellMembsInfoMaxLoc,nCellMembsInfoMin, &
               nCellMembsInfoMinLoc,nLayersInfoMax,nLayersInfoMaxLoc, &
               nLayersInfoMin,nLayersInfoMinLoc,nLayersMax,nRows,order, &
               orderNominal,sCount,stencilSizeMax,stencilSizeMin,nCols, &
               iRow,iCol
    INTEGER :: bf2v(4)
    INTEGER, DIMENSION(:), ALLOCATABLE :: f2cs
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: bFaceMembs,layerInfo      
    REAL(RFREAL) :: dx,dy,dz,term
    REAL(RFREAL) :: colMax(4)     
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: a,aInv  
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid            

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildBF2CStencil',&
  'RFLU_ModStencilsBFaces.F90')

    IF ( global%myProcid == MASTERPROC ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                                 'Building boundary-face-to-cell stencil...'
                               
        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A,1X,I4)') SOLVER_NAME,'Patch:', &
                                         pPatch%iPatchLocal 
        END IF ! global%verbLevel
      END IF ! global%verbLevel                                                                
    END IF ! global%myProcid

! ******************************************************************************
!   Set grid pointer 
! ******************************************************************************

    pGrid => pRegion%grid
    
! ******************************************************************************
!   For non-virtual patches, build stencil
! ******************************************************************************

    IF ( pPatch%bcType /= BC_VIRTUAL ) THEN 
      orderNominal   = pPatch%bf2csInfo%orderNominal 
      nLayersMax     = pPatch%bf2csInfo%nLayersMax     
      nBFaceMembsMax = pPatch%bf2csInfo%nBFaceMembsMax  
      stencilSizeMax = pPatch%bf2csInfo%nCellMembsMax           
      stencilSizeMin = pPatch%bf2csInfo%nCellMembsMin

      nCellMembsInfoMax = 0
      nCellMembsInfoMin = HUGE(1)

      nLayersInfoMax = 0
      nLayersInfoMin = HUGE(1) 

      nBFaceMembsMaxTemp = 2*nBFaceMembsMax 

      IF ( pPatch%renumFlag .EQV. .FALSE. ) THEN 
        CALL ErrorStop(global,ERR_PATCH_RENUMFLAG,__LINE__)
      END IF ! pPatch%renumFlag 

! ==============================================================================             
!     Allocate temporary memory
! ==============================================================================             

      ALLOCATE(f2cs(stencilSizeMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'f2cs')
      END IF ! global%error 

      ALLOCATE(bFaceMembs(2,stencilSizeMax),STAT=errorFlag)
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

! ==============================================================================             
!     Loop over faces
! ==============================================================================             

      DO ifg = 1,pPatch%nBFaces

! ------------------------------------------------------------------------------
!       Initialize
! ------------------------------------------------------------------------------

        DO isl = 1,stencilSizeMax
          f2cs(isl) = 0
        END DO ! isl           

        DO iLayer = 1,nLayersMax
          layerInfo(X2CS_LAYER_BEG,iLayer) = 0
          layerInfo(X2CS_LAYER_END,iLayer) = 0          
        END DO ! iLayer              

! ------------------------------------------------------------------------------
!       Build basic stencil
! ------------------------------------------------------------------------------

        degr = 0

        IF ( pPatch%bcType == BC_SYMMETRY ) THEN
                
! ------- Build basic stencil consisting of cells meeting at vertices of face.
!         NOTE this will lead to larger stencils, depending on the specified
!         minimum size, because first layer will already include string of
!         cells which are on boundary. Enlarging the stencil by an additional
!         layer will increase support in direction away from boundary, but
!         also along boundary, so stencils using this approach are large along
!         the boundary. Need this kind of construction for symmetry boundaries
!         so that stencils themselves are also symmetric.

          DO ivl = 1,4                     
            IF ( pPatch%bf2v(ivl,ifg) /= VERT_NONE ) THEN
              bf2v(ivl) = pPatch%bv(pPatch%bf2v(ivl,ifg))
            ELSE                                           
              bf2v(ivl) = VERT_NONE
            END IF ! pPatch%bf2v
          END DO ! ivl                                                          

          CALL RFLU_AddFaceVertNeighbs(global,pGrid,stencilSizeMax,bf2v,degr, & 
                                       f2cs)

          pPatch%bf2cs(ifg)%nLayers = 1

          layerInfo(X2CS_LAYER_BEG,1) = 1                  
          layerInfo(X2CS_LAYER_END,1) = degr
        ELSE

! ------- Build basic stencil consisting of cell abutting face. NOTE this will
!         lead to smaller stencils than approach above.

          degr = 1

          pPatch%bf2cs(ifg)%nLayers = 1

          f2cs(1) = pPatch%bf2c(ifg)

          layerInfo(X2CS_LAYER_BEG,1) = 1
          layerInfo(X2CS_LAYER_END,1) = degr
        END IF ! bcType

! ------------------------------------------------------------------------------
!       Extend basic stencil
! ------------------------------------------------------------------------------

        DO iLayer = 2,nLayersMax
          order  = orderNominal
          sCount = 0 
              
! ------- Check whether stencil weights are singular ---------------------------

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
                  icg = f2cs(isl)

                  dx = pGrid%cofg(XCOORD,icg) - pPatch%fc(XCOORD,ifg)
                  dy = pGrid%cofg(YCOORD,icg) - pPatch%fc(YCOORD,ifg)

                  term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)

                  a(isl,1) = term
                  a(isl,2) = term*dx
                  a(isl,3) = term*dy                                
                END DO ! isl             
              CASE ( 3 )
                DO isl = 1,degr
                  icg = f2cs(isl)

                  dx = pGrid%cofg(XCOORD,icg) - pPatch%fc(XCOORD,ifg)
                  dy = pGrid%cofg(YCOORD,icg) - pPatch%fc(YCOORD,ifg)
                  dz = pGrid%cofg(ZCOORD,icg) - pPatch%fc(ZCOORD,ifg)

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

! ------- Check whether to reject or accept stencil. If singular or too small, 
!         add layer of cells. Pass 0 instead of ifg because want to prevent  
!         present cell from being added, not present face. If pass present   
!         face, could actually prevent a proper cell from being added in     
!         rare cases...                                                      

          IF ( sCount /= 0 .OR. degr < stencilSizeMin ) THEN
            f2csBeg = layerInfo(X2CS_LAYER_BEG,iLayer-1)
            f2csEnd = layerInfo(X2CS_LAYER_END,iLayer-1)              

            CALL RFLU_AddCellLayer(global,pGrid,stencilSizeMax,0,degr, &
                                   f2csBeg,f2csEnd,f2cs)

            pPatch%bf2cs(ifg)%nLayers = pPatch%bf2cs(ifg)%nLayers + 1              

            layerInfo(X2CS_LAYER_BEG,iLayer) = & 
              layerInfo(X2CS_LAYER_END,iLayer-1) + 1 
            layerInfo(X2CS_LAYER_END,iLayer) = degr                                               
          ELSE 
            EXIT
          END IF ! sCount                                
        END DO ! iLayer 

! ----- Store stencil ----------------------------------------------------------

        pPatch%bf2cs(ifg)%nCellMembs = degr

        ALLOCATE(pPatch%bf2cs(ifg)%cellMembs(pPatch%bf2cs(ifg)%nCellMembs), &
                 STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, & 
                         'pPatch%bf2cs%cellMembs')
        END IF ! global%error

        DO isl = 1,pPatch%bf2cs(ifg)%nCellMembs
          pPatch%bf2cs(ifg)%cellMembs(isl) = f2cs(isl)
        END DO ! isl

        ALLOCATE(pPatch%bf2cs(ifg)%layerInfo(X2CS_LAYER_BEG:X2CS_LAYER_END, &
                 pPatch%bf2cs(ifg)%nLayers),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, &
                         'pPatch%bf2cs%layerInfo')
        END IF ! global%error        

        DO iLayer = 1,pPatch%bf2cs(ifg)%nLayers
          pPatch%bf2cs(ifg)%layerInfo(X2CS_LAYER_BEG,iLayer) = &
            layerInfo(X2CS_LAYER_BEG,iLayer)
          pPatch%bf2cs(ifg)%layerInfo(X2CS_LAYER_END,iLayer) = &
            layerInfo(X2CS_LAYER_END,iLayer)            
        END DO ! iLayer 

! ----- Add boundary faces to stencil ------------------------------------------

        nBFaceMembs = 0      

        IF ( (pPatch%bcType == BC_NOSLIPWALL_HFLUX) .OR. & 
             (pPatch%bcType == BC_NOSLIPWALL_TEMP ) .OR. & 
             (pPatch%bcType == BC_INJECTION       ) ) THEN                           
          CALL RFLU_AddBFaces(pRegion,nBFaceMembsMaxTemp, & 
                 pPatch%bf2cs(ifg)%nCellMembs,& 
                 pPatch%bf2cs(ifg)%cellMembs(1:pPatch%bf2cs(ifg)%nCellMembs), & 
                 nBFaceMembs,bFaceMembs)
        END IF ! pPatch%bcType

! ----- Check whether boundary faces were added --------------------------------

        IF ( nBFaceMembs > 0 ) THEN

! ------- Sort boundary faces by distance         

          CALL RFLU_SortBFaces(pRegion,pPatch%fc(XCOORD:ZCOORD,ifg), &
                               nBFaceMembs,bFaceMembs(1:2,1:nBFaceMembs))         

! ------- Remove first face if it is the same as present face in loop

          IF ( bFaceMembs(1,1) == pPatch%iPatchLocal .AND. & 
               bFaceMembs(2,1) == ifg ) THEN 
            DO isl = 2,nBFaceMembs
              bFaceMembs(1,isl-1) = bFaceMembs(1,isl)
              bFaceMembs(2,isl-1) = bFaceMembs(2,isl) 
            END DO ! isl

            nBFaceMembs = nBFaceMembs - 1
          END IF ! bFaceMembs          

          pPatch%bf2cs(ifg)%nBFaceMembs = MIN(nBFaceMembs,nBFaceMembsMax)

          IF ( pPatch%bf2cs(ifg)%nBFaceMembs > 0 ) THEN 
            ALLOCATE(pPatch%bf2cs(ifg)%bFaceMembs(2, & 
                     pPatch%bf2cs(ifg)%nBFaceMembs),STAT=errorFlag)
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) THEN 
              CALL ErrorStop(global,ERR_ALLOCATE,__LINE__, &
                             'pPatch%bf2cs%bFaceMembs')
            END IF ! global%error
          ELSE 
            NULLIFY(pPatch%bf2cs(ifg)%bFaceMembs)
          END IF ! pPatch%bf2cs(ifg)%nBFaceMembs

          DO isl = 1,pPatch%bf2cs(ifg)%nBFaceMembs
            pPatch%bf2cs(ifg)%bFaceMembs(1,isl) = bFaceMembs(1,isl)
            pPatch%bf2cs(ifg)%bFaceMembs(2,isl) = bFaceMembs(2,isl)
          END DO ! isl
        ELSE 
          pPatch%bf2cs(ifg)%nBFaceMembs = 0

          NULLIFY(pPatch%bf2cs(ifg)%bFaceMembs)
        END IF ! nBFaceMembs

! ----- Extract information for later printing ---------------------------------

        IF ( pPatch%bf2cs(ifg)%nLayers < nLayersInfoMin ) THEN 
          nLayersInfoMin    = pPatch%bf2cs(ifg)%nLayers
          nLayersInfoMinLoc = ifg
        END IF ! pPatch%bf2cs(ifg)%nLayers                   

        IF ( pPatch%bf2cs(ifg)%nLayers > nLayersInfoMax ) THEN 
          nLayersInfoMax    = pPatch%bf2cs(ifg)%nLayers
          nLayersInfoMaxLoc = ifg
        END IF ! pPatch%bf2cs(ifg)%nLayers

        IF ( pPatch%bf2cs(ifg)%nCellMembs < nCellMembsInfoMin ) THEN 
          nCellMembsInfoMin    = pPatch%bf2cs(ifg)%nCellMembs
          nCellMembsInfoMinLoc = ifg
        END IF ! pPatch%bf2cs(ifg)%nCellMembs                   

        IF ( pPatch%bf2cs(ifg)%nCellMembs > nCellMembsInfoMax ) THEN 
          nCellMembsInfoMax    = pPatch%bf2cs(ifg)%nCellMembs
          nCellMembsInfoMaxLoc = ifg        
        END IF ! pPatch%bf2cs(ifg)%nCellMembs                                                                  
      END DO ! ifg  

! ==============================================================================             
!     Write out information on stencils
! ==============================================================================             

      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel > VERBOSE_LOW ) THEN 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Statistics:'         
        WRITE(STDOUT,'(A,7X,A,2(1X,I3),2(1X,I9))') SOLVER_NAME, &
              'Minimum/maximum number of cell layers: ',nLayersInfoMin, &
              nLayersInfoMax,nLayersInfoMinLoc,nLayersInfoMaxLoc 
        WRITE(STDOUT,'(A,7X,A,2(1X,I3),2(1X,I9))') SOLVER_NAME, &
              'Minimum/maximum number of cell members:',nCellMembsInfoMin, &
              nCellMembsInfoMax,nCellMembsInfoMinLoc,nCellMembsInfoMaxLoc         
      END IF ! global%myProcid

! ==============================================================================             
!     Deallocate temporary memory
! ==============================================================================             

      DEALLOCATE(f2cs,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'f2cs')
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

#ifdef CHECK_DATASTRUCT
! ==============================================================================             
!     Data structure output for checking
! ==============================================================================             

      WRITE(STDOUT,'(A)') SOLVER_NAME
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Boundary face-to-cell stencils'
      WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Maximum number of layers:', & 
                                     pPatch%bf2csInfo%nLayersMax      
      WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Minimum stencil size:', & 
                                     pPatch%bf2csInfo%nCellMembsMin

      DO ifg = 1,pPatch%nBFaces
        WRITE(STDOUT,'(A,1X,I6,2(1X,I3),3X,20(1X,I6))') SOLVER_NAME,ifg, & 
              pPatch%bf2cs(ifg)%nLayers,pPatch%bf2cs(ifg)%nCellMembs, & 
              pPatch%bf2cs(ifg)%cellMembs(1:pPatch%bf2cs(ifg)%nCellMembs)
      END DO ! ifg

      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'    
      WRITE(STDOUT,'(A)') SOLVER_NAME
#endif         
    END IF ! pPatch%bcType

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Building boundary-face-to-cell stencil done.'            
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildBF2CStencil

 
 
  
  
  


! *******************************************************************************
!
! Purpose: Wrapper for building boundary face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPatch              Pointer to patch
!   constrInput         Flag indicating whether have constrained reconstruction
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildBF2CStencilWrapper(pRegion,pPatch,constrInput)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN), OPTIONAL :: constrInput
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    LOGICAL :: addBFaces
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildBF2CStencilWrapper',&
  'RFLU_ModStencilsBFaces.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

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

    SELECT CASE ( pRegion%mixtInput%stencilDimensBFaces )
      CASE ( 1 ) 
        CALL RFLU_BuildBF2CStencil_1D(pRegion,pPatch)
      CASE ( 2,3 ) 
        CALL RFLU_BuildBF2CStencil(pRegion,pPatch)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%stencilDimensBFaces      

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildBF2CStencilWrapper

  
  
  
  
! *******************************************************************************
!
! Purpose: Create 1D boundary-face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPatch              Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_CreateBF2CStencil_1D(pRegion,pPatch)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ifl
    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateBF2CStencil_1D',&
  'RFLU_ModStencilsBFaces.F90')

    IF ( global%myProcid == MASTERPROC ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                                 'Creating 1D boundary-face-to-cell stencil...'
                               
        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A,1X,I4)') SOLVER_NAME,'Patch:', &
                                         pPatch%iPatchLocal 
        END IF ! global%verbLevel
      END IF ! global%verbLevel                                                                
    END IF ! global%myProcid

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyBF2CStencil_1D(pRegion,pPatch)

! ******************************************************************************
!   Allocate memory and initialize
! ******************************************************************************

    IF ( pPatch%bcType /= BC_VIRTUAL ) THEN      
      ALLOCATE(pPatch%bf2cs1D(pPatch%nBFaces),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2cs1D')
      END IF ! global%error

      DO ifl = 1,pPatch%nBFaces
        pPatch%bf2cs1D(ifl)%nCellMembs  = 0
        pPatch%bf2cs1D(ifl)%nBFaceMembs = 0
      END DO ! ifl
    END IF ! pPatch%bcType          

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Creating 1D boundary-face-to-cell stencil done.'            
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_CreateBF2CStencil_1D     
  
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Create boundary-face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPatch              Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_CreateBF2CStencil(pRegion,pPatch)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ifl
    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateBF2CStencil',&
  'RFLU_ModStencilsBFaces.F90')

    IF ( global%myProcid == MASTERPROC ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                                 'Creating boundary-face-to-cell stencil...'
                               
        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A,1X,I4)') SOLVER_NAME,'Patch:', &
                                         pPatch%iPatchLocal 
        END IF ! global%verbLevel
      END IF ! global%verbLevel                                                                
    END IF ! global%myProcid

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyBF2CStencil(pRegion,pPatch)

! ******************************************************************************
!   Allocate memory and initialize
! ******************************************************************************

    IF ( pPatch%bcType /= BC_VIRTUAL ) THEN      
      ALLOCATE(pPatch%bf2cs(pPatch%nBFaces),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2cs')
      END IF ! global%error

      DO ifl = 1,pPatch%nBFaces
        pPatch%bf2cs(ifl)%nCellMembs  = 0
        pPatch%bf2cs(ifl)%nBFaceMembs = 0
      END DO ! ifl
    END IF ! pPatch%bcType          

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Creating boundary-face-to-cell stencil done.'            
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_CreateBF2CStencil    
  
  
  
  
  
  

! *******************************************************************************
!
! Purpose: Wrapper routine for creating boundary face-to-cell stencils.
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

  SUBROUTINE RFLU_CreateBF2CStencilWrapper(pRegion,pPatch)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
      
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateBF2CStencilWrapper',&
  'RFLU_ModStencilsBFaces.F90')

! ******************************************************************************
!   Call routines to create stencils
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%stencilDimensBFaces )
      CASE ( 1 )     
        CALL RFLU_CreateBF2CStencil_1D(pRegion,pPatch) 
      CASE ( 2,3 ) 
        CALL RFLU_CreateBF2CStencil(pRegion,pPatch)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%stencilDimensBFaces

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_CreateBF2CStencilWrapper

  
  



! *******************************************************************************
!
! Purpose: Destroy 1D boundary-face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPatch              Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_DestroyBF2CStencil_1D(pRegion,pPatch)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ifg
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global 

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyBF2CStencil_1D',&
  'RFLU_ModStencilsBFaces.F90')

    IF ( global%myProcid == MASTERPROC ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                                 'Destroying 1D boundary-face-to-cell stencil...'
                               
        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A,1X,I4)') SOLVER_NAME,'Patch:', &
                                         pPatch%iPatchLocal 
        END IF ! global%verbLevel
      END IF ! global%verbLevel                                                                
    END IF ! global%myProcid

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    IF ( pPatch%bcType /= BC_VIRTUAL ) THEN  
      DO ifg = 1,pPatch%nBFaces
        DEALLOCATE(pPatch%bf2cs1D(ifg)%cellMembs,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__, & 
                         'pPatch%bf2cs1D%cellMembs')
        END IF ! global%error

        IF ( pPatch%bf2cs1D(ifg)%nBFaceMembs > 0 ) THEN
          DEALLOCATE(pPatch%bf2cs1D(ifg)%bFaceMembs,STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__, & 
                           'pPatch%bf2cs1D%bFaceMembs')
          END IF ! global%error 
        END IF ! pPatch%bf2cs1D%nBFaceMembs                   
      END DO ! ifg

      DEALLOCATE(pPatch%bf2cs1D,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bf2cs1D')
      END IF ! global%error
    END IF ! pPatch       

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyBF2CStencil_1D(pRegion,pPatch)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying 1D '// &
                               'boundary-face-to-cell stencil done.'
    END IF ! global%verbLevel   

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyBF2CStencil_1D 

  
  


! *******************************************************************************
!
! Purpose: Destroy boundary-face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPatch              Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_DestroyBF2CStencil(pRegion,pPatch)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ifg
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global 

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyBF2CStencil',&
  'RFLU_ModStencilsBFaces.F90')

    IF ( global%myProcid == MASTERPROC ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                                 'Destroying boundary-face-to-cell stencil...'
                               
        IF ( global%verbLevel > VERBOSE_LOW ) THEN
          WRITE(STDOUT,'(A,3X,A,1X,I4)') SOLVER_NAME,'Patch:', &
                                         pPatch%iPatchLocal 
        END IF ! global%verbLevel
      END IF ! global%verbLevel                                                                
    END IF ! global%myProcid

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    IF ( pPatch%bcType /= BC_VIRTUAL ) THEN  
      DO ifg = 1,pPatch%nBFaces
        DEALLOCATE(pPatch%bf2cs(ifg)%cellMembs,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__, & 
                         'pPatch%bf2cs%cellMembs')
        END IF ! global%error

        IF ( pPatch%bf2cs(ifg)%nBFaceMembs > 0 ) THEN
          DEALLOCATE(pPatch%bf2cs(ifg)%bFaceMembs,STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__, & 
                           'pPatch%bf2cs%bFaceMembs')
          END IF ! global%error 
        END IF ! pPatch%bf2cs%nBFaceMembs                   
      END DO ! ifg

      DEALLOCATE(pPatch%bf2cs,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bf2cs')
      END IF ! global%error
    END IF ! pPatch       

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyBF2CStencil(pRegion,pPatch)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Destroying boundary-face-to-cell stencil done.'
    END IF ! global%verbLevel   

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyBF2CStencil 







! *******************************************************************************
!
! Purpose: Wrapper routine for destroying boundary face-to-cell stencils.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pPatch      Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_DestroyBF2CStencilWrapper(pRegion,pPatch)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
      
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyBF2CStencilWrapper',&
  'RFLU_ModStencilsBFaces.F90')

! ******************************************************************************
!   Call routines to destroy stencils
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%stencilDimensBFaces )
      CASE ( 1 )    
        CALL RFLU_DestroyBF2CStencil_1D(pRegion,pPatch)
      CASE ( 2,3 ) 
        CALL RFLU_DestroyBF2CStencil(pRegion,pPatch)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%stencilDimensBFaces

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyBF2CStencilWrapper







! *******************************************************************************
!
! Purpose: Nullify 1D boundary-face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPatch              Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_NullifyBF2CStencil_1D(pRegion,pPatch)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_NullifyBF2CStencil_1D',&
  'RFLU_ModStencilsBFaces.F90')

! ******************************************************************************
!   Nullify memory 
! ******************************************************************************

    IF ( pPatch%bcType /= BC_VIRTUAL ) THEN 
      NULLIFY(pPatch%bf2cs1D)     
    END IF ! pPatch%bcType

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_NullifyBF2CStencil_1D


  


! *******************************************************************************
!
! Purpose: Nullify boundary-face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPatch              Pointer to patch
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_NullifyBF2CStencil(pRegion,pPatch)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_NullifyBF2CStencil',&
  'RFLU_ModStencilsBFaces.F90')

! ******************************************************************************
!   Nullify memory 
! ******************************************************************************

    IF ( pPatch%bcType /= BC_VIRTUAL ) THEN 
      NULLIFY(pPatch%bf2cs)     
    END IF ! pPatch%bcType

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_NullifyBF2CStencil     
  
  
  
  
  
  

! *******************************************************************************
!
! Purpose: Set boundary-face-to-cell stencil information.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPatch              Pointer to patch
!   orderNominalInput   Nominal order of accuracy
!
! Output: None.
!
! Notes: 
!   1. NOTE need to guard against orderInput being zero if running with 
!      first-order scheme.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetInfoBF2CStencil_1D(pRegion,pPatch,orderNominalInput)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: orderNominalInput
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: nBFaceMembsMax,nLayersMax,orderNominal,stencilSizeMax, &
               stencilSizeMin    
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_SetInfoBF2CStencil_1D',&
  'RFLU_ModStencilsBFaces.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
            'Setting 1D boundary-face-to-cell stencil information...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set stencil information. NOTE nBFaceMembsMax must be one less than the  
!   number of unknowns (or columns). NOTE orderNominal must be at least 2 
!   so can get stencil of at least 2 cells.
! ******************************************************************************

    orderNominal = MAX(orderNominalInput,2)

    nLayersMax     = orderNominal 
    nBFaceMembsMax = 0            ! TEMPORARY
    stencilSizeMin = orderNominal ! No difference between min and max value
    stencilSizeMax = orderNominal     
      
    pPatch%bf2cs1DInfo%orderNominal   = orderNominal
    pPatch%bf2cs1DInfo%nLayersMax     = nLayersMax      
    pPatch%bf2cs1DInfo%nBFaceMembsMax = nBFaceMembsMax 
    pPatch%bf2cs1DInfo%nCellMembsMax  = stencilSizeMax    
    pPatch%bf2cs1DInfo%nCellMembsMin  = stencilSizeMin      

! ******************************************************************************
!   Print stencil information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
            'Maximum allowed number of cell layers:  ',nLayersMax
      WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
            'Minimum required number of cell members:',stencilSizeMin                          
      WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
            'Maximum allowed number of cell members: ',stencilSizeMax   
    END IF ! global%myProcid
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
            'Setting 1D boundary-face-to-cell stencil information done.'
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_SetInfoBF2CStencil_1D  
  
  
  
  



! *******************************************************************************
!
! Purpose: Set boundary-face-to-cell stencil information.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPatch              Pointer to patch
!   orderNominalInput   Nominal order of accuracy
!
! Output: None.
!
! Notes: 
!   1. NOTE need to guard against orderInput being zero if running with 
!      first-order scheme.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetInfoBF2CStencil(pRegion,pPatch,orderNominalInput)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: orderNominalInput
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: nBFaceMembsMax,nLayersMax,orderNominal,stencilSizeMax, &
               stencilSizeMin    
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_SetInfoBF2CStencil',&
  'RFLU_ModStencilsBFaces.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
            'Setting boundary-face-to-cell stencil information...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set stencil information. NOTE nBFaceMembsMax must be one less than the  
!   number of unknowns (or columns), otherwise LAPACK routine used for 
!   constrained least-squares problem always gives trivial solution.
! ******************************************************************************

    orderNominal = MAX(orderNominalInput,1)

    nLayersMax     =  6 
    nBFaceMembsMax = 12
    stencilSizeMin = RFLU_ComputeStencilSize(global,pRegion%mixtInput%dimens, &
                                             1,orderNominal)
    stencilSizeMax = 10*stencilSizeMin       
      
    pPatch%bf2csInfo%orderNominal   = orderNominal
    pPatch%bf2csInfo%nLayersMax     = nLayersMax      
    pPatch%bf2csInfo%nBFaceMembsMax = nBFaceMembsMax 
    pPatch%bf2csInfo%nCellMembsMax  = stencilSizeMax    
    pPatch%bf2csInfo%nCellMembsMin  = stencilSizeMin      

! ******************************************************************************
!   Print stencil information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
            'Maximum allowed number of cell layers:  ',nLayersMax
      WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
            'Minimum required number of cell members:',stencilSizeMin                          
      WRITE(STDOUT,'(A,3X,A,1X,I3)') SOLVER_NAME, &
            'Maximum allowed number of cell members: ',stencilSizeMax   
    END IF ! global%myProcid
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
            'Setting boundary-face-to-cell stencil information done.'
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_SetInfoBF2CStencil








! *******************************************************************************
!
! Purpose: Wrapper routine for setting info for boundary face-to-cell stencils.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPatch              Pointer to patch
!   orderNominal        Nominal order of accuracy
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetInfoBF2CStencilWrapper(pRegion,pPatch,orderNominal)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: orderNominal
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================
      
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_SetInfoBF2CStencilWrapper',&
  'RFLU_ModStencilsBFaces.F90')

! ******************************************************************************
!   Call routines to set info for stencils
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%stencilDimensBFaces )
      CASE ( 1 ) 
        CALL RFLU_SetInfoBF2CStencil_1D(pRegion,pPatch,orderNominal)
      CASE ( 2,3 ) 
        CALL RFLU_SetInfoBF2CStencil(pRegion,pPatch,orderNominal)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%stencilDimensBFaces

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_SetInfoBF2CStencilWrapper 








! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModStencilsBFaces


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModStencilsBFaces.F90,v $
! Revision 1.10  2008/12/06 08:44:24  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:35  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2007/07/08 21:45:03  gzheng
! changed the PRESENT is used for PGI compiler
!
! Revision 1.7  2006/12/15 13:41:41  haselbac
! Changed stencil construction for sy patches so get symmetric stencils
!
! Revision 1.6  2006/12/15 13:26:36  haselbac
! Fixed bug in format statement, found by ifort
!
! Revision 1.5  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.4  2006/04/07 14:50:59  haselbac
! Added wrapper funcs, 1D stencil capability
!
! Revision 1.3  2005/10/27 19:19:35  haselbac
! Changed names, clean-up
!
! Revision 1.2  2005/10/18 03:00:56  haselbac
! Increased nBFaceMembsMax
!
! Revision 1.1  2005/10/05 14:33:44  haselbac
! Initial revision
!
! ******************************************************************************




















