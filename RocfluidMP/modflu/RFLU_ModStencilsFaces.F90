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
! Purpose: Suite of routines to construct face stencils.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModStencilsFaces.F90,v 1.11 2008/12/06 08:44:24 mtcampbe Exp $
!
! Copyright: (c) 2005-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModStencilsFaces

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
  PUBLIC :: RFLU_BuildListCF2CStencil, & 
            RFLU_BuildF2CStencilWrapper, &
            RFLU_CreateF2CStencilWrapper, &
            RFLU_DestroyF2CStencilWrapper, &
            RFLU_DestroyListCF2CStencil, & 
            RFLU_SetInfoF2CStencilWrapper
            
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModStencilsFaces.F90,v $ $Revision: 1.11 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
 






! *******************************************************************************
!
! Purpose: Build 1D face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Restricted to hexahedra.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildF2CStencil_1D(pRegion)

    USE ModTools, ONLY: FloatEqual
 
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

    INTEGER :: fnDir,f2cs1DBeg,f2cs1DEnd,degr,errorFlag,ifg,iLayer,iloc,iPatch,isl, &
               nCellMembsInfoMax,nCellMembsInfoMaxLoc,nCellMembsInfoMin, &
               nCellMembsInfoMinLoc,nFaces,nLayersInfoMax,nLayersInfoMaxLoc, &
               nLayersInfoMin,nLayersInfoMinLoc,nLayersMax,stencilSizeMax, &
               stencilSizeMin
    INTEGER, DIMENSION(:), ALLOCATABLE :: f2cs1D
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: layerInfo
    REAL(RFREAL) :: nx,ny,nz,nm
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildF2CStencil_1D',&
  'RFLU_ModStencilsFaces.F90')

    IF ( (global%myProcid == MASTERPROC) .AND. &
         (global%verbLevel > VERBOSE_NONE) ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building 1D face-to-cell stencil...'           
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

    nLayersMax     = pGrid%f2csInfo%nLayersMax
    stencilSizeMax = pGrid%f2csInfo%nCellMembsMax
    stencilSizeMin = pGrid%f2csInfo%nCellMembsMin

    nCellMembsInfoMax = 0
    nCellMembsInfoMin = HUGE(1)

    nLayersInfoMax = 0
    nLayersInfoMin = HUGE(1)

! ******************************************************************************
!   Allocate temporary memory
! ******************************************************************************

    ALLOCATE(f2cs1D(stencilSizeMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'f2cs1D')
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

    DO ifg = 1,pGrid%nFaces

! ==============================================================================
!     Initialize variables
! ==============================================================================

      nx = ABS(pGrid%fn(XCOORD,ifg))
      ny = ABS(pGrid%fn(YCOORD,ifg))      
      nz = ABS(pGrid%fn(ZCOORD,ifg))
            
      nm = MAX(nx,ny,nz)       

      IF ( FloatEqual(nm,1.0_RFREAL,1.0E-6_RFREAL) .EQV. .TRUE. ) THEN  
        IF ( nx > ny .AND. nx > nz ) THEN 
          fnDir = XCOORD
        ELSE IF ( ny > nx .AND. ny > nz ) THEN 
          fnDir = YCOORD
        ELSE IF ( nz > nx .AND. nz > ny ) THEN 
          fnDir = ZCOORD
        END IF ! nx
      ELSE 
! TEMPORARY
        WRITE(*,*) 'ERROR - Face not aligned with coordinate axes!'
        STOP
! END TEMPORARY      
      END IF ! FloatEqual

      degr = 0

      DO isl = 1,stencilSizeMax
        f2cs1D(isl) = 0
      END DO ! isl

      DO iLayer = 1,nLayersMax
        layerInfo(X2CS_LAYER_BEG,iLayer) = 0
        layerInfo(X2CS_LAYER_END,iLayer) = 0          
      END DO ! iLayer 

! ==============================================================================             
!     Build basic stencil consisting of cells meeting at vertices of face
! ==============================================================================             

      degr = 2            
      
      f2cs1D(1) = pGrid%f2c(1,ifg)
      f2cs1D(2) = pGrid%f2c(2,ifg)

      pGrid%f2cs1D(ifg)%nLayers = 1        

      layerInfo(X2CS_LAYER_BEG,1) = 1
      layerInfo(X2CS_LAYER_END,1) = degr
   
! ==============================================================================             
!     Extend basic stencil. NOTE for 1D stencil do not have to check weight 
!     singularity
! ==============================================================================             

      DO iLayer = 2,nLayersMax
        IF ( degr < stencilSizeMin ) THEN 
          f2cs1DBeg = layerInfo(X2CS_LAYER_BEG,iLayer-1)
          f2cs1DEnd = layerInfo(X2CS_LAYER_END,iLayer-1)              

          CALL RFLU_AddCellLayer_1D(global,pGrid,stencilSizeMax,0,degr, &
                                    f2cs1DBeg,f2cs1DEnd,f2cs1D,fnDir)

          pGrid%f2cs1D(ifg)%nLayers = pGrid%f2cs1D(ifg)%nLayers + 1

          layerInfo(X2CS_LAYER_BEG,iLayer) = &
            layerInfo(X2CS_LAYER_END,iLayer-1) + 1 
          layerInfo(X2CS_LAYER_END,iLayer) = degr
        ELSE 
          EXIT
        END IF ! degr       
      END DO ! iLayer 

! ==============================================================================        
!     Store stencil and layer info
! ==============================================================================

      pGrid%f2cs1D(ifg)%nCellMembs = degr

      ALLOCATE(pGrid%f2cs1D(ifg)%cellMembs(degr),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%f2cs1D%cellMembs')
      END IF ! global%error

      DO isl = 1,degr
        pGrid%f2cs1D(ifg)%cellMembs(isl) = f2cs1D(isl)
      END DO ! isl  

      ALLOCATE(pGrid%f2cs1D(ifg)%layerInfo(X2CS_LAYER_BEG:X2CS_LAYER_END, &
               pGrid%f2cs1D(ifg)%nLayers),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%f2cs1D%layerInfo')
      END IF ! global%error        

      DO iLayer = 1,pGrid%f2cs1D(ifg)%nLayers
        pGrid%f2cs1D(ifg)%layerInfo(X2CS_LAYER_BEG,iLayer) = &
          layerInfo(X2CS_LAYER_BEG,iLayer)
        pGrid%f2cs1D(ifg)%layerInfo(X2CS_LAYER_END,iLayer) = &
          layerInfo(X2CS_LAYER_END,iLayer)            
      END DO ! iLayer 

! ==============================================================================
!     Extract information for later printing 
! ==============================================================================      
                      
      IF ( pGrid%f2cs1D(ifg)%nLayers < nLayersInfoMin ) THEN 
        nLayersInfoMin    = pGrid%f2cs1D(ifg)%nLayers
        nLayersInfoMinLoc = ifg
      END IF ! pGrid%f2cs(ifg)%nLayers                   

      IF ( pGrid%f2cs1D(ifg)%nLayers > nLayersInfoMax ) THEN 
        nLayersInfoMax    = pGrid%f2cs1D(ifg)%nLayers
        nLayersInfoMaxLoc = ifg
      END IF ! pGrid%f2cs(ifg)%nLayers                      
                        
      IF ( pGrid%f2cs1D(ifg)%nCellMembs < nCellMembsInfoMin ) THEN 
        nCellMembsInfoMin    = pGrid%f2cs1D(ifg)%nCellMembs
        nCellMembsInfoMinLoc = ifg
      END IF ! pGrid%f2cs(ifg)%nCellMembs                   

      IF ( pGrid%f2cs1D(ifg)%nCellMembs > nCellMembsInfoMax ) THEN 
        nCellMembsInfoMax    = pGrid%f2cs1D(ifg)%nCellMembs
        nCellMembsInfoMaxLoc = ifg        
      END IF ! pGrid%f2cs1D(ifg)%nCellMembs
    END DO ! ifg     

! ******************************************************************************
!   Write out information on stencils
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Statistics:' 
      WRITE(STDOUT,'(A,5X,A,2(1X,I3),2(1X,I9))') SOLVER_NAME, &
            'Minimum/maximum number of cell layers: ',nLayersInfoMin, &
            nLayersInfoMax,nLayersInfoMinLoc,nLayersInfoMaxLoc               
      WRITE(STDOUT,'(A,5X,A,2(1X,I3),2(1X,I9))') SOLVER_NAME, &
            'Minimum/maximum number of cell members:',nCellMembsInfoMin, &
            nCellMembsInfoMax,nCellMembsInfoMinLoc,nCellMembsInfoMaxLoc       
    END IF ! global%myProcid
    
! ******************************************************************************
!   Deallocate temporary memory
! ******************************************************************************

    DEALLOCATE(f2cs1D,STAT=errorFlag)
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
         (global%verbLevel > VERBOSE_NONE) ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Building 1D face-to-cell stencil done.'            
    END IF ! global%myProcid   

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildF2CStencil_1D







! *******************************************************************************
!
! Purpose: Build face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   addBFaces           Flag indicating whether should add boundary faces
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildF2CStencil(pRegion,addBFaces)

    USE RFLU_ModFaceList, ONLY: RFLU_GetOpposingFaces
    
    USE ModInterfaces, ONLY: RFLU_InvertMatrixSVD

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    LOGICAL, INTENT(IN) :: addBFaces
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: degr,errorFlag,f2csBeg,f2csEnd,icg,ifg,ifg2,ifl,iLayer,iloc, &
               iPatch,isl,ivl,iv2c,nBFaceMembs,nBFaceMembsMax, &
               nBFaceMembsMaxTemp,nCellMembsInfoMax,nCellMembsInfoMaxLoc, &
               nCellMembsInfoMin,nCellMembsInfoMinLoc,nFacesOpp, &
               nLayersInfoMax,nLayersInfoMaxLoc,nLayersInfoMin, &
               nLayersInfoMinLoc,nLayersMax,nRows,order,orderNominal,sCount, &
               stencilSizeMax,stencilSizeMin,nCols,iRow,iCol
    INTEGER :: bf2v(4)
    INTEGER :: faceOppInfo(2,2)
    INTEGER, DIMENSION(:), ALLOCATABLE :: f2cs
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: bFaceMembs,layerInfo 
    REAL(RFREAL) :: dx,dy,dz,term
    REAL(RFREAL) :: colMax(4)     
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: a,aInv 
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_patch), POINTER :: pPatch      
    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildF2CStencil',&
  'RFLU_ModStencilsFaces.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building face-to-cell stencil...'
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set grid pointer 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Get variables
! ******************************************************************************

    orderNominal   = pGrid%f2csInfo%orderNominal 
    nLayersMax     = pGrid%f2csInfo%nLayersMax     
    nBFaceMembsMax = pGrid%f2csInfo%nBFaceMembsMax  
    stencilSizeMax = pGrid%f2csInfo%nCellMembsMax           
    stencilSizeMin = pGrid%f2csInfo%nCellMembsMin                     

    nCellMembsInfoMax = 0
    nCellMembsInfoMin = HUGE(1)
    
    nLayersInfoMax = 0
    nLayersInfoMin = HUGE(1)    
    
    nBFaceMembsMaxTemp = 2*nBFaceMembsMax
                    
    IF ( (global%myProcid == MASTERPROC) .AND. &
         (global%verbLevel > VERBOSE_LOW) ) THEN
      WRITE(STDOUT,'(A,3X,A,1X,L1)') SOLVER_NAME,'Adding boundary faces:', &
                                     addBFaces            
    END IF ! global%verbLevel      
    
! ******************************************************************************
!   Allocate temporary memory
! ******************************************************************************

    ALLOCATE(f2cs(stencilSizeMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'f2cs')
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
!   Loop over interior faces
! ******************************************************************************

    DO ifg = 1,pGrid%nFaces

! ==============================================================================                       
!    Initialize
! ==============================================================================             

      DO isl = 1,stencilSizeMax
        f2cs(isl) = 0
      END DO ! isl           

      DO iLayer = 1,nLayersMax
        layerInfo(X2CS_LAYER_BEG,iLayer) = 0
        layerInfo(X2CS_LAYER_END,iLayer) = 0          
      END DO ! iLayer                 

! ==============================================================================             
!     Build basic stencil consisting of cells meeting at vertices of face
! ==============================================================================             

      degr = 0            
      
      CALL RFLU_AddFaceVertNeighbs(global,pGrid,stencilSizeMax, & 
                                   pGrid%f2v(1:4,ifg),degr,f2cs)

      pGrid%f2cs(ifg)%nLayers = 1            

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
                icg = f2cs(isl)

                dx = pGrid%cofg(XCOORD,icg) - pGrid%fc(XCOORD,ifg)
                dy = pGrid%cofg(YCOORD,icg) - pGrid%fc(YCOORD,ifg)
                
                term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)
                
                a(isl,1) = term
                a(isl,2) = term*dx
                a(isl,3) = term*dy                                
              END DO ! isl             
            CASE ( 3 )
              DO isl = 1,degr
                icg = f2cs(isl)

                dx = pGrid%cofg(XCOORD,icg) - pGrid%fc(XCOORD,ifg)
                dy = pGrid%cofg(YCOORD,icg) - pGrid%fc(YCOORD,ifg)
                dz = pGrid%cofg(ZCOORD,icg) - pGrid%fc(ZCOORD,ifg)
                
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
!       add layer of cells. Pass 0 instead of ifg because want to prevent
!       present cell from being added, not present face. If pass present
!       face, could actually prevent a proper cell from being added in 
!       rare cases...
! ------------------------------------------------------------------------------

        IF ( sCount /= 0 .OR. degr < stencilSizeMin ) THEN             
          f2csBeg = layerInfo(X2CS_LAYER_BEG,iLayer-1)
          f2csEnd = layerInfo(X2CS_LAYER_END,iLayer-1)              

          CALL RFLU_AddCellLayer(global,pGrid,stencilSizeMax,0,degr, &
                                 f2csBeg,f2csEnd,f2cs)

          pGrid%f2cs(ifg)%nLayers = pGrid%f2cs(ifg)%nLayers + 1

          layerInfo(X2CS_LAYER_BEG,iLayer) = & 
            layerInfo(X2CS_LAYER_END,iLayer-1) + 1 
          layerInfo(X2CS_LAYER_END,iLayer) = degr 
        ELSE 
          EXIT
        END IF ! sCount       
      END DO ! iLayer        

! ==============================================================================        
!     Store stencil and layer info
! ==============================================================================

      pGrid%f2cs(ifg)%nCellMembs = degr

      ALLOCATE(pGrid%f2cs(ifg)%cellMembs(degr),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%f2cs%cellMembs')
      END IF ! global%error

      DO isl = 1,degr
        pGrid%f2cs(ifg)%cellMembs(isl) = f2cs(isl)
      END DO ! isl  

      ALLOCATE(pGrid%f2cs(ifg)%layerInfo(X2CS_LAYER_BEG:X2CS_LAYER_END, &
               pGrid%f2cs(ifg)%nLayers),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%f2cs%layerInfo')
      END IF ! global%error        

      DO iLayer = 1,pGrid%f2cs(ifg)%nLayers
        pGrid%f2cs(ifg)%layerInfo(X2CS_LAYER_BEG,iLayer) = &
          layerInfo(X2CS_LAYER_BEG,iLayer)
        pGrid%f2cs(ifg)%layerInfo(X2CS_LAYER_END,iLayer) = &
          layerInfo(X2CS_LAYER_END,iLayer)            
      END DO ! iLayer 

! ==============================================================================        
!     Add boundary faces to stencil. If the stencil contains boundary faces,
!     sort them by distance and pick <nBFaceMembsMax> closest ones.
! ==============================================================================        

      nBFaceMembs = 0

      IF ( addBFaces .EQV. .TRUE. ) THEN 
        CALL RFLU_AddBFaces(pRegion,nBFaceMembsMaxTemp, & 
                            pGrid%f2cs(ifg)%nCellMembs, & 
                            pGrid%f2cs(ifg)%cellMembs(1:pGrid%f2cs(ifg)%nCellMembs), &
                            nBFaceMembs,bFaceMembs)
      END IF ! addBFaces

      IF ( nBFaceMembs > 0 ) THEN 
        CALL RFLU_SortBFaces(pRegion,pGrid%fc(XCOORD:ZCOORD,icg), &
                             nBFaceMembs,bFaceMembs(1:2,1:nBFaceMembs)) 

        pGrid%f2cs(ifg)%nBFaceMembs = MIN(nBFaceMembs,nBFaceMembsMax)

        ALLOCATE(pGrid%f2cs(ifg)%bFaceMembs(2,pGrid%f2cs(ifg)%nBFaceMembs), & 
                 STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%f2cs%bFaceMembs')
        END IF ! global%error

        DO isl = 1,pGrid%f2cs(ifg)%nBFaceMembs
          pGrid%f2cs(ifg)%bFaceMembs(1,isl) = bFaceMembs(1,isl)
          pGrid%f2cs(ifg)%bFaceMembs(2,isl) = bFaceMembs(2,isl)            
        END DO ! isl
      ELSE 
        pGrid%f2cs(ifg)%nBFaceMembs = 0
        
        NULLIFY(pGrid%f2cs(ifg)%bFaceMembs)
      END IF ! nBFaceMembs 
      
! ==============================================================================
!     Extract information for later printing 
! ==============================================================================      
                      
      IF ( pGrid%f2cs(ifg)%nLayers < nLayersInfoMin ) THEN 
        nLayersInfoMin    = pGrid%f2cs(ifg)%nLayers
        nLayersInfoMinLoc = ifg
      END IF ! pGrid%f2cs(ifg)%nLayers                   

      IF ( pGrid%f2cs(ifg)%nLayers > nLayersInfoMax ) THEN 
        nLayersInfoMax    = pGrid%f2cs(ifg)%nLayers
        nLayersInfoMaxLoc = ifg
      END IF ! pGrid%f2cs(ifg)%nLayers                      
                        
      IF ( pGrid%f2cs(ifg)%nCellMembs < nCellMembsInfoMin ) THEN 
        nCellMembsInfoMin    = pGrid%f2cs(ifg)%nCellMembs
        nCellMembsInfoMinLoc = ifg
      END IF ! pGrid%f2cs(ifg)%nCellMembs                   

      IF ( pGrid%f2cs(ifg)%nCellMembs > nCellMembsInfoMax ) THEN 
        nCellMembsInfoMax    = pGrid%f2cs(ifg)%nCellMembs
        nCellMembsInfoMaxLoc = ifg        
      END IF ! pGrid%f2cs(ifg)%nCellMembs
    END DO ! ifg             

! ******************************************************************************
!   Write out information on stencils
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_LOW ) THEN 
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
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Face-to-cell stencils'
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Maximum number of layers:', & 
                                   pGrid%f2csInfo%nLayersMax        
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Minimum stencil size:', & 
                                   pGrid%f2csInfo%nCellMembsMin   
    DO ifg = 1,pGrid%nFaces
      WRITE(STDOUT,'(A,1X,I6,2(1X,I3),3X,20(1X,I6))') SOLVER_NAME,ifg, & 
            pGrid%f2cs(ifg)%nLayers,pGrid%f2cs(ifg)%nCellMembs, &
            pGrid%f2cs(ifg)%cellMembs(1:pGrid%f2cs(ifg)%nCellMembs)
    END DO ! ifg

    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'    
    WRITE(STDOUT,'(A)') SOLVER_NAME
#endif          

! ******************************************************************************
!   Deallocate temporary memory
! ******************************************************************************

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

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building face-to-cell stencil done.'            
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildF2CStencil

  
  
  
  

! *******************************************************************************
!
! Purpose: Wrapper for building face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   constrInput         Flag indicating whether have constrained reconstruction
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildF2CStencilWrapper(pRegion,constrInput)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN), OPTIONAL :: constrInput
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

    CALL RegisterFunction(global,'RFLU_BuildF2CStencilWrapper',&
  'RFLU_ModStencilsFaces.F90')

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

    SELECT CASE ( pRegion%mixtInput%stencilDimensFaces )
      CASE ( 1 ) 
        CALL RFLU_BuildF2CStencil_1D(pRegion)
      CASE ( 2,3 ) 
        CALL RFLU_BuildF2CStencil(pRegion,addBFaces)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%stencilDimensFaces
      

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildF2CStencilWrapper

  
  
  
  

! *******************************************************************************
!
! Purpose: Build list of face-to-cell stencils which are constrained.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildListCF2CStencil(pRegion)

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

    INTEGER :: errorFlag,ifg
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildListCF2CStencil',&
  'RFLU_ModStencilsFaces.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Building list of constrained ', & 
                                 'face-to-cell stencil...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Count and build list of constrained cell-to-cell stencils
! ******************************************************************************

    pGrid%nFacesConstr = 0

    IF ( pRegion%mixtInput%cReconstFaces > CONSTR_NONE ) THEN 
      DO ifg = 1,pGrid%nFaces
        IF ( pGrid%f2cs(ifg)%nBFaceMembs > 0 ) THEN
          pGrid%nFacesConstr = pGrid%nFacesConstr + 1
        END IF ! pGrid%f2cs(ifg)%nBFaceMembs
      END DO ! ifg
      
      IF ( pGrid%nFacesConstr > 0 ) THEN
        ALLOCATE(pGrid%ifgConstr(pGrid%nFacesConstr),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%ifgConstr')
        END IF ! global%error
        
        pGrid%nFacesConstr = 0
        
        DO ifg = 1,pGrid%nFaces
          IF ( pGrid%f2cs(ifg)%nBFaceMembs > 0 ) THEN
            pGrid%nFacesConstr = pGrid%nFacesConstr + 1

            pGrid%ifgConstr(pGrid%nFacesConstr) = ifg
          END IF ! pGrid%f2cs(ifg)%nBFaceMembs
        END DO ! ifg
      ELSE 
        NULLIFY(pGrid%ifgConstr)        
      END IF ! pGrid%nFacesConstr
    ELSE 
      NULLIFY(pGrid%ifgConstr)
    END IF ! pRegion%mixtInput%cReconstFaces
        
! ******************************************************************************
!   Print info
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,3X,A,A,1X,I5)') SOLVER_NAME,'Number of constrained ', & 
                                       'face-to-cell stencils:', &
                                       pGrid%nFacesConstr           
    END IF ! global%verbLevel     
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Building list of constrained ', & 
                                 'face-to-cell stencil done.'   
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildListCF2CStencil   
  
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Create 1D face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_CreateF2CStencil_1D(pRegion)

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

    INTEGER :: errorFlag,ifg
    TYPE(t_grid), POINTER :: pGrid       
    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateF2CStencil_1D',&
  'RFLU_ModStencilsFaces.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating 1D face-to-cell stencil...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyF2CStencil_1D(pRegion)

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory and initialize
! ******************************************************************************

    ALLOCATE(pGrid%f2cs1D(pGrid%nFaces),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%f2cs1D')
    END IF ! global%error               

    DO ifg = 1,pGrid%nFaces
      pGrid%f2cs1D(ifg)%nCellMembs  = 0
      pGrid%f2cs1D(ifg)%nBFaceMembs = 0
    END DO ! ifg

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Creating 1D face-to-cell stencil done.'            
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_CreateF2CStencil_1D  
  
  
  
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Create face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_CreateF2CStencil(pRegion)

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

    INTEGER :: errorFlag,ifg
    TYPE(t_grid), POINTER :: pGrid       
    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateF2CStencil',&
  'RFLU_ModStencilsFaces.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating face-to-cell stencil...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyF2CStencil(pRegion)

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory and initialize
! ******************************************************************************

    ALLOCATE(pGrid%f2cs(pGrid%nFaces),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%f2cs')
    END IF ! global%error               

    DO ifg = 1,pGrid%nFaces
      pGrid%f2cs(ifg)%nCellMembs  = 0
      pGrid%f2cs(ifg)%nBFaceMembs = 0
    END DO ! ifg

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Creating face-to-cell stencil done.'            
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_CreateF2CStencil    
     






! *******************************************************************************
!
! Purpose: Wrapper routine for creating face-to-cell stencils.
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

  SUBROUTINE RFLU_CreateF2CStencilWrapper(pRegion)

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

    CALL RegisterFunction(global,'RFLU_CreateF2CStencilWrapper',&
  'RFLU_ModStencilsFaces.F90')

! ******************************************************************************
!   Call routines to create stencils
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%stencilDimensFaces )
      CASE ( 1 )     
        CALL RFLU_CreateF2CStencil_1D(pRegion) 
      CASE ( 2,3 ) 
        CALL RFLU_CreateF2CStencil(pRegion)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%stencilDimensFaces

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_CreateF2CStencilWrapper



  




! *******************************************************************************
!
! Purpose: Destroy 1D face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_DestroyF2CStencil_1D(pRegion)

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

    INTEGER :: errorFlag,ifg
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global 

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyF2CStencil_1D',&
  'RFLU_ModStencilsFaces.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Destroying 1D face-to-cell stencil...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    DO ifg = 1,pGrid%nFaces
      DEALLOCATE(pGrid%f2cs1D(ifg)%cellMembs,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%f2cs1D%cellMembs')
      END IF ! global%error 

      IF ( pGrid%f2cs1D(ifg)%nBFaceMembs > 0 ) THEN        
        DEALLOCATE(pGrid%f2cs1D(ifg)%bFaceMembs,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__, & 
                         'pGrid%f2cs1D%bFaceMembs')
        END IF ! global%error 
      END IF ! pGrid%f2cs1D%nBFaceMembs                       
    END DO ! ifg

    DEALLOCATE(pGrid%f2cs1D,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%f2cs1D')
    END IF ! global%error
       
! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyF2CStencil_1D(pRegion)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Destroying 1D face-to-cell stencil done.'
    END IF ! global%verbLevel   

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyF2CStencil_1D






! *******************************************************************************
!
! Purpose: Destroy face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_DestroyF2CStencil(pRegion)

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

    INTEGER :: errorFlag,ifg
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global 

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyF2CStencil',&
  'RFLU_ModStencilsFaces.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Destroying face-to-cell stencil...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    DO ifg = 1,pGrid%nFaces
      DEALLOCATE(pGrid%f2cs(ifg)%cellMembs,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%f2cs%cellMembs')
      END IF ! global%error 

      IF ( pGrid%f2cs(ifg)%nBFaceMembs > 0 ) THEN        
        DEALLOCATE(pGrid%f2cs(ifg)%bFaceMembs,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__, & 
                         'pGrid%f2cs%bFaceMembs')
        END IF ! global%error 
      END IF ! pGrid%f2cs%nBFaceMembs                       
    END DO ! ifg

    DEALLOCATE(pGrid%f2cs,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%f2cs')
    END IF ! global%error
       
! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyF2CStencil(pRegion)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Destroying face-to-cell stencil done.'
    END IF ! global%verbLevel   

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyF2CStencil 
  
    
  
  
  
  

! *******************************************************************************
!
! Purpose: Wrapper routine for destroying face-to-cell stencils.
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

  SUBROUTINE RFLU_DestroyF2CStencilWrapper(pRegion)

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

    CALL RegisterFunction(global,'RFLU_DestroyF2CStencilWrapper',&
  'RFLU_ModStencilsFaces.F90')

! ******************************************************************************
!   Call routines to destroy stencils
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%stencilDimensFaces )
      CASE ( 1 )    
        CALL RFLU_DestroyF2CStencil_1D(pRegion)
      CASE ( 2,3 ) 
        CALL RFLU_DestroyF2CStencil(pRegion)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%stencilDimensFaces

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyF2CStencilWrapper    
  
  
  



! *******************************************************************************
!
! Purpose: Destroy list of face-to-cell stencils which are constrained.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_DestroyListCF2CStencil(pRegion)

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

    CALL RegisterFunction(global,'RFLU_DestroyListCF2CStencil',&
  'RFLU_ModStencilsFaces.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Destroying list of ', & 
                                 'constrained face-to-cell stencil...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Destroy list of constrained cell-to-cell stencils
! ******************************************************************************

    IF ( pRegion%mixtInput%cReconstFaces > CONSTR_NONE ) THEN 
      IF ( pGrid%nFacesConstr > 0 ) THEN
        DEALLOCATE(pGrid%ifgConstr,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%ifgConstr')
        END IF ! global%error
        
        pGrid%nFacesConstr = 0    
      END IF ! pGrid%nFacesConstr
    END IF ! pRegion%mixtInput%cReconstFaces
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Destroying list of ', & 
                                 'constrained face-to-cell stencil done.'   
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyListCF2CStencil 





! *******************************************************************************
!
! Purpose: Nullify face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_NullifyF2CStencil_1D(pRegion)

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

    INTEGER :: iPatch
    TYPE(t_grid), POINTER :: pGrid 
    TYPE(t_patch), POINTER :: pPatch       
    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_NullifyF2CStencil_1D',&
  'RFLU_ModStencilsFaces.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Nullifying 1D face-to-cell stencil...'            
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory 
! ******************************************************************************

    NULLIFY(pGrid%f2cs1D)         

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Nullifying 1D face-to-cell stencil done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_NullifyF2CStencil_1D


  
  
  
! *******************************************************************************
!
! Purpose: Nullify face-to-cell stencil.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_NullifyF2CStencil(pRegion)

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

    INTEGER :: iPatch
    TYPE(t_grid), POINTER :: pGrid 
    TYPE(t_patch), POINTER :: pPatch       
    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_NullifyF2CStencil',&
  'RFLU_ModStencilsFaces.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Nullifying face-to-cell stencil...'            
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory 
! ******************************************************************************

    NULLIFY(pGrid%f2cs)         

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Nullifying face-to-cell stencil done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_NullifyF2CStencil    
     






! *******************************************************************************
!
! Purpose: Set 1D face-to-cell stencil information.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   orderNominalInput   Nominal order of accuracy
!
! Output: None.
!
! Notes: 
!   1. NOTE need to guard against orderInput being zero if running with 
!      first-order scheme.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetInfoF2CStencil_1D(pRegion,orderNominalInput)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: orderNominalInput
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: nBFaceMembsMax,nLayersMax,orderNominal,stencilSizeMax, &
               stencilSizeMin
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_SetInfoF2CStencil_1D',&
  'RFLU_ModStencilsFaces.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Setting 1D face-to-cell stencil information...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

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

    pGrid%f2csInfo%orderNominal   = orderNominal
    pGrid%f2csInfo%nLayersMax     = nLayersMax      
    pGrid%f2csInfo%nBFaceMembsMax = nBFaceMembsMax 
    pGrid%f2csInfo%nCellMembsMax  = stencilSizeMax    
    pGrid%f2csInfo%nCellMembsMin  = stencilSizeMin

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
            'Setting 1D face-to-cell stencil information done.'
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_SetInfoF2CStencil_1D






! *******************************************************************************
!
! Purpose: Set face-to-cell stencil information.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   orderNominalInput   Nominal order of accuracy
!
! Output: None.
!
! Notes: 
!   1. NOTE need to guard against orderInput being zero if running with 
!      first-order scheme.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetInfoF2CStencil(pRegion,orderNominalInput)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: orderNominalInput
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: nBFaceMembsMax,nLayersMax,orderNominal,stencilSizeMax, &
               stencilSizeMin
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_SetInfoF2CStencil',&
  'RFLU_ModStencilsFaces.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Setting face-to-cell stencil information...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Set stencil information. NOTE nBFaceMembsMax must be one less than the  
!   number of unknowns (or columns), otherwise LAPACK routine used for 
!   constrained least-squares problem always gives trivial solution.
! ******************************************************************************

    orderNominal = MAX(orderNominalInput,1)

    nLayersMax     = 6
!    nBFaceMembsMax = RFLU_ComputeStencilSize(global,pRegion%mixtInput%dimens, &
!                                             1,orderNominal) - 1  
    nBFaceMembsMax = 3 
    stencilSizeMin = RFLU_ComputeStencilSize(global,pRegion%mixtInput%dimens, &
                                             1,orderNominal)
    stencilSizeMax = 10*stencilSizeMin       

    pGrid%f2csInfo%orderNominal   = orderNominal
    pGrid%f2csInfo%nLayersMax     = nLayersMax      
    pGrid%f2csInfo%nBFaceMembsMax = nBFaceMembsMax 
    pGrid%f2csInfo%nCellMembsMax  = stencilSizeMax    
    pGrid%f2csInfo%nCellMembsMin  = stencilSizeMin

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
                               'Setting face-to-cell stencil information done.'
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_SetInfoF2CStencil







! *******************************************************************************
!
! Purpose: Wrapper routine for setting info for face-to-cell stencils.
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

  SUBROUTINE RFLU_SetInfoF2CStencilWrapper(pRegion,orderNominal)

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

    CALL RegisterFunction(global,'RFLU_SetInfoF2CStencilWrapper',&
  'RFLU_ModStencilsFaces.F90')

! ******************************************************************************
!   Call routines to set info for stencils
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%stencilDimensFaces )
      CASE ( 1 ) 
        CALL RFLU_SetInfoF2CStencil_1D(pRegion,orderNominal)
      CASE ( 2,3 ) 
        CALL RFLU_SetInfoF2CStencil(pRegion,orderNominal)
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%stencilDimensFaces

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_SetInfoF2CStencilWrapper 






! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModStencilsFaces


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModStencilsFaces.F90,v $
! Revision 1.11  2008/12/06 08:44:24  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:35  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2007/07/08 21:45:03  gzheng
! changed the PRESENT is used for PGI compiler
!
! Revision 1.8  2006/12/15 13:26:36  haselbac
! Fixed bug in format statement, found by ifort
!
! Revision 1.7  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.6  2006/04/07 14:52:04  haselbac
! Adapted to new stencilDimens param, bug fixes for 1D stencils
!
! Revision 1.5  2006/03/09 15:04:51  haselbac
! Bug fix and put routines in right order
!
! Revision 1.4  2006/03/09 14:08:59  haselbac
! Wrapperified module bcos of 1D routines, removed CF2C list routines
!
! Revision 1.3  2005/12/25 15:33:58  haselbac
! Added face-specific constraint flag
!
! Revision 1.2  2005/10/27 19:19:36  haselbac
! Changed names, clean-up
!
! Revision 1.1  2005/10/05 14:33:44  haselbac
! Initial revision
!
! ******************************************************************************






















