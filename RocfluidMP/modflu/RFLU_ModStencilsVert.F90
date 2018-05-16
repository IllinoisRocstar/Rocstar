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
! Purpose: Suite of routines to construct vertex stencils.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModStencilsVert.F90,v 1.5 2008/12/06 08:44:24 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModStencilsVert

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
  PUBLIC :: RFLU_BuildStencilVert2Cell, &
            RFLU_CreateStencilVert2Cell, & 
            RFLU_DestroyStencilVert2Cell, &
            RFLU_NullifyStencilVert2Cell, &              
            RFLU_SetInfoStencilVert2Cell
            
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModStencilsVert.F90,v $ $Revision: 1.5 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS





! *******************************************************************************
!
! Purpose: Build vertex-to-cell stencil.
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

  SUBROUTINE RFLU_BuildStencilVert2Cell(pRegion)

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

    INTEGER :: degr,errorFlag,icg,icl,ict,iLayer,iloc,isl,ivg,ivl,iv2c, &
               nBFaceMembsMax,nBFaceMembsMaxTemp,nCellMembsInfoMax, &
               nCellMembsInfoMaxLoc,nCellMembsInfoMin,nCellMembsInfoMinLoc, &
               nLayersInfoMax,nLayersInfoMaxLoc,nLayersInfoMin, &
               nLayersInfoMinLoc,nLayersMax,nRows,order,orderNominal,sCount, &
               stencilSizeMax,stencilSizeMin,v2csBeg,v2csEnd
    INTEGER, DIMENSION(:), ALLOCATABLE :: v2cs 
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: layerInfo       
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: dr,wts   
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildStencilVert2Cell',&
  'RFLU_ModStencilsVert.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Building vertex-to-cell stencil...'            
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set grid pointer 
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Set variables
! ******************************************************************************

    orderNominal   = pGrid%v2csInfo%orderNominal 
    nLayersMax     = pGrid%v2csInfo%nLayersMax     
    nBFaceMembsMax = pGrid%v2csInfo%nBFaceMembsMax  
    stencilSizeMax = pGrid%v2csInfo%nCellMembsMax           
    stencilSizeMin = pGrid%v2csInfo%nCellMembsMin                     
                  
    nCellMembsInfoMax = 0
    nCellMembsInfoMin = HUGE(1)

    nLayersInfoMax = 0
    nLayersInfoMin = HUGE(1) 
    
    nBFaceMembsMaxTemp = 2*nBFaceMembsMax 
    
! ******************************************************************************
!   Allocate temporary memory
! ******************************************************************************

    ALLOCATE(v2cs(stencilSizeMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'v2cs')
    END IF ! global%error              

    ALLOCATE(layerInfo(X2CS_LAYER_BEG:X2CS_LAYER_END,nLayersMax), &
             STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'layerInfo')
    END IF ! global%error               

! ******************************************************************************
!   Loop over vertices
! ******************************************************************************

    DO ivg = 1,pGrid%nVertTot            

! ==============================================================================
!     Initialize
! ==============================================================================

      degr = 0

      DO isl = 1,stencilSizeMax
        v2cs(isl) = 0
      END DO ! isl

      DO iLayer = 1,nLayersMax
        layerInfo(X2CS_LAYER_BEG,iLayer) = 0
        layerInfo(X2CS_LAYER_END,iLayer) = 0          
      END DO ! iLayer    

      pGrid%v2cs(ivg)%nLayers = 1 

! ==============================================================================
!     Build basic stencil
! ==============================================================================

      DO iv2c = pGrid%v2cInfo(V2C_BEG,ivg),pGrid%v2cInfo(V2C_END,ivg)
        degr = degr + 1
        v2cs(degr) = pGrid%v2c(iv2c)
      END DO ! iv2c

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

        IF ( degr >= 4 ) THEN
          nRows = degr

          ALLOCATE(dr(XCOORD:ZCOORD,nRows),STAT=errorFlag)
          global%error = errorFlag   
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'dr')
          END IF ! global%error 

          ALLOCATE(wts(1,nRows),STAT=errorFlag)
          global%error = errorFlag   
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'wts')
          END IF ! global%error          

          DO isl = 1,degr
            icg = v2cs(isl)

            dr(XCOORD,isl) = pGrid%cofg(XCOORD,icg) - pGrid%xyz(XCOORD,ivg)
            dr(YCOORD,isl) = pGrid%cofg(YCOORD,icg) - pGrid%xyz(YCOORD,ivg)
            dr(ZCOORD,isl) = pGrid%cofg(ZCOORD,icg) - pGrid%xyz(ZCOORD,ivg) 
          END DO ! isl

          CALL RFLU_ComputeStencilWeights(global,pRegion%mixtInput%dimens, & 
                                          COMPWTS_MODE_FIXED, & 
                                          COMPWTS_SCAL_INVDIST,DERIV_DEGREE_0, &
                                          order,nRows,dr,wts,sCount)

          DEALLOCATE(dr,STAT=errorFlag)
          global%error = errorFlag   
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'dr')
          END IF ! global%error 

          DEALLOCATE(wts,STAT=errorFlag)
          global%error = errorFlag   
          IF ( global%error /= ERR_NONE ) THEN 
            CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'wts')
          END IF ! global%error 
        END IF ! degr

! ------------------------------------------------------------------------------
!       Check whether to reject or accept stencil. If singular or too small, 
!       add layer of cells. Pass 0 instead of ivg because want to prevent
!       present cell from being added, not present vertex. If pass present
!       vertex, could actually prevent a proper cell from being added in 
!       rare cases...
! ------------------------------------------------------------------------------

        IF ( sCount /= 0 .OR. degr <= stencilSizeMin ) THEN
          v2csBeg = layerInfo(X2CS_LAYER_BEG,iLayer-1)
          v2csEnd = layerInfo(X2CS_LAYER_END,iLayer-1)            

          CALL RFLU_AddCellLayer(global,pGrid,stencilSizeMax,0,degr, &
                                 v2csBeg,v2csEnd,v2cs)

          pGrid%v2cs(ivg)%nLayers = pGrid%v2cs(ivg)%nLayers + 1     

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

      pGrid%v2cs(ivg)%nCellMembs = degr

      ALLOCATE(pGrid%v2cs(ivg)%cellMembs(pGrid%v2cs(ivg)%nCellMembs), & 
               STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%v2cs%cellMembs')
      END IF ! global%error

      DO isl = 1,pGrid%v2cs(ivg)%nCellMembs
        pGrid%v2cs(ivg)%cellMembs(isl) = v2cs(isl)
      END DO ! isl  

      ALLOCATE(pGrid%v2cs(ivg)%layerInfo(X2CS_LAYER_BEG:X2CS_LAYER_END, &
               pGrid%v2cs(ivg)%nLayers),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%v2cs%layerInfo')
      END IF ! global%error        

      DO iLayer = 1,pGrid%v2cs(ivg)%nLayers
        pGrid%v2cs(ivg)%layerInfo(X2CS_LAYER_BEG,iLayer) = &
          layerInfo(X2CS_LAYER_BEG,iLayer)
        pGrid%v2cs(ivg)%layerInfo(X2CS_LAYER_END,iLayer) = &
          layerInfo(X2CS_LAYER_END,iLayer)            
      END DO ! iLayer   
      
! ==============================================================================
!     Extract information for later printing 
! ==============================================================================

      IF ( pGrid%v2cs(ivg)%nLayers < nLayersInfoMin ) THEN 
        nLayersInfoMin    = pGrid%v2cs(ivg)%nLayers
        nLayersInfoMinLoc = ivg
      END IF ! pGrid%v2cs(ivg)%nLayers                   

      IF ( pGrid%v2cs(ivg)%nLayers > nLayersInfoMax ) THEN 
        nLayersInfoMax    = pGrid%v2cs(ivg)%nLayers
        nLayersInfoMaxLoc = ivg
      END IF ! pGrid%v2cs(ivg)%nLayers

      IF ( pGrid%v2cs(ivg)%nCellMembs < nCellMembsInfoMin ) THEN 
        nCellMembsInfoMin    = pGrid%v2cs(ivg)%nCellMembs
        nCellMembsInfoMinLoc = ivg
      END IF ! pGrid%v2cs(ivg)%nCellMembs                   

      IF ( pGrid%v2cs(ivg)%nCellMembs > nCellMembsInfoMax ) THEN 
        nCellMembsInfoMax    = pGrid%v2cs(ivg)%nCellMembs
        nCellMembsInfoMaxLoc = ivg        
      END IF ! pGrid%v2cs(ivg)%nCellMembs             
    END DO ! ivg      

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

    DEALLOCATE(v2cs,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'v2cs')
    END IF ! global%error                   

    DEALLOCATE(layerInfo,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'layerInfo')
    END IF ! global%error            

#ifdef CHECK_DATASTRUCT
! ******************************************************************************
!   Data structure output for checking
! ******************************************************************************

    WRITE(STDOUT,'(A)') SOLVER_NAME
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Vertex-to-cell stencils'
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Maximum number of layers:', & 
                                   pGrid%v2csInfo%nLayersMax      
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Minimum stencil size:', & 
                                   pGrid%v2csInfo%nCellMembsMin   

    DO ivg = 1,pGrid%nVertTot
      WRITE(STDOUT,'(A,1X,I6,2(1X,I3),3X,20(1X,I6))') SOLVER_NAME,ivg, & 
            pGrid%v2cs(ivg)%nLayers,pGrid%v2cs(ivg)%nCellMembs, & 
            pGrid%v2cs(ivg)%cellMembs(1:pGrid%v2cs(ivg)%nCellMembs)
    END DO ! ivg

    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'    
    WRITE(STDOUT,'(A)') SOLVER_NAME
#endif          

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Building vertex-to-cell stencil done.'
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BuildStencilVert2Cell

     
     
  
  
  
  
! *******************************************************************************
!
! Purpose: Create vertex-to-cell stencil.
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

  SUBROUTINE RFLU_CreateStencilVert2Cell(pRegion)

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

    INTEGER :: errorFlag,ivg
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildStencilVert2Cell',&
  'RFLU_ModStencilsVert.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Creating vertex-to-cell stencil...'            
    END IF ! global%verbLevel

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyStencilVert2Cell(pRegion)

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Allocate memory and initialize
! ******************************************************************************

    ALLOCATE(pGrid%v2cs(pGrid%nVertTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%v2cs')
    END IF ! global%error               

    DO ivg = 1,pGrid%nVertTot
      pGrid%v2cs(ivg)%nCellMembs  = 0
      pGrid%v2cs(ivg)%nBFaceMembs = 0
    END DO ! ivg

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Creating vertex-to-cell stencil done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_CreateStencilVert2Cell    
 
 
 





! *******************************************************************************
!
! Purpose: Destroy vertex-to-cell stencil.
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

  SUBROUTINE RFLU_DestroyStencilVert2Cell(pRegion)

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

    INTEGER :: errorFlag,ivg
    TYPE(t_grid), POINTER :: pGrid        
    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyStencilVert2Cell',&
  'RFLU_ModStencilsVert.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Destroying vertex-to-cell stencil...'            
    END IF ! global%verbLevel 

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    DO ivg = 1,pGrid%nVertTot
      DEALLOCATE(pGrid%v2cs(ivg)%cellMembs,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%v2cs%cellMembs')
      END IF ! global%error          
    END DO ! ivg      

    DEALLOCATE(pGrid%v2cs,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%v2cs')
    END IF ! global%error

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyStencilVert2Cell(pRegion)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Destroying vertex-to-cell stencil done.'
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyStencilVert2Cell


     
     
  
  
! *******************************************************************************
!
! Purpose: Nullify vertex-to-cell stencil.
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

  SUBROUTINE RFLU_NullifyStencilVert2Cell(pRegion)

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

    CALL RegisterFunction(global,'RFLU_NullifyStencilVert2Cell',&
  'RFLU_ModStencilsVert.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Nullifying vertex-to-cell stencil...'            
    END IF ! global%verbLevel

! ******************************************************************************
!   Set grid pointer
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    NULLIFY(pGrid%v2cs)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Nullifying vertex-to-cell stencil done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_NullifyStencilVert2Cell     





! *******************************************************************************
!
! Purpose: Set vertex-to-cell stencil information.
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

  SUBROUTINE RFLU_SetInfoStencilVert2Cell(pRegion,orderNominal)

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

    CALL RegisterFunction(global,'RFLU_SetInfoStencilVert2Cell',&
  'RFLU_ModStencilsVert.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Setting vert-to-cell stencil information...'            
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
    nBFaceMembsMax = RFLU_ComputeStencilSize(global,pRegion%mixtInput%dimens, &
                                             1,orderNominal) - 1  
    stencilSizeMin = RFLU_ComputeStencilSize(global,pRegion%mixtInput%dimens, &
                                             2,orderNominal)
    stencilSizeMax = 10*stencilSizeMin       

    pGrid%v2csInfo%orderNominal   = orderNominal
    pGrid%v2csInfo%nLayersMax     = nLayersMax      
    pGrid%v2csInfo%nBFaceMembsMax = nBFaceMembsMax 
    pGrid%v2csInfo%nCellMembsMax  = stencilSizeMax    
    pGrid%v2csInfo%nCellMembsMin  = stencilSizeMin

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
                               'Setting vertex-to-cell stencil information done.'
    END IF ! global%verbLevel    

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_SetInfoStencilVert2Cell





! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModStencilsVert


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModStencilsVert.F90,v $
! Revision 1.5  2008/12/06 08:44:24  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:35  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/12/15 13:26:36  haselbac
! Fixed bug in format statement, found by ifort
!
! Revision 1.2  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.1  2005/10/05 14:33:44  haselbac
! Initial revision
!
! ******************************************************************************











