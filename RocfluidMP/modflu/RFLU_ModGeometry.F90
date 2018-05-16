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
! Purpose: Suite of routines to compute geometry.
!
! Description: None.
!
! Notes: 
!   1. Need special destruction modes because cofg needed for interpolation from
!      cell centers to vertices, so once have the entire geometry constructed, 
!      can deallocate everything except cofg.
!   2. Initialize geometry also in creation routine, because cofg is used in
!      RFLU_PrintLocInfo, and although it is allocated before being used, the
!      fact that it is not initialized appears to cause a floating-point
!      exception (INVALID) on modi4. 
!
! ******************************************************************************
!
! $Id: RFLU_ModGeometry.F90,v 1.35 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModGeometry

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  IMPLICIT NONE
   
  PRIVATE
  PUBLIC :: RFLU_BuildBVertexNormals, & 
            RFLU_BuildGeometry, &             
            RFLU_ComputeApproxCentroids, & 
            RFLU_ComputeFaceDist, &             
            RFLU_CreateApproxCentroids, & 
            RFLU_CreateFaceDist, & 
            RFLU_CreateGeometry, & 
            RFLU_DestroyApproxCentroids, &
            RFLU_DestroyFaceDist, &  
            RFLU_DestroyGeometry, & 
            RFLU_NullifyApproxCentroids, &
            RFLU_NullifyGeometry 
                                   
  SAVE    
     
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
     
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModGeometry.F90,v $ $Revision: 1.35 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
 





! ******************************************************************************
!
! Purpose: Build boundary-vertex normals list 
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
    
  SUBROUTINE RFLU_BuildBVertexNormals(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: errorFlag,ibv,ic,ifc,iPatch
    REAL(RFREAL) :: term
    REAL(RFREAL), DIMENSION(:,:), POINTER :: pBvn
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildBVertexNormals',&
  'RFLU_ModGeometry.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building boundary-vertex '// & 
                               'normals...' 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

    pGrid => pRegion%grid  

! ******************************************************************************
!   Build list of boundary vertex normals for each patch
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)                    
      pBvn   => pPatch%bvn

! ==============================================================================
!     Accumulate face normals to vertices
! ==============================================================================

      DO ifc = 1,pPatch%nBFacesTot
        DO ic = 1,4
          ibv = pPatch%bf2v(ic,ifc)

          IF ( ibv /= VERT_NONE ) THEN             
            pBvn(XCOORD,ibv) = pBvn(XCOORD,ibv) + pPatch%fn(XCOORD,ifc)
            pBvn(YCOORD,ibv) = pBvn(YCOORD,ibv) + pPatch%fn(YCOORD,ifc)
            pBvn(ZCOORD,ibv) = pBvn(ZCOORD,ibv) + pPatch%fn(ZCOORD,ifc) 
          END IF ! ibv
        END DO ! ic         
      END DO ! ifc

! ==============================================================================
!     Normalize accumulated face normals
! ==============================================================================

      DO ibv = 1,pPatch%nBVertTot        
        term = 1.0_RFREAL/(SQRT(pBvn(XCOORD,ibv)*pBvn(XCOORD,ibv) & 
                              + pBvn(YCOORD,ibv)*pBvn(YCOORD,ibv) &
                              + pBvn(ZCOORD,ibv)*pBvn(ZCOORD,ibv)))        

        pBvn(XCOORD,ibv) = term*pBvn(XCOORD,ibv)
        pBvn(YCOORD,ibv) = term*pBvn(YCOORD,ibv)
        pBvn(ZCOORD,ibv) = term*pBvn(ZCOORD,ibv)
      END DO ! ibv               

      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        IF ( pPatch%nBVert > 0 ) THEN 
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Normal component extrema:'         
          WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'Patch:',iPatch             
          WRITE(STDOUT,'(A,7X,A,2(1X,E23.16))') SOLVER_NAME,'x-direction:', & 
                MINVAL(pBvn(XCOORD,1:pPatch%nBVert)), & 
                MAXVAL(pBvn(XCOORD,1:pPatch%nBVert))
          WRITE(STDOUT,'(A,7X,A,2(1X,E23.16))') SOLVER_NAME,'y-direction:', & 
                MINVAL(pBvn(YCOORD,1:pPatch%nBVert)), & 
                MAXVAL(pBvn(YCOORD,1:pPatch%nBVert)) 
          WRITE(STDOUT,'(A,7X,A,2(1X,E23.16))') SOLVER_NAME,'z-direction:', & 
                MINVAL(pBvn(ZCOORD,1:pPatch%nBVert)), & 
                MAXVAL(pBvn(ZCOORD,1:pPatch%nBVert))  
        END IF ! pPatch%nBVert                               
      END IF ! global%verbLevel      
    END DO ! iPatch

#ifdef CHECK_DATASTRUCT
! ==============================================================================
!   Data structure output for checking
! ==============================================================================

    IF ( pGrid%nPatches > 0 ) THEN
      WRITE(STDOUT,'(A)') SOLVER_NAME
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Boundary vertex normals:'

      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch) 
        WRITE(STDOUT,'(A,1X,A,1X,I3,3X,A)') SOLVER_NAME,'Patch:', & 
                                            iPatch,pPatch%bcName
        WRITE(STDOUT,'(A,1X,A,1X,A,1X,I7))') SOLVER_NAME,'Number of actual', & 
                                             'vertices:',pPatch%nBVert
        WRITE(STDOUT,'(A,1X,A,1X,A,1X,I7))') SOLVER_NAME,'Number of total', & 
                                             'vertices:',pPatch%nBVertTot

        DO ibv = 1,pPatch%nBVertTot
          WRITE(STDOUT,'(A,1X,I7,3(1X,E18.9))') SOLVER_NAME,ibv, & 
                                                pPatch%bvn(XCOORD:ZCOORD,ibv)
        END DO ! ibv

      END DO ! iPatch
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
      WRITE(STDOUT,'(A)') SOLVER_NAME
    END IF ! pGrid%nPatches       
#endif

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building boundary-vertex '// & 
                               'normals done.'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME 
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_BuildBVertexNormals








! ******************************************************************************
!
! Purpose: Build geometry - compute face areas and centroids and cell volumes.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   sypeFaceFlag        Flag indicating whether faces on symmetry and periodic
!                       patches should contribute to computation of volumes 
!                       of cells adjacent to patches
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_BuildGeometry(pRegion,sypeFaceFlag)

    USE ModSortSearch

    USE RFLU_ModCellFaceEdgeInfo     

    USE ModInterfaces, ONLY: FaceCentroidQuad, & 
                             FaceCentroidTria, & 
                             FaceVectorQuad, & 
                             FaceVectorTria, &
                             RFLU_PrintLocInfo

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    LOGICAL, OPTIONAL :: sypeFaceFlag
    TYPE (t_region), POINTER :: pRegion    

! ==============================================================================
!   Locals     
! ==============================================================================

    LOGICAL :: ignoreSypeFaces
    CHARACTER(CHRLEN) :: errorString
    INTEGER :: c1,c1k,c1t,c2,c2k,c2t,errorFlag,ibv,ic,icl,icg,ifc,ifk, & 
               iPatch,iv,v1,v2,v3,v4
    INTEGER, DIMENSION(:) :: dummyLoc(1),volLoc(2,MIN_VAL:MAX_VAL)
    REAL(RFREAL) :: faceSumMax,fCenX,fCenY,fCenZ,fSumLimit,fVecM,fVecMSum, & 
                    fVecX,fVecY,fVecZ,patchFaceSum,term,volErr,volSum1,volSum2
    REAL(RFREAL), PARAMETER :: THRD = 1.0_RFREAL/3.0_RFREAL, & 
                               VOL_ERR_LIMIT = 1.0E-10_RFREAL
    REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: fnDummy,volDummy
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: faceSum
    REAL(RFREAL) :: xyzAvg(XCOORD:ZCOORD),xyzNodes(XCOORD:ZCOORD,4) 
    TYPE(t_grid), POINTER :: pGrid,pGridOld      
    TYPE(t_patch), POINTER :: pPatch            
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_BuildGeometry',&
  'RFLU_ModGeometry.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building geometry...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

    pGrid => pRegion%grid

    IF ( .NOT. PRESENT(sypeFaceFlag) ) THEN 
      ignoreSypeFaces = .TRUE.
    ELSE 
      ignoreSypeFaces = sypeFaceFlag
    END IF ! PRESENT(ignoreSypeFaces)

! ******************************************************************************
!   Initialize geometry - do here because of mesh motion
! ******************************************************************************

    DO ic = 1,pGrid%nCellsTot ! Explicit loop because of ASCI White problem
      pGrid%vol(ic)         = 0.0_RFREAL  
      pGrid%cofg(XCOORD,ic) = 0.0_RFREAL 
      pGrid%cofg(YCOORD,ic) = 0.0_RFREAL 
      pGrid%cofg(ZCOORD,ic) = 0.0_RFREAL                               
    END DO ! ic

    DO ifc = 1,pGrid%nFacesTot 
      pGrid%fn(XCOORD,ifc) = 0.0_RFREAL
      pGrid%fn(YCOORD,ifc) = 0.0_RFREAL
      pGrid%fn(ZCOORD,ifc) = 0.0_RFREAL
      pGrid%fn(XYZMAG,ifc) = 0.0_RFREAL                                       
      pGrid%fc(XCOORD,ifc) = 0.0_RFREAL
      pGrid%fc(YCOORD,ifc) = 0.0_RFREAL
      pGrid%fc(ZCOORD,ifc) = 0.0_RFREAL                
    END DO ! ifc      

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)  

      pPatch%pc(XCOORD) = 0.0_RFREAL
      pPatch%pc(YCOORD) = 0.0_RFREAL
      pPatch%pc(ZCOORD) = 0.0_RFREAL

      DO ifc = 1,pPatch%nBFacesTot 
        pPatch%fn(XCOORD,ifc) = 0.0_RFREAL
        pPatch%fn(YCOORD,ifc) = 0.0_RFREAL
        pPatch%fn(ZCOORD,ifc) = 0.0_RFREAL
        pPatch%fn(XYZMAG,ifc) = 0.0_RFREAL

        pPatch%fc(XCOORD,ifc) = 0.0_RFREAL 
        pPatch%fc(YCOORD,ifc) = 0.0_RFREAL 
        pPatch%fc(ZCOORD,ifc) = 0.0_RFREAL                               
      END DO ! ifc 

      IF ( pRegion%mixtInput%moveGrid .EQV. .TRUE. ) THEN       
        DO ibv = 1,pPatch%nBVertTot
          pPatch%bvn(XCOORD,ibv) = 0.0_RFREAL
          pPatch%bvn(YCOORD,ibv) = 0.0_RFREAL
          pPatch%bvn(ZCOORD,ibv) = 0.0_RFREAL                        
        END DO ! ibv
      END IF ! pRegion%mixtInput                
    END DO ! iPatch                  

! ******************************************************************************
!   Create and build approximate cell centroids (for later use)
! ******************************************************************************

    CALL RFLU_CreateApproxCentroids(pRegion)
    CALL RFLU_ComputeApproxCentroids(pRegion)

! ******************************************************************************
!   Interior faces
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Non-boundary faces...'
    END IF ! global%verbLevel

    DO ifc = 1,pGrid%nFacesTot
      v1 = pGrid%f2v(1,ifc)
      v2 = pGrid%f2v(2,ifc)
      v3 = pGrid%f2v(3,ifc)

      c1 = pGrid%f2c(1,ifc)
      c2 = pGrid%f2c(2,ifc)

! ==============================================================================
!     Compute face vector and centroid 
! ==============================================================================

      IF ( pGrid%f2v(4,ifc) == VERT_NONE ) THEN ! triangular face
        xyzNodes(XCOORD:ZCOORD,1) = pGrid%xyz(XCOORD:ZCOORD,v1)
        xyzNodes(XCOORD:ZCOORD,2) = pGrid%xyz(XCOORD:ZCOORD,v2)
        xyzNodes(XCOORD:ZCOORD,3) = pGrid%xyz(XCOORD:ZCOORD,v3)

        CALL FaceVectorTria(xyzNodes(XCOORD:ZCOORD,1:3),fVecX,fVecY,fVecZ)
        CALL FaceCentroidTria(xyzNodes(XCOORD:ZCOORD,1:3),fCenX,fCenY,fCenZ)
      ELSE ! quadrilateral face
        v4 = pGrid%f2v(4,ifc)

        xyzNodes(XCOORD:ZCOORD,1) = pGrid%xyz(XCOORD:ZCOORD,v1)
        xyzNodes(XCOORD:ZCOORD,2) = pGrid%xyz(XCOORD:ZCOORD,v2)
        xyzNodes(XCOORD:ZCOORD,3) = pGrid%xyz(XCOORD:ZCOORD,v3)
        xyzNodes(XCOORD:ZCOORD,4) = pGrid%xyz(XCOORD:ZCOORD,v4)

        CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4),fVecX,fVecY,fVecZ) 
        CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4),fCenX,fCenY,fCenZ)
      END IF ! quadFaces

      fVecM = SQRT(fVecX*fVecX + fVecY*fVecY + fVecZ*fVecZ)
      term  = 1.0_RFREAL/fVecM

      pGrid%fn(XCOORD,ifc) = fVecX*term
      pGrid%fn(YCOORD,ifc) = fVecY*term
      pGrid%fn(ZCOORD,ifc) = fVecZ*term
      pGrid%fn(XYZMAG,ifc) = fVecM  

      pGrid%fc(XCOORD,ifc) = fCenX
      pGrid%fc(YCOORD,ifc) = fCenY
      pGrid%fc(ZCOORD,ifc) = fCenZ  

! ==============================================================================
!     Compute contribution to volumes and centroids 
! ==============================================================================

      IF ( c1 /= CELL_TYPE_EXT .AND. c1 /= CELL_TYPE_BND ) THEN
        term = (fCenX - pGrid%cofgApp(XCOORD,c1))*fVecX & 
             + (fCenY - pGrid%cofgApp(YCOORD,c1))*fVecY & 
             + (fCenZ - pGrid%cofgApp(ZCOORD,c1))*fVecZ        

        pGrid%vol(c1) = pGrid%vol(c1) + term        

        term = fCenX*fVecX + fCenY*fVecY + fCenZ*fVecZ         

        pGrid%cofg(XCOORD,c1) = pGrid%cofg(XCOORD,c1) + term*fCenX
        pGrid%cofg(YCOORD,c1) = pGrid%cofg(YCOORD,c1) + term*fCenY
        pGrid%cofg(ZCOORD,c1) = pGrid%cofg(ZCOORD,c1) + term*fCenZ    
      END IF ! c1

      IF ( c2 /= CELL_TYPE_EXT .AND. c2 /= CELL_TYPE_BND ) THEN 
        term = (fCenX - pGrid%cofgApp(XCOORD,c2))*fVecX & 
             + (fCenY - pGrid%cofgApp(YCOORD,c2))*fVecY & 
             + (fCenZ - pGrid%cofgApp(ZCOORD,c2))*fVecZ        

        pGrid%vol(c2) = pGrid%vol(c2) - term        

        term = fCenX*fVecX + fCenY*fVecY + fCenZ*fVecZ

        pGrid%cofg(XCOORD,c2) = pGrid%cofg(XCOORD,c2) - term*fCenX
        pGrid%cofg(YCOORD,c2) = pGrid%cofg(YCOORD,c2) - term*fCenY
        pGrid%cofg(ZCOORD,c2) = pGrid%cofg(ZCOORD,c2) - term*fCenZ
      END IF ! c2        
    END DO ! ifc

! ******************************************************************************
!   Boundary faces, loop over patches
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary faces...'
    END IF ! global%verbLevel

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN
        WRITE(STDOUT,'(A,5X,A,I3)') SOLVER_NAME,'Patch: ',iPatch
      END IF ! global%verbLevel       

      fVecMSum = 0.0_RFREAL

! ==============================================================================
!     Loop over faces
! ==============================================================================

      DO ifc = 1,pPatch%nBFacesTot 
        v1 = pPatch%bv(pPatch%bf2v(1,ifc))
        v2 = pPatch%bv(pPatch%bf2v(2,ifc))
        v3 = pPatch%bv(pPatch%bf2v(3,ifc))

        c1 = pPatch%bf2c(ifc)

! ------------------------------------------------------------------------------
!       Compute face vector and centroid
! ------------------------------------------------------------------------------

        IF ( pPatch%bf2v(4,ifc) == VERT_NONE ) THEN ! triangular face
          xyzNodes(XCOORD:ZCOORD,1) = pGrid%xyz(XCOORD:ZCOORD,v1)
          xyzNodes(XCOORD:ZCOORD,2) = pGrid%xyz(XCOORD:ZCOORD,v2)
          xyzNodes(XCOORD:ZCOORD,3) = pGrid%xyz(XCOORD:ZCOORD,v3)

          CALL FaceVectorTria(xyzNodes(XCOORD:ZCOORD,1:3),fVecX,fVecY,fVecZ)
          CALL FaceCentroidTria(xyzNodes(XCOORD:ZCOORD,1:3),fCenX,fCenY,fCenZ)
        ELSE ! quadrilateral face
          v4 = pPatch%bv(pPatch%bf2v(4,ifc))

          xyzNodes(XCOORD:ZCOORD,1) = pGrid%xyz(XCOORD:ZCOORD,v1)
          xyzNodes(XCOORD:ZCOORD,2) = pGrid%xyz(XCOORD:ZCOORD,v2)
          xyzNodes(XCOORD:ZCOORD,3) = pGrid%xyz(XCOORD:ZCOORD,v3)
          xyzNodes(XCOORD:ZCOORD,4) = pGrid%xyz(XCOORD:ZCOORD,v4)

          CALL FaceVectorQuad(xyzNodes(XCOORD:ZCOORD,1:4),fVecX,fVecY,fVecZ) 
          CALL FaceCentroidQuad(xyzNodes(XCOORD:ZCOORD,1:4),fCenX,fCenY,fCenZ) 
        END IF ! quadFaces

        fVecM = SQRT(fVecX*fVecX + fVecY*fVecY + fVecZ*fVecZ)
        term  = 1.0_RFREAL/fVecM

        pPatch%fn(XCOORD,ifc) = fVecX*term
        pPatch%fn(YCOORD,ifc) = fVecY*term
        pPatch%fn(ZCOORD,ifc) = fVecZ*term
        pPatch%fn(XYZMAG,ifc) = fVecM  

        pPatch%fc(XCOORD,ifc) = fCenX
        pPatch%fc(YCOORD,ifc) = fCenY
        pPatch%fc(ZCOORD,ifc) = fCenZ
           
! ------------------------------------------------------------------------------
!       Accumulate patch centroid 
! ------------------------------------------------------------------------------

        fVecMSum = fVecMSum + fVecM

        pPatch%pc(XCOORD) = pPatch%pc(XCOORD) + fCenX*fVecM
        pPatch%pc(YCOORD) = pPatch%pc(YCOORD) + fCenY*fVecM
        pPatch%pc(ZCOORD) = pPatch%pc(ZCOORD) + fCenZ*fVecM

! ------------------------------------------------------------------------------
!       Add contribution to volume and centroids. NOTE add contribution only if 
!       face is not on symmetry or periodic patch.
! ------------------------------------------------------------------------------

        IF ( ignoreSypeFaces .EQV. .TRUE. ) THEN  
          IF ( pPatch%bcType /= BC_SYMMETRY .AND. & 
               pPatch%bcType /= BC_PERIODIC ) THEN
            term = (fCenX - pGrid%cofgApp(XCOORD,c1))*fVecX & 
                 + (fCenY - pGrid%cofgApp(YCOORD,c1))*fVecY & 
                 + (fCenZ - pGrid%cofgApp(ZCOORD,c1))*fVecZ        

            pGrid%vol(c1) = pGrid%vol(c1) + term 

            term = fCenX*fVecX + fCenY*fVecY + fCenZ*fVecZ  

            pGrid%cofg(XCOORD,c1) = pGrid%cofg(XCOORD,c1) + term*fCenX
            pGrid%cofg(YCOORD,c1) = pGrid%cofg(YCOORD,c1) + term*fCenY
            pGrid%cofg(ZCOORD,c1) = pGrid%cofg(ZCOORD,c1) + term*fCenZ
          END IF ! pPatch%bcType
        ELSE 
         term = (fCenX - pGrid%cofgApp(XCOORD,c1))*fVecX & 
              + (fCenY - pGrid%cofgApp(YCOORD,c1))*fVecY & 
              + (fCenZ - pGrid%cofgApp(ZCOORD,c1))*fVecZ        

          pGrid%vol(c1) = pGrid%vol(c1) + term 

          term = fCenX*fVecX + fCenY*fVecY + fCenZ*fVecZ  

          pGrid%cofg(XCOORD,c1) = pGrid%cofg(XCOORD,c1) + term*fCenX
          pGrid%cofg(YCOORD,c1) = pGrid%cofg(YCOORD,c1) + term*fCenY
          pGrid%cofg(ZCOORD,c1) = pGrid%cofg(ZCOORD,c1) + term*fCenZ          
        END IF ! ignoreSypeFaces
      END DO ! ifc

! =============================================================================
!     Finalize patch centroid
! =============================================================================

      term = 1.0_RFREAL/fVecMSum

      pPatch%pc(XCOORD) = term*pPatch%pc(XCOORD)
      pPatch%pc(YCOORD) = term*pPatch%pc(YCOORD)
      pPatch%pc(ZCOORD) = term*pPatch%pc(ZCOORD)
    END DO ! iPatch

! ******************************************************************************
!   Finalize volumes and cell centroids
! ******************************************************************************

    DO ic = 1,pGrid%nCellsTot
      pGrid%vol(ic)  = THRD*pGrid%vol(ic)   

      term = 1.0_RFREAL/(4.0_RFREAL*pGrid%vol(ic))  

      pGrid%cofg(XCOORD,ic) = term*pGrid%cofg(XCOORD,ic)
      pGrid%cofg(YCOORD,ic) = term*pGrid%cofg(YCOORD,ic)
      pGrid%cofg(ZCOORD,ic) = term*pGrid%cofg(ZCOORD,ic) 
    END DO ! ic

! ******************************************************************************
!   Destroy approximate cell centroids
! ******************************************************************************

    CALL RFLU_DestroyApproxCentroids(pRegion)

! ******************************************************************************
!   Print out some information for diagnostic purposes
! ******************************************************************************

! ==============================================================================
!   Volume information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      dummyLoc  = MINLOC(pGrid%vol(1:pGrid%nCells))
      volLoc(1,MIN_VAL) = dummyLoc(1)
      volLoc(1,MAX_VAL) = dummyLoc(1)        

      dummyLoc  = MAXLOC(pGrid%vol(1:pGrid%nCells))                        
      volLoc(2,MIN_VAL) = dummyLoc(1)
      volLoc(2,MAX_VAL) = dummyLoc(1)                    

      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Statistics (actual cells only):'
      WRITE(STDOUT,'(A,5X,A,7X,E23.16,1X,I8)') SOLVER_NAME, &
        'Minimum volume:',MINVAL(pGrid%vol(1:pGrid%nCells)),volLoc(1,MIN_VAL)
      WRITE(STDOUT,'(A,5X,A,7X,E23.16,1X,I8)') SOLVER_NAME, &
        'Maximum volume:',MAXVAL(pGrid%vol(1:pGrid%nCells)),volLoc(2,MAX_VAL)

      CALL RFLU_PrintLocInfo(pRegion,volLoc,2,LOCINFO_MODE_VERBOSE, & 
                             OUTPUT_MODE_MASTER_ONLY)            
    END IF ! global%verbLevel

    IF ( MINVAL(pGrid%vol(1:pGrid%nCells)) <= 0.0_RFREAL ) THEN 
      CALL ErrorStop(global,ERR_VOLUME_NEGATIVE,__LINE__)
    END IF ! MINVAL

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      IF ( ASSOCIATED(pRegion%gridOld%vol) .EQV. .TRUE. ) THEN 
        IF ( MINVAL(pRegion%gridOld%vol(1:pGrid%nCells)) > 0.0_RFREAL ) THEN
          pGridOld => pRegion%gridOld

          WRITE(STDOUT,'(A,5X,A,1X,E23.16,2(1X,I8))') SOLVER_NAME, & 
            'Minimum volume ratio:', & 
            MINVAL(pGrid%vol(1:pGrid%nCells))/ &
            MINVAL(pGridOld%vol(1:pGrid%nCells)), &             
            MINLOC(pGrid%vol(1:pGrid%nCells)), &
            MINLOC(pGridOld%vol(1:pGrid%nCells))
          WRITE(STDOUT,'(A,5X,A,1X,E23.16,2(1X,I8))') SOLVER_NAME, & 
            'Maximum volume ratio:',& 
            MAXVAL(pGrid%vol(1:pGrid%nCells))/ &
            MAXVAL(pGridOld%vol(1:pGrid%nCells)), &             
            MAXLOC(pGrid%vol(1:pGrid%nCells)), &
            MAXLOC(pGridOld%vol(1:pGrid%nCells))
        END IF ! MINVAL 
      END IF ! ASSOCIATED
    END IF ! global%myProcid 

! ==============================================================================
!   Patch areas (only actual ones - because of parallel runs)
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,5X,A,1X,A)') SOLVER_NAME,'Boundary patch areas', & 
                                    '(actual faces only):'

      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        ALLOCATE(fnDummy(pPatch%nBFaces),STAT=errorFlag)
        global%error = errorFlag   
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'fnDummy')
        END IF ! global%error    

        DO ifc = 1,pPatch%nBFaces
          fnDummy(ifc) = pPatch%fn(XYZMAG,ifc)
        END DO ! ifc

        IF ( pPatch%nBFaces > 0 ) THEN ! Might have 0 actual faces...
          CALL QuickSortRFREAL(fnDummy,pPatch%nBFaces) ! Put in ascending order
        END IF ! pPatch%nBFaces

        patchFaceSum = 0.0_RFREAL

        DO ifc = 1,pPatch%nBFaces
          patchFaceSum = patchFaceSum + fnDummy(ifc)
        END DO ! ifc

        DEALLOCATE(fnDummy,STAT=errorFlag)
        global%error = errorFlag   
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'fnDummy')
        END IF ! global%error              

        WRITE(STDOUT,'(A,7X,A,1X,I3,3X,E23.16)') SOLVER_NAME,'Patch:', & 
                                                 iPatch,patchFaceSum
      END DO ! iPatch
    END IF ! global%myProcid

! ******************************************************************************
!   Check volume computation: total volume
! ******************************************************************************

    IF ( global%checkLevel > CHECK_NONE ) THEN
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN  
        WRITE(STDOUT,'(A,3X,A,1X,A)') SOLVER_NAME,'Check total volume', & 
                                      '(actual cells only):'
      END IF ! global%verbLevel

! ==============================================================================
!     Part 1: Sum individual volumes. NOTE want only sum interior volumes, so 
!     need to be careful with summation because sorting means that interior
!     and dummy cells are jumbled up. Therefore copy only interior cells into
!     volDummy array and then sort   
! ==============================================================================

      volSum1 = 0.0_RFREAL  

      ALLOCATE(volDummy(pGrid%nCells),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'volDummy')
      END IF ! global%error    

      IF ( pGrid%nTets > 0 ) THEN 
        DO icl = 1,pGrid%nTets
          icg = pGrid%tet2CellGlob(icl)
          volDummy(icg) = pGrid%vol(icg)
        END DO ! icl
      END IF ! pGrid%nTets

      IF ( pGrid%nHexs > 0 ) THEN
        DO icl = 1,pGrid%nHexs
          icg = pGrid%hex2CellGlob(icl)
          volDummy(icg) = pGrid%vol(icg)
        END DO ! icl
      END IF ! pGrid%nHexs        

      IF ( pGrid%nPris > 0 ) THEN
        DO icl = 1,pGrid%nPris
          icg = pGrid%pri2CellGlob(icl)
          volDummy(icg) = pGrid%vol(icg)
        END DO ! icl
      END IF ! pGrid%nPris 

      IF ( pGrid%nPyrs > 0 ) THEN
        DO icl = 1,pGrid%nPyrs
          icg = pGrid%pyr2CellGlob(icl)
          volDummy(icg) = pGrid%vol(icg)
        END DO ! icl
      END IF ! pGrid%nPyrs 

      CALL QuickSortRFREAL(volDummy,pGrid%nCells) ! Put in ascending order

      DO ic = 1,pGrid%nCells  
        volSum1 = volSum1 + volDummy(ic)
      END DO ! ic

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN  
        WRITE(STDOUT,'(A,5X,A,1X,E23.16)') SOLVER_NAME,& 
          'Total volume from sum of control volumes:',volSum1            
      END IF ! global%verbLevel

      DEALLOCATE(volDummy,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'volDummy')
      END IF ! global%error 

! ==============================================================================
!     Part 2: Compute total volume from boundary faces
! ==============================================================================

! ------------------------------------------------------------------------------
!     Compute approximate centroid for region 
! ------------------------------------------------------------------------------

      xyzAvg(XCOORD) = 0.0_RFREAL
      xyzAvg(YCOORD) = 0.0_RFREAL
      xyzAvg(ZCOORD) = 0.0_RFREAL

      DO iv = 1,pGrid%nVertTot
        xyzAvg(XCOORD) = xyzAvg(XCOORD) + pGrid%xyz(XCOORD,iv)
        xyzAvg(YCOORD) = xyzAvg(YCOORD) + pGrid%xyz(YCOORD,iv)
        xyzAvg(ZCOORD) = xyzAvg(ZCOORD) + pGrid%xyz(ZCOORD,iv)
      END DO ! iv

      term = 1.0_RFREAL/REAL(pGrid%nVertTot,KIND=RFREAL)

      xyzAvg(XCOORD) = term*xyzAvg(XCOORD)
      xyzAvg(YCOORD) = term*xyzAvg(YCOORD)
      xyzAvg(ZCOORD) = term*xyzAvg(ZCOORD)        

! ------------------------------------------------------------------------------
!     Compute volume. NOTE need to leave out faces on symmetry and periodic 
!     patches because taken into account through interpartition faces below.
! ------------------------------------------------------------------------------

      volSum2 = 0.0_RFREAL

      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        IF ( ignoreSypeFaces .EQV. .TRUE. ) THEN  
          IF ( pPatch%bcType /= BC_SYMMETRY .AND. & 
               pPatch%bcType /= BC_PERIODIC ) THEN
            DO ifc = 1,pPatch%nBFaces
              fVecX = pPatch%fn(XCOORD,ifc)
              fVecY = pPatch%fn(YCOORD,ifc)
              fVecZ = pPatch%fn(ZCOORD,ifc)
              fVecM = pPatch%fn(XYZMAG,ifc)

              fCenX = pPatch%fc(XCOORD,ifc)
              fCenY = pPatch%fc(YCOORD,ifc)
              fCenZ = pPatch%fc(ZCOORD,ifc)

              volSum2 = volSum2 + ((fCenX - xyzAvg(XCOORD))*fVecX & 
                                +  (fCenY - xyzAvg(YCOORD))*fVecY & 
                                +  (fCenZ - xyzAvg(ZCOORD))*fVecZ)*fVecM
            END DO ! ifc
          END IF ! pPatch%bcType
        ELSE 
          DO ifc = 1,pPatch%nBFaces
            fVecX = pPatch%fn(XCOORD,ifc)
            fVecY = pPatch%fn(YCOORD,ifc)
            fVecZ = pPatch%fn(ZCOORD,ifc)
            fVecM = pPatch%fn(XYZMAG,ifc)

            fCenX = pPatch%fc(XCOORD,ifc)
            fCenY = pPatch%fc(YCOORD,ifc)
            fCenZ = pPatch%fc(ZCOORD,ifc)

            volSum2 = volSum2 + ((fCenX - xyzAvg(XCOORD))*fVecX & 
                              +  (fCenY - xyzAvg(YCOORD))*fVecY & 
                              +  (fCenZ - xyzAvg(ZCOORD))*fVecZ)*fVecM
          END DO ! ifc                      
        END IF ! ignoreSypefaces
      END DO ! iPatch

! ==============================================================================
!     Part 3: Include contribution of interpartition faces. NOTE need to take 
!     into account direction of face. NOTE the second comparison in the if 
!     statement (ifk == FACE_KIND_AB) is retained so that should no dummy  
!     cells be included in a debugging run, one still gets the correct
!     computation of the total volume. In normal runs, where dummy cells 
!     exist, this statement should never be active
! ==============================================================================

      DO ifc = 1,pGrid%nFacesTot      
        c1 = pGrid%f2c(1,ifc)
        c2 = pGrid%f2c(2,ifc)

        c1k = RFLU_GetGlobalCellKind(global,pGrid,c1)        
        c2k = RFLU_GetGlobalCellKind(global,pGrid,c2)  
        ifk = RFLU_GetFaceKind(global,c1k,c2k)        

        IF ( (ifk == FACE_KIND_AV) .OR. (ifk == FACE_KIND_AB) ) THEN 
          fVecX = pGrid%fn(XCOORD,ifc)
          fVecY = pGrid%fn(YCOORD,ifc)
          fVecZ = pGrid%fn(ZCOORD,ifc)
          fVecM = pGrid%fn(XYZMAG,ifc)

          fCenX = pGrid%fc(XCOORD,ifc)
          fCenY = pGrid%fc(YCOORD,ifc)
          fCenZ = pGrid%fc(ZCOORD,ifc)

          term = ((fCenX - xyzAvg(XCOORD))*fVecX & 
               +  (fCenY - xyzAvg(YCOORD))*fVecY & 
               +  (fCenZ - xyzAvg(ZCOORD))*fVecZ)*fVecM

          IF ( c1k == CELL_KIND_ACTUAL ) THEN 
            volSum2 = volSum2 + term
          ELSE        
            volSum2 = volSum2 - term
          END IF ! c1k                                 
        END IF ! ifk
      END DO ! ifc

      volSum2 = THRD*volSum2
      volErr  = (volSum2-volSum1)/(0.5_RFREAL*(volSum1+volSum2)*100.0_RFREAL)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN
        WRITE(STDOUT,'(A,5X,A,4X,E23.16)') SOLVER_NAME,& 
          'Total volume from boundary polyhedron:',volSum2
        WRITE(STDOUT,'(A,5X,A,4X,E23.16)') SOLVER_NAME,& 
          'Error (in % of average total volume): ',volErr
      END IF ! global%verbLevel

      IF ( ABS(volErr) > VOL_ERR_LIMIT ) THEN
        WRITE(errorString,'(A,1X,E13.6)') 'Error:',volErr                  
        CALL ErrorStop(global,ERR_VOLUME_DIFF,__LINE__,TRIM(errorString))
      END IF ! volErr
    END IF ! global%checkLevel

! ******************************************************************************
!   Check volume computation: closed-ness of control volumes. NOTE need to 
!   check this for ALL volumes so that can make sure that geometry is computed
!   correctly for virtual cells.
! ******************************************************************************

    IF ( global%checkLevel > CHECK_NONE ) THEN 
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Check closedness of '// & 
                                 'control volumes (all cells):'
      END IF ! global%verbLevel

      ALLOCATE(faceSum(XCOORD:ZCOORD,pGrid%nCellsTot),STAT=errorFlag)
      global%error = errorFlag    
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'faceSum')
      END IF ! global%error

      DO ic = 1,pGrid%nCellsTot
        faceSum(XCOORD,ic) = 0.0_RFREAL
        faceSum(YCOORD,ic) = 0.0_RFREAL
        faceSum(ZCOORD,ic) = 0.0_RFREAL                
      END DO ! ic

! ==============================================================================
!     Part 1: Contribution from interior faces  
! ==============================================================================

      DO ifc = 1,pGrid%nFacesTot
        c1 = pGrid%f2c(1,ifc)
        c2 = pGrid%f2c(2,ifc)

        fVecX = pGrid%fn(XCOORD,ifc)
        fVecY = pGrid%fn(YCOORD,ifc)
        fVecZ = pGrid%fn(ZCOORD,ifc)       
        fVecM = pGrid%fn(XYZMAG,ifc)

        IF ( c1 /= CELL_TYPE_EXT .AND. c1 /= CELL_TYPE_BND ) THEN     
          faceSum(XCOORD,c1) = faceSum(XCOORD,c1) + fVecX*fVecM
          faceSum(YCOORD,c1) = faceSum(YCOORD,c1) + fVecY*fVecM
          faceSum(ZCOORD,c1) = faceSum(ZCOORD,c1) + fVecZ*fVecM 
        END IF ! c1  

        IF ( c2 /= CELL_TYPE_EXT .AND. c2 /= CELL_TYPE_BND ) THEN       
          faceSum(XCOORD,c2) = faceSum(XCOORD,c2) - fVecX*fVecM
          faceSum(YCOORD,c2) = faceSum(YCOORD,c2) - fVecY*fVecM
          faceSum(ZCOORD,c2) = faceSum(ZCOORD,c2) - fVecZ*fVecM          
        END IF ! c2
      END DO ! ifc

! ==============================================================================
!     Part 2: Contribution from boundary faces
! ==============================================================================

      DO iPatch = 1,pGrid%nPatches
        pPatch => pRegion%patches(iPatch)

        IF ( ignoreSypeFaces .EQV. .TRUE. ) THEN  
          IF ( pPatch%bcType /= BC_SYMMETRY .AND. & 
               pPatch%bcType /= BC_PERIODIC ) THEN
            DO ifc = 1,pPatch%nBFacesTot
              c1 = pPatch%bf2c(ifc)

              fVecX = pPatch%fn(XCOORD,ifc)
              fVecY = pPatch%fn(YCOORD,ifc)
              fVecZ = pPatch%fn(ZCOORD,ifc)
              fVecM = pPatch%fn(XYZMAG,ifc)

              faceSum(XCOORD,c1) = faceSum(XCOORD,c1) + fVecX*fVecM
              faceSum(YCOORD,c1) = faceSum(YCOORD,c1) + fVecY*fVecM
              faceSum(ZCOORD,c1) = faceSum(ZCOORD,c1) + fVecZ*fVecM
            END DO ! ifc
          END IF ! pPatch%bcType
        ELSE 
          DO ifc = 1,pPatch%nBFacesTot
            c1 = pPatch%bf2c(ifc)

            fVecX = pPatch%fn(XCOORD,ifc)
            fVecY = pPatch%fn(YCOORD,ifc)
            fVecZ = pPatch%fn(ZCOORD,ifc)
            fVecM = pPatch%fn(XYZMAG,ifc)

            faceSum(XCOORD,c1) = faceSum(XCOORD,c1) + fVecX*fVecM
            faceSum(YCOORD,c1) = faceSum(YCOORD,c1) + fVecY*fVecM
            faceSum(ZCOORD,c1) = faceSum(ZCOORD,c1) + fVecZ*fVecM
          END DO ! ifc        
        END IF ! ignoreSypeFaces
      END DO ! iPatch

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN 
        WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Minimum/maximum value of '// & 
              'sum of face vectors:'
        WRITE(STDOUT,'(A,7X,A,2(1X,E23.16),2(1X,I8))') SOLVER_NAME, & 
              'x-direction:', & 
              MINVAL(faceSum(XCOORD,1:pGrid%nCellsTot)), & 
              MAXVAL(faceSum(XCOORD,1:pGrid%nCellsTot)), & 
              MINLOC(faceSum(XCOORD,1:pGrid%nCellsTot)), & 
              MAXLOC(faceSum(XCOORD,1:pGrid%nCellsTot))
        WRITE(STDOUT,'(A,7X,A,2(1X,E23.16),2(1X,I8))') SOLVER_NAME, & 
              'y-direction:', & 
              MINVAL(faceSum(YCOORD,1:pGrid%nCellsTot)), & 
              MAXVAL(faceSum(YCOORD,1:pGrid%nCellsTot)), & 
              MINLOC(faceSum(YCOORD,1:pGrid%nCellsTot)), & 
              MAXLOC(faceSum(YCOORD,1:pGrid%nCellsTot))
        WRITE(STDOUT,'(A,7X,A,2(1X,E23.16),2(1X,I8))') SOLVER_NAME, & 
              'z-direction:', & 
              MINVAL(faceSum(ZCOORD,1:pGrid%nCellsTot)), & 
              MAXVAL(faceSum(ZCOORD,1:pGrid%nCellsTot)), & 
              MINLOC(faceSum(ZCOORD,1:pGrid%nCellsTot)), & 
              MAXLOC(faceSum(ZCOORD,1:pGrid%nCellsTot))
      END IF ! global%verbLevel  

      faceSumMax = MAX(ABS(MINVAL(faceSum(XCOORD,1:pGrid%nCellsTot))), & 
                       ABS(MAXVAL(faceSum(XCOORD,1:pGrid%nCellsTot))), & 
                       ABS(MINVAL(faceSum(YCOORD,1:pGrid%nCellsTot))), & 
                       ABS(MAXVAL(faceSum(YCOORD,1:pGrid%nCellsTot))), & 
                       ABS(MINVAL(faceSum(ZCOORD,1:pGrid%nCellsTot))), & 
                       ABS(MAXVAL(faceSum(ZCOORD,1:pGrid%nCellsTot))))  

      IF ( faceSumMax > MINVAL(pGrid%fn(XYZMAG,1:pGrid%nFacesTot)) ) THEN
        WRITE(errorString,'(A,1X,E13.6)') 'Error:',faceSumMax 
        CALL ErrorStop(global,ERR_FACESUM,__LINE__,TRIM(errorString))
      END IF ! faceSumMax  

      DEALLOCATE(faceSum,STAT=errorFlag)
      global%error = errorFlag    
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'faceSum')
      END IF ! global%error
    END IF ! global%checkLevel

#ifdef CHECK_DATASTRUCT
! ******************************************************************************
!   Data structure output for checking
! ******************************************************************************

    WRITE(STDOUT,'(A)') SOLVER_NAME
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Cell centroid locations'
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Number of actual cells: ', & 
          pGrid%nCells
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Number of virtual cells:', & 
          pGrid%nCellsTot-pGrid%nCells                               

    DO ic = 1,pGrid%nCellsTot
      WRITE(STDOUT,'(A,1X,I6,3(1X,E18.9))') SOLVER_NAME,ic, & 
                                            pGrid%cofg(XCOORD,ic), & 
                                            pGrid%cofg(YCOORD,ic), & 
                                            pGrid%cofg(ZCOORD,ic)
    END DO ! ic

    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
    WRITE(STDOUT,'(A)') SOLVER_NAME
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Cell volumes'
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Number of actual cells: ', & 
          pGrid%nCells
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Number of virtual cells:', & 
          pGrid%nCellsTot-pGrid%nCells       

    DO ic = 1,pGrid%nCellsTot
      WRITE(STDOUT,'(A,1X,I6,1X,E18.9)') SOLVER_NAME,ic,pGrid%vol(ic)
    END DO ! ic

    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
    WRITE(STDOUT,'(A)') SOLVER_NAME   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###' 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Face centroid locations'

    DO ifc = 1,pGrid%nFacesTot
      WRITE(STDOUT,'(A,1X,I6,3(1X,E18.9))') SOLVER_NAME,ifc, & 
                                            pGrid%fc(XCOORD,ifc), & 
                                            pGrid%fc(YCOORD,ifc), & 
                                            pGrid%fc(ZCOORD,ifc)
    END DO ! ifc

    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
    WRITE(STDOUT,'(A)') SOLVER_NAME   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###' 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Face normals'

    DO ifc = 1,pGrid%nFacesTot
      WRITE(STDOUT,'(A,1X,I6,3(1X,E18.9))') SOLVER_NAME,ifc, & 
                                            pGrid%fn(XCOORD,ifc), & 
                                            pGrid%fn(YCOORD,ifc), & 
                                            pGrid%fn(ZCOORD,ifc)
    END DO ! ifc

    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
    WRITE(STDOUT,'(A)') SOLVER_NAME                            
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###' 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Boundary face centroid locations'

    DO iPatch = 1,pGrid%nPatches 
      pPatch => pRegion%patches(iPatch)
      WRITE(STDOUT,'(A,3X,A,1X,I7)') SOLVER_NAME,'Patch:',iPatch  
      WRITE(STDOUT,'(A,3X,A,1X,I7)') SOLVER_NAME, & 
                                     'Actual number of faces:', & 
                                     pPatch%nBFaces
      WRITE(STDOUT,'(A,3X,A,1X,I7)') SOLVER_NAME, & 
                                     'Total number of faces: ', & 
                                     pPatch%nBFacesTot

      DO ifc = 1,pPatch%nBFacesTot
        WRITE(STDOUT,'(A,1X,I6,3(1X,E18.9))') SOLVER_NAME,ifc, & 
                                              pPatch%fc(XCOORD,ifc), & 
                                              pPatch%fc(YCOORD,ifc), & 
                                              pPatch%fc(ZCOORD,ifc)
      END DO ! ifc
    END DO ! iPatch

    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
    WRITE(STDOUT,'(A)') SOLVER_NAME                            
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###' 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Boundary face normals'

    DO iPatch = 1,pGrid%nPatches 
      pPatch => pRegion%patches(iPatch)
      WRITE(STDOUT,'(A,3X,A,1X,I7)') SOLVER_NAME,'Patch:',iPatch  
      WRITE(STDOUT,'(A,3X,A,1X,I7)') SOLVER_NAME, & 
                                     'Actual number of faces:', & 
                                     pPatch%nBFaces
      WRITE(STDOUT,'(A,3X,A,1X,I7)') SOLVER_NAME, & 
                                     'Total number of faces: ', & 
                                     pPatch%nBFacesTot

      DO ifc = 1,pPatch%nBFacesTot
        WRITE(STDOUT,'(A,1X,I6,3(1X,E18.9))') SOLVER_NAME,ifc, & 
                                              pPatch%fn(XCOORD,ifc), & 
                                              pPatch%fn(YCOORD,ifc), & 
                                              pPatch%fn(ZCOORD,ifc)
      END DO ! ifc
    END DO ! iPatch

    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
    WRITE(STDOUT,'(A)') SOLVER_NAME 
#endif

! ******************************************************************************
!   Compute boundary vertex normals for moving grids
! ******************************************************************************

    IF ( pRegion%mixtInput%moveGrid .EQV. .TRUE. ) THEN 
      CALL RFLU_BuildBVertexNormals(pRegion)
    END IF ! pRegion%mixtInput

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Building geometry done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_BuildGeometry







! ******************************************************************************
!
! Purpose: Compute approximate cell centroids from simple vertex-coordinate 
!  averages
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

  SUBROUTINE RFLU_ComputeApproxCentroids(pRegion)

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals     
! ==============================================================================

    INTEGER :: icl,icg,v1,v2,v3,v4,v5,v6,v7,v8
    REAL(RFREAL) :: term,x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4,y5,y6,y7,y8, & 
                    z1,z2,z3,z4,z5,z6,z7,z8
    TYPE(t_grid), POINTER :: pGrid              
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeApproxCentroids',&
  'RFLU_ModGeometry.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Computing approximate ', & 
                                 'centroids...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Loop over cell types and compute approximate centroids
! ******************************************************************************

! ==============================================================================
!   Tetrahedra
! ==============================================================================

    term = 1.0_RFREAL/4.0_RFREAL 

    DO icl = 1,pGrid%nTetsTot
      v1 = pGrid%tet2v(1,icl)
      v2 = pGrid%tet2v(2,icl)
      v3 = pGrid%tet2v(3,icl)
      v4 = pGrid%tet2v(4,icl)                        

      x1 = pGrid%xyz(XCOORD,v1)
      x2 = pGrid%xyz(XCOORD,v2)        
      x3 = pGrid%xyz(XCOORD,v3)
      x4 = pGrid%xyz(XCOORD,v4)

      y1 = pGrid%xyz(YCOORD,v1)
      y2 = pGrid%xyz(YCOORD,v2)        
      y3 = pGrid%xyz(YCOORD,v3)
      y4 = pGrid%xyz(YCOORD,v4)

      z1 = pGrid%xyz(ZCOORD,v1)
      z2 = pGrid%xyz(ZCOORD,v2)        
      z3 = pGrid%xyz(ZCOORD,v3)
      z4 = pGrid%xyz(ZCOORD,v4)

      icg = pGrid%tet2CellGlob(icl)

      pGrid%cofgApp(XCOORD,icg) = term*(x1 + x2 + x3 + x4)
      pGrid%cofgApp(YCOORD,icg) = term*(y1 + y2 + y3 + y4)
      pGrid%cofgApp(ZCOORD,icg) = term*(z1 + z2 + z3 + z4)               
    END DO ! icl

! ==============================================================================
!   Hexahedra
! ==============================================================================

    term = 1.0_RFREAL/8.0_RFREAL 

    DO icl = 1,pGrid%nHexsTot
      v1 = pGrid%hex2v(1,icl)
      v2 = pGrid%hex2v(2,icl)
      v3 = pGrid%hex2v(3,icl)
      v4 = pGrid%hex2v(4,icl)                        
      v5 = pGrid%hex2v(5,icl)
      v6 = pGrid%hex2v(6,icl)
      v7 = pGrid%hex2v(7,icl)
      v8 = pGrid%hex2v(8,icl)                        

      x1 = pGrid%xyz(XCOORD,v1)
      x2 = pGrid%xyz(XCOORD,v2)        
      x3 = pGrid%xyz(XCOORD,v3)
      x4 = pGrid%xyz(XCOORD,v4)
      x5 = pGrid%xyz(XCOORD,v5)
      x6 = pGrid%xyz(XCOORD,v6)        
      x7 = pGrid%xyz(XCOORD,v7)
      x8 = pGrid%xyz(XCOORD,v8)

      y1 = pGrid%xyz(YCOORD,v1)
      y2 = pGrid%xyz(YCOORD,v2)        
      y3 = pGrid%xyz(YCOORD,v3)
      y4 = pGrid%xyz(YCOORD,v4)
      y5 = pGrid%xyz(YCOORD,v5)
      y6 = pGrid%xyz(YCOORD,v6)        
      y7 = pGrid%xyz(YCOORD,v7)
      y8 = pGrid%xyz(YCOORD,v8)

      z1 = pGrid%xyz(ZCOORD,v1)
      z2 = pGrid%xyz(ZCOORD,v2)        
      z3 = pGrid%xyz(ZCOORD,v3)
      z4 = pGrid%xyz(ZCOORD,v4)
      z5 = pGrid%xyz(ZCOORD,v5)
      z6 = pGrid%xyz(ZCOORD,v6)        
      z7 = pGrid%xyz(ZCOORD,v7)
      z8 = pGrid%xyz(ZCOORD,v8)

      icg = pGrid%hex2CellGlob(icl)

      pGrid%cofgApp(XCOORD,icg) = term*(x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8)
      pGrid%cofgApp(YCOORD,icg) = term*(y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8)
      pGrid%cofgApp(ZCOORD,icg) = term*(z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8)
    END DO ! icl

! ==============================================================================
!   Prisms
! ==============================================================================

    term = 1.0_RFREAL/6.0_RFREAL 

    DO icl = 1,pGrid%nPrisTot
      v1 = pGrid%pri2v(1,icl)
      v2 = pGrid%pri2v(2,icl)
      v3 = pGrid%pri2v(3,icl)
      v4 = pGrid%pri2v(4,icl)                        
      v5 = pGrid%pri2v(5,icl)
      v6 = pGrid%pri2v(6,icl)

      x1 = pGrid%xyz(XCOORD,v1)
      x2 = pGrid%xyz(XCOORD,v2)        
      x3 = pGrid%xyz(XCOORD,v3)
      x4 = pGrid%xyz(XCOORD,v4)
      x5 = pGrid%xyz(XCOORD,v5)
      x6 = pGrid%xyz(XCOORD,v6)        

      y1 = pGrid%xyz(YCOORD,v1)
      y2 = pGrid%xyz(YCOORD,v2)        
      y3 = pGrid%xyz(YCOORD,v3)
      y4 = pGrid%xyz(YCOORD,v4)
      y5 = pGrid%xyz(YCOORD,v5)
      y6 = pGrid%xyz(YCOORD,v6)        

      z1 = pGrid%xyz(ZCOORD,v1)
      z2 = pGrid%xyz(ZCOORD,v2)        
      z3 = pGrid%xyz(ZCOORD,v3)
      z4 = pGrid%xyz(ZCOORD,v4)
      z5 = pGrid%xyz(ZCOORD,v5)
      z6 = pGrid%xyz(ZCOORD,v6)        

      icg = pGrid%pri2CellGlob(icl)

      pGrid%cofgApp(XCOORD,icg) = term*(x1 + x2 + x3 + x4 + x5 + x6)
      pGrid%cofgApp(YCOORD,icg) = term*(y1 + y2 + y3 + y4 + y5 + y6)
      pGrid%cofgApp(ZCOORD,icg) = term*(z1 + z2 + z3 + z4 + z5 + z6)
    END DO ! icl

! ==============================================================================
!   Pyramids
! ==============================================================================

    term = 1.0_RFREAL/5.0_RFREAL 

    DO icl = 1,pGrid%nPyrsTot
      v1 = pGrid%pyr2v(1,icl)
      v2 = pGrid%pyr2v(2,icl)
      v3 = pGrid%pyr2v(3,icl)
      v4 = pGrid%pyr2v(4,icl)                        
      v5 = pGrid%pyr2v(5,icl)

      x1 = pGrid%xyz(XCOORD,v1)
      x2 = pGrid%xyz(XCOORD,v2)        
      x3 = pGrid%xyz(XCOORD,v3)
      x4 = pGrid%xyz(XCOORD,v4)
      x5 = pGrid%xyz(XCOORD,v5)

      y1 = pGrid%xyz(YCOORD,v1)
      y2 = pGrid%xyz(YCOORD,v2)        
      y3 = pGrid%xyz(YCOORD,v3)
      y4 = pGrid%xyz(YCOORD,v4)
      y5 = pGrid%xyz(YCOORD,v5)

      z1 = pGrid%xyz(ZCOORD,v1)
      z2 = pGrid%xyz(ZCOORD,v2)        
      z3 = pGrid%xyz(ZCOORD,v3)
      z4 = pGrid%xyz(ZCOORD,v4)
      z5 = pGrid%xyz(ZCOORD,v5)

      icg = pGrid%pyr2CellGlob(icl)

      pGrid%cofgApp(XCOORD,icg) = term*(x1 + x2 + x3 + x4 + x5)
      pGrid%cofgApp(YCOORD,icg) = term*(y1 + y2 + y3 + y4 + y5)
      pGrid%cofgApp(ZCOORD,icg) = term*(z1 + z2 + z3 + z4 + z5)               
    END DO ! icl

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Computing approximate '// & 
                                 'centroids done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeApproxCentroids






! ******************************************************************************
!
! Purpose: Compute distance from face to cell centroids. 
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
    
  SUBROUTINE RFLU_ComputeFaceDist(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: c1,c2,errorFlag,ifc,iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ComputeFaceDist',&
  'RFLU_ModGeometry.F90')

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing face distance...' 
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

    pGrid => pRegion%grid  

! ******************************************************************************
!   Compute distances
! ******************************************************************************

! ==============================================================================
!   Interior faces      
! ==============================================================================
    
    DO ifc = 1,pGrid%nFaces
      c1 = pGrid%f2c(1,ifc)
      c2 = pGrid%f2c(2,ifc)

      pGrid%cofgDist(1,ifc) = DOT_PRODUCT(pGrid%fc(:,ifc)-pGrid%cofg(:,c1), & 
                                          pGrid%fn(1:3,ifc))
      pGrid%cofgDist(2,ifc) = DOT_PRODUCT(pGrid%cofg(:,c2)-pGrid%fc(:,ifc), &
                                          pGrid%fn(1:3,ifc))
    END DO ! ifc

! ==============================================================================
!   Boundary faces         
! ==============================================================================

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      DO ifc = 1,pPatch%nBFaces
        c1 = pPatch%bf2c(ifc)

        pPatch%cofgDist(ifc) = DOT_PRODUCT(pPatch%fc(:,ifc)-pGrid%cofg(:,c1), &
                                           pPatch%fn(1:3,ifc))
      END DO ! ifc
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel >= VERBOSE_HIGH ) THEN  
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Computing face distance done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ComputeFaceDist






! ******************************************************************************
!
! Purpose: Create approximate cell centroids.
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

  SUBROUTINE RFLU_CreateApproxCentroids(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
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

    CALL RegisterFunction(global,'RFLU_CreateApproxCentroids',&
  'RFLU_ModGeometry.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating approximate centroids...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyApproxCentroids(pRegion)

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    ALLOCATE(pGrid%cofgApp(XCOORD:ZCOORD,pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%cofgApp')
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Creating approximate '// & 
                                 'centroids done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_CreateApproxCentroids






! ******************************************************************************
!
! Purpose: Create distance from cell centroid to face centroid.
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

  SUBROUTINE RFLU_CreateFaceDist(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: errorFlag,iPatch
    TYPE(t_grid), POINTER :: pGrid 
    TYPE(t_patch), POINTER :: pPatch             
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateFaceDist',&
  'RFLU_ModGeometry.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating face distance...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyFaceDist(pRegion)

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    ALLOCATE(pGrid%cofgDist(2,pGrid%nFaces),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%cofgDist')
    END IF ! global%error

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)  

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN
        WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'Patch:',iPatch
      END IF ! global%verbLevel  

      ALLOCATE(pPatch%cofgDist(pPatch%nBFaces),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%cofgDist')
      END IF ! global%error
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Creating face distance done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_CreateFaceDist
  
  




  
! ******************************************************************************
!
! Purpose: Create geometry.
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

  SUBROUTINE RFLU_CreateGeometry(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: errorFlag,ibv,ic,ifc,iPatch
    TYPE(t_grid), POINTER :: pGrid      
    TYPE(t_patch), POINTER :: pPatch            
    TYPE(t_global), POINTER :: global      

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_CreateGeometry',&
  'RFLU_ModGeometry.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating geometry...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyGeometry(pRegion)

! ******************************************************************************
!   Allocate memory for interior grid geometry
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Interior geometry...'
    END IF ! global%verbLevel

! ==============================================================================
!   Volume and volume centroids
! ==============================================================================

    ALLOCATE(pGrid%vol(pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%vol')
    END IF ! global%error

    ALLOCATE(pGrid%cofg(XCOORD:ZCOORD,pGrid%nCellsTot),STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%cofg')
    END IF ! global%error  

    DO ic = 1,pGrid%nCellsTot ! Explicit loop because of ASCI White problem
      pGrid%vol(ic)         = 0.0_RFREAL  
      pGrid%cofg(XCOORD,ic) = 0.0_RFREAL 
      pGrid%cofg(YCOORD,ic) = 0.0_RFREAL 
      pGrid%cofg(ZCOORD,ic) = 0.0_RFREAL                               
    END DO ! ic

! ==============================================================================
!   Face normals and centroids
! ==============================================================================

    ALLOCATE(pGrid%fn(XCOORD:XYZMAG,pGrid%nFacesTot),STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%fn')
    END IF ! global%error

    ALLOCATE(pGrid%fc(XCOORD:ZCOORD,pGrid%nFacesTot),STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%fc')
    END IF ! global%error

    DO ifc = 1,pGrid%nFacesTot 
      pGrid%fn(XCOORD,ifc) = 0.0_RFREAL
      pGrid%fn(YCOORD,ifc) = 0.0_RFREAL
      pGrid%fn(ZCOORD,ifc) = 0.0_RFREAL
      pGrid%fn(XYZMAG,ifc) = 0.0_RFREAL

      pGrid%fc(XCOORD,ifc) = 0.0_RFREAL
      pGrid%fc(YCOORD,ifc) = 0.0_RFREAL
      pGrid%fc(ZCOORD,ifc) = 0.0_RFREAL                
    END DO ! ifc 

! ******************************************************************************
!   Allocate memory for patch geometry
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Patch geometry...'
    END IF ! global%verbLevel  

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)  

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN
        WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'Patch:',iPatch
      END IF ! global%verbLevel  

! ==============================================================================
!     Face normal and centroid  
! ==============================================================================

      ALLOCATE(pPatch%fn(XCOORD:XYZMAG,pPatch%nBFacesTot),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%patches%fn')
      END IF ! global%error

      ALLOCATE(pPatch%fc(XCOORD:ZCOORD,pPatch%nBFacesTot),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%patches%fc')
      END IF ! global%error 

      DO ifc = 1,pPatch%nBFacesTot 
        pPatch%fn(XCOORD,ifc) = 0.0_RFREAL
        pPatch%fn(YCOORD,ifc) = 0.0_RFREAL
        pPatch%fn(ZCOORD,ifc) = 0.0_RFREAL
        pPatch%fn(XYZMAG,ifc) = 0.0_RFREAL

        pPatch%fc(XCOORD,ifc) = 0.0_RFREAL 
        pPatch%fc(YCOORD,ifc) = 0.0_RFREAL 
        pPatch%fc(ZCOORD,ifc) = 0.0_RFREAL                               
      END DO ! ifc         

! ==============================================================================
!     Vertex normals 
! ==============================================================================

      IF ( pRegion%mixtInput%moveGrid .EQV. .TRUE. ) THEN 
        ALLOCATE(pPatch%bvn(XCOORD:ZCOORD,pPatch%nBVertTot),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bvn')
        END IF ! global%error

        DO ibv = 1,pPatch%nBVertTot
          pPatch%bvn(XCOORD,ibv) = 0.0_RFREAL
          pPatch%bvn(YCOORD,ibv) = 0.0_RFREAL
          pPatch%bvn(ZCOORD,ibv) = 0.0_RFREAL                        
        END DO ! ibv
      ELSE
        NULLIFY(pPatch%bvn)
      END IF ! pRegion%mixtInput%moveGrid        
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************  

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating geometry done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_CreateGeometry








! ******************************************************************************
!
! Purpose: Destroy approximate cell centroids
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

  SUBROUTINE RFLU_DestroyApproxCentroids(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
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

    CALL RegisterFunction(global,'RFLU_DestroyApproxCentroids',&
  'RFLU_ModGeometry.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Destroying approximate ', & 
                                 'centroids...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel
    
    pGrid => pRegion%grid

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    DEALLOCATE(pGrid%cofgApp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%cofgApp')
    END IF ! global%error

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyApproxCentroids(pRegion)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Destroying approximate '// & 
                                 'centroids done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_DestroyApproxCentroids






! ******************************************************************************
!
! Purpose: Destroy distance from cell centroid to face centroid.
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

  SUBROUTINE RFLU_DestroyFaceDist(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: errorFlag,iPatch
    TYPE(t_grid), POINTER :: pGrid 
    TYPE(t_patch), POINTER :: pPatch             
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyFaceDist',&
  'RFLU_ModGeometry.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying face distance...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Deallocate memory
! ******************************************************************************

    DEALLOCATE(pGrid%cofgDist,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%cofgDist')
    END IF ! global%error

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)  

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel >= VERBOSE_HIGH ) THEN
        WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'Patch:',iPatch
      END IF ! global%verbLevel  

      DEALLOCATE(pPatch%cofgDist,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%cofgDist')
      END IF ! global%error
    END DO ! iPatch
    
! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyFaceDist(pRegion)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Destroying face distance done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_DestroyFaceDist





! ******************************************************************************
!
! Purpose: Destroy geometry.
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

  SUBROUTINE RFLU_DestroyGeometry(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: errorFlag,iPatch
    TYPE(t_grid), POINTER :: pGrid      
    TYPE(t_patch), POINTER :: pPatch            
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_DestroyGeometry',&
  'RFLU_ModGeometry.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying geometry...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

    pGrid => pRegion%grid   

! ******************************************************************************
!   Deallocate memory for interior grid geometry
! ******************************************************************************

    DEALLOCATE(pGrid%vol,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%grid%vol')
    END IF ! global%error

    DEALLOCATE(pGrid%fn,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%grid%fn')
    END IF ! global%error

    IF ( ASSOCIATED(pGrid%fc) .EQV. .TRUE. ) THEN 
      DEALLOCATE(pGrid%fc,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%grid%fc')
      END IF ! global%error
    END IF ! ASSOCIATED

    IF ( ASSOCIATED(pGrid%cofg) .EQV. .TRUE. ) THEN 
      DEALLOCATE(pGrid%cofg,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%grid%cofg')
      END IF ! global%error
    END IF ! ASSOCIATED 

! ******************************************************************************
!   Deallocate memory for patch geometry
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)  

      DEALLOCATE(pPatch%fn,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%patches%fn')
      END IF ! global%error

      IF ( ASSOCIATED(pPatch%fc) .EQV. .TRUE. ) THEN
        DEALLOCATE(pPatch%fc,STAT=errorFlag)
        global%error = errorFlag   
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'region%patches%fc')
        END IF ! global%error
      END IF ! ASSOCIATED

      IF ( ASSOCIATED(pPatch%bvn) .EQV. .TRUE. ) THEN 
        DEALLOCATE(pPatch%bvn,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%bvn')
        END IF ! global%error 
      END IF ! ASSOCIATED                                     
    END DO ! iPatch

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    CALL RFLU_NullifyGeometry(pRegion)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying geometry done.'
    END IF ! global%verbLevel 

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_DestroyGeometry





! ******************************************************************************
!
! Purpose: Nullify approximate cell centroids.
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

  SUBROUTINE RFLU_NullifyApproxCentroids(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
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

    CALL RegisterFunction(global,'RFLU_NullifyApproxCentroids',&
  'RFLU_ModGeometry.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Nullifying approximate centroids...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    NULLIFY(pGrid%cofgApp)

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME, & 
                                 'Nullifying approximate centroids done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_NullifyApproxCentroids






! ******************************************************************************
!
! Purpose: Nullify face distance.
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

  SUBROUTINE RFLU_NullifyFaceDist(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
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

    CALL RegisterFunction(global,'RFLU_NullifyFaceDist',&
  'RFLU_ModGeometry.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Nullifying face distance...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory
! ******************************************************************************

    NULLIFY(pGrid%cofgDist)

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)  

      NULLIFY(pPatch%cofgDist)
    END DO ! iPatch
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME, & 
                                 'Nullifying face distance done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_NullifyFaceDist






! ******************************************************************************
!
! Purpose: Nullify geometry.
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

  SUBROUTINE RFLU_NullifyGeometry(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
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

    CALL RegisterFunction(global,'RFLU_NullifyGeometry',&
  'RFLU_ModGeometry.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Nullifying geometry...'
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', &
                                       pRegion%iRegionGlobal
    END IF ! global%verbLevel

    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify memory for interior grid geometry
! ******************************************************************************

    NULLIFY(pGrid%vol)
    NULLIFY(pGrid%cofg)

    NULLIFY(pGrid%fn)
    NULLIFY(pGrid%fc)

! ******************************************************************************
!   Nullify memory for patch geometry
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)  

      NULLIFY(pPatch%fn)
      NULLIFY(pPatch%fc)
      NULLIFY(pPatch%bvn)
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************  

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Nullifying geometry done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_NullifyGeometry



! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModGeometry


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModGeometry.F90,v $
! Revision 1.35  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.34  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.33  2007/07/08 21:45:03  gzheng
! changed the PRESENT is used for PGI compiler
!
! Revision 1.32  2007/04/14 15:54:36  mtcampbe
! Sudden death fix
!
! Revision 1.31  2006/04/13 18:08:17  haselbac
! Cosmetics only
!
! Revision 1.30  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.29  2006/03/25 21:52:57  haselbac
! Substantial changes because of sype patches
!
! Revision 1.28  2005/01/13 21:42:19  haselbac
! Bug fix in computation of face distances
!
! Revision 1.27  2004/12/19 15:47:36  haselbac
! Added routines for face distance
!
! Revision 1.26  2004/10/19 19:27:58  haselbac
! Substantial clean-up
!
! Revision 1.25  2004/01/22 16:03:59  haselbac
! Made contents of modules PRIVATE, only procs PUBLIC, to avoid errors on ALC 
! and titan
!
! Revision 1.24  2003/11/25 21:03:24  haselbac
! Added additional output when checking data structure
!
! Revision 1.23  2003/09/02 02:52:32  haselbac
! Added value of error to error messages
!
! Revision 1.22  2003/08/20 02:09:58  haselbac
! Changed verbosity conditions to reduce solver output in GENx runs
!
! Revision 1.21  2003/06/04 22:08:30  haselbac
! Added Nullify routines, some cosmetics
!
! Revision 1.20  2003/05/06 01:04:02  haselbac
! Added NULLIFY for pPatch%bvn
!
! Revision 1.19  2003/05/02 02:34:41  haselbac
! Changed if-statements around volume ratio output (f90 quirk)
!
! Revision 1.18  2003/04/16 18:34:23  haselbac
! Removed commented-out statement
!
! Revision 1.17  2003/04/07 14:25:14  haselbac
! Changed volErr definition and limit
!
! Revision 1.16  2003/03/22 00:01:35  haselbac
! Fixed bug: Init in wrong places
!
! Revision 1.15  2003/03/21 23:09:23  haselbac
! Fixed subtle bug: Added init to RFLU_CreateGeometry
!
! Revision 1.14  2003/03/15 18:08:54  haselbac
! Bug fixes for centroid, || comps, other changes
!
! Revision 1.13  2003/01/28 16:31:18  haselbac
! Now contains BuildBVertexNormals, print out locs of extrema in volumes, 
! bug fixes, scaled coordinates
!
! Revision 1.12  2002/12/20 23:18:37  haselbac
! Fixed output bug: no output for verbosity=0
!
! Revision 1.11  2002/11/08 21:27:17  haselbac
! Write out volume ratio for moving-grid cases
!
! Revision 1.10  2002/10/27 19:05:57  haselbac
! Bug fixes (removed HACK_PERIODIC distinctions), changed creation of 
! geometry
!
! Revision 1.9  2002/10/16 21:13:34  haselbac
! Increased VOL_ERR_LIMIT to 1.0 (X29-new)
!
! Revision 1.8  2002/10/08 15:49:21  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.7  2002/10/05 19:04:39  haselbac
! Increased volume tolerance to 0.5, more CHECK output
!
! Revision 1.6  2002/09/09 15:07:22  haselbac
! global now under regions, bug fix for periodic hack
!
! Revision 1.5  2002/07/25 15:01:35  haselbac
! Deallocation changed for OLES, more checking output and error checks, 
! verbosity levels changed
!
! Revision 1.4  2002/06/27 15:50:34  haselbac
! Modifications for parallelization, added CHECK_DATASTRUCT flag
!
! Revision 1.3  2002/06/17 13:39:45  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.2  2002/05/04 16:39:51  haselbac
! New modules, more checking, deallocate fc for 1st order
!
! Revision 1.1  2002/04/11 18:48:48  haselbac
! Initial revision
!
! ******************************************************************************



















