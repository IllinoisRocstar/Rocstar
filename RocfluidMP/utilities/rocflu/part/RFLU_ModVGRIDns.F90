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
! Purpose: Collection of routines to read and convert MESH3D grids.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModVGRIDns.F90,v 1.5 2008/12/06 08:45:03 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModVGRIDns

  USE ModParameters
  USE ModDataTypes  
  USE ModError  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModMPI

  USE ModBuildFileNames, ONLY: BuildFileNamePlain
  
  USE RFLU_ModDimensions, ONLY: RFLU_SetMaxDimension  

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  PRIVATE

! ==============================================================================
! Data
! ==============================================================================

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModVGRIDns.F90,v $ $Revision: 1.5 $'        

  TYPE t_gridVGRIDns
    INTEGER :: nBQuads,nBTris,nMappings,nPatches
    INTEGER, DIMENSION(:), POINTER :: bTri2p
    INTEGER, DIMENSION(:,:), POINTER :: bTri2v,patch2bc
  END TYPE t_gridVGRIDns
  
  TYPE(t_gridVGRIDns) :: gridVGRIDns

! ==============================================================================
! Public procedures
! ==============================================================================

  PUBLIC :: RFLU_ConvVGRIDns2ROCFLU, & 
            RFLU_ReadGridVGRIDns

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
  


! ******************************************************************************
!
! Purpose: Convert grid format from VGRIDns to ROCFLU.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Read triangles as 132 to reverse normal vector (inwards in VGRIDns)
!
! ******************************************************************************

  SUBROUTINE RFLU_ConvVGRIDns2ROCFLU(pRegion)
 
    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
 
! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: errorFlag,i,iBegMax,iBegMin,iBeg1,iBeg2,ic,iEndMax, &
               iEndMin,iEnd1,iEnd2,iFile,ifl,iMap,iMap2,iPatch,it,term
    INTEGER, DIMENSION(:), ALLOCATABLE :: cntr
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch  
    TYPE(t_global), POINTER :: global

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion 
  
! ******************************************************************************
!   Start
! ******************************************************************************  

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ConvVGRIDns2ROCFLU', &
                          'RFLU_ModVGRIDns.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Converting from VGRIDns to ROCFLU format...'
    END IF ! global%verbLevel

! ==============================================================================
!   Set grid pointer and initialize variables
! ==============================================================================

    pGrid => pRegion%grid

    pGrid%nEdges    = 0
    pGrid%nEdgesTot = 0

    pGrid%nFaces    = 0
    pGrid%nFacesTot = 0
  
! ==============================================================================
!   Read additional info on grid, required for mapping of patches
! ==============================================================================

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Reading VGRIDns information file...'
    END IF ! global%verbLevel

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.vgi',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag )
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error  

    READ(iFile,*) pGrid%nPatches

    ALLOCATE(pRegion%patches(pGrid%nPatches),STAT=errorFlag)
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%patches')
    END IF ! global%error     

    READ(iFile,*) gridVGRIDns%nMappings

    ALLOCATE(gridVGRIDns%patch2bc(3,gridVGRIDns%nMappings), STAT=errorFlag )
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridVGRIDns%patch2bc')
    END IF ! global%error   

    DO iMap = 1,gridVGRIDns%nMappings
      READ(iFile,*) (gridVGRIDns%patch2bc(i,iMap),i=1,3)
    END DO ! iMap

! ==============================================================================
!   Check for consistent input - somewhat complicated...
! ============================================================================== 

    IF ( global%checkLevel > CHECK_NONE ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN  
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Checking patch mapping entries...'
      END IF ! global%verbLevel

      IF ( MINVAL(gridVGRIDns%patch2bc(1,:)) /= 1 .OR. & 
           MAXVAL(gridVGRIDns%patch2bc(2,:)) /= gridVGRIDns%nPatches ) THEN
        IF ( global%verbLevel > VERBOSE_NONE ) THEN           
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel    
        CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
      END IF ! MINVAL    

      DO iMap = 1,gridVGRIDns%nMappings
        IF ( gridVGRIDns%patch2bc(2,iMap) < gridVGRIDns%patch2bc(1,iMap) ) THEN 
          IF ( global%verbLevel > VERBOSE_NONE ) THEN
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
          END IF ! global%verbLevel                        
          CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
        END IF ! gridVGRIDns
      END DO ! iMap   

      DO iMap = 1,gridVGRIDns%nMappings
        IF ( MINVAL(gridVGRIDns%patch2bc(3,:)) /= 1 .OR. & 
             MAXVAL(gridVGRIDns%patch2bc(3,:)) /= pGrid%nPatches ) THEN 
          IF ( global%verbLevel > VERBOSE_NONE ) THEN
            WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
          END IF ! global%verbLevel     
          CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
        END IF ! gridVGRIDns
      END DO ! iMap       

      DO iMap = 1,gridVGRIDns%nMappings
        DO iMap2 = 1,gridVGRIDns%nMappings

          IF ( iMap /= iMap2 ) THEN 
            iBeg1 = gridVGRIDns%patch2bc(1,iMap)
            iEnd1 = gridVGRIDns%patch2bc(2,iMap)

            iBeg2 = gridVGRIDns%patch2bc(1,iMap2)
            iEnd2 = gridVGRIDns%patch2bc(2,iMap2)        

            IF ( iBeg1 < iBeg2 ) THEN 
              iBegMin = iBeg1
              iEndMin = iEnd1
              iBegMax = iBeg2
              iEndMax = iEnd2 
            ELSE IF ( iBeg1 > iBeg2 ) THEN 
              iBegMin = iBeg2
              iEndMin = iEnd2
              iBegMax = iBeg1
              iEndMax = iEnd1         
            ELSE ! iBeg1 and iBeg2 have the same value
              IF ( global%verbLevel > VERBOSE_NONE ) THEN
                WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
              END IF ! global%verbLevel  
              CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
            END IF ! iBeg1

            IF ( iEndMax <= iEndMin .OR. iEndMin >= iBegMax ) THEN
              IF ( global%verbLevel > VERBOSE_NONE ) THEN
                WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
              END IF ! global%verbLevel  
              CALL ErrorStop(global,ERR_PATCH_NUMBERING,__LINE__)
            END IF ! iEndMax
          END IF ! iMap

        END DO ! iMap2
      END DO ! iMap

      IF ( global%verbLevel > VERBOSE_NONE ) THEN  
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
                                 'Checking patch mapping entries done.'
      END IF ! global%verbLevel
    END IF ! global%checkLevel

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
                               'Reading VGRIDns information file done.'
    END IF ! global%verbLevel
  
! ******************************************************************************
!   Convert to ROCFLU format
! ******************************************************************************
 
    global%nPatches = pGrid%nPatches 
 
! ==============================================================================
!   Generate boundary triangle lists
! ============================================================================== 

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
                               'Generating boundary triangle lists...'  
    END IF ! global%verbLevel  

! ------------------------------------------------------------------------------
!   Count number of boundary faces
! ------------------------------------------------------------------------------

    ALLOCATE(cntr(pGrid%nPatches),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'cntr')
    END IF ! global%error   

    cntr(:) = 0 

    DO it = 1,gridVGRIDns%nBTris
      DO iMap = 1,gridVGRIDns%nMappings
        IF ( gridVGRIDns%bTri2p(it) >= gridVGRIDns%patch2bc(1,iMap) .AND. & 
             gridVGRIDns%bTri2p(it) <= gridVGRIDns%patch2bc(2,iMap) ) THEN 
          iPatch = gridVGRIDns%patch2bc(3,iMap)

          cntr(iPatch) = cntr(iPatch) + 1
          EXIT
        END IF ! gridVGRIDns
      END DO ! iMap
    END DO ! it

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)  

      pPatch%nBTris  = cntr(iPatch)
      pPatch%nBQuads = 0 
      pPatch%nBVert  = 0  

      pPatch%iPatchGlobal = iPatch
      pPatch%iBorder      = PATCH_IBORDER_DEFAULT     
      pPatch%renumFlag    = .FALSE.

      IF ( iPatch > 1 ) THEN 
        cntr(iPatch) = cntr(iPatch) + cntr(iPatch-1)
      END IF ! iPatch
    END DO ! iPatch

    IF ( cntr(pGrid%nPatches) /= gridVGRIDns%nBTris ) THEN 
      CALL ErrorStop(global,ERR_NBFACES_WRONG,__LINE__)   
    END IF ! grid  

! ------------------------------------------------------------------------------
!   Set total boundary patch quantities and number of boundary faces
! ------------------------------------------------------------------------------

    pGrid%nBFaces = 0 
    
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      pPatch%nBFaces = pPatch%nBTris + pPatch%nBQuads
      pGrid%nBFaces  = pGrid%nBFaces + pPatch%nBFaces

      pPatch%nBFacesTot = pPatch%nBFaces
      pPatch%nBQuadsTot = pPatch%nBQuads  
      pPatch%nBTrisTot  = pPatch%nBTris    
      pPatch%nBVertTot  = pPatch%nBVert  
      
      pPatch%nBTrisMax  = RFLU_SetMaxDimension(global,pPatch%nBTrisTot)
      pPatch%nBQuadsMax = RFLU_SetMaxDimension(global,pPatch%nBQuadsTot) 
      pPatch%nBFacesMax = RFLU_SetMaxDimension(global,pPatch%nBFacesTot)             
      pPatch%nBVertMax  = RFLU_SetMaxDimension(global,pPatch%nBVertTot)
      
      pPatch%nBCellsVirt = 0                    
    END DO ! iPatch  
    
    pGrid%nBFacesTot = pGrid%nBFaces    

! ------------------------------------------------------------------------------
!   Build face lists
! ------------------------------------------------------------------------------

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)    

      ALLOCATE(pPatch%bTri2v(3,pPatch%nBTrisMax),STAT=errorFlag)
      global%error = errorFlag 
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bTri2v')
      END IF ! global%error 
    END DO ! iPatch  

    cntr(:) = 0  

    DO it = 1,gridVGRIDns%nBTris
      DO iMap = 1,gridVGRIDns%nMappings
        IF ( gridVGRIDns%bTri2p(it) >= gridVGRIDns%patch2bc(1,iMap) .AND. & 
             gridVGRIDns%bTri2p(it) <= gridVGRIDns%patch2bc(2,iMap) ) THEN   
          iPatch = gridVGRIDns%patch2bc(3,iMap)
          pPatch => pRegion%patches(iPatch)           

          cntr(iPatch) = cntr(iPatch) + 1

          pPatch%bTri2v(1,cntr(iPatch)) = gridVGRIDns%bTri2v(1,it)
          pPatch%bTri2v(2,cntr(iPatch)) = gridVGRIDns%bTri2v(3,it)
          pPatch%bTri2v(3,cntr(iPatch)) = gridVGRIDns%bTri2v(2,it)        

          EXIT
        END IF ! gridVGRIDns
      END DO ! iMap        
    END DO ! it    

    DEALLOCATE(cntr,STAT=errorFlag)
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'cntr')
    END IF ! global%error

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, & 
                               'Generating boundary triangle lists done.'  
    END IF ! global%verbLevel  

! ******************************************************************************
!   Convert tetrahedra to ROCFLU format - NOTE renumbering!
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Renumbering tetrahedra...'  
    END IF ! global%verbLevel 

    DO ic = 1,pGrid%nTetsTot
      term              = pGrid%tet2v(2,ic)
      pGrid%tet2v(2,ic) = pGrid%tet2v(3,ic) 
      pGrid%tet2v(3,ic) = term
    END DO ! iv

    IF ( global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Renumbering tetrahedra done.'  
    END IF ! global%verbLevel 

! ==============================================================================
!   Initialize other cell types
! ==============================================================================

    pGrid%nHexsTot = 0
    pGrid%nPrisTot = 0
    pGrid%nPyrsTot = 0

    pGrid%nHexs = pGrid%nHexsTot
    pGrid%nPris = pGrid%nPrisTot
    pGrid%nPyrs = pGrid%nPyrsTot
  
! ******************************************************************************
!   Deallocate VGRIDns memory
! ******************************************************************************

    DEALLOCATE(gridVGRIDns%bTri2v,STAT=errorFlag)
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gridVGRIDns%bTri2v')
    END IF ! global%error

    DEALLOCATE(gridVGRIDns%bTri2p,STAT=errorFlag)
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gridVGRIDns%bTri2p')
    END IF ! global%error

! ******************************************************************************
!   Allocate memory for boundary face lists bf2c and bf2v
! ******************************************************************************

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      ALLOCATE(pPatch%bf2c(pPatch%nBFacesMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2c')
      END IF ! global%error

      ALLOCATE(pPatch%bf2v(4,pPatch%nBFacesMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bf2v')
      END IF ! global%error       

      DO ifl = 1,pPatch%nBFacesMax
        pPatch%bf2v(1,ifl) = VERT_NONE 
        pPatch%bf2v(2,ifl) = VERT_NONE
        pPatch%bf2v(3,ifl) = VERT_NONE
        pPatch%bf2v(4,ifl) = VERT_NONE
      END DO ! ifl  
    END DO ! iPatch   

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Converting from VGRIDns to ROCFLU format done.'
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ConvVGRIDns2ROCFLU  
  
  
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Read grid file from VGRIDns.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. .mapbc file not read - information contained in .mapbc file is 
!      found in more convenient way in user-provided .vgi file. See routine
!      RFLU_ConvVGRIDns2ROCFLU.F90.
!   2. Only deals with tetrahedra!
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadGridVGRIDns(pRegion)
  
  USE ModBuildFileNames, ONLY: BuildFileNamePlain  
  
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
  
    CHARACTER(CHRLEN) :: dummyString,iFileName  
    INTEGER :: cvmax,cvmin,dummyInteger,errorFlag,i,ic,iFile,it,iv
    INTEGER, DIMENSION(:), ALLOCATABLE :: cntr   
    REAL(RFREAL) :: dummyReal
    TYPE(t_grid), POINTER :: pGrid    
    TYPE(t_global), POINTER :: global  
        
! ******************************************************************************
!   Start
! ******************************************************************************
  
    gridVGRIDns%nBTris  = 0
    gridVGRIDns%nBQuads = 0

! ******************************************************************************
!   Read bc file
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadGridVGRIDns', &
                          'RFLU_ModVGRIDns.F90')

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading VGRIDns grid file...'
    END IF ! global%verbLevel

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.vbc',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error  

! ==============================================================================
!   General information
! ==============================================================================

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'General information...'
    END IF ! global%verbLevel

    READ(iFile,*) gridVGRIDns%nBTris,dummyString,gridVGRIDns%nPatches, & 
                  dummyInteger
    READ(iFile,'(A)') dummyString

! ==============================================================================
!   Boundary faces
! ==============================================================================

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary faces...'
    END IF ! global%verbLevel

    ALLOCATE(gridVGRIDns%bTri2p(gridVGRIDns%nBTris),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridVGRIDns%bTri2p')
    END IF ! global%error

    ALLOCATE(gridVGRIDns%bTri2v(3,gridVGRIDns%nBTris),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gridVGRIDns%bTri2v')
    END IF ! global%error

    DO it = 1,gridVGRIDns%nBTris
      READ(iFile,*) dummyInteger,gridVGRIDns%bTri2p(it), & 
                    (gridVGRIDns%bTri2v(iv,it),iv=1,3)
    END DO ! it

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'File read successfully.'  
    END IF ! global%verbLevel

! ******************************************************************************
!   Read connectivity and coordinate file
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading VGRIDns .cogsg file...'
    END IF ! global%verbLevel

    iFile = IF_GRID

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.cogsg',iFileName)

    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Connectivity
! ==============================================================================

    pGrid => pRegion%grid

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header...'
    END IF ! global%verbLevel

    READ(iFile) dummyInteger,pGrid%nTetsTot,pGrid%nVertTot
    REWIND(iFile)

    pGrid%nTetsMax = RFLU_SetMaxDimension(global,pGrid%nTetsTot)
    pGrid%nVertMax = RFLU_SetMaxDimension(global,pGrid%nVertTot)

! ==============================================================================
!   Connectivity
! ==============================================================================

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Connectivity...'
    END IF ! global%verbLevel

    ALLOCATE(pGrid%tet2v(4,pGrid%nTetsMax),STAT=errorFlag)
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%tet2v')
    END IF ! global%error

    READ(iFile) dummyInteger,dummyInteger,dummyInteger, & 
                dummyInteger,dummyInteger,dummyInteger, & 
                dummyReal,((pGrid%tet2v(iv,ic),ic=1,pGrid%nTetsTot),iv=1,4)

    pGrid%nCellsTot = pGrid%nTetsTot
    pGrid%nCellsMax = RFLU_SetMaxDimension(global,pGrid%nCellsTot)    

! ------------------------------------------------------------------------------
!   Checking - only valid for tetrahedral grids, cf. RFLU_ReadGridCENTAUR.F90
! ------------------------------------------------------------------------------

    IF ( global%checkLevel > CHECK_NONE ) THEN
      IF ( global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                                 'Checking face connectivity array entries...'
      END IF ! global%verbLevel

      cvmin = MINVAL(gridVGRIDns%bTri2v(1:3,1:gridVGRIDns%nBTris))
      cvmax = MAXVAL(gridVGRIDns%bTri2v(1:3,1:gridVGRIDns%nBTris))

      IF ( cvmin < 1 .OR. cvmax > pGrid%nVertTot ) THEN
        IF ( global%verbLevel > VERBOSE_NONE ) THEN       
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Check failed.'
        END IF ! global%verbLevel        
        CALL ErrorStop(global,ERR_VERTEX_NUMBER,__LINE__)
      END IF ! cvmin
            
      IF ( global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
              'Checking face connectivity array entries done.'
      END IF ! global%verbLevel

      IF ( global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
                                 'Checking cell connectivity array entries...'
      END IF ! global%verbLevel 

      cvmin = MINVAL(pGrid%tet2v(1:4,1:pGrid%nTetsTot))
      cvmax = MAXVAL(pGrid%tet2v(1:4,1:pGrid%nTetsTot))    

      IF ( cvmin /= 1 .OR. cvmax /= pGrid%nVertTot ) THEN
        IF ( global%verbLevel > VERBOSE_NONE ) THEN       
          WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'CHECK failed.'
        END IF ! global%verbLevel        
        CALL ErrorStop(global,ERR_VERTEX_NUMBER,__LINE__)
      END IF ! cvmin

      IF ( global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME, &
              'Checking cell connectivity array entries done.'
      END IF ! global%verbLevel 
    END IF ! global%checkLevel

! ==============================================================================
!   Coordinates
! ==============================================================================

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Coordinates...'
    END IF ! global%verbLevel

    ALLOCATE(pGrid%xyz(3,pGrid%nVertMax),STAT=errorFlag)
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%xyz')
    END IF ! global%error

    READ(iFile) ((pGrid%xyz(i,iv),iv=1,pGrid%nVertTot),i=1,3) 

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error
  
! ******************************************************************************
!   Set grid size variables
! ******************************************************************************

    pGrid%nVert  = pGrid%nVertTot
    pGrid%nCells = pGrid%nCellsTot
    pGrid%nTets  = pGrid%nTetsTot
  
! ******************************************************************************
!   Print grid statistics
! ******************************************************************************

    IF ( global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,3X,A)')       SOLVER_NAME,'Grid Statistics:'
      WRITE(STDOUT,'(A,5X,A,2X,I9)') SOLVER_NAME,'Vertices:           ', & 
                                     pGrid%nVertTot
      WRITE(STDOUT,'(A,5X,A,2X,I9)') SOLVER_NAME,'Tetrahedra:         ', & 
                                     pGrid%nTetsTot  
      WRITE(STDOUT,'(A,5X,A,2X,I9)') SOLVER_NAME,'Boundary triangles: ', & 
                                     gridVGRIDns%nBTris
    END IF ! global%verbLevel      
      
! ******************************************************************************
!   End
! ******************************************************************************
  
    IF ( global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading VGRIDns grid file done.' 
    END IF ! global%verbLevel            

    CALL DeregisterFunction(global)    
  
  END SUBROUTINE RFLU_ReadGridVGRIDns

  
  
  
  
  
  
! ******************************************************************************
! End
! ******************************************************************************


END MODULE RFLU_ModVGRIDns

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModVGRIDns.F90,v $
! Revision 1.5  2008/12/06 08:45:03  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:18:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/18 21:02:38  haselbac
! Fixed checks; bug arose bcos of max-dimensioned arrays for sype
!
! Revision 1.2  2006/03/25 22:04:29  haselbac
! Changes because of sype patches
!
! Revision 1.1  2005/04/15 15:09:14  haselbac
! Initial revision
!
! Revision 1.4  2005/01/20 14:54:56  haselbac
! Added setting of nBFaces and nBFacesTot
!
! Revision 1.3  2004/11/03 17:09:24  haselbac
! Removed setting of vertex and cell flags
!
! Revision 1.2  2004/10/19 19:31:15  haselbac
! Removed renumbering of bface lists
!
! Revision 1.1  2004/07/06 15:15:50  haselbac
! Initial revision
!
! ******************************************************************************
  

  








