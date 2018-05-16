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
! Purpose: Collection of utility routines for writing ENSIGHT files.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModENSIGHTUtils.F90,v 1.4 2008/12/06 08:45:06 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModENSIGHTUtils

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PARAMETER, PRIVATE :: &
    RCSIdentString = '$RCSfile: RFLU_ModENSIGHTUtils.F90,v $ $Revision: 1.4 $'

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: RFLU_ENS_WriteGrid, & 
            RFLU_ENS_WriteScalar, & 
            RFLU_ENS_WriteVector


! ******************************************************************************
! Subroutines and functions
! ******************************************************************************

  CONTAINS





! ******************************************************************************
!
! Purpose: Write grid to ENSIGHT geometry file.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   emptyPartFlag       Flag indicating whether part should be empty
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_WriteGrid(pRegion,emptyPartFlag)

  USE RFLU_ModRenumberList, ONLY: RFLU_DenumberList, & 
                                  RFLU_RenumberList
  USE RFLU_ModTopologyUtils, ONLY: RFLU_BuildConnVertList

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  LOGICAL, INTENT(IN) :: emptyPartFlag
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(80) :: dummyString
  INTEGER :: errorFlag,i,icl,ifl,iPatch,ivg,ivl,nBVert,nBVertEst
  INTEGER, DIMENSION(:), ALLOCATABLE :: vList
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ENS_WriteGrid', &
                        'RFLU_ModENSIGHTUtils.F90')

  pGrid => pRegion%grid

  global%postPartNumber = global%postPartNumber + 1

! ******************************************************************************
! Volume grid
! ******************************************************************************

  dummyString = 'part'
  WRITE(IF_ENS_GEOMETRY) dummyString
  WRITE(IF_ENS_GEOMETRY) global%postPartNumber
  
  WRITE(dummyString,'(A,I5.5)') 'VOL_',pRegion%iRegionGlobal
  WRITE(IF_ENS_GEOMETRY) dummyString     

  IF ( emptyPartFlag .EQV. .FALSE. ) THEN

! ==============================================================================
!   Coordinates
! ==============================================================================

    dummyString = 'coordinates'
    WRITE(IF_ENS_GEOMETRY) dummyString  
    WRITE(IF_ENS_GEOMETRY) pGrid%nVertTot
    WRITE(IF_ENS_GEOMETRY) (REAL(pGrid%xyz(XCOORD,ivg),KIND=SPREAL), &
                           ivg=1,pGrid%nVertTot)
    WRITE(IF_ENS_GEOMETRY) (REAL(pGrid%xyz(YCOORD,ivg),KIND=SPREAL), &
                           ivg=1,pGrid%nVertTot)
    WRITE(IF_ENS_GEOMETRY) (REAL(pGrid%xyz(ZCOORD,ivg),KIND=SPREAL), &
                           ivg=1,pGrid%nVertTot)         

! ==============================================================================
!   Cell connectivity
! ==============================================================================

! ------------------------------------------------------------------------------
!   Tetrahedra
! ------------------------------------------------------------------------------
  
    IF ( pGrid%nTets > 0 ) THEN 
      dummyString = 'tetra4'
      WRITE(IF_ENS_GEOMETRY) dummyString 
      WRITE(IF_ENS_GEOMETRY) pGrid%nTets
      WRITE(IF_ENS_GEOMETRY) ((pGrid%tet2v(ivl,icl),ivl=1,4),icl=1,pGrid%nTets) 
    END IF ! pGrid%nTets
    
    IF ( pGrid%nTetsTot > pGrid%nTets ) THEN 
      dummyString = 'g_tetra4'
      WRITE(IF_ENS_GEOMETRY) dummyString 
      WRITE(IF_ENS_GEOMETRY) pGrid%nTetsTot-pGrid%nTets
      WRITE(IF_ENS_GEOMETRY) ((pGrid%tet2v(ivl,icl),ivl=1,4), &
                             icl=pGrid%nTets+1,pGrid%nTetsTot)
    END IF ! pGrid%nTetsTot    

! ------------------------------------------------------------------------------
!   Hexahedra
! ------------------------------------------------------------------------------
    
    IF ( pGrid%nHexs > 0 ) THEN 
      dummyString = 'hexa8'
      WRITE(IF_ENS_GEOMETRY) dummyString
      WRITE(IF_ENS_GEOMETRY) pGrid%nHexs
      WRITE(IF_ENS_GEOMETRY) ((pGrid%hex2v(ivl,icl),ivl=1,8),icl=1,pGrid%nHexs) 
    END IF ! pGrid%nHexs
    
    IF ( pGrid%nHexsTot > pGrid%nHexs ) THEN 
      dummyString = 'g_hexa8'
      WRITE(IF_ENS_GEOMETRY) dummyString 
      WRITE(IF_ENS_GEOMETRY) pGrid%nHexsTot-pGrid%nHexs
      WRITE(IF_ENS_GEOMETRY) ((pGrid%hex2v(ivl,icl),ivl=1,8), &
                             icl=pGrid%nHexs+1,pGrid%nHexsTot)
    END IF ! pGrid%nHexsTot    
  
! ------------------------------------------------------------------------------
!   Prisms
! ------------------------------------------------------------------------------
  
    IF ( pGrid%nPris > 0 ) THEN 
      dummyString = 'penta6'
      WRITE(IF_ENS_GEOMETRY) dummyString 
      WRITE(IF_ENS_GEOMETRY) pGrid%nPris
      WRITE(IF_ENS_GEOMETRY) ((pGrid%pri2v(ivl,icl),ivl=1,6),icl=1,pGrid%nPris) 
    END IF ! pGrid%nPris 
    
    IF ( pGrid%nPrisTot > pGrid%nPris ) THEN 
      dummyString = 'g_penta6'
      WRITE(IF_ENS_GEOMETRY) dummyString 
      WRITE(IF_ENS_GEOMETRY) pGrid%nPrisTot-pGrid%nPris
      WRITE(IF_ENS_GEOMETRY) ((pGrid%pri2v(ivl,icl),ivl=1,6) , &
                             icl=pGrid%nPris+1,pGrid%nPrisTot)
    END IF ! pGrid%nPrisTot       

! ------------------------------------------------------------------------------
!   Pyramids
! ------------------------------------------------------------------------------
  
    IF ( pGrid%nPyrs > 0 ) THEN 
      dummyString = 'pyramid5'
      WRITE(IF_ENS_GEOMETRY) dummyString
      WRITE(IF_ENS_GEOMETRY) pGrid%nPyrs
      WRITE(IF_ENS_GEOMETRY) ((pGrid%pyr2v(ivl,icl),ivl=1,5),icl=1,pGrid%nPyrs) 
    END IF ! pGrid%nPyrs
    
    IF ( pGrid%nPyrsTot > pGrid%nPyrs ) THEN 
      dummyString = 'g_pyramid5'
      WRITE(IF_ENS_GEOMETRY) dummyString 
      WRITE(IF_ENS_GEOMETRY) pGrid%nPyrsTot-pGrid%nPyrs
      WRITE(IF_ENS_GEOMETRY) ((pGrid%pyr2v(ivl,icl),ivl=1,5), &
                             icl=pGrid%nPyrs+1,pGrid%nPyrsTot)
    END IF ! pGrid%nPyrsTot    
  END IF ! emptyPartFlag

! ******************************************************************************
! Surface grid
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
    
! ==============================================================================
!   Triangles
! ==============================================================================
   
! ------------------------------------------------------------------------------
!   Actual triangles 
! ------------------------------------------------------------------------------   
   
    IF ( pPatch%nBTris > 0 ) THEN 
      global%postPartNumber = global%postPartNumber + 1    
    
      dummyString = 'part'
      WRITE(IF_ENS_GEOMETRY) dummyString
      WRITE(IF_ENS_GEOMETRY) global%postPartNumber
      WRITE(dummyString,'(A,I3.3,A,I5.5)') 'PAT_',iPatch,'_TRI-A_', &
                                            pRegion%iRegionGlobal
      WRITE(IF_ENS_GEOMETRY) dummyString  

      IF ( emptyPartFlag .EQV. .FALSE. ) THEN
      
! ----- Build list of vertices -------------------------------------------------     
      
        nBVertEst = 3*pPatch%nBTris 

        ALLOCATE(vList(nBVertEst),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vList')
        END IF ! global%error

        CALL RFLU_BuildConnVertList(global,pPatch%bTri2v(1:3,1:pPatch%nBTris), & 
                                    3,pPatch%nBTris,vList,nBVertEst,nBVert)

! ----- Write coordinates ------------------------------------------------------     

        dummyString = 'coordinates'
        WRITE(IF_ENS_GEOMETRY) dummyString  
        WRITE(IF_ENS_GEOMETRY) nBVert
        WRITE(IF_ENS_GEOMETRY) (REAL(pGrid%xyz(XCOORD,vList(ivl)), & 
                               KIND=SPREAL),ivl=1,nBVert) 
        WRITE(IF_ENS_GEOMETRY) (REAL(pGrid%xyz(YCOORD,vList(ivl)), &
                               KIND=SPREAL),ivl=1,nBVert) 
        WRITE(IF_ENS_GEOMETRY) (REAL(pGrid%xyz(ZCOORD,vList(ivl)), &
                               KIND=SPREAL),ivl=1,nBVert)                 

! ----- Write renumbered connectivity ------------------------------------------

        dummyString = 'tria3'
        WRITE(IF_ENS_GEOMETRY) dummyString 
        WRITE(IF_ENS_GEOMETRY) pPatch%nBTris

        CALL RFLU_RenumberList(global,3,pPatch%nBTris, & 
                               pPatch%bTri2v(1:3,1:pPatch%nBTris),nBVert, & 
                               vList(1:nBVert))      

        WRITE(IF_ENS_GEOMETRY) ((pPatch%bTri2v(ivl,ifl),ivl=1,3), & 
                               ifl=1,pPatch%nBTris) 

        CALL RFLU_DenumberList(global,3,pPatch%nBTris, & 
                               pPatch%bTri2v(1:3,1:pPatch%nBTris),nBVert, & 
                               vList(1:nBVert))  

! ----- Destroy list of vertices -----------------------------------------------      

        DEALLOCATE(vList,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vList')
        END IF ! global%error
      END IF ! emptyPartFlag
    END IF ! pPatch%nBTris
    
! ------------------------------------------------------------------------------
!   Virtual triangles 
! ------------------------------------------------------------------------------   
    
    IF ( pPatch%nBTrisTot > pPatch%nBTris ) THEN 
      global%postPartNumber = global%postPartNumber + 1    
    
      dummyString = 'part'
      WRITE(IF_ENS_GEOMETRY) dummyString
      WRITE(IF_ENS_GEOMETRY) global%postPartNumber
      WRITE(dummyString,'(A,I3.3,A,I5.5)') 'PAT_',iPatch,'_TRI-V_', &
                                            pRegion%iRegionGlobal
      WRITE(IF_ENS_GEOMETRY) dummyString  

      IF ( emptyPartFlag .EQV. .FALSE. ) THEN

! ----- Build list of vertices -------------------------------------------------     

        nBVertEst = 3*(pPatch%nBTrisTot-pPatch%nBTris) 

        ALLOCATE(vList(nBVertEst),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vList')
        END IF ! global%error

        CALL RFLU_BuildConnVertList(global, &
             pPatch%bTri2v(1:3,pPatch%nBTris+1:pPatch%nBTrisTot), & 
             3,pPatch%nBTrisTot-pPatch%nBTris,vList,nBVertEst,nBVert)

! ----- Write coordinates ------------------------------------------------------    

        dummyString = 'coordinates'
        WRITE(IF_ENS_GEOMETRY) dummyString  
        WRITE(IF_ENS_GEOMETRY) nBVert
        WRITE(IF_ENS_GEOMETRY) (REAL(pGrid%xyz(XCOORD,vList(ivl)), & 
                               KIND=SPREAL),ivl=1,nBVert) 
        WRITE(IF_ENS_GEOMETRY) (REAL(pGrid%xyz(YCOORD,vList(ivl)), &
                               KIND=SPREAL),ivl=1,nBVert) 
        WRITE(IF_ENS_GEOMETRY) (REAL(pGrid%xyz(ZCOORD,vList(ivl)), &
                               KIND=SPREAL),ivl=1,nBVert) 

! ----- Write coordinates ------------------------------------------------------      

        dummyString = 'g_tria3'
        WRITE(IF_ENS_GEOMETRY) dummyString 
        WRITE(IF_ENS_GEOMETRY) pPatch%nBTrisTot-pPatch%nBTris

        CALL RFLU_RenumberList(global,3,pPatch%nBTrisTot-pPatch%nBTris, & 
             pPatch%bTri2v(1:3,pPatch%nBTris+1:pPatch%nBTrisTot),nBVert, & 
             vList(1:nBVert))      

        WRITE(IF_ENS_GEOMETRY) ((pPatch%bTri2v(ivl,ifl),ivl=1,3), & 
                               ifl=pPatch%nBTris+1,pPatch%nBTrisTot) 

        CALL RFLU_DenumberList(global,3,pPatch%nBTrisTot-pPatch%nBTris, & 
             pPatch%bTri2v(1:3,pPatch%nBTris+1:pPatch%nBTrisTot),nBVert, & 
             vList(1:nBVert))  

! ----- Destroy list of vertices -----------------------------------------------      

        DEALLOCATE(vList,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vList')
        END IF ! global%error
      END IF ! emptyPartFlag
    END IF ! pPatch%nBTrisTot   
    
! ==============================================================================
!   Quadrilaterals
! ==============================================================================
  
! ------------------------------------------------------------------------------
!   Actual quadrilaterals   
! ------------------------------------------------------------------------------
  
    IF ( pPatch%nBQuads > 0 ) THEN
      global%postPartNumber = global%postPartNumber + 1    
     
      dummyString = 'part'
      WRITE(IF_ENS_GEOMETRY) dummyString
      WRITE(IF_ENS_GEOMETRY) global%postPartNumber
      WRITE(dummyString,'(A,I3.3,A,I5.5)') 'PAT_',iPatch,'_QUAD-A_', &
                                            pRegion%iRegionGlobal
      WRITE(IF_ENS_GEOMETRY) dummyString   

      IF ( emptyPartFlag .EQV. .FALSE. ) THEN

! ----- Build list of vertices -------------------------------------------------     

        nBVertEst = 4*pPatch%nBQuads 

        ALLOCATE(vList(nBVertEst),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vList')
        END IF ! global%error

        CALL RFLU_BuildConnVertList(global, & 
                                    pPatch%bQuad2v(1:4,1:pPatch%nBQuads), & 
                                    4,pPatch%nBQuads,vList,nBVertEst,nBVert)

! ----- Write coordinates ------------------------------------------------------       

        dummyString = 'coordinates'
        WRITE(IF_ENS_GEOMETRY) dummyString  
        WRITE(IF_ENS_GEOMETRY) nBVert
        WRITE(IF_ENS_GEOMETRY) (REAL(pGrid%xyz(XCOORD,vList(ivl)), & 
                               KIND=SPREAL),ivl=1,nBVert) 
        WRITE(IF_ENS_GEOMETRY) (REAL(pGrid%xyz(YCOORD,vList(ivl)), &
                               KIND=SPREAL),ivl=1,nBVert) 
        WRITE(IF_ENS_GEOMETRY) (REAL(pGrid%xyz(ZCOORD,vList(ivl)), &
                               KIND=SPREAL),ivl=1,nBVert)     

! ----- Write coordinates       

        dummyString = 'quad4'
        WRITE(IF_ENS_GEOMETRY) dummyString
        WRITE(IF_ENS_GEOMETRY) pPatch%nBQuads

        CALL RFLU_RenumberList(global,4,pPatch%nBQuads, & 
                               pPatch%bQuad2v(1:4,1:pPatch%nBQuads),nBVert, & 
                               vList(1:nBVert)) 

        WRITE(IF_ENS_GEOMETRY) ((pPatch%bQuad2v(ivl,ifl),ivl=1,4) , &
                               ifl=1,pPatch%nBQuads)

        CALL RFLU_DenumberList(global,4,pPatch%nBQuads, & 
                               pPatch%bQuad2v(1:4,1:pPatch%nBQuads),nBVert, & 
                               vList(1:nBVert)) 

! ----- Destroy list of vertices -----------------------------------------------      

        DEALLOCATE(vList,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vList')
        END IF ! global%error
      END IF ! emptyPartFlag              
    END IF ! pPatch%nBQuads 
    
! ------------------------------------------------------------------------------
!   Virtual quadrilaterals    
! ------------------------------------------------------------------------------
    
    IF ( pPatch%nBQuadsTot > pPatch%nBQuads ) THEN
      global%postPartNumber = global%postPartNumber + 1    
     
      dummyString = 'part'
      WRITE(IF_ENS_GEOMETRY) dummyString
      WRITE(IF_ENS_GEOMETRY) global%postPartNumber
      WRITE(dummyString,'(A,I3.3,A,I5.5)') 'PAT_',iPatch,'_QUAD-V_', &
                                            pRegion%iRegionGlobal
      WRITE(IF_ENS_GEOMETRY) dummyString   

      IF ( emptyPartFlag .EQV. .FALSE. ) THEN

! ----- Build list of vertices -------------------------------------------------     

        nBVertEst = 4*(pPatch%nBQuadsTot-pPatch%nBQuads) 

        ALLOCATE(vList(nBVertEst),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vList')
        END IF ! global%error

        CALL RFLU_BuildConnVertList(global, & 
             pPatch%bQuad2v(1:4,pPatch%nBQuads+1:pPatch%nBQuadsTot), & 
             4,pPatch%nBQuadsTot-pPatch%nBQuads,vList,nBVertEst,nBVert)

! ----- Write coordinates ------------------------------------------------------      

        dummyString = 'coordinates'
        WRITE(IF_ENS_GEOMETRY) dummyString  
        WRITE(IF_ENS_GEOMETRY) nBVert
        WRITE(IF_ENS_GEOMETRY) (REAL(pGrid%xyz(XCOORD,vList(ivl)), & 
                               KIND=SPREAL),ivl=1,nBVert) 
        WRITE(IF_ENS_GEOMETRY) (REAL(pGrid%xyz(YCOORD,vList(ivl)), &
                               KIND=SPREAL),ivl=1,nBVert) 
        WRITE(IF_ENS_GEOMETRY) (REAL(pGrid%xyz(ZCOORD,vList(ivl)), &
                               KIND=SPREAL),ivl=1,nBVert)    

! ----- Write coordinates      

        dummyString = 'g_quad4'
        WRITE(IF_ENS_GEOMETRY) dummyString 
        WRITE(IF_ENS_GEOMETRY) pPatch%nBQuadsTot-pPatch%nBQuads

        CALL RFLU_RenumberList(global,4,pPatch%nBQuadsTot-pPatch%nBQuads, & 
             pPatch%bQuad2v(1:4,pPatch%nBQuads+1:pPatch%nBQuadsTot), &
             nBVert,vList(1:nBVert)) 

        WRITE(IF_ENS_GEOMETRY) ((pPatch%bQuad2v(ivl,ifl),ivl=1,4), &
                               ifl=pPatch%nBQuads+1,pPatch%nBQuadsTot) 

        CALL RFLU_DenumberList(global,4,pPatch%nBQuadsTot-pPatch%nBQuads, & 
             pPatch%bQuad2v(1:4,pPatch%nBQuads+1:pPatch%nBQuadsTot), &
             nBVert,vList(1:nBVert)) 

! ----- Destroy list of vertices -----------------------------------------------      

        DEALLOCATE(vList,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vList')
        END IF ! global%error
      END IF ! emptyPartFlag              
    END IF ! pPatch%nBQuadsTot
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ENS_WriteGrid








! ******************************************************************************
!
! Purpose: Write scalar variable to ENSIGHT file.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   var                 Pointer to scalar
!   iFile               File index
!   emptyPartFlag       Flag indicating whether part should be empty
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_WriteScalar(pRegion,var,iFile,emptyPartFlag)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  LOGICAL, INTENT(IN) :: emptyPartFlag
  INTEGER, INTENT(IN) :: iFile
  REAL(RFREAL), DIMENSION(:), POINTER :: var
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(80) :: dummyString
  INTEGER :: errorFlag,icl,ifl,iPatch,offs
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ENS_WriteScalar', &
                        'RFLU_ModENSIGHTUtils.F90')

  pGrid => pRegion%grid

  global%postPartNumber = global%postPartNumber + 1
  
! ******************************************************************************
! Volume data
! ******************************************************************************

  dummyString = 'part'
  WRITE(iFile) dummyString
  WRITE(iFile) global%postPartNumber

  IF ( emptyPartFlag .EQV. .FALSE. ) THEN 
  
! ==============================================================================
!   Tetrahedra
! ==============================================================================

    IF ( pGrid%nTets > 0 ) THEN 
      dummyString = 'tetra4'
      WRITE(iFile) dummyString     
      WRITE(iFile) (REAL(var(pGrid%tet2CellGlob(icl)),KIND=SPREAL), &
                   icl=1,pGrid%nTets)
    END IF ! pGrid%nTets
    
    IF ( pGrid%nTetsTot > pGrid%nTets ) THEN 
      dummyString = 'g_tetra4'
      WRITE(iFile) dummyString 
      WRITE(iFile) (REAL(var(pGrid%tet2CellGlob(icl)),KIND=SPREAL), & 
                   icl=pGrid%nTets+1,pGrid%nTetsTot)
    END IF ! pGrid%nTetsTot     

! ==============================================================================
!   Hexahedra
! ==============================================================================

    IF ( pGrid%nHexs > 0 ) THEN 
      dummyString = 'hexa8'
      WRITE(iFile) dummyString     
      WRITE(iFile) (REAL(var(pGrid%hex2CellGlob(icl)),KIND=SPREAL), &
                   icl=1,pGrid%nHexs)
    END IF ! pGrid%nHexs
    
    IF ( pGrid%nHexsTot > pGrid%nHexs ) THEN 
      dummyString = 'g_hexa8'
      WRITE(iFile) dummyString     
      WRITE(iFile) (REAL(var(pGrid%hex2CellGlob(icl)),KIND=SPREAL), &
                   icl=pGrid%nHexs+1,pGrid%nHexsTot)
    END IF ! pGrid%nHexsTot    

! ==============================================================================
!   Prisms
! ==============================================================================

    IF ( pGrid%nPris > 0 ) THEN 
      dummyString = 'penta6'
      WRITE(iFile) dummyString     
      WRITE(iFile) (REAL(var(pGrid%pri2CellGlob(icl)),KIND=SPREAL), &
                   icl=1,pGrid%nPris)
    END IF ! pGrid%nPris

    IF ( pGrid%nPrisTot > pGrid%nPris ) THEN 
      dummyString = 'g_penta6'
      WRITE(iFile) dummyString     
      WRITE(iFile) (REAL(var(pGrid%pri2CellGlob(icl)),KIND=SPREAL), &
                   icl=pGrid%nPris+1,pGrid%nPrisTot)
    END IF ! pGrid%nPrisTot  

! ==============================================================================
!   Pyramids
! ==============================================================================

    IF ( pGrid%nPyrs > 0 ) THEN 
      dummyString = 'pyramid5'
      WRITE(iFile) dummyString     
      WRITE(iFile) (REAL(var(pGrid%pyr2CellGlob(icl)),KIND=SPREAL), & 
                   icl=1,pGrid%nPyrs)
    END IF ! pGrid%nPyrs
    
    IF ( pGrid%nPyrsTot > pGrid%nPyrs ) THEN 
      dummyString = 'g_pyramid5'
      WRITE(iFile) dummyString     
      WRITE(iFile) (REAL(var(pGrid%pyr2CellGlob(icl)),KIND=SPREAL), &
                   icl=pGrid%nPyrs+1,pGrid%nPyrsTot)
    END IF ! pGrid%nPyrsTot    
  END IF ! emptyPartFlag

! ******************************************************************************
! Patch data
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
    
! ==============================================================================
!   Triangles
! ==============================================================================

! ------------------------------------------------------------------------------
!   Actual triangles 
! ------------------------------------------------------------------------------
   
    IF ( pPatch%nBTris > 0 ) THEN 
      global%postPartNumber = global%postPartNumber + 1    
    
      dummyString = 'part'
      WRITE(iFile) dummyString
      WRITE(iFile) global%postPartNumber

      IF ( emptyPartFlag .EQV. .FALSE. ) THEN 
        dummyString = 'tria3'
        WRITE(iFile) dummyString 
        WRITE(iFile) (REAL(var(pPatch%bf2c(ifl)),KIND=SPREAL), & 
                     ifl=1,pPatch%nBTris) 
      END IF ! emptyPartFlag        
    END IF ! pPatch%nBTris

! ------------------------------------------------------------------------------
!   Virtual triangles 
! ------------------------------------------------------------------------------
    
    IF ( pPatch%nBTrisTot > pPatch%nBTris ) THEN 
      global%postPartNumber = global%postPartNumber + 1    
    
      dummyString = 'part'
      WRITE(iFile) dummyString
      WRITE(iFile) global%postPartNumber

      offs = pPatch%nBTris + pPatch%nBQuads

      IF ( emptyPartFlag .EQV. .FALSE. ) THEN 
        dummyString = 'g_tria3'
        WRITE(iFile) dummyString 
        WRITE(iFile) (REAL(var(pPatch%bf2c(ifl+offs)),KIND=SPREAL), &
                     ifl=1,pPatch%nBTrisTot-pPatch%nBTris) 
      END IF ! emptyPartFlag        
    END IF ! pPatch%nBTrisTot    
    
! ==============================================================================
!   Quadrilaterals
! ==============================================================================
  
! ------------------------------------------------------------------------------
!   Actual quadrilaterals
! ------------------------------------------------------------------------------
    
    IF ( pPatch%nBQuads > 0 ) THEN
      global%postPartNumber = global%postPartNumber + 1    
     
      dummyString = 'part'
      WRITE(iFile) dummyString
      WRITE(iFile) global%postPartNumber

      offs = pPatch%nBTris
              
      IF ( emptyPartFlag .EQV. .FALSE. ) THEN 
        dummyString = 'quad4'
        WRITE(iFile) dummyString 
        WRITE(iFile) (REAL(var(pPatch%bf2c(ifl+offs)),KIND=SPREAL), & 
                     ifl=1,pPatch%nBQuads) 
      END IF ! emptyPartFlag 
    END IF ! pPatch%nBQuads 
    
! ------------------------------------------------------------------------------
!   Virtual quadrilaterals
! ------------------------------------------------------------------------------
    
    IF ( pPatch%nBQuadsTot > pPatch%nBQuads ) THEN
      global%postPartNumber = global%postPartNumber + 1    
     
      dummyString = 'part'
      WRITE(iFile) dummyString
      WRITE(iFile) global%postPartNumber
      
      offs = pPatch%nBTrisTot + pPatch%nBQuads
      
      IF ( emptyPartFlag .EQV. .FALSE. ) THEN 
        dummyString = 'g_quad4'
        WRITE(iFile) dummyString 
        WRITE(iFile) (REAL(var(pPatch%bf2c(ifl+offs)),KIND=SPREAL), &
                     ifl=1,pPatch%nBQuadsTot-pPatch%nBQuads) 
      END IF ! emptyPartFlag 
    END IF ! pPatch%nBQuads        
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ENS_WriteScalar








! ******************************************************************************
!
! Purpose: Write vector variable to ENSIGHT file.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   var                 Pointer to vector
!   iFile               File index
!   emptyPartFlag       Flag indicating whether part should be empty
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_ENS_WriteVector(pRegion,var,iFile,emptyPartFlag)

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  LOGICAL, INTENT(IN) :: emptyPartFlag
  INTEGER, INTENT(IN) :: iFile
  REAL(RFREAL), DIMENSION(:,:), POINTER :: var
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(80) :: dummyString
  INTEGER :: errorFlag,icl,ifl,iPatch,offs
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ENS_WriteVector', &
                        'RFLU_ModENSIGHTUtils.F90')

  pGrid => pRegion%grid

  global%postPartNumber = global%postPartNumber + 1

! ******************************************************************************
! Volume data
! ******************************************************************************

  dummyString = 'part'
  WRITE(iFile) dummyString
  WRITE(iFile) global%postPartNumber

  IF ( emptyPartFlag .EQV. .FALSE. ) THEN 
  
! ==============================================================================
!   Tetrahedra
! ==============================================================================

    IF ( pGrid%nTets > 0 ) THEN 
      dummyString = 'tetra4'
      WRITE(iFile) dummyString     
      WRITE(iFile) (REAL(var(1,pGrid%tet2CellGlob(icl)),KIND=SPREAL), &
                   icl=1,pGrid%nTets)
      WRITE(iFile) (REAL(var(2,pGrid%tet2CellGlob(icl)),KIND=SPREAL), &
                   icl=1,pGrid%nTets)
      WRITE(iFile) (REAL(var(3,pGrid%tet2CellGlob(icl)),KIND=SPREAL), &
                   icl=1,pGrid%nTets)      
    END IF ! pGrid%nTets
    
    IF ( pGrid%nTetsTot > pGrid%nTets ) THEN 
      dummyString = 'g_tetra4'
      WRITE(iFile) dummyString     
      WRITE(iFile) (REAL(var(1,pGrid%tet2CellGlob(icl)),KIND=SPREAL), & 
                   icl=pGrid%nTets+1,pGrid%nTetsTot)
      WRITE(iFile) (REAL(var(2,pGrid%tet2CellGlob(icl)),KIND=SPREAL), &
                   icl=pGrid%nTets+1,pGrid%nTetsTot)
      WRITE(iFile) (REAL(var(3,pGrid%tet2CellGlob(icl)),KIND=SPREAL), &
                   icl=pGrid%nTets+1,pGrid%nTetsTot)
    END IF ! pGrid%nTetsTot    

! ==============================================================================
!   Hexahedra
! ==============================================================================

    IF ( pGrid%nHexs > 0 ) THEN 
      dummyString = 'hexa8'
      WRITE(iFile) dummyString     
      WRITE(iFile) (REAL(var(1,pGrid%hex2CellGlob(icl)),KIND=SPREAL), & 
                   icl=1,pGrid%nHexs)
      WRITE(iFile) (REAL(var(2,pGrid%hex2CellGlob(icl)),KIND=SPREAL), &
                   icl=1,pGrid%nHexs)
      WRITE(iFile) (REAL(var(3,pGrid%hex2CellGlob(icl)),KIND=SPREAL), &
                   icl=1,pGrid%nHexs)       
    END IF ! pGrid%nHexs
    
    IF ( pGrid%nHexsTot > pGrid%nHexs ) THEN 
      dummyString = 'g_hexa8'
      WRITE(iFile) dummyString     
      WRITE(iFile) (REAL(var(1,pGrid%hex2CellGlob(icl)),KIND=SPREAL), & 
                   icl=pGrid%nHexs+1,pGrid%nHexsTot)
      WRITE(iFile) (REAL(var(2,pGrid%hex2CellGlob(icl)),KIND=SPREAL), &
                   icl=pGrid%nHexs+1,pGrid%nHexsTot)
      WRITE(iFile) (REAL(var(3,pGrid%hex2CellGlob(icl)),KIND=SPREAL), &
                   icl=pGrid%nHexs+1,pGrid%nHexsTot)
    END IF ! pGrid%nHexsTot    

! ==============================================================================
!   Prisms
! ==============================================================================

    IF ( pGrid%nPris > 0 ) THEN 
      dummyString = 'penta6'
      WRITE(iFile) dummyString     
      WRITE(iFile) (REAL(var(1,pGrid%pri2CellGlob(icl)),KIND=SPREAL), &
                   icl=1,pGrid%nPris)
      WRITE(iFile) (REAL(var(2,pGrid%pri2CellGlob(icl)),KIND=SPREAL), &
                   icl=1,pGrid%nPris)
      WRITE(iFile) (REAL(var(3,pGrid%pri2CellGlob(icl)),KIND=SPREAL), &
                   icl=1,pGrid%nPris)       
    END IF ! pGrid%nPris

    IF ( pGrid%nPrisTot > pGrid%nPris ) THEN 
      dummyString = 'g_penta6'
      WRITE(iFile) dummyString      
      WRITE(iFile) (REAL(var(1,pGrid%pri2CellGlob(icl)),KIND=SPREAL), & 
                   icl=pGrid%nPris+1,pGrid%nPrisTot)
      WRITE(iFile) (REAL(var(2,pGrid%pri2CellGlob(icl)),KIND=SPREAL), &
                   icl=pGrid%nPris+1,pGrid%nPrisTot)
      WRITE(iFile) (REAL(var(3,pGrid%pri2CellGlob(icl)),KIND=SPREAL), &
                   icl=pGrid%nPris+1,pGrid%nPrisTot)  
    END IF ! pGrid%nPrisTot

! ==============================================================================
!   Pyramids
! ==============================================================================

    IF ( pGrid%nPyrs > 0 ) THEN 
      dummyString = 'pyramid5'
      WRITE(iFile) dummyString     
      WRITE(iFile) (REAL(var(1,pGrid%pyr2CellGlob(icl)),KIND=SPREAL), &
                   icl=1,pGrid%nPyrs)
      WRITE(iFile) (REAL(var(2,pGrid%pyr2CellGlob(icl)),KIND=SPREAL), &
                   icl=1,pGrid%nPyrs)
      WRITE(iFile) (REAL(var(3,pGrid%pyr2CellGlob(icl)),KIND=SPREAL), &
                   icl=1,pGrid%nPyrs)       
    END IF ! pGrid%nPyrs
    
    IF ( pGrid%nPyrsTot > pGrid%nPyrs ) THEN 
      dummyString = 'g_pyramid5'
      WRITE(iFile) dummyString      
      WRITE(iFile) (REAL(var(1,pGrid%pyr2CellGlob(icl)),KIND=SPREAL), & 
                   icl=pGrid%nPyrs+1,pGrid%nPyrsTot)
      WRITE(iFile) (REAL(var(2,pGrid%pyr2CellGlob(icl)),KIND=SPREAL), &
                   icl=pGrid%nPyrs+1,pGrid%nPyrsTot)
      WRITE(iFile) (REAL(var(3,pGrid%pyr2CellGlob(icl)),KIND=SPREAL), &
                   icl=pGrid%nPyrs+1,pGrid%nPyrsTot)  
    END IF ! pGrid%nPyrsTot    
  END IF ! emptyPatchFlag

! ******************************************************************************
! Patch data
! ******************************************************************************

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
    
! ==============================================================================
!   Triangles
! ==============================================================================
   
! ------------------------------------------------------------------------------   
!   Actual triangles
! ------------------------------------------------------------------------------   
   
    IF ( pPatch%nBTris > 0 ) THEN 
      global%postPartNumber = global%postPartNumber + 1    
    
      dummyString = 'part'
      WRITE(iFile) dummyString
      WRITE(iFile) global%postPartNumber

      IF ( emptyPartFlag .EQV. .FALSE. ) THEN 
        dummyString = 'tria3'
        WRITE(iFile) dummyString
        WRITE(iFile) (REAL(var(1,pPatch%bf2c(ifl)),KIND=SPREAL), &
                     ifl=1,pPatch%nBTris) 
        WRITE(iFile) (REAL(var(2,pPatch%bf2c(ifl)),KIND=SPREAL), &
                     ifl=1,pPatch%nBTris)
        WRITE(iFile) (REAL(var(3,pPatch%bf2c(ifl)),KIND=SPREAL), &
                     ifl=1,pPatch%nBTris)        
      END IF ! emptyPartFlag                     
    END IF ! pPatch%nBTris

! ------------------------------------------------------------------------------   
!   Virtual triangles   
! ------------------------------------------------------------------------------   
    
    IF ( pPatch%nBTrisTot > pPatch%nBTris ) THEN 
      global%postPartNumber = global%postPartNumber + 1    
    
      dummyString = 'part'
      WRITE(iFile) dummyString
      WRITE(iFile) global%postPartNumber

      offs = pPatch%nBTris + pPatch%nBQuads

      IF ( emptyPartFlag .EQV. .FALSE. ) THEN 
        dummyString = 'g_tria3'
        WRITE(iFile) dummyString 
        WRITE(iFile) (REAL(var(1,pPatch%bf2c(ifl+offs)),KIND=SPREAL), &
                     ifl=1,pPatch%nBTrisTot-pPatch%nBTris) 
        WRITE(iFile) (REAL(var(2,pPatch%bf2c(ifl+offs)),KIND=SPREAL), &
                     ifl=1,pPatch%nBTrisTot-pPatch%nBTris)
        WRITE(iFile) (REAL(var(3,pPatch%bf2c(ifl+offs)),KIND=SPREAL), &
                     ifl=1,pPatch%nBTrisTot-pPatch%nBTris)  
      END IF ! emptyPartFlag                     
    END IF ! pPatch%nBTrisTot    
    
! ==============================================================================
!   Quadrilaterals
! ==============================================================================

! ------------------------------------------------------------------------------   
!   Actual quadrilaterals   
! ------------------------------------------------------------------------------   
  
    IF ( pPatch%nBQuads > 0 ) THEN
      global%postPartNumber = global%postPartNumber + 1    
     
      dummyString = 'part'
      WRITE(iFile) dummyString
      WRITE(iFile) global%postPartNumber
      
      offs = pPatch%nBTris      
      
      IF ( emptyPartFlag .EQV. .FALSE. ) THEN 
        dummyString = 'quad4'
        WRITE(iFile) dummyString 
        WRITE(iFile) (REAL(var(1,pPatch%bf2c(ifl+offs)),KIND=SPREAL), &
                     ifl=1,pPatch%nBQuads)
        WRITE(iFile) (REAL(var(2,pPatch%bf2c(ifl+offs)),KIND=SPREAL), &
                     ifl=1,pPatch%nBQuads)
        WRITE(iFile) (REAL(var(3,pPatch%bf2c(ifl+offs)),KIND=SPREAL), &
                     ifl=1,pPatch%nBQuads) 
      END IF ! emptyPartFlag             
    END IF ! pPatch%nBQuads 

! ------------------------------------------------------------------------------   
!   Virtual quadrilaterals   
! ------------------------------------------------------------------------------   
    
    IF ( pPatch%nBQuadsTot > pPatch%nBQuads ) THEN
      global%postPartNumber = global%postPartNumber + 1    
     
      dummyString = 'part'
      WRITE(iFile) dummyString
      WRITE(iFile) global%postPartNumber
      
      offs = pPatch%nBTrisTot + pPatch%nBQuads
      
      IF ( emptyPartFlag .EQV. .FALSE. ) THEN 
        dummyString = 'g_quad4'
        WRITE(iFile) dummyString 
        WRITE(iFile) (REAL(var(1,pPatch%bf2c(ifl+offs)),KIND=SPREAL), &
                     ifl=1,pPatch%nBQuadsTot-pPatch%nBQuads) 
        WRITE(iFile) (REAL(var(2,pPatch%bf2c(ifl+offs)),KIND=SPREAL), &
                     ifl=1,pPatch%nBQuadsTot-pPatch%nBQuads)  
        WRITE(iFile) (REAL(var(3,pPatch%bf2c(ifl+offs)),KIND=SPREAL), &
                     ifl=1,pPatch%nBQuadsTot-pPatch%nBQuads)                                             
      END IF ! emptyPartFlag             
    END IF ! pPatch%nBQuads               
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ENS_WriteVector







END MODULE RFLU_ModENSIGHTUtils

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModENSIGHTUtils.F90,v $
! Revision 1.4  2008/12/06 08:45:06  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:16  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.1  2005/10/05 20:23:34  haselbac
! Initial revision
!
! ******************************************************************************









