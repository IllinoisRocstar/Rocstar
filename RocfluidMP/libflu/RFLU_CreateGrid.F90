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
! Purpose: Create grid.
!
! Description: None.
!
! Input:
!   pRegion             Region pointer
!
! Output: None.
!
! Notes: 
!   1. This routine creates the basic grid arrays only, i.e., those read in
!      by RFLU_ReadGrid{ASCII/Binary}.F90. Other grid-associated arrays
!      are created in other routines. One exception is the local boundary 
!      face connectivity arrays required for GENX runs. These arrays are 
!      created by the call to RFLU_CreateBFaceLocLists. 
!   2. The dimensions of the various arrays MUST have been previously read 
!      by a call to RFLU_ReadDimensions.F90.
!
! ******************************************************************************
!
! $Id: RFLU_CreateGrid.F90,v 1.25 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CreateGrid(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModBndPatch, ONLY: t_patch
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModMPI

  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,iPatch
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_CreateGrid.F90,v $ $Revision: 1.25 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_CreateGrid',&
  'RFLU_CreateGrid.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating grid...'
  END IF ! global%verbLevel

! ******************************************************************************
! Set grid pointer
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Coordinates and vertex flags
! ******************************************************************************

  ALLOCATE(pGrid%xyz(3,pGrid%nVertMax),STAT=errorFlag)
  global%error = errorFlag         
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%xyz')
  END IF ! global%error
  
! ******************************************************************************
! Connectivity
! ******************************************************************************

! ==============================================================================
! Tetrahedra
! ==============================================================================

  IF ( pGrid%nTetsMax > 0 ) THEN 
    ALLOCATE(pGrid%tet2v(4,pGrid%nTetsMax),STAT=errorFlag) 
    global%error = errorFlag         
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%tet2v')
    END IF ! global%error  
  ELSE 
    NULLIFY(pGrid%tet2v) 
  END IF ! pGrid%nTetsMax

! ==============================================================================
! Hexahedra
! ==============================================================================

  IF ( pGrid%nHexsMax > 0 ) THEN 
    ALLOCATE(pGrid%hex2v(8,pGrid%nHexsMax),STAT=errorFlag) 
    global%error = errorFlag 
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%hex2v')
    END IF ! global%error  
  ELSE 
    NULLIFY(pGrid%hex2v) 
  END IF ! pGrid%nHexsTot

! ==============================================================================
! Prisms
! ==============================================================================

  IF ( pGrid%nPrisMax > 0 ) THEN 
    ALLOCATE(pGrid%pri2v(6,pGrid%nPrisMax),STAT=errorFlag)
    global%error = errorFlag          
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%pri2v')
    END IF ! global%error 
  ELSE 
    NULLIFY(pGrid%pri2v)  
  END IF ! pGrid%nPrisTot               
        
! ==============================================================================
! Pyramids 
! ==============================================================================

  IF ( pGrid%nPyrsMax > 0 ) THEN 
    ALLOCATE(pGrid%pyr2v(5,pGrid%nPyrsMax),STAT=errorFlag)
    global%error = errorFlag          
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%grid%pyr2v')
    END IF ! global%error
  ELSE 
    NULLIFY(pGrid%pyr2v)  
  END IF ! pGrid%nPyrsTot

! ==============================================================================
! Initialize some dimensions
! ==============================================================================

  pGrid%nCellsSpecial = 0
  pGrid%nFacesSpecial = 0

  pGrid%nBFaces    = 0 
  pGrid%nBFacesTot = 0 
  
  pGrid%nCellsConstr = 0
  pGrid%nFacesConstr = 0  

! ******************************************************************************
! Patches
! ******************************************************************************

  IF ( pGrid%nPatches > 0 ) THEN 
    ALLOCATE(pRegion%patches(pGrid%nPatches),STAT=errorFlag)  
    global%error = errorFlag         
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%patches')
    END IF ! global%error
  ELSE 
    NULLIFY(pRegion%patches)
  END IF ! pGrid%nPatches

! ==============================================================================
! Loop over patches
! ==============================================================================

  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)

! ------------------------------------------------------------------------------
!   Copy data read in from dimension file into patch type
! ------------------------------------------------------------------------------

    pPatch%iPatchLocal  = iPatch
    pPatch%iPatchGlobal = pGrid%patchDimens(PATCH_DIMENS_IPGLOBAL,iPatch)
    
    pPatch%nBTris      = pGrid%patchDimens(PATCH_DIMENS_NBTRIS     ,iPatch)
    pPatch%nBTrisTot   = pGrid%patchDimens(PATCH_DIMENS_NBTRISTOT  ,iPatch)
    pPatch%nBTrisMax   = pGrid%patchDimens(PATCH_DIMENS_NBTRISMAX  ,iPatch) 
       
    pPatch%nBQuads     = pGrid%patchDimens(PATCH_DIMENS_NBQUADS    ,iPatch)    
    pPatch%nBQuadsTot  = pGrid%patchDimens(PATCH_DIMENS_NBQUADSTOT ,iPatch)
    pPatch%nBQuadsMax  = pGrid%patchDimens(PATCH_DIMENS_NBQUADSMAX ,iPatch)    

    pPatch%nBCellsVirt = pGrid%patchDimens(PATCH_DIMENS_NBCELLSVIRT,iPatch) 

    pPatch%nBFaces    = pPatch%nBTris    + pPatch%nBQuads 
    pPatch%nBFacesTot = pPatch%nBTrisTot + pPatch%nBQuadsTot
    pPatch%nBFacesMax = pPatch%nBTrisMax + pPatch%nBQuadsMax    

    pGrid%nBFaces     = pGrid%nBFaces    + pPatch%nBFaces
    pGrid%nBFacesTot  = pGrid%nBFacesTot + pPatch%nBFacesTot 

    pPatch%nBVert     = 0     
    pPatch%nBVertTot  = 0 
                
! ------------------------------------------------------------------------------
!   Set variables
! ------------------------------------------------------------------------------

    pPatch%plotFlag      = .TRUE. 
    pPatch%renumFlag     = .FALSE.
    pPatch%flatFlag      = .FALSE.
    pPatch%transformFlag = .FALSE.

    pPatch%movePatchDir = 0

! ------------------------------------------------------------------------------
!   Allocate remaining arrays (see also note above)  
! ------------------------------------------------------------------------------

! - Face arrays ----------------------------------------------------------------

    IF ( pPatch%nBTrisMax > 0 ) THEN 
      ALLOCATE(pPatch%bTri2v(3,pPatch%nBTrisMax),STAT=errorFlag)
      global%error = errorFlag             
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%patches%bTri2v')
      END IF ! global%error
    ELSE 
      NULLIFY(pPatch%bTri2v) 
    END IF ! pPatch%nBTrisTot

    IF ( pPatch%nBQuadsMax > 0 ) THEN 
      ALLOCATE(pPatch%bQuad2v(4,pPatch%nBQuadsMax),STAT=errorFlag)
      global%error = errorFlag             
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'region%patch%bQuad2v')
      END IF ! global%error    
    ELSE 
      NULLIFY(pPatch%bQuad2v)
    END IF ! pPatch%nBQuadsTot

! - Virtual cell array ---------------------------------------------------------
    
    IF ( pPatch%nBCellsVirt > 0 ) THEN                                                                          
      ALLOCATE(pPatch%bvc(pPatch%nBCellsVirt),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%bvc')
      END IF ! global%error
    ELSE 
      NULLIFY(pPatch%bvc)        
    END IF ! pPatch%nBCellsVirt
    
! ------------------------------------------------------------------------------
!   Initialize transformation matrix  
! ------------------------------------------------------------------------------
    
    pPatch%tm(XCOORD,XCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    pPatch%tm(XCOORD,YCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    pPatch%tm(XCOORD,ZCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    pPatch%tm(XCOORD,XYZMAG) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    pPatch%tm(YCOORD,XCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    pPatch%tm(YCOORD,YCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    pPatch%tm(YCOORD,ZCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    pPatch%tm(YCOORD,XYZMAG) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    pPatch%tm(ZCOORD,XCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    pPatch%tm(ZCOORD,YCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    pPatch%tm(ZCOORD,ZCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    pPatch%tm(ZCOORD,XYZMAG) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    pPatch%tm(XYZMAG,XCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    pPatch%tm(XYZMAG,YCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    pPatch%tm(XYZMAG,ZCOORD) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    pPatch%tm(XYZMAG,XYZMAG) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************  

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating grid done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CreateGrid

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CreateGrid.F90,v $
! Revision 1.25  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.24  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.23  2006/08/18 13:59:01  haselbac
! Added init of transform flag and tm
!
! Revision 1.22  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.21  2006/04/07 14:42:00  haselbac
! Added setting of pPatch%iPatchLocal
!
! Revision 1.20  2006/03/25 21:41:48  haselbac
! Changes made bcos of sype boundaries
!
! Revision 1.19  2005/10/27 18:55:36  haselbac
! Added init of nCellsConstr and nFacesConstr
!
! Revision 1.18  2005/08/09 00:53:33  haselbac
! Added init of plotFlag for patches
!
! Revision 1.17  2005/06/13 22:42:40  haselbac
! Removed setting of cnstrType
!
! Revision 1.16  2005/06/10 19:44:12  haselbac
! Bug fix: Wrong module
!
! Revision 1.15  2005/06/10 18:04:38  haselbac
! Bug fix for backward-compatibility in GENx with cnstr_type
!
! Revision 1.14  2005/01/14 21:09:18  haselbac
! Removed init of nBorders, otherwise cannot read comm maps in solver
!
! Revision 1.13  2004/12/29 21:01:14  haselbac
! Added setting of pGrid%nBFaces and pGrid%nBFacesTot
!
! Revision 1.12  2004/12/04 03:21:37  haselbac
! Added initialization of nBorders
!
! Revision 1.11  2004/11/09 00:27:04  haselbac
! Bug fix: Used nTetsTot instead of nTetsMax
!
! Revision 1.10  2004/11/03 16:58:39  haselbac
! Removed allocation of vertex and cell flags
!
! Revision 1.9  2004/10/19 19:24:13  haselbac
! Introduced Max dims, removed bf arrays and derived dims
!
! Revision 1.8  2004/10/08 21:08:00  fnajjar
! ACH: Bug fix: Added initialization of nFacesSpecial
!
! Revision 1.7  2004/07/06 15:14:13  haselbac
! Bug fix: patchDimens now under grid type, cosmetics
!
! Revision 1.6  2003/11/03 03:48:53  haselbac
! Cosmetic changes only
!
! Revision 1.5  2003/06/04 22:01:16  haselbac
! Added NULLIFY for patches
!
! Revision 1.4  2003/04/01 19:37:14  haselbac
! Added initialization of nCellsSpecial
!
! Revision 1.3  2003/03/25 19:11:52  haselbac
! Added setting of patch renumbering flag
!
! Revision 1.2  2003/03/15 16:52:04  haselbac
! Added vertex arrays, allocate based on *Tot
!
! Revision 1.1  2003/01/28 15:53:31  haselbac
! Initial revision
!
! ******************************************************************************







