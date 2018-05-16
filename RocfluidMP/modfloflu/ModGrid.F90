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
! Purpose: define the derived data type related to grid.
!
! Description: none
!
! Notes: none
!
! ******************************************************************************
!
! $Id: ModGrid.F90,v 1.70 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2001-2006 by the University of Illinois
!
! ******************************************************************************

MODULE ModGrid

  USE ModDataTypes
  USE ModParameters
  
#ifdef RFLU  
  USE ModColoring, ONLY: t_soc
  USE ModStencil, ONLY: t_stencil,t_stencilInfo
  USE ModBorder, ONLY: t_border
#endif  
  
  IMPLICIT NONE

! ******************************************************************************
! Type definition
! ******************************************************************************

  TYPE t_grid
    INTEGER :: indSvel
#ifdef RFLO
! ==============================================================================
!   Basic grid quantities
! ==============================================================================

    INTEGER :: ipc, jpc, kpc
    LOGICAL :: boundMoved(6), allExternal(6), edgeMoved(12)
    LOGICAL :: boundFlat(6), edgeStraight(12)
    REAL(RFREAL) :: skewness, minVol
    INTEGER, POINTER :: ijkDgen(:)                           ! degeneration flag
    REAL(RFREAL), POINTER :: c2fCoI(:,:) ,c2fCoJ(:,:) ,c2fCoK(:,:)
    REAL(RFREAL), POINTER :: c2eCoI(:,:) ,c2eCoJ(:,:) ,c2eCoK(:,:)
    REAL(RFREAL), POINTER :: si(:,:), sj(:,:), sk(:,:)
    REAL(RFREAL), POINTER :: siVel(:), sjVel(:), skVel(:)
    REAL(RFREAL), POINTER :: arcLen12(:,:), arcLen34(:,:), arcLen56(:,:)
    REAL(RFREAL), POINTER :: xyz(:,:), xyzOld(:,:)
    REAL(RFREAL), POINTER :: vol(:), cofg(:,:)
    REAL(RFREAL), POINTER :: cfcI(:,:), cfcJ(:,:), cfcK(:,:) ! face-cntrs (FOMS)

! ==============================================================================
!   Moving grid
! ==============================================================================

    INTEGER :: remesh                            ! block remeshing
    INTEGER,      POINTER :: nCorns(:)           ! number of motion points
    INTEGER,      POINTER :: ijkCorn(:,:)        ! ijk value of corn(iCorn,iReg)
    INTEGER,      POINTER :: nghbor(:,:,:)       ! for block corners motion
    INTEGER,      POINTER :: nShared(:)          ! # shared corners (iCorn,iReg)
    INTEGER,      POINTER :: cshared(:,:,:)      ! crnrs shared wother ptch/blck
    REAL(RFREAL), POINTER :: regCorn(:,:,:)      ! current corner motion
    REAL(RFREAL), POINTER :: regCornOld(:,:,:)   ! old corner motion
    REAL(RFREAL), POINTER :: regCornOrig(:,:,:)  ! previous tstep corner motion
    REAL(RFREAL), POINTER :: regCornBuff(:,:,:)  ! buffer store latest cornCoord
    REAL(RFREAL), POINTER :: regCornOrth(:,:,:)  ! orthogonally moved corners

    REAL(RFREAL), POINTER :: xyzOrth(:,:), xyzTemp(:,:)      ! work-space (FOMS)

    REAL(RFREAL), POINTER :: stu(:,:), stuOld(:,:)           ! param grid (EPDE)
    REAL(RFREAL), POINTER :: stui(:,:), stuj(:,:), stuk(:,:) ! 1st deriv. stu
    REAL(RFREAL), POINTER :: stuii(:,:) ,stujj(:,:)          ! 2nd deriv. stu
    REAL(RFREAL), POINTER :: stukk(:,:) ,stuij(:,:)          ! 2nd deriv. stu
    REAL(RFREAL), POINTER :: stuik(:,:) ,stujk(:,:)          ! 2nd deriv. stu
    REAL(RFREAL), POINTER :: pmat(:,:,:)                     ! control function

    REAL(RFREAL), POINTER :: aijk(:),aimjk(:),aipjk(:)       ! matrix coeffs.
    REAL(RFREAL), POINTER :: aijmk(:),aijpk(:),aijkm(:)      ! in Elliptic PDE
    REAL(RFREAL), POINTER :: aijkp(:),aimjmk(:),aipjmk(:)    ! grid motion
    REAL(RFREAL), POINTER :: aimjpk(:),aipjpk(:),aimjkm(:)
    REAL(RFREAL), POINTER :: aipjkm(:),aimjkp(:),aipjkp(:)
    REAL(RFREAL), POINTER :: aijmkm(:),aijpkm(:),aijmkp(:)
    REAL(RFREAL), POINTER :: aijpkp(:)

! ==============================================================================
!   Conversion to unstructured grid
! ==============================================================================

    INTEGER, POINTER :: tofluLoc2g(:,:,:) 
#endif

#ifdef RFLU
! ==============================================================================
!   Basic grid quantities
! ==============================================================================

! ------------------------------------------------------------------------------
!   Dimensions
! ------------------------------------------------------------------------------

    INTEGER :: nCells,nCellsTot,nCellsMax
    INTEGER :: nEdges,nEdgesTot,nEdgesEst
    INTEGER :: nFaces,nFacesTot,nFacesEst,nFacesCut,nFacesAV,nFacesVV
    INTEGER :: nHexs,nHexsTot,nHexsMax
    INTEGER :: nPris,nPrisTot,nPrisMax    
    INTEGER :: nPyrs,nPyrsTot,nPyrsMax    
    INTEGER :: nTets,nTetsTot,nTetsMax
    INTEGER :: nVert,nVertTot,nVertMax,nVertInt

    INTEGER :: nBCells,nBCellsTot
    INTEGER :: nBFaces,nBFacesTot
    INTEGER :: nBVert,nBVertEst
    INTEGER :: nPatches,nPatchesMax

! ------------------------------------------------------------------------------
!   Cell-to-vertex connectivity
! ------------------------------------------------------------------------------

    INTEGER, DIMENSION(:,:), POINTER :: hex2v,pri2v,pyr2v,tet2v

! ------------------------------------------------------------------------------
!   Local-to-global and global-to-local cell indexing
! ------------------------------------------------------------------------------
    
    INTEGER, DIMENSION(:), POINTER :: hex2CellGlob,pri2CellGlob, &
                                      pyr2CellGlob,tet2CellGlob
    INTEGER, DIMENSION(:,:), POINTER :: cellGlob2Loc

! ------------------------------------------------------------------------------
!   Data related to partitioning
! ------------------------------------------------------------------------------
     
    INTEGER, DIMENSION(:), POINTER :: pc2sc,pbf2sbfCSR,pbf2sbfCSRInfo, &
                                      pv2sv,r2pcCSR,r2pbcCSRInfo,r2pcCSRInfo, &
                                      sc2r
    INTEGER, DIMENSION(:), POINTER :: vertKind
    INTEGER, DIMENSION(:), POINTER :: avfCSR,avfCSRInfo
    INTEGER, DIMENSION(:,:), POINTER :: avf,r2pbcCSR,sbc2pc,sc2pc,sv2pv 
    INTEGER, DIMENSION(:,:,:), POINTER :: patchCounter

! ------------------------------------------------------------------------------
!   Edge lists
! ------------------------------------------------------------------------------
                                          
    INTEGER, DIMENSION(:), POINTER :: e2c,e2cDegr,e2cStrt,e2rDegr
    INTEGER, DIMENSION(:,:), POINTER :: e2v,e2vTemp                                         

! ------------------------------------------------------------------------------
!   Face lists
! ------------------------------------------------------------------------------
                                         
    INTEGER, DIMENSION(:,:), POINTER :: f2c,f2cTemp
    INTEGER, DIMENSION(:,:), POINTER :: f2v,f2vTemp      

! ------------------------------------------------------------------------------
!   Vertex-to-cell list
! ------------------------------------------------------------------------------

    INTEGER, DIMENSION(:), POINTER :: v2c
    INTEGER, DIMENSION(:,:), POINTER :: v2cInfo

! ------------------------------------------------------------------------------
!   Cell-to-face lists
! ------------------------------------------------------------------------------
                                        
    INTEGER, DIMENSION(:), POINTER :: hex2fCntr,pri2fCntr,pyr2fCntr,tet2fCntr
    INTEGER, DIMENSION(:,:,:), POINTER :: hex2f,pri2f,pyr2f,tet2f                                              

! ------------------------------------------------------------------------------
!   Actual-virtual-face-to-border and patch lists
! ------------------------------------------------------------------------------
                                        
    INTEGER, DIMENSION(:), POINTER :: avf2p                                        
    INTEGER, DIMENSION(:,:), POINTER :: avf2b

! ==============================================================================
!   Other
! ==============================================================================
    
    INTEGER :: indGs  
    INTEGER, DIMENSION(:), POINTER :: bcm,bvm
                
! ==============================================================================
!   Patch dimensions
! ==============================================================================
    
    INTEGER, DIMENSION(PATCH_DIMENS_BEG:PATCH_DIMENS_END, & 
                       PATCH_DIMENS_NPATCHMAX) :: patchDimens    

! ==============================================================================
!   Borders
! ==============================================================================
        
    INTEGER :: nBorders    
    INTEGER, DIMENSION(BORDER_INFO_BEG:BORDER_INFO_END, & 
                       BORDER_INFO_MAX) :: borderInfo
    INTEGER, DIMENSION(:), POINTER :: borderCntr                        
    TYPE(t_border), DIMENSION(:), POINTER :: borders    
        
! ==============================================================================
!   Coloring
! ==============================================================================
        
    INTEGER :: nSoc,nSocMax   
    INTEGER, DIMENSION(:), POINTER:: col
    TYPE(t_soc), DIMENSION(:), POINTER :: soc           
        
! ==============================================================================
!   Geometry
! ==============================================================================

    REAL(RFREAL), DIMENSION(:), POINTER :: vol
    REAL(RFREAL), DIMENSION(:,:), POINTER :: cofg,cofgApp,cofgDist,fc,fn,xyz, &
                                             xyzOld
        
! ==============================================================================
!   Stencils 
! ==============================================================================

    INTEGER :: nCellsConstr,nFacesConstr
    INTEGER, DIMENSION(:), POINTER :: c2csKeyCSR,c2csCSR,icgConstr,ifgConstr
    TYPE(t_stencil), DIMENSION(:), POINTER :: c2cs,f2cs,f2cs1D,v2cs
    TYPE(t_stencil), DIMENSION(:,:), POINTER :: c2cs1D            
    TYPE(t_stencilInfo) :: c2csInfo,f2csInfo,f2cs1DInfo,v2csInfo

! ==============================================================================
!   Special cells and faces
! ==============================================================================
    
    INTEGER :: nCellsSpecial,nFacesSpecial
    INTEGER, DIMENSION(NCELLS_SPECIAL_MAX) :: cellsSpecial 
    INTEGER, DIMENSION(2,NFACES_SPECIAL_MAX) :: facesSpecial

! ==============================================================================
!   Optimal LES 
! ==============================================================================

    INTEGER, DIMENSION(:), POINTER :: f2fpOLES,fp2fOLES
    INTEGER, DIMENSION(:,:), POINTER :: fsOLES
    REAL(RFREAL) :: deltaOLES
    REAL(RFREAL), DIMENSION(:), POINTER :: rhoOLES,rhsOLES
    REAL(RFREAL), DIMENSION(:,:), POINTER :: lhsOLES,lhsInvOLES
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: intLimOLES,int1OLES, &
                                               int20OLES,int21OLES,int22OLES, &
                                               int31OLES,int32OLES, & 
                                               int40OLES,int41OLES,int42OLES, &
                                               int50OLES,int51OLES,int52OLES
    REAL(RFREAL), DIMENSION(:,:,:,:), POINTER :: wtLinOLES
    REAL(RFREAL), DIMENSION(:,:,:,:,:,:), POINTER :: wtQuadOLES 
    
! ==============================================================================
!   Poisson matrix
! ==============================================================================    
    
    INTEGER, DIMENSION(:), POINTER :: poissonNNZ
    INTEGER, DIMENSION(:), POINTER :: poissonCol,poissonFirst,poissonLast
    REAL(RFREAL), DIMENSION(:), POINTER :: poissonA
    
! ==============================================================================
!   Grid Motion 
! ==============================================================================

    INTEGER, DIMENSION(:), POINTER :: degr
    REAL(RFREAL) :: fsScaleFactor
    REAL(RFREAL), DIMENSION(:), POINTER :: gmEdgeWght,gmVertWght,gs,volMin
    REAL(RFREAL), DIMENSION(:,:), POINTER :: disp,rhs

! ==============================================================================
!   Charm 
! ==============================================================================

    INTEGER :: charmCellType,charmVertType

! ==============================================================================
!   GENx 
! ==============================================================================

    INTEGER :: pconnSizeTot,pconnSizeGhost
    INTEGER, DIMENSION(:), POINTER :: pconn

! ==============================================================================
!   Periodic boundaries - Hack
! ==============================================================================

    INTEGER, DIMENSION(:), POINTER :: cellMapPeriod
#endif

  END TYPE t_grid

END MODULE ModGrid

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModGrid.F90,v $
! Revision 1.70  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.69  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.68  2006/08/18 13:59:53  haselbac
! Added avf2p
!
! Revision 1.67  2006/03/25 21:46:35  haselbac
! Added nFacesVV for sype patches
!
! Revision 1.66  2006/03/18 11:05:39  wasistho
! added minvol
!
! Revision 1.65  2006/03/15 06:41:37  wasistho
! added region skewness
!
! Revision 1.64  2006/03/14 04:33:08  wasistho
! added edgeStraight
!
! Revision 1.63  2006/03/12 10:32:30  wasistho
! added boundFlat
!
! Revision 1.62  2006/03/09 14:06:05  haselbac
! Added f2cs1D
!
! Revision 1.61  2006/01/06 22:08:12  haselbac
! Added c2cs1D
!
! Revision 1.60  2005/12/03 09:34:59  wasistho
! added arrays for movegrid EPDE
!
! Revision 1.59  2005/11/29 04:49:28  wasistho
! added new variables for elliptic PDE gridmotion
!
! Revision 1.58  2005/11/11 07:17:01  wasistho
! added ijkDgen
!
! Revision 1.57  2005/10/27 18:57:41  haselbac
! Added variables for constraints
!
! Revision 1.56  2005/10/27 05:10:50  wasistho
! added xyzOrth
!
! Revision 1.55  2005/10/20 06:48:48  wasistho
! added cfcI,J,K
!
! Revision 1.54  2005/09/19 18:39:03  haselbac
! Added borderCntr
!
! Revision 1.53  2005/08/25 23:09:15  wasistho
! added moving grid variables for block orthogonality
!
! Revision 1.52  2005/08/19 18:24:56  haselbac
! Fixed bug in definition of col
!
! Revision 1.51  2005/08/19 02:32:35  haselbac
! Renamed cColors to col
!
! Revision 1.50  2005/08/17 20:17:37  hdewey2
! Added coloring variables
!
! Revision 1.49  2005/06/10 19:32:05  wasistho
! move memory allocation of movegridFrame variables to modflo/RFLO_ModMoveGridFrame
!
! Revision 1.48  2005/06/05 23:02:47  wasistho
! distinguish external boundary to be fully and partly external
!
! Revision 1.47  2005/05/29 20:44:21  wasistho
! added regCornOrig for rflo
!
! Revision 1.46  2005/05/28 08:05:43  wasistho
! added data structure for moveGridFrame
!
! Revision 1.45  2005/05/27 01:52:25  wasistho
! added rflo remeshing
!
! Revision 1.44  2005/05/18 22:06:47  fnajjar
! ACH: Added nFacesAV for parallel PLAG
!
! Revision 1.43  2005/04/29 23:00:32  haselbac
! Added avf2b array
!
! Revision 1.42  2005/01/20 14:48:28  haselbac
! Added arrays for partitioning
!
! Revision 1.41  2005/01/13 21:37:01  haselbac
! Bug fix: poissonA now declared as REAL
!
! Revision 1.40  2004/12/29 21:03:04  haselbac
! Added arrays for parallelization
!
! Revision 1.39  2004/12/19 15:44:57  haselbac
! Added arrays for incompressible solver
!
! Revision 1.38  2004/12/04 03:24:31  haselbac
! Added further variables for partitioning
!
! Revision 1.37  2004/11/09 00:28:17  haselbac
! Added data related to partitioning
!
! Revision 1.36  2004/11/03 16:59:57  haselbac
! Removed vertex and cell flag arrays, cosmetics
!
! Revision 1.35  2004/10/19 19:28:53  haselbac
! Clean-up
!
! Revision 1.34  2004/09/27 01:36:42  haselbac
! Added counter and array for special faces
!
! Revision 1.33  2004/08/17 00:55:02  wasistho
! prepared for utilities/rocflo/toflu
!
! Revision 1.32  2004/08/03 00:33:51  wasistho
! added c2eConI,J,K in RFLO
!
! Revision 1.31  2004/08/02 19:33:19  wasistho
! changed grid%avgCo to grid%c2fCo
!
! Revision 1.30  2004/07/30 17:23:40  wasistho
! provide cell2face averaging coefficients
!
! Revision 1.29  2004/07/06 15:14:14  haselbac
! Bug fix: Removed patchDimens, now under grid type
!
! Revision 1.28  2004/06/16 20:00:53  haselbac
! Added nBFaces, cosmetics
!
! Revision 1.27  2004/04/14 02:08:12  haselbac
! Added fsScaleFactor
!
! Revision 1.26  2003/12/05 16:54:11  haselbac
! Added c2csKeyCSR and c2csCSR
!
! Revision 1.25  2003/12/04 03:28:24  haselbac
! Added new stencil arrays, deleted old ones
!
! Revision 1.24  2003/08/11 21:50:00  jblazek
! Splitted ModInterfaces into 4 sections.
!
! Revision 1.23  2003/07/22 02:01:21  haselbac
! Added vertex-interpolation weights and stencil
!
! Revision 1.22  2003/04/28 22:42:02  haselbac
! Added vars for merged bv list for RFLU
!
! Revision 1.21  2003/04/07 14:23:24  haselbac
! Added arrays for cell-to-face lists
!
! Revision 1.20  2003/04/01 16:38:26  haselbac
! Added cellsSpecial and nCellsSpecial
!
! Revision 1.19  2003/03/31 16:13:45  haselbac
! Added variables for grid motion based on disp
!
! Revision 1.18  2003/03/15 17:47:51  haselbac
! Added some new variables
!
! Revision 1.17  2003/03/14 22:05:11  jblazek
! Improved mesh motion algorithm - node movement exchaged between blocks.
!
! Revision 1.16  2003/01/28 16:41:46  haselbac
! Added approx cofg, removed boundary faces numbers
!
! Revision 1.15  2002/11/08 21:24:48  haselbac
! Some clean-up
!
! Revision 1.14  2002/10/27 19:02:12  haselbac
! Added several variables for edge list and grid motion
!
! Revision 1.13  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.12  2002/09/09 14:54:49  haselbac
! Added nBFaces, face gradient weights, and OLES arrays
!
! Revision 1.11  2002/08/29 22:36:57  jblazek
! Changed name of index pointer for grid speeds.
!
! Revision 1.10  2002/07/25 15:11:33  haselbac
! Added OLES stuff
!
! Revision 1.9  2002/06/27 15:55:17  haselbac
! Added arrays for parallelization
!
! Revision 1.8  2002/06/14 20:14:53  haselbac
! Added grid speed stuff and vertex and cell flags
!
! Revision 1.7  2002/04/11 18:51:50  haselbac
! Added arrays for v2c data structure
!
! Revision 1.6  2002/03/14 19:05:43  haselbac
! Added entries for face centroid and face normal
!
! Revision 1.5  2002/03/01 16:48:35  haselbac
! Added some arrays for data structure and its generation
!
! Revision 1.4  2002/02/27 18:38:20  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
! Revision 1.3  2001/12/21 23:00:21  haselbac
! Added ROCFLU boundary variables
!
! Revision 1.2  2001/12/04 17:15:45  haselbac
! Added ROCFLU variables
!
! Revision 1.1  2001/12/04 00:07:01  jblazek
! Modules BndPatch, Global and Grid moved to modfloflu directory.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
! ******************************************************************************






