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
! Purpose: Move grid by smoothing coordinates.
!
! Description: Simple pseudo-Laplacian smoothing of grid given boundary 
!   distortion. 
!
! Input: 
!   regions     Region data
!   context     Context in which grid smoothing is to be done 
!
! Output: None.
!
! Notes: 
!   1. The input variable context allows the routine to distinguish between 
!      two kinds of grid smoothing operations: The first is an operation on 
!      a grid with the goal of smoothing the interior grid with no boundary
!      motion. This could be done to smooth an original grid. The second is 
!      the smoothing of the grid with boundary motion. Since the operations 
!      to be carried out are slighly different, it is convenient to use the
!      context variable to distinguish between them. The context variable can
!      take the two (integer parameter) values MOVEGRID_CONTEXT_ONLYSMOOTH 
!      and MOVEGRID_CONTEXT_MOVESMOOTH.
!   2. The two parameters influencing the grid motion, i.e., the number of
!      iterations and the smoothing coefficient, may need some tuning as more
!      complicated geometries are attempted...
!   3. The integrity of patches with no imposed motion is ensured by one of 
!      two methods. Both may be cast in terms of a boundary condition on the 
!      residual of the Laplacian operator. The first is a von Neumann bc on 
!      the residual and can be used for flat boundaries which are not aligned
!      the normal to which is not aligned with a coordinate direction. The 
!      second is a Dirichlet boundary condition which can be used for flat 
!      boundaries the normal to which is aligned with a coordinate direction.
!      Of course, the first method subsumes the second, but the latter may be
!      preferable because the first can introduce errors on the order of the 
!      machine precision into the constant coordinate component. Which of the 
!      two methods is more suitable is determined by a simple algorithm in
!      the routine RFLU_SetMoveGridOptions.F90. 
!
! ******************************************************************************
!
! $Id: RFLU_MoveGridXyz.F90,v 1.9 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_MoveGridXyz(regions,context)

  USE ModDataTypes
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: context
  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,ic,ie,iPatch,iReg,it,iv,ivg,nIter,patchCheckCounter, & 
             v1,v2,v3,v4
  REAL(RFREAL) :: nx,ny,nz,sFact,term,volMaxGlob,volMaxLoc,volMinGlob, &
                  volMinLoc
  REAL(RFREAL) :: dXyz(XCOORD:ZCOORD)
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pRhs,pXyz,pXyzOld
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid,pGridOld
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_MoveGridXyz.F90,v $ $Revision: 1.9 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_MoveGridXyz',&
  'RFLU_MoveGridXyz.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel /= VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Moving grid based on coordinates...'
  END IF ! global%myProcid
  
! ******************************************************************************
! Set variables
! ******************************************************************************
    
  nIter = regions(1)%mixtInput%moveGridNIter
  sFact = regions(1)%mixtInput%moveGridSFact 

! TEMPORARY Disable smoothing for MPI version
  nIter = 0 
! END TEMPORARY
  
! ******************************************************************************
! Copy grid data to old grid
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal 
    pGrid    => regions(iReg)%grid
    pGridOld => regions(iReg)%gridOld    
    
    IF ( context == MOVEGRID_CONTEXT_MOVESMOOTH ) THEN 
      DO ic = 1,pGrid%nCellsTot ! Explicit copy to avoid ASCI White problem
        pGridOld%vol(ic) = pGrid%vol(ic)
      END DO ! ic  
    END IF ! context   
         
    DO iv = 1,pGrid%nVertTot ! Explicit copy to avoid ASCI White problem    
      pGridOld%xyz(XCOORD,iv) = pGrid%xyz(XCOORD,iv)
      pGridOld%xyz(YCOORD,iv) = pGrid%xyz(YCOORD,iv)
      pGridOld%xyz(ZCOORD,iv) = pGrid%xyz(ZCOORD,iv)            
    END DO ! iv 
  END DO ! iReg
  
! ******************************************************************************
! Allocate vertex degree array
! ******************************************************************************  
  
  DO iReg = 1,global%nRegionsLocal
    pGrid => regions(iReg)%grid
       
    ALLOCATE(pGrid%degr(pGrid%nVertTot),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%degr')
    END IF ! global%error  
        
    IF ( pGrid%nTetsTot == pGrid%nCellsTot ) THEN ! Only for tetrahedra
      ALLOCATE(pGrid%volMin(pGrid%nVertTot),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%volMin')
      END IF ! global%error 
    END IF ! pGrid%nTetsTot    
  END DO ! iReg  
    
! ******************************************************************************  
! Compute modified edge degree. NOTE the modified vertex degree depends on the
! geometric properties of the grid; in this particular version on the volume. 
! Because the volume does not get updated during the iterations to smooth the
! grid, the modified edge degree does not have to be recomputed. 
! ******************************************************************************  
  
! ==============================================================================  
! Original vertex degree. NOTE must take into account virtual edges.
! ==============================================================================  

  DO iReg = 1,global%nRegionsLocal
    pGrid => regions(iReg)%grid

    DO iv = 1,pGrid%nVertTot ! Explicit loop to avoid ASCI White problem 
      pGrid%degr(iv) = 0
    END DO ! iv 

    DO ie = 1,pGrid%nEdgesTot
      v1 = pGrid%e2v(1,ie)
      v2 = pGrid%e2v(2,ie)

      pGrid%degr(v1) = pGrid%degr(v1) + 1
      pGrid%degr(v2) = pGrid%degr(v2) + 1      
    END DO ! ie
  END DO ! iReg     

! ==============================================================================  
! Modified vertex degree (for tetrahedral grids only). NOTE the modification 
! is usually beneficial if the grid contains very small and very large cells 
! and the small cells are of poor quality, which can lead to negative volumes.
! ==============================================================================  
        
  DO iReg = 1,global%nRegionsLocal 
    pGrid => regions(iReg)%grid

    IF ( pGrid%nTetsTot == pGrid%nCellsTot ) THEN 
      volMinLoc = MINVAL(pGrid%vol(1:pGrid%nTets))
      volMaxLoc = MAXVAL(pGrid%vol(1:pGrid%nTets))

! TO DO 
! Need to use reductions to find volume extrema
      volMinGlob = volMinLoc
      volMaxGlob = volMaxLoc  
! END TO DO 
   
      DO iv = 1,pGrid%nVertTot ! Explicit loop to avoid ASCI White problem 
        pGrid%volMin(iv) = HUGE(1.0_RFREAL)
      END DO ! iv

      DO ic = 1,pGrid%nTetsTot
        v1 = pGrid%tet2v(1,ic)
        v2 = pGrid%tet2v(2,ic)
        v3 = pGrid%tet2v(3,ic)
        v4 = pGrid%tet2v(4,ic) 

        pGrid%volMin(v1) = MIN(pGrid%volMin(v1),pGrid%vol(ic))
        pGrid%volMin(v2) = MIN(pGrid%volMin(v2),pGrid%vol(ic))
        pGrid%volMin(v3) = MIN(pGrid%volMin(v3),pGrid%vol(ic))            
        pGrid%volMin(v4) = MIN(pGrid%volMin(v4),pGrid%vol(ic))
      END DO ! ic

      DO iv = 1,pGrid%nVert
        pGrid%degr(iv) = pGrid%degr(iv) & 
                       + (volMaxGlob - volMinGlob)/pGrid%volMin(iv)
      END DO ! iv        
    END IF ! pGrid%nTetsTot
  END DO ! iReg       
    
! ******************************************************************************
! Initialize coordinates to old coordinates plus patch deformation
! ******************************************************************************
  
  DO iReg = 1,global%nRegionsLocal 
    pGrid    => regions(iReg)%grid
    pGridOld => regions(iReg)%gridOld     
  
    pXyz    => pGrid%xyz
    pXyzOld => pGridOld%xyz    
    
    IF ( context == MOVEGRID_CONTEXT_MOVESMOOTH ) THEN 
      DO iPatch = 1,pGrid%nPatches
        pPatch => regions(iReg)%patches(iPatch)

        IF ( pPatch%movePatchDir /= MOVEPATCH_DIR_NONE ) THEN       
          DO iv = 1,pPatch%nBVert
            ivg = pPatch%bv(iv)            

            pXyz(XCOORD,ivg) = pXyzOld(XCOORD,ivg) + pPatch%dXyz(XCOORD,iv)
            pXyz(YCOORD,ivg) = pXyzOld(YCOORD,ivg) + pPatch%dXyz(YCOORD,iv)
            pXyz(ZCOORD,ivg) = pXyzOld(ZCOORD,ivg) + pPatch%dXyz(ZCOORD,iv)
          END DO ! iv
        END IF ! pPatch%movePatchDir
      END DO ! iPatch
    ELSE IF ( context == MOVEGRID_CONTEXT_ONLYSMOOTH ) THEN 
      DO iPatch = 1,pGrid%nPatches
        pPatch => regions(iReg)%patches(iPatch)

        IF ( pPatch%movePatchDir /= MOVEPATCH_DIR_NONE ) THEN       
          DO iv = 1,pPatch%nBVert
            ivg = pPatch%bv(iv)

            pXyz(XCOORD,ivg) = pXyzOld(XCOORD,ivg)
            pXyz(YCOORD,ivg) = pXyzOld(YCOORD,ivg)
            pXyz(ZCOORD,ivg) = pXyzOld(ZCOORD,ivg)
          END DO ! iv
        END IF ! pPatch%movePatchDir
      END DO ! iPatch    
    END IF ! context
  END DO ! iReg      
     
! ******************************************************************************
! Apply pseudo-Laplacian smoothing to coordinates of all vertices
! ******************************************************************************

  DO it = 1,nIter
             
! ==============================================================================
!   Compute right-hand side with boundary restrictions already imposed
! ==============================================================================  
   
    DO iReg = 1,global%nRegionsLocal  
      pGrid => regions(iReg)%grid
      pXyz  => pGrid%xyz
      pRhs  => pGrid%rhs     
      
! ------------------------------------------------------------------------------
!     Initialize right-hand side
! ------------------------------------------------------------------------------
      
      DO iv = 1,pGrid%nVertTot ! Explicit loop to avoid ASCI White problem 
        pRhs(XCOORD,iv) = 0.0_RFREAL
        pRhs(YCOORD,iv) = 0.0_RFREAL
        pRhs(ZCOORD,iv) = 0.0_RFREAL                
      END DO ! iv        
      
! ------------------------------------------------------------------------------
!     Compute right-hand side. NOTE loop only over actual edges and take into
!     accont the region degree of each edge to deal with multiple contributions
!     of edges on region boundaries.
! ------------------------------------------------------------------------------
  
      DO ie = 1,pGrid%nEdges 
        v1 = pGrid%e2v(1,ie)
        v2 = pGrid%e2v(2,ie)

        term = 1.0_RFREAL/(REAL(pGrid%e2rDegr(ie),KIND=RFREAL))

        dXyz(XCOORD) = term*(pXyz(XCOORD,v2) - pXyz(XCOORD,v1))
        dXyz(YCOORD) = term*(pXyz(YCOORD,v2) - pXyz(YCOORD,v1))
        dXyz(ZCOORD) = term*(pXyz(ZCOORD,v2) - pXyz(ZCOORD,v1)) 

        pRhs(XCOORD,v1) = pRhs(XCOORD,v1) + dXyz(XCOORD)
        pRhs(YCOORD,v1) = pRhs(YCOORD,v1) + dXyz(YCOORD)
        pRhs(ZCOORD,v1) = pRhs(ZCOORD,v1) + dXyz(ZCOORD)     

        pRhs(XCOORD,v2) = pRhs(XCOORD,v2) - dXyz(XCOORD)
        pRhs(YCOORD,v2) = pRhs(YCOORD,v2) - dXyz(YCOORD)
        pRhs(ZCOORD,v2) = pRhs(ZCOORD,v2) - dXyz(ZCOORD)         
      END DO ! ie

! ------------------------------------------------------------------------------
!     Apply boundary deformation restrictions
! ------------------------------------------------------------------------------  
  
      patchCheckCounter = 0
  
      IF ( context == MOVEGRID_CONTEXT_MOVESMOOTH ) THEN 
      
! ----- Enforce imposed deformation --------------------------------------------
      
        DO iPatch = 1,pGrid%nPatches
          pPatch => regions(iReg)%patches(iPatch)

          IF ( pPatch%movePatchDir /= MOVEPATCH_DIR_NONE ) THEN                  
            patchCheckCounter = patchCheckCounter + 1   
                   
            DO iv = 1,pPatch%nBVert
              ivg = pPatch%bv(iv)

              pRhs(XCOORD,ivg) = 0.0_RFREAL
              pRhs(YCOORD,ivg) = 0.0_RFREAL
              pRhs(ZCOORD,ivg) = 0.0_RFREAL
            END DO ! iv
          END IF ! pPatch%movePatch
        END DO ! iPatch                   
        
! ----- Ensure boundary integrity (no imposed motion or smoothing) -------------   
        
        DO iPatch = 1,pGrid%nPatches
          pPatch => regions(iReg)%patches(iPatch)

          IF ( pPatch%movePatchDir == MOVEPATCH_DIR_NONE ) THEN                 
            patchCheckCounter = patchCheckCounter + 1                  
                 
            SELECT CASE ( pPatch%moveBcType ) 
            
! ----------- Normal displacement set to zero with Neumann bc              
            
              CASE ( MOVEGRID_BCTYPE_NEUMANN )
                DO iv = 1,pPatch%nBVert
                  ivg = pPatch%bv(iv)

                  nx = pPatch%bvn(XCOORD,iv)
                  ny = pPatch%bvn(YCOORD,iv)
                  nz = pPatch%bvn(ZCOORD,iv)                    

                  term = pRhs(XCOORD,ivg)*nx & 
                       + pRhs(YCOORD,ivg)*ny & 
                       + pRhs(ZCOORD,ivg)*nz

                  pRhs(XCOORD,ivg) = pRhs(XCOORD,ivg) - term*nx
                  pRhs(YCOORD,ivg) = pRhs(YCOORD,ivg) - term*ny
                  pRhs(ZCOORD,ivg) = pRhs(ZCOORD,ivg) - term*nz            
                END DO ! iv
                
! ----------- Normal displacement set to zero with Dirichlet bc                 
                
              CASE ( MOVEGRID_BCTYPE_DIRICHX, & 
                     MOVEGRID_BCTYPE_DIRICHY, & 
                     MOVEGRID_BCTYPE_DIRICHZ )
                DO iv = 1,pPatch%nBVert
                  ivg = pPatch%bv(iv)
                  
                  pRhs(pPatch%moveBcType,ivg) = 0.0_RFREAL               
                END DO ! iv

! ----------- No bc (must only occur with no movement and smoothing)

              CASE ( MOVEGRID_BCTYPE_NONE ) 
                IF ( pPatch%smoothGrid .EQV. .FALSE. ) THEN 
                  DO iv = 1,pPatch%nBVert
                    ivg = pPatch%bv(iv)

                    pRhs(XCOORD,ivg) = 0.0_RFREAL
                    pRhs(YCOORD,ivg) = 0.0_RFREAL
                    pRhs(ZCOORD,ivg) = 0.0_RFREAL
                  END DO ! iv                  
                ELSE ! Defensive programming
                  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
                END IF ! 

              CASE DEFAULT ! Defensive programming
                CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
            END SELECT ! pPatch%moveBcType 
                      
          END IF ! pPatch
        END DO ! iPatch
        
      ELSE IF ( context == MOVEGRID_CONTEXT_ONLYSMOOTH ) THEN 
        DO iPatch = 1,pGrid%nPatches
          pPatch => regions(iReg)%patches(iPatch)
       
          patchCheckCounter = patchCheckCounter + 1 
       
          DO iv = 1,pPatch%nBVert
            ivg = pPatch%bv(iv)

            pRhs(XCOORD,ivg) = 0.0_RFREAL
            pRhs(YCOORD,ivg) = 0.0_RFREAL
            pRhs(ZCOORD,ivg) = 0.0_RFREAL
          END DO ! iv
        END DO ! iPatch      
      END IF ! context
    END DO ! iReg

! ------------------------------------------------------------------------------
!   Check that all patches have mesh motion boundary conditions set. This 
!   should never happen - if it does, there is a bug in the above code 
!   segment (defensive coding)
! ------------------------------------------------------------------------------  
            
    IF ( patchCheckCounter /= pGrid%nPatches ) THEN 
      CALL ErrorStop(global,ERR_MOVEPATCH_BC_NOTSET,__LINE__)
    END IF ! patchCheckCounter   

! TO DO     
! Need to use reduction to add contributions at interface vertices 
! END TO DO 
              
! ==============================================================================
!   Update coordinates of actual vertices
! ==============================================================================

    DO iReg = 1,global%nRegionsLocal
      pGrid => regions(iReg)%grid
      pXyz  => pGrid%xyz
      pRhs  => pGrid%rhs    
    
      DO iv = 1,pGrid%nVert
        term = sFact/REAL(pGrid%degr(iv),RFREAL)     
      
        pXyz(XCOORD,iv) = pXyz(XCOORD,iv) + term*pRhs(XCOORD,iv)
        pXyz(YCOORD,iv) = pXyz(YCOORD,iv) + term*pRhs(YCOORD,iv)
        pXyz(ZCOORD,iv) = pXyz(ZCOORD,iv) + term*pRhs(ZCOORD,iv)
      END DO ! iv
    END DO ! iReg

! TO DO 
! Update virtual vertices
! END TO DO 
  END DO ! it

! ******************************************************************************
! Deallocate residual and degree arrays
! ******************************************************************************  
  
  DO iReg = 1,global%nRegionsLocal 
    pGrid => regions(iReg)%grid
        
    DEALLOCATE(pGrid%degr,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%degr')
    END IF ! global%error 
    
    IF ( pGrid%nTetsTot == pGrid%nCellsTot ) THEN ! Only for tetrahedra
      DEALLOCATE(pGrid%volMin,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%volMin')
      END IF ! global%error 
    END IF ! pGrid%nTetsTot             
  END DO ! iReg  

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel /= VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Moving grid based on coordinates done.'
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_MoveGridXyz

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_MoveGridXyz.F90,v $
! Revision 1.9  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.6  2005/06/09 20:22:58  haselbac
! Replaced movePatch by movePatchDir
!
! Revision 1.5  2005/04/15 15:07:21  haselbac
! Converted to MPI, grid motion can only be used for serial runs
!
! Revision 1.4  2003/08/20 02:09:58  haselbac
! Changed verbosity conditions to reduce solver output in GENx runs
!
! Revision 1.3  2003/07/09 22:36:53  haselbac
! Added IF statement to avoid probs with single-proc charm jobs
!
! Revision 1.2  2003/05/01 14:12:03  haselbac
! Added logic to deal with non-moving and non-smoothed patches
!
! Revision 1.1  2003/03/31 16:05:39  haselbac
! Initial revision
!
! Revision 1.8  2003/03/27 19:31:36  haselbac
! Fixed bug in initializing pGrid%degr
!
! Revision 1.7  2003/03/15 18:51:03  haselbac
! Completed || of grid motion
!
! Revision 1.6  2003/02/20 19:43:51  haselbac
! Some simplification, added additional check
!
! Revision 1.5  2003/01/28 14:43:19  haselbac
! Now get user input, enforce no smoothing if necessary, small area modification activated
!
! Revision 1.4  2002/11/26 15:28:35  haselbac
! Added Dirichlet bc for grid motion
!
! Revision 1.3  2002/11/19 23:46:41  haselbac
! Cleaned up boundary integrity enforcement of grid motion
!
! Revision 1.2  2002/11/08 21:33:11  haselbac
! Several modifications, mainly cosmetic apart from degree-modification
!
! Revision 1.1  2002/10/27 19:15:20  haselbac
! Initial revision
!
! Revision 1.1  2002/10/16 21:15:11  haselbac
! Initial revision
!
! ******************************************************************************







