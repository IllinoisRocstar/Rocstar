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
! Purpose: Move grid by smoothing displacements.
!
! Description: Simple pseudo-Laplacian smoothing of grid given boundary 
!   distortion. 
!
! Input: 
!   regions     Region data
!
! Output: None.
!
! Notes: 
!   1. The two parameters influencing the grid motion, i.e., the number of
!      iterations and the smoothing coefficient, may need some tuning as more
!      complicated geometries are attempted...
!   2. The integrity of patches with no imposed motion is ensured by one of 
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
! $Id: RFLU_MoveGridDisp.F90,v 1.11 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_MoveGridDisp(regions)

  USE ModDataTypes
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters

  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,ic,ie,iPatch,iReg,it,iv,ivg,nIter,patchCheckCounter, & 
             v1,v2
  REAL(RFREAL) :: nx,ny,nz,sFact,term,x1,x2,y1,y2,z1,z2
  REAL(RFREAL) :: dDisp(XCOORD:ZCOORD)
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pDisp,pRhs,pXyz
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid,pGridOld
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_MoveGridDisp.F90,v $ $Revision: 1.11 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_MoveGridDisp',&
  'RFLU_MoveGridDisp.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel /= VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Moving grid based on displacements...'
  END IF ! global%myProcid
  
! *****************************************************************************
! Set variables
! *****************************************************************************
    
  nIter = regions(1)%mixtInput%moveGridNIter
  sFact = regions(1)%mixtInput%moveGridSFact 

! TEMPORARY disable smoothing for MPI version
  nIter = 0
! END TEMPORARY
  
! *****************************************************************************
! Copy grid data to old grid. Must be done even if have zero smoothing 
! iterations. 
! *****************************************************************************

  DO iReg = 1,global%nRegionsLocal 
    pGrid    => regions(iReg)%grid
    pGridOld => regions(iReg)%gridOld    
    
    DO ic = 1,pGrid%nCellsTot ! Explicit copy to avoid ASCI White problem
      pGridOld%vol(ic) = pGrid%vol(ic)
    END DO ! ic  
         
    DO iv = 1,pGrid%nVertTot ! Explicit copy to avoid ASCI White problem    
      pGridOld%xyz(XCOORD,iv) = pGrid%xyz(XCOORD,iv)
      pGridOld%xyz(YCOORD,iv) = pGrid%xyz(YCOORD,iv)
      pGridOld%xyz(ZCOORD,iv) = pGrid%xyz(ZCOORD,iv)            
    END DO ! iv 
  END DO ! iReg
  
! *****************************************************************************
! Allocate temporary memory for and compute edge weights 
! *****************************************************************************  
    
  IF ( nIter > 0 ) THEN 
    DO iReg = 1,global%nRegionsLocal
      pGrid => regions(iReg)%grid
       
      ALLOCATE(pGrid%gmEdgeWght(pGrid%nEdgesTot),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%gmEdgeWght')
      END IF ! global%error     )   
       
      ALLOCATE(pGrid%gmVertWght(pGrid%nVertTot),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pGrid%gmVertWght')
      END IF ! global%error     
    END DO ! iReg  
    
    DO iReg = 1,global%nRegionsLocal
      pGrid => regions(iReg)%grid
  
      DO iv = 1,pGrid%nVertTot
        pGrid%gmVertWght(iv) = 0.0_RFREAL
      END DO ! iv

      DO ie = 1,pGrid%nEdgesTot
        v1 = pGrid%e2v(1,ie)
        v2 = pGrid%e2v(2,ie)

        x1 = pGrid%xyz(XCOORD,v1)
        y1 = pGrid%xyz(YCOORD,v1)
        z1 = pGrid%xyz(ZCOORD,v1)
      
        x2 = pGrid%xyz(XCOORD,v2)
        y2 = pGrid%xyz(YCOORD,v2)
        z2 = pGrid%xyz(ZCOORD,v2)                  

        term = 1.0_RFREAL/SQRT((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
  
        pGrid%gmEdgeWght(ie) = term
      
        pGrid%gmVertWght(v1) = pGrid%gmVertWght(v1) + term
        pGrid%gmVertWght(v2) = pGrid%gmVertWght(v2) + term    
      END DO ! ie    
    END DO ! iReg     
  END IF ! nIter
    
! ******************************************************************************
! Initialize displacements to zero and impose patch displacements. NOTE need to
! impose patch displacements so that get boundary deformation even if have zero
! smoothing iterations. 
! ******************************************************************************    
  
  DO iReg = 1,global%nRegionsLocal 
    pGrid => regions(iReg)%grid
    pDisp => pGrid%disp  
    
    DO iv = 1,pGrid%nVertTot
      pDisp(XCOORD,iv) = 0.0_RFREAL 
      pDisp(YCOORD,iv) = 0.0_RFREAL       
      pDisp(ZCOORD,iv) = 0.0_RFREAL          
    END DO ! iv 
    
    DO iPatch = 1,pGrid%nPatches
      pPatch => regions(iReg)%patches(iPatch)
 
      IF ( pPatch%movePatchDir /= MOVEPATCH_DIR_NONE ) THEN       
        DO iv = 1,pPatch%nBVert
          ivg = pPatch%bv(iv)            
 
          pDisp(XCOORD,ivg) = pPatch%dXyz(XCOORD,iv)
          pDisp(YCOORD,ivg) = pPatch%dXyz(YCOORD,iv)
          pDisp(ZCOORD,ivg) = pPatch%dXyz(ZCOORD,iv)
        END DO ! iv
      END IF ! pPatch%movePatchDir
    END DO ! iPatch         
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
      pDisp => pGrid%disp 
      pRhs  => pGrid%rhs           
      
      DO iv = 1,pGrid%nVert
        pRhs(XCOORD,iv) = 0.0_RFREAL  
        pRhs(YCOORD,iv) = 0.0_RFREAL              
        pRhs(ZCOORD,iv) = 0.0_RFREAL  
      END DO ! iv    
      
! ------------------------------------------------------------------------------
!     Enforce patch deformation. NOTE that strictly speaking this is not
!     needed during first smoothing iteration because initialized outside 
!     iteration loop. 
! ------------------------------------------------------------------------------      
            
      DO iPatch = 1,pGrid%nPatches
        pPatch => regions(iReg)%patches(iPatch)

        IF ( pPatch%movePatchDir /= MOVEPATCH_DIR_NONE ) THEN       
          DO iv = 1,pPatch%nBVert
            ivg = pPatch%bv(iv)            

            pDisp(XCOORD,ivg) = pPatch%dXyz(XCOORD,iv)
            pDisp(YCOORD,ivg) = pPatch%dXyz(YCOORD,iv)
            pDisp(ZCOORD,ivg) = pPatch%dXyz(ZCOORD,iv)
          END DO ! iv
        END IF ! pPatch%movePatchDir
      END DO ! iPatch      
      
! ------------------------------------------------------------------------------
!     Compute right-hand side. NOTE loop only over actual edges and take into
!     account the region degree of each edge to deal with multiple contributions
!     of edges on region boundaries.
! ------------------------------------------------------------------------------
  
      DO ie = 1,pGrid%nEdges 
        v1 = pGrid%e2v(1,ie)
        v2 = pGrid%e2v(2,ie)

        term = sFact*pGrid%gmEdgeWght(ie)/(REAL(pGrid%e2rDegr(ie),KIND=RFREAL))

        dDisp(XCOORD) = term*(pDisp(XCOORD,v2) - pDisp(XCOORD,v1))
        dDisp(YCOORD) = term*(pDisp(YCOORD,v2) - pDisp(YCOORD,v1))
        dDisp(ZCOORD) = term*(pDisp(ZCOORD,v2) - pDisp(ZCOORD,v1)) 

        pRhs(XCOORD,v1) = pRhs(XCOORD,v1) + dDisp(XCOORD)
        pRhs(YCOORD,v1) = pRhs(YCOORD,v1) + dDisp(YCOORD)
        pRhs(ZCOORD,v1) = pRhs(ZCOORD,v1) + dDisp(ZCOORD)     

        pRhs(XCOORD,v2) = pRhs(XCOORD,v2) - dDisp(XCOORD)
        pRhs(YCOORD,v2) = pRhs(YCOORD,v2) - dDisp(YCOORD)
        pRhs(ZCOORD,v2) = pRhs(ZCOORD,v2) - dDisp(ZCOORD)         
      END DO ! ie

! ------------------------------------------------------------------------------
!     Apply boundary deformation restrictions
! ------------------------------------------------------------------------------  
  
      patchCheckCounter = 0
        
! --- Enforce imposed deformation ----------------------------------------------     
      
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
        END IF ! pPatch%movePatchDir
      END DO ! iPatch                   

! --- Ensure boundary integrity (no imposed motion or smoothing) ---------------   

      DO iPatch = 1,pGrid%nPatches
        pPatch => regions(iReg)%patches(iPatch)

        IF ( pPatch%movePatchDir == MOVEPATCH_DIR_NONE ) THEN                 
          patchCheckCounter = patchCheckCounter + 1                  

          SELECT CASE ( pPatch%moveBcType ) 

! --------- Normal displacement set to zero with Neumann bc              

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

! --------- Normal displacement set to zero with Dirichlet bc                 

            CASE ( MOVEGRID_BCTYPE_DIRICHX, & 
                   MOVEGRID_BCTYPE_DIRICHY, & 
                   MOVEGRID_BCTYPE_DIRICHZ )                                             
              DO iv = 1,pPatch%nBVert
                ivg = pPatch%bv(iv)

                pRhs(pPatch%moveBcType,ivg) = 0.0_RFREAL               
              END DO ! iv 

! --------- No bc (must only occur with no movement and smoothing)

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

        END IF ! pPatch%movePatchDir
      END DO ! iPatch                   
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
! Reduction of rhs 
! END TO DO 
                        
! ==============================================================================
!   Update displacements
! ==============================================================================

    DO iReg = 1,global%nRegionsLocal
      pGrid => regions(iReg)%grid
      pDisp => pGrid%disp
      pRhs  => pGrid%rhs    
    
      DO iv = 1,pGrid%nVertTot
        term = 1.0_RFREAL/REAL(pGrid%gmVertWght(iv),RFREAL)     
      
        pDisp(XCOORD,iv) = pDisp(XCOORD,iv) + term*pRhs(XCOORD,iv)
        pDisp(YCOORD,iv) = pDisp(YCOORD,iv) + term*pRhs(YCOORD,iv)
        pDisp(ZCOORD,iv) = pDisp(ZCOORD,iv) + term*pRhs(ZCOORD,iv)
      END DO ! iv      
    END DO ! iReg

! TO DO 
! Update virtual vertices
! END TO DO 
  END DO ! it

! ******************************************************************************
! Update coordinates
! ******************************************************************************

  IF ( nIter > 0 ) THEN 
    DO iReg = 1,global%nRegionsLocal 
      pGrid => regions(iReg)%grid
      pDisp => pGrid%disp
      pXyz  => pGrid%xyz
    
      DO iv = 1,pGrid%nVertTot
        pXyz(XCOORD,iv) = pXyz(XCOORD,iv) + pDisp(XCOORD,iv)
        pXyz(YCOORD,iv) = pXyz(YCOORD,iv) + pDisp(YCOORD,iv)
        pXyz(ZCOORD,iv) = pXyz(ZCOORD,iv) + pDisp(ZCOORD,iv)            
      END DO ! iv
    END DO ! iReg
  END IF ! nIter

! ******************************************************************************
! Deallocate temporary arrays
! ******************************************************************************  
  
  IF ( nIter > 0 ) THEN 
    DO iReg = 1,global%nRegionsLocal 
      pGrid => regions(iReg)%grid
          
      DEALLOCATE(pGrid%gmEdgeWght,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%gmEdgeWght')
      END IF ! global%error        
          
      DEALLOCATE(pGrid%gmVertWght,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pGrid%gmVertWght')
      END IF ! global%error               
    END DO ! iReg  
  END IF ! nIter

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel /= VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                             'Moving grid based on displacements done.'
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_MoveGridDisp

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_MoveGridDisp.F90,v $
! Revision 1.11  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2005/06/09 20:22:58  haselbac
! Replaced movePatch by movePatchDir
!
! Revision 1.8  2005/04/15 15:07:20  haselbac
! Converted to MPI, grid motion can only be used for serial runs
!
! Revision 1.7  2004/09/17 19:59:58  haselbac
! Added IF statements to prevent operations for zero iterations
!
! Revision 1.6  2004/04/14 03:57:11  haselbac
! Bug fix: Now get correct boundary deformation for zero smoothing iterations
!
! Revision 1.5  2003/08/20 02:09:58  haselbac
! Changed verbosity conditions to reduce solver output in GENx runs
!
! Revision 1.4  2003/07/09 22:36:53  haselbac
! Added IF statement to avoid probs with single-proc charm jobs
!
! Revision 1.3  2003/05/01 14:12:10  haselbac
! Added logic to deal with non-moving and non-smoothed patches
!
! Revision 1.2  2003/04/17 00:16:19  haselbac
! Bug fix: Weight should be inverse sqrt, little effect
!
! Revision 1.1  2003/03/31 16:05:39  haselbac
! Initial revision
!
! ******************************************************************************







