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
! Purpose: Suite of routines related to WENO method.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModWENO.F90,v 1.4 2008/12/06 08:44:25 mtcampbe Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModWENO

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid  
  USE ModBndPatch, ONLY: t_patch
  USE ModMPI
    
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_WENOGradCellsWrapper, &
            RFLU_WENOGradCellsXYZWrapper
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModWENO.F90,v $ $Revision: 1.4 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  








! ******************************************************************************
!
! Purpose: Compute 2d WENO gradients of any vector or scalar at cell centers. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!
! Output:
!   grad        Gradients of variables at cell centers
!
! Notes: None.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_WENOGradCells_2D(pRegion,iBegGrad,iEndGrad,grad)
      
    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: iBegGrad,iEndGrad
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg,icg2,iGrad,isl
    REAL(RFREAL) :: smooIndSum,term
    REAL(RFREAL) :: gradLocal(XCOORD:YCOORD)
    REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: smooInd  
    REAL(RFREAL), DIMENSION(:,:,:), ALLOCATABLE :: gradENO
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WENOGradCells_2D',&
  'RFLU_ModWENO.F90' )

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::WENOGradCells_2D")
#endif

! ******************************************************************************    
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Compute ENO gradients
! ******************************************************************************

! ==============================================================================
!   Allocate memory         
! ==============================================================================

    ALLOCATE(gradENO(XCOORD:YCOORD,iBegGrad:iEndGrad,pGrid%nCellsTot), &
             STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gradENO')
    END IF ! global%error 

    ALLOCATE(smooInd(0:pGrid%c2csInfo%nCellMembsMax),STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'smooInd')
    END IF ! global%error

! ==============================================================================
!   Loop over cells 
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot

! ------------------------------------------------------------------------------
!     Compute weighted gradient based on smoothness indicator
! ------------------------------------------------------------------------------

      DO iGrad = iBegGrad,iEndGrad

! ----- Compute smoothness indicator -------------------------------------------

        term = ABS(grad(XCOORD,iGrad,icg)) + ABS(grad(YCOORD,iGrad,icg)) 

        smooInd(0) = 1.0_RFREAL/(term*term + 1.0E-15_RFREAL)
        smooIndSum = smooInd(0)

        DO isl = 1,pGrid%c2cs(icg)%nCellMembs                   
          icg2 = pGrid%c2cs(icg)%cellMembs(isl)        

          term = ABS(grad(XCOORD,iGrad,icg2)) + ABS(grad(YCOORD,iGrad,icg2))                                                                      

          smooInd(isl) = 1.0_RFREAL/(term*term + 1.0E-15_RFREAL)          
          smooIndSum   = smooIndSum + smooInd(isl)                             
        END DO ! isl

! ----- Compute weighted gradient ----------------------------------------------

        term = smooInd(0)/smooIndSum  

        gradLocal(XCOORD) = term*grad(XCOORD,iGrad,icg)  
        gradLocal(YCOORD) = term*grad(YCOORD,iGrad,icg)

        DO isl = 1,pGrid%c2cs(icg)%nCellMembs
          icg2 = pGrid%c2cs(icg)%cellMembs(isl) 
          term = smooInd(isl)/smooIndSum

          gradLocal(XCOORD) = gradLocal(XCOORD) + term*grad(XCOORD,iGrad,icg2)
          gradLocal(YCOORD) = gradLocal(YCOORD) + term*grad(YCOORD,iGrad,icg2)       
        END DO ! isl

        gradENO(XCOORD,iGrad,icg) = gradLocal(XCOORD)
        gradENO(YCOORD,iGrad,icg) = gradLocal(YCOORD)
      END DO ! iGrad                                                       
    END DO ! icg

! ==============================================================================
!   Reassign gradients
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot
      DO iGrad = iBegGrad,iEndGrad
        grad(XCOORD,iGrad,icg) = gradENO(XCOORD,iGrad,icg)
        grad(YCOORD,iGrad,icg) = gradENO(YCOORD,iGrad,icg)       
      END DO ! iGrad  
    END DO ! icg

! ==============================================================================
!   Deallocate memory         
! ==============================================================================

    DEALLOCATE(smooInd,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'smooInd')
    END IF ! global%error
      
    DEALLOCATE(gradENO,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gradENO')
    END IF ! global%error 

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::WENOGradCells_2D")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WENOGradCells_2D







! ******************************************************************************
!
! Purpose: Compute 3d WENO gradients of any vector or scalar at cell centers. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!
! Output:
!   grad        Gradients of variables at cell centers
!
! Notes: None.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_WENOGradCells_3D(pRegion,iBegGrad,iEndGrad,grad)
      
    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: iBegGrad,iEndGrad
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg,icg2,iGrad,isl
    REAL(RFREAL) :: smooIndSum,term
    REAL(RFREAL) :: gradLocal(XCOORD:ZCOORD)
    REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: smooInd  
    REAL(RFREAL), DIMENSION(:,:,:), ALLOCATABLE :: gradENO
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WENOGradCells_3D',&
  'RFLU_ModWENO.F90' )

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::WENOGradCells_3D")
#endif

! ******************************************************************************    
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Compute ENO gradients
! ******************************************************************************

! ==============================================================================
!   Allocate memory         
! ==============================================================================

    ALLOCATE(gradENO(XCOORD:ZCOORD,iBegGrad:iEndGrad,pGrid%nCellsTot), &
             STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gradENO')
    END IF ! global%error 

    ALLOCATE(smooInd(0:pGrid%c2csInfo%nCellMembsMax),STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'smooInd')
    END IF ! global%error

! ==============================================================================
!   Loop over cells 
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot

! ------------------------------------------------------------------------------
!     Compute weighted gradient based on smoothness indicator
! ------------------------------------------------------------------------------

      DO iGrad = iBegGrad,iEndGrad

! ----- Compute smoothness indicator -------------------------------------------

        term = ABS(grad(XCOORD,iGrad,icg)) & 
             + ABS(grad(YCOORD,iGrad,icg)) & 
             + ABS(grad(ZCOORD,iGrad,icg)) 

        smooInd(0) = 1.0_RFREAL/(term*term + 1.0E-15_RFREAL)
        smooIndSum = smooInd(0)

        DO isl = 1,pGrid%c2cs(icg)%nCellMembs                   
          icg2 = pGrid%c2cs(icg)%cellMembs(isl)        

          term = ABS(grad(XCOORD,iGrad,icg2)) & 
               + ABS(grad(YCOORD,iGrad,icg2)) & 
               + ABS(grad(ZCOORD,iGrad,icg2))                                                                      

          smooInd(isl) = 1.0_RFREAL/(term*term + 1.0E-15_RFREAL)          
          smooIndSum   = smooIndSum + smooInd(isl)                             
        END DO ! isl

! ----- Compute weighted gradient ----------------------------------------------

        term = smooInd(0)/smooIndSum  

        gradLocal(XCOORD) = term*grad(XCOORD,iGrad,icg)  
        gradLocal(YCOORD) = term*grad(YCOORD,iGrad,icg)
        gradLocal(ZCOORD) = term*grad(ZCOORD,iGrad,icg)

        DO isl = 1,pGrid%c2cs(icg)%nCellMembs
          icg2 = pGrid%c2cs(icg)%cellMembs(isl) 
          term = smooInd(isl)/smooIndSum

          gradLocal(XCOORD) = gradLocal(XCOORD) + term*grad(XCOORD,iGrad,icg2)
          gradLocal(YCOORD) = gradLocal(YCOORD) + term*grad(YCOORD,iGrad,icg2)
          gradLocal(ZCOORD) = gradLocal(ZCOORD) + term*grad(ZCOORD,iGrad,icg2)        
        END DO ! isl

        gradENO(XCOORD,iGrad,icg) = gradLocal(XCOORD)
        gradENO(YCOORD,iGrad,icg) = gradLocal(YCOORD)
        gradENO(ZCOORD,iGrad,icg) = gradLocal(ZCOORD) 
      END DO ! iGrad                                                       
    END DO ! icg

! ==============================================================================
!   Reassign gradients
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot
      DO iGrad = iBegGrad,iEndGrad
        grad(XCOORD,iGrad,icg) = gradENO(XCOORD,iGrad,icg)
        grad(YCOORD,iGrad,icg) = gradENO(YCOORD,iGrad,icg)
        grad(ZCOORD,iGrad,icg) = gradENO(ZCOORD,iGrad,icg)            
      END DO ! iGrad  
    END DO ! icg

! ==============================================================================
!   Deallocate memory         
! ==============================================================================

    DEALLOCATE(smooInd,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'smooInd')
    END IF ! global%error
      
    DEALLOCATE(gradENO,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gradENO')
    END IF ! global%error 

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::WENOGradCells_3D")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WENOGradCells_3D






! ******************************************************************************
!
! Purpose: Compute WENO gradients of any vector or scalar at cell centers. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!   grad        Gradients of variables at cell centers
!
! Output:
!   grad        Weighted gradients of variables at cell centers
!
! Notes: 
!   1. If stencil dimension is 1, calling xyz routine, anything else does not
!      make sense.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_WENOGradCellsWrapper(pRegion,iBegGrad,iEndGrad,grad)

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iBegGrad,iEndGrad
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_global), POINTER :: global

! *****************************************************************************
!   Start
! *****************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WENOGradCellsWrapper',&
  'RFLU_ModWENO.F90' )

! ******************************************************************************    
!   Call weighting routines
! ******************************************************************************    

    SELECT CASE ( pRegion%mixtInput%stencilDimensCells )
      CASE ( 1 ) ! NOTE calling xyz routine, anything else does not make sense
        CALL RFLU_WENOGradCellsXYZ_1D(pRegion,iBegGrad,iEndGrad,grad)               
      CASE ( 2 ) 
        CALL RFLU_WENOGradCells_2D(pRegion,iBegGrad,iEndGrad,grad)  
      CASE ( 3 ) 
        CALL RFLU_WENOGradCells_3D(pRegion,iBegGrad,iEndGrad,grad)          
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pMixtInput%stencilDimensCells

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WENOGradCellsWrapper






! ******************************************************************************
!
! Purpose: Compute 1d WENO gradients of any vector or scalar at cell centers 
!   component-wise. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!
! Output:
!   grad        Gradients of variables at cell centers
!
! Notes: None.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_WENOGradCellsXYZ_1D(pRegion,iBegGrad,iEndGrad,grad)
      
    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: iBegGrad,iEndGrad
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg,icg2,iDir,iDirEnd,iGrad,isl
    REAL(RFREAL) :: smooIndSum,term
    REAL(RFREAL), DIMENSION(XCOORD:ZCOORD) :: gradLocal
    REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: smooInd  
    REAL(RFREAL), DIMENSION(:,:,:), ALLOCATABLE :: gradENO
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WENOGradCellsXYZ_1D',&
  'RFLU_ModWENO.F90' )

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::WENOGradCellsXYZ_1D")
#endif

! ******************************************************************************    
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

    SELECT CASE ( pRegion%mixtInput%dimens ) 
      CASE ( 1 ) 
        iDirEnd = XCOORD
      CASE ( 2 )
        iDirEnd = YCOORD
      CASE ( 3 )       
        iDirEnd = ZCOORD        
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%dimens

! ******************************************************************************
!   Compute ENO gradients
! ******************************************************************************

! ==============================================================================
!   Allocate memory         
! ==============================================================================

    ALLOCATE(gradENO(XCOORD:iDirEnd,iBegGrad:iEndGrad,pGrid%nCellsTot), &
             STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gradENO')
    END IF ! global%error 

    ALLOCATE(smooInd(0:pGrid%c2csInfo%nCellMembsMax),STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'smooInd')
    END IF ! global%error

! ==============================================================================
!   Loop over cells 
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot

! ------------------------------------------------------------------------------
!     Compute weighted gradient based on smoothness indicator
! ------------------------------------------------------------------------------

      DO iGrad = iBegGrad,iEndGrad
      
! ----- Loop over directions ---------------------------------------------------      
      
        DO iDir = XCOORD,iDirEnd     

! ------- Compute smoothness indicator 

          term = grad(iDir,iGrad,icg)

          smooInd(0) = 1.0_RFREAL/(term*term + 1.0E-15_RFREAL)
          smooIndSum = smooInd(0)

          DO isl = 1,pGrid%c2cs1D(iDir,icg)%nCellMembs                   
            icg2 = pGrid%c2cs1D(iDir,icg)%cellMembs(isl)        

            term = grad(iDir,iGrad,icg2)

            smooInd(isl) = 1.0_RFREAL/(term*term + 1.0E-15_RFREAL)
            smooIndSum   = smooIndSum + smooInd(isl)
          END DO ! isl

! ------- Compute weighted gradient 

          term = smooInd(0)/smooIndSum

          gradLocal(iDir) = term*grad(iDir,iGrad,icg)  

          DO isl = 1,pGrid%c2cs1D(iDir,icg)%nCellMembs
            icg2 = pGrid%c2cs1D(iDir,icg)%cellMembs(isl) 

            term = smooInd(isl)/smooIndSum

            gradLocal(iDir) = gradLocal(iDir) + term*grad(iDir,iGrad,icg2)       
          END DO ! isl

          gradENO(iDir,iGrad,icg) = gradLocal(iDir)
        END DO ! iDir
      END DO ! iGrad                                                       
    END DO ! icg

! ==============================================================================
!   Reassign gradients
! ==============================================================================

    SELECT CASE ( pRegion%mixtInput%dimens ) 
      CASE ( 1 ) 
        DO icg = 1,pGrid%nCellsTot
          DO iGrad = iBegGrad,iEndGrad
            grad(XCOORD,iGrad,icg) = gradENO(XCOORD,iGrad,icg)
            grad(YCOORD,iGrad,icg) = 0.0_RFREAL  
            grad(ZCOORD,iGrad,icg) = 0.0_RFREAL           
          END DO ! iGrad  
        END DO ! icg
      CASE ( 2 )
        DO icg = 1,pGrid%nCellsTot
          DO iGrad = iBegGrad,iEndGrad
            grad(XCOORD,iGrad,icg) = gradENO(XCOORD,iGrad,icg)
            grad(YCOORD,iGrad,icg) = gradENO(YCOORD,iGrad,icg)
            grad(ZCOORD,iGrad,icg) = 0.0_RFREAL           
          END DO ! iGrad  
        END DO ! icg
      CASE ( 3 )       
        DO icg = 1,pGrid%nCellsTot
          DO iGrad = iBegGrad,iEndGrad
            grad(XCOORD,iGrad,icg) = gradENO(XCOORD,iGrad,icg)
            grad(YCOORD,iGrad,icg) = gradENO(YCOORD,iGrad,icg)
            grad(ZCOORD,iGrad,icg) = gradENO(ZCOORD,iGrad,icg)            
          END DO ! iGrad  
        END DO ! icg  
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%dimens

! ==============================================================================
!   Deallocate memory         
! ==============================================================================

    DEALLOCATE(smooInd,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'smooInd')
    END IF ! global%error
      
    DEALLOCATE(gradENO,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gradENO')
    END IF ! global%error 

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::WENOGradCellsXYZ_1D")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WENOGradCellsXYZ_1D






! ******************************************************************************
!
! Purpose: Compute 2d WENO gradients of any vector or scalar at cell centers 
!   component-wise. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!
! Output:
!   grad        Gradients of variables at cell centers
!
! Notes: None.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_WENOGradCellsXYZ_2D(pRegion,iBegGrad,iEndGrad,grad)
      
    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: iBegGrad,iEndGrad
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg,icg2,iGrad,isl
    REAL(RFREAL) :: termX,termY
    REAL(RFREAL), DIMENSION(XCOORD:YCOORD) :: gradLocal,smooIndSum
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: smooInd  
    REAL(RFREAL), DIMENSION(:,:,:), ALLOCATABLE :: gradENO
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WENOGradCellsXYZ_2D',&
  'RFLU_ModWENO.F90' )

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::WENOGradCellsXYZ_2D")
#endif

! ******************************************************************************    
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Compute ENO gradients
! ******************************************************************************

! ==============================================================================
!   Allocate memory         
! ==============================================================================

    ALLOCATE(gradENO(XCOORD:YCOORD,iBegGrad:iEndGrad,pGrid%nCellsTot), &
             STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gradENO')
    END IF ! global%error 

    ALLOCATE(smooInd(XCOORD:YCOORD,0:pGrid%c2csInfo%nCellMembsMax), & 
             STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'smooInd')
    END IF ! global%error

! ==============================================================================
!   Loop over cells 
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot

! ------------------------------------------------------------------------------
!     Compute weighted gradient based on smoothness indicator
! ------------------------------------------------------------------------------

      DO iGrad = iBegGrad,iEndGrad

! ----- Compute smoothness indicator -------------------------------------------

        termX = grad(XCOORD,iGrad,icg)
        termY = grad(YCOORD,iGrad,icg)

        smooInd(XCOORD,0) = 1.0_RFREAL/(termX*termX + 1.0E-15_RFREAL)
        smooInd(YCOORD,0) = 1.0_RFREAL/(termY*termY + 1.0E-15_RFREAL)

        smooIndSum(XCOORD) = smooInd(XCOORD,0)
        smooIndSum(YCOORD) = smooInd(YCOORD,0)

        DO isl = 1,pGrid%c2cs(icg)%nCellMembs                   
          icg2 = pGrid%c2cs(icg)%cellMembs(isl)        

          termX = grad(XCOORD,iGrad,icg2)
          termY = grad(YCOORD,iGrad,icg2)

          smooInd(XCOORD,isl) = 1.0_RFREAL/(termX*termX + 1.0E-15_RFREAL)
          smooInd(YCOORD,isl) = 1.0_RFREAL/(termY*termY + 1.0E-15_RFREAL)

          smooIndSum(XCOORD) = smooIndSum(XCOORD) + smooInd(XCOORD,isl)
          smooIndSum(YCOORD) = smooIndSum(YCOORD) + smooInd(YCOORD,isl)
        END DO ! isl

! ----- Compute weighted gradient ----------------------------------------------

        termX = smooInd(XCOORD,0)/smooIndSum(XCOORD)  
        termY = smooInd(YCOORD,0)/smooIndSum(YCOORD)  

        gradLocal(XCOORD) = termX*grad(XCOORD,iGrad,icg)  
        gradLocal(YCOORD) = termY*grad(YCOORD,iGrad,icg)

        DO isl = 1,pGrid%c2cs(icg)%nCellMembs
          icg2 = pGrid%c2cs(icg)%cellMembs(isl) 

          termX = smooInd(XCOORD,isl)/smooIndSum(XCOORD)
          termY = smooInd(YCOORD,isl)/smooIndSum(YCOORD)

          gradLocal(XCOORD) = gradLocal(XCOORD) + termX*grad(XCOORD,iGrad,icg2)
          gradLocal(YCOORD) = gradLocal(YCOORD) + termY*grad(YCOORD,iGrad,icg2)       
        END DO ! isl

        gradENO(XCOORD,iGrad,icg) = gradLocal(XCOORD)
        gradENO(YCOORD,iGrad,icg) = gradLocal(YCOORD)
      END DO ! iGrad                                                       
    END DO ! icg

! ==============================================================================
!   Reassign gradients
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot
      DO iGrad = iBegGrad,iEndGrad
        grad(XCOORD,iGrad,icg) = gradENO(XCOORD,iGrad,icg)
        grad(YCOORD,iGrad,icg) = gradENO(YCOORD,iGrad,icg)      
      END DO ! iGrad  
    END DO ! icg

! DEBUG
!    DO iGrad = iBegGrad,iEndGrad
!      WRITE(*,*) iGrad,MINVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		       MINVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		       MAXVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		       MAXVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot))
!    END DO ! iGrad
!    
!    STOP
! END DEBUG

! ==============================================================================
!   Deallocate memory         
! ==============================================================================

    DEALLOCATE(smooInd,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'smooInd')
    END IF ! global%error
      
    DEALLOCATE(gradENO,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gradENO')
    END IF ! global%error 

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::WENOGradCellsXYZ_2D")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WENOGradCellsXYZ_2D








! ******************************************************************************
!
! Purpose: Compute 3d WENO gradients of any vector or scalar at cell centers 
!   component-wise. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!
! Output:
!   grad        Gradients of variables at cell centers
!
! Notes: None.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_WENOGradCellsXYZ_3D(pRegion,iBegGrad,iEndGrad,grad)
      
    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: iBegGrad,iEndGrad
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg,icg2,iGrad,isl
    REAL(RFREAL) :: termX,termY,termZ
    REAL(RFREAL), DIMENSION(XCOORD:ZCOORD) :: gradLocal,smooIndSum
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: smooInd  
    REAL(RFREAL), DIMENSION(:,:,:), ALLOCATABLE :: gradENO
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WENOGradCellsXYZ_3D',&
  'RFLU_ModWENO.F90' )

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::WENOGradCellsXYZ_3D")
#endif

! ******************************************************************************    
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Compute ENO gradients
! ******************************************************************************

! ==============================================================================
!   Allocate memory         
! ==============================================================================

    ALLOCATE(gradENO(XCOORD:ZCOORD,iBegGrad:iEndGrad,pGrid%nCellsTot), &
             STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gradENO')
    END IF ! global%error 

    ALLOCATE(smooInd(XCOORD:ZCOORD,0:pGrid%c2csInfo%nCellMembsMax), & 
             STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'smooInd')
    END IF ! global%error

! ==============================================================================
!   Loop over cells 
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot

! ------------------------------------------------------------------------------
!     Compute weighted gradient based on smoothness indicator
! ------------------------------------------------------------------------------

      DO iGrad = iBegGrad,iEndGrad

! ----- Compute smoothness indicator -------------------------------------------

        termX = grad(XCOORD,iGrad,icg)
        termY = grad(YCOORD,iGrad,icg)
        termZ = grad(ZCOORD,iGrad,icg)

        smooInd(XCOORD,0) = 1.0_RFREAL/(termX*termX + 1.0E-15_RFREAL)
        smooInd(YCOORD,0) = 1.0_RFREAL/(termY*termY + 1.0E-15_RFREAL)
        smooInd(ZCOORD,0) = 1.0_RFREAL/(termZ*termZ + 1.0E-15_RFREAL)

        smooIndSum(XCOORD) = smooInd(XCOORD,0)
        smooIndSum(YCOORD) = smooInd(YCOORD,0)
        smooIndSum(ZCOORD) = smooInd(ZCOORD,0)

        DO isl = 1,pGrid%c2cs(icg)%nCellMembs                   
          icg2 = pGrid%c2cs(icg)%cellMembs(isl)        

          termX = grad(XCOORD,iGrad,icg2)
          termY = grad(YCOORD,iGrad,icg2)
          termZ = grad(ZCOORD,iGrad,icg2)

          smooInd(XCOORD,isl) = 1.0_RFREAL/(termX*termX + 1.0E-15_RFREAL)
          smooInd(YCOORD,isl) = 1.0_RFREAL/(termY*termY + 1.0E-15_RFREAL)
          smooInd(ZCOORD,isl) = 1.0_RFREAL/(termZ*termZ + 1.0E-15_RFREAL)

          smooIndSum(XCOORD) = smooIndSum(XCOORD) + smooInd(XCOORD,isl)
          smooIndSum(YCOORD) = smooIndSum(YCOORD) + smooInd(YCOORD,isl)
          smooIndSum(ZCOORD) = smooIndSum(ZCOORD) + smooInd(ZCOORD,isl)
        END DO ! isl

! ----- Compute weighted gradient ----------------------------------------------

        termX = smooInd(XCOORD,0)/smooIndSum(XCOORD)  
        termY = smooInd(YCOORD,0)/smooIndSum(YCOORD)  
        termZ = smooInd(ZCOORD,0)/smooIndSum(ZCOORD)  

        gradLocal(XCOORD) = termX*grad(XCOORD,iGrad,icg)  
        gradLocal(YCOORD) = termY*grad(YCOORD,iGrad,icg)
        gradLocal(ZCOORD) = termZ*grad(ZCOORD,iGrad,icg)

        DO isl = 1,pGrid%c2cs(icg)%nCellMembs
          icg2 = pGrid%c2cs(icg)%cellMembs(isl) 

          termX = smooInd(XCOORD,isl)/smooIndSum(XCOORD)
          termY = smooInd(YCOORD,isl)/smooIndSum(YCOORD)
          termZ = smooInd(ZCOORD,isl)/smooIndSum(ZCOORD)

          gradLocal(XCOORD) = gradLocal(XCOORD) + termX*grad(XCOORD,iGrad,icg2)
          gradLocal(YCOORD) = gradLocal(YCOORD) + termY*grad(YCOORD,iGrad,icg2)
          gradLocal(ZCOORD) = gradLocal(ZCOORD) + termZ*grad(ZCOORD,iGrad,icg2)        
        END DO ! isl

        gradENO(XCOORD,iGrad,icg) = gradLocal(XCOORD)
        gradENO(YCOORD,iGrad,icg) = gradLocal(YCOORD)
        gradENO(ZCOORD,iGrad,icg) = gradLocal(ZCOORD) 
      END DO ! iGrad                                                       
    END DO ! icg

! ==============================================================================
!   Reassign gradients
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot
      DO iGrad = iBegGrad,iEndGrad
        grad(XCOORD,iGrad,icg) = gradENO(XCOORD,iGrad,icg)
        grad(YCOORD,iGrad,icg) = gradENO(YCOORD,iGrad,icg)
        grad(ZCOORD,iGrad,icg) = gradENO(ZCOORD,iGrad,icg)            
      END DO ! iGrad  
    END DO ! icg

! DEBUG
!    DO iGrad = iBegGrad,iEndGrad
!      WRITE(*,*) iGrad,MINVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		       MINVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot)), &
!		       MINVAL(grad(ZCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		       MAXVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		       MAXVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot)), &
!		       MAXVAL(grad(ZCOORD,iGrad,1:pGrid%nCellsTot))
!    END DO ! iGrad
!    
!    STOP
! END DEBUG

! ==============================================================================
!   Deallocate memory         
! ==============================================================================

    DEALLOCATE(smooInd,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'smooInd')
    END IF ! global%error
      
    DEALLOCATE(gradENO,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gradENO')
    END IF ! global%error 

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::WENOGradCellsXYZ_3D")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WENOGradCellsXYZ_3D








! ******************************************************************************
!
! Purpose: Compute 2d ENO gradients of any vector or scalar at cell centers 
!   component wise. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!
! Output:
!   grad        Gradients of variables at cell centers
!
! Notes: 
!   1. Optimized by Adam Moody and Charles Shereda, LLNL.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_WENOGradCellsXYZFast_2D(pRegion,iBegGrad,iEndGrad,grad)
      
    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: iBegGrad,iEndGrad
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg,icg2,iGrad,isl,nCellMembs
    INTEGER, DIMENSION(:), ALLOCATABLE :: icg_ary    
    REAL(RFREAL) :: gradLocalX,gradLocalY,nextX,nextX2,nextY,nextY2, &
                    smooIndSumX,smooIndX,smooIndSumY,smooIndY,termX,termY
    REAL(RFREAL), DIMENSION(XCOORD:YCOORD) :: gradLocal,smooIndSum
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: smooInd  
    REAL(RFREAL), DIMENSION(:,:,:), ALLOCATABLE :: gradENO
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WENOGradCellsXYZFast_2D',&
  'RFLU_ModWENO.F90' )

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::WENOGradCellsXYZFast_2D")
#endif

! ******************************************************************************    
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Compute ENO gradients
! ******************************************************************************

! ==============================================================================
!   Allocate memory         
! ==============================================================================

    ALLOCATE(gradENO(XCOORD:YCOORD,iBegGrad:iEndGrad,pGrid%nCellsTot), &
             STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gradENO')
    END IF ! global%error 

    ALLOCATE(icg_ary(0:pGrid%c2csInfo%nCellMembsMax+2),STAT=errorFlag )
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'icg_ary')
    END IF ! global%error

! ==============================================================================
!   Loop over cells 
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot
      nCellMembs = pGrid%c2cs(icg)%nCellMembs
      
      DO isl = 1,nCellMembs
        icg_ary(isl) = pGrid%c2cs(icg)%cellMembs(isl)
      END DO ! isl
      
      icg_ary(nCellMembs+1) = icg
      icg_ary(nCellMembs+2) = icg

      DO iGrad = iBegGrad,iEndGrad
        termX = grad(XCOORD,iGrad,icg)
        termY = grad(YCOORD,iGrad,icg)
        
        nextX = grad(XCOORD,iGrad,icg_ary(1))
        nextY = grad(YCOORD,iGrad,icg_ary(1))
        
        nextX2 = grad(XCOORD,iGrad,icg_ary(2))
        nextY2 = grad(YCOORD,iGrad,icg_ary(2))

        smooIndSumX = 1.0_RFREAL/(termX*termX + 1.0E-15_RFREAL)
        smooIndSumY = 1.0_RFREAL/(termY*termY + 1.0E-15_RFREAL)

        gradLocalX = smooIndSumX*termX
        gradLocalY = smooIndSumY*termY

        DO isl = 1,nCellMembs                  
          icg2 = icg_ary(isl+2)

          termX = nextX
          termY = nextY

          nextX = nextX2
          nextY = nextY2

          nextX2 = grad(XCOORD,iGrad,icg2)
          nextY2 = grad(YCOORD,iGrad,icg2)

          smooIndX = 1.0_RFREAL/(termX*termX + 1.0E-15_RFREAL)
          smooIndY = 1.0_RFREAL/(termY*termY + 1.0E-15_RFREAL)

          smooIndSumX = smooIndSumX + smooIndX
          smooIndSumY = smooIndSumY + smooIndY

          gradLocalX = gradLocalX + smooIndX*termX
          gradLocalY = gradLocalY + smooIndY*termY
        END DO ! isl

        gradENO(XCOORD,iGrad,icg) = gradLocalX/smooIndSumX
        gradENO(YCOORD,iGrad,icg) = gradLocalY/smooIndSumY
      END DO ! iGrad                                                       
    END DO ! icg

    DO icg = 1,pGrid%nCellsTot
      DO iGrad = iBegGrad,iEndGrad
        grad(XCOORD,iGrad,icg) = gradENO(XCOORD,iGrad,icg)
        grad(YCOORD,iGrad,icg) = gradENO(YCOORD,iGrad,icg)           
      END DO ! iGrad  
    END DO ! icg

! DEBUG
!    DO iGrad = iBegGrad,iEndGrad
!      WRITE(*,*) iGrad,MINVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		       MINVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		       MAXVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		       MAXVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot))
!    END DO ! iGrad
!    
!    STOP
! END DEBUG

! ==============================================================================
!   Deallocate memory         
! ==============================================================================

    DEALLOCATE(icg_ary,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'smooInd')
    END IF ! global%error
      
    DEALLOCATE(gradENO,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gradENO')
    END IF ! global%error 

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::WENOGradCellsXYZFast_2D")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WENOGradCellsXYZFast_2D







! ******************************************************************************
!
! Purpose: Compute 3d ENO gradients of any vector or scalar at cell centers 
!   component wise. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!
! Output:
!   grad        Gradients of variables at cell centers
!
! Notes: 
!   1. Optimized by Adam Moody and Charles Shereda, LLNL.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_WENOGradCellsXYZFast_3D(pRegion,iBegGrad,iEndGrad,grad)
      
    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER :: iBegGrad,iEndGrad
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,icg,icg2,iGrad,isl,nCellMembs
    INTEGER, DIMENSION(:), ALLOCATABLE :: icg_ary    
    REAL(RFREAL) :: gradLocalX,gradLocalY,gradLocalZ,nextX,nextX2,nextY, &
                    nextY2,nextZ,nextZ2,smooIndSumX,smooIndX,smooIndSumY, &
                    smooIndY,smooIndSumZ,smooIndZ,termX,termY,termZ
    REAL(RFREAL), DIMENSION(XCOORD:ZCOORD) :: gradLocal,smooIndSum
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: smooInd  
    REAL(RFREAL), DIMENSION(:,:,:), ALLOCATABLE :: gradENO
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WENOGradCellsXYZFast_3D',&
  'RFLU_ModWENO.F90' )

#ifdef ROCPROF
    CALL FPROFILER_BEGINS("RFLU::WENOGradCellsXYZFast_3D")
#endif

! ******************************************************************************    
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Compute ENO gradients
! ******************************************************************************

! ==============================================================================
!   Allocate memory         
! ==============================================================================

    ALLOCATE(gradENO(XCOORD:ZCOORD,iBegGrad:iEndGrad,pGrid%nCellsTot), &
             STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'gradENO')
    END IF ! global%error 

    ALLOCATE(icg_ary(0:pGrid%c2csInfo%nCellMembsMax+2),STAT=errorFlag )
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'icg_ary')
    END IF ! global%error

! ==============================================================================
!   Loop over cells 
! ==============================================================================

    DO icg = 1,pGrid%nCellsTot
      nCellMembs = pGrid%c2cs(icg)%nCellMembs
      
      DO isl = 1,nCellMembs
        icg_ary(isl) = pGrid%c2cs(icg)%cellMembs(isl)
      END DO ! isl
      
      icg_ary(nCellMembs+1) = icg
      icg_ary(nCellMembs+2) = icg

      DO iGrad = iBegGrad,iEndGrad
	termX = grad(XCOORD,iGrad,icg)
	termY = grad(YCOORD,iGrad,icg)
	termZ = grad(ZCOORD,iGrad,icg)
                     
	nextX = grad(XCOORD,iGrad,icg_ary(1))
	nextY = grad(YCOORD,iGrad,icg_ary(1))
	nextZ = grad(ZCOORD,iGrad,icg_ary(1))
       
	nextX2 = grad(XCOORD,iGrad,icg_ary(2))
	nextY2 = grad(YCOORD,iGrad,icg_ary(2))
	nextZ2 = grad(ZCOORD,iGrad,icg_ary(2))

        smooIndSumX = 1.0_RFREAL/(termX*termX + 1.0E-15_RFREAL)
        smooIndSumY = 1.0_RFREAL/(termY*termY + 1.0E-15_RFREAL)
        smooIndSumZ = 1.0_RFREAL/(termZ*termZ + 1.0E-15_RFREAL)

        gradLocalX = smooIndSumX*termX
        gradLocalY = smooIndSumY*termY
        gradLocalZ = smooIndSumZ*termZ

        DO isl = 1,nCellMembs                  
          icg2 = icg_ary(isl+2)

          termX = nextX
          termY = nextY
          termZ = nextZ

          nextX = nextX2
          nextY = nextY2
          nextZ = nextZ2

          nextX2 = grad(XCOORD,iGrad,icg2)
          nextY2 = grad(YCOORD,iGrad,icg2)
          nextZ2 = grad(ZCOORD,iGrad,icg2)

          smooIndX = 1.0_RFREAL/(termX*termX + 1.0E-15_RFREAL)
          smooIndY = 1.0_RFREAL/(termY*termY + 1.0E-15_RFREAL)
          smooIndZ = 1.0_RFREAL/(termZ*termZ + 1.0E-15_RFREAL)

          smooIndSumX = smooIndSumX + smooIndX
          smooIndSumY = smooIndSumY + smooIndY
          smooIndSumZ = smooIndSumZ + smooIndZ

          gradLocalX = gradLocalX + smooIndX*termX
          gradLocalY = gradLocalY + smooIndY*termY
          gradLocalZ = gradLocalZ + smooIndZ*termZ
        END DO ! isl

        gradENO(XCOORD,iGrad,icg) = gradLocalX/smooIndSumX
        gradENO(YCOORD,iGrad,icg) = gradLocalY/smooIndSumY
        gradENO(ZCOORD,iGrad,icg) = gradLocalZ/smooIndSumZ
      END DO ! iGrad                                                       
    END DO ! icg

    DO icg = 1,pGrid%nCellsTot
      DO iGrad = iBegGrad,iEndGrad
        grad(XCOORD,iGrad,icg) = gradENO(XCOORD,iGrad,icg)
        grad(YCOORD,iGrad,icg) = gradENO(YCOORD,iGrad,icg)
        grad(ZCOORD,iGrad,icg) = gradENO(ZCOORD,iGrad,icg)            
      END DO ! iGrad  
    END DO ! icg

! DEBUG
!    DO iGrad = iBegGrad,iEndGrad
!      WRITE(*,*) iGrad,MINVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		       MINVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot)), &
!		       MINVAL(grad(ZCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		       MAXVAL(grad(XCOORD,iGrad,1:pGrid%nCellsTot)), & 
!		       MAXVAL(grad(YCOORD,iGrad,1:pGrid%nCellsTot)), &
!		       MAXVAL(grad(ZCOORD,iGrad,1:pGrid%nCellsTot))
!    END DO ! iGrad
!    
!    STOP
! END DEBUG

! ==============================================================================
!   Deallocate memory         
! ==============================================================================

    DEALLOCATE(icg_ary,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'smooInd')
    END IF ! global%error
      
    DEALLOCATE(gradENO,STAT=errorFlag)
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'gradENO')
    END IF ! global%error 

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF
    CALL FPROFILER_ENDS("RFLU::WENOGradCellsXYZFast_3D")
#endif

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WENOGradCellsXYZFast_3D
  





! ******************************************************************************
!
! Purpose: Compute WENO-gradients of any vector or scalar at cell centers 
!  component-wise. 
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region data
!   iBegGrad    Beginning index of data in grad
!   iEndGrad    Ending index of data in grad
!   grad        Gradients of variables at cell centers
!
! Output:
!   grad        Weighted gradients of variables at cell centers
!
! Notes: None.
!
! ******************************************************************************
    
  SUBROUTINE RFLU_WENOGradCellsXYZWrapper(pRegion,iBegGrad,iEndGrad,grad)

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iBegGrad,iEndGrad
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: grad
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    TYPE(t_global), POINTER :: global

! *****************************************************************************
!   Start
! *****************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WENOGradCellsXYZWrapper',&
  'RFLU_ModWENO.F90' )

! ******************************************************************************    
!   Call weighting routines
! ******************************************************************************    

    SELECT CASE ( pRegion%mixtInput%stencilDimensCells )
      CASE ( 1 ) 
        CALL RFLU_WENOGradCellsXYZ_1D(pRegion,iBegGrad,iEndGrad,grad)               
      CASE ( 2 ) 
!        CALL RFLU_WENOGradCellsXYZ_2D(pRegion,iBegGrad,iEndGrad,grad)
        CALL RFLU_WENOGradCellsXYZFast_2D(pRegion,iBegGrad,iEndGrad,grad)	  
      CASE ( 3 ) 
!        CALL RFLU_WENOGradCellsXYZ_3D(pRegion,iBegGrad,iEndGrad,grad)
        CALL RFLU_WENOGradCellsXYZFast_3D(pRegion,iBegGrad,iEndGrad,grad)                    
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pMixtInput%stencilDimensCells

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WENOGradCellsXYZWrapper





! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModWENO


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModWENO.F90,v $
! Revision 1.4  2008/12/06 08:44:25  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:36  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/19 19:43:17  haselbac
! Added tuned routines
!
! Revision 1.1  2006/04/07 14:36:19  haselbac
! Initial revision
!
! ******************************************************************************















