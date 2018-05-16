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
! Purpose: Suite of routines for MPI interaction.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModMPI.F90,v 1.17 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModMPI

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt_input
  USE ModBorder, ONLY: t_border,t_border_data
  USE ModMPI

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_MPI_ClearRequestWrapper, &
            RFLU_MPI_CopyWrapper, & 
            RFLU_MPI_CreateBufferIPclSend, &        
            RFLU_MPI_CreateBuffersWrapper, &    
            RFLU_MPI_DestroyBufferIPclSend, &           
            RFLU_MPI_DestroyBuffersWrapper, &        
            RFLU_MPI_ISendWrapper, & 
            RFLU_MPI_RecreateBufferIPclSend, &
            RFLU_MPI_RecvWrapper, & 
            RFLU_MPI_SetTagsWrapper
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
      
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModMPI.F90,v $ $Revision: 1.17 $' 
                      
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS






! ******************************************************************************
!
! Purpose: Clearing send requests.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   request     Request (to be cleared)
!
! Output: 
!   request     Request (cleared)
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_ClearRequest(global,request)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(INOUT) :: request
    TYPE(t_global), POINTER :: global  
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag
    INTEGER :: status(MPI_STATUS_SIZE)          
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_ClearRequest',&
  'RFLU_ModMPI.F90')

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    CALL MPI_Wait(request,status,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_ClearRequest









! ******************************************************************************
!
! Purpose: Wrapper for clearing send requests.
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

  SUBROUTINE RFLU_MPI_ClearRequestWrapper(pRegion)

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

    INTEGER :: errorFlag,iBorder
    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid             
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_ClearRequestWrapper',&
  'RFLU_ModMPI.F90')

#ifdef ROCPROF 
    CALL FPROFILER_BEGINS("RFLU::ClearRequestWrapper")
#endif

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
    
! ==============================================================================
!     Clear request if not on same process
! ==============================================================================    
    
      IF ( pBorder%iProc /= global%myProcid ) THEN 
        IF ( pBorder%nCellsSend > 0 ) THEN 

! ------------------------------------------------------------------------------
!         Mixture
! ------------------------------------------------------------------------------      

          CALL RFLU_MPI_ClearRequest(global,pBorder%mixt%sendRequest)

! ------------------------------------------------------------------------------
!         Physical modules
! ------------------------------------------------------------------------------      

#ifdef SPEC
          IF ( global%specUsed .EQV. .TRUE. ) THEN 
            CALL RFLU_MPI_ClearRequest(global,pBorder%spec%sendRequest)   
          END IF ! global%specUsed
#endif  
        END IF ! pBorder%nCellsSend
      END IF ! pBorder
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF 
    CALL FPROFILER_ENDS("RFLU::ClearRequestWrapper")
#endif

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_ClearRequestWrapper






! ******************************************************************************
!
! Purpose: Copy cell data.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pBorder     Pointer to border
!   cellData    Array with cell data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_CopyCellData(global,pBorder,pBorder2,cellData,cellData2)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cellData
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: cellData2    
    TYPE(t_border), POINTER :: pBorder,pBorder2
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icg2,icl,iVar,nVars        
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_CopyCellData',&
  'RFLU_ModMPI.F90')

! ******************************************************************************
!   Set variables
! ******************************************************************************

    nVars = SIZE(cellData,1)
    
    IF ( nVars /= SIZE(cellData2,1) ) THEN 
      CALL ErrorStop(global,ERR_DATADIM_MISMATCH,__LINE__)
    END IF ! nVars 
    
! ******************************************************************************
!   Copy data 
! ******************************************************************************

    DO icl = 1,pBorder%nCellsSend
      icg  = pBorder%icgSend(icl)
      icg2 = pBorder2%icgRecv(icl)

      DO iVar = 1,nVars
        cellData2(iVar,icg2) = cellData(iVar,icg) 
      END DO ! iVar      
    END DO ! icl
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_CopyCellData









! ******************************************************************************
!
! Purpose: Wrapper for copying data.
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

  SUBROUTINE RFLU_MPI_CopyWrapper(regions)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), DIMENSION(:), POINTER :: regions
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iBorder,iBorder2,iReg,iReg2
    TYPE(t_border), POINTER :: pBorder,pBorder2 
    TYPE(t_global), POINTER :: global 
    TYPE(t_grid), POINTER :: pGrid     
    TYPE(t_region), POINTER :: pRegion,pRegion2
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => regions(0)%global

    CALL RegisterFunction(global,'RFLU_MPI_CopyWrapper',&
  'RFLU_ModMPI.F90')

#ifdef ROCPROF 
    CALL FPROFILER_BEGINS("RFLU::CopyWrapper")
#endif

! ******************************************************************************
!   Set pointers
! ******************************************************************************
    
! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)
      pGrid   => pRegion%grid      

      DO iBorder = 1,pGrid%nBorders
        pBorder => pGrid%borders(iBorder)

! ==============================================================================
!       Copy data if on same process
! ==============================================================================    

        IF ( pBorder%iProc == global%myProcid ) THEN
          iReg2    = pBorder%iRegionLocal
          iBorder2 = pBorder%iBorder

          pRegion2 => regions(iReg2)     
          pBorder2 => pRegion2%grid%borders(iBorder2) 

! ------------------------------------------------------------------------------
!         Check dimensions
! ------------------------------------------------------------------------------

          IF ( pBorder%nCellsSend /= pBorder2%nCellsRecv ) THEN 
            CALL ErrorStop(global,ERR_BUFFERDIM_MISMATCH,__LINE__)
          END IF ! pBorder
          
! ------------------------------------------------------------------------------
!         Mixture
! ------------------------------------------------------------------------------      
      
          CALL RFLU_MPI_CopyCellData(global,pBorder,pBorder2, & 
                                     pRegion%mixt%cv,pRegion2%mixt%cv)      
                                     
! ------------------------------------------------------------------------------
!         Physical modules
! ------------------------------------------------------------------------------      

#ifdef SPEC
          IF ( global%specUsed .EQV. .TRUE. ) THEN 
            CALL RFLU_MPI_CopyCellData(global,pBorder,pBorder2, & 
                                       pRegion%spec%cv,pRegion2%spec%cv)     
          END IF ! global%specUsed
#endif           
                          
        END IF ! pBorder
      END DO ! iBorder
    END DO ! iReg

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF 
    CALL FPROFILER_ENDS("RFLU::CopyWrapper")
#endif

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_CopyWrapper








! ******************************************************************************
!
! Purpose: Create iPclSend buffers.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pBorder     Pointer to border (Optional)
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_CreateBufferIPclSend(pRegion,pBorder)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_border), POINTER, OPTIONAL :: pBorder
    TYPE(t_region), POINTER :: pRegion 
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iBorder,nVars
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid     
    TYPE(t_border), POINTER :: pBorder2
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_CreateBufferIPclSend',&
  'RFLU_ModMPI.F90')

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders when pBorder is not present
!    else allocate only for select pBorder
! ******************************************************************************

#ifdef PLAG
    IF ( PRESENT(pBorder) ) THEN
      nVars = SIZE(pBorder%iPclSend,1)

      ALLOCATE(pBorder%iPclSend(nVars,pBorder%nPclsSendMax),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%iPclSend')
      END IF ! global%error

    ELSE
      DO iBorder = 1,pGrid%nBorders
        pBorder2 => pGrid%borders(iBorder)

        nVars = 2
        pBorder2%nPclsSendMax = 1000

        ALLOCATE(pBorder2%iPclSend(nVars,pBorder2%nPclsSendMax),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pBorder%iPclSend')
        END IF ! global%error
      END DO ! iBorder
    END IF ! PRESENT(pBorder)
#endif

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
      
  END SUBROUTINE RFLU_MPI_CreateBufferIPclSend








! ******************************************************************************
!
! Purpose: Create buffers.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pBorder     Pointer to border
!   borderData  Border data
!   nVars       Number of variables in borderData
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_CreateBuffers(global,pBorder,borderData,nVars)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: nVars
    TYPE(t_global), POINTER :: global
    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_border_data) :: borderData
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag       
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_CreateBuffers',&
  'RFLU_ModMPI.F90')

! ******************************************************************************
!   Allocate memory
! ******************************************************************************
       
    ALLOCATE(borderData%sendBuff(nVars,pBorder%nCellsSend),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'borderData%sendBuff')
    END IF ! global%error

    ALLOCATE(borderData%recvBuff(nVars,pBorder%nCellsRecv),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'borderData%recvBuff')
    END IF ! global%error   

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
      
  END SUBROUTINE RFLU_MPI_CreateBuffers







! ******************************************************************************
!
! Purpose: Wrapper for creating buffers.
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

  SUBROUTINE RFLU_MPI_CreateBuffersWrapper(pRegion)

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

    INTEGER :: errorFlag,iBorder   
    TYPE(t_border), POINTER :: pBorder        
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid 
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_CreateBuffersWrapper',&
  'RFLU_ModMPI.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating buffers...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ==============================================================================
!     Create buffers if not on same process
! ==============================================================================    
      
      IF ( pBorder%iProc /= global%myProcid ) THEN
      
! ------------------------------------------------------------------------------
!       Mixture
! ------------------------------------------------------------------------------      
 
        CALL RFLU_MPI_CreateBuffers(global,pBorder,pBorder%mixt, &
                                    pRegion%mixtInput%nCv)

! ------------------------------------------------------------------------------
!       Physical modules
! ------------------------------------------------------------------------------      

#ifdef SPEC
        IF ( global%specUsed .EQV. .TRUE. ) THEN 
          CALL RFLU_MPI_CreateBuffers(global,pBorder,pBorder%spec, & 
                                      pRegion%specInput%nSpecies)   

        END IF ! global%specUsed        
#endif                           
      END IF ! pBorder%iProc
    END DO ! iBorder
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating buffers done.'
    END IF ! global%verbLevel
      
  END SUBROUTINE RFLU_MPI_CreateBuffersWrapper








! ******************************************************************************
!
! Purpose: Destroy iPclSend buffers.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pBorder     Pointer to border (optional)
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_DestroyBufferIPclSend(pRegion,pBorder)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_border), POINTER, OPTIONAL :: pBorder
    TYPE(t_region), POINTER :: pRegion 
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iBorder      
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid     
    TYPE(t_border), POINTER :: pBorder2
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_DestroyBufferIPclSend',&
  'RFLU_ModMPI.F90')

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders when pBorder is not present
!    else deallocate only for select pBorder
! ******************************************************************************

#ifdef PLAG
    IF ( PRESENT(pBorder) ) THEN
      DEALLOCATE(pBorder%iPclSend,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pBorder%iPclSend')
      END IF ! global%error

    ELSE
      DO iBorder = 1,pGrid%nBorders
        pBorder2 => pGrid%borders(iBorder)

        DEALLOCATE(pBorder2%iPclSend,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pBorder%iPclSend')
        END IF ! global%error
      END DO ! iBorder
    END IF ! PRESENT(pBorder)
#endif

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
      
  END SUBROUTINE RFLU_MPI_DestroyBufferIPclSend








! ******************************************************************************
!
! Purpose: Destroy buffers.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   pBorder     Pointer to border
!   borderData  Border data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_DestroyBuffers(global,pBorder,borderData)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_global), POINTER :: global
    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_border_data) :: borderData
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag       
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_DestroyBuffers',&
  'RFLU_ModMPI.F90')

! ******************************************************************************
!   Allocate memory
! ******************************************************************************
       
    DEALLOCATE(borderData%sendBuff,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'borderData%sendBuff')
    END IF ! global%error

    DEALLOCATE(borderData%recvBuff,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'borderData%recvBuff')
    END IF ! global%error   

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
      
  END SUBROUTINE RFLU_MPI_DestroyBuffers







! ******************************************************************************
!
! Purpose: Wrapper for destroying buffers.
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

  SUBROUTINE RFLU_MPI_DestroyBuffersWrapper(pRegion)

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

    INTEGER :: errorFlag,iBorder   
    TYPE(t_border), POINTER :: pBorder        
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid 
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_DestroyBuffersWrapper',&
  'RFLU_ModMPI.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying buffers...'
    END IF ! global%verbLevel
    
! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ==============================================================================
!     Create buffers if not on same process
! ==============================================================================    
      
      IF ( pBorder%iProc /= global%myProcid ) THEN
      
! ------------------------------------------------------------------------------
!       Mixture
! ------------------------------------------------------------------------------      

        CALL RFLU_MPI_DestroyBuffers(global,pBorder,pBorder%mixt)      

! ------------------------------------------------------------------------------
!       Physical modules
! ------------------------------------------------------------------------------      

#ifdef SPEC
        IF ( global%specUsed .EQV. .TRUE. ) THEN 
          CALL RFLU_MPI_DestroyBuffers(global,pBorder,pBorder%spec)   
        END IF ! global%specUsed        
#endif       
      END IF ! pBorder%iProc
    END DO ! iBorder
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_HIGH ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Destroying buffers done.'
    END IF ! global%verbLevel
      
  END SUBROUTINE RFLU_MPI_DestroyBuffersWrapper









! ******************************************************************************
!
! Purpose: Send cell data.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   pBorder             Pointer to border
!   cellDataBuff        Buffer array
!   cellData            Data array
!   tag                 Tag
!
! Output: 
!   request             Request
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_ISendCellData(global,pBorder,cellDataBuff,cellData,tag, & 
                                    request)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: tag
    INTEGER, INTENT(OUT) :: request
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cellData
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: cellDataBuff
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icl,iVar,nVars        
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_ISendCellData',&
  'RFLU_ModMPI.F90')
    
! ******************************************************************************
!   Set variables
! ******************************************************************************

    nVars = SIZE(cellData,1)

! ******************************************************************************
!   Pack data into buffer
! ******************************************************************************

    DO icl = 1,pBorder%nCellsSend
      icg = pBorder%icgSend(icl)

      DO iVar = 1,nVars
        cellDataBuff(iVar,icl) = cellData(iVar,icg)
      END DO ! iVar      
    END DO ! icl

! ******************************************************************************
!   Send data  
! ******************************************************************************

    IF ( pBorder%nCellsSend > 0 ) THEN 
      CALL MPI_ISend(cellDataBuff,pBorder%nCellsSend*nVars,MPI_RFREAL, & 
                     pBorder%iProc,tag,global%mpiComm,request,errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
      END IF ! global%error                                      
    END IF ! pBorder%nCellsSend

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_ISendCellData








! ******************************************************************************
!
! Purpose: Recreate iPclSend buffers.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!   pBorder     Pointer to border (Optional)
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_RecreateBufferIPclSend(pRegion,pBorder)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_border), POINTER, OPTIONAL :: pBorder
    TYPE(t_region), POINTER :: pRegion 
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,iPcl,iVar,nPclsSendMax,nPclsSendMaxOld,nVars
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iPclSendTemp     
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid     
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_RecreateBufferIPclSend',&
  'RFLU_ModMPI.F90')

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

#ifdef PLAG
! ******************************************************************************
!   Set variables 
! ******************************************************************************

    nVars           = SIZE(pBorder%iPclSend,1)
    nPclsSendMaxOld = SIZE(pBorder%iPclSend,2)

    pBorder%nPclsSendMax = & 
      NINT(1.20_RFREAL*REAL(pBorder%nPclsSend,KIND=RFREAL))

! ******************************************************************************
!   Allocate temporary array
! ******************************************************************************

    ALLOCATE(iPclSendTemp(nVars,pBorder%nPclsSendMax),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'iPclSendTemp')
    END IF ! global%error
    
! ******************************************************************************
!   Copy array
! ******************************************************************************
    
    DO iPcl = 1,nPclsSendMaxOld
      DO iVar = 1,nVars
        iPclSendTemp(iVar,iPcl) = pBorder%iPclSend(iVar,iPcl)
      END DO ! iVar
    END DO ! iPcl 

! ******************************************************************************
!   Deallocate iPclSend array
! ******************************************************************************      

    CALL RFLU_MPI_DestroyBufferIPclSend(pRegion,pBorder)

! ******************************************************************************
!   Reallocate iPclSend array
! ******************************************************************************      

    CALL RFLU_MPI_CreateBufferIPclSend(pRegion,pBorder)   

! ******************************************************************************
!   Copy array
! ******************************************************************************
    
    DO iPcl = 1,pBorder%nPclsSend
      DO iVar = 1,nVars
        pBorder%iPclSend(iVar,iPcl) = iPclSendTemp(iVar,iPcl)
      END DO ! iVar
    END DO ! iPcl     

! ******************************************************************************
!   Deallocate temporary array
! ******************************************************************************

    DEALLOCATE(iPclSendTemp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'iPclSendTemp')
    END IF ! global%error
#endif

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
      
  END SUBROUTINE RFLU_MPI_RecreateBufferIPclSend







! ******************************************************************************
!
! Purpose: Wrapper for sending data.
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

  SUBROUTINE RFLU_MPI_ISendWrapper(pRegion)

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

    INTEGER :: errorFlag,iBorder
    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid             
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_ISendWrapper',&
  'RFLU_ModMPI.F90')

#ifdef ROCPROF 
    CALL FPROFILER_BEGINS("RFLU::ISendWrapper")
#endif

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
    
! ==============================================================================
!     Send data if not on same process
! ==============================================================================    
    
      IF ( pBorder%iProc /= global%myProcid ) THEN 
      
! ------------------------------------------------------------------------------
!       Mixture
! ------------------------------------------------------------------------------      

        CALL RFLU_MPI_ISendCellData(global,pBorder,pBorder%mixt%sendBuff, & 
                                    pRegion%mixt%cv,pBorder%mixt%tag, &
                                    pBorder%mixt%sendRequest)

! ------------------------------------------------------------------------------
!       Physical modules
! ------------------------------------------------------------------------------      

#ifdef SPEC
        IF ( global%specUsed .EQV. .TRUE. ) THEN 
          CALL RFLU_MPI_ISendCellData(global,pBorder,pBorder%spec%sendBuff, & 
                                      pRegion%spec%cv,pBorder%spec%tag, &
                                      pBorder%spec%sendRequest)        
        END IF ! global%specUsed
#endif  
      END IF ! pBorder
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF 
    CALL FPROFILER_ENDS("RFLU::ISendWrapper")
#endif

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_ISendWrapper








! ******************************************************************************
!
! Purpose: Receive cell data.
!
! Description: None.
!
! Input:
!   global              Pointer to global data
!   pBorder             Pointer to border
!   cellDataBuff        Buffer array
!   cellData            Data array
!   tag                 Tag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_MPI_RecvCellData(global,pBorder,cellDataBuff,cellData,tag)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: tag
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: cellDataBuff
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: cellData
    TYPE(t_border), POINTER :: pBorder
    TYPE(t_global), POINTER :: global         
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: errorFlag,icg,icl,iVar,nVars
    INTEGER :: status(MPI_STATUS_SIZE)        
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_RecvCellData',&
  'RFLU_ModMPI.F90')

! ******************************************************************************
!   Set variables
! ******************************************************************************

    nVars = SIZE(cellData,1)
    
! ******************************************************************************
!   Recv data 
! ******************************************************************************

    IF ( pBorder%nCellsRecv > 0 ) THEN 
      CALL MPI_Recv(cellDataBuff,pBorder%nCellsRecv*nVars,MPI_RFREAL, & 
                    pBorder%iProc,tag,global%mpiComm,status,errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
      END IF ! global%error  
    END IF ! pBorder%nCellsRecv

! ******************************************************************************
!   Unpack data from buffer
! ******************************************************************************

    DO icl = 1,pBorder%nCellsRecv
      icg = pBorder%icgRecv(icl)

      DO iVar = 1,nVars
        cellData(iVar,icg) = cellDataBuff(iVar,icl) 
      END DO ! iVar      
    END DO ! icl
        
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_RecvCellData







! ******************************************************************************
!
! Purpose: Wrapper for receiving data.
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

  SUBROUTINE RFLU_MPI_RecvWrapper(pRegion)

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

    INTEGER :: errorFlag,iBorder
    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_global), POINTER :: global 
    TYPE(t_grid), POINTER :: pGrid     
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_RecvWrapper',&
  'RFLU_ModMPI.F90')

#ifdef ROCPROF 
    CALL FPROFILER_BEGINS("RFLU::RecvWrapper")
#endif

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)

! ==============================================================================
!     Receive data if not on same process
! ==============================================================================    

      IF ( pBorder%iProc /= global%myProcid ) THEN 

! ------------------------------------------------------------------------------
!       Mixture
! ------------------------------------------------------------------------------      

        CALL RFLU_MPI_RecvCellData(global,pBorder,pBorder%mixt%recvBuff, & 
                                   pRegion%mixt%cv,pBorder%mixt%tag)      

! ------------------------------------------------------------------------------
!       Physical modules
! ------------------------------------------------------------------------------      

#ifdef SPEC
        IF ( global%specUsed .EQV. .TRUE. ) THEN 
          CALL RFLU_MPI_RecvCellData(global,pBorder,pBorder%spec%recvBuff, & 
                                     pRegion%spec%cv,pBorder%spec%tag)           
        END IF ! global%specUsed
#endif         
      END IF ! pBorder
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

#ifdef ROCPROF 
    CALL FPROFILER_ENDS("RFLU::RecvWrapper")
#endif

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_RecvWrapper





! ******************************************************************************
!
! Purpose: Set tag.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   iReg1       Index of first region
!   iReg2       Index of second region
!   iMsg        Index of message
!   tagMax      Maximum allowed value of tag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  INTEGER FUNCTION RFLU_MPI_SetTag(global,iReg1,iReg2,iMsg,tagMax)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: iMsg,iReg1,iReg2,tagMax
    TYPE(t_global), POINTER :: global 
       
! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER :: iRegMax,iRegMin  
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_MPI_SetTag',&
  'RFLU_ModMPI.F90')

! ******************************************************************************
!   Set tag
! ******************************************************************************

    iRegMax = MAX(iReg1,iReg2)
    iRegMin = MIN(iReg1,iReg2)

    RFLU_MPI_SetTag = iRegMin + (iRegMax-1)*global%nRegions &
                    + (iRegMax-1)*(iMsg-1)*global%nRegions*global%nRegions 

    IF ( RFLU_MPI_SetTag > tagMax ) THEN 
      CALL ErrorStop(global,ERR_MPI_TAGMAX,__LINE__)
    END IF ! RFLU_MPI_SetTag

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END FUNCTION RFLU_MPI_SetTag










! ******************************************************************************
!
! Purpose: Wrapper for setting tags.
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

  SUBROUTINE RFLU_MPI_SetTagsWrapper(pRegion)

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

    LOGICAL :: dummyLogical
    INTEGER :: errorFlag,iBorder,iMsg,tagMax
    TYPE(t_border), POINTER :: pBorder 
    TYPE(t_global), POINTER :: global  
    TYPE(t_grid), POINTER :: pGrid             
  
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_MPI_SetTagsWrapper',&
  'RFLU_ModMPI.F90')

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid

! ******************************************************************************
!   Get maximum allowed value of tag. NOTE must always use MPI_COMM_WORLD - 
!   cannot use global%mpiComm because may be split off in GENx computations
!   and get zero tagMax as a result.
! ******************************************************************************

    CALL MPI_Attr_get(MPI_COMM_WORLD,MPI_TAG_UB,tagMax,dummyLogical,errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
    END IF ! global%error

! ******************************************************************************
!   Loop over borders
! ******************************************************************************

    DO iBorder = 1,pGrid%nBorders
      pBorder => pGrid%borders(iBorder)
    
! ==============================================================================
!     Mixture
! ==============================================================================
            
      iMsg = 1

      pBorder%mixt%tag = RFLU_MPI_SetTag(global,pRegion%iRegionGlobal, & 
                                         pBorder%iRegionGlobal,iMsg,tagMax)

! ==============================================================================
!     Physical modules
! ==============================================================================

#ifdef SPEC
      iMsg = iMsg + 1

      pBorder%spec%tag = RFLU_MPI_SetTag(global,pRegion%iRegionGlobal, & 
                                         pBorder%iRegionGlobal,iMsg,tagMax)            
#endif        

#ifdef PLAG
      iMsg = iMsg + 1

      pBorder%plag%tagCount = RFLU_MPI_SetTag(global,pRegion%iRegionGlobal, & 
                                              pBorder%iRegionGlobal,iMsg,tagMax)            
      iMsg = iMsg + 1

      pBorder%plag%tagInt = RFLU_MPI_SetTag(global,pRegion%iRegionGlobal, & 
                                            pBorder%iRegionGlobal,iMsg,tagMax)            
      iMsg = iMsg + 1

      pBorder%plag%tag = RFLU_MPI_SetTag(global,pRegion%iRegionGlobal, & 
                                         pBorder%iRegionGlobal,iMsg,tagMax)            
#endif        
    END DO ! iBorder

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_MPI_SetTagsWrapper









! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModMPI


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModMPI.F90,v $
! Revision 1.17  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.16  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.15  2007/07/08 21:45:03  gzheng
! changed the PRESENT is used for PGI compiler
!
! Revision 1.14  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.13  2006/02/09 03:37:44  haselbac
! Bug fix: Must always use MPI_COMM_WORLD to get max tag
!
! Revision 1.12  2005/12/14 21:50:32  haselbac
! Cosmetics
!
! Revision 1.11  2005/12/14 21:20:28  fnajjar
! Added subroutine and made changes for dynamic allocation of iPclsSend
!
! Revision 1.10  2005/12/13 23:30:53  haselbac
! Cosmetics
!
! Revision 1.9  2005/12/13 23:06:56  fnajjar
! Added defs of tag for PLAG
!
! Revision 1.8  2005/12/08 03:01:01  haselbac
! Major bug fix: spec tag was not set
!
! Revision 1.7  2005/12/03 19:48:23  haselbac
! Bug fix: Only clear request if have indeed sent message, cosmetics
!
! Revision 1.6  2005/09/19 18:40:37  haselbac
! Added IFs for border sizes before send and recv
!
! Revision 1.5  2005/07/08 15:01:29  haselbac
! Added profiling calls
!
! Revision 1.4  2005/05/26 22:01:01  haselbac
! Fixed two serious bugs: status is array in recv and wait
!
! Revision 1.3  2005/05/18 22:12:04  fnajjar
! ACH: Added routines to create and destroy iPclSend buffers
!
! Revision 1.2  2005/04/29 00:05:47  haselbac
! Added MPI_Wait routines
!
! Revision 1.1  2005/04/15 15:06:43  haselbac
! Initial revision
!
! ******************************************************************************























