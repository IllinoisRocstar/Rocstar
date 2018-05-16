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
!******************************************************************************
!
! Purpose: Compute integral 2 of optimal LES approach.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to current region
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_ComputeIntegral2OLES.F90,v 1.7 2008/12/06 08:44:29 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ComputeIntegral2OLES(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModMPI

  USE ModInterfaces, ONLY: RFLU_ComputeDCUHREInfo

  USE RFLU_ModOLES

  IMPLICIT NONE

! --- parameters      

  TYPE(t_region), POINTER :: pRegion

! --- local variables

  INTEGER :: c1g,c2g,errorFlag,hLoc,iInt,ic1l,ic2l,ifc,ifcp,iFun,iFunNZ,j, & 
             key,l,loopCounter,m,n,nCells,nFun,nFunNZ,restartFlag,vLoc
  REAL(RFREAL) :: normFact,normFactTerm
  
  TYPE(t_grid), POINTER :: pGrid    
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Start
! ==============================================================================

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeIntegral2OLES',&
  'RFLU_ComputeIntegral2OLES.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Computing integral 2...'                                              
  END IF ! global%verbLevel

! ==============================================================================
! Set grid pointer
! ==============================================================================

  pGrid => pRegion%grid

! ==============================================================================
! Set various quantities
! ==============================================================================

! ------------------------------------------------------------------------------
! Set normalization factor term, modified inside face loop 
! ------------------------------------------------------------------------------

  normFactTerm = 1.0_RFREAL/pGrid%deltaOLES**6

! ------------------------------------------------------------------------------
! Set DCUHRE parameters
! ------------------------------------------------------------------------------

  nDim   = 6
  nFun   = 9
  nFunNZ = 3
  key    = 0 ! Use accurate integration

  errAbsReq = 10.0_RFREAL*EPSILON(1.0_RFREAL)
  errRelReq = 2.0E-3_RFREAL 

! ------------------------------------------------------------------------------
! Compute DCUHRE information and allocate arrays 
! ------------------------------------------------------------------------------

  maxCalls = MAX_CALLS_LIMIT

  CALL RFLU_ComputeDCUHREInfo(global,nDim,nFun,key,maxCalls,workArraySize)
  CALL RFLU_AllocateDCUHREArrays(global,nDim,nFunNZ,nFun)

  ALLOCATE(workArray(workArraySize),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'workArray')
  END IF ! global%error

  ALLOCATE(integral(nFun),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'integral')
  END IF ! global%error

  integral(:)   = 0.0_RFREAL
  integralNZ(:) = 0.0_RFREAL

! ------------------------------------------------------------------------------
! Set non-zero functions 
! ------------------------------------------------------------------------------

  CALL RFLU_SetMapFunNZ2FunCorr22(nFunNZ)

! ==============================================================================
! Loop over protoype faces
! ==============================================================================

  DO ifcp = 1,3 
    ifc = pGrid%fp2fOLES(ifcp)

    nCells = SIZE(pGrid%fsOLES,1)

! ------------------------------------------------------------------------------
!   Loop over cells
! ------------------------------------------------------------------------------        

    DO ic1l = 1,nCells
      c1g = pGrid%fsOLES(ic1l,ifc)

      lowLim(1) = pGrid%intLimOLES(INT_LIM_LOW,XCOORD,c1g) 
      lowLim(2) = pGrid%intLimOLES(INT_LIM_LOW,YCOORD,c1g)
      lowLim(3) = pGrid%intLimOLES(INT_LIM_LOW,ZCOORD,c1g)                

      uppLim(1) = pGrid%intLimOLES(INT_LIM_UPP,XCOORD,c1g)
      uppLim(2) = pGrid%intLimOLES(INT_LIM_UPP,YCOORD,c1g)
      uppLim(3) = pGrid%intLimOLES(INT_LIM_UPP,ZCOORD,c1g)              

      DO ic2l = 1,nCells
        c2g = pGrid%fsOLES(ic2l,ifc)

        lowLim(4) = pGrid%intLimOLES(INT_LIM_LOW,XCOORD,c2g)
        lowLim(5) = pGrid%intLimOLES(INT_LIM_LOW,YCOORD,c2g)
        lowLim(6) = pGrid%intLimOLES(INT_LIM_LOW,ZCOORD,c2g)              

        uppLim(4) = pGrid%intLimOLES(INT_LIM_UPP,XCOORD,c2g)
        uppLim(5) = pGrid%intLimOLES(INT_LIM_UPP,YCOORD,c2g)
        uppLim(6) = pGrid%intLimOLES(INT_LIM_UPP,ZCOORD,c2g)            

! ----- Loop over integrals --------------------------------------------------

        DO iInt = 1,2        
          errorFlag   = 0
          maxCalls    = MAX_CALLS_START
          restartFlag = 0 

          DO loopCounter = 1,DCUHRE_LOOP_LIMIT           
            IF ( iInt == 1 ) THEN 
              CALL DCUHRE(nDim,nFunNZ,lowLim,uppLim,MIN_CALLS,maxCalls, & 
                          RFLU_DefineCorrelation220,errAbsReq,errRelReq, & 
                          key,workArraySize,restartFlag,integralNZ, & 
                          errAbsEst,nEval,errorFlag,workArray)        
            ELSE 
              CALL DCUHRE(nDim,nFunNZ,lowLim,uppLim,MIN_CALLS,maxCalls, & 
                          RFLU_DefineCorrelation221,errAbsReq,errRelReq, & 
                          key,workArraySize,restartFlag,integralNZ, & 
                          errAbsEst,nEval,errorFlag,workArray)
            END IF ! iInt

            IF ( errorFlag == 0 .OR. loopCounter == DCUHRE_LOOP_LIMIT ) THEN
              EXIT
            ELSE IF ( errorFlag == 1 ) THEN 
              restartFlag = 1
              maxCalls    = MAX_CALLS_FACTOR*maxCalls

              CALL RFLU_ComputeDCUHREInfo(global,nDim,nFunNZ,key,maxCalls, & 
                                          workArraySizeNew)

              IF ( workArraySizeNew > workArraySize ) THEN 
                EXIT       
              END IF ! workArraySizeNew
            ELSE 
              CALL ErrorStop(global,ERR_DCUHRE_OUTPUT,__LINE__)                            
            END IF ! errorFlag
          END DO ! loopCounter

! ------- Normalize integral ---------------------------------------------------

          IF ( iInt == 1 ) THEN 
            normFact = normFactTerm
          ELSE 
            normFact = normFactTerm*CONST_KOLMOGOROV/ & 
                       (6.0_RFREAL*pGrid%rhoOLES(ifc)**(2.0_RFREAL/3.0_RFREAL))     
          END IF ! iInt

          DO iFunNZ = 1,nFunNZ 
            iFun = mapFunNZ2FunCorr22(iFunNZ)                       
            integral(iFun) = normFact*integralNZ(iFunNZ)
          END DO ! iFunNZ

! ------- Store integral in array ----------------------------------------------
 
          DO iFun = 1,nFun              
            CALL RFLU_MapK2IJ(iFun,l,j)  
             
            vLoc = RFLU_GetI1PosOLES(l,ic1l)
            hLoc = RFLU_GetLPosOLES(j,ic2l)
 
            IF ( iInt == 1 ) THEN
              pGrid%int20OLES(ifcp,vLoc,hLoc) = integral(iFun) 
            ELSE
              pGrid%int21OLES(ifcp,vLoc,hLoc) = integral(iFun)             
            END IF ! iInt

          END DO ! iFun     
        END DO ! iInt

      END DO ! icl2
    END DO ! icl1                 
  END DO ! ifc

#ifdef CHECK_DATASTRUCT
! --- Data structure output for checking
      WRITE(STDOUT,'(A)') SOLVER_NAME
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Optimal LES I20 integral matrix'
      DO ifcp = 1,3 ! loop over prototype faces
        WRITE(STDOUT,'(A,3X,A,1X,I2)') SOLVER_NAME,'Face:',ifcp
        DO vLoc = 1,3*nCells ! loop over components
          WRITE(STDOUT,'(A,1X,I6,6(1X,E11.4))') SOLVER_NAME,vLoc, & 
                                           pGrid%int20OLES(ifcp,vLoc,1:3*nCells)
        END DO ! vLoc
      END DO ! ifcp    
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Optimal LES I21 integral matrix'
      DO ifcp = 1,3 ! loop over prototype faces
        WRITE(STDOUT,'(A,3X,A,1X,I2)') SOLVER_NAME,'Face:',ifcp
        DO vLoc = 1,3*nCells ! loop over components
          WRITE(STDOUT,'(A,1X,I6,6(1X,E11.4))') SOLVER_NAME,vLoc, & 
                                           pGrid%int21OLES(ifcp,vLoc,1:3*nCells)
        END DO ! vLoc
      END DO ! ifcp     
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
      WRITE(STDOUT,'(A)') SOLVER_NAME   
#endif 

! ==============================================================================
! Deallocate arrays
! ==============================================================================

  CALL RFLU_DeallocateDCUHREArrays(global)

  DEALLOCATE(integral,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'integral')
  END IF ! global%error

  DEALLOCATE(workArray,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'workArray')
  END IF ! global%error      

! ==============================================================================
! End
! ============================================================================== 

  CALL DeregisterFunction(global)     

END SUBROUTINE RFLU_ComputeIntegral2OLES

!*******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeIntegral2OLES.F90,v $
! Revision 1.7  2008/12/06 08:44:29  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2003/03/15 18:30:58  haselbac
! Added footer
!
!*******************************************************************************







