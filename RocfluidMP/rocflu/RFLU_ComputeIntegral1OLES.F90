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
! Purpose: Compute integral 1 of optimal LES approach.
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
! $Id: RFLU_ComputeIntegral1OLES.F90,v 1.8 2008/12/06 08:44:29 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ComputeIntegral1OLES(pRegion)

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

  INTEGER :: c1g,errorFlag,i,ic1l,ifc,ifcp,iFun,iFunNZ,key,l, & 
             loopCounter,m,n,nCells,nFun,nFunNZ,restartFlag,vLoc
  REAL(RFREAL) :: normFact,normFactTerm
  
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Start
! ==============================================================================

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_ComputeIntegral1OLES',&
  'RFLU_ComputeIntegral1OLES.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Computing integral 1...'        
  END IF ! global%verbLevel

! ==============================================================================
! Set grid pointer
! ==============================================================================

  pGrid => pRegion%grid

! ==============================================================================
! Set DCUHRE parameters
! ==============================================================================

  nDim   = 5
  nFun   = 9
  nFunNZ = 3 ! Change this to nFun if you want all integrals evaluated
  key    = 0 ! Use accurate integration

  errAbsReq = 10.0_RFREAL*EPSILON(1.0_RFREAL)
  errRelReq = 5.0E-3_RFREAL 

! ==============================================================================
! Compute DCUHRE information and allocate arrays
! ==============================================================================

  maxCalls = MAX_CALLS_LIMIT

  CALL RFLU_ComputeDCUHREInfo(global,nDim,nFunNZ,key,maxCalls,workArraySize)
  CALL RFLU_AllocateDCUHREArrays(global,nDim,nFunNZ,nFun)

  ALLOCATE(integral(nFun),STAT=errorFlag)
  global%error = errorFlag  
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'integral')
  END IF ! global%error

  integral(:) = 0.0_RFREAL

  ALLOCATE(workArray(workArraySize),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'workArray')
  END IF ! global%error

! ==============================================================================
! Set normalization factor
! ==============================================================================

  normFact = 1.0_RFREAL/(15.0_RFREAL*pGrid%deltaOLES**6)

! ==============================================================================
! Set non-zero functions
! ==============================================================================

  CALL RFLU_SetMapFunNZ2FunCorr32(nFunNZ)

! ==============================================================================
! Loop over prototype faces
! ==============================================================================

  DO ifcp = 1,3 
    ifc = pGrid%fp2fOLES(ifcp)

! ------------------------------------------------------------------------------
!   Get face limits from normal face list
! ------------------------------------------------------------------------------

    c1g = pGrid%f2c(1,ifc)

    dummy = MAXLOC(ABS(pGrid%fn(1:3,ifc)))    
    nzLoc = dummy(1) ! Used in RFLU_DefineCorrelation32 (next two also)     
    nzVal = pGrid%fc(nzLoc,ifc)       
    nzSgn = NINT(pGrid%fn(nzLoc,ifc)) 

    n = nDim-2

    DO m = 1,3
      IF ( m /= nzLoc ) THEN
        n = n + 1

        lowLim(n) = pGrid%intLimOLES(INT_LIM_LOW,m,c1g) 
        uppLim(n) = pGrid%intLimOLES(INT_LIM_UPP,m,c1g)
      END IF ! m 
    END DO ! m

! ------------------------------------------------------------------------------
!   Loop over cells
! ------------------------------------------------------------------------------        

    nCells = SIZE(pGrid%fsOLES,1)

    DO ic1l = 1,nCells
      c1g = pGrid%fsOLES(ic1l,ifc)

      lowLim(1) = pGrid%intLimOLES(INT_LIM_LOW,XCOORD,c1g) 
      lowLim(2) = pGrid%intLimOLES(INT_LIM_LOW,YCOORD,c1g)
      lowLim(3) = pGrid%intLimOLES(INT_LIM_LOW,ZCOORD,c1g)                

      uppLim(1) = pGrid%intLimOLES(INT_LIM_UPP,XCOORD,c1g)
      uppLim(2) = pGrid%intLimOLES(INT_LIM_UPP,YCOORD,c1g)
      uppLim(3) = pGrid%intLimOLES(INT_LIM_UPP,ZCOORD,c1g)        


      errorFlag   = 0
      maxCalls    = MAX_CALLS_START
      restartFlag = 0 

! --- Compute integrals --------------------------------------------------------

      DO loopCounter = 1,DCUHRE_LOOP_LIMIT 
        CALL DCUHRE(nDim,nFunNZ,lowLim,uppLim,MIN_CALLS,maxCalls, & 
                    RFLU_DefineCorrelation32,errAbsReq,errRelReq, & 
                    key,workArraySize,restartFlag,integralNZ, & 
                    errAbsEst,nEval,errorFlag,workArray)

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

! --- Normalize integral -------------------------------------------------------

      DO iFunNZ = 1,nFunNZ 
        iFun = mapFunNZ2FunCorr32(iFunNZ) 
        integral(iFun) = normFact*integralNZ(iFunNZ) 
      END DO ! iFunNZ

! --- Store integral in array --------------------------------------------------

      DO iFun = 1,nFun
        CALL RFLU_MapK2IJ(iFun,l,i)        
        vLoc = RFLU_GetI1PosOLES(l,ic1l)     
        pGrid%int1OLES(i,ifcp,vLoc) = integral(iFun)
      END DO ! iFun             
    END DO ! ic1l       
  END DO ! ifcp


#ifdef CHECK_DATASTRUCT
! --- Data structure output for checking
      WRITE(STDOUT,'(A)') SOLVER_NAME
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Optimal LES I1 integral vector'
      DO i = 1,3 ! loop over components
        WRITE(STDOUT,'(2(A,1X),I1)') SOLVER_NAME,'Component:',i
        DO ifcp = 1,3 ! loop over prototype faces
          WRITE(STDOUT,'(A,3X,A,1X,I2)') SOLVER_NAME,'Face:',ifcp
          DO vLoc = 1,3*nCells ! loop over components
            WRITE(STDOUT,'(A,1X,I6,1X,E11.4)') SOLVER_NAME,vLoc, & 
                                               pGrid%int1OLES(i,ifcp,vLoc)
          END DO ! vLoc
        END DO ! ifcp
      END DO ! i 
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

END SUBROUTINE RFLU_ComputeIntegral1OLES

!*******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeIntegral1OLES.F90,v $
! Revision 1.8  2008/12/06 08:44:29  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2003/05/16 02:27:44  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.5  2003/03/15 18:29:42  haselbac
! Added KIND qualifyer to NINT, added footer
!
!*******************************************************************************







