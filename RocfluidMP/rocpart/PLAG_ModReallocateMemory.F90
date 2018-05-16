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
!*******************************************************************************
!
! Purpose: Suite of routines to dynamically reallocate memory.
!
! Description: None.
!
! Notes: None.
!
!*******************************************************************************
!
! $Id: PLAG_ModReallocateMemory.F90,v 1.9 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!*******************************************************************************

MODULE PLAG_ModReallocateMemory

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal,     ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModPartLag,    ONLY: t_plag
  USE ModMPI 
  USE PLAG_ModParameters
  USE ModInterfacesLagrangian, ONLY: PLAG_INRT_AllocMemTStep,   &
                                     PLAG_INRT_DeallocMemTStep, &  
                                     PLAG_RFLU_AllocMemSol,     &
                                     PLAG_RFLU_AllocMemTStep,   &                                    
                                     PLAG_RFLU_DeallocMemSol,   &
                                     PLAG_RFLU_DeallocMemTStep
 
  USE PLAG_ModDimensions,      ONLY: PLAG_SetMaxDimensions

  IMPLICIT NONE

  REAL(RFREAL), PARAMETER :: PLAG_EXPAND_RATIO =   0.90_RFREAL
  REAL(RFREAL), PARAMETER :: PLAG_SHRINK_RATIO =   0.25_RFREAL

  PRIVATE
  PUBLIC :: PLAG_ReallocMemWrapper

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: PLAG_ModReallocateMemory.F90,v $ $Revision: 1.9 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS




! ******************************************************************************
!
! Purpose: Copy dimensions for Lagrangian particle solution.
!
! Description: None.
!
! Input:
!   global       Global pointer
!   pPlag        Plag pointer
!   pPlagCopy    Plag pointer to copy
!
! Output: None.
!
! Notes: Memory dimensions operation is performed at the last RK-stage following
!        particle relocation, communication and injection.
!
! ******************************************************************************

SUBROUTINE PLAG_CopyDimensions(global,pPlag,pPlagCopy)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_plag),   POINTER :: pPlag,pPlagCopy

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ModReallocateMemory.F90,v $ $Revision: 1.9 $'

  CALL RegisterFunction(global,'PLAG_CopyDimensions',&
  'PLAG_ModReallocateMemory.F90')

! *****************************************************************************
! Copy particle dimensions from original datastructure to temporary datastructure
! *****************************************************************************

  pPlagCopy%nAiv = pPlag%nAiv
  pPlagCopy%nArv = pPlag%nArv
  pPlagCopy%nCv  = pPlag%nCv
  pPlagCopy%nDv  = pPlag%nDv
  pPlagCopy%nTv  = pPlag%nTv
  pPlagCopy%nPcls = pPlag%nPcls

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_CopyDimensions



! ******************************************************************************
!
! Purpose: Copy memory for Lagrangian particle solution.
!
! Description: None.
!
! Input:
!   global       Global pointer
!   pPlag        Plag pointer
!   pPlagCopy    Plag pointer to copy
!
! Output: None.
!
! Notes: 
!   1.  Memory copy operation is performed at the last RK-stage following
!       particle relocation, communication and injection.
!   2.  Copying of dv and tv data is not required given the memory swap
!       occurs at the last RK stage; but it is kept for consistency.
!
! ******************************************************************************

SUBROUTINE PLAG_CopyMemory(global,pPlag,pPlagCopy)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global
  TYPE(t_plag),   POINTER :: pPlag,pPlagCopy

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iPcl,iVar,nAiv,nArv,nCont,nCv,nDv,nPcls,nTv

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ModReallocateMemory.F90,v $ $Revision: 1.9 $'

  CALL RegisterFunction(global,'PLAG_CopyMemory',&
  'PLAG_ModReallocateMemory.F90')

! ******************************************************************************
! Set variables 
! ******************************************************************************

  nPcls = pPlag%nPcls

  nAiv = pPlag%nAiv
  nArv = pPlag%nArv

  nCv  = pPlag%nCv
  nDv  = pPlag%nDv
  nTv  = pPlag%nTv 

! *****************************************************************************
! Copy particle data from original datastructure to temporary datastructure
! *****************************************************************************

  DO iPcl = 1,nPcls
    DO iVar = 1,nCv
      pPlagCopy%cv(iVar,iPcl)     = pPlag%cv(iVar,iPcl)
      pPlagCopy%cvOld(iVar,iPcl)  = pPlag%cvOld(iVar,iPcl)
      pPlagCopy%rhsSum(iVar,iPcl) = pPlag%rhsSum(iVar,iPcl)
    END DO ! iVar
    
    DO iVar = 1,nDv
      pPlagCopy%dv(iVar,iPcl) = pPlag%dv(iVar,iPcl)
    END DO ! iVar
   
    DO iVar = 1,nTv
      pPlagCopy%tv(iVar,iPcl) = pPlag%tv(iVar,iPcl)
    END DO ! iVar
    
    DO iVar = 1,nAiv
      pPlagCopy%aiv(iVar,iPcl)    = pPlag%aiv(iVar,iPcl)
      pPlagCopy%aivOld(iVar,iPcl) = pPlag%aivOld(iVar,iPcl)
    END DO ! iVar
    
    DO iVar = 1,nArv   
      pPlagCopy%arv(iVar,iPcl)    = pPlag%arv(iVar,iPcl)
      pPlagCopy%arvOld(iVar,iPcl) = pPlag%arvOld(iVar,iPcl)
    END DO ! iVar
  END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_CopyMemory








! ******************************************************************************
!
! Purpose: Decide to reallocate memory for particle datastructure.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
! 
! Output: 
!  pDecideReallocMem Logical flag for decision
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE PLAG_DecideReallocMem(pRegion,pDecideReallocMem)

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion
  LOGICAL, INTENT(INOUT)  :: pDecideReallocMem

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: errorFlag,nPcls,nPclsMax
 
  REAL(RFREAL) :: pclRatio

  TYPE(t_global), POINTER :: global 

! ******************************************************************************
! Start
! ******************************************************************************
 
  RCSIdentString = '$RCSfile: PLAG_ModReallocateMemory.F90,v $ $Revision: 1.9 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_DecideReallocMem',&
  'PLAG_ModReallocateMemory.F90')

! ******************************************************************************
! Set variables 
! ******************************************************************************

  nPclsMax    = pRegion%plag%nPclsMax
  nPcls       = pRegion%plag%nPcls

  pDecideReallocMem = .FALSE.

! ******************************************************************************
! Set logical flag depending on ratio of nPcls to maximum size of datastructure
! ******************************************************************************

  pclRatio = REAL(nPcls,KIND=RFREAL)/REAL(nPclsMax,KIND=RFREAL)

  IF ( pclRatio >= PLAG_EXPAND_RATIO   .OR.  &
       ( pclRatio <= PLAG_SHRINK_RATIO .AND. & 
         nPclsMax > NPCLS_TOT_MIN )          ) pDecideReallocMem = .TRUE.

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_DecideReallocMem







! ******************************************************************************
!
! Purpose: Realllocate memory for Lagrangian particle datastructure.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE PLAG_ReallocMem(pRegion)

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

  INTEGER :: errorFlag,nPcls,nPlagPclsTot

  TYPE(t_global), POINTER :: global 
  TYPE(t_plag),   POINTER :: pPlag,pPlagTemp

! ******************************************************************************
! Start
! ******************************************************************************
 
  RCSIdentString = '$RCSfile: PLAG_ModReallocateMemory.F90,v $ $Revision: 1.9 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_ReallocMem',&
  'PLAG_ModReallocateMemory.F90')

! ******************************************************************************
! Set pointers and variables 
! ******************************************************************************

  pPlag => pRegion%plag
  pPlagTemp => pRegion%plagTemp

! ******************************************************************************
! Copy particle dimensions from original datastructure 
!   to temporary datastructure
! ******************************************************************************

  CALL PLAG_CopyDimensions( global,pPlag,pPlagTemp )

! ******************************************************************************
! Allocate memory for intermediate datastructure
! ******************************************************************************

  CALL PLAG_RFLU_AllocMemSol( pRegion,pPlagTemp )
  CALL PLAG_RFLU_AllocMemTStep( pRegion,pPlagTemp )

! *****************************************************************************
! Copy particle data from original datastructure to temporary datastructure
! *****************************************************************************

  CALL PLAG_CopyMemory( global,pPlag,pPlagTemp )

! ******************************************************************************
! Destroy memory of original datastructure
! ******************************************************************************

  CALL PLAG_RFLU_DeallocMemSol( pRegion,pPlag )
  CALL PLAG_RFLU_DeallocMemTStep( pRegion,pPlag )
  CALL PLAG_INRT_DeallocMemTStep( pRegion,pPlag )

! ******************************************************************************
! Allocate memory for original datastructure with updated maximum value
! ******************************************************************************

  CALL PLAG_RFLU_AllocMemSol( pRegion,pPlag )
  CALL PLAG_RFLU_AllocMemTStep( pRegion,pPlag )
  CALL PLAG_INRT_AllocMemTStep( pRegion,pPlag )

! *****************************************************************************
! Copy particle data from temporary datastructure to original datastructre
! *****************************************************************************

  CALL PLAG_CopyMemory( global,pPlagTemp,pPlag )

! ******************************************************************************
! Destroy memory of temporary datastructure
! ******************************************************************************

  CALL PLAG_RFLU_DeallocMemSol( pRegion,pPlagTemp )
  CALL PLAG_RFLU_DeallocMemTStep( pRegion,pPlagTemp )

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_ReallocMem







! ******************************************************************************
!
! Purpose: Wrapper routine to realllocate memory for Lagrangian 
!particle datastructure.
!
! Description: None.
!
! Input:
!   pRegion        Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE PLAG_ReallocMemWrapper(pRegion)

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

  LOGICAL :: pDecideReallocMem

  TYPE(t_global), POINTER :: global 
  TYPE(t_plag),   POINTER :: pPlag

! ******************************************************************************
! Start
! ******************************************************************************
 
  RCSIdentString = '$RCSfile: PLAG_ModReallocateMemory.F90,v $ $Revision: 1.9 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_ReallocMem',&
  'PLAG_ModReallocateMemory.F90')

! ******************************************************************************
! Set pointers  
! ******************************************************************************

  pPlag => pRegion%plag

! ******************************************************************************
! Driver for memory reallocation
! ******************************************************************************

  CALL PLAG_DecideReallocMem( pRegion, pDecideReallocMem )

  IF ( pDecideReallocMem .EQV. .TRUE. ) THEN
    CALL PLAG_SetMaxDimensions( pRegion )
    CALL PLAG_ReallocMem( pRegion )

! TEMPORARY
    WRITE(STDOUT,'(A,I2,2X,I10,2X,I10)') ' PLAG_ReallocMem-iReg: nPcls, nPclsMax = ',&
    pRegion%iRegionGlobal,pRegion%plag%nPcls,pRegion%plag%nPclsMax
! END TEMPORARY

  END IF ! pDecideReallocMem

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_ReallocMemWrapper







! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE PLAG_ModReallocateMemory

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModReallocateMemory.F90,v $
! Revision 1.9  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2007/03/20 17:37:25  fnajjar
! Moved USE call for PLAG_SetMaxDimensions to new module PLAG_ModDimensions
!
! Revision 1.6  2007/03/06 23:15:32  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.5  2006/05/02 17:47:19  fnajjar
! Increased integer size in format statement
!
! Revision 1.4  2005/12/13 23:09:09  fnajjar
! Added memory allocation, copying and deallocation of cvOld,rhsSum,aivOld and arvOld to fix bug for parallel PLAG w RFLU
!
! Revision 1.3  2004/07/29 20:03:56  fnajjar
! Bug fix in PLAG_DecideReallocMem
!
! Revision 1.2  2004/07/29 16:32:50  fnajjar
! Included temporary io statement for monitoring
!
! Revision 1.1  2004/07/28 18:59:50  fnajjar
! Initial import for dynamic memory reallocation
!
!******************************************************************************











