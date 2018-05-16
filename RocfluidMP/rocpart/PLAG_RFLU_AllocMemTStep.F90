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
! Purpose: Allocate memory for Lagrangian particles related to time stepping.
!
! Description: None.
!
! Input: 
!   pRegion     Region pointer
!   pPlag       Pointer to particle data structure
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_RFLU_AllocMemTStep.F90,v 1.8 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_AllocMemTStep(pRegion,pPlag)

  USE ModDataTypes
  USE ModError
  USE ModParameters  
  USE ModGlobal, ONLY: t_global 
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModPartLag, ONLY: t_plag  
  USE ModMPI
   
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion
  TYPE(t_plag), POINTER :: pPlag

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,iPcl,iVar,nAiv,nArv,nCv,nPclsMax
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_AllocMemTStep.F90,v $ $Revision: 1.8 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_AllocMemTStep',&
  'PLAG_RFLU_AllocMemTStep.F90')

! ******************************************************************************
! Set pointers and variables 
! ******************************************************************************

  pGrid  => pRegion%grid

  nPclsMax = pRegion%plag%nPclsMax

  nAiv = pPlag%nAiv
  nArv = pPlag%nArv

  nCv  = pPlag%nCv

! ******************************************************************************
! Allocate memory
! ******************************************************************************

! ==============================================================================
! Old state vector
! ==============================================================================

  ALLOCATE(pPlag%cvOld(nCv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%cvOld')
  END IF ! global%error 

! ==============================================================================
! Additional variables  
! ==============================================================================

  ALLOCATE(pPlag%aivOld(nAiv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%aivOld')
  END IF ! global%error 

  ALLOCATE(pPlag%arvOld(nArv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%arvOld')
  END IF ! global%error

! ==============================================================================     
! Residuals 
! ==============================================================================
  
  ALLOCATE(pPlag%rhs(nCv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%rhs')
  END IF ! global%error  
  
  ALLOCATE(pPlag%rhsSum(nCv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%rhsSum')
  END IF ! global%error  

! ******************************************************************************
! Initialize memory
! ******************************************************************************

  DO iPcl = 1,nPclsMax
    DO iVar = 1,nCv
      pPlag%cvOld(iVar,iPcl)  = 0.0_RFREAL
      pPlag%rhs(iVar,iPcl)    = 0.0_RFREAL
      pPlag%rhsSum(iVar,iPcl) = 0.0_RFREAL
    END DO ! iVar

    DO iVar = 1,nAiv
      pPlag%aivOld(iVar,iPcl) = 0
    END DO ! iVar
    
    DO iVar = 1,nArv   
      pPlag%arvOld(iVar,iPcl) = 0
    END DO ! iVar
  END DO ! iPcl
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_AllocMemTStep

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_AllocMemTStep.F90,v $
! Revision 1.8  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2007/03/06 23:15:32  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.5  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.4  2004/07/28 18:58:13  fnajjar
! Included new definition for nPclsTot from dynamic memory reallocation
!
! Revision 1.3  2004/07/16 20:09:39  fnajjar
! Bug fix to initialize cvOld instead of cv
!
! Revision 1.2  2004/03/08 23:02:28  fnajjar
! Added initialization section within DO-loop construct
!
! Revision 1.1  2004/02/26 21:00:39  haselbac
! Initial revision
!
!******************************************************************************







