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
! Purpose: Allocate memory for Lagrangian particle solution.
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
! $Id: PLAG_RFLU_AllocMemSol.F90,v 1.7 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_AllocMemSol(pRegion,pPlag)

  USE ModDataTypes
  USE ModError
  USE ModParameters  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModPartLag, ONLY: t_plag
  USE ModMPI
   
  USE PLAG_ModParameters  
   
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
  INTEGER :: errorFlag,iCont,iPcl,iVar,nAiv,nArv,nCont,nCv,nDv,nPclsMax,nTv
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_AllocMemSol.F90,v $ $Revision: 1.7 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_AllocMemSol',&
  'PLAG_RFLU_AllocMemSol.F90')

! ******************************************************************************
! Set pointers and variables 
! ******************************************************************************

  pGrid => pRegion%grid

  nPclsMax = pRegion%plag%nPclsMax
  nCont    = pRegion%plagInput%nCont

  nAiv = pPlag%nAiv
  nArv = pPlag%nArv

  nCv  = pPlag%nCv
  nDv  = pPlag%nDv
  nTv  = pPlag%nTv

! ******************************************************************************
! Allocate memory
! ******************************************************************************

! ==============================================================================
! State vector
! ==============================================================================

  ALLOCATE(pPlag%cv(nCv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%cv')
  END IF ! global%error 

! ==============================================================================
! Dependent variables  
! ==============================================================================
    
  ALLOCATE(pPlag%dv(nDv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%dv')
  END IF ! global%error 

! ==============================================================================
! Transport variables  
! ==============================================================================
    
  ALLOCATE(pPlag%tv(nTv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%tv')
  END IF ! global%error 

! ==============================================================================
! Additional variables  
! ==============================================================================

  ALLOCATE(pPlag%aiv(nAiv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%aiv')
  END IF ! global%error 

  ALLOCATE(pPlag%arv(nArv,nPclsMax),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%arv')
  END IF ! global%error

! ==============================================================================
! Lagrangian particles mass and volume indices
! ==============================================================================

  ALLOCATE(pPlag%cvPlagMass(nCont),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global, ERR_ALLOCATE,__LINE__ ,'pPlag%cvPlagMass')
  END IF ! global%error

  ALLOCATE(pPlag%dvPlagVolu(nCont),STAT=errorFlag)
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop(global, ERR_ALLOCATE,__LINE__ ,'pPlag%dvPlagVolu')
  END IF ! global%error

  DO iCont = 1,nCont
    pPlag%cvPlagMass(iCont) = CV_PLAG_LAST + iCont
    pPlag%dvPlagVolu(iCont) = DV_PLAG_LAST + iCont
  END DO ! iCont

! ******************************************************************************
! Initialize memory
! ******************************************************************************

  DO iPcl = 1,nPclsMax
    DO iVar = 1,nCv
      pPlag%cv(iVar,iPcl) = 0.0_RFREAL
    END DO ! iVar
    
    DO iVar = 1,nDv
      pPlag%dv(iVar,iPcl) = 0.0_RFREAL
    END DO ! iVar
   
    DO iVar = 1,nTv
      pPlag%tv(iVar,iPcl) = 0.0_RFREAL
    END DO ! iVar
    
    DO iVar = 1,nAiv
      pPlag%aiv(iVar,iPcl) = 0
    END DO ! iVar
    
    DO iVar = 1,nArv   
      pPlag%arv(iVar,iPcl) = 0
    END DO ! iVar
  END DO ! iPcl

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_AllocMemSol

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_AllocMemSol.F90,v $
! Revision 1.7  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2007/03/06 23:15:32  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.4  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.3  2004/07/28 18:58:13  fnajjar
! Included new definition for nPclsTot from dynamic memory reallocation
!
! Revision 1.2  2004/03/08 23:02:28  fnajjar
! Added initialization section within DO-loop construct
!
! Revision 1.1  2004/02/26 21:00:36  haselbac
! Initial revision
!
!******************************************************************************







