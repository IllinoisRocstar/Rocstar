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
! Purpose: Allocate mixture memory for vertex variables.
!
! Description: None.
!
! Input: 
!   pRegion     Region pointer
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_AllocateMemoryVert.F90,v 1.8 2008/12/06 08:45:05 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_AllocateMemoryVert(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModGrid, ONLY: t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModMixture, ONLY: t_mixt_input
  USE ModMPI
      
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Parameters
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_mixt_input), POINTER :: pMixtInput
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: RFLU_AllocateMemoryVert.F90,v $ $Revision: 1.8 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_AllocateMemoryVert', &
                        'RFLU_AllocateMemoryVert.F90')

! ******************************************************************************
! Set pointers 
! ******************************************************************************

  pGrid      => pRegion%grid
  pMixtInput => pRegion%mixtInput

! ******************************************************************************
! Allocate memory
! ******************************************************************************

! ==============================================================================  
! Conserved variables
! ==============================================================================  
  
  ALLOCATE(pRegion%mixt%cvVert(CV_MIXT_DENS:CV_MIXT_ENER,pGrid%nVertTot), & 
           STAT=errorFlag)
  global%error = errorFlag           
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mixt%cvVert')
  END IF ! global%error  

! ==============================================================================  
! Dependent variables  
! ==============================================================================  
 
  ALLOCATE(pRegion%mixt%dvVert(pMixtInput%nDv,pGrid%nVertTot),STAT=errorFlag)
  global%error = errorFlag           
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mixt%dvVert')
  END IF ! global%error  
    
! ==============================================================================  
! Gas variables 
! ==============================================================================  

  IF ( pMixtInput%nGvAct == 0 ) THEN
    ALLOCATE(pRegion%mixt%gvVert(pMixtInput%nGv,0:1),STAT=errorFlag)
  ELSE
    ALLOCATE(pRegion%mixt%gvVert(pMixtInput%nGv,pGrid%nVertTot),STAT=errorFlag)
  END IF ! pMixtInput%nGvAct
  global%error = errorFlag           
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pRegion%mixt%gvVert')
  END IF ! global%error   

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_AllocateMemoryVert

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_AllocateMemoryVert.F90,v $
! Revision 1.8  2008/12/06 08:45:05  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:18:16  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.5  2006/03/13 15:39:21  haselbac
! Bug fix: Incorrect allocation of gvVert
!
! Revision 1.4  2005/11/14 17:04:43  haselbac
! Generalized to support pseudo-gas model
!
! Revision 1.3  2005/11/10 02:46:22  haselbac
! Cosmetics only
!
! Revision 1.2  2005/10/31 21:09:39  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/02/26 21:01:21  haselbac
! Initial revision
!
! ******************************************************************************







