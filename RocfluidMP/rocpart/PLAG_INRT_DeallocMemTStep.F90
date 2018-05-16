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
! Purpose: Allocate memory for Lagrangian particles sources 
!          related to time stepping.
!
! Description: none.

! Input: 
!   pRegion     Region pointer
!   pPlag       Pointer to particle data structure
!
! Output: None.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_INRT_DeallocMemTStep.F90,v 1.4 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_INRT_DeallocMemTStep( pRegion,pPlag )

  USE ModDataTypes
  USE ModError
  USE ModParameters  
  USE ModGlobal, ONLY: t_global 
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
  TYPE(t_plag),   POINTER :: pPlag

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: maxDisEdges,errorFlag

  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_INRT_DeallocMemTStep.F90,v $ $Revision: 1.4 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_INRT_DeallocMemTStep',&
  'PLAG_INRT_DeallocMemTStep.F90')

! ******************************************************************************
! Set variables 
! ******************************************************************************

  maxDisEdges = pRegion%inrtInput%maxDisEdges

! ******************************************************************************
! Deallocate memory
! ******************************************************************************

  IF (maxDisEdges > 0) THEN
    DEALLOCATE( pPlag%inrtSources,stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= ERR_NONE) THEN 
      CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__,'pPlag%inrtSources')
    END IF ! global%error 

  ELSE
    NULLIFY( pPlag%inrtSources )

  END IF ! maxDisEdges
 
! ******************************************************************************
! End
! ******************************************************************************

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_INRT_DeallocMemTStep

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_INRT_DeallocMemTStep.F90,v $
! Revision 1.4  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.1  2004/07/26 19:01:20  fnajjar
! Initial import
!
!******************************************************************************







