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
! Purpose: Integrate source terms using operator-splitting method.
!
! Description: None.
!
! Input: 
!   regions             Data for all grid regions
!
! Output: None.
!
! Notes: 
!   1. This routine is work in progress... 
!
! ******************************************************************************
!
! $Id: IntegrateSourceTermsMP.F90,v 1.6 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE IntegrateSourceTermsMP(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters

#ifdef RFLU
#ifdef SPEC
  USE SPEC_RFLU_ModChemistry, ONLY: SPEC_RFLU_IntegrateChemSrcTerm
  USE ModInterfaces, ONLY: RFLU_SetVarsWrapper
#endif  
#endif

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

  TYPE(t_region), POINTER :: regions(:)
#ifdef RFLU
  TYPE(t_region), POINTER :: pRegion  
#endif

! ==============================================================================
! Arguments
! ==============================================================================

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
#ifdef RFLU
  INTEGER :: iReg
#endif  
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: IntegrateSourceTermsMP.F90,v $ $Revision: 1.6 $'

! ******************************************************************************
! Set pointers and variables
! ****************************************************************************** 

  global => regions(1)%global

  CALL RegisterFunction(global,'IntegrateSourceTermsMP',&
  'IntegrateSourceTermsMP.F90')

#ifdef RFLU
#ifdef SPEC  
! ******************************************************************************
! Integrate chemistry source term and update dependent variables
! ******************************************************************************  

! TEMPORARY - At present, do not use operator-split integration of chemistry
!             source terms based on modifications by Luca. This means that 
!             call source term residual function directly in SourceTermsMP.F90
!  DO iReg = 1,global%nRegionsLocal
!    pRegion => regions(iReg)
!      
!    IF ( pRegion%specInput%sourceFlag .EQV. .TRUE. ) THEN 
!      CALL SPEC_RFLU_IntegrateChemSrcTerm(pRegion)
!      CALL RFLU_SetVarsWrapper(pRegion,1,pRegion%grid%nCellsTot)
!    END IF ! pRegion%specInput%sourceType
!  END DO ! iReg
! END TEMPORARY
#endif
#endif

! ******************************************************************************
! End
! ****************************************************************************** 

  CALL DeregisterFunction(global)

END SUBROUTINE IntegrateSourceTermsMP

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: IntegrateSourceTermsMP.F90,v $
! Revision 1.6  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.3  2005/06/06 14:22:09  haselbac
! Adapted to Lucas changes
!
! Revision 1.2  2005/04/15 15:06:02  haselbac
! Adapted call to RFLU_SetVarsWrapper
!
! Revision 1.1  2004/12/01 16:48:43  haselbac
! Initial revision after changing case
!
! Revision 1.3  2004/11/14 20:30:23  haselbac
! Bug fix: Moved USE inside ifdef RFLU
!
! Revision 1.2  2004/11/14 19:35:06  haselbac
! Replaced call to UpdateDependentVarsMP by RFLU_SetVarsWrapper, cosmetics
!
! Revision 1.1  2004/04/01 21:22:17  haselbac
! Initial revision
!
! ******************************************************************************







