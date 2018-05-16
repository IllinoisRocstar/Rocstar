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
! Purpose: Update dependent variables.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data
!
! Output: None.
!
! Notes: 
!   1. This body of this routine may eventually be split up into separate 
!      parts like the routine MixtureProperties.
!
! ******************************************************************************
!
! $Id: SPEC_UpdateDependentVars.F90,v 1.6 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_UpdateDependentVars(pRegion)

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region

  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: icg,iSpec
  REAL(RFREAL) :: mum
  REAL(RFREAL), DIMENSION(:,:), POINTER :: tvm,tvs
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_UpdateDependentVars.F90,v $ $Revision: 1.6 $'

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_UpdateDependentVars',&
  'SPEC_UpdateDependentVars.F90')

! ******************************************************************************
! Update transport variables
! ******************************************************************************

  IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN 
    tvm => pRegion%mixt%tv
    tvs => pRegion%spec%tv

    DO icg = 1,pRegion%grid%nCellsTot
      mum = tvm(TV_MIXT_MUEL,icg)

      DO iSpec = 1,pRegion%specInput%nSpecies 
        tvs(iSpec,icg) = mum/pRegion%specInput%specType(iSpec)%schmidtNumber
      END DO ! iSpec
    END DO ! icg
  END IF ! pRegion%mixtInput%flowModel

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_UpdateDependentVars

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_UpdateDependentVars.F90,v $
! Revision 1.6  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/04/07 15:19:25  haselbac
! Removed tabs
!
! Revision 1.3  2004/12/22 00:40:18  haselbac
! Bug fix: Added missing POINTER statement
!
! Revision 1.2  2004/11/14 19:48:46  haselbac
! Changed interface, clean-up
!
! Revision 1.1  2004/01/29 22:59:32  haselbac
! Initial revision
!
! ******************************************************************************







