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
! Purpose: Read in user input related to mixture.
!
! Description: None.
!
! Input: 
!   regions     Data for regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReadMixtureSection.F90,v 1.5 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReadMixtureSection(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  USE ModDataStruct, ONLY: t_region  
  
  USE ModInterfaces, ONLY: ReadSection  
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Local variables
! ==============================================================================

  INTEGER :: iReg,nVals

#ifdef RFLU
  INTEGER, PARAMETER :: NVALS_MAX = 2
#endif

  CHARACTER(10) :: keys(NVALS_MAX)
  LOGICAL :: defined(NVALS_MAX)
  REAL(RFREAL) :: vals(NVALS_MAX)
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction(global,'ReadMixtureSection',&
  'ReadMixtureSection.F90')

! specify keywords and search for them

  nVals = NVALS_MAX

#ifdef RFLU
  keys(1) = 'FROZENFLAG'
  keys(2) = 'GASMODEL'     
#endif

#ifdef RFLU
  CALL ReadSection(global,IF_INPUT,nVals,keys,vals,defined)

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)  
    IF ( defined(1) .EQV. .TRUE. ) THEN 
      IF ( NINT(vals(1)) == 1 ) THEN 
        regions(iReg)%mixtInput%frozenFlag = .TRUE.
      ELSE 
        regions(iReg)%mixtInput%frozenFlag = .FALSE.
      END IF ! NINT(vals(1))
    END IF ! defined               
  END DO ! iReg

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)  
    IF ( defined(2) .EQV. .TRUE. ) THEN 
      regions(iReg)%mixtInput%gasModel = NINT(vals(2))
    END IF ! defined               
  END DO ! iReg
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE ReadMixtureSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadMixtureSection.F90,v $
! Revision 1.5  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.2  2005/11/10 22:20:06  fnajjar
! ACH: Added frozenFlag
!
! Revision 1.1  2005/10/31 19:23:53  haselbac
! Initial revision
!
! ******************************************************************************







