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
! Purpose: allocate memory for time-dependent boundary conditions
!
! Description: none.
!
! Input: none.
!
! Output: none.
!
! Notes: this cannot be called from RFLU_AllocateMemoryWrapper because it
!        requires that the BC input file has been read already
!
! ******************************************************************************
!
! $Id: RFLU_AllocateMemoryTbc.F90,v 1.11 2008/12/06 08:44:11 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_AllocateMemoryTbc(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModMixture, ONLY: t_mixt_input
  USE ModBndPatch, ONLY: t_patch,t_bcvalues
  USE ModError
  USE ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: errorFlag,fType,iReg,iPatch,iLev
  TYPE(t_global),     POINTER :: global
  TYPE(t_mixt_input), POINTER :: input
  TYPE(t_patch),      POINTER :: pPatch
  TYPE(t_bcvalues),   POINTER :: bc

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction( global,'RFLU_AllocateMemoryTbc',&
  'RFLU_AllocateMemoryTbc.F90' )

! ******************************************************************************
! Allocate TBCs
! ******************************************************************************

  input => pRegion%mixtInput

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

    DO fType = 1,FTYPE_MAX
      SELECT CASE ( fType )
        CASE ( FTYPE_MIXT )
          bc => pPatch%mixt
        CASE ( FTYPE_TURB )
          CYCLE
        CASE ( FTYPE_PLAG )
          CYCLE
        CASE ( FTYPE_SPEC )
          IF ( global%specUsed .EQV. .FALSE. ) THEN 
	    CYCLE
	  END IF ! global%specUsed
	  
          bc => pPatch%spec
        CASE ( FTYPE_RADI )
          IF ( input%radiUsed .EQV. .FALSE. ) THEN
            CYCLE
          END IF ! input%radiUsed
          
          bc => pPatch%valRadi
        CASE DEFAULT ! Required because FTYPE_MAX includes PEUL
          CYCLE
      END SELECT ! fType

      IF ( bc%nData > 0 ) THEN
        ALLOCATE(bc%tbcs(bc%nData),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
        END IF ! global%error

        bc%tbcs(:)%tbcType = TBC_NONE
      END IF ! bc%nData
    END DO ! fType
  END DO ! iPatch

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE RFLU_AllocateMemoryTbc

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_AllocateMemoryTbc.F90,v $
! Revision 1.11  2008/12/06 08:44:11  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:24  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2006/08/19 15:38:41  mparmar
! Renamed patch variables
!
! Revision 1.8  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.7  2006/01/06 06:58:37  wasistho
! cycled tbc for turbulence
!
! Revision 1.6  2004/10/19 19:24:04  haselbac
! Bug fix: Added CASE DEFAULT to deal with PEUL case
!
! Revision 1.5  2004/07/28 15:29:19  jferry
! created global variable for spec use
!
! Revision 1.4  2004/07/23 22:43:15  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.3  2004/03/05 22:09:01  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.2  2004/01/29 23:02:07  haselbac
! Fixed bug: Use specUsed instead of specModel
!
! Revision 1.1  2003/06/04 20:05:54  jferry
! re-worked implementation of TBCs in unstructured code
!
! ******************************************************************************







