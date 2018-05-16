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
! Purpose: Read in user input related to injection boundary condition for 
!   species.
!
! Description: None.
!
! Input: 
!   pRegion     Region pointer
!
! Output: None.
!
! Notes: 
!   1. Define additional keyword, SPEC_, which can be used to set a default 
!      for all species. Individual species can be overridden by specifying
!      the appropriate keyword SPECn. 
!
!******************************************************************************
!
! $Id: SPEC_RFLU_ReadBcInjectSection.F90,v 1.5 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE SPEC_RFLU_ReadBcInjectSection(pRegion)

  USE ModDataTypes
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModError
  USE ModParameters
  
  USE ModInterfaces, ONLY: MakeNumberedKeys,ReadPatchSection 

  IMPLICIT NONE

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region), POINTER :: pRegion

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: bcName,RCSIdentString
  CHARACTER(10), DIMENSION(:), ALLOCATABLE :: keys
  CHARACTER(256) :: fileName
  LOGICAL, DIMENSION(:), ALLOCATABLE :: defined
  INTEGER :: checkSum,distrib,errorFlag,ifc,iKey,iPatch,iPatchBeg,iPatchEnd, &
             iReg,iVal,nBFacesTot,nKeys
  REAL(RFREAL), DIMENSION(:), ALLOCATABLE :: vals
  TYPE(t_grid) :: grid
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: SPEC_RFLU_ReadBcInjectSection.F90,v $ $Revision: 1.5 $'

  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_ReadBcInjectSection',&
  'SPEC_RFLU_ReadBcInjectSection.F90')

! *****************************************************************************
! Allocate memory 
! *****************************************************************************

  nKeys = pRegion%specInput%nSpecies + 1 ! Add one because of default key

  ALLOCATE(keys(nKeys),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'keys')
  END IF ! global%error

  ALLOCATE(vals(nKeys),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'vals')
  END IF ! global%error  

  ALLOCATE(defined(nKeys),STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'defined')
  END IF ! global%error 

! *****************************************************************************
! Generate keys. NOTE first key is default key.
! *****************************************************************************

  keys(1) = 'SPEC_'
  
  CALL MakeNumberedKeys(keys,2,'SPEC',1,nKeys,1)

! *****************************************************************************
! Read section
! *****************************************************************************

  CALL ReadPatchSection(global,IF_INPUT,nKeys,keys,vals,iPatchBeg,iPatchEnd, &
                        distrib,fileName,bcName,defined)

! *****************************************************************************
! Check if specified number of patches exceeds available ones
! *****************************************************************************

  IF ( iPatchEnd > global%nPatches ) THEN 
    CALL ErrorStop(global,ERR_PATCH_RANGE,__LINE__)
  END IF ! iPatchEnd

! *****************************************************************************
! Set options and check if all necessary values defined
! *****************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

! =============================================================================
!   Check whether this global patch exists in this region
! =============================================================================

    IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. & 
         pPatch%iPatchGlobal <= iPatchEnd ) THEN

! -----------------------------------------------------------------------------
!     Set options. NOTE at this stage, have already read information for 
!     mixture, so for GENx runs, pPatch%bcCoupled has been set and checked. 
!     Hence can use it directly to determine whether need distribution.
! -----------------------------------------------------------------------------

      pPatch%spec%nData     = pRegion%specInput%nSpecies
      pPatch%spec%nSwitches = 0
      
      IF ( pPatch%bcCoupled == BC_BURNING ) THEN 
        pPatch%spec%distrib = BCDAT_DISTRIB ! MUST have distribution
      ELSE 
        pPatch%spec%distrib = distrib 
      END IF ! pPatch%bcCoupled

! -----------------------------------------------------------------------------
!     Check whether all values defined (if no default defined)
! -----------------------------------------------------------------------------

      IF ( defined(1) .EQV. .FALSE. ) THEN ! No default defined
        checkSum = 0
      
        DO iKey = 2,nKeys
          IF ( defined(iKey) .EQV. .TRUE. ) THEN 
            checkSum = checkSum + 1
          END IF ! defined
        END DO ! iKey
        
        IF ( checkSum /= pRegion%specInput%nSpecies ) THEN 
          CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)
        END IF ! checkSum 
      END IF ! defined(1)

    END IF ! pPatch%iPatchGlobal
  END DO ! iPatch

! *****************************************************************************
! Copy values/distribution to variables
! *****************************************************************************

  DO iPatch = 1,pRegion%grid%nPatches
    pPatch => pRegion%patches(iPatch)

! =============================================================================
!   Check whether this global patch exists in this region
! =============================================================================

    IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. & 
         pPatch%iPatchGlobal <= iPatchEnd ) THEN

! -----------------------------------------------------------------------------
!     Distribution 
! -----------------------------------------------------------------------------

      IF ( pPatch%spec%distrib == BCDAT_DISTRIB ) THEN
        nBFacesTot = pPatch%nBTrisTot + pPatch%nBQuadsTot
        
        ALLOCATE(pPatch%spec%vals(pPatch%spec%nData,nBFacesTot), & 
                 STAT=errorFlag)                              
        global%error = errorFlag         
        IF ( global%error /= 0 ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%spec%vals')
        END IF ! global        

! ----- If not coupled, get boundary data from file ---------------------------

        IF ( pPatch%bcCoupled == BC_NOT_COUPLED ) THEN          
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        ELSE ! initialize - important for GENX runs from scratch
          DO ifc = 1,nBFacesTot
            DO iVal = 1,pPatch%spec%nData
              pPatch%spec%vals(iVal,ifc) = 0.0_RFREAL
            END DO ! iVal
          END DO ! ifc
        END IF ! pPatch%bcCoupled

! -----------------------------------------------------------------------------
!     Constant value
! -----------------------------------------------------------------------------

      ELSE
        ALLOCATE(pPatch%spec%vals(pPatch%spec%nData,0:1), & 
                 STAT=errorFlag)
        global%error = errorFlag         
        IF ( global%error /= 0 ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%spec%vals')
        END IF ! global%error

        DO iVal = 1,pPatch%spec%nData
          IF ( defined(1+iVal) .EQV. .TRUE. ) THEN ! Set to input
            pPatch%spec%vals(iVal,0:1) = vals(1+iVal)                      
          ELSE ! Set to default
            pPatch%spec%vals(iVal,0:1) = vals(1)
          END IF ! defined          
        END DO ! iVal
      END IF  ! pPatch%spec%distrib
    END IF ! pPatch%iPatchGlobal
  END DO ! iPatch

! *****************************************************************************
! Deallocate memory 
! *****************************************************************************

  DEALLOCATE(keys,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'keys')
  END IF ! global%error

  DEALLOCATE(vals,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'vals')
  END IF ! global%error  

  DEALLOCATE(defined,STAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'defined')
  END IF ! global%error 
  
! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_RFLU_ReadBcInjectSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_ReadBcInjectSection.F90,v $
! Revision 1.5  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/08/19 15:40:32  mparmar
! Renamed patch variables
!
! Revision 1.2  2006/04/07 15:19:25  haselbac
! Removed tabs
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
!******************************************************************************







