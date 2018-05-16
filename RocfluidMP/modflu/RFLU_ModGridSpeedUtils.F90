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
! Purpose: Collection of utility routines for manipulating grid speeds.
!
! Description: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModGridSpeedUtils.F90,v 1.7 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModGridSpeedUtils

  USE ModDataTypes
  USE ModParameters
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModError

  IMPLICIT NONE
    
  PRIVATE
  PUBLIC :: RFLU_DecideNeedGridSpeeds, & 
            RFLU_DescaleGridSpeed, &
            RFLU_DescaleGridSpeeds, &
            RFLU_InitGridSpeedScaleFactor, &
            RFLU_ScaleGridSpeed, &
            RFLU_ScaleGridSpeeds, &     
            RFLU_SetGridSpeedScaleFactor
                
  CONTAINS
  






! ******************************************************************************
!
! Purpose: Decide whether need grid speeds.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  LOGICAL FUNCTION RFLU_DecideNeedGridSpeeds(pRegion)
 
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion
  
! ******************************************************************************
!   Decide whether need grid speeds
! ******************************************************************************

    RFLU_DecideNeedGridSpeeds = .FALSE. 
    
    IF ( pRegion%mixtInput%moveGrid .EQV. .TRUE. ) THEN 
      RFLU_DecideNeedGridSpeeds = .TRUE.
    END IF ! pRegion%mixtInput%moveGrid
 
! ******************************************************************************
!   End
! ******************************************************************************

  END FUNCTION RFLU_DecideNeedGridSpeeds









! ******************************************************************************
!
! Purpose: Descale grid speed for single face.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   fs                  Scaled grid speed
!
! Output: 
!   fs                  Descaled grid speed
!
! Notes: None.
!
! ******************************************************************************

  FUNCTION RFLU_DescaleGridSpeed(pRegion,fs)
 
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

    REAL(RFREAL) :: RFLU_DescaleGridSpeed

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    REAL(RFREAL), INTENT(INOUT) :: fs
    TYPE(t_region), POINTER :: pRegion
  
! ******************************************************************************
!   Scale grid speed
! ******************************************************************************

    RFLU_DescaleGridSpeed = fs/pRegion%grid%fsScaleFactor
 
! ******************************************************************************
!   End
! ******************************************************************************

  END FUNCTION RFLU_DescaleGridSpeed






! ******************************************************************************
!
! Purpose: Descale grid speeds for faces and boundary patches.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_DescaleGridSpeeds(pRegion)
 
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  

    INTEGER :: ifg,ifl,iPatch
    REAL(RFREAL) :: scaleFactor
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! *****************************************************************************
!   Set pointers and variables
! *****************************************************************************

    pGrid => pRegion%grid

    scaleFactor = 1.0_RFREAL/pGrid%fsScaleFactor
  
! ******************************************************************************
!   Scale grid speeds
! ******************************************************************************

    DO ifg = LBOUND(pGrid%gs,1),UBOUND(pGrid%gs,1)
      pGrid%gs(ifg) = scaleFactor*pGrid%gs(ifg)
    END DO ! ifg
 
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      DO ifl = LBOUND(pPatch%gs,1),UBOUND(pPatch%gs,1)
        pPatch%gs(ifl) = scaleFactor*pPatch%gs(ifl)
      END DO ! ifl
    END DO ! iPatch
 
! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_DescaleGridSpeeds






! ******************************************************************************
!
! Purpose: Initialize grid speed scale factor.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!
! Output: None.
!
! Notes: 
!   1. Needed for steady flows without grid motion so that calls to scaling 
!      routines work properly.
!
! ******************************************************************************

  SUBROUTINE RFLU_InitGridSpeedScaleFactor(pRegion)
 
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

          
! ******************************************************************************
!   Initialize scaling factor
! ******************************************************************************

    pRegion%grid%fsScaleFactor = 1.0_RFREAL
 
! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_InitGridSpeedScaleFactor







! ******************************************************************************
!
! Purpose: Scale grid speed for single face.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   fs                  Descaled grid speed
!
! Output: 
!   fs                  Scaled grid speed
!
! Notes: None.
!
! ******************************************************************************

  FUNCTION RFLU_ScaleGridSpeed(pRegion,fs)
 
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

    REAL(RFREAL) :: RFLU_ScaleGridSpeed

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    REAL(RFREAL), INTENT(INOUT) :: fs
    TYPE(t_region), POINTER :: pRegion
  
! ******************************************************************************
!   Scale grid speed
! ******************************************************************************

    RFLU_ScaleGridSpeed = pRegion%grid%fsScaleFactor*fs
 
! ******************************************************************************
!   End
! ******************************************************************************

  END FUNCTION RFLU_ScaleGridSpeed






! ******************************************************************************
!
! Purpose: Scale grid speeds for faces and boundary patches.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ScaleGridSpeeds(pRegion)
 
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  

    INTEGER :: ifg,ifl,iPatch
    REAL(RFREAL) :: scaleFactor
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! *****************************************************************************
!   Set pointers and variables
! *****************************************************************************

    pGrid => pRegion%grid

    scaleFactor = pGrid%fsScaleFactor
  
! ******************************************************************************
!   Scale grid speeds
! ******************************************************************************

    DO ifg = LBOUND(pGrid%gs,1),UBOUND(pGrid%gs,1) 
      pGrid%gs(ifg) = scaleFactor*pGrid%gs(ifg)
    END DO ! ifg
 
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      DO ifl = LBOUND(pPatch%gs,1),UBOUND(pPatch%gs,1) 
        pPatch%gs(ifl) = scaleFactor*pPatch%gs(ifl)
      END DO ! ifl
    END DO ! iPatch
 
! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_ScaleGridSpeeds






! ******************************************************************************
!
! Purpose: Set grid speed scale factor.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_SetGridSpeedScaleFactor(pRegion)
 
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  

    INTEGER :: iRk,iRkStep,nRkSteps
    REAL(RFREAL) :: term
    REAL(RFREAL) :: ark(5),grk(5)
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
        
! *****************************************************************************
!   Set pointers and variables
! *****************************************************************************

    global => pRegion%global
    pGrid  => pRegion%grid

    iRkStep  = pRegion%irkStep
    nRkSteps = pRegion%global%nrkSteps

    ark(:) = pRegion%mixtInput%ark(:)
    grk(:) = pRegion%mixtInput%grk(:)
  
! ******************************************************************************
!   Determine scaling factor
! ******************************************************************************

    IF ( iRkStep > 1 .AND. iRkStep < nRkSteps ) THEN 
      pGrid%fsScaleFactor   = ark(iRkStep-1)/ark(iRkStep)
    ELSE IF ( iRkStep == 1 ) THEN
      pGrid%fsScaleFactor   = 1.0_RFREAL/ark(iRkStep)
    ELSE IF ( iRkStep == nRkSteps ) THEN
      term = 0.0_RFREAL

      DO iRk = 1,nRkSteps-1
        term = term + grk(iRk)/ark(iRk)
      END DO ! iRk

      pGrid%fsScaleFactor   = (1.0_RFREAL/ark(nRkSteps) - term)/ark(iRkStep-1)
    ELSE ! Defensive programming
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! iRkStep
 
! ******************************************************************************
!   End
! ******************************************************************************

  END SUBROUTINE RFLU_SetGridSpeedScaleFactor




END MODULE RFLU_ModGridSpeedUtils

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModGridSpeedUtils.F90,v $
! Revision 1.7  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.4  2004/10/19 19:28:03  haselbac
! Added routine to decide whether need grid speeds
!
! Revision 1.3  2004/06/16 20:01:02  haselbac
! Modification of loop limits to increase efficiency
!
! Revision 1.2  2004/04/19 20:21:07  haselbac
! Bug fix: Missing indGs added
!
! Revision 1.1  2004/04/14 02:05:11  haselbac
! Initial revision
!
! ******************************************************************************






