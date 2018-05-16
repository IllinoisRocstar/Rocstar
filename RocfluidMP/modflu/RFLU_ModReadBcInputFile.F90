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
! Purpose: Collection of routines to read boundary condition input file.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModReadBcInputFile.F90,v 1.27 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModReadBcInputFile

  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_ReadBcInputFileWrapper

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  CHARACTER(CHRLEN) :: RCSIdentString = &
    '$RCSfile: RFLU_ModReadBcInputFile.F90,v $ $Revision: 1.27 $'

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS







! ******************************************************************************
!
! Purpose: Read in user input related to farfield boundary condition.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadBcFarfSection(pRegion)

    USE ModInterfaces, ONLY: ReadPatchSection

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER, PARAMETER :: NVALS_MAX = 13

    LOGICAL, DIMENSION(NVALS_MAX) :: defined
    CHARACTER(10) :: keys(NVALS_MAX)
    CHARACTER(256) :: iFileName
    CHARACTER(CHRLEN) :: bcName
    INTEGER :: checkSum,distrib,errorFlag,i,iData,ifl,iPatch,iPatchBeg, &
               iPatchEnd,nBFacesTot,nVals
    REAL(RFREAL), DIMENSION(NVALS_MAX) :: vals
    TYPE(t_grid) :: grid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadBcFarfSection',&
  'RFLU_ModReadBcInputFile.F90')

! ******************************************************************************
!   Specify keywords and search for them
! ******************************************************************************

    nVals = NVALS_MAX

    keys(1) = 'CORR'
    keys(2) = 'MACH'
    keys(3) = 'ATTACK'
    keys(4) = 'SLIP'
    keys(5) = 'PRESS'
    keys(6) = 'TEMP'
    keys(7) = 'MVPATCH'
    keys(8) = 'SMGRID'
    keys(9) = 'ORDER'
    keys(10) = 'MOVEDIR'
    keys(11) = 'COUPLED'
    keys(12) = 'KIND'
    keys(13) = 'THRUSTFLAG'

    CALL ReadPatchSection(global,IF_INPUT,nVals,keys,vals,iPatchBeg,iPatchEnd, &
                          distrib,iFileName,bcName,defined )

! ******************************************************************************
!   Check if specified number of patches exceeds available ones
! ******************************************************************************

    IF ( iPatchEnd > global%nPatches ) THEN
      CALL ErrorStop(global,ERR_PATCH_RANGE,__LINE__)
    END IF ! iPatchEnd

! ******************************************************************************
!   Get switches and check that all necessary values defined
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN
        pPatch%bcType = BC_FARFIELD
        pPatch%bcName = bcName
        
!        pPatch%bcCoupled    = BC_NOT_COUPLED
!        pPatch%movePatchDir = MOVEPATCH_DIR_NONE
        
        pPatch%cReconst      = CONSTR_NONE        
        pPatch%plotStatsFlag = .FALSE.

! TEMPORARY - No longer used, keep for backward compatibility 
        pPatch%movePatch  = .FALSE. 
        pPatch%smoothGrid = .FALSE. 
! END TEMPORARY         

! ------------------------------------------------------------------------------
!       initialize Boundary Condition kind
! ------------------------------------------------------------------------------

        IF ( defined(12) .EQV. .TRUE. ) THEN
          IF ( (NINT(vals(12)) >= BC_KIND_MIN) .AND. & 
               (NINT(vals(12)) <= BC_KIND_MAX) ) THEN
            pPatch%bcKind = NINT(vals(12))
          ELSE
! TEMPORARY : issue a warning here ...
            pPatch%bcKind = BC_KIND_SIMPLE ! Initialize with Default BC Kind
          END IF ! checking range of vals(12)
        ELSE
          pPatch%bcKind = BC_KIND_SIMPLE ! Default BC Kind
        END IF ! defined(12)

! ------------------------------------------------------------------------------
!       initialize patch thrustFlag
! ------------------------------------------------------------------------------

        IF ( defined(13) .EQV. .TRUE. ) THEN
          IF ( (vals(13) > 0.5_RFREAL) .AND. & 
               (vals(13) < 1.5_RFREAL) ) THEN
            pPatch%thrustFlag = .TRUE.
          ELSE
            pPatch%thrustFlag = .FALSE.
          END IF ! 
        ELSE
          pPatch%thrustFlag = .FALSE.
        END IF ! defined(13)

! ------------------------------------------------------------------------------
!       Set switches
! ------------------------------------------------------------------------------

        pPatch%mixt%nSwitches = 1

        ALLOCATE(pPatch%mixt%switches(pPatch%mixt%nSwitches), &
                 STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
        END IF ! global
        
! ------------------------------------------------------------------------------
!       Check if switches defined
! ------------------------------------------------------------------------------

        IF ( defined(1) .EQV. .TRUE. ) THEN
          IF ( (NINT(vals(1)) == BCOPT_CORR_YES) .OR. & 
               (NINT(vals(1)) == BCOPT_CORR_NO ) ) THEN
            pPatch%mixt%switches(BCSWI_FARF_CORR) = NINT(vals(1))
          ELSE
            CALL ErrorStop(global,ERR_VAL_BCSWITCH,__LINE__,'(farfield type).')
          END IF ! NINT
        ELSE
          CALL ErrorStop(global,ERR_NO_BCSWITCH,__LINE__,'(farfield type).')
        END IF ! defined       
        
! ------------------------------------------------------------------------------
!       Check whether appropriate values specified
! ------------------------------------------------------------------------------

        pPatch%mixt%nData   = 5
        pPatch%mixt%distrib = distrib

        IF ( pPatch%mixt%distrib == BCDAT_CONSTANT ) THEN
          checkSum = 0

          DO i = 0,pPatch%mixt%nData-1
            IF ( defined(2+i) .EQV. .TRUE. ) THEN
              checkSum = checkSum + 1
            END IF ! defined
          END DO ! i

          IF ( checkSum /= pPatch%mixt%nData ) THEN
            CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)
          END IF ! checkSum
        END IF ! pPatch%mixt%distrib
        
! ------------------------------------------------------------------------------
!       Set patch spatial order
! ------------------------------------------------------------------------------

        IF ( defined(9) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(9)) == 2 ) THEN  
            pPatch%spaceOrder = 2
          ELSE 
            pPatch%spaceOrder = 1
          END IF ! NINT(vals(9))  
        ELSE 
          pPatch%spaceOrder = pRegion%mixtInput%spaceOrderBFaces
        END IF ! defined                              

! ------------------------------------------------------------------------------
!       Set patch motion variable
! ------------------------------------------------------------------------------

        IF ( defined(10) .EQV. .TRUE. ) THEN
          pPatch%movePatchDir = vals(10)
        ELSE
          pPatch%movePatchDir = MOVEPATCH_DIR_NONE
        END IF ! defined

! ------------------------------------------------------------------------------
!       Set coupling variable
! ------------------------------------------------------------------------------

        IF ( defined(11) .EQV. .TRUE. ) THEN
          pPatch%bcCoupled = vals(11)
        ELSE
          pPatch%bcCoupled = BC_NOT_COUPLED
        END IF ! defined

      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! ******************************************************************************
!   Copy values/distribution to variables 
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!   Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN

! ------------------------------------------------------------------------------
!       Distribution from file: Allocate and initialize, actual values read in
!       at later stage
! ------------------------------------------------------------------------------

        IF ( pPatch%mixt%distrib == BCDAT_DISTRIB ) THEN
	  nBFacesTot = pPatch%nBTrisTot + pPatch%nBQuadsTot
	
          ALLOCATE(pPatch%mixt%vals(pPatch%mixt%nData,nBFacesTot), &
                   STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%vals')
          END IF ! global

          DO ifl = 1,nBFacesTot
            DO iData = 1,pPatch%mixt%nData
              pPatch%mixt%vals(iData,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
            END DO ! iData
          END DO ! ifl

! ------------------------------------------------------------------------------
!       Constant value
! ------------------------------------------------------------------------------

        ELSE
          ALLOCATE(pPatch%mixt%vals(pPatch%mixt%nData,0:1), &
                   STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%vals')
          END IF ! global

          pPatch%mixt%vals(BCDAT_FARF_MACH  ,0:1) = vals(2)
          pPatch%mixt%vals(BCDAT_FARF_ATTACK,0:1) = vals(3)*global%deg2rad
          pPatch%mixt%vals(BCDAT_FARF_SLIP  ,0:1) = vals(4)*global%deg2rad
          pPatch%mixt%vals(BCDAT_FARF_PRESS ,0:1) = vals(5)
          pPatch%mixt%vals(BCDAT_FARF_TEMP  ,0:1) = vals(6)
        END IF ! pPatch%mixt%distrib

      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadBcFarfSection







! ******************************************************************************
!
! Purpose: Read in user input related to inflow boundary condition based on 
!   total conditions and flow angles.
!
! Description: None.
!
! Input:
!   pRegion        Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadBcInflowTotAngSection(pRegion)

    USE ModInterfaces, ONLY: ReadPatchSection

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER, PARAMETER :: NVALS_MAX = 15

    LOGICAL :: iFileExists
    LOGICAL, DIMENSION(NVALS_MAX) :: defined
    CHARACTER(10) :: keys(NVALS_MAX)
    CHARACTER(256) :: iFileName
    CHARACTER(CHRLEN) :: bcName
    INTEGER :: checkSum,distrib,errorFlag,i,iData,ifl,iPatch,iReg,nVals, &
               iPatchBeg,iPatchEnd,nBFacesTot,switch
    REAL(RFREAL), DIMENSION(NVALS_MAX) :: vals
    TYPE(t_grid) :: grid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadBcInflowTotAngSection',&
  'RFLU_ModReadBcInputFile.F90')

! ******************************************************************************
!   Specify keywords and search for them
! ******************************************************************************

    nVals = NVALS_MAX

    keys(1) = 'TYPE'
    keys(2) = 'FIXED'
    keys(3) = 'PTOT'
    keys(4) = 'TTOT'
    keys(5) = 'BETAH'
    keys(6) = 'BETAV'
    keys(7) = 'MACH'
    keys(8) = 'MVPATCH'
    keys(9) = 'SMGRID'
    keys(10) = 'MOVEDIR'
    keys(11) = 'COUPLED'
    keys(12) = 'KIND'
    keys(13) = 'REFLECT'
    keys(14) = 'THRUSTFLAG'
    keys(15) = 'ORDER'

    CALL ReadPatchSection(global,IF_INPUT,nVals,keys,vals,iPatchBeg,iPatchEnd, &
                          distrib,iFileName,bcName,defined)

! ******************************************************************************
!   Check if specified number of patches exceeds available ones
! ******************************************************************************

    IF ( iPatchEnd > global%nPatches ) THEN
      CALL ErrorStop(global,ERR_PATCH_RANGE,__LINE__)
    END IF ! iPatchEnd

! ******************************************************************************
!   Get switches and check that all necessary values defined
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN
        pPatch%bcType = BC_INFLOW_TOTANG
        pPatch%bcName = bcName
        
!        pPatch%bcCoupled    = BC_NOT_COUPLED
!        pPatch%movePatchDir = MOVEPATCH_DIR_NONE
        
        pPatch%cReconst      = CONSTR_NONE           
        pPatch%plotStatsFlag = .FALSE.  

! TEMPORARY - No longer used, keep for backward compatibility 
        pPatch%movePatch  = .FALSE. 
        pPatch%smoothGrid = .FALSE. 
! END TEMPORARY 

        pPatch%mixt%distrib = distrib

! ------------------------------------------------------------------------------
!       Set switches
! ------------------------------------------------------------------------------

        pPatch%mixt%nSwitches = 2

        ALLOCATE(pPatch%mixt%switches(pPatch%mixt%nSwitches), &
                 STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
        END IF ! global

! ------------------------------------------------------------------------------
!       Check if switches defined
! ------------------------------------------------------------------------------

        IF ( defined(1) .EQV. .TRUE. ) THEN
          IF ( (NINT(vals(1)) >= BCOPT_SUPERSONIC) .AND. &
               (NINT(vals(1)) <= BCOPT_MIXED     ) ) THEN
            pPatch%mixt%switches(BCSWI_INFLOW_TYPE) = NINT(vals(1))
          ELSE
            CALL ErrorStop(global,ERR_VAL_BCSWITCH,__LINE__,'(inflow type).')
          END IF ! NINT
        ELSE
          CALL ErrorStop(global,ERR_NO_BCSWITCH,__LINE__,'(inflow type).')
        END IF ! defined

        IF ( defined(2) .EQV. .TRUE. ) THEN
          IF ( (NINT(vals(2)) >= BCOPT_FIXED_NO ) .AND. &
               (NINT(vals(2)) <= BCOPT_FIXED_YES) ) THEN
            pPatch%mixt%switches(BCSWI_INFLOW_FIXED) = NINT(vals(2))
          ELSE
            CALL ErrorStop(global,ERR_VAL_BCSWITCH,__LINE__,'(fixed).')
          END IF ! NINT
        ELSE
          pPatch%mixt%switches(BCSWI_INFLOW_FIXED) = BCOPT_FIXED_NO
        END IF ! defined

! ------------------------------------------------------------------------------
!       initialize Boundary Condition kind
! ------------------------------------------------------------------------------

        IF ( defined(12) .EQV. .TRUE. ) THEN
          IF ( (NINT(vals(12)) >= BC_KIND_MIN) .AND. & 
               (NINT(vals(12)) <= BC_KIND_MAX) ) THEN
            pPatch%bcKind = NINT(vals(12))
          ELSE
! TEMPORARY : issue a warning here ...
            pPatch%bcKind = BC_KIND_SIMPLE ! Initialize with Default BC Kind
          END IF ! checking range of vals(12)
        ELSE
          pPatch%bcKind = BC_KIND_SIMPLE ! Default BC Kind
        END IF ! defined(12)

! ------------------------------------------------------------------------------
!       initialize if inflow BC is reflecting or non reflecting
! ------------------------------------------------------------------------------

        IF ( defined(13) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(13)) == BC_REFLECTING ) THEN 
            pPatch%reflect = BC_REFLECTING
          ELSE
            pPatch%reflect = BC_NONREFLECTING
          END IF ! 
        ELSE
          pPatch%reflect = BC_REFLECTING
        END IF ! defined(13)

! ------------------------------------------------------------------------------
!       initialize patch thrustFlag
! ------------------------------------------------------------------------------

        IF ( defined(14) .EQV. .TRUE. ) THEN
          IF ( (vals(14) > 0.5_RFREAL) .AND. & 
               (vals(14) < 1.5_RFREAL) ) THEN
            pPatch%thrustFlag = .TRUE.
          ELSE
            pPatch%thrustFlag = .FALSE.
          END IF ! 
        ELSE
          pPatch%thrustFlag = .FALSE.
        END IF ! defined(14)

! ------------------------------------------------------------------------------
!       Set patch spatial order
! ------------------------------------------------------------------------------

        IF ( defined(15) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(15)) == 2 ) THEN
            pPatch%spaceOrder = 2
          ELSE
            pPatch%spaceOrder = 1
          END IF ! NINT(vals(15))
        ELSE
          pPatch%spaceOrder = pRegion%mixtInput%spaceOrderBFaces
        END IF ! defined

! ------------------------------------------------------------------------------
!       Check whether appropriate values specified: for subsonic inflow need 4
!       quantities, for supersonic and mixed inflow need 5 quantities
! ------------------------------------------------------------------------------

        IF ( pPatch%mixt%switches(BCSWI_INFLOW_TYPE) /= BCOPT_SUBSONIC ) THEN
          pPatch%mixt%nData = 5
        ELSE
          pPatch%mixt%nData = 4
        END IF ! pPatch

        IF ( pPatch%mixt%distrib == BCDAT_CONSTANT ) THEN
          checkSum = 0

          DO i = 0,pPatch%mixt%nData-1
            IF ( defined(3+i) .EQV. .TRUE. ) THEN
              checkSum = checkSum + 1
            END IF ! defined
          END DO ! i

          IF ( checkSum /= pPatch%mixt%nData ) THEN
            CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)
          END IF ! checkSum
        END IF ! pPatch

! ------------------------------------------------------------------------------
!       Set patch motion variable
! ------------------------------------------------------------------------------

        IF ( defined(10) .EQV. .TRUE. ) THEN
          pPatch%movePatchDir = vals(10)
        ELSE
          pPatch%movePatchDir = MOVEPATCH_DIR_NONE
        END IF ! defined

! ------------------------------------------------------------------------------
!       Set coupling variable
! ------------------------------------------------------------------------------

        IF ( defined(11) .EQV. .TRUE. ) THEN
          pPatch%bcCoupled = vals(11)
        ELSE
          pPatch%bcCoupled = BC_NOT_COUPLED

        END IF ! defined
      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! ******************************************************************************
!   Copy values/distribution to variables for non-adiabatic walls
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN
        switch = pPatch%mixt%switches(BCSWI_INFLOW_TYPE)

! ------------------------------------------------------------------------------
!       Distribution from file: Allocate and initialize, actual values read in
!       at later stage
! ------------------------------------------------------------------------------

        IF ( pPatch%mixt%distrib == BCDAT_DISTRIB ) THEN
	  nBFacesTot = pPatch%nBTrisTot + pPatch%nBQuadsTot
	
          ALLOCATE(pPatch%mixt%vals(pPatch%mixt%nData,nBFacesTot), &
                   STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%vals')
          END IF ! global

          DO ifl = 1,nBFacesTot
            DO iData = 1,pPatch%mixt%nData
              pPatch%mixt%vals(iData,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
            END DO ! iData
          END DO ! ifl

! ------------------------------------------------------------------------------
!       Constant value
! ------------------------------------------------------------------------------

        ELSE
          ALLOCATE(pPatch%mixt%vals(pPatch%mixt%nData,0:1),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%vals')
          END IF ! global

          pPatch%mixt%vals(BCDAT_INFLOW_PTOT, 0:1) = vals(3)
          pPatch%mixt%vals(BCDAT_INFLOW_TTOT, 0:1) = vals(4)
          pPatch%mixt%vals(BCDAT_INFLOW_BETAH,0:1) = vals(5)*global%deg2rad
          pPatch%mixt%vals(BCDAT_INFLOW_BETAV,0:1) = vals(6)*global%deg2rad

          IF ( switch /= BCOPT_SUBSONIC ) THEN
            pPatch%mixt%vals(BCDAT_INFLOW_MACH,0:1) = vals(7)
          END IF ! switch
        END IF  ! pPatch%mixt%distrib

      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadBcInflowTotAngSection









! ******************************************************************************
!
! Purpose: Read in user input related to inflow boundary condition based on 
!   velocities and temperature.
!
! Description: None.
!
! Input:
!   pRegion        Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadBcInflowVelTempSection(pRegion)

    USE ModInterfaces, ONLY: ReadPatchSection

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER, PARAMETER :: NVALS_MAX = 14

    LOGICAL :: iFileExists
    LOGICAL, DIMENSION(NVALS_MAX) :: defined
    CHARACTER(10) :: keys(NVALS_MAX)
    CHARACTER(256) :: iFileName
    CHARACTER(CHRLEN) :: bcName
    INTEGER :: checkSum,distrib,errorFlag,i,iData,ifl,iPatch,iReg,nBFacesTot, & 
               nVals,iPatchBeg,iPatchEnd,switch
    REAL(RFREAL), DIMENSION(NVALS_MAX) :: vals
    TYPE(t_grid) :: grid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadBcInflowVelTempSection',&
  'RFLU_ModReadBcInputFile.F90')

! ******************************************************************************
!   Specify keywords and search for them
! ******************************************************************************

    nVals = NVALS_MAX

    keys(1)  = 'TYPE'
    keys(2)  = 'VELX'
    keys(3)  = 'VELY'
    keys(4)  = 'VELZ'
    keys(5)  = 'TEMP'
    keys(6)  = 'PRESS'
    keys(7)  = 'MVPATCH'
    keys(8)  = 'SMGRID'
    keys(9)  = 'MOVEDIR'
    keys(10) = 'COUPLED'
    keys(11) = 'KIND' 
    keys(12) = 'REFLECT'
    keys(13) = 'THRUSTFLAG'
    keys(14) = 'ORDER'

    CALL ReadPatchSection(global,IF_INPUT,nVals,keys,vals,iPatchBeg,iPatchEnd, &
                          distrib,iFileName,bcName,defined)

! ******************************************************************************
!   Check if specified number of patches exceeds available ones
! ******************************************************************************

    IF ( iPatchEnd > global%nPatches ) THEN
      CALL ErrorStop(global,ERR_PATCH_RANGE,__LINE__)
    END IF ! iPatchEnd

! ******************************************************************************
!   Get switches and check that all necessary values defined
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN
        pPatch%bcType = BC_INFLOW_VELTEMP
        pPatch%bcName = bcName
        
!        pPatch%bcCoupled    = BC_NOT_COUPLED
!        pPatch%movePatchDir = MOVEPATCH_DIR_NONE
        
        pPatch%cReconst      = CONSTR_NONE           
        pPatch%plotStatsFlag = .FALSE.

! TEMPORARY - No longer used, keep for backward compatibility 
        pPatch%movePatch  = .FALSE. 
        pPatch%smoothGrid = .FALSE. 
! END TEMPORARY 

        pPatch%mixt%distrib = distrib

! ------------------------------------------------------------------------------
!       Set switches
! ------------------------------------------------------------------------------

        pPatch%mixt%nSwitches = 1

        ALLOCATE(pPatch%mixt%switches(pPatch%mixt%nSwitches), &
                 STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
        END IF ! global

! ------------------------------------------------------------------------------
!       Check if switches defined
! ------------------------------------------------------------------------------

        IF ( defined(1) .EQV. .TRUE. ) THEN
          IF ( (NINT(vals(1)) >= BCOPT_SUPERSONIC) .AND. &
               (NINT(vals(1)) <= BCOPT_MIXED     ) ) THEN
            pPatch%mixt%switches(BCSWI_INFLOW_TYPE) = NINT(vals(1))
          ELSE
            CALL ErrorStop(global,ERR_VAL_BCSWITCH,__LINE__,'(inflow type).')
          END IF ! NINT
        ELSE
          CALL ErrorStop(global,ERR_NO_BCSWITCH,__LINE__,'(inflow type).')
        END IF ! defined

! ------------------------------------------------------------------------------
!       initialize Boundary Condition kind
! ------------------------------------------------------------------------------

        IF ( defined(11) .EQV. .TRUE. ) THEN
          IF ( (NINT(vals(11)) >= BC_KIND_MIN) .AND. & 
               (NINT(vals(11)) <= BC_KIND_MAX) ) THEN
            pPatch%bcKind = NINT(vals(11))
          ELSE
! TEMPORARY : issue a warning here ...
            pPatch%bcKind = BC_KIND_SIMPLE ! Initialize with Default BC Kind
          END IF ! checking range of vals(11)
        ELSE
          pPatch%bcKind = BC_KIND_SIMPLE ! Default BC Kind
        END IF ! defined(11)

! ------------------------------------------------------------------------------
!       initialize if inflow BC is reflecting or non reflecting
! ------------------------------------------------------------------------------

        IF ( defined(12) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(12)) == BC_REFLECTING ) THEN 
            pPatch%reflect = BC_REFLECTING
          ELSE
            pPatch%reflect = BC_NONREFLECTING
          END IF ! 
        ELSE
          pPatch%reflect = BC_REFLECTING
        END IF ! defined(12)

! ------------------------------------------------------------------------------
!       initialize patch thrustFlag
! ------------------------------------------------------------------------------

        IF ( defined(13) .EQV. .TRUE. ) THEN
          IF ( (vals(13) > 0.5_RFREAL) .AND. & 
               (vals(13) < 1.5_RFREAL) ) THEN
            pPatch%thrustFlag = .TRUE.
          ELSE
            pPatch%thrustFlag = .FALSE.
          END IF ! 
        ELSE
          pPatch%thrustFlag = .FALSE.
        END IF ! defined(13)

! ------------------------------------------------------------------------------
!       Set patch spatial order
! ------------------------------------------------------------------------------

        IF ( defined(14) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(14)) == 2 ) THEN
            pPatch%spaceOrder = 2
          ELSE
            pPatch%spaceOrder = 1
          END IF ! NINT(vals(14))
        ELSE
          pPatch%spaceOrder = pRegion%mixtInput%spaceOrderBFaces
        END IF ! defined

! ------------------------------------------------------------------------------
!       Check whether appropriate values specified: for subsonic inflow need 4
!       quantities, for supersonic and mixed inflow need 5 quantities
! ------------------------------------------------------------------------------

        IF ( pPatch%mixt%switches(BCSWI_INFLOW_TYPE) /= BCOPT_SUBSONIC ) THEN
          pPatch%mixt%nData = 5
        ELSE
          pPatch%mixt%nData = 4
        END IF ! pPatch

        IF ( pPatch%mixt%distrib == BCDAT_CONSTANT ) THEN
          checkSum = 0

          DO i = 0,pPatch%mixt%nData-1
            IF ( defined(2+i) .EQV. .TRUE. ) THEN
              checkSum = checkSum + 1
            END IF ! defined
          END DO ! i

          IF ( checkSum /= pPatch%mixt%nData ) THEN
            CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)
          END IF ! checkSum
        END IF ! pPatch

! ------------------------------------------------------------------------------
!       Set patch motion variable
! ------------------------------------------------------------------------------

        IF ( defined(9) .EQV. .TRUE. ) THEN
          pPatch%movePatchDir = vals(9)
        ELSE
          pPatch%movePatchDir = MOVEPATCH_DIR_NONE
        END IF ! defined

! ------------------------------------------------------------------------------
!       Set coupling variable
! ------------------------------------------------------------------------------

        IF ( defined(10) .EQV. .TRUE. ) THEN
          pPatch%bcCoupled = vals(10)
        ELSE
          pPatch%bcCoupled = BC_NOT_COUPLED
        END IF ! defined

      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! ******************************************************************************
!   Copy values/distribution to variables for non-adiabatic walls
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN
        switch = pPatch%mixt%switches(BCSWI_INFLOW_TYPE)

! ------------------------------------------------------------------------------
!       Distribution from file: Allocate and initialize, actual values read in
!       at later stage
! ------------------------------------------------------------------------------

        IF ( pPatch%mixt%distrib == BCDAT_DISTRIB ) THEN
	  nBFacesTot = pPatch%nBTrisTot + pPatch%nBQuadsTot
	
          ALLOCATE(pPatch%mixt%vals(pPatch%mixt%nData,nBFacesTot), &
                   STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%vals')
          END IF ! global

          DO ifl = 1,nBFacesTot
            DO iData = 1,pPatch%mixt%nData
              pPatch%mixt%vals(iData,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
            END DO ! iData
          END DO ! ifl

! ------------------------------------------------------------------------------
!       Constant value
! ------------------------------------------------------------------------------

        ELSE
          ALLOCATE(pPatch%mixt%vals(pPatch%mixt%nData,0:1),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%vals')
          END IF ! global

          pPatch%mixt%vals(BCDAT_INFLOW_U,0:1) = vals(2)
          pPatch%mixt%vals(BCDAT_INFLOW_V,0:1) = vals(3)
          pPatch%mixt%vals(BCDAT_INFLOW_W,0:1) = vals(4)
          pPatch%mixt%vals(BCDAT_INFLOW_T,0:1) = vals(5)

          IF ( switch /= BCOPT_SUBSONIC ) THEN
            pPatch%mixt%vals(BCDAT_INFLOW_P,0:1) = vals(6)
          END IF ! switch
        END IF  ! pPatch%mixt%distrib

      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadBcInflowVelTempSection










! ******************************************************************************
!
! Purpose: Read in user input related to injection boundary condition.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadBcInjectSection(pRegion)

    USE ModInterfaces, ONLY: ReadPatchSection

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER, PARAMETER :: NVALS_MAX = 12

    LOGICAL, DIMENSION(NVALS_MAX) :: defined
    CHARACTER(10) :: keys(NVALS_MAX)
    CHARACTER(256) :: fname
    CHARACTER(CHRLEN) :: bcName
    INTEGER :: checkSum,distrib,errorFlag,i,iData,ifl,iPatch,iPatchBeg, &
               iPatchEnd,iReg,nBFacesTot,nVals
    REAL(RFREAL), DIMENSION(NVALS_MAX) :: vals
    TYPE(t_grid) :: grid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadBcInjectSection',&
  'RFLU_ModReadBcInputFile.F90')

! ******************************************************************************
!   Specify keywords and search for them
! ******************************************************************************

    nVals = NVALS_MAX

    keys(1) = 'MFRATE'
    keys(2) = 'TEMP'
    keys(3) = 'COUPLED'
    keys(4) = 'BFLAG'
    keys(5) = 'MVPATCH'
    keys(6) = 'SMGRID'
    keys(7) = 'MOVEDIR'
    keys(8) = 'STATS'
    keys(9) = 'CRECONST'
    keys(10)= 'KIND'
    keys(11)= 'THRUSTFLAG'
    keys(12)= 'ORDER'

    CALL ReadPatchSection(global,IF_INPUT,nVals,keys,vals,iPatchBeg,iPatchEnd, &
                          distrib,fname,bcName,defined )

! ******************************************************************************
!   Check if specified number of patches exceeds available ones
! ******************************************************************************

    IF ( iPatchEnd > global%nPatches ) THEN
      CALL ErrorStop(global,ERR_PATCH_RANGE,__LINE__)
    END IF ! iPatchEnd

! ******************************************************************************
!   Get switches and check that all necessary values defined
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN
        pPatch%bcType = BC_INJECTION
        pPatch%bcName = bcName

! TEMPORARY - No longer used, keep for backward compatibility 
        pPatch%movePatch  = .FALSE. 
        pPatch%smoothGrid = .FALSE. 
! END TEMPORARY 

        pPatch%mixt%nData     = 2
        pPatch%mixt%nSwitches = 0

! ------------------------------------------------------------------------------
!       initialize Boundary Condition kind
! ------------------------------------------------------------------------------

        IF ( defined(10) .EQV. .TRUE. ) THEN
          IF ( (NINT(vals(10)) >= BC_KIND_MIN) .AND. & 
               (NINT(vals(10)) <= BC_KIND_MAX) ) THEN
            pPatch%bcKind = NINT(vals(10))
          ELSE
! TEMPORARY : issue a warning here ...
            pPatch%bcKind = BC_KIND_SIMPLE ! Initialize with Default BC Kind
          END IF ! checking range of vals(10)
        ELSE
          pPatch%bcKind = BC_KIND_SIMPLE ! Default BC Kind
        END IF ! defined(10)

! ------------------------------------------------------------------------------
!       initialize patch thrustFlag
! ------------------------------------------------------------------------------

        IF ( defined(11) .EQV. .TRUE. ) THEN
          IF ( (vals(11) > 0.5_RFREAL) .AND. & 
               (vals(11) < 1.5_RFREAL) ) THEN
            pPatch%thrustFlag = .TRUE.
          ELSE
            pPatch%thrustFlag = .FALSE.
          END IF ! 
        ELSE
          pPatch%thrustFlag = .FALSE.
        END IF ! defined(11)

! ------------------------------------------------------------------------------
!       Set patch spatial order
! ------------------------------------------------------------------------------

        IF ( defined(12) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(12)) == 2 ) THEN
            pPatch%spaceOrder = 2
          ELSE
            pPatch%spaceOrder = 1
          END IF ! NINT(vals(12))
        ELSE
          pPatch%spaceOrder = pRegion%mixtInput%spaceOrderBFaces
        END IF ! defined

! ------------------------------------------------------------------------------
!       Check whether coupling flag defined
! ------------------------------------------------------------------------------

#ifdef GENX
        IF ( defined(3) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(3)) == BC_BURNING ) THEN
            pPatch%bcCoupled       = NINT(vals(3))
            pPatch%mixt%distrib = BCDAT_DISTRIB ! MUST have distribution
          ELSE IF ( NINT(vals(3)) == BC_NOT_BURNING ) THEN 
            global%warnCounter = global%warnCounter + 1

            WRITE(STDOUT,'(A,3X,A,1X,A,I3,A)') SOLVER_NAME,'*** WARNING ***', &
                  'Invalid coupling input for patch ',iPatch, &
                  '. Overriding coupling input.'
		  
            pPatch%bcCoupled       = BC_BURNING	  
            pPatch%mixt%distrib = BCDAT_DISTRIB ! MUST have distribution
	  ELSE IF ( NINT(vals(3)) == BC_NOT_COUPLED ) THEN 
            pPatch%bcCoupled       = NINT(vals(3))	    
            pPatch%mixt%distrib = distrib
	  ELSE 
            global%warnCounter = global%warnCounter + 1

            WRITE(STDOUT,'(A,3X,A,1X,A,I3,A)') SOLVER_NAME,'*** WARNING ***', &
                  'Invalid coupling input for patch ',iPatch, &
                  '. Overriding coupling input.'
                  
            pPatch%bcCoupled       = BC_NOT_COUPLED 
            pPatch%mixt%distrib = distrib	    
          END IF ! NINT(vals)
        ELSE
          pPatch%bcCoupled       = BC_NOT_COUPLED
          pPatch%mixt%distrib = distrib
        END IF ! defined
#else   
        IF ( defined(3) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(3)) /= BC_NOT_COUPLED ) THEN
            global%warnCounter = global%warnCounter + 1
          
            WRITE(STDOUT,'(A,3X,A,1X,A,I3,A)') SOLVER_NAME,'*** WARNING ***', &
                  'Invalid coupling input for patch ',iPatch, &
                  '. Overriding user input.'                      

            pPatch%bcCoupled = BC_NOT_COUPLED
          ELSE
            pPatch%bcCoupled = BC_NOT_COUPLED
	  END IF ! NINT(vals(3))
	  
          pPatch%mixt%distrib = distrib        
        ELSE 
          pPatch%bcCoupled       = BC_NOT_COUPLED
          pPatch%mixt%distrib = BCDAT_CONSTANT
	END IF ! defined(3)			
#endif

#ifdef GENX
! ------------------------------------------------------------------------------
!       Check whether bFlagInit defined
! ------------------------------------------------------------------------------

        IF ( defined(4) .EQV. .TRUE. ) THEN
          pPatch%bFlagInit = vals(4)

          IF ( pPatch%bFlagInit < 0 .OR. pPatch%bFlagInit > 1 ) THEN
            global%warnCounter = global%warnCounter + 1

            WRITE(STDOUT,'(A,3X,A,1X,A,I3,A)') SOLVER_NAME,'*** WARNING ***', &
                  'Invalid burning flag input for patch ',iPatch, &
                  '. Overriding burning flag input.'
            pPatch%bFlagInit = 1
          END IF ! pPatch%bFlagInit
        ELSE
          global%warnCounter = global%warnCounter + 1

          WRITE(STDOUT,'(A,3X,A,1X,A,I3,A)') SOLVER_NAME,'*** WARNING ***', &
                'Missing burning flag input for patch ',iPatch, &
                '. Setting burning flag input.'
          pPatch%bFlagInit = 1
        END IF ! defined
#endif
 
! ------------------------------------------------------------------------------
!       Check that all flags were set
! ------------------------------------------------------------------------------

        IF ( pPatch%mixt%distrib == BCDAT_CONSTANT ) THEN
          checkSum = 0

          DO i = 0,pPatch%mixt%nData-1
            IF ( defined(1+i) .EQV. .TRUE. ) THEN
              checkSum = checkSum + 1
            END IF ! defined
          END DO ! i

          IF ( checkSum /= pPatch%mixt%nData ) THEN
            CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)
          END IF ! checkSum
        END IF ! patch

! ------------------------------------------------------------------------------
!       Set patch motion variable
! ------------------------------------------------------------------------------

        IF ( defined(7) .EQV. .TRUE. ) THEN 
          pPatch%movePatchDir = vals(7)
        ELSE 
          pPatch%movePatchDir = MOVEPATCH_DIR_NONE
        END IF ! defined
        
! ------------------------------------------------------------------------------
!       Set patch statistics plotting flag
! ------------------------------------------------------------------------------

        IF ( defined(8) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(8)) == 1 ) THEN  
            pPatch%plotStatsFlag = .TRUE.
          ELSE
            pPatch%plotStatsFlag = .FALSE.
          END IF ! NINT
        ELSE 
          pPatch%plotStatsFlag = .FALSE.
        END IF ! defined 
        
! ------------------------------------------------------------------------------
!       Set constraint variable
! ------------------------------------------------------------------------------

        IF ( defined(9) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(9)) == 1 ) THEN  
            pPatch%cReconst = CONSTR_WEIGHTED
          ELSE 
            pPatch%cReconst = CONSTR_NONE
          END IF ! NINT(vals(9))  
        ELSE 
          pPatch%cReconst = CONSTR_NONE
        END IF ! defined                               
      END IF ! pPatch%iPatchGlobal            
    END DO ! iPatch

! ******************************************************************************
!   Copy values/distribution to variables 
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN

! ------------------------------------------------------------------------------
!       Distribution from file: Allocate and initialize, actual values read or
!       obtained at later stage
! ------------------------------------------------------------------------------

        IF ( pPatch%mixt%distrib == BCDAT_DISTRIB ) THEN
          nBFacesTot = pPatch%nBTrisTot + pPatch%nBQuadsTot

          ALLOCATE(pPatch%mixt%vals(pPatch%mixt%nData,nBFacesTot), &
                   STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%vals')
          END IF ! global

! ------- If not coupled, get boundary data from file --------------------------

          IF ( pPatch%bcCoupled == BC_NOT_COUPLED ) THEN
            DO ifl = 1,nBFacesTot
              DO iData = 1,pPatch%mixt%nData
                pPatch%mixt%vals(iData,ifl) = &
                  REAL(CRAZY_VALUE_INT,KIND=RFREAL)
              END DO ! iData
            END DO ! ifl

! ------- If coupled, get boundary data from Roccom ----------------------------

          ELSE ! Initialize - important for GENX runs from scratch
            DO ifl = 1,nBFacesTot
              DO iData = 1,pPatch%mixt%nData
                pPatch%mixt%vals(iData,ifl) = 0.0_RFREAL
              END DO ! iData
            END DO ! ifl
          END IF ! pPatch%bcCoupled

! ------------------------------------------------------------------------------
!       Constant value
! ------------------------------------------------------------------------------

        ELSE
          ALLOCATE(pPatch%mixt%vals(pPatch%mixt%nData,0:1), &
                   STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%vals')
          END IF ! global

          pPatch%mixt%vals(BCDAT_INJECT_MFRATE,0:1) = vals(1)
          pPatch%mixt%vals(BCDAT_INJECT_TEMP  ,0:1) = vals(2)
        END IF ! pPatch%mixt%distrib

      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadBcInjectSection








! ******************************************************************************
!
! Purpose: Read in user input related to boundary conditions (done on all
!   processors).
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadBcInputFile(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNamePlain

    IMPLICIT NONE

! ******************************************************************************
!   Definitions and declarations
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Locals
! ==============================================================================

    CHARACTER(CHRLEN) :: bcCoupledString,iFileName,moveString,smoothString
    CHARACTER(256) :: line
    INTEGER :: errorFlag,iPatch,loopCounter
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_grid) :: grid
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadBcInputFile',&
  'RFLU_ModReadBcInputFile.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_LOW ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Reading Rocflu boundary condition file...'
    END IF ! global%verbLevel

! ******************************************************************************
!   Open file
! ******************************************************************************

    CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.bc',iFileName)

    OPEN(IF_INPUT,FILE=iFileName,FORM='FORMATTED',STATUS='OLD',IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error

! ******************************************************************************
!   Read file looking for keywords
! ******************************************************************************

    loopCounter = 0

    keyWordLoop: DO
      READ(IF_INPUT,'(A256)',IOSTAT=errorFlag) line

      IF ( errorFlag > 0 ) THEN ! Error occurred
        CALL ErrorStop(global,ERR_FILE_READ,__LINE__,'File: '//TRIM(iFileName))
      ELSE IF ( errorFlag < 0 ) THEN ! Encountered end of file
        EXIT keyWordLoop
      END IF ! errorFlag

      SELECT CASE( TRIM(line) )
        CASE ('# BC_FARF')
          CALL RFLU_ReadBcFarfSection(pRegion)
! TEMPORARY - Keep this for backward compatibility
        CASE ('# BC_INFLOW')
          CALL RFLU_ReadBcInflowTotAngSection(pRegion)
! END TEMPORARY              
        CASE ('# BC_INFLOW_TOTANG')
          CALL RFLU_ReadBcInflowTotAngSection(pRegion)
        CASE ('# BC_INFLOW_VELTEMP')
          CALL RFLU_ReadBcInflowVelTempSection(pRegion)
        CASE ('# BC_INJECT')
          CALL RFLU_ReadBcInjectSection(pRegion)                                
! TEMPORARY - Keep this for backward compatibility
        CASE ('# BC_NOSLIP')
          CALL RFLU_ReadBcNoSlipWallHeatSect(pRegion)
! END TEMPORARY    
        CASE ('# BC_NOSLIP_HFLUX')
          CALL RFLU_ReadBcNoSlipWallHeatSect(pRegion)
        CASE ('# BC_NOSLIP_TEMP')
          CALL RFLU_ReadBcNoSlipWallTempSect(pRegion) 
        CASE ('# BC_OUTFLOW')
          CALL RFLU_ReadBcOutflowSection(pRegion)
        CASE ('# BC_PERIODIC')
          CALL RFLU_ReadBcPeriodicSection(pRegion)
        CASE ('# BC_SLIPW')
          CALL RFLU_ReadBcSlipWallSection(pRegion)
        CASE ('# BC_SYMMETRY')
          CALL RFLU_ReadBcSymmetrySection(pRegion)          
        CASE ('# BC_VIRTUAL')
          CALL RFLU_ReadBcVirtualSection(pRegion)          
        CASE ('# END')
          EXIT keyWordLoop
      END SELECT ! TRIM(line)

      loopCounter = loopCounter + 1 
      
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN ! Prevent infinite loop
        CALL ErrorStop(global,ERR_INFINITE_LOOP ,__LINE__)
      END IF ! loopCounter
    END DO keyWordLoop

! ******************************************************************************
!   Close file
! ******************************************************************************

    CLOSE(IF_INPUT,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(iFileName))
    END IF ! global%error

! ******************************************************************************
!   Write out information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary condition information:'
      WRITE(STDOUT,'(A,5X,A,2X,A,2X,A,2X,A,12X,A,2X,A,3X,A,15X,A)') &
            SOLVER_NAME,'Local','Global','Type','Name','Order','Constr', &
            'Coupling','Motion'

      DO iPatch = 1,pRegion%grid%nPatches
        pPatch => pRegion%patches(iPatch)

        IF ( pPatch%bcCoupled == BC_NOT_COUPLED ) THEN
          bcCoupledString = 'Not coupled'
        ELSE IF ( pPatch%bcCoupled == BC_NOT_BURNING ) THEN
          bcCoupledString = 'Coupled, not burning'
        ELSE IF ( pPatch%bcCoupled == BC_BURNING ) THEN
          bcCoupledString = 'Coupled, burning'
        ELSE ! Defensive programming
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! pPatch%bcCoupled

        WRITE(STDOUT,'(A,4X,I4,4X,I4,4X,I3,3X,A,2X,I2,6X,I2,5X,A,5X,I2)') &
              SOLVER_NAME,iPatch,pPatch%iPatchGlobal,pPatch%bcType, &
              pPatch%bcName(1:15),pPatch%spaceOrder,pPatch%cReconst, &
              bcCoupledString(1:20),pPatch%movePatchDir
      END DO ! iPatch
    END IF ! global%verbLevel

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel >= VERBOSE_MED ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Reading Rocflu boundary condition file done.'
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadBcInputFile






! ******************************************************************************
!
! Purpose: Read in user input related to boundary conditions (done on all
!   processors).
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadBcInputFileWrapper(pRegion)

#ifdef SPEC
    USE ModInterfacesSpecies, ONLY: SPEC_RFLU_ReadBcInputFile
#endif

    IMPLICIT NONE

! *****************************************************************************
!   Definitions and declarations
! *****************************************************************************

! =============================================================================
!   Arguments
! =============================================================================

    TYPE(t_region), POINTER :: pRegion

! =============================================================================
!   Locals
! =============================================================================

    TYPE(t_global), POINTER :: global

! *****************************************************************************
!   Start
! *****************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadBcInputFileWrapper',&
  'RFLU_ModReadBcInputFile.F90')

! *****************************************************************************
!   Call routines to read boundary condition information
! *****************************************************************************

! =============================================================================
!   Mixture
! =============================================================================

    CALL RFLU_ReadBcInputFile(pRegion)

! =============================================================================
!   Physical modules
! =============================================================================

#ifdef SPEC
    IF ( global%specUsed .EQV. .TRUE. ) THEN
      CALL SPEC_RFLU_ReadBcInputFile(pRegion)
    END IF ! global%specUsed
#endif

! *****************************************************************************
!   End
! *****************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadBcInputFileWrapper






! ******************************************************************************
!
! Purpose: Read in user input related to no slip-wall boundary condition with 
!   imposed heat flux.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadBcNoSlipWallHeatSect(pRegion)

    USE ModInterfaces, ONLY: ReadPatchSection

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER, PARAMETER :: NVALS_MAX = 10
    LOGICAL, DIMENSION(NVALS_MAX) :: defined
    CHARACTER(10) :: keys(NVALS_MAX)
    CHARACTER(256) :: iFileName
    CHARACTER(CHRLEN) :: bcName
    INTEGER :: checkSum,distrib,errorFlag,i,iData,ifl,iPatch,nVals,iPatchBeg, &
               iPatchEnd,nBFacesTot
    REAL(RFREAL), DIMENSION(NVALS_MAX) :: vals
    TYPE(t_grid) :: grid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadBcNoSlipWallHeatSect',&
  'RFLU_ModReadBcInputFile.F90')

! ******************************************************************************
!   Specify keywords and search for them
! ******************************************************************************

    nVals = NVALS_MAX

    keys(1) = 'HFLUX'
    keys(2) = 'COUPLED'
    keys(3) = 'MVPATCH'
    keys(4) = 'SMGRID'
    keys(5) = 'MOVEDIR'
    keys(6) = 'STATS'
    keys(7) = 'CRECONST'
    keys(8) = 'KIND'
    keys(9) = 'THRUSTFLAG'
    keys(10)= 'ORDER'

    CALL ReadPatchSection(global,IF_INPUT,nVals,keys,vals,iPatchBeg,iPatchEnd, &
                          distrib,iFileName,bcName,defined)

! ******************************************************************************
!   Check if specified number of patches exceeds available ones
! ******************************************************************************

    IF ( iPatchEnd > global%nPatches ) THEN
      CALL ErrorStop(global,ERR_PATCH_RANGE,__LINE__)
    END IF ! iPatchEnd

! ******************************************************************************
!   Get switches and check that all necessary values defined
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN
        pPatch%bcType = BC_NOSLIPWALL_HFLUX
        pPatch%bcName = bcName

! TEMPORARY - No longer used, keep for backward compatibility 
        pPatch%movePatch  = .FALSE. 
        pPatch%smoothGrid = .FALSE. 
! END TEMPORARY 

        pPatch%mixt%nData     = 1
        pPatch%mixt%nSwitches = 0
        
        pPatch%mixt%distrib = BCDAT_CONSTANT ! TEMPORARY restriction

! ------------------------------------------------------------------------------
!       Check that necessary variables defined
! ------------------------------------------------------------------------------

        IF ( pPatch%mixt%distrib == BCDAT_CONSTANT ) THEN
          checkSum = 0

          DO i = 0,pPatch%mixt%nData-1
            IF ( defined(1+i) .EQV. .TRUE. ) THEN
              checkSum = checkSum + 1
            END IF ! defined
          END DO ! i

          IF ( checkSum /= pPatch%mixt%nData ) THEN
            CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)
          END IF ! checkSum
        END IF ! pPatch
        
! ------------------------------------------------------------------------------
!       initialize Boundary Condition kind
! ------------------------------------------------------------------------------

        IF ( defined(8) .EQV. .TRUE. ) THEN
          IF ( (NINT(vals(8)) >= BC_KIND_MIN) .AND. & 
               (NINT(vals(8)) <= BC_KIND_MAX) ) THEN
            pPatch%bcKind = NINT(vals(8))
          ELSE
! TEMPORARY : issue a warning here ...
            pPatch%bcKind = BC_KIND_SIMPLE ! Initialize with Default BC Kind
          END IF ! checking range of vals(8)
        ELSE
          pPatch%bcKind = BC_KIND_SIMPLE ! Default BC Kind
        END IF ! defined(8)

! ------------------------------------------------------------------------------
!       initialize patch thrustFlag
! ------------------------------------------------------------------------------

        IF ( defined(9) .EQV. .TRUE. ) THEN
          IF ( (vals(9) > 0.5_RFREAL) .AND. & 
               (vals(9) < 1.5_RFREAL) ) THEN
            pPatch%thrustFlag = .TRUE.
          ELSE
            pPatch%thrustFlag = .FALSE.
          END IF ! 
        ELSE
          pPatch%thrustFlag = .FALSE.
        END IF ! defined(9)

! ------------------------------------------------------------------------------
!       Set patch spatial order
! ------------------------------------------------------------------------------

        IF ( defined(10) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(10)) == 2 ) THEN
            pPatch%spaceOrder = 2
          ELSE
            pPatch%spaceOrder = 1
          END IF ! NINT(vals(10))
        ELSE
          pPatch%spaceOrder = pRegion%mixtInput%spaceOrderBFaces
        END IF ! defined

! ------------------------------------------------------------------------------
!       Check whether coupling flag defined
! ------------------------------------------------------------------------------

#ifdef GENX
        IF ( defined(2) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(2)) == BC_NOT_BURNING ) THEN
            pPatch%bcCoupled       = NINT(vals(2))
            pPatch%mixt%distrib = BCDAT_DISTRIB ! MUST have distribution
          ELSE IF ( NINT(vals(2)) == BC_BURNING ) THEN 
            global%warnCounter = global%warnCounter + 1

            WRITE(STDOUT,'(A,3X,A,1X,A,I3,A)') SOLVER_NAME,'*** WARNING ***', &
                  'Invalid coupling input for patch ',iPatch, &
                  '. Overriding user input.'
                  
            pPatch%bcCoupled       = BC_NOT_BURNING
            pPatch%mixt%distrib = BCDAT_DISTRIB ! MUST have distribution	    	    
          ELSE IF ( NINT(vals(2)) == BC_NOT_COUPLED ) THEN
            pPatch%bcCoupled       = NINT(vals(2))
            pPatch%mixt%distrib = distrib	    	  
	  ELSE 
            global%warnCounter = global%warnCounter + 1

            WRITE(STDOUT,'(A,3X,A,1X,A,I3,A)') SOLVER_NAME,'*** WARNING ***', &
                  'Invalid coupling input for patch ',iPatch, &
                  '. Overriding user input.'

            pPatch%bcCoupled       = BC_NOT_COUPLED
            pPatch%mixt%distrib = distrib
          END IF ! NINT(vals)
        ELSE
          pPatch%bcCoupled       = BC_NOT_COUPLED
          pPatch%mixt%distrib = distrib
        END IF ! defined
#else
        IF ( defined(2) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(2)) /= BC_NOT_COUPLED ) THEN
            global%warnCounter = global%warnCounter + 1
          
            WRITE(STDOUT,'(A,3X,A,1X,A,I3,A)') SOLVER_NAME,'*** WARNING ***', &
                  'Invalid coupling input for patch ',iPatch, &
                  '. Overriding user input.'                      

            pPatch%bcCoupled = BC_NOT_COUPLED
          ELSE
            pPatch%bcCoupled = BC_NOT_COUPLED
	  END IF ! NINT(vals(2))
	  
          pPatch%mixt%distrib = distrib        
        ELSE
          pPatch%bcCoupled       = BC_NOT_COUPLED
          pPatch%mixt%distrib = BCDAT_CONSTANT
	END IF ! defined(2)
#endif

! ------------------------------------------------------------------------------
!       Set patch motion variable
! ------------------------------------------------------------------------------

        IF ( defined(5) .EQV. .TRUE. ) THEN 
          pPatch%movePatchDir = vals(5)
        ELSE 
          pPatch%movePatchDir = MOVEPATCH_DIR_NONE
        END IF ! defined
        
! ------------------------------------------------------------------------------
!       Set patch statistics plotting flag
! ------------------------------------------------------------------------------

        IF ( defined(6) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(6)) == 1 ) THEN  
            pPatch%plotStatsFlag = .TRUE.
          ELSE
            pPatch%plotStatsFlag = .FALSE.
          END IF ! NINT
        ELSE 
          pPatch%plotStatsFlag = .FALSE.
        END IF ! defined    
        
! ------------------------------------------------------------------------------
!       Set constraint variable
! ------------------------------------------------------------------------------

        IF ( defined(7) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(7)) == 1 ) THEN  
            pPatch%cReconst = CONSTR_WEIGHTED
          ELSE 
            pPatch%cReconst = CONSTR_NONE
          END IF ! NINT(vals(7))  
        ELSE 
          pPatch%cReconst = CONSTR_NONE
        END IF ! defined                                
      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! ******************************************************************************
!   Copy values/distribution to variables 
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN

! ------------------------------------------------------------------------------
!       Distribution from file: Allocate and initialize, actual values read in
!       at later stage
! ------------------------------------------------------------------------------

        IF ( pPatch%mixt%distrib == BCDAT_DISTRIB ) THEN
	  nBFacesTot = pPatch%nBTrisTot + pPatch%nBQuadsTot	
	
          ALLOCATE(pPatch%mixt%vals(pPatch%mixt%nData,nBFacesTot), &
                   STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%vals')
          END IF ! global

          DO ifl = 1,nBFacesTot
            DO iData = 1,pPatch%mixt%nData
              pPatch%mixt%vals(iData,ifl) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
            END DO ! iData
          END DO ! ifl

! ------------------------------------------------------------------------------
!       Constant value
! ------------------------------------------------------------------------------

        ELSE
          ALLOCATE(pPatch%mixt%vals(pPatch%mixt%nData,0:1), &
                   STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%vals')
          END IF ! global

          pPatch%mixt%vals(BCDAT_NOSLIP_Q,0:1) = vals(1)
        END IF ! pPatch%mixt%distrib
      END IF ! pPatch%mixt%switches
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadBcNoSlipWallHeatSect








! ******************************************************************************
!
! Purpose: Read in user input related to no slip-wall boundary condition with 
!   imposed temperature.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadBcNoSlipWallTempSect(pRegion)

    USE ModInterfaces, ONLY: ReadPatchSection

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER, PARAMETER :: NVALS_MAX = 10
    LOGICAL, DIMENSION(NVALS_MAX) :: defined
    CHARACTER(10) :: keys(NVALS_MAX)
    CHARACTER(256) :: iFileName
    CHARACTER(CHRLEN) :: bcName
    INTEGER :: checkSum,distrib,errorFlag,i,iData,ifl,iPatch,nBFacesTot, & 
               nVals,iPatchBeg,iPatchEnd
    REAL(RFREAL), DIMENSION(NVALS_MAX) :: vals
    TYPE(t_grid) :: grid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadBcNoSlipWallTempSect',&
  'RFLU_ModReadBcInputFile.F90')

! ******************************************************************************
!   Specify keywords and search for them
! ******************************************************************************

    nVals = NVALS_MAX

    keys(1) = 'TEMP'
    keys(2) = 'COUPLED'
    keys(3) = 'MVPATCH'
    keys(4) = 'SMGRID'
    keys(5) = 'MOVEDIR'
    keys(6) = 'STATS'
    keys(7) = 'CRECONST'
    keys(8) = 'KIND'
    keys(9) = 'THRUSTFLAG'
    keys(10)= 'ORDER'

    CALL ReadPatchSection(global,IF_INPUT,nVals,keys,vals,iPatchBeg,iPatchEnd, &
                          distrib,iFileName,bcName,defined)

! ******************************************************************************
!   Check if specified number of patches exceeds available ones
! ******************************************************************************

    IF ( iPatchEnd > global%nPatches ) THEN
      CALL ErrorStop(global,ERR_PATCH_RANGE,__LINE__)
    END IF ! iPatchEnd

! ******************************************************************************
!   Get switches and check that all necessary values defined
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN
        pPatch%bcType = BC_NOSLIPWALL_TEMP
        pPatch%bcName = bcName

! TEMPORARY - No longer used, keep for backward compatibility 
        pPatch%movePatch  = .FALSE. 
        pPatch%smoothGrid = .FALSE. 
! END TEMPORARY 

        pPatch%mixt%nData     = 1
        pPatch%mixt%nSwitches = 0
        
! ------------------------------------------------------------------------------
!       initialize Boundary Condition kind
! ------------------------------------------------------------------------------

        IF ( defined(8) .EQV. .TRUE. ) THEN
          IF ( (NINT(vals(8)) >= BC_KIND_MIN) .AND. & 
               (NINT(vals(8)) <= BC_KIND_MAX) ) THEN
            pPatch%bcKind = NINT(vals(8))
          ELSE
! TEMPORARY : issue a warning here ...
            pPatch%bcKind = BC_KIND_SIMPLE ! Initialize with Default BC Kind
          END IF ! checking range of vals(8)
        ELSE
          pPatch%bcKind = BC_KIND_SIMPLE ! Default BC Kind
        END IF ! defined(8)

! ------------------------------------------------------------------------------
!       initialize patch thrustFlag
! ------------------------------------------------------------------------------

        IF ( defined(9) .EQV. .TRUE. ) THEN
          IF ( (vals(9) > 0.5_RFREAL) .AND. &
               (vals(9) < 1.5_RFREAL) ) THEN
            pPatch%thrustFlag = .TRUE.
          ELSE
            pPatch%thrustFlag = .FALSE.
          END IF !
        ELSE
          pPatch%thrustFlag = .FALSE.
        END IF ! defined(9)

! ------------------------------------------------------------------------------
!       Set patch spatial order
! ------------------------------------------------------------------------------

        IF ( defined(10) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(10)) == 2 ) THEN
            pPatch%spaceOrder = 2
          ELSE
            pPatch%spaceOrder = 1
          END IF ! NINT(vals(10))
        ELSE
          pPatch%spaceOrder = pRegion%mixtInput%spaceOrderBFaces
        END IF ! defined

! ------------------------------------------------------------------------------
!       Check whether coupling flag defined
! ------------------------------------------------------------------------------

#ifdef GENX
        IF ( defined(2) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(2)) == BC_NOT_BURNING ) THEN
            pPatch%bcCoupled       = NINT(vals(2))
            pPatch%mixt%distrib = BCDAT_DISTRIB ! MUST have distribution
          ELSE IF ( NINT(vals(2)) == BC_BURNING ) THEN 
            global%warnCounter = global%warnCounter + 1

            WRITE(STDOUT,'(A,3X,A,1X,A,I3,A)') SOLVER_NAME,'*** WARNING ***', &
                  'Invalid coupling input for patch ',iPatch, &
                  '. Overriding user input.'
                  
            pPatch%bcCoupled       = BC_NOT_BURNING
            pPatch%mixt%distrib = BCDAT_DISTRIB ! MUST have distribution	    	    
          ELSE IF ( NINT(vals(2)) == BC_NOT_COUPLED ) THEN
            pPatch%bcCoupled       = NINT(vals(2))
            pPatch%mixt%distrib = distrib	    	  
	  ELSE 
            global%warnCounter = global%warnCounter + 1

            WRITE(STDOUT,'(A,3X,A,1X,A,I3,A)') SOLVER_NAME,'*** WARNING ***', &
                  'Invalid coupling input for patch ',iPatch, &
                  '. Overriding user input.'

            pPatch%bcCoupled       = BC_NOT_COUPLED
            pPatch%mixt%distrib = distrib
          END IF ! NINT(vals)
        ELSE
          pPatch%bcCoupled       = BC_NOT_COUPLED
          pPatch%mixt%distrib = distrib
        END IF ! defined
#else
        IF ( defined(2) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(2)) /= BC_NOT_COUPLED ) THEN
            global%warnCounter = global%warnCounter + 1
          
            WRITE(STDOUT,'(A,3X,A,1X,A,I3,A)') SOLVER_NAME,'*** WARNING ***', &
                  'Invalid coupling input for patch ',iPatch, &
                  '. Overriding user input.'                      

            pPatch%bcCoupled = BC_NOT_COUPLED
          ELSE
            pPatch%bcCoupled = BC_NOT_COUPLED
	  END IF ! NINT(vals(2))
	  
          pPatch%mixt%distrib = distrib        
        ELSE
          pPatch%bcCoupled       = BC_NOT_COUPLED
          pPatch%mixt%distrib = BCDAT_CONSTANT
	END IF ! defined(2)
#endif

! ------------------------------------------------------------------------------
!       Check that necessary variables defined
! ------------------------------------------------------------------------------

        IF ( pPatch%mixt%distrib == BCDAT_CONSTANT ) THEN
          checkSum = 0

          DO i = 0,pPatch%mixt%nData-1
            IF ( defined(1+i) .EQV. .TRUE. ) THEN
              checkSum = checkSum + 1
            END IF ! defined
          END DO ! i

          IF ( checkSum /= pPatch%mixt%nData ) THEN
            CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)
          END IF ! checkSum
        END IF ! pPatch

! ------------------------------------------------------------------------------
!       Set patch motion variable
! ------------------------------------------------------------------------------

        IF ( defined(5) .EQV. .TRUE. ) THEN 
          pPatch%movePatchDir = vals(5)
        ELSE 
          pPatch%movePatchDir = MOVEPATCH_DIR_NONE
        END IF ! defined
        
! ------------------------------------------------------------------------------
!       Set patch statistics plotting flag
! ------------------------------------------------------------------------------

        IF ( defined(6) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(6)) == 1 ) THEN  
            pPatch%plotStatsFlag = .TRUE.
          ELSE
            pPatch%plotStatsFlag = .FALSE.
          END IF ! NINT
        ELSE 
          pPatch%plotStatsFlag = .FALSE.
        END IF ! defined   
        
! ------------------------------------------------------------------------------
!       Set constraint variable
! ------------------------------------------------------------------------------

        IF ( defined(7) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(7)) == 1 ) THEN  
            pPatch%cReconst = CONSTR_WEIGHTED
          ELSE 
            pPatch%cReconst = CONSTR_NONE
          END IF ! NINT(vals(7))  
        ELSE 
          pPatch%cReconst = CONSTR_NONE
        END IF ! defined                                   
      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! ******************************************************************************
!   Copy values/distribution to variables 
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN

! ------------------------------------------------------------------------------
!       Distribution from file: Allocate and initialize, actual values read in
!       at later stage
! ------------------------------------------------------------------------------

        IF ( pPatch%mixt%distrib == BCDAT_DISTRIB ) THEN
	  nBFacesTot = pPatch%nBTrisTot + pPatch%nBQuadsTot
	
          ALLOCATE(pPatch%mixt%vals(pPatch%mixt%nData,nBFacesTot), &
                   STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%vals')
          END IF ! global

! ------- If not coupled, get boundary data from file --------------------------

          IF ( pPatch%bcCoupled == BC_NOT_COUPLED ) THEN
            DO ifl = 1,nBFacesTot
              DO iData = 1,pPatch%mixt%nData
                pPatch%mixt%vals(iData,ifl) = &
                  REAL(CRAZY_VALUE_INT,KIND=RFREAL)
              END DO ! iData
            END DO ! ifl

! ------- If coupled, get boundary data from Roccom ----------------------------

          ELSE ! Initialize - important for GENX runs from scratch
            DO ifl = 1,nBFacesTot
              DO iData = 1,pPatch%mixt%nData
                pPatch%mixt%vals(iData,ifl) = 0.0_RFREAL
              END DO ! iData
            END DO ! ifl
          END IF ! pPatch%bcCoupled

! ------------------------------------------------------------------------------
!       Constant value
! ------------------------------------------------------------------------------

        ELSE
          ALLOCATE(pPatch%mixt%vals(pPatch%mixt%nData,0:1), &
                   STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%vals')
          END IF ! global

          pPatch%mixt%vals(BCDAT_NOSLIP_T,0:1) = vals(1)
        END IF ! pPatch%mixt%distrib
      END IF ! pPatch%mixt%switches
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadBcNoSlipWallTempSect









! ******************************************************************************
!
! Purpose: Read in user input related to outflow boundary condition.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadBcOutflowSection(pRegion)

    USE ModInterfaces, ONLY: ReadPatchSection

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER, PARAMETER :: NVALS_MAX = 12

    LOGICAL, DIMENSION(NVALS_MAX) :: defined
    CHARACTER(10) :: keys(NVALS_MAX)
    CHARACTER(256) :: iFileName
    CHARACTER(CHRLEN) :: bcName
    INTEGER :: checkSum,distrib,errorFlag,i,iData,ifl,iPatch,iReg,iPatchBeg, &
               iPatchEnd,nVals
    REAL(RFREAL), DIMENSION(NVALS_MAX) :: vals
    TYPE(t_grid) :: grid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadBcOutflowSection',&
  'RFLU_ModReadBcInputFile.F90')

! ******************************************************************************
!   Specify keywords and search for them
! ******************************************************************************

    nVals = NVALS_MAX

    keys(1) = 'TYPE'
    keys(2) = 'PRESS'
    keys(3) = 'MVPATCH'
    keys(4) = 'SMGRID'
    keys(5) = 'ORDER'
    keys(6) = 'STATS'
    keys(7) = 'MOVEDIR'
    keys(8) = 'COUPLED'
    keys(9) = 'KIND'
    keys(10) = 'REFLECT'
    keys(11) = 'NSCBCK'
    keys(12) = 'THRUSTFLAG'

    CALL ReadPatchSection(global,IF_INPUT,nVals,keys,vals,iPatchBeg,iPatchEnd, &
                          distrib,iFileName,bcName,defined )

! ******************************************************************************
!   Check if specified number of patches exceeds available ones
! ******************************************************************************

    IF ( iPatchEnd > global%nPatches ) THEN
      CALL ErrorStop(global,ERR_PATCH_RANGE,__LINE__)
    END IF ! iPatchEnd

! ******************************************************************************
!   Get switches and check that all necessary values defined
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN
        pPatch%bcType = BC_OUTFLOW
        pPatch%bcName = bcName
        
!        pPatch%bcCoupled    = BC_NOT_COUPLED
!        pPatch%movePatchDir = MOVEPATCH_DIR_NONE
        
        pPatch%cReconst      = CONSTR_NONE
        pPatch%plotStatsFlag = .FALSE.

! TEMPORARY - No longer used, keep for backward compatibility 
        pPatch%movePatch  = .FALSE. 
        pPatch%smoothGrid = .FALSE. 
! END TEMPORARY 

! ------------------------------------------------------------------------------
!       Set switches
! ------------------------------------------------------------------------------

        pPatch%mixt%nSwitches = 1
        pPatch%mixt%distrib   = distrib

        ALLOCATE(pPatch%mixt%switches(pPatch%mixt%nSwitches), &
                 STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%switches')
        END IF ! global

! ------------------------------------------------------------------------------
!       initialize Boundary Condition kind
! ------------------------------------------------------------------------------

        IF ( defined(9) .EQV. .TRUE. ) THEN
          IF ( (NINT(vals(9)) >= BC_KIND_MIN) .AND. & 
               (NINT(vals(9)) <= BC_KIND_MAX) ) THEN
            pPatch%bcKind = NINT(vals(9))
          ELSE
! TEMPORARY : issue a warning here ...
            pPatch%bcKind = BC_KIND_SIMPLE ! Initialize with Default BC Kind
          END IF ! checking range of vals(9)
        ELSE
          pPatch%bcKind = BC_KIND_SIMPLE ! Default BC Kind
        END IF ! defined(9)

! ------------------------------------------------------------------------------
!       initialize if far field is reflecting or non reflecting
! ------------------------------------------------------------------------------

        IF ( defined(10) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(10)) == BC_REFLECTING ) THEN 
            pPatch%reflect = BC_REFLECTING
          ELSE
            pPatch%reflect = BC_NONREFLECTING
          END IF ! 
        ELSE
          pPatch%reflect = BC_REFLECTING
        END IF ! defined(10)

        IF ( defined(11) .EQV. .TRUE. ) THEN
          pPatch%nscbcK = vals(11) 
        END IF ! defined(11)

! ------------------------------------------------------------------------------
!       initialize patch thrustFlag
! ------------------------------------------------------------------------------

        IF ( defined(12) .EQV. .TRUE. ) THEN
          IF ( (vals(12) > 0.5_RFREAL) .AND. &
               (vals(12) < 1.5_RFREAL) ) THEN
            pPatch%thrustFlag = .TRUE.
          ELSE
            pPatch%thrustFlag = .FALSE.
          END IF !
        ELSE
          pPatch%thrustFlag = .FALSE.
        END IF ! defined(12)

! ------------------------------------------------------------------------------
!       Check if switch defined
! ------------------------------------------------------------------------------

        IF ( defined(1) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(1)) == 0 ) THEN
            pPatch%mixt%switches(BCSWI_OUTFLOW_TYPE) = BCOPT_SUPERSONIC
          ELSE IF ( NINT(vals(1)) == 1 ) THEN
            pPatch%mixt%switches(BCSWI_OUTFLOW_TYPE) = BCOPT_SUBSONIC
          ELSE IF ( NINT(vals(1)) == 2 ) THEN
            pPatch%mixt%switches(BCSWI_OUTFLOW_TYPE) = BCOPT_MIXED
          ELSE
            CALL ErrorStop(global,ERR_VAL_BCSWITCH,__LINE__,'(outflow type).')
          END IF ! NINT
        ELSE
          CALL ErrorStop(global,ERR_NO_BCSWITCH,__LINE__,'(outflow type).')
        END IF ! defined

! ------------------------------------------------------------------------------
!       Check whether appropriate values specified
! ------------------------------------------------------------------------------

        IF ( pPatch%mixt%switches(BCSWI_OUTFLOW_TYPE) /= &
             BCOPT_SUPERSONIC ) THEN
          pPatch%mixt%nData = 1

          IF ( pPatch%mixt%distrib == BCDAT_CONSTANT ) THEN
            IF ( defined(2) .EQV. .FALSE. ) THEN
              CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)
            END IF ! defined
          END IF ! patch
        ELSE
          pPatch%mixt%nData = 0
        END IF ! pPatch
        
! ------------------------------------------------------------------------------
!       Set patch spatial order
! ------------------------------------------------------------------------------

        IF ( defined(5) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(5)) == 2 ) THEN  
            pPatch%spaceOrder = 2
          ELSE 
            pPatch%spaceOrder = 1
          END IF ! NINT(vals(5))  
        ELSE 
          pPatch%spaceOrder = pRegion%mixtInput%spaceOrderBFaces
        END IF ! defined                    

! ------------------------------------------------------------------------------
!       Set patch statistics plotting flag
! ------------------------------------------------------------------------------

        IF ( defined(6) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(6)) == 1 ) THEN
            pPatch%plotStatsFlag = .TRUE.
          ELSE
            pPatch%plotStatsFlag = .FALSE.
          END IF ! NINT
        ELSE
          pPatch%plotStatsFlag = .FALSE.
        END IF ! defined

! ------------------------------------------------------------------------------
!       Set patch motion variable
! ------------------------------------------------------------------------------

        IF ( defined(7) .EQV. .TRUE. ) THEN
          pPatch%movePatchDir = vals(7)
        ELSE
          pPatch%movePatchDir = MOVEPATCH_DIR_NONE
        END IF ! defined

! ------------------------------------------------------------------------------
!       Set coupling variable
! ------------------------------------------------------------------------------

        IF ( defined(8) .EQV. .TRUE. ) THEN
          pPatch%bcCoupled = vals(8)
        ELSE
          pPatch%bcCoupled = BC_NOT_COUPLED
        END IF ! defined

      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! ******************************************************************************
!   Copy values/distribution to variables for non-adiabatic walls
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN

        IF ( pPatch%mixt%switches(BCSWI_OUTFLOW_TYPE) /= &
             BCOPT_SUPERSONIC ) THEN

! ------------------------------------------------------------------------------
!         Distribution from file: Allocate and initialize, actual values read in
!         at later stage
! ------------------------------------------------------------------------------

          IF ( pPatch%mixt%distrib == BCDAT_DISTRIB ) THEN
            ALLOCATE(pPatch%mixt%vals(pPatch%mixt%nData,pPatch%nBFaces), &
                     STAT=errorFlag)
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) THEN
              CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
            END IF ! global

            DO ifl = 1,pPatch%nBFaces
              DO iData = 1,pPatch%mixt%nData
                pPatch%mixt%vals(iData,ifl) = &
                  REAL(CRAZY_VALUE_INT,KIND=RFREAL)
              END DO ! iData
            END DO ! ifl

! ------------------------------------------------------------------------------
!         Constant value
! ------------------------------------------------------------------------------

          ELSE
            ALLOCATE(pPatch%mixt%vals(pPatch%mixt%nData,0:1), &
                     STAT=errorFlag)
            global%error = errorFlag
            IF ( global%error /= ERR_NONE ) THEN
              CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%mixt%vals')
            END IF ! global

            pPatch%mixt%vals(BCDAT_OUTFLOW_PRESS,0:1) = vals(2)
          END IF ! pPatch%mixt%distrib
        ELSE
          NULLIFY(pPatch%mixt%vals)
        END IF ! pPatch%mixt%switches

      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadBcOutflowSection









! ******************************************************************************
!
! Purpose: Read in user input related to periodic boundary condition.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadBcPeriodicSection(pRegion)

    USE ModInterfaces, ONLY: ReadPatchSection

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER, PARAMETER :: NVALS_MAX = 5

    LOGICAL, DIMENSION(NVALS_MAX) :: defined
    CHARACTER(10) :: keys(NVALS_MAX)
    CHARACTER(256) :: iFileName
    CHARACTER(CHRLEN) :: bcName
    INTEGER :: checkSum,distrib,errorFlag,i,iData,ifl,iPatch,iReg,iPatchBeg, &
               iPatchEnd,nVals
    REAL(RFREAL), DIMENSION(NVALS_MAX) :: vals
    TYPE(t_grid) :: grid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadBcPeriodicSection',&
  'RFLU_ModReadBcInputFile.F90')

! ******************************************************************************
!   Specify keywords and search for them
! ******************************************************************************

    nVals = NVALS_MAX

    keys(1) = 'RELPATCH'
    keys(2) = 'ANGLE'
    keys(3) = 'AXIS'
    keys(4) = 'MVPATCH'
    keys(5) = 'SMGRID'

    CALL ReadPatchSection(global,IF_INPUT,nVals,keys,vals,iPatchBeg,iPatchEnd, &
                          distrib,iFileName,bcName,defined )

! ******************************************************************************
!   Check if specified number of patches exceeds available ones
! ******************************************************************************

    IF ( iPatchEnd > global%nPatches ) THEN
      CALL ErrorStop(global,ERR_PATCH_RANGE,__LINE__)
    END IF ! iPatchEnd

! ******************************************************************************
!   Get switches and check that all necessary values defined
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN
        pPatch%bcType    = BC_PERIODIC
        pPatch%bcName    = bcName
        pPatch%bcCoupled = BC_NOT_COUPLED

        pPatch%mixt%nData     = 0
        pPatch%mixt%nSwitches = 0
        pPatch%mixt%distrib   = BCDAT_CONSTANT
       
        pPatch%bcKind        = BC_KIND_SIMPLE ! Value immaterial 
        pPatch%thrustFlag    = .FALSE. ! Value immaterial

        pPatch%spaceOrder = pRegion%mixtInput%spaceOrderBFaces ! Value immaterial
        pPatch%cReconst   = CONSTR_NONE
        
! TEMPORARY - No longer used, keep for backward compatibility 
        pPatch%movePatch  = .FALSE. 
        pPatch%smoothGrid = .FALSE. 
! END TEMPORARY         
        
        IF ( defined(1) .EQV. .FALSE. ) THEN
          CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)
        ELSE
          pPatch%iPatchRelated = NINT(vals(1))
        END IF ! defined
        
        IF ( defined(2) .EQV. .FALSE. ) THEN
          CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)
        ELSE
          pPatch%angleRelated = vals(2)*global%deg2rad
        END IF ! defined   
        
        IF ( defined(3) .EQV. .FALSE. ) THEN
          CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)
        ELSE
          pPatch%axisRelated = NINT(vals(3))
        END IF ! defined              
      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadBcPeriodicSection








! ******************************************************************************
!
! Purpose: Read in user input related to slip-wall boundary condition.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadBcSlipWallSection(pRegion)

    USE ModInterfaces, ONLY: ReadPatchSection

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER, PARAMETER :: NVALS_MAX = 9

    LOGICAL, DIMENSION(NVALS_MAX) :: defined
    CHARACTER(10) :: keys(NVALS_MAX)
    CHARACTER(256) :: iFileName
    CHARACTER(CHRLEN) :: bcName
    INTEGER :: distrib,errorFlag,iPatch,nVals,iPatchBeg,iPatchEnd
    REAL(RFREAL), DIMENSION(NVALS_MAX) :: vals
    TYPE(t_grid) :: grid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadBcSlipWallSection',&
  'RFLU_ModReadBcInputFile.F90')

! ******************************************************************************
!   Specify keywords and search for them
! ******************************************************************************

    nVals = NVALS_MAX

    keys(1) = 'COUPLED'
    keys(2) = 'MVPATCH'
    keys(3) = 'SMGRID'
    keys(4) = 'MOVEDIR'
    keys(5) = 'STATS'
    keys(6) = 'CRECONST'
    keys(7) = 'KIND'
    keys(8) = 'THRUSTFLAG'
    keys(9) = 'ORDER'

    CALL ReadPatchSection(global,IF_INPUT,nVals,keys,vals,iPatchBeg,iPatchEnd, &
                          distrib,iFileName,bcName,defined)

! ******************************************************************************
!   Check if specified number of patches exceeds available ones
! ******************************************************************************

    IF ( iPatchEnd > global%nPatches ) THEN
      CALL ErrorStop(global,ERR_PATCH_RANGE,__LINE__)
    END IF ! iPatchEnd

! ******************************************************************************
!   Loop over patches and copy values/distribution to variables
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN
        pPatch%bcType = BC_SLIPWALL
        pPatch%bcName = bcName

! TEMPORARY - No longer used, keep for backward compatibility 
        pPatch%movePatch  = .FALSE. 
        pPatch%smoothGrid = .FALSE. 
! END TEMPORARY 

        pPatch%mixt%nData     = 0
        pPatch%mixt%nSwitches = 0
        pPatch%mixt%distrib   = BCDAT_CONSTANT

        NULLIFY(pPatch%mixt%switches)
        NULLIFY(pPatch%mixt%vals)

! ------------------------------------------------------------------------------
!       initialize Boundary Condition kind
! ------------------------------------------------------------------------------

        IF ( defined(7) .EQV. .TRUE. ) THEN
          IF ( (NINT(vals(7)) >= BC_KIND_MIN) .AND. & 
               (NINT(vals(7)) <= BC_KIND_MAX) ) THEN
            pPatch%bcKind = NINT(vals(7))
          ELSE
! TEMPORARY : issue a warning here ...
            pPatch%bcKind = BC_KIND_SIMPLE ! Initialize with Default BC Kind
          END IF ! checking range of vals(7)
        ELSE
          pPatch%bcKind = BC_KIND_SIMPLE ! Default BC Kind
        END IF ! defined(7)

! ------------------------------------------------------------------------------
!       initialize patch thrustFlag
! ------------------------------------------------------------------------------

        IF ( defined(8) .EQV. .TRUE. ) THEN
          IF ( (vals(8) > 0.5_RFREAL) .AND. &
               (vals(8) < 1.5_RFREAL) ) THEN
            pPatch%thrustFlag = .TRUE.
          ELSE
            pPatch%thrustFlag = .FALSE.
          END IF !
        ELSE
          pPatch%thrustFlag = .FALSE.
        END IF ! defined(8)

! ------------------------------------------------------------------------------
!       Set patch spatial order
! ------------------------------------------------------------------------------

        IF ( defined(9) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(9)) == 2 ) THEN
            pPatch%spaceOrder = 2
          ELSE
            pPatch%spaceOrder = 1
          END IF ! NINT(vals(9))
        ELSE
          pPatch%spaceOrder = pRegion%mixtInput%spaceOrderBFaces
        END IF ! defined

! ------------------------------------------------------------------------------
!       Check whether coupling flag defined
! ------------------------------------------------------------------------------

#ifdef GENX
        IF ( defined(1) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(1)) == BC_NOT_BURNING ) THEN 
            pPatch%bcCoupled = NINT(vals(1))
          ELSE IF ( NINT(vals(1)) == BC_BURNING ) THEN 
            global%warnCounter = global%warnCounter + 1

            WRITE(STDOUT,'(A,3X,A,1X,A,I3,A)') SOLVER_NAME,'*** WARNING ***', &
                  'Invalid coupling input for patch ',iPatch, &
                  '. Overriding user input.'
                  
            pPatch%bcCoupled = BC_NOT_BURNING     
          ELSE IF ( NINT(vals(1)) == BC_NOT_COUPLED ) THEN 
            pPatch%bcCoupled = NINT(vals(1))
          ELSE 
            global%warnCounter = global%warnCounter + 1

            WRITE(STDOUT,'(A,3X,A,1X,A,I3,A)') SOLVER_NAME,'*** WARNING ***', &
                  'Invalid coupling input for patch ',iPatch, &
                  '. Overriding user input.'
                  
            pPatch%bcCoupled = BC_NOT_COUPLED     
          END IF ! NINT(vals(1))
        ELSE
          pPatch%bcCoupled = BC_NOT_COUPLED
        END IF ! defined
#else
        IF ( defined(1) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(1)) /= BC_NOT_COUPLED ) THEN
            global%warnCounter = global%warnCounter + 1
          
            WRITE(STDOUT,'(A,3X,A,1X,A,I3,A)') SOLVER_NAME,'*** WARNING ***', &
                  'Invalid coupling input for patch ',iPatch, &
                  '. Overriding user input.'                      

            pPatch%bcCoupled = BC_NOT_COUPLED
          ELSE
            pPatch%bcCoupled = BC_NOT_COUPLED
          END IF ! NINT(vals(1))   
        ELSE
          pPatch%bcCoupled = BC_NOT_COUPLED
        END IF ! defined(1)
#endif

! ------------------------------------------------------------------------------
!       Set patch motion variable
! ------------------------------------------------------------------------------

        IF ( defined(4) .EQV. .TRUE. ) THEN 
          pPatch%movePatchDir = vals(4)
        ELSE 
          pPatch%movePatchDir = MOVEPATCH_DIR_NONE
        END IF ! defined
        
! ------------------------------------------------------------------------------
!       Set patch statistics plotting flag
! ------------------------------------------------------------------------------

        IF ( defined(5) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(5)) == 1 ) THEN  
            pPatch%plotStatsFlag = .TRUE.
          ELSE
            pPatch%plotStatsFlag = .FALSE.
          END IF ! NINT
        ELSE 
          pPatch%plotStatsFlag = .FALSE.
        END IF ! defined              
        
! ------------------------------------------------------------------------------
!       Set constraint variable
! ------------------------------------------------------------------------------

        IF ( defined(6) .EQV. .TRUE. ) THEN
          IF ( NINT(vals(6)) == 1 ) THEN  
            pPatch%cReconst = CONSTR_WEIGHTED
          ELSE 
            pPatch%cReconst = CONSTR_NONE
          END IF ! NINT(vals(6))  
        ELSE 
          pPatch%cReconst = CONSTR_NONE
        END IF ! defined                                
      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadBcSlipWallSection










! ******************************************************************************
!
! Purpose: Read in user input related to symmetry boundary condition.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadBcSymmetrySection(pRegion)

    USE ModInterfaces, ONLY: ReadPatchSection

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER, PARAMETER :: NVALS_MAX = 2

    LOGICAL, DIMENSION(NVALS_MAX) :: defined
    CHARACTER(10) :: keys(NVALS_MAX)
    CHARACTER(256) :: iFileName
    CHARACTER(CHRLEN) :: bcName
    INTEGER :: distrib,errorFlag,iPatch,nVals,iPatchBeg,iPatchEnd
    REAL(RFREAL), DIMENSION(NVALS_MAX) :: vals
    TYPE(t_grid) :: grid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadBcSymmetrySection',&
  'RFLU_ModReadBcInputFile.F90')

! ******************************************************************************
!   Specify keywords and search for them
! ******************************************************************************

    nVals = NVALS_MAX

    keys(1) = 'MVPATCH'
    keys(2) = 'SMGRID'

    CALL ReadPatchSection(global,IF_INPUT,nVals,keys,vals,iPatchBeg,iPatchEnd, &
                          distrib,iFileName,bcName,defined)

! ******************************************************************************
!   Check if specified number of patches exceeds available ones
! ******************************************************************************

    IF ( iPatchEnd > global%nPatches ) THEN
      CALL ErrorStop(global,ERR_PATCH_RANGE,__LINE__)
    END IF ! iPatchEnd

! ******************************************************************************
!   Loop over patches and copy values/distribution to variables
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN

        pPatch%bcType    = BC_SYMMETRY
        pPatch%bcName    = bcName
        pPatch%bcCoupled = BC_NOT_COUPLED

        pPatch%iPatchRelated = iPatch

        pPatch%mixt%nData     = 0
        pPatch%mixt%nSwitches = 0
        pPatch%mixt%distrib   = BCDAT_CONSTANT

        pPatch%bcKind        = BC_KIND_SIMPLE ! Value immaterial
        pPatch%thrustFlag    = .FALSE. ! Value immaterial

        pPatch%spaceOrder = pRegion%mixtInput%spaceOrderBFaces ! Value immaterial
        pPatch%cReconst   = CONSTR_NONE

! TEMPORARY - No longer used, keep for backward compatibility 
        pPatch%movePatch  = .FALSE. 
        pPatch%smoothGrid = .FALSE. 
! END TEMPORARY 

        ALLOCATE(pPatch%mixt%switches(pPatch%mixt%nSwitches), &
                 STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
        END IF ! global

        NULLIFY(pPatch%mixt%vals)
      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadBcSymmetrySection









! ******************************************************************************
!
! Purpose: Read in user input related to virtual boundary condition.
!
! Description: None.
!
! Input:
!   pRegion                Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadBcVirtualSection(pRegion)

    USE ModInterfaces, ONLY: ReadPatchSection

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion

! ==============================================================================
!   Local variables
! ==============================================================================

    INTEGER, PARAMETER :: NVALS_MAX = 1

    LOGICAL, DIMENSION(NVALS_MAX) :: defined
    CHARACTER(10) :: keys(NVALS_MAX)
    CHARACTER(256) :: iFileName
    CHARACTER(CHRLEN) :: bcName
    INTEGER :: distrib,errorFlag,iPatch,nVals,iPatchBeg,iPatchEnd
    REAL(RFREAL), DIMENSION(NVALS_MAX) :: vals
    TYPE(t_grid) :: grid
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global

! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadBcVirtualSection',&
  'RFLU_ModReadBcInputFile.F90')

! ******************************************************************************
!   Specify keywords and search for them
! ******************************************************************************

    nVals = NVALS_MAX

    keys(1) = 'SMGRID'

    CALL ReadPatchSection(global,IF_INPUT,nVals,keys,vals,iPatchBeg,iPatchEnd, &
                          distrib,iFileName,bcName,defined)

! ******************************************************************************
!   Check if specified number of patches exceeds available ones
! ******************************************************************************

    IF ( iPatchEnd > global%nPatches ) THEN
      CALL ErrorStop(global,ERR_PATCH_RANGE,__LINE__)
    END IF ! iPatchEnd

! ******************************************************************************
!   Loop over patches and copy values/distribution to variables
! ******************************************************************************

    DO iPatch = 1,pRegion%grid%nPatches
      pPatch => pRegion%patches(iPatch)

! ==============================================================================
!     Check whether this global patch exists in this region
! ==============================================================================

      IF ( pPatch%iPatchGlobal >= iPatchBeg .AND. &
           pPatch%iPatchGlobal <= iPatchEnd ) THEN
        pPatch%bcType = BC_VIRTUAL
        pPatch%bcName = bcName

        pPatch%bcCoupled    = BC_NOT_COUPLED 
        pPatch%movePatchDir = MOVEPATCH_DIR_NONE
	
        pPatch%mixt%nData     = 0
        pPatch%mixt%nSwitches = 0
        pPatch%mixt%distrib   = BCDAT_CONSTANT

        pPatch%bcKind        = BC_KIND_SIMPLE ! Value immaterial
        pPatch%thrustFlag    = .FALSE. ! Value immaterial

        pPatch%spaceOrder    = pRegion%mixtInput%spaceOrderBFaces ! Value immaterial
        pPatch%cReconst      = CONSTR_NONE
        pPatch%plotStatsFlag = .FALSE.

! TEMPORARY - No longer used, keep for backward compatibility 
        pPatch%movePatch  = .FALSE. 
        pPatch%smoothGrid = .FALSE. 
! END TEMPORARY 

        NULLIFY(pPatch%mixt%switches)
        NULLIFY(pPatch%mixt%vals)
      END IF ! pPatch%iPatchGlobal
    END DO ! iPatch

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_ReadBcVirtualSection







! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLU_ModReadBcInputFile


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModReadBcInputFile.F90,v $
! Revision 1.27  2008/12/06 08:44:23  mtcampbe
! Updated license.
!
! Revision 1.26  2008/11/19 22:17:34  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.25  2006/10/20 21:19:05  mparmar
! added reading of THRUSTFLAG and cleaned up reading of NSCBC related keywords
!
! Revision 1.24  2006/08/19 15:39:12  mparmar
! Renamed patch variables, added reading of KIND,REFLECT,NSCBCK
!
! Revision 1.23  2006/08/10 17:19:11  rfiedler
! Corrected max array size for outflow BC.
!
! Revision 1.22  2006/08/09 19:19:05  rfiedler
! Allow COUPLED and MOVEDIR to be specified for far field, inflow, and outflow.
!
! Revision 1.21  2006/08/08 17:23:35  rfiedler
! Use MOVEDIR from *.bc to get cnstr_type, not the HDF values.
!
! Revision 1.20  2006/05/02 17:56:42  haselbac
! Cosmetics
!
! Revision 1.19  2006/05/02 17:40:25  fnajjar
! Added STATS key for outflow bc
!
! Revision 1.18  2006/04/17 19:55:44  haselbac
! Bug fix: Added setting of spaceOrder for symmetry and virtual patches
!
! Revision 1.17  2006/04/15 17:01:39  haselbac
! Added ORDER and RECONST params
!
! Revision 1.16  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.15  2006/03/25 21:54:38  haselbac
! Added routines to read input for sype patches
!
! Revision 1.14  2005/11/10 02:27:58  haselbac
! Cosmetics only
!
! Revision 1.13  2005/10/18 02:56:54  haselbac
! Bug fix: Incorrect setting of bcCoupled for non-GENX computations
!
! Revision 1.12  2005/10/14 14:07:50  haselbac
! Significant changes to checks for GENX sims, added GENX support for noslip walls
!
! Revision 1.11  2005/10/05 16:19:23  haselbac
! Fixed problems with names longer than 31 chars
!
! Revision 1.10  2005/10/05 14:04:53  haselbac
! Split no-slip wall sections, cosmetics
!
! Revision 1.9  2005/09/23 18:58:03  haselbac
! Added setting of plotStatsFlag
!
! Revision 1.8  2005/06/09 20:22:13  haselbac
! Removed calls to RFLU_CheckMoveGridInput, changed to MOVEDIR keyword
!
! Revision 1.7  2005/05/04 03:35:18  haselbac
! Removed setting of pPatch%writeGrid, was not used for some time
!
! Revision 1.6  2005/04/27 02:11:50  haselbac
! Added routine to read INFLOW_VELTEMP, made most routines private
!
! Revision 1.5  2005/03/09 15:05:31  haselbac
! Added BC_VIRTUAL, some clean-up
!
! Revision 1.4  2004/12/27 23:29:19  haselbac
! Added parameter for farfield bc
!
! Revision 1.3  2004/10/19 19:28:21  haselbac
! Changed reading of injection boundary input
!
! Revision 1.2  2004/07/28 15:29:20  jferry
! created global variable for spec use
!
! Revision 1.1  2004/07/06 15:14:28  haselbac
! Initial revision
!
! ******************************************************************************



















