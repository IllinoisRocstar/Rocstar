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
! Purpose: Initialize input parameters for turbulence to default values.
!
! Description: none.
!
! Input: none.
!
! Output: regions = initial input values.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_InitInputValues.F90,v 1.15 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_InitInputValues( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModTurbulence, ONLY : t_turb_input
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  USE TURB_ModParameters
#ifdef RFLU
  USE ModInterfaces, ONLY           : RFLU_CreateGrid
  USE RFLU_ModDimensions, ONLY      : RFLU_ReadDimensions
  USE RFLU_ModReadBcInputFile, ONLY : RFLU_ReadBcInputFileWrapper
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iPatch
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN)           :: RCSIdentString
  TYPE(t_patch), POINTER      :: patch1
  TYPE(t_global), POINTER     :: global
  TYPE(t_turb_input), POINTER :: input
#ifdef RFLU
  TYPE(t_region), POINTER     :: pRegion
#endif
  INTEGER :: errorFlag

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_InitInputValues.F90,v $'

  global => regions(1)%global
  CALL RegisterFunction( global,'TURB_InitInputValues',&
  'TURB_InitInputValues.F90' )

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Entering TURB_InitInputValues...'
  END IF ! global%verbLevel

! global values ---------------------------------------------------------------

  global%turbActive        = .FALSE.
  global%turbCalcWDist     = .FALSE.
  global%turbWorkUnused    = .TRUE.
  global%turbWallDim       = 0
  global%turbCalcWDistFreq = CALCWDIST_INI

! region related values -------------------------------------------------------

#ifdef RFLO
  DO iReg=1,global%nRegions
#endif
#ifdef RFLU
!$ Start Temporary Procedure
  IF (global%moduleType == MODULE_TYPE_SOLVER) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      pRegion => regions(iReg)
      CALL RFLU_ReadDimensions( pRegion )         ! needed to get grid%nPatches
      CALL RFLU_CreateGrid( pRegion )             ! needed to allocate patches
      CALL RFLU_ReadBcInputFileWrapper( pRegion ) ! needed to get patch%bcType
    ENDDO
  ELSE
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%grid%nPatches = 0
    ENDDO
  ENDIF   ! moduleType
!$ End Temporary

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
#endif

    input => regions(iReg)%turbInput

! - general input parameters
 
    input%nOutField      = 1
    input%modelClass     = MODEL_NONE
    input%nZof           = 0

! - RaNS/DES input parameters

    input%wDistMethod    = WDIST_DIRECT
    input%cDes           =  0.65_RFREAL
    input%smoocf         = -1._RFREAL
    input%spaceDiscr     = RANS_DISCR_UPW
    input%vis2           =  0.50_RFREAL
    input%vis4           =  0._RFREAL
    input%spaceOrder     = RANS_DISCR_ORD1
    input%functV1        = SA_FV1_POW3

! - LES input parameters

    input%cSmag          = 0.1_RFREAL
    input%xyzSmag(:)     = -HUGE( 1.0_RFREAL )
    input%filterType     = FILTYPE_UNIFORM
    input%deltaType      = DELTYPE_CBRT
    input%filterWidth(:) = FILWIDTH_ONE
    input%homDir(:)      = OFF
    input%engModel       = ACTIVE
    input%calcVort       = CALCVORT_NO

! - wlm input parameters

    input%wallModel      = WLM_MODEL_NOMODEL ! These r initial single vals per  
    input%wlmRefPoint    = 1                 ! region, and copied to all patches 
                                             ! below. Later they are assigned
                                             ! (TURB_CoWlmReadBcSectionFlo/Flu)
                                             ! max. values over all patches and
                                             ! presented in TURB_PrinUserInput
#ifdef RFLO
    DO iPatch=1,regions(iReg)%nPatches

! --- only finest level defined yet

      patch1 => regions(iReg)%levels(1)%patches(iPatch)

      patch1%turb%nData = 0                              ! used by TBC
      patch1%turb%nSwitches = 0
      patch1%turb%distrib = BCDAT_CONSTANT 

      patch1%valBola%bcSet = .false.   

      IF ((patch1%bcType>=BC_NOSLIPWALL .AND. &
           patch1%bcType<=BC_NOSLIPWALL+BC_RANGE) .AND. &   ! my boundary type,
           regions(iReg)%procid==global%myProcid  .AND. &   ! region active and
           regions(iReg)%active==ACTIVE) THEN               ! on my processor
#endif
#ifdef RFLU
    DO iPatch=1,regions(iReg)%grid%nPatches     ! to getrid of above RFLU calls
                                                ! this should move to TURB_Read
! --- defined in current level                  ! BcInputFile to be called from 
                                                ! RFLU_ReadBcInputFileWrapper.
      patch1 => regions(iReg)%patches(iPatch)   ! CheckParamInput for FLU patch
                                                ! moves after TURB_ReadBcInputF

      patch1%turb%nData = 0                         ! used by TBC
      patch1%turb%nSwitches = 0                     ! RFLU turb TBC is
      patch1%turb%distrib = BCDAT_CONSTANT          ! cycled for now (RFLU_
                                                       ! AllocateMemoryTbc)
      IF (patch1%bcType>=BC_NOSLIPWALL .AND. &
          patch1%bcType<=BC_NOSLIPWALL+BC_RANGE) THEN  ! my boundary type
#endif

! ----- initiate parameters in all no-slip patches

        patch1%valBola%nData = 0
        patch1%valBola%nSwitches = 0
        patch1%valBola%distrib = BCDAT_CONSTANT 
        patch1%valBola%nSwitches = patch1%valBola%nSwitches + WLM_NSWITCH

        ALLOCATE( patch1%valBola%switches(patch1%valBola%nSwitches), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! ----- initiate wlm switches in all no-slip patches

        patch1%valBola%switches(WLM_INPUT_MODEL)    = input%wallModel
        patch1%valBola%switches(WLM_INPUT_REFPOINT) = input%wlmRefPoint

! ----- allocation and valuation of patch1%valBola%vals, which roughness data
!       is part of, starts in TURB_CoWlmReadBcSectionFlo/Flu

      ENDIF ! bcType
    ENDDO   ! iPatch
  ENDDO  ! iReg

! finalize --------------------------------------------------------------------

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Leaving TURB_InitInputValues.'
  END IF ! global%verbLevel

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_InitInputValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_InitInputValues.F90,v $
! Revision 1.15  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.14  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.13  2006/08/19 15:40:35  mparmar
! Renamed patch variables
!
! Revision 1.12  2006/03/03 23:03:19  wasistho
! added moduleType condition for Rocflu
!
! Revision 1.11  2006/02/11 02:46:06  wasistho
! early obtaining nPatches, patches and bcType from rocflu
!
! Revision 1.10  2006/02/05 04:43:05  wasistho
! skip loop over patches in Rocflu
!
! Revision 1.9  2006/02/04 04:59:04  wasistho
! added enter and leave statements
!
! Revision 1.8  2006/01/12 09:48:01  wasistho
! enabled tripping fixed Smagorinsky
!
! Revision 1.7  2006/01/06 07:02:44  wasistho
! set valturb%nData to 0
!
! Revision 1.6  2005/01/12 01:13:23  wasistho
! removed single quote signs since SUN has trouble with it
!
! Revision 1.5  2004/04/20 20:45:57  wasistho
! added user option for frequency of computing wall distance
!
! Revision 1.4  2004/03/20 03:28:29  wasistho
! prepared for RFLU
!
! Revision 1.3  2004/03/19 02:47:08  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/08 23:30:56  wasistho
! changed turb nomenclature
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.12  2004/02/19 04:02:58  wasistho
! added new rans/SA parameter VISCFUNCTION
!
! Revision 1.11  2004/02/14 03:42:48  wasistho
! added new WLM parameter: reference point
!
! Revision 1.10  2004/02/11 03:24:20  wasistho
! added feature: variable number of turbulence output fields
!
! Revision 1.9  2003/10/26 00:10:44  wasistho
! added multiple discr.types and order
!
! Revision 1.8  2003/10/15 03:40:49  wasistho
! added 2nd order dissipation coeff. k2
!
! Revision 1.7  2003/10/09 20:48:53  wasistho
! added DES lengthscale coefficient CDES
!
! Revision 1.6  2003/10/07 02:05:49  wasistho
! initial installation of RaNS-SA and DES
!
! Revision 1.5  2003/08/02 00:19:13  wasistho
! set initial calcVort to zero
!
! Revision 1.4  2003/08/01 22:17:37  wasistho
! prepared rocturb for Genx
!
! Revision 1.3  2003/07/22 02:58:58  wasistho
! prepare more accurate rocturb restart
!
! Revision 1.2  2003/05/31 01:46:43  wasistho
! installed turb. wall layer model
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************







