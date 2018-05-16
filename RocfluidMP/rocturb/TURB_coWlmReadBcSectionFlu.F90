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
! Purpose: read in user input related to wall layer model on noslip wall bc.
!
! Description: none.
!
! Input: boundary condition file.
!
! Output: regions = BC data pertinent to wlm.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_coWlmReadBcSectionFlu.F90,v 1.4 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_CoWlmReadBcSection( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadPatchSection
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  CHARACTER(15)     :: keys(3)
  CHARACTER(256)    :: fname
  CHARACTER(CHRLEN) :: bcName

  INTEGER :: brbeg, brend, prbeg, prend, distrib, switch
  INTEGER :: ijBeg, ijEnd, errorFlag

  LOGICAL :: defined(3)

  REAL(RFREAL) :: vals(3)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_coWlmReadBcSectionFlu.F90,v $ $Revision: 1.4 $'

  global => regions(1)%global
  CALL RegisterFunction( global,'TURB_CoWlmReadBcSection',&
  'TURB_coWlmReadBcSectionFlu.F90' )

! specify keywords and search for them ----------------------------------------

  keys(1) = 'MODEL'
  keys(2) = 'REFPOINT'
  keys(3) = 'ROUGHNESS'

  CALL ReadPatchSection( global,IF_INPUT,3,keys,vals, &
                         prbeg,prend,distrib,fname,bcName,defined )

  IF ( prend > global%nPatches ) THEN 
    CALL ErrorStop(global,ERR_PATCH_RANGE,__LINE__)
  END IF

! get switches & check if all necessary values defined ------------------------

  brbeg = LBOUND(regions,1)   ! temporary for now before zonal modeling apply
  brend = UBOUND(regions,1)

  DO iReg=brbeg,brend
    IF (regions(iReg)%mixtInput%turbModel<=TURB_MODEL_NONE) THEN
      CALL ErrorStop( global,ERR_TURB_REGION,__LINE__,'Wall model unapplicable.' )
    ENDIF

    regions(iReg)%turbInput%wallModel = 0
    regions(iReg)%turbInput%wallRough = 0._RFREAL
    DO iPatch=prbeg,MIN(prend,regions(iReg)%grid%nPatches)

      patch => regions(iReg)%patches(iPatch)

      IF (patch%bcType>=BC_NOSLIPWALL .AND. &
          patch%bcType<=BC_NOSLIPWALL+BC_RANGE) THEN   ! my boundary type,

        IF (patch%bcCoupled == BC_EXTERNAL) THEN  ! data from outside
          patch%valBola%distrib = BCDAT_DISTRIB    ! => always distribution
        ELSE
          patch%valBola%distrib = distrib
        ENDIF

! ----- get input switches

        IF (defined(1)) THEN
          patch%valBola%switches(WLM_INPUT_MODEL)=MAX( 0,INT(vals(1)+0.5_RFREAL) )
        ELSE
          CALL ErrorStop( global,ERR_NO_BCSWITCH,__LINE__,'WLM model missing.' )
        ENDIF

        IF (patch%valBola%switches(WLM_INPUT_MODEL) == WLM_MODEL_NOMODEL) goto 999

        IF (defined(2)) THEN
          patch%valBola%switches(WLM_INPUT_REFPOINT) = INT(vals(2)+0.5_RFREAL)
        ELSE
          CALL ErrorStop( global,ERR_NO_BCSWITCH,__LINE__, &
                          'WLM reference point missing.' )
        ENDIF

! ----- allocate memory for the roughness distribution

        patch%valBola%nData = patch%valBola%nData + 2*TENSOR_ALL_NELM + 3 + 9
        ijBeg = 1
        ijEnd = patch%nBFaces

        ALLOCATE( patch%valBola%vals(ijBeg:ijEnd,patch%valBola%nData), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

        IF (patch%valBola%distrib==BCDAT_DISTRIB) THEN
! ------- roughness distribution from file
!         CALL TURB_WlmReadRoughness( global,fname,patch )
          CALL ErrorStop( global,ERR_VAL_BCVAL,__LINE__,'No variable roughness yet'  )
        ELSE
! ------- distribution has constant value
          IF (defined(3)) THEN
            patch%valBola%vals(:,WLM_VALS_ROUGH) = vals(3)
          ELSE
            patch%valBola%vals(:,WLM_VALS_ROUGH) = 0._RFREAL
          ENDIF
        ENDIF  ! distribution?

        IF (patch%valBola%switches(WLM_INPUT_MODEL) == WLM_MODEL_EXTERN) THEN
! ------- wall stress distribution from file
!         CALL TURB_WlmReadWallStress( global,fname,patch )
        ENDIF

999 CONTINUE

      ENDIF   ! my BC

! --- copy max wlm param values to input param for screen-print
!      IF (patch%bcType>=BC_NOSLIPWALL .AND. & 
!          patch%bcType<=BC_NOSLIPWALL+BC_RANGE) THEN

!        regions(iReg)%turbInput%wallModel= &
!          MAX( regions(iReg)%turbInput%wallModel,INT(vals(1)+0.5_RFREAL) )

!        regions(iReg)%turbInput%wlmRefPoint= &
!          MAX( regions(iReg)%turbInput%wlmRefPoint,INT(vals(2)+0.5_RFREAL) )

!        regions(iReg)%turbInput%wallRough= &
!          MAX( regions(iReg)%turbInput%wallRough,vals(3) )
                                          
!      ENDIF  ! my boundary type 
    ENDDO  ! iPatch
  ENDDO  ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_CoWlmReadBcSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_coWlmReadBcSectionFlu.F90,v $
! Revision 1.4  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.1  2004/03/25 04:42:58  wasistho
! prepared for RFLU
!
!
!
!******************************************************************************







