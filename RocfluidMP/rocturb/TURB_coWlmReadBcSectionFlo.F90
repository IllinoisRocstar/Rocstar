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
! $Id: TURB_coWlmReadBcSectionFlo.F90,v 1.7 2009/08/26 12:28:53 mtcampbe Exp $
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

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  CHARACTER(15)     :: keys(3)
  CHARACTER(256)    :: fname

  INTEGER :: brbeg, brend, prbeg, prend, distrib, switch
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag

  LOGICAL :: defined(3)

  REAL(RFREAL) :: vals(3)

  TYPE(t_patch), POINTER  :: patch1
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_coWlmReadBcSectionFlo.F90,v $ $Revision: 1.7 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'TURB_CoWlmReadBcSection',&
  'TURB_coWlmReadBcSectionFlo.F90' )

! specify keywords and search for them ----------------------------------------

  keys(1) = 'MODEL'
  keys(2) = 'REFPOINT'
  keys(3) = 'ROUGHNESS'

  CALL ReadPatchSection( global,IF_INPUT,3,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,fname,defined )

! get switches & check if all necessary values defined ------------------------

  DO iReg=brbeg,brend
    IF (regions(iReg)%mixtInput%turbModel<=TURB_MODEL_NONE) THEN
      CALL ErrorStop( global,ERR_TURB_REGION,__LINE__,'Wall model unapplicable.' )
    ENDIF

    regions(iReg)%turbInput%wallModel = 0
    regions(iReg)%turbInput%wallRough = 0._RFREAL
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch1 => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch1%bcType>=BC_NOSLIPWALL .AND. &
           patch1%bcType<=BC_NOSLIPWALL+BC_RANGE) .AND. &   ! my boundary type,
           regions(iReg)%procid==global%myProcid  .AND. &   ! region active and
           regions(iReg)%active==ACTIVE) THEN               ! on my processor

        IF (patch1%valBola%bcSet .eqv. .true.) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,__LINE__,'Wall layer model.' )

        IF (patch1%bcCoupled == BC_EXTERNAL) THEN   ! data from outside
          patch1%valBola%distrib = BCDAT_DISTRIB    ! => always distribution
        ELSE
          patch1%valBola%distrib = distrib
        ENDIF

! ----- get input switches

        IF (defined(1) .eqv. .true.) THEN
          patch1%valBola%switches(WLM_INPUT_MODEL)=MAX( 0,INT(vals(1)+0.5_RFREAL) )
        ELSE
          CALL ErrorStop( global,ERR_NO_BCSWITCH,__LINE__,'WLM model missing.' )
        ENDIF

        IF (patch1%valBola%switches(WLM_INPUT_MODEL) == WLM_MODEL_NOMODEL) goto 999

        IF (defined(2) .eqv. .true.) THEN
          patch1%valBola%switches(WLM_INPUT_REFPOINT) = INT(vals(2)+0.5_RFREAL)
        ELSE
          CALL ErrorStop( global,ERR_NO_BCSWITCH,__LINE__, &
                          'WLM reference point missing.' )
        ENDIF

!        IF (defined(4) .eqv. .true.) THEN
!          patch1%valBola%switches(WLM_INPUT_HOMDIR) = INT(vals(4)+0.5_RFREAL)
!        ELSE
!          CALL ErrorStop( global,ERR_NO_BCSWITCH,__LINE__,'WLM homdir missing.' )
!        ENDIF

! ----- allocate memory for the roughness distribution

        patch1%valBola%nData = patch1%valBola%nData + 2*TENSOR_ALL_NELM + 3 + 9
        n1    = ABS(patch1%l1end-patch1%l1beg)
        n2    = ABS(patch1%l2end-patch1%l2beg)
        iOff  = n1 + 1
        ijBeg = IndIJ( 0, 0,iOff)
        ijEnd = IndIJ(n1,n2,iOff)

        ALLOCATE( patch1%valBola%vals(ijBeg:ijEnd,patch1%valBola%nData), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

        IF (patch1%valBola%distrib==BCDAT_DISTRIB) THEN
! ------- roughness distribution from file
!         CALL TURB_WlmReadRoughness( global,fname,patch )
          CALL ErrorStop( global,ERR_VAL_BCVAL,__LINE__,'No variable roughness yet'  )
        ELSE
! ------- distribution has constant value
          IF (defined(3) .eqv. .true.) THEN
            patch1%valBola%vals(:,WLM_VALS_ROUGH) = vals(3)
          ELSE
            patch1%valBola%vals(:,WLM_VALS_ROUGH) = 0._RFREAL
          ENDIF
        ENDIF  ! distribution?

        IF (patch1%valBola%switches(WLM_INPUT_MODEL) == WLM_MODEL_EXTERN) THEN
! ------- wall stress distribution from file
!         CALL TURB_WlmReadWallStress( global,fname,patch )
        ENDIF

999 CONTINUE

! ----- set flag to BC specified
        patch1%valBola%bcSet = .true.

      ENDIF   ! my BC & processor, active

! --- copy max wlm param values to input param for screen-print
!      IF (patch1%bcType>=BC_NOSLIPWALL .AND. & 
!          patch1%bcType<=BC_NOSLIPWALL+BC_RANGE) THEN

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
! $Log: TURB_coWlmReadBcSectionFlo.F90,v $
! Revision 1.7  2009/08/26 12:28:53  mtcampbe
! Ported to Hera.   Fixed logical expression syntax errors.  Replaced all
! IF (logical_variable)  with IF (logical_variable .eqv. .true.) as
! consistent with the specification.  Also changed: IF( ASSOCIATED(expr) )
! to IF ( ASSOCIATED(expr) .eqv. .true. ).   Intel compilers produce code
! which silently fails for some mal-formed expressions, so these changes
! are a net which should ensure that they are evaluated as intended.
!
! Revision 1.6  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2004/03/25 04:40:41  wasistho
! prepared for RFLU
!
! Revision 1.3  2004/03/23 03:35:00  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/13 03:13:56  wasistho
! get rid of flo/flu identifier in TURB_Co.. routines
!
! Revision 1.1  2004/03/08 23:33:31  wasistho
! changed turb nomenclature
!
! Revision 1.4  2004/03/02 03:50:43  wasistho
! bug fixed nvals=2 to nvals=3, forgot colon after Id and Log
!
!
!******************************************************************************







