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
! Purpose: read in user input related to APN injection boundary condition.
!
! Description: none.
!
! Input: boundary condition file.
!
! Output: regions = BC data.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ReadBcInjectAPNSection.F90,v 1.5 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReadBcInjectAPNSection( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadPatchSection, RFLO_ReadBcFromFile
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  CHARACTER(10)  :: keys(6)
  CHARACTER(256) :: fname

  INTEGER :: brbeg, brend, prbeg, prend, distrib
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag

  LOGICAL :: defined(6)

  REAL(RFREAL) :: vals(6)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ReadBcInjectAPNSection',&
  'RFLO_ReadBcInjectAPNSection.F90' )

! specify keywords and search for them ----------------------------------------

  keys(1) = 'SDENS'
  keys(2) = 'ACOEFF'
  keys(3) = 'NPOWER'
  keys(4) = 'TEMP'
  keys(5) = 'EXTRAPOL'
  keys(6) = 'MAXCHANGE'

  CALL ReadPatchSection( global,IF_INPUT,6,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,fname,defined )
  
! check if all values defined -------------------------------------------------

  IF (distrib==BCDAT_CONSTANT .AND. &
      (.NOT. (defined(1).eqv..true.) .OR. &
       .NOT. (defined(2).eqv..true.) .OR. &
       .NOT. (defined(3).eqv..true.) .OR. &
       .NOT. (defined(4).eqv..true.))) CALL ErrorStop( global,ERR_BCVAL_MISSING,&
       __LINE__ )

  IF (.NOT. (defined(5).eqv..true.) .OR. &
      .NOT. (defined(6).eqv..true.))   CALL ErrorStop( global,ERR_BCVAL_MISSING,&
      __LINE__ )

! copy values/distribution to variables ---------------------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_INJECTION .AND. &
           patch%bcType<=BC_INJECTION+BC_RANGE) .AND. &   ! my boundary type
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        patch%bcType = BC_INJECTION_APN

        IF (patch%mixt%bcSet.eqv..true.) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,&
          __LINE__,'Injection boundary.' )

        patch%mixt%nData     = 8
        patch%mixt%nSwitches = 1
        patch%mixt%bcSet     = .true.
        IF (patch%bcCoupled == BC_EXTERNAL) THEN   ! data from outside
          patch%mixt%distrib = BCDAT_DISTRIB    ! => always distribution
        ELSE
          patch%mixt%distrib = distrib
        ENDIF

! ----- get value of switch

        ALLOCATE( patch%mixt%switches(patch%mixt%nSwitches), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
        __LINE__ )

          patch%mixt%switches(BCSWI_INJECT_EXTRAP) = EXTRAPOL_CONST
        IF (vals(5) > 0.1) &
          patch%mixt%switches(BCSWI_INJECT_EXTRAP) = EXTRAPOL_LINEAR

        patch%mixt%maxChange = vals(6)

! ----- allocate memory for the values

        n1    = ABS(patch%l1end-patch%l1beg)
        n2    = ABS(patch%l2end-patch%l2beg)
        iOff  = n1 + 1
        ijBeg = IndIJ( 0, 0,iOff)
        ijEnd = IndIJ(n1,n2,iOff)

        ALLOCATE( patch%mixt%vals(patch%mixt%nData,ijBeg:ijEnd), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
        __LINE__ )

! ----- distribution from file

        IF (patch%mixt%distrib==BCDAT_DISTRIB .AND. &
            patch%bcCoupled      /=BC_EXTERNAL  ) THEN
          patch%mixt%nData = 4
          CALL RFLO_ReadBcFromFile( global,fname,patch )
          patch%mixt%nData = 8

! ----- distribution from external source / constant value

        ELSE
          patch%mixt%vals(BCDAT_INJECT_SDENS ,:) = ABS( vals(1) )
          patch%mixt%vals(BCDAT_INJECT_ACOEFF,:) = ABS( vals(2) )
          patch%mixt%vals(BCDAT_INJECT_NPOWER,:) = ABS( vals(3) )
          patch%mixt%vals(BCDAT_INJECT_TEMP  ,:) = ABS( vals(4) )
        ENDIF  ! distribution?

! ----- patch%mixt%distrib for APN always non-constant for later use

        patch%mixt%distrib = BCDAT_DISTRIB

      ENDIF    ! bcType, active region on my processor

    ENDDO      ! iPatch
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ReadBcInjectAPNSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReadBcInjectAPNSection.F90,v $
! Revision 1.5  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2008/10/23 18:20:53  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.2  2006/08/19 15:38:14  mparmar
! Renamed patch variables
!
! Revision 1.1  2006/01/20 06:17:54  wasistho
! initial import
!
!
!******************************************************************************







