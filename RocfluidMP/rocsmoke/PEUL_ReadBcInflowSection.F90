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
! Purpose: read in user input related to inflow boundary condition for
!          Eulerian particles
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
! $Id: PEUL_ReadBcInflowSection.F90,v 1.4 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_ReadBcInflowSection( regions )

  USE ModDataTypes
  USE ModBndPatch,   ONLY : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters

  USE ModInterfaces, ONLY : ReadPatchSection, MakeNumberedKeys
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch, ipt

! ... local variables
  INTEGER, PARAMETER :: NKEYS_MAX  = 20
  INTEGER, PARAMETER :: NPEUL_KEYS = 10

  CHARACTER(CHRLEN) :: RCSIdentString
  CHARACTER(10)     :: keys(NKEYS_MAX)
  CHARACTER(256)    :: fname

  INTEGER :: nKeys, brbeg, brend, prbeg, prend, distrib
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag
  INTEGER :: iKeyDens,iKeyDens0

  LOGICAL :: defined(NKEYS_MAX), alldef

  REAL(RFREAL) :: vals(NKEYS_MAX)

  TYPE(t_patch),  POINTER :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PEUL_ReadBcInflowSection.F90,v $ $Revision: 1.4 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PEUL_ReadBcInflowSection',&
  'PEUL_ReadBcInflowSection.F90' )

! begin -----------------------------------------------------------------------

! define keys

  iKeyDens  = 1
  iKeyDens0 = iKeyDens

  keys(iKeyDens) = 'DENS_'
  CALL MakeNumberedKeys(keys,iKeyDens0+1,'DENS',1,NPEUL_KEYS,1)

  nKeys = iKeyDens0 + NPEUL_KEYS

  IF (nKeys > NKEYS_MAX) CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

! Read smoke BC section from BC input file

  CALL ReadPatchSection( global,IF_INPUT,nKeys,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,fname,defined )

! check if all necessary values defined ---------------------------------------

  DO iReg=brbeg,brend

    IF (regions(iReg)%peulInput%nPtypes > NPEUL_KEYS) &
      CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_INFLOW .AND. &
           patch%bcType<=BC_INFLOW+BC_RANGE)    .AND. &   ! my boundary type,
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        IF (patch%peul%bcSet) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,__LINE__, &
            'PEUL Inflow boundary.' )

!        IF (patch%bcCoupled == BC_EXTERNAL) THEN   ! data from outside
!          patch%peul%distrib = BCDAT_DISTRIB    ! => always distribution
!          CALL ErrorStop( global,ERR_PEUL_EXTERNAL,__LINE__ )
!        ELSE
!          patch%peul%distrib = distrib
!        END IF ! bcCoupled
          patch%peul%distrib = distrib

! ----- check if appropriate values specified
        IF ( patch%peul%distrib == BCDAT_CONSTANT ) THEN
          IF (.NOT. defined(iKeyDens)) THEN
            alldef = .TRUE.
            DO ipt = 1, regions(iReg)%peulInput%nPtypes
              alldef = alldef .AND. defined(iKeyDens0+ipt)
            END DO ! ipt
            IF (.NOT. alldef) CALL ErrorStop(global,ERR_BCVAL_MISSING,__LINE__)
          END IF
        ELSE
          CALL ErrorStop( global,ERR_PEUL_DISTRIB,__LINE__ )
        END IF

! ----- check if extra values specified

        DO ipt = 1, NPEUL_KEYS
          IF (defined(iKeyDens0+ipt)) THEN
            IF (ipt > regions(iReg)%peulInput%nPtypes) THEN
              CALL ErrorStop(global,ERR_PEUL_BCVAL_EXTRA,__LINE__)
            END IF ! ipt
          END IF   ! defined
        END DO     ! ipt

! ----- set flag to BC specified
        patch%peul%bcSet = .TRUE.

      END IF   ! my BC & processor, active
    END DO     ! iPatch
  END DO       ! iReg

! copy values/distribution to variables ---------------------------------------

  DO iReg=brbeg,brend

    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_INFLOW .AND. &
           patch%bcType<=BC_INFLOW+BC_RANGE)    .AND. &   ! my boundary type,
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

          patch%peul%nData = regions(iReg)%peulInput%nPtypes

! ----- allocate memory for the values

        IF (patch%peul%distrib == BCDAT_DISTRIB) THEN
          n1    = ABS(patch%l1end-patch%l1beg)
          n2    = ABS(patch%l2end-patch%l2beg)
          iOff  = n1 + 1
          ijBeg = IndIJ( 0, 0,iOff)
          ijEnd = IndIJ(n1,n2,iOff)
        ELSE
          ijBeg = 0
          ijEnd = 1
        END IF

        NULLIFY(patch%peul%vals)
        IF (patch%peul%nData > 0) THEN
          ALLOCATE( patch%peul%vals(patch%peul%nData,ijBeg:ijEnd), &
                    stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
        END IF

! ----- distribution from file

        IF (patch%peul%distrib==BCDAT_DISTRIB .AND. &
            patch%bcCoupled      /=BC_EXTERNAL  ) THEN

          CALL ErrorStop( global,ERR_PEUL_DISTRIB,__LINE__ )

! ----- distribution from external source / constant value

        ELSE

          IF (ASSOCIATED(patch%peul%vals)) THEN
            IF (defined(iKeyDens)) patch%peul%vals(:,:) = vals(iKeyDens)
            DO ipt = 1, patch%peul%nData
              IF (defined(iKeyDens0+ipt)) patch%peul%vals(ipt,:) = &
                vals(iKeyDens0+ipt)
            END DO ! ipt
          END IF

        END IF  ! distribution?

      END IF    ! bcType, active region on my processor

    END DO      ! iPatch
  END DO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_ReadBcInflowSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_ReadBcInflowSection.F90,v $
! Revision 1.4  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:40:23  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/12/01 21:09:41  haselbac
! Initial revision after changing case
!
! Revision 1.9  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.8  2004/03/02 21:45:12  jferry
! Added check on number of keys
!
! Revision 1.7  2003/11/21 23:20:03  jferry
! Turned off error trap for BC_EXTERNAL
!
! Revision 1.6  2003/09/19 20:34:44  jferry
! Added underscore character to default keys to make key set prefix-free
!
! Revision 1.5  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.4  2003/03/24 23:30:53  jferry
! overhauled rocinteract to allow interaction design to use user input
!
! Revision 1.3  2003/03/04 19:26:47  jferry
! Cleaned up routines that read sections of input files
!
! Revision 1.2  2003/02/12 23:34:48  jferry
! Replaced [io]stat=global%error with local errorFlag
!
! Revision 1.1  2003/02/11 22:52:51  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







