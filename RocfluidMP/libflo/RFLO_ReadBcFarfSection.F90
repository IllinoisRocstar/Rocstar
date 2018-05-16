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
! Purpose: read in user input related to far field boundary condition.
!
! Description: none.
!
! Input: boundary condition file.
!
! Output: regions = BC data
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ReadBcFarfSection.F90,v 1.5 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReadBcFarfSection( regions )

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
  CHARACTER(10)  :: keys(5)
  CHARACTER(256) :: fname

  INTEGER :: brbeg, brend, prbeg, prend, distrib
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag

  LOGICAL :: defined(5)

  REAL(RFREAL) :: vals(5)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ReadBcFarfSection',&
  'RFLO_ReadBcFarfSection.F90' )

! specify keywords and search for them ----------------------------------------

  keys(1) = 'MACH'
  keys(2) = 'ATTACK'
  keys(3) = 'SLIP'
  keys(4) = 'PRESS'
  keys(5) = 'TEMP'

  CALL ReadPatchSection( global,IF_INPUT,5,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,fname,defined )

! check if all values defined -------------------------------------------------

  IF (distrib==BCDAT_CONSTANT .AND. &
      (.NOT. (defined(1).eqv..true.) .OR. &
       .NOT. (defined(2).eqv..true.) .OR. &
       .NOT. (defined(3).eqv..true.) .OR. &
       .NOT. (defined(4).eqv..true.) .OR. &
       .NOT. (defined(5).eqv..true.))) CALL ErrorStop( global,ERR_BCVAL_MISSING,&
       __LINE__ )

! copy values/distribution to variables ---------------------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_FARFIELD .AND. &
           patch%bcType<=BC_FARFIELD+BC_RANGE)  .AND. &   ! my boundary type
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        IF (patch%mixt%bcSet.eqv..true.) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,&
          __LINE__,'Farfield boundary.' )

        patch%mixt%nData     = 5
        patch%mixt%nSwitches = 0
        patch%mixt%bcSet     = .true.
        IF (patch%bcCoupled == BC_EXTERNAL) THEN   ! data from outside
          patch%mixt%distrib = BCDAT_DISTRIB    ! => always distribution
        ELSE
          patch%mixt%distrib = distrib
        ENDIF

! ----- allocate memory for the values

        IF (patch%mixt%distrib == BCDAT_DISTRIB) THEN
          n1    = ABS(patch%l1end-patch%l1beg)
          n2    = ABS(patch%l2end-patch%l2beg)
          iOff  = n1 + 1
          ijBeg = IndIJ( 0, 0,iOff)
          ijEnd = IndIJ(n1,n2,iOff)
        ELSE
          ijBeg = 0
          ijEnd = 1
        ENDIF
        ALLOCATE( patch%mixt%vals(patch%mixt%nData,ijBeg:ijEnd), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
        __LINE__ )

! ----- distribution from file

        IF (patch%mixt%distrib==BCDAT_DISTRIB .AND. &
            patch%bcCoupled      /=BC_EXTERNAL  ) THEN
          CALL RFLO_ReadBcFromFile( global,fname,patch )

          patch%mixt%vals(BCDAT_FARF_ATTACK,:) = &
            patch%mixt%vals(BCDAT_FARF_ATTACK,:)*global%rad
          patch%mixt%vals(BCDAT_FARF_SLIP  ,:) = &
            patch%mixt%vals(BCDAT_FARF_SLIP  ,:)*global%rad

! ----- distribution from external source / constant value

        ELSE
          patch%mixt%vals(BCDAT_FARF_MACH  ,:) = vals(1)
          patch%mixt%vals(BCDAT_FARF_ATTACK,:) = vals(2)*global%rad
          patch%mixt%vals(BCDAT_FARF_SLIP  ,:) = vals(3)*global%rad
          patch%mixt%vals(BCDAT_FARF_PRESS ,:) = vals(4)
          patch%mixt%vals(BCDAT_FARF_TEMP  ,:) = vals(5)
        ENDIF  ! distribution?

      ENDIF    ! bcType, active region on my processor

    ENDDO      ! iPatch
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ReadBcFarfSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReadBcFarfSection.F90,v $
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
! Revision 1.2  2006/08/19 15:38:07  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.6  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.3  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.2  2003/02/11 22:30:21  jferry
! Re-worked BC and TBC input routines to add multi-physics capability
!
! Revision 1.1  2002/10/19 00:40:30  jblazek
! Added utility (rflosurf) to write out surface grids for GenX.
!
! Revision 1.12  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.11  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.10  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.9  2002/09/17 13:43:00  jferry
! Added Time-dependent boundary conditions
!
! Revision 1.8  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.7  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.6  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.5  2002/02/04 15:30:25  jblazek
! Added injection boundary condition.
!
! Revision 1.4  2002/01/31 20:56:30  jblazek
! Added basic boundary conditions.
!
! Revision 1.3  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.2  2001/12/22 00:09:39  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1  2001/12/08 00:18:42  jblazek
! Added routines to read BC input file.
!
!******************************************************************************







