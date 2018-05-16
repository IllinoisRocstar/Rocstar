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
! Purpose: read in user input related to MRate injection boundary condition.
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
! $Id: RFLO_ReadBcInjectMrateSection.F90,v 1.6 2010/02/18 21:47:38 juzhang Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReadBcInjectMrateSection( regions )

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
  CHARACTER(10)  :: keys(4)
  CHARACTER(256) :: fname

  INTEGER :: brbeg, brend, prbeg, prend, distrib
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag

  LOGICAL :: defined(4)

  REAL(RFREAL) :: vals(4)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ReadBcInjectMrateSection',&
  'RFLO_ReadBcInjectMrateSection.F90' )

! specify keywords and search for them ----------------------------------------

  keys(1) = 'MFRATE'
  keys(2) = 'TEMP'
  keys(3) = 'EXTRAPOL'
  keys(4) = 'MAXCHANGE'

  CALL ReadPatchSection( global,IF_INPUT,4,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,fname,defined )
  
! check if all values defined -------------------------------------------------

  IF (distrib==BCDAT_CONSTANT .AND. &
      (.NOT. (defined(1).eqv..true.) .OR. &
       .NOT. (defined(2).eqv..true.))) CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )

  IF (.NOT. (defined(3).eqv..true.) .OR. &
      .NOT. (defined(4).eqv..true.))   CALL ErrorStop( global,ERR_BCVAL_MISSING,__LINE__ )

! copy values/distribution to variables ---------------------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_INJECTION .AND. &
           patch%bcType<=BC_INJECTION+BC_RANGE) .AND. &   ! my boundary type
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

!!JZ comment out temporarily..
!!        patch%bcType = BC_INJECTION_MRATE

        IF (patch%mixt%bcSet.eqv..true.) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,__LINE__,'Injection boundary.' )

        patch%mixt%nData     = 5
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
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

          patch%mixt%switches(BCSWI_INJECT_EXTRAP) = EXTRAPOL_CONST
        IF (vals(3) > 0.1) &
          patch%mixt%switches(BCSWI_INJECT_EXTRAP) = EXTRAPOL_LINEAR

        patch%mixt%maxChange = vals(4)

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
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! ----- distribution from file

        IF (patch%mixt%distrib==BCDAT_DISTRIB .AND. &
            patch%bcCoupled      /=BC_EXTERNAL  ) THEN
          patch%mixt%nData = 2
          CALL RFLO_ReadBcFromFile( global,fname,patch )
          patch%mixt%nData = 5

! ----- distribution from external source / constant value

        ELSE
          patch%mixt%vals(BCDAT_INJECT_MFRATE,:) = vals(1)
          patch%mixt%vals(BCDAT_INJECT_TEMP  ,:) = vals(2)
        ENDIF  ! distribution?

      ENDIF    ! bcType, active region on my processor

    ENDDO      ! iPatch
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ReadBcInjectMrateSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReadBcInjectMrateSection.F90,v $
! Revision 1.6  2010/02/18 21:47:38  juzhang
! Heat transfer bc for non-propellant surface documented in Rocburn_PY_HT.pdf in Rocburn_PY directory is implemented within Rocburn_PY. Major changes were made to Rocburn, Rocman3, RocfluidMP/genx, RocfluidMP/modflo directories.
!
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
! Revision 1.2  2006/08/19 15:38:15  mparmar
! Renamed patch variables
!
! Revision 1.1  2006/01/20 06:17:54  wasistho
! initial import
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.7  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.4  2003/05/19 21:18:21  jblazek
! Automated switch to 0th-order extrapolation at slip walls and injection boundaries.
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
! Revision 1.13  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.12  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.11  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.10  2002/09/17 13:43:00  jferry
! Added Time-dependent boundary conditions
!
! Revision 1.9  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.8  2002/02/27 18:38:20  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
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







