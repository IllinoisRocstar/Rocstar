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
! Purpose: read in user input related to slip-wall boundary condition.
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
! $Id: RFLO_ReadBcSlipWallSection.F90,v 1.9 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReadBcSlipWallSection( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadPatchSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  CHARACTER(10)  :: keys(6)
  CHARACTER(256) :: fname

  INTEGER :: brbeg, brend, prbeg, prend, distrib, errorFlag

  LOGICAL :: defined(6)

  REAL(RFREAL) :: vals(6)

  TYPE(t_patch), POINTER :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ReadBcSlipWallSection',&
  'RFLO_ReadBcSlipWallSection.F90' )

! specify keywords and search for them ----------------------------------------

  keys(1) = 'EXTRAPOL'
  keys(2) = 'MAXCHANGE'
  keys(3) = 'MOTION'
  keys(4) = 'XVEL'
  keys(5) = 'YVEL'
  keys(6) = 'ZVEL'

  vals(3:6) = 0._RFREAL   ! initialize vars not read by all regs and patches

  CALL ReadPatchSection( global,IF_INPUT,6,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,fname,defined )
  
! check if all values defined -------------------------------------------------

  IF (.NOT. (defined(1).eqv..true.) .OR. &
      .NOT. (defined(2).eqv..true.)) CALL ErrorStop( global,ERR_BCVAL_MISSING,&
      __LINE__ )

! copy values/distribution to variables ---------------------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_SLIPWALL .AND. &
           patch%bcType<=BC_SLIPWALL+BC_RANGE) .AND. &    ! my boundary type
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        IF (patch%mixt%bcSet.eqv..true.) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,__LINE__,'Slip-wall boundary.' )

        patch%mixt%nData     = 0
        patch%mixt%nSwitches = 1
        patch%mixt%distrib   = BCDAT_CONSTANT
        patch%mixt%bcSet     = .true.

! ----- get value of switch

        ALLOCATE( patch%mixt%switches(patch%mixt%nSwitches), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

        patch%mixt%switches(BCSWI_SLIPW_EXTRAP)   = EXTRAPOL_CONST
        IF (vals(1) > 0.1) &
          patch%mixt%switches(BCSWI_SLIPW_EXTRAP) = EXTRAPOL_LINEAR

        patch%mixt%maxChange = vals(2)

        IF (vals(3) > 0.1) THEN
          patch%mixt%setMotion      = .true.
          patch%mixt%bndVel(XCOORD) = vals(4)
          patch%mixt%bndVel(YCOORD) = vals(5)
          patch%mixt%bndVel(ZCOORD) = vals(6)
        ENDIF

      ENDIF    ! bcType, active region on my processor

      IF ((defined(3).eqv..true.) .AND. (vals(3) > 0.1)) global%internDeform = 1

    ENDDO      ! iPatch
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ReadBcSlipWallSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReadBcSlipWallSection.F90,v $
! Revision 1.9  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:21  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2008/10/23 18:20:53  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.6  2006/08/19 15:38:21  mparmar
! Renamed patch variables
!
! Revision 1.5  2005/09/30 03:31:00  wasistho
! initialize vals(3:6)
!
! Revision 1.4  2005/09/30 01:31:03  wasistho
! move internDeform assignment to all procs
!
! Revision 1.3  2005/09/30 01:03:27  wasistho
! assigned global%internDeform
!
! Revision 1.2  2005/09/09 03:22:34  wasistho
! added boundary motion parameters
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.6  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.3  2003/05/19 21:18:21  jblazek
! Automated switch to 0th-order extrapolation at slip walls and injection boundaries.
!
! Revision 1.2  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/19 00:40:30  jblazek
! Added utility (rflosurf) to write out surface grids for GenX.
!
! Revision 1.7  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.6  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.5  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.4  2002/09/17 13:43:00  jferry
! Added Time-dependent boundary conditions
!
! Revision 1.3  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/03/30 00:50:49  jblazek
! Cleaned up with flint.
!
! Revision 1.1  2002/02/27 18:38:20  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
!******************************************************************************







