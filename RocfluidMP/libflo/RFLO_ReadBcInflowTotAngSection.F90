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
! Purpose: read in user input related to inflow boundary condition.
!
! Description: present inflow bc is based on total pressure, total temperature
!              and flow angle.
!
! Input: boundary condition file.
!
! Output: regions = BC data.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ReadBcInflowTotAngSection.F90,v 1.5 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReadBcInflowTotAngSection( regions )

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
  CHARACTER(10)  :: keys(7)
  CHARACTER(256) :: fname

  INTEGER :: brbeg, brend, prbeg, prend, distrib, switch
  INTEGER :: n1, n2, iOff, ijBeg, ijEnd, errorFlag

  LOGICAL :: defined(7)

  REAL(RFREAL) :: vals(7)

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ReadBcInflowTotAngSection',&
  'RFLO_ReadBcInflowTotAngSection.F90' )

! specify keywords and search for them ----------------------------------------

  keys(1) = 'TYPE'
  keys(2) = 'FIXED'
  keys(3) = 'PTOT'
  keys(4) = 'TTOT'
  keys(5) = 'BETAH'
  keys(6) = 'BETAV'
  keys(7) = 'MACH'

  CALL ReadPatchSection( global,IF_INPUT,7,keys,vals,brbeg,brend, &
                         prbeg,prend,distrib,fname,defined )

! get switches & check if all necessary values defined ------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_INFLOW .AND. &
           patch%bcType<=BC_INFLOW+BC_RANGE)    .AND. &   ! my boundary type,
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        patch%bcType = BC_INFLOW_TOTANG

        IF (patch%mixt%bcSet.eqv..true.) &
          CALL ErrorStop( global,ERR_PATCH_OVERSPEC,&
          __LINE__,'Inflow boundary.' )

        patch%mixt%nSwitches = 2
        IF (patch%bcCoupled == BC_EXTERNAL) THEN   ! data from outside
          patch%mixt%distrib = BCDAT_DISTRIB    ! => always distribution
        ELSE
          patch%mixt%distrib = distrib
        ENDIF

        ALLOCATE( patch%mixt%switches(patch%mixt%nSwitches), &
                  stat=errorFlag )
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
        __LINE__ )

! ----- check if switch defined
        IF (defined(1).eqv..true.) THEN
          patch%mixt%switches(BCSWI_INFLOW_TYPE)   = BCOPT_SUBSONIC
          IF (vals(1) < 0.1) &
            patch%mixt%switches(BCSWI_INFLOW_TYPE) = BCOPT_SUPERSONIC
          IF (vals(1) > 1.9) &
            patch%mixt%switches(BCSWI_INFLOW_TYPE) = BCOPT_MIXED
        ELSE
          CALL ErrorStop( global,ERR_NO_BCSWITCH,&
          __LINE__,'(inflow type).' )
        ENDIF

        IF (defined(2).eqv..true.) THEN
          IF (vals(2) < 0.1) &
            patch%mixt%switches(BCSWI_INFLOW_FIXED) = BCOPT_FIXED_NO
          IF (vals(2) > 0.9) &
            patch%mixt%switches(BCSWI_INFLOW_FIXED) = BCOPT_FIXED_YES
        ELSE
          patch%mixt%switches(BCSWI_INFLOW_FIXED) = BCOPT_FIXED_NO
        ENDIF

! ----- check if appropriate values specified
        IF (patch%mixt%switches(BCSWI_INFLOW_TYPE) == BCOPT_SUBSONIC) THEN
          IF (patch%mixt%distrib==BCDAT_CONSTANT .AND. &
              (.NOT. (defined(3).eqv..true.) .OR. &
               .NOT. (defined(4).eqv..true.) .OR. &
               .NOT. (defined(5).eqv..true.) .OR. &
               .NOT. (defined(6).eqv..true.))) CALL ErrorStop( global,ERR_BCVAL_MISSING,&
               __LINE__ )
        ENDIF
        IF (patch%mixt%switches(BCSWI_INFLOW_TYPE) == BCOPT_SUPERSONIC .OR. &
            patch%mixt%switches(BCSWI_INFLOW_TYPE) == BCOPT_MIXED) THEN
          IF (patch%mixt%distrib==BCDAT_CONSTANT .AND. &
              (.NOT. (defined(3).eqv..true.) .OR. &
               .NOT. (defined(4).eqv..true.) .OR. &
               .NOT. (defined(5).eqv..true.) .OR. &
               .NOT. (defined(6).eqv..true.) .OR. &
               .NOT. (defined(7).eqv..true.))) CALL ErrorStop( global,ERR_BCVAL_MISSING,&
               __LINE__ )
        ENDIF

! ----- set flag to BC specified
        patch%mixt%bcSet = .true.

      ENDIF   ! my BC & processor, active
    ENDDO
  ENDDO

! copy values/distribution to variables ---------------------------------------

  DO iReg=brbeg,brend
    DO iPatch=prbeg,MIN(prend,regions(iReg)%nPatches)

      patch => regions(iReg)%levels(1)%patches(iPatch)

      IF ((patch%bcType>=BC_INFLOW .AND. &
           patch%bcType<=BC_INFLOW+BC_RANGE)    .AND. &   ! my boundary type,
          regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor

        switch = patch%mixt%switches(BCSWI_INFLOW_TYPE)
        IF (switch == BCOPT_SUBSONIC) THEN
          patch%mixt%nData = 4
        ELSE
          patch%mixt%nData = 5
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

          IF (switch == BCOPT_SUBSONIC) THEN
            patch%mixt%vals(BCDAT_INFLOW_BETAH,:) = &
              patch%mixt%vals(BCDAT_INFLOW_BETAH,:)*global%rad
            patch%mixt%vals(BCDAT_INFLOW_BETAV,:) = &
              patch%mixt%vals(BCDAT_INFLOW_BETAV,:)*global%rad
          ENDIF

! ----- distribution from external source / constant value

        ELSE
          patch%mixt%vals(BCDAT_INFLOW_PTOT ,:) = vals(3)
          patch%mixt%vals(BCDAT_INFLOW_TTOT ,:) = vals(4)
          patch%mixt%vals(BCDAT_INFLOW_BETAH,:) = vals(5)*global%rad
          patch%mixt%vals(BCDAT_INFLOW_BETAV,:) = vals(6)*global%rad
          IF (switch /= BCOPT_SUBSONIC) THEN
            patch%mixt%vals(BCDAT_INFLOW_MACH,:) = vals(7)
          ENDIF
        ENDIF  ! distribution?

      ENDIF    ! bcType, active region on my processor

    ENDDO      ! iPatch
  ENDDO        ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ReadBcInflowTotAngSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReadBcInflowTotAngSection.F90,v $
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
! Revision 1.2  2006/08/19 15:38:10  mparmar
! Renamed patch variables
!
! Revision 1.1  2005/04/28 05:48:28  wasistho
! added velocity based inflow BC
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.7  2004/01/29 22:55:30  haselbac
! Read new argument
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
! Revision 1.8  2002/06/22 01:13:38  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.7  2002/03/29 23:15:22  jblazek
! Corrected bug in MPI send.
!
! Revision 1.6  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.5  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
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







