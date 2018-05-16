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
! Purpose: read user input and store it in the data structure.
!          Copy region topology to coarser grids.
!          Check user input.
!
! Description: none.
!
! Input: regions = dimensions and topology (finest grid).
!
! Output: regions = dimensions, topology and user input (finest grid).
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_UserInput.F90,v 1.7 2008/12/06 08:44:28 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_UserInput( regions )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_InitInputValues, RFLO_ReadRegionmapSection, &
        RFLO_MapRegionsProcessors, ReadInputFile, RFLO_DerivedInputValues, &
        RFLO_CheckUserInput, RFLO_ReadBcInputFile
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  CHARACTER(CHRLEN) :: msg

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_UserInput',&
  'RFLO_UserInput.F90' )

! initialize parameters

  CALL RFLO_InitInputValues( regions )

! read & do region mapping

  CALL RFLO_ReadRegionmapSection( global )

  CALL RFLO_MapRegionsProcessors( regions )

! read user input

  CALL ReadInputFile( regions )


! change noslip to slip wall for Euler

  DO iReg=1,global%nRegions
    IF (regions(iReg)%mixtInput%flowModel == FLOW_EULER) THEN
      DO iPatch=1,regions(iReg)%nPatches
        patch => regions(iReg)%levels(1)%patches(iPatch)   ! only finest level
        IF (patch%bcType>=BC_NOSLIPWALL .AND. &            ! defined yet
            patch%bcType<=BC_NOSLIPWALL+BC_RANGE) THEN
          patch%bcType        = BC_SLIPWALL
          patch%mixt%bcSet = .false.
        ENDIF
      ENDDO
    ENDIF
  ENDDO

! read boundary conditions

  CALL RFLO_ReadBcInputFile( regions )

! set model & numerical parameters from user input

  CALL RFLO_DerivedInputValues( regions )

! check input (other than BCs, for which all modules are checked later)

  CALL RFLO_CheckUserInput( regions )

! set start and current grid levels

  DO iReg=1,global%nRegions
    IF (regions(iReg)%nGridLevels < global%startLevel) THEN
      WRITE(msg,1000) iReg,regions(iReg)%nGridLevels
      CALL ErrorStop( global,ERR_GRID_LEVEL,__LINE__,msg )
    ELSE
      regions(iReg)%startLevel = global%startLevel
      regions(iReg)%currLevel  = global%startLevel
    ENDIF
  ENDDO

! finalize

  CALL DeregisterFunction( global )

1000 FORMAT('Region ',I5,', grid level= ',I2,'.')

END SUBROUTINE RFLO_UserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_UserInput.F90,v $
! Revision 1.7  2008/12/06 08:44:28  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:39  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2008/10/23 18:20:57  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.4  2006/08/19 15:39:49  mparmar
! Renamed patch variables
!
! Revision 1.3  2005/10/17 22:34:30  wasistho
! moved calcCellCtr test to from UserInput to DerivedInputValues
!
! Revision 1.2  2005/05/21 05:42:34  wasistho
! set calcCellCtr
!
! Revision 1.1  2004/11/29 20:51:40  wasistho
! lower to upper case
!
! Revision 1.15  2003/11/20 16:40:40  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.11  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.10  2003/02/11 22:49:53  jferry
! Re-worked BC and TBC input routines to add multi-physics capability
!
! Revision 1.1  2003/02/11 22:30:21  jferry
! Re-worked BC and TBC input routines to add multi-physics capability
!
! Revision 1.9  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.8  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.7  2002/02/27 18:38:20  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
! Revision 1.6  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.5  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.4  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.3  2002/01/02 16:04:20  jblazek
! Added routines to generate geometry for dummy cells.
!
! Revision 1.2  2001/12/08 00:18:42  jblazek
! Added routines to read BC input file.
!
! Revision 1.1.1.1  2001/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************







