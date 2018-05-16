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
! Purpose: read in user input related to boundary conditions
!          (done on all processors).
!
! Description: none.
!
! Input: boundary condition file.
!
! Output: regions = BC data for all regions.
!
! Notes: 
!
!******************************************************************************
!
! $Id: RFLO_ReadBcInputFile.F90,v 1.7 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReadBcInputFile( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_ReadBcNoslipSection, &
        RFLO_ReadBcInflowTotAngSection, RFLO_ReadBcInflowVelSection, &
        RFLO_ReadBcOutflowSection, RFLO_ReadBcFarfSection, &
        RFLO_ReadBcInjectMrateSection, RFLO_ReadBcInjectAPNSection, &
        RFLO_ReadBcSlipWallSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  CHARACTER(2*CHRLEN+9) :: fname
  CHARACTER(256)        :: line

  INTEGER :: errorFlag

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ReadBcInputFile',&
  'RFLO_ReadBcInputFile.F90' )

! - open file

  WRITE(fname,'(A)') TRIM(global%inDir)//TRIM(global%casename)//'.bc'
  OPEN(IF_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! - read file looking for keywords

  DO
    READ(IF_INPUT,'(A256)',err=10,end=86) line
    SELECT CASE(TRIM(line))

    CASE ('# BC_SLIPW')
      CALL RFLO_ReadBcSlipWallSection( regions )

    CASE ('# BC_NOSLIP')
      CALL RFLO_ReadBcNoslipSection( regions )

! TEMPORARY - Keep this for backward compatibility
    CASE ('# BC_INFLOW')
      CALL RFLO_ReadBcInflowTotAngSection( regions )
! END TEMPORARY
    
    CASE ('# BC_INFLOW_TOTANG')
      CALL RFLO_ReadBcInflowTotAngSection( regions )

    CASE ('# BC_INFLOW_VELTEMP')
      CALL RFLO_ReadBcInflowVelSection( regions,BC_INFLOW_VELTEMP )                    

    CASE ('# BC_INFLOW_VELPRESS')
      CALL RFLO_ReadBcInflowVelSection( regions,BC_INFLOW_VELPRESS )                    

    CASE ('# BC_OUTFLOW')
      CALL RFLO_ReadBcOutflowSection( regions )
    
    CASE ('# BC_FARF')
      CALL RFLO_ReadBcFarfSection( regions )
    
! TEMPORARY - Keep this for backward compatibility
    CASE ('# BC_INJECT')
      CALL RFLO_ReadBcInjectMrateSection( regions )
! END TEMPORARY

   CASE ('# BC_INJECT_MRATE')
      CALL RFLO_ReadBcInjectMrateSection( regions )

    CASE ('# BC_INJECT_APN')
      CALL RFLO_ReadBcInjectAPNSection( regions )

    END SELECT
  ENDDO

86   CONTINUE

! close file ------------------------------------------------------------------
  CLOSE(IF_INPUT,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )

! finalization & error handling -----------------------------------------------

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) )

999  CONTINUE

END SUBROUTINE RFLO_ReadBcInputFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReadBcInputFile.F90,v $
! Revision 1.7  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2008/10/23 18:20:53  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.4  2006/01/20 06:14:37  wasistho
! added ReadBcInjectMrate and ReadBcInjectAPN
!
! Revision 1.3  2005/04/28 22:06:56  wasistho
! fixed ModInterfaces list
!
! Revision 1.2  2005/04/28 05:45:37  wasistho
! added velocity based inflow BC
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.7  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.4  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.3  2003/02/26 23:38:30  jferry
! eliminated end=999 convention to ensure that files get closed
!
! Revision 1.2  2003/02/11 22:30:21  jferry
! Re-worked BC and TBC input routines to add multi-physics capability
!
! Revision 1.1  2002/10/19 00:40:30  jblazek
! Added utility (rflosurf) to write out surface grids for GenX.
!
! Revision 1.9  2002/10/14 22:11:53  jblazek
! No more number of regions in the name of the BC file.
!
! Revision 1.8  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.7  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.6  2002/09/17 13:43:00  jferry
! Added Time-dependent boundary conditions
!
! Revision 1.5  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/02/27 18:38:20  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
! Revision 1.3  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.2  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.1  2001/12/08 00:18:42  jblazek
! Added routines to read BC input file.
!
!******************************************************************************







