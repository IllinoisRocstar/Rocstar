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
! Purpose: read in user input related to region mapping and activation.
!
! Description: none.
!
! Input: user input file.
!
! Output: global = type of mapping regions to processors.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ReadRegionMapSection.F90,v 1.5 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReadRegionMapSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  CHARACTER(2*CHRLEN+4) :: fname
  CHARACTER(10)         :: keys(1)
  CHARACTER(256)        :: line

  INTEGER :: errorFlag

  LOGICAL :: defined(1)

  REAL(RFREAL) :: vals(1)

!******************************************************************************

  CALL RegisterFunction( global,'RFLO_ReadRegionMapSection',&
  'RFLO_ReadRegionMapSection.F90' )

! open file and search for keywords

  keys(1) = 'NBLOCKS'

  fname = TRIM(global%inDir)//TRIM(global%casename)//'.inp'
  OPEN(IF_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,fname )

  DO
    READ(IF_INPUT,'(A256)',err=10,end=86) line
    IF (TRIM(line) == '# BLOCKMAP') THEN
      CALL ReadSection( global,IF_INPUT,1,keys,vals,defined )
      IF (defined(1).eqv. .true.) global%nRegionsProc = NINT(vals(1))
      EXIT
    ENDIF
  ENDDO

86   CONTINUE

! close file ------------------------------------------------------------------

  CLOSE(IF_INPUT,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,fname )

! finalization & error handling -----------------------------------------------

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__,fname )

999  CONTINUE

END SUBROUTINE RFLO_ReadRegionMapSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReadRegionMapSection.F90,v $
! Revision 1.5  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2008/10/23 18:20:57  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.2  2005/01/12 04:19:11  wasistho
! applied NINT to vals(1)
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.10  2003/11/20 16:40:40  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.5  2003/02/26 23:38:30  jferry
! eliminated end=999 convention to ensure that files get closed
!
! Revision 1.4  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.3  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.2  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.3  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.2  2001/12/22 00:09:39  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1.1.1  2001/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************







