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
! Purpose: read in user input for PEUL information (done on all processors).
!
! Description: none.
!
! Input: user input file.
!
! Output: regions = PEUL information
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_ReadInputFile.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_ReadInputFile( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters

  USE PEUL_ModInterfaces, ONLY : PEUL_ReadConPartSection, &
                                 PEUL_ReadConPartPtypeSection
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iRead, iReg

! ... local variables
  CHARACTER(CHRLEN)   :: RCSIdentString
  CHARACTER(CHRLEN+4) :: fname
  CHARACTER(256)      :: line

  LOGICAL :: usedSomewhere, unusedSomewhere

  INTEGER :: errorFlag, iPtype, nPtypes, nSecs, brbeg, brend, readStatus

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_ReadInputFile.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PEUL_ReadInputFile',&
  'PEUL_ReadInputFile.F90' )

! begin -----------------------------------------------------------------------

! search for CONPART and CONPART_PTYPE sections

  fname = TRIM(global%inDir)//TRIM(global%casename)//'.inp'

  DO iRead = 1,2

! - open file

    OPEN(IF_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! - read file looking for keywords

    SELECT CASE (iRead)

    CASE (1) ! on first pass, find # of sections and Ptypes

      nSecs = 0 ! initialize counter for multiply referenced sections

      DO
        READ(IF_INPUT,'(A256)',err=10,end=86) line

        SELECT CASE(TRIM(line))

        CASE ('# CONPART')
          IF (nSecs == 1) THEN
            nPtypes = iPtype
          ELSE IF (nSecs > 1) THEN
            IF (iPtype /= nPtypes) & ! nPtypes same in each sec?
              CALL ErrorStop( global,ERR_PEUL_NPTYPES,__LINE__ )
          END IF ! nSecs
          nSecs = nSecs + 1
          iPtype = 0

        CASE ('# CONPART_PTYPE')
          IF (nSecs == 0) &
            CALL ErrorStop( global,ERR_PEUL_PTYPE,__LINE__ )
          iPtype = iPtype + 1

        END SELECT ! line
      END DO

86    CONTINUE

      IF (nSecs == 1) THEN
        nPtypes = iPtype
      ELSE IF (nSecs > 1) THEN
        IF (iPtype /= nPtypes) & ! nPtypes same in each sec?
          CALL ErrorStop( global,ERR_PEUL_NPTYPES,__LINE__ )
      END IF ! nSecs

    CASE (2) ! on second pass, search for real- and string-valued keys

      DO
        READ(IF_INPUT,'(A256)',err=10,end=87) line

        SELECT CASE(TRIM(line))

        CASE ('# CONPART')
          CALL PEUL_ReadConPartSection( regions,nPtypes,brbeg,brend )
          iPtype = 0

        CASE ('# CONPART_PTYPE')
          iPtype = iPtype + 1
          CALL PEUL_ReadConPartPtypeSection( regions,brbeg,brend,iPtype )

        END SELECT ! line
      END DO

87    CONTINUE

    END SELECT ! iRead

! close file ------------------------------------------------------------------

    CLOSE(IF_INPUT,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )

    IF (nSecs == 0) EXIT

  END DO ! iRead

! set global%peulUsed -------------------------------------------------------

  usedSomewhere   = .FALSE.
  unusedSomewhere = .FALSE.

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
    readStatus = regions(iReg)%peulInput%readStatus
    usedSomewhere   = usedSomewhere  .OR.(readStatus == 1)
    unusedSomewhere = unusedSomewhere.OR.(readStatus /= 1) ! == 0 or -1
  END DO ! iReg

  IF (usedSomewhere.AND.unusedSomewhere) THEN
    CALL ErrorStop( global,ERR_MP_ALLORNONE,__LINE__ )
  END IF ! usedSomewhere.AND.unusedSomewhere

  global%peulUsed = usedSomewhere

! finalization & error handling -----------------------------------------------

  CALL DeregisterFunction( global )
  GOTO 999

10  CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) )

999 CONTINUE

END SUBROUTINE PEUL_ReadInputFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_ReadInputFile.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:47  haselbac
! Initial revision after changing case
!
! Revision 1.7  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2003/04/14 16:32:24  jferry
! minor edits
!
! Revision 1.5  2003/03/17 18:42:45  jblazek
! Added inDir to path of the input file.
!
! Revision 1.4  2003/03/11 16:04:57  jferry
! Created data type for material properties
!
! Revision 1.3  2003/02/26 23:38:30  jferry
! eliminated end=999 convention to ensure that files get closed
!
! Revision 1.2  2003/02/12 23:34:48  jferry
! Replaced [io]stat=global%error with local errorFlag
!
! Revision 1.1  2003/02/11 22:52:51  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







