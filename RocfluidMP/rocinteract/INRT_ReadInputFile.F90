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
! Purpose: read in user input for INRT information (done on all processors).
!
! Description: none.
!
! Input: user input file.
!
! Output: regions = INRT information
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_ReadInputFile.F90,v 1.3 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_ReadInputFile( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModParameters

  USE INRT_ModInterfaces, ONLY : INRT_ReadDefaultSection,INRT_ReadDrag, &
    INRT_ReadHeatTransferNonBurn,INRT_ReadScouring, INRT_ReadBurning,   &
    INRT_ReadBoilingRegulation
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iRead,iReg

! ... local variables
  CHARACTER(CHRLEN)   :: RCSIdentString
  CHARACTER(CHRLEN+4) :: fname
  CHARACTER(256)      :: line

  LOGICAL :: usedSomewhere, unusedSomewhere

  INTEGER :: errorFlag

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_ReadInputFile.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'INRT_ReadInputFile',&
  'INRT_ReadInputFile.F90' )

! begin -----------------------------------------------------------------------

  fname = TRIM(global%inDir)//TRIM(global%casename)//'.inp'

  DO iRead = 1,2

! - open file

    OPEN(IF_INPUT,file=fname,form='formatted',status='old',iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! - read file looking for keywords

    SELECT CASE (iRead)

    CASE (1) ! on first pass, look for INRT_DEFAULT sections

      DO
        READ(IF_INPUT,'(A256)',err=10,end=86) line

        SELECT CASE(TRIM(line))

        CASE ('# INRT_DEFAULT')
          CALL INRT_ReadDefaultSection( regions )

        END SELECT ! line
      END DO

86    CONTINUE

    CASE (2) ! on second pass, look for all other INRT sections

      DO
        READ(IF_INPUT,'(A256)',err=10,end=87) line

        SELECT CASE(TRIM(line))

        CASE ('# INRT_DRAG')
          CALL INRT_ReadDrag( regions )

        CASE ('# INRT_HEAT_TRANSFER_NONBURN')
          CALL INRT_ReadHeatTransferNonBurn( regions )

        CASE ('# INRT_SCOURING')
          CALL INRT_ReadScouring( regions )

        CASE ('# INRT_BURNING')
          CALL INRT_ReadBurning( regions )

        CASE ('# INRT_BOILING_REGULATION')
          CALL INRT_ReadBoilingRegulation( regions )

        END SELECT ! line
      END DO

87    CONTINUE

    CASE DEFAULT
      CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

    END SELECT ! iRead

! - close file ----------------------------------------------------------------

    CLOSE(IF_INPUT,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )

  END DO ! iRead

! set global%inrtUsed -------------------------------------------------------

  usedSomewhere   = .FALSE.
  unusedSomewhere = .FALSE.

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
    usedSomewhere   = usedSomewhere   .OR. &
                             regions(iReg)%inrtInput%defaultRead
    unusedSomewhere = unusedSomewhere .OR. &
                        .NOT.regions(iReg)%inrtInput%defaultRead
  END DO ! iReg

  IF (usedSomewhere.AND.unusedSomewhere) THEN
    CALL ErrorStop( global,ERR_MP_ALLORNONE,__LINE__ )
  END IF ! usedSomewhere.AND.unusedSomewhere

  global%inrtUsed = usedSomewhere

! finalization & error handling -----------------------------------------------

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname) )

999  CONTINUE

END SUBROUTINE INRT_ReadInputFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_ReadInputFile.F90,v $
! Revision 1.3  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:56:37  fnajjar
! Initial revision after changing case
!
! Revision 1.9  2004/07/28 15:42:13  jferry
! deleted defunct constructs: useDetangle, useSmokeDrag, useSmokeHeatTransfer
!
! Revision 1.8  2004/07/23 22:43:17  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.7  2004/04/15 16:04:21  jferry
! minor formatting (removed trailing spaces)
!
! Revision 1.6  2004/03/02 21:48:09  jferry
! First phase of replacing Detangle interaction
!
! Revision 1.5  2003/09/25 15:48:43  jferry
! implemented Boiling Regulation interaction
!
! Revision 1.4  2003/04/02 22:32:04  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.3  2003/03/24 23:30:52  jferry
! overhauled rocinteract to allow interaction design to use user input
!
! Revision 1.2  2003/03/11 16:09:39  jferry
! Added comments
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************







