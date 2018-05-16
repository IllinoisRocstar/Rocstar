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
! Purpose: read state of random number generator for structured code
!
! Description: the following solution formats are supported:
!              - RocfloMP ASCII
!              - RocfloMP binary.
!
! Input: regions = dimensions of all regions.
!
! Output: regions(:)%randData = state of RNG for each region
!         global%currentTime  = physical time
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ReadRandomState.F90,v 1.3 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ReadRandomState( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  USE ModRandom

  USE ModInterfaces, ONLY : RFLO_ReadDataFileInt, RFLO_ReadDataFileReal
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, i

! ... local variables
  INTEGER, PARAMETER     :: N_RVAR = 1
  INTEGER, PARAMETER     :: N_IVAR_EXTRA = 1
  INTEGER, PARAMETER     :: RANDSTATE_OFFSET = N_IVAR_EXTRA+1
  INTEGER, PARAMETER     :: N_RANDSTATE = RAND_TOTAL_SIZE+N_IVAR_EXTRA

  CHARACTER(2*CHRLEN+18) :: fname, fhead
  CHARACTER(CHRLEN)      :: RCSIdentString, msg

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: errorFlag, mti, randState(N_RANDSTATE,1)

  REAL(RFREAL) :: rvar(N_RVAR,1)

  LOGICAL :: fileExists

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ReadRandomState',&
  'RFLO_ReadRandomState.F90' )

! begin -----------------------------------------------------------------------

! determine if RNG state file exists ------------------------------------------

  fhead = TRIM(global%inDir)//TRIM(global%casename)

  SELECT CASE(global%solutFormat)
  CASE (FORMAT_ASCII)
    fhead = TRIM(fhead)//'.randa_'
  CASE (FORMAT_BINARY)
    fhead = TRIM(fhead)//'.randb_'
  CASE DEFAULT
    CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
  END SELECT

  IF (global%flowType == FLOW_UNSTEADY) THEN
    WRITE(fname,'(A,ES11.5)') TRIM(fhead), global%timeStamp
  ELSE
    WRITE(fname,'(A,I6.6)')   TRIM(fhead), global%currentIter
  ENDIF

  INQUIRE(file=fname,exist=fileExists)

  IF (.NOT.fileExists) THEN
    IF (global%myProcid == MASTERPROC .AND. &
        global%currentTime > 0._RFREAL) THEN
      WRITE(STDOUT,'(A)') SOLVER_NAME// &
        '### WARNING: No Random State file found: '// &
        'Random Number Generator re-initialized'
    ENDIF ! MASTERPROC
    GO TO 999
  ENDIF ! fileExists

! open random state file (only master proc.) ----------------------------------

  IF (global%myProcid == MASTERPROC) THEN

    SELECT CASE(global%solutFormat)
    CASE (FORMAT_ASCII)
      OPEN(IF_RAND_STATE,file=fname,form=  'formatted',status='old', &
        iostat=errorFlag)
    CASE (FORMAT_BINARY)
      OPEN(IF_RAND_STATE,file=fname,form='unformatted',status='old', &
        iostat=errorFlag)
    CASE DEFAULT
      CALL ErrorStop( global,ERR_UNKNOWN_FORMAT,__LINE__ )
    END SELECT

    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! - Read time stamp of RNG state, and yield an error if not the current time --

    CALL RFLO_ReadDataFileReal( global,IF_RAND_STATE,global%solutFormat, &
                                N_RVAR,1,rvar )

    IF (global%flowType==FLOW_UNSTEADY.AND.global%currentTime/=rvar(1,1)) THEN
      WRITE(msg,1000) rvar(1,1),global%currentTime
      CALL ErrorStop( global,ERR_TIME_SOLUTION,__LINE__,msg//' File: '// &
        TRIM(fname) )
    ENDIF

  ENDIF ! MASTERPROC

! Receive data from master ----------------------------------------------------

  DO iReg=1,global%nRegions

! - Read the RNG state for each region (only master) --------------------------

    IF (global%myProcid == MASTERPROC) THEN

      CALL RFLO_ReadDataFileInt( global,IF_RAND_STATE,global%solutFormat, &
                                 N_RANDSTATE,1,randState )

      IF (randState(1,1) /= iReg) &
        CALL ErrorStop( global,ERR_REGION_NUMBER,__LINE__,'File: '// &
          TRIM(fname) )

! --- Allowable values of mti are within array bounds 0:RAND_BUFFER_SIZE-1,
! --- or the special case mti = RAND_BUFFER_SIZE (generates array refill)

      mti = randState(RAND_MTI_INDEX+RANDSTATE_OFFSET,1)

      IF (mti < 0 .OR. mti > RAND_BUFFER_SIZE) &
        CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

#ifdef MPI
      IF (regions(iReg)%procid /= MASTERPROC) THEN

        CALL MPI_Send( randState,N_RANDSTATE,MPI_INTEGER, &
                       regions(iReg)%procid,iReg, &
                       global%mpiComm,global%mpierr )
        IF (global%mpierr /=0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

      ENDIF

    ELSE ! not the master

      IF (regions(iReg)%procid == global%myProcid) THEN

        CALL MPI_Recv( randState,N_RANDSTATE,MPI_INTEGER,MASTERPROC,iReg, &
                       global%mpiComm,status,global%mpierr )
        IF (global%mpierr /=0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

      ENDIF
#endif

    ENDIF ! MASTERPROC

! - copy temporary array into randData

    IF (regions(iReg)%procid == global%myProcid) THEN
      DO i=0,RAND_TOTAL_SIZE-1
        regions(iReg)%randData%mt(i) = randState(i+RANDSTATE_OFFSET,1)
      ENDDO ! i
    ENDIF ! global%myProcid

  ENDDO ! iReg

! finalize --------------------------------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CLOSE(IF_RAND_STATE,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname) )
  ENDIF

999  CONTINUE
  CALL DeregisterFunction( global )

1000 FORMAT('Time in file is= ',ES12.5,' but it should be= ',ES12.5,'.')

END SUBROUTINE RFLO_ReadRandomState

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ReadRandomState.F90,v $
! Revision 1.3  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:21  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.1  2003/11/21 22:28:42  fnajjar
! Added Read capability for Random Number Generator
!
!******************************************************************************







