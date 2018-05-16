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
! Purpose: shut down Rocflo-MP.
!
! Description: none.
!
! Input: regions = all grid regions.
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_EndFlowSolver.F90,v 1.3 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_EndFlowSolver( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iprobe

! ... local variables
  INTEGER :: errorFlag

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_EndFlowSolver',&
  'RFLO_EndFlowSolver.F90' )

! close convergence, probe & thrust file --------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CLOSE(IF_CONVER,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__ )
    IF (global%thrustType /= THRUST_NONE) THEN
      CLOSE(IF_THRUST,iostat=errorFlag)
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__ )
    ENDIF
  ENDIF

  IF (global%nProbes > 0) THEN
    DO iprobe=1,global%nProbes
      IF (regions(global%probePos(iprobe,1))%procid==global%myProcid .AND. &
          regions(global%probePos(iprobe,1))%active==ACTIVE) THEN
        CLOSE(IF_PROBE+iprobe-1,iostat=errorFlag)
        global%error = errorFlag
        IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__ )
      ENDIF
    ENDDO
  ENDIF

  IF (global%myProcid==MASTERPROC .AND. global%verbLevel/=VERBOSE_NONE) &
    WRITE(STDOUT,'(/,A,/)') SOLVER_NAME//' Program finished.'

! finish MPI ------------------------------------------------------------------

#ifndef GENX
#ifdef MPI
  CALL MPI_Barrier( global%mpiComm,global%mpierr )
  CALL MPI_Finalize( global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif
#endif

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_EndFlowSolver

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_EndFlowSolver.F90,v $
! Revision 1.3  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.9  2003/06/02 17:12:00  jblazek
! Added computation of thrust.
!
! Revision 1.8  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.7  2002/10/25 18:36:47  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.6  2002/09/25 17:57:09  jblazek
! Added call to EndFlowSolver from genx/fluid_finalize.
!
! Revision 1.5  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.4  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.3  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.2  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.1  2002/02/25 22:36:52  jblazek
! Simplified solver initialization routine.
!
!******************************************************************************







