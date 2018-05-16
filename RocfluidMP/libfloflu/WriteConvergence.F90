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
! Purpose: write convergence history to screen and to file.
!
! Description: none.
!
! Input: global%forceX/Y/Z = forces in x-, y- and z-direction
!        global%massIn/Out = mass flow in and out of domain
!        global%residual   = actual density residual
!        global%resInit    = initial residual (for normalization purposes)
!
! Output: to screen and file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: WriteConvergence.F90,v 1.4 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE WriteConvergence( global )

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY : t_global
  USE ModTools, ONLY  : MakeNonZero
  USE ModMPI
  USE ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  REAL(RFREAL) :: resNorm
#ifdef MPI
  REAL(RFREAL) :: localVals(6), globalVals(6)
#endif

!******************************************************************************

  CALL RegisterFunction( global,'WriteConvergence',&
  'WriteConvergence.F90' )

! sum up data from other processors -------------------------------------------

#ifdef MPI
  localVals(1) = global%forceX
  localVals(2) = global%forceY
  localVals(3) = global%forceZ
  localVals(4) = global%massIn
  localVals(5) = global%massOut
  localVals(6) = global%stopRun

!  CALL MPI_Reduce( localVals,globalVals,6,MPI_RFREAL,MPI_SUM,MASTERPROC, &
!                   global%mpiComm,global%mpierr )
  CALL MPI_Allreduce( localVals,globalVals,6,MPI_RFREAL,MPI_SUM, &
                      global%mpiComm,global%mpierr )
  IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

  global%forceX  = globalVals(1)
  global%forceY  = globalVals(2)
  global%forceZ  = globalVals(3)
  global%massIn  = globalVals(4)
  global%massOut = globalVals(5)
  global%stopRun = globalVals(6)
#endif

! steady flow -----------------------------------------------------------------

  IF (global%flowType==FLOW_STEADY .AND. global%myProcid==MASTERPROC) THEN
    IF (global%currentIter == 1) THEN
      resNorm = 1._RFREAL
    ELSE
      resNorm = global%residual/MakeNonZero(global%resInit)
      resNorm = MakeNonZero(resNorm)
    ENDIF

    IF (global%verbLevel /= VERBOSE_NONE) &
      WRITE(STDOUT,1000) SOLVER_NAME,global%currentIter,LOG10(resNorm), &
                         global%forceX,global%forceY,global%forceZ, &
                         global%massIn,global%massOut

    WRITE(IF_CONVER,1010,err=10) global%currentIter,LOG10(resNorm), &
                                 global%forceX,global%forceY,global%forceZ, &
                                 global%massIn,global%massOut

! unsteady flow ---------------------------------------------------------------

  ELSE IF (global%flowType==FLOW_UNSTEADY .AND. global%myProcid==MASTERPROC) THEN
    IF (global%verbLevel /= VERBOSE_NONE) &
      WRITE(STDOUT,2000) SOLVER_NAME,global%currentTime,global%dtMin, &
                         global%forceX,global%forceY,global%forceZ, &
                         global%massIn,global%massOut

    WRITE(IF_CONVER,2010,err=10) global%currentTime,global%dtMin, &
                                 global%forceX,global%forceY,global%forceZ, &
                                 global%massIn,global%massOut
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_WRITE,__LINE__,'Convergence file.' )

1000 FORMAT(A,1X,I6,1PE13.4,5E13.4)
1010 FORMAT(     I6,1PE13.4,5E13.4)
2000 FORMAT(A,1X,1PE12.5,6E13.4)
2010 FORMAT(     1PE12.5,6E13.4)

999  CONTINUE

END SUBROUTINE WriteConvergence

!******************************************************************************
!
! RCS Revision history:
!
! $Log: WriteConvergence.F90,v $
! Revision 1.4  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:24  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/04/15 15:06:07  haselbac
! Removed Charm/FEM stuff, routine no longer called by RFLU
!
! Revision 1.1  2004/12/01 16:52:21  haselbac
! Initial revision after changing case
!
! Revision 1.16  2003/05/19 21:18:21  jblazek
! Automated switch to 0th-order extrapolation at slip walls and injection boundaries.
!
! Revision 1.15  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.14  2003/02/06 01:33:38  jblazek
! Taken resNorm out of conditional compilation.
!
! Revision 1.13  2003/02/05 21:07:30  jblazek
! Coordinated stop of a run works now for MPI.
!
! Revision 1.12  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.11  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.10  2002/08/30 01:47:58  jblazek
! Added support for moving grids.
!
! Revision 1.9  2002/07/25 14:49:31  haselbac
! Added FEM call to get force and mass flows for parallel runs
!
! Revision 1.8  2002/06/17 13:32:55  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.7  2002/05/04 16:19:23  haselbac
! Added MakeNonZero function and ModTools
!
! Revision 1.6  2002/03/18 22:25:45  jblazek
! Finished multiblock and MPI.
!
! Revision 1.5  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/02/16 07:16:00  jblazek
! Added implicit residual smoothing.
!
! Revision 1.3  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.2  2002/01/31 20:56:30  jblazek
! Added basic boundary conditions.
!
! Revision 1.1  2002/01/30 02:16:24  jblazek
! Convergence output moved to common library.
!
!******************************************************************************







