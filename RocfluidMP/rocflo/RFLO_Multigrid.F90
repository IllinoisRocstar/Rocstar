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
! Purpose: integrate the governing equations in time using
!          the multigrid methodology.
!
! Description: none.
!
! Input: regions = data of all grid regions.
!
! Output: regions = flow variables for all grid regions.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_Multigrid.F90,v 1.3 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_Multigrid( dIterSystem,regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : 
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  INTEGER :: dIterSystem

  TYPE(t_region), POINTER :: regions(:)

! ... loop variables


! ... local variables
  INTEGER :: iter

  LOGICAL :: stopExists

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_Multigrid',&
  'RFLO_Multigrid.F90' )

  iter = 0

! iterations ------------------------------------------------------------------

  DO

! - update time / iteration number

    global%currentIter = global%currentIter + 1
    iter = iter + 1

! - check for stop file 

    INQUIRE(file="STOP",exist=stopExists)
    IF (stopExists) THEN 
      IF (global%myProcid==MASTERPROC .AND. global%verbLevel/=VERBOSE_NONE) &
        WRITE(STDOUT,1000) SOLVER_NAME
    ENDIF

! - check for end of time stepping

    IF (iter==dIterSystem .OR. global%residual<=global%resTol .OR. &
        stopExists) EXIT

  ENDDO  ! loop over iter

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1000  FORMAT(/,A,' File ''STOP'' detected !')

END SUBROUTINE RFLO_Multigrid

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_Multigrid.F90,v $
! Revision 1.3  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.10  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.5  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.3  2002/06/22 01:13:38  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.2  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/16 22:03:35  jblazek
! Added time-stepping routines.
!
!******************************************************************************







