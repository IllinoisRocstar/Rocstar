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
! Purpose: calculate the minimum time step for all regions on all processors,
!          compare it to the imposed time step and store the minimum of both.
!
! Description: none.
!
! Input: regions           = data of all grid regions
!        global%dtImposed = user prescribed time step.
!
! Output: global%dtMin = minimum time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_MinimumTimeStep.F90,v 1.3 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_MinimumTimeStep( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  REAL(RFREAL) :: cfl, dtMin, dtMinTot

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_MinimumTimeStep',&
  'RFLO_MinimumTimeStep.F90' )

! compare imposed and max. stable time step

  dtMin = global%dtImposed

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &     ! region active and
        regions(iReg)%active==ACTIVE) THEN                ! on my processor
      cfl   = regions(iReg)%mixtInput%cfl
      dtMin = MIN(dtMin,global%dtMin*cfl)
    ENDIF   ! active
  ENDDO     ! iReg

  global%dtMin = dtMin

! exhange time step between processors

#ifdef MPI
  CALL MPI_Allreduce( global%dtMin,dtMinTot,1,MPI_RFREAL,MPI_MIN, &
                      global%mpiComm,global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

  global%dtMin = dtMinTot
#endif

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_MinimumTimeStep

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_MinimumTimeStep.F90,v $
! Revision 1.3  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.7  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.6  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.5  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.3  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.2  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.1  2002/01/16 22:03:35  jblazek
! Added time-stepping routines.
!
!******************************************************************************







