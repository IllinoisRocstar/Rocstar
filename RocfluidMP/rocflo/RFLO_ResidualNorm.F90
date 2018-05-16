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
! Purpose: compute the 2-norm of the residual.
!
! Description: none.
!
! Input: blck = data of current region.
!
! Output: global%residual = norm (all regions, all processors).
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ResidualNorm.F90,v 1.3 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ResidualNorm( regions )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, i, j, k

! ... local variables
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, ijkC, iCOff, ijCOff

  REAL(RFREAL)          :: dr, drho, drhoTot
  REAL(RFREAL), POINTER :: cv(:,:), cvOld(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_ResidualNorm',&
  'RFLO_ResidualNorm.F90' )

! sum drho for all regions on this processor ----------------------------------

  drho = 0._RFREAL

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &     ! region active and
        regions(iReg)%active==ACTIVE) THEN                ! on my processor

! --- get dimensions and pointers

      iLev = regions(iReg)%currLevel

      CALL RFLO_GetDimensPhys( regions(iReg),iLev,ipcbeg,ipcend, &
                               jpcbeg,jpcend,kpcbeg,kpcend )
      CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )

      cv    => regions(iReg)%levels(iLev)%mixt%cv
      cvOld => regions(iReg)%levels(iLev)%mixt%cvOld

! --- add to the sum

      DO k=kpcbeg,kpcend
        DO j=jpcbeg,jpcend
          DO i=ipcbeg,ipcend
            ijkC = IndIJK(i,j,k,iCOff,ijCOff)
            dr   = cv(CV_MIXT_DENS,ijkC) - cvOld(CV_MIXT_DENS,ijkC)
            drho = drho + dr*dr
          ENDDO
        ENDDO
      ENDDO

    ENDIF   ! active
  ENDDO     ! iReg

! exchange between processors -------------------------------------------------

#ifdef MPI
  CALL MPI_Allreduce( drho,drhoTot,1,MPI_RFREAL,MPI_SUM, &
                      global%mpiComm,global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#else
  drhoTot = drho
#endif

  global%residual = SQRT(drhoTot)

  IF (global%currentIter == 1) THEN
    global%resInit = global%residual
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_ResidualNorm

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ResidualNorm.F90,v $
! Revision 1.3  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:40  wasistho
! lower to upper case
!
! Revision 1.10  2003/11/20 16:40:40  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.5  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.4  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.3  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.2  2002/01/31 20:56:30  jblazek
! Added basic boundary conditions.
!
! Revision 1.1  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
!******************************************************************************







