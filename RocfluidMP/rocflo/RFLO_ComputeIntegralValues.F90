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
! Purpose: Compute integrals for GENx checking.
!
! Description: None.
!
! Input:
!   regions     Region data
!
! Output:
!   integ       Vector of integrals (for output by Rocman)
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLO_ComputeIntegralValues.F90,v 1.4 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

#ifdef GENX
SUBROUTINE RFLO_ComputeIntegralValues(regions,integ)
#else
SUBROUTINE RFLO_ComputeIntegralValues(regions)
#endif

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  USE ModMPI
  
  IMPLICIT NONE
#include "Indexing.h"

#ifdef GENX
  INCLUDE 'rocmanf90.h'
#endif

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, i, j, k

! ... local variables
  INTEGER :: iLev, ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iCOff, ijCOff, ijkC

  REAL(RFREAL) :: enerLocal,ibAreaLocal,inbAreaLocal,massLocal,xMomLocal, & 
                  yMomLocal,zMomLocal,volLocal
#ifdef GENX
  DOUBLE PRECISION, DIMENSION(MAN_INTEG_SIZE) :: integ
  REAL(RFREAL), DIMENSION(MAN_INTEG_SIZE) :: globalVals,localVals
#else
  REAL(RFREAL), DIMENSION(2) :: globalVals,localVals
#endif
  REAL(RFREAL), POINTER :: cv(:,:), vol(:)

  TYPE(t_global), POINTER :: global

! *****************************************************************************

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLO_ComputeIntegralValues',&
  'RFLO_ComputeIntegralValues.F90')

! get dimensions and pointers -------------------------------------------------

! Compute total volume and total mass -----------------------------------------

  volLocal  = 0.0_RFREAL
  massLocal = 0.0_RFREAL

  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &  ! region active and
        regions(iReg)%active==ACTIVE) THEN             ! on my processor

      iLev = regions(iReg)%currLevel

      CALL RFLO_GetDimensPhys( regions(iReg),iLev,ipcbeg,ipcend, &
                               jpcbeg,jpcend,kpcbeg,kpcend )
      CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )

      cv  => regions(iReg)%levels(iLev)%mixt%cv
      vol => regions(iReg)%levels(iLev)%grid%vol
    
! --- loop over cells ---------------------------------------------------------

      DO k=kpcbeg,kpcend
        DO j=jpcbeg,jpcend
          DO i=ipcbeg,ipcend
            ijkC = IndIJK(i,j,k,iCOff,ijCOff)
            volLocal  = volLocal + vol(ijkC)
            massLocal = massLocal + cv(CV_MIXT_DENS,ijkC)*vol(ijkC)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDDO ! iReg

#ifdef GENX
! Compute momenta components and energy ---------------------------------------

  xMomLocal = 0.0_RFREAL ! Not computed at present
  yMomLocal = 0.0_RFREAL
  zMomLocal = 0.0_RFREAL
  enerLocal = 0.0_RFREAL 

! Compute interacting surface areas -------------------------------------------

  inbAreaLocal = 0.0_RFREAL ! Not computed at present 
  ibAreaLocal  = 0.0_RFREAL
#endif 

! *****************************************************************************
! Gather data
! *****************************************************************************

#ifdef GENX
  localVals(MAN_INTEG_VOL    ) = volLocal
  localVals(MAN_INTEG_MASS   ) = massLocal
  localVals(MAN_INTEG_XMOM   ) = xMomLocal
  localVals(MAN_INTEG_YMOM   ) = yMomLocal
  localVals(MAN_INTEG_ZMOM   ) = zMomLocal
  localVals(MAN_INTEG_ENER   ) = enerLocal
  localVals(MAN_INTEG_IBAREA ) = ibAreaLocal  
  localVals(MAN_INTEG_INBAREA) = inbAreaLocal      
#else
  localVals(1) = volLocal
  localVals(2) = massLocal  
#endif

! Perform reduction operation -------------------------------------------------

#ifdef MPI
  CALL MPI_AllReduce( localVals,globalVals,SIZE(localVals),MPI_RFREAL,MPI_SUM,&
                      global%mpiComm,global%mpierr )
  IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#else
  DO i = 1,SIZE(localVals)
    globalVals(i) = localVals(i)
  END DO ! i
#endif   

! Assign data -----------------------------------------------------------------

#ifdef GENX
  DO i = 1,MAN_INTEG_SIZE
    integ(i) = globalVals(i)
  END DO ! i
#else
  global%totalVol  = globalVals(1)
  global%totalMass = globalVals(2)  
#endif

! Finalize --------------------------------------------------------------------

  CALL DeregisterFunction(global)

END SUBROUTINE RFLO_ComputeIntegralValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ComputeIntegralValues.F90,v $
! Revision 1.4  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:21  haselbac
! Removed tabs
!
! Revision 1.1  2005/02/26 04:07:17  wasistho
! added RFLO_ComputeIntegralValues
!
!
!******************************************************************************







