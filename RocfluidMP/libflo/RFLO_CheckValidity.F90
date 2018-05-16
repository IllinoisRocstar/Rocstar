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
! Purpose: Check validity of variables
!
! Description: It detects NaN solution and/or negative pressure or density.
!
! Input:  region = Region data
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_CheckValidity.F90,v 1.6 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CheckValidity(region)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModParameters
  USE ModMPI
  USE ModTools, ONLY: IsNan
  USE ModInterfaces, ONLY: MixtPerf_G_CpR, MixtPerf_P_DEoGVm2, MixtPerf_R_M, &
                           MixtPerf_T_DPR, &
                           RFLO_GetDimensPhys, RFLO_GetCellOffset
  
  IMPLICIT NONE
#include "Indexing.h"

! ... parameters
  TYPE(t_region), TARGET :: region

! ... loop variables
  INTEGER :: i, j, k, iPatch

! ... local variables
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patches(:)

  INTEGER, PARAMETER :: MAX_INVALID_LOCS = 10
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, iCOff, ijCOff, ijkC
  INTEGER :: indCp,indMol,nLocs
#ifdef MPI
  INTEGER :: MPIerrCode
#endif
  REAL(RFREAL) :: Eo,gamma,p,rgas,rho,rrho,t,u,v,Vm2,w,rmin,pmin
  REAL(RFREAL), DIMENSION(:,:), POINTER :: cv,gv
  LOGICAL :: foundNan

!******************************************************************************

  global => region%global

  CALL RegisterFunction(global,'RFLO_CheckValidity',&
  'RFLO_CheckValidity.F90')

  nLocs = 0

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  cv      => region%levels(iLev)%mixt%cv
  gv      => region%levels(iLev)%mixt%gv
  patches => region%levels(iLev)%patches

  indCp  = region%levels(iLev)%mixt%indCp
  indMol = region%levels(iLev)%mixt%indMol 

! loop over cells and check for positivity ------------------------------------

  rmin = 1.E+13_RFREAL
  pmin = 1.E+13_RFREAL
  foundNan = .FALSE.

  DO k=kpcbeg,kpcend
    DO j=jpcbeg,jpcend
      DO i=ipcbeg,ipcend
        ijkC = IndIJK(i,j,k,iCOff,ijCOff)

        rho  = cv(CV_MIXT_DENS,ijkC)
        rrho = 1.0_RFREAL/rho
        u    = rrho*cv(CV_MIXT_XMOM,ijkC)
        v    = rrho*cv(CV_MIXT_YMOM,ijkC)
        w    = rrho*cv(CV_MIXT_ZMOM,ijkC)        
        Eo   = rrho*cv(CV_MIXT_ENER,ijkC)

        rgas  = MixtPerf_R_M(gv(GV_MIXT_MOL,ijkC*indMol))
        gamma = MixtPerf_G_CpR(gv(GV_MIXT_CP,ijkC*indCp),rgas)
        Vm2   = u*u + v*v + w*w

        p = MixtPerf_P_DEoGVm2(rho,Eo,gamma,Vm2)

        rmin  = MIN( rmin,rho )
        pmin  = MIN( pmin,p )

        foundNan = IsNan(rho)
        IF (foundNan) GOTO 888

      ENDDO      ! i
    ENDDO        ! j
  ENDDO          ! k

888 CONTINUE

  IF (foundNan .OR. rmin < 0._RFREAL .OR. pmin < 0._RFREAL) THEN

    DO k=kpcbeg,kpcend
     DO j=jpcbeg,jpcend
      DO i=ipcbeg,ipcend
        ijkC = IndIJK(i,j,k,iCOff,ijCOff)

        rho  = cv(CV_MIXT_DENS,ijkC)
        rrho = 1.0_RFREAL/rho
        u    = rrho*cv(CV_MIXT_XMOM,ijkC)
        v    = rrho*cv(CV_MIXT_YMOM,ijkC)
        w    = rrho*cv(CV_MIXT_ZMOM,ijkC)        
        Eo   = rrho*cv(CV_MIXT_ENER,ijkC)

        rgas  = MixtPerf_R_M(gv(GV_MIXT_MOL,ijkC*indMol))
        gamma = MixtPerf_G_CpR(gv(GV_MIXT_CP,ijkC*indCp),rgas)
        Vm2   = u*u + v*v + w*w

        p = MixtPerf_P_DEoGVm2(rho,Eo,gamma,Vm2)
        t = MixtPerf_T_DPR(rho,p,rgas)

        IF ( (IsNan(rho) .EQV. .TRUE.) .OR. & 
             (IsNan(u)   .EQV. .TRUE.) .OR. &
             (IsNan(v)   .EQV. .TRUE.) .OR. & 
             (IsNan(w)   .EQV. .TRUE.) .OR. & 
             (IsNan(p)   .EQV. .TRUE.) .OR. &
             (IsNan(t)   .EQV. .TRUE.) .OR. &
             (      rho  <  0._RFREAL) .OR. &
             (      p    <  0._RFREAL) ) THEN
          nLocs = nLocs + 1   

          IF ( nLocs == 1 ) THEN 
            WRITE(STDOUT,'(A,1X,A,1X,I9)') SOLVER_NAME, & 
                  'Invalid variables detected!'
              
            IF ( global%flowType == FLOW_UNSTEADY ) THEN 
              WRITE(STDOUT,'(A,3X,A,1X,1PE12.5)') SOLVER_NAME, &
                                                  'Current time:', &
                                                  global%currentTime              
            ELSE 
              WRITE(STDOUT,'(A,3X,A,1X,I6)') SOLVER_NAME, &
                                             'Current iteration number:', &
                                             global%currentIter           
            END IF ! global%flowType                 
                                            
            WRITE(STDOUT,1000) SOLVER_NAME,'Region:',region%iRegionGlobal, & 
                        ', bc-types:', (patches(iPatch)%bcType, &
                          iPatch = 1,region%nPatches)
                                  
            WRITE(STDOUT,'(A,6X,A,7(1X,A))') SOLVER_NAME,'#', &
                                             '   Density   ', &
                                             '  x-velocity ', &
                                             '  y-velocity ', &
                                             '  z-velocity ', &
                                             '   Pressure  ', &
                                             ' Temperature ', &
                                             '   i   ,   j   ,   k'       
          END IF ! nLocs

          IF ( nLocs <= MAX_INVALID_LOCS ) THEN 
            WRITE(STDOUT,'(A,4X,I3,6(1X,E13.6),2X,3I7)') SOLVER_NAME,nLocs, & 
                                                  rho,u,v,w,p,t,i,j,k
          END IF ! nLocs
        END IF   ! dv    
      ENDDO      ! i
     ENDDO       ! j
    ENDDO        ! k

#ifdef MPI
    CALL MPI_Abort( global%mpiComm,MPIerrCode,global%mpierr )
#endif
    STOP
  ENDIF          ! foundNan...

! finalize --------------------------------------------------------------------

1000 FORMAT( A,3X,A,I5,A,60I5 )

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_CheckValidity

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CheckValidity.F90,v $
! Revision 1.6  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/05/12 21:31:14  wasistho
! write i,j,k in the same line as previous vars
!
! Revision 1.3  2005/05/12 20:53:51  wasistho
! added bctypes in error msg
!
! Revision 1.2  2005/01/11 00:26:15  wasistho
! changed mpi_finalize to mpi_abort
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.6  2004/08/24 01:02:23  wasistho
! search nan per cell i.o. the sum value
!
! Revision 1.5  2004/08/16 17:05:52  wasistho
! moved MPI finalize within IF statement
!
! Revision 1.4  2004/08/04 00:29:35  wasistho
! speedup check validity
!
! Revision 1.3  2004/07/26 20:03:13  wasistho
! added i,j,k locations
!
! Revision 1.2  2004/07/26 19:30:44  wasistho
! changed POINTER to TARGET for region parameter
!
! Revision 1.1  2004/07/26 19:10:06  wasistho
! initial import RFLO_CheckValidity
!
!
!******************************************************************************







