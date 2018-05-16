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
! Purpose: Compute the fourth energy subgrid contribution, alpha_4.
!
! Description: The contribution of the energy subgrid alpha_4 is added to 
!              the viscous flux residual. 
!              Alpha_1 is the energy transfer term from resolved to sgs,
!                      has been added in VisFluxEddy of VFluxHybrid. 
!              Alpha_4 is the present turbulence dissipation rate term.
!              Alpha_4 = C_eps.bar[rho].k^(3/2)/delta, with
!                      k = 1/2.tau_ii and 
!                      C_eps = int{alpha_1}dx/int{bar[rho].k^(3/2)/delta)}dx
!              Alpha_2 (effect on the heat conductivity) and
!              Alpha_3 (pressure dilatation term) have been modeled together
!                      in the above viscous flux routines using the 
!                      eddy-diffusivity model (turb. heat flux). 
!
! Input: region  = data of current region 
!
! Output: mixt%diss, containing the alpha_4 contribution
!
! Notes: alpha_1,2,3 has been added to mixt%diss in a previous treatment.
!
!******************************************************************************
!
! $Id: TURB_LesEsgModel4.F90,v 1.7 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_LesEsgModel4( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset

#include "Indexing.h"
#endif
  USE TURB_ModInterfaces, ONLY : TURB_StatCCollector
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region)    :: region

! ... loop variables
  INTEGER :: i, j, k, ijkC

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  REAL(RFREAL) :: oo3,oo12,beta,delFac,delta
  REAL(RFREAL) :: esg1Glo,esg4Glo,cEpsilon,esg4
  REAL(RFREAL), POINTER :: cv(:,:),diss(:,:),trace(:),vol(:)
  REAL(RFREAL), POINTER :: colVar(:,:)

  INTEGER :: ibn,ien
#ifdef RFLO
  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend
  INTEGER :: ipcbeg,ipcend,jpcbeg,jpcend,kpcbeg,kpcend
  INTEGER :: iLev,iCOff,ijCOff
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_LesEsgModel4.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_LesEsgModel4',&
  'TURB_LesEsgModel4.F90' )

! get constants, indices and pointers -----------------------------------------

  oo3    = 1.0_RFREAL/3.0_RFREAL
  oo12   = 1.0_RFREAL/12.0_RFREAL
  beta   = region%mixtInput%betrk(region%irkStep)
  delFac = SQRT(region%turbInput%delFac2)

#ifdef RFLO
  iLev   =  region%currLevel
  cv     => region%levels(ilev)%mixt%cv
  diss   => region%levels(ilev)%mixt%diss
  trace  => region%levels(iLev)%turb%trace
  vol    => region%levels(ilev)%grid%vol

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  ibeg = ipcbeg
  iend = ipcend
  jbeg = jpcbeg
  jend = jpcend
  kbeg = kpcbeg
  kend = kpcend

! global sum of rho*k^(3/2)/delta contribution (global%esg4Sum) has been 
! initialized in TURB_LesRkInit

! proceed computation of the model constant C_epsilon
  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,iend
        ijkC = IndIJK(i ,j ,k ,iCOff,ijCOff)
#endif
#ifdef RFLU
  cv    => region%mixt%cv
  diss  => region%mixt%diss
  trace => region%turb%trace
  vol   => region%grid%vol
  
  DO ijkC = 1,region%grid%nCells
#endif

! ----- obtain the subgrid kinetic energy k == 0.5*rho*tau_ii in trace at cell 
!       centers; rho*tau_ii is the trace of rho*tau_ij tensor; note that in 
!       the formal formulation k == 0.5*tau_ii.

        trace(ijkC) = oo12*trace(ijkC)

! ----- get global weighted sum of alpha_4; store rho*k^(3/2)/delta in trace
        delta = delFac*vol(ijkC)**oo3
        trace(ijkC) = SQRT(trace(ijkC)**3/cv(CV_MIXT_DENS,ijkC))/delta
        global%esg4Sum = global%esg4Sum + trace(ijkC)*vol(ijkC)

#ifdef RFLO
      ENDDO ! i
    ENDDO   ! j
  ENDDO     ! k
#endif
#ifdef RFLU
  ENDDO     ! ijkC
#endif

  esg1Glo = global%esg1Psum
  esg4Glo = global%esg4Psum

  cEpsilon = esg1Glo/(esg4Glo+1.E-16_RFREAL)

! add contribution of the 4th ESG model into diss

#ifdef RFLO
  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,iend
        ijkC = IndIJK(i ,j ,k ,iCOff,ijCOff)
#endif
#ifdef RFLU
  DO ijkC = 1,region%grid%nCells
#endif
        esg4 = cEpsilon*trace(ijkC)
        diss(CV_MIXT_ENER,ijkC) = diss(CV_MIXT_ENER,ijkC) + esg4*beta 
#ifdef RFLO
      ENDDO ! i
    ENDDO   ! j
  ENDDO     ! k
#endif
#ifdef RFLU
  ENDDO     ! ijkC
#endif

#ifdef STATS
! if desired, collect quantities of interest at cell centers for statistics

  IF (region%turbInput%nSt > 0) THEN
    ibn = LBOUND( trace,1 )
    ien = UBOUND( trace,1 )
    ALLOCATE( colVar(1,ibn:ien) )
    colVar(1,:) = cEpsilon*trace(:)
    CALL TURB_StatCCollector( region,3,3,colVar )
    DEALLOCATE( colVar)
  ENDIF
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_LesEsgModel4

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_LesEsgModel4.F90,v $
! Revision 1.7  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2004/10/22 23:17:46  wasistho
! collect esg4 into statistics
!
! Revision 1.4  2004/04/08 20:22:02  wasistho
! replaced REAL_SMALL by actual value
!
! Revision 1.3  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.2  2004/03/19 02:49:45  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.6  2004/02/26 21:23:16  wasistho
! changed turb%esg.. to global%esg..
!
! Revision 1.5  2004/02/24 21:03:39  wasistho
! used local (region) production-dissipation balance for non-uniform MPI load
!
! Revision 1.4  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.3  2003/04/04 22:38:00  wasistho
! remove warning to inactivate energy model
!
! Revision 1.2  2002/12/12 03:28:09  wasistho
! Change mpi_comm_world to global%mpiComm
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************







