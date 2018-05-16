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
! Purpose: update pressure gradient based on mass flux balance
!
! Description: when computed mass flux is lower than reference mass flux
!              pressure gradient is adjusted in such away that the bulk flow 
!              accelerate, in other words mass flux is conserved
!
! Input: region = data of current region.
!
! Output: region%levels%mixt%rhs = updated pressure gradient added to residual.
!
! Notes: If MPI, nProcAlloc should be equal to nRegions, otherwise the wall
!        normal extent should be covered in one region.
!
!******************************************************************************
!
! $Id: PERI_coPgradUpdateFlu.F90,v 1.9 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_CoPgradUpdate( region )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModBndPatch, ONLY   : t_patch
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PERI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: ijkC, ifc, iPatch

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch

  INTEGER :: n, ijkC0, ijkCB, istage

  REAL(RFREAL)  :: rhoc, rhouc, rhovc, rhowc, rhoec, pc
  REAL(RFREAL)  :: refMassFlux, rescaler, alpha, massFlux
  REAL(RFREAL)  :: cfl, dt, rucm, delta, deltaMassFlux, averSurf
  REAL(RFREAL)  :: sndMassFlux, rcvMassFlux, sndAverSurf,  rcvAverSurf
  REAL(RFREAL)  :: sndDt, rcvDt, avgDt, dtFactor, refMeanPgrad, refTauWall
  REAL(RFREAL)  :: tauWall, avgTauWall, sndTauWall, rcvTauWall
  REAL(RFREAL)  :: rn, sndRn, rcvRn, dy, fnMag, volTot
  REAL(RFREAL), POINTER :: xyz(:,:), vol(:), cofg(:,:)
  REAL(RFREAL), POINTER :: cv(:,:), tv(:,:), rdt(:)

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_coPgradUpdateFlu.F90,v $ $Revision: 1.9 $'

  global => region%global
  CALL RegisterFunction( global,'PERI_CoPgradUpdate',&
  'PERI_coPgradUpdateFlu.F90' )

! get dimensions, parameters and pointers -------------------------------------
 
  xyz  => region%grid%xyz
  vol  => region%grid%vol
  cofg => region%grid%cofg
  cv   => region%mixt%cv
  tv   => region%mixt%tv
  rdt  => region%dt

  istage      = MOD(region%irkStep,global%nrkSteps) + 1
  alpha       = region%mixtInput%ark(istage)
  refMassFlux = region%periInput%bulkmFlux
  delta       = global%refLength

  IF (region%periInput%flowKind == PERI_FLOW_CPR) THEN
    rescaler = region%periInput%cprEpsilon
    IF (global%flowType == FLOW_UNSTEADY .AND. &
        global%solverType == SOLV_EXPLICIT ) alpha = 1._RFREAL
  ELSEIF (region%periInput%flowKind == PERI_FLOW_CHANNEL) THEN
    rescaler = 1._RFREAL
  ENDIF

  IF (global%flowType==FLOW_STEADY) THEN
    cfl   = region%mixtInput%cfl
    avgDt = 0._RFREAL
    IF (global%solverType == SOLV_IMPLICIT ) THEN ! Dual Tst
      dtFactor = 10._RFREAL
    ELSE
      dtFactor = 0.9_RFREAL
    ENDIF
    n = 0
    DO ijkC0=1,region%grid%nCellsTot
      n = n+1
      avgDt = avgDt + rdt(ijkC0)
    ENDDO
    rn    = REAL(n)
#ifdef MPI
    sndDt = avgDt
    CALL MPI_ALLREDUCE( sndDt,rcvDt,1,MPI_DOUBLE_PRECISION, MPI_SUM, &
                        global%mpiComm, global%mpierr )
    avgDt = rcvDt
    sndRn = rn
    CALL MPI_ALLREDUCE( sndRn,rcvRn,1,MPI_DOUBLE_PRECISION, MPI_SUM, &
                        global%mpiComm, global%mpierr )
    rn = rcvRn
#endif
    avgDt = dtFactor*avgDt/rn
    dt    = cfl*avgDt
!    IF (global%solverType == SOLV_IMPLICIT ) dt = global%dtMin ! Dual Tst
  ELSE
    dt   = global%dtMin
  ENDIF

! check if pressure gradient constant in all regions

  refMeanPgrad = region%periInput%meanPgrad

#ifdef MPI
  IF (refMeanPgrad /= global%moduleVar(1) ) THEN
    CALL ErrorStop( global,ERR_PERI_PHYSPARAM,__LINE__, &
         'mean pressure gradient vary between regions' )
  ENDIF
#endif 

! compute mass flux prior rho.u update -------------------------------------

  rucm   = 0._RFREAL
  volTot = 0._RFREAL
  DO ijkC = 1, region%grid%nCellsTot
    rucm   = rucm + cv(CV_MIXT_XMOM,ijkC0)*vol(ijkC)
    volTot = volTot + vol(ijkC)
  END DO
  massFlux = rucm
  averSurf = 0.5_RFREAL*volTot/delta

! this part only for channel flow

  avgTauWall = 0._RFREAL

  IF (region%periInput%flowKind == PERI_FLOW_CHANNEL .AND. &
      region%periInput%split(JCOORD) /= OFF ) THEN

    DO iPatch=1,region%grid%nPatches
      patch => region%patches(iPatch)
      IF (patch%bcType>=BC_NOSLIPWALL .AND. &
          patch%bcType<=BC_NOSLIPWALL+BC_RANGE) THEN   ! my boundary type

        DO ifc=1,patch%nBFaces
          ijkCB = patch%bf2c(ifc)

          dy    = ABS( cofg(YCOORD,ijkCB)-patch%fc(YCOORD,ifc) )
          fnMag = patch%fn(XYZMAG,ifc)

          tauWall = 0.5_RFREAL*tv(TV_MIXT_MUEL,ijkCB)* &
                    cv(CV_MIXT_XMOM,ijkCB)/cv(CV_MIXT_DENS,ijkCB)/dy
          avgTauWall = avgTauWall + tauWall*fnMag
        ENDDO ! ifc
      ENDIF   ! NOSLIPWALL
    ENDDO     ! iPatch  
  ENDIF       ! flowKind=channel

! get average mass flux over all processors

#ifdef MPI
  sndMassFlux = massFlux
  CALL MPI_ALLREDUCE( sndMassFlux,rcvMassFlux,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                      global%mpiComm, global%mpierr )
  massFlux = rcvMassFlux

  IF (region%periInput%flowKind == PERI_FLOW_CHANNEL) THEN
    sndTauWall = avgTauWall
    CALL MPI_ALLREDUCE( sndTauWall,rcvTauWall,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                        global%mpiComm, global%mpierr )
    avgTauWall = rcvTauWall
  ENDIF

  IF (region%periInput%split(JCOORD) == OFF) THEN
    sndAverSurf = averSurf
    CALL MPI_ALLREDUCE( sndAverSurf,rcvAverSurf,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                        global%mpiComm, global%mpierr )
    averSurf = rcvAverSurf
  ENDIF
#endif
  IF (global%myProcId==MASTERPROC) write(*,*)massFlux,averSurf

! channel half width needed to update pressure gradient

  massFlux   = massFlux/averSurf
  avgTauWall = 0.5_RFREAL*avgTauWall/averSurf

! update pressure gradient

  deltaMassFlux = refMassFlux - massFlux
  IF (global%myProcId==MASTERPROC) WRITE(*,*) region%iRegionGlobal, &
    'refMassFlux-massFlux prior updating u',refMassFlux,massFlux,deltaMassFlux, &
                                            avgTauWall/delta

  refMeanPgrad  = refMeanPgrad - 0.5_RFREAL*deltaMassFlux/ &
                  (delta*rescaler*alpha*dt)

  IF ((region%periInput%flowKind == PERI_FLOW_CHANNEL) .AND. &
      (region%periInput%pgradType == CNL_PGRAD_TAUWALL)) THEN
!    refTauWall   = global%refDensity*region%periInput%cnlUtau**2
!    refMeanPgrad = -refTauWall/global%refLength
    refMeanPgrad = MIN( refMeanPgrad, -avgTauWall/global%refLength )
  ENDIF

  IF (global%myProcId==MASTERPROC) WRITE(*,*) region%iRegionGlobal, &
    ' meanPgrad ',refMeanPgrad,dt,alpha

! store pressure gradient in data structure
  region%periInput%meanPgrad = refMeanPgrad

! update x-momentum due to updated pressure gradient/ mass flux ---------------

  IF (region%periInput%flowKind == PERI_FLOW_CPR) THEN

    rucm   = 0._RFREAL
    DO ijkC = 1,region%grid%nCellsTot
      rucm   = rucm + cv(CV_MIXT_XMOM,ijkC0)*vol(ijkC)
    END DO
    massFlux = rucm

#ifdef MPI
    CALL MPI_BARRIER( global%mpiComm, global%mpierr )
    sndMassFlux = massFlux
    CALL MPI_ALLREDUCE( sndMassFlux,rcvMassFlux,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                        global%mpiComm, global%mpierr )
    massFlux = rcvMassFlux
#endif

    massFlux = massFlux/averSurf

    deltaMassFlux = refMassFlux - massFlux
    IF (global%myProcId==MASTERPROC) WRITE(*,*) region%iRegionGlobal, &
      'refMassFlux-massFlux AFTER updating u',refMassFlux,massFlux,deltaMassFlux

  ENDIF ! flowKind=CPR

! save pressure gradient in a global variable for restart purpose

  global%moduleVar(1) = refMeanPgrad

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PERI_CoPgradUpdate

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_coPgradUpdateFlu.F90,v $
! Revision 1.9  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2004/12/11 03:49:52  wasistho
! resereve dt=global%dtMin for dual tst
!
! Revision 1.6  2004/12/08 19:42:03  wasistho
! used dtFactor
!
! Revision 1.5  2004/12/08 01:47:50  wasistho
! bug fix, MPI_MAX to MPI_SUM
!
! Revision 1.4  2004/12/04 06:51:32  wasistho
! take global%dtMin for dt regardless steady or unsteady
!
! Revision 1.3  2004/06/17 20:04:02  wasistho
! compiled with RFLU
!
! Revision 1.2  2004/06/17 00:50:53  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/06/09 01:10:26  wasistho
! changed nomenclature
!
! Revision 1.15  2004/01/23 22:42:49  wasistho
! set avgDt to 0.9 of the true value for stability of steady state computations
!
! Revision 1.14  2004/01/23 03:11:22  wasistho
! modified dt, added computed tauwall, unst.impl.condition, no cnl massflx correction
!
! Revision 1.13  2003/10/21 20:26:37  wasistho
! global max for dt of unsteady flow
!
! Revision 1.12  2003/10/20 00:32:07  wasistho
! modified istage definition
!
! Revision 1.11  2003/10/18 00:47:50  wasistho
! modified dt for unsteady flow
!
! Revision 1.10  2003/09/18 01:56:26  wasistho
! added ijksplit and pgradType in PERI_PgradUpdate
!
! Revision 1.9  2003/09/12 20:57:45  wasistho
! update pressure gradient only for CPR, fixed for channel
!
! Revision 1.8  2003/09/07 23:32:02  wasistho
! shift ISTAGE 1 stage forward to be consistent
!
! Revision 1.7  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.6  2003/04/25 23:17:49  wasistho
! update cpr pgrad per time-step for unsteady
!
! Revision 1.5  2003/04/05 04:06:39  wasistho
! shift forward time stage
!
! Revision 1.4  2003/04/04 22:27:19  wasistho
! correct averSize for multi processors
!
! Revision 1.3  2003/04/03 20:57:02  wasistho
! replace regions to region in pgradUpdate
!
! Revision 1.2  2003/04/03 00:30:07  wasistho
! uniform pgrad multiblock/proc
!
! Revision 1.1.1.1  2003/03/29 03:36:30  wasistho
! install ROCPERI
!
!
!******************************************************************************







