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
! ******************************************************************************
!
! Purpose: calculate solution at a new time level/iteration.
!
! Description: the governing equations are integrated in time using
!              an explicit, multistage (Runge-Kutta type) temporal scheme.
!              Spatial and temporal discretizations are independent of
!              each other based on the method of lines.
!
! Input: regions    = data of all regions
!        ftermNew   = compute new forcing term (true/false)
!        residFterm = add forcing term to residual (true/false).
!
! Output: regions%levels%mixt = new solution (and forcing term) after
!                               one time step/iteration.
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: ExplicitMultistage.F90,v 1.15 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ExplicitMultistage( regions,ftermNew,residFterm )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI
  USE ModParameters

  USE RFLU_ModConvertCv
  USE RFLU_ModDifferentiationCells
  USE RFLU_ModLimiters, ONLY: RFLU_LimitGradCellsSimple
  USE RFLU_ModMPI
  USE ModInterfaces, ONLY: ConvectiveFluxes, &
                           MixtureProperties, &
                           NumericalDissipation, &
                           RFLU_CheckPositivityWrapper, &
                           RFLU_CheckValidityWrapper, &
                           RFLU_SetVars, & 
                           RFLU_TimeStepInviscid, &
                           RFLU_TimeStepViscous, &
                           RFLU_ZeroVirtualCellVars, &
                           SourceTerms, &
                           UpdateTbc, &
                           ViscousFluxes
#ifdef PEUL
  USE ModInterfacesEulerian, ONLY : PEUL_SpectralRadii, &
                                    PEUL_ResidualSmoothingCoeffs, &
                                    PEUL_ResidualSmoothing
#endif
#ifdef INRT
  USE ModInterfacesInteract, ONLY : INRT_SourceTerms
#endif
#ifdef RADI
  USE ModInterfacesRadiation, ONLY : RADI_SourceTerms
#endif
#ifdef PERI
  USE ModInterfacesPeriodic, ONLY : PERI_SourceTerms, PERI_SolutionUpdate
#endif
#ifdef TURB
  USE ModInterfacesTurbulence, ONLY : TURB_EmsInit, TURB_RansConvectiveFluxes, &
                TURB_RansNumericalDissipation,      TURB_RansSourceTerms, &
                TURB_RansZeroDummyCells,            TURB_SolutionUpdate
#endif

  IMPLICIT NONE

#ifdef GENX
  INCLUDE 'comf90.h'
#endif

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  LOGICAL :: ftermNew, residFterm

! ... loop variables
  INTEGER :: iReg, iRegLocal, ic, istage

! ... local variables
  INTEGER :: ibc, iec, iecTot, ldiss(5), flowModel, gasModel

  LOGICAL :: moveGrid

  REAL(RFREAL)          :: cfl, ark(5), betrk(5), blend1, fac, adtv
  REAL(RFREAL)          :: alpha, time
  REAL(RFREAL), POINTER :: cv(:,:), cvOld(:,:), dt(:), diss(:,:), rhs(:,:)
  REAL(RFREAL), POINTER :: vol(:), fterm(:,:)

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'ExplicitMultistage',&
  'ExplicitMultistage.F90' )

! loop over stages and regions ------------------------------------------------

  DO istage=1,regions(1)%global%nrkSteps
    DO iRegLocal = 1,global%nRegionsLocal
      iReg = iRegLocal
        regions(iReg)%irkStep = istage

! ----- get dimensions and pointers

        ldiss(:)  = regions(iReg)%mixtInput%ldiss(:)
        cfl       = regions(iReg)%mixtInput%cfl
        ark(:)    = regions(iReg)%mixtInput%ark(:)
        betrk(:)  = regions(iReg)%mixtInput%betrk(:)
        flowModel = regions(iReg)%mixtInput%flowModel
        gasModel  = regions(iReg)%mixtInput%gasModel
        moveGrid  = regions(iReg)%mixtInput%moveGrid

        pRegion => regions(iReg)

        ibc    = 1
        iec    = pRegion%grid%nCells
        iecTot = pRegion%grid%nCellsTot

        cv    => pRegion%mixt%cv
        cvOld => pRegion%mixt%cvOld
        diss  => pRegion%mixt%diss
        rhs   => pRegion%mixt%rhs
        vol   => pRegion%grid%vol
        dt    => pRegion%dt
        IF (residFterm) fterm => pRegion%mixt%fterm

! ----- store previous solution; set dissipation to zero (first stage)

        IF (istage == 1) THEN
          DO ic=ibc,iecTot
            cvOld(CV_MIXT_DENS,ic) = cv(CV_MIXT_DENS,ic)
            cvOld(CV_MIXT_XMOM,ic) = cv(CV_MIXT_XMOM,ic)
            cvOld(CV_MIXT_YMOM,ic) = cv(CV_MIXT_YMOM,ic)
            cvOld(CV_MIXT_ZMOM,ic) = cv(CV_MIXT_ZMOM,ic)
            cvOld(CV_MIXT_ENER,ic) = cv(CV_MIXT_ENER,ic)
            diss(CV_MIXT_DENS,ic)  = 0._RFREAL
            diss(CV_MIXT_XMOM,ic)  = 0._RFREAL
            diss(CV_MIXT_YMOM,ic)  = 0._RFREAL
            diss(CV_MIXT_ZMOM,ic)  = 0._RFREAL
            diss(CV_MIXT_ENER,ic)  = 0._RFREAL
          ENDDO
        ENDIF

! ----- initialize dissipation (later stages)

        IF (istage>1 .AND. ldiss(istage)/=0) THEN
          blend1 = 1._RFREAL - betrk(istage)
          DO ic=ibc,iecTot
            diss(CV_MIXT_DENS,ic) = blend1*diss(CV_MIXT_DENS,ic)
            diss(CV_MIXT_XMOM,ic) = blend1*diss(CV_MIXT_XMOM,ic)
            diss(CV_MIXT_YMOM,ic) = blend1*diss(CV_MIXT_YMOM,ic)
            diss(CV_MIXT_ZMOM,ic) = blend1*diss(CV_MIXT_ZMOM,ic)
            diss(CV_MIXT_ENER,ic) = blend1*diss(CV_MIXT_ENER,ic)
          ENDDO
        ENDIF

! ----- Compute cell gradients for higher-order scheme

        IF ( regions(iReg)%mixtInput%spaceOrder > 1 ) THEN
          pRegion => regions(iReg)
          CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)
          CALL RFLU_ComputeGradCells(pRegion,CV_MIXT_DENS,CV_MIXT_PRES, &
                                     GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                     pRegion%mixt%cv,pRegion%mixt%gradCell)
          CALL RFLU_ComputeGradCellsENO(pRegion,GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                        pRegion%mixt%gradCell)
          CALL RFLU_LimitGradCellsSimple(pRegion,CV_MIXT_DENS,CV_MIXT_PRES, &
                                         GRC_MIXT_DENS,GRC_MIXT_PRES, &
                                         pRegion%mixt%cv,pRegion%mixt%cvInfo, &
                                         pRegion%mixt%gradCell)
          CALL RFLU_ConvertCvPrim2Cons(pRegion,CV_MIXT_STATE_CONS)
        END IF ! regions

#ifdef TURB
          IF (flowModel == FLOW_NAVST .AND. &
              regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) &
            CALL TURB_EmsInit( regions(iReg),istage )
#endif

! ----- compute numerical dissipation

        IF (ldiss(istage) /= 0) THEN
          CALL NumericalDissipation( regions(iReg) )
        ENDIF

#ifdef TURB
        IF (ldiss(istage) /= 0 .AND. flowModel == FLOW_NAVST .AND. &
            regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) THEN
          CALL TURB_RansNumericalDissipation( regions(iReg) )
        ENDIF
#endif

! ----- compute viscous fluxes

        IF (flowModel==FLOW_NAVST .AND. ldiss(istage)/=0) THEN
          CALL ViscousFluxes( regions(iReg) )
        ENDIF

! ----- compute convective fluxes; form residual

        CALL ConvectiveFluxes( regions(iReg) )

#ifdef TURB
        IF (flowModel == FLOW_NAVST .AND. &
            regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) THEN
          CALL TURB_RansConvectiveFluxes( regions(iReg) )
        ENDIF
#endif

! ----- add source terms

        CALL SourceTerms( regions(iReg) )

#ifdef INRT
        IF (global%inrtUsed) THEN
          CALL INRT_SourceTerms( regions(iReg) )
        ENDIF
#endif
#ifdef RADI
        CALL RADI_SourceTerms( regions(iReg) )
#endif
#ifdef PERI
        CALL PERI_SourceTerms( regions(iReg) )
#endif
#ifdef TURB
       IF (flowModel == FLOW_NAVST .AND. &
           regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) THEN
         CALL TURB_RansSourceTerms( regions(iReg) )
       ENDIF
#endif

! ----- zero out residuals in dummy cells

        pRegion => regions(iReg)
        CALL RFLU_ZeroVirtualCellVars(pRegion,rhs)
#ifdef TURB
        IF (flowModel == FLOW_NAVST .AND. &
            regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) THEN
          CALL TURB_RansZeroDummyCells( regions(iReg) )
        ENDIF
#endif

! ----- compute new forcing term

        IF (ftermNew) THEN
          DO ic=ibc,iecTot
            fterm(CV_MIXT_DENS,ic) = fterm(CV_MIXT_DENS,ic) - &
                                       rhs(CV_MIXT_DENS,ic)
            fterm(CV_MIXT_XMOM,ic) = fterm(CV_MIXT_XMOM,ic) - &
                                       rhs(CV_MIXT_XMOM,ic)
            fterm(CV_MIXT_YMOM,ic) = fterm(CV_MIXT_YMOM,ic) - &
                                       rhs(CV_MIXT_YMOM,ic)
            fterm(CV_MIXT_ZMOM,ic) = fterm(CV_MIXT_ZMOM,ic) - &
                                       rhs(CV_MIXT_ZMOM,ic)
            fterm(CV_MIXT_ENER,ic) = fterm(CV_MIXT_ENER,ic) - &
                                       rhs(CV_MIXT_ENER,ic)
          ENDDO
        ENDIF

        ftermNew = .false.

! ----- residual (+ forcing term) * time step / volume

        fac = ark(istage)*cfl

        IF (residFterm) THEN
          DO ic=ibc,iec
            adtv = fac*dt(ic)/vol(ic)
            rhs(CV_MIXT_DENS,ic) = adtv*(rhs(CV_MIXT_DENS,ic)+ &
                                         fterm(CV_MIXT_DENS,ic))
            rhs(CV_MIXT_XMOM,ic) = adtv*(rhs(CV_MIXT_XMOM,ic)+ &
                                         fterm(CV_MIXT_XMOM,ic))
            rhs(CV_MIXT_YMOM,ic) = adtv*(rhs(CV_MIXT_YMOM,ic)+ &
                                         fterm(CV_MIXT_YMOM,ic))
            rhs(CV_MIXT_ZMOM,ic) = adtv*(rhs(CV_MIXT_ZMOM,ic)+ &
                                         fterm(CV_MIXT_ZMOM,ic))
            rhs(CV_MIXT_ENER,ic) = adtv*(rhs(CV_MIXT_ENER,ic)+ &
                                         fterm(CV_MIXT_ENER,ic))
          ENDDO
        ELSE
          DO ic=ibc,iec
            adtv = fac*dt(ic)/vol(ic)
            rhs(CV_MIXT_DENS,ic) = adtv*rhs(CV_MIXT_DENS,ic)
            rhs(CV_MIXT_XMOM,ic) = adtv*rhs(CV_MIXT_XMOM,ic)
            rhs(CV_MIXT_YMOM,ic) = adtv*rhs(CV_MIXT_YMOM,ic)
            rhs(CV_MIXT_ZMOM,ic) = adtv*rhs(CV_MIXT_ZMOM,ic)
            rhs(CV_MIXT_ENER,ic) = adtv*rhs(CV_MIXT_ENER,ic)
          ENDDO
        ENDIF

! ----- update solution

        IF (global%solverType == SOLV_IMPLICIT) THEN
          fac = 1.5_RFREAL*ark(istage)*cfl/global%dtMin
          DO ic=ibc,iec
            adtv = 1._RFREAL/(1._RFREAL+fac*dt(ic))
            cv(CV_MIXT_DENS,ic) = cvOld(CV_MIXT_DENS,ic) - &
                                  adtv*rhs(CV_MIXT_DENS,ic)
            cv(CV_MIXT_XMOM,ic) = cvOld(CV_MIXT_XMOM,ic) - &
                                  adtv*rhs(CV_MIXT_XMOM,ic)
            cv(CV_MIXT_YMOM,ic) = cvOld(CV_MIXT_YMOM,ic) - &
                                  adtv*rhs(CV_MIXT_YMOM,ic)
            cv(CV_MIXT_ZMOM,ic) = cvOld(CV_MIXT_ZMOM,ic) - &
                                  adtv*rhs(CV_MIXT_ZMOM,ic)
            cv(CV_MIXT_ENER,ic) = cvOld(CV_MIXT_ENER,ic) - &
                                  adtv*rhs(CV_MIXT_ENER,ic)
          ENDDO
        ELSE
          DO ic=ibc,iec
            cv(CV_MIXT_DENS,ic) = cvOld(CV_MIXT_DENS,ic) - rhs(CV_MIXT_DENS,ic)
            cv(CV_MIXT_XMOM,ic) = cvOld(CV_MIXT_XMOM,ic) - rhs(CV_MIXT_XMOM,ic)
            cv(CV_MIXT_YMOM,ic) = cvOld(CV_MIXT_YMOM,ic) - rhs(CV_MIXT_YMOM,ic)
            cv(CV_MIXT_ZMOM,ic) = cvOld(CV_MIXT_ZMOM,ic) - rhs(CV_MIXT_ZMOM,ic)
            cv(CV_MIXT_ENER,ic) = cvOld(CV_MIXT_ENER,ic) - rhs(CV_MIXT_ENER,ic)
          ENDDO
        ENDIF

#ifdef PERI
        IF (regions(iReg)%periInput%flowKind /= OFF) THEN
           CALL PERI_SolutionUpdate( regions(iReg) )
        ENDIF
#endif
#ifdef TURB
        IF (flowModel == FLOW_NAVST .AND. &
            regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE) THEN
          CALL TURB_SolutionUpdate( regions(iReg),istage,ibc,iec )
        ENDIF
#endif

! ----- check for positivity/validity -----------------------------------------

        pRegion => regions(iReg)
        CALL RFLU_CheckValidityWrapper(pRegion)
        CALL RFLU_CheckPositivityWrapper(pRegion)

! ----- exchange conservative variables with other processors

        pRegion => regions(iReg)
        CALL RFLU_MPI_ISendWrapper(pRegion)
        CALL RFLU_SetVars(pRegion,1,pRegion%grid%nCells)

! ----- update dependent variables

    ENDDO    ! iReg

    CALL RFLU_MPI_CopyWrapper(regions)

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)

      CALL RFLU_MPI_RecvWrapper(pRegion)    
      CALL RFLU_SetVars(pRegion,pRegion%grid%nCells+1,pRegion%grid%nCellsTot) 
    END DO ! iReg

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)

      CALL RFLU_MPI_ClearRequestWrapper(pRegion)
    END DO ! iReg

  ENDDO   ! istage

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE ExplicitMultistage

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ExplicitMultistage.F90,v $
! Revision 1.15  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.14  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.13  2006/03/26 20:21:14  haselbac
! Changed to wrappers bcos of GL model
!
! Revision 1.12  2005/12/03 19:45:05  haselbac
! Apparent bug fix: Separated call to RFLU_MPI_ClearRequestWrapper into separate loop
!
! Revision 1.11  2005/10/31 21:09:34  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.10  2005/10/05 13:47:23  haselbac
! Adapted to new module for cell grads
!
! Revision 1.9  2005/07/13 17:26:20  haselbac
! Bug fix: Adapted to changes to RFLU_LimitGradCells
!
! Revision 1.8  2005/05/16 20:38:42  haselbac
! Renamed RFLU_ZeroDummyCells, moved computation of dt into RFLU_TimeStepping
!
! Revision 1.7  2005/04/29 00:06:09  haselbac
! Added routines to clear send requests
!
! Revision 1.6  2005/04/15 15:06:01  haselbac
! Converted to MPI
!
! Revision 1.5  2005/04/06 01:59:42  wasistho
! mv call to PERI_CoMeanCorrection to after BC treatment, commented for now
!
! Revision 1.4  2005/03/10 02:08:25  wasistho
! commented PERI_coMeanCorrection temporarily for testing
!
! Revision 1.3  2005/03/07 05:05:18  wasistho
! install hybrid DESSA turbulence model
!
! Revision 1.2  2004/12/28 22:49:15  wasistho
! moved RFLO_Bcond* and RFLO_BoundaryCond* routines into RFLO_ModBoundaryConditions
!
! Revision 1.1  2004/12/01 16:48:30  haselbac
! Initial revision after changing case
!
! Revision 1.68  2004/11/17 23:44:33  wasistho
! used generic RK-update for rocturb
!
! Revision 1.67  2004/11/14 19:34:26  haselbac
! Replaced call to mixtureProperties by RFLU_SetVars
!
! Revision 1.66  2004/10/19 19:25:34  haselbac
! Adapted calls to time step routines
!
! Revision 1.65  2004/07/26 19:08:59  wasistho
! add RFLO_CheckValidity
!
! Revision 1.64  2004/03/27 03:02:15  wasistho
! added ifdef RFLO within ifdef TURB
!
! Revision 1.63  2004/03/20 00:26:35  wasistho
! set turb_rflo_ransNumericalDiss to turb_ransNumerical..
!
! Revision 1.62  2004/03/19 02:39:55  wasistho
! renamed TURB_RFLO_RansZeroDummyCells to TURB_RansZeroDummyCells
!
! Revision 1.61  2004/03/11 03:32:13  wasistho
! changed rocturb nomenclature
!
! Revision 1.60  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.59  2004/03/02 21:49:21  jferry
! Added inrtUsed flag to mixture data structure
!
! Revision 1.58  2004/02/26 21:13:55  wasistho
! changed TURB_ransEmsInit to TURB_emsInit
!
! Revision 1.57  2004/01/29 22:52:41  haselbac
! Added info argument to cell-gradient routine call
!
! Revision 1.56  2003/12/04 03:22:59  haselbac
! Added second-order scheme, viscous fluxes, validity check
!
! Revision 1.55  2003/11/25 21:01:39  haselbac
! Added routine to update dummy cells
!
! Revision 1.54  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.51  2003/10/27 04:49:29  wasistho
! replace ransCentralDiss by ransNumericalDiss.
!
! Revision 1.50  2003/10/21 03:57:42  wasistho
! put turb and peul spectralradii before its smoothingcoef
!
! Revision 1.49  2003/10/16 20:12:27  wasistho
! completed incorporation of RaNS/DES
!
! Revision 1.48  2003/10/03 20:12:11  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.47  2003/09/26 21:45:34  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.46  2003/08/28 20:32:48  wasistho
! excluced ModInterfacesTurbulence,Radiation,Periodic from ModInterfaces
!
! Revision 1.45  2003/08/13 02:19:42  wasistho
! added call to radiation source term routine
!
! Revision 1.44  2003/07/22 01:52:40  haselbac
! Bug fix: dual-time stepping only for RFLO, lead to core dump
!
! Revision 1.43  2003/07/03 21:48:44  jblazek
! Implemented dual-time stepping.
!
! Revision 1.42  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.41  2003/05/09 17:01:03  jiao
! Renamed the COM_call_function_handlers to COM_call_function.
!
! Revision 1.40  2003/04/05 02:01:08  wasistho
! regions to region in PERI_solutionUpdate
!
! Revision 1.39  2003/03/29 03:26:21  wasistho
! install ROCPERI
!
! Revision 1.38  2003/03/15 16:17:58  haselbac
! Added zeroing of dummies (?) and changed FEM to IDXL
!
! Revision 1.37  2003/03/05 20:41:21  jiao
! ACH: Split update inbuff calls into get correct dependency of rhof on mdot
!
! Revision 1.36  2003/03/04 22:12:34  jferry
! Initial import of Rocinteract
!
! Revision 1.35  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
! Revision 1.34  2002/12/06 22:29:25  jblazek
! Corrected bug for geometry exchange between minimal patches.
!
! Revision 1.33  2002/10/27 18:47:32  haselbac
! Added call to RFLU_CheckPositivity
!
! Revision 1.32  2002/10/17 06:50:42  jiao
! Changed 0. to 0._RFREAL.
!
! Revision 1.31  2002/10/12 20:21:25  jblazek
! Rearranged UpdateTbc (after getting values from GenX);
! corrected bug in externalBc.
!
! Revision 1.30  2002/10/05 18:34:08  haselbac
! Moved ifdef CHARM inside ifdef RFLU
!
! Revision 1.29  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.27  2002/09/09 13:57:49  haselbac
! added viscous routines, bug fix for iec, mixtInput under regions
!
! Revision 1.26  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.25  2002/08/29 21:53:18  jblazek
! Added support for moving grids.
!
! Revision 1.24  2002/08/16 21:33:47  jblazek
! Changed interface to MixtureProperties.
!
! Revision 1.23  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
! Revision 1.22  2002/07/25 14:50:52  haselbac
! Added RFLU_ModFEM and call for update of dummy cells
!
! Revision 1.21  2002/07/23 20:28:19  wasistho
! #ifdef RFLO around ViscousFluxes removed
!
! Revision 1.20  2002/07/22 17:01:23  jblazek
! Removed MPI_Barrier at the end of stage loop.
!
! Revision 1.19  2002/07/05 23:20:46  jblazek
! Corrected bug in perfgasDependentVars.F90; did some cosmetics.
!
! Revision 1.18  2002/06/27 15:46:27  haselbac
! Added FEM calls for communication - not activated yet
!
! Revision 1.17  2002/06/14 20:11:48  haselbac
! Deleted ModLocal, renamed local%nRegions to global%nRegionsLocal
!
! Revision 1.16  2002/06/07 16:40:36  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.15  2002/05/28 14:14:39  haselbac
! Enclosed viscous fluxes within RFLO conditional statement
!
! Revision 1.14  2002/05/28 13:44:05  haselbac
! Cosmetic changes only?
!
! Revision 1.13  2002/05/21 01:51:56  wasistho
! add viscous terms
!
! Revision 1.12  2002/05/04 16:31:49  haselbac
! Added RFLU statements
!
! Revision 1.11  2002/04/12 17:36:23  jblazek
! Added timer.
!
! Revision 1.10  2002/04/01 19:36:07  jblazek
! Added routine to clear send requests.
!
! Revision 1.9  2002/03/30 00:50:49  jblazek
! Cleaned up with flint.
!
! Revision 1.8  2002/03/29 23:15:22  jblazek
! Corrected bug in MPI send.
!
! Revision 1.7  2002/03/18 22:25:45  jblazek
! Finished multiblock and MPI.
!
! Revision 1.6  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.5  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.4  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.3  2002/01/23 23:37:32  jblazek
! All blocks passed to time integration routines.
!
! Revision 1.2  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.1  2002/01/16 22:03:34  jblazek
! Added time-stepping routines.
!
!******************************************************************************







