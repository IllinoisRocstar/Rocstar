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
! Purpose: calculate solution at a new time level.
!
! Description: the governing equations are integrated in time using
!              the classical 4-stage Runge-Kutta method (4th-order in
!              time) in low-storage formulation.
!
! Input: regions = data of all regions.
!
! Output: regions%levels%mixt = new solution after one time step.
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: RungeKuttaMP.F90,v 1.20 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RungeKuttaMP( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters

#ifdef RFLU
  USE RFLU_ModBoundXvUtils
  USE RFLU_ModGridSpeedUtils, ONLY: RFLU_DescaleGridSpeeds, &
                                    RFLU_ScaleGridSpeeds, &
                                    RFLU_SetGridSpeedScaleFactor
  USE RFLU_ModTimeZoom, ONLY: RFLU_TimeZoomDriver
  USE RFLU_ModMPI
  USE RFLU_ModNSCBC
  USE RFLU_ModRelatedPatches, ONLY: RFLU_RELP_TransformWrapper
  USE RFLU_ModTime, ONLY: RFLU_SetTimeRK
#endif

  USE ModInterfaces, ONLY : AfterUpdateMP, &
                            CellGradientsMP, ConvectiveFluxesMP, &
                            GlobalCommunicationMP, InitCommunicationMP, &
                            NumericalDissipationMP, RKInitMP, RKUpdateMP, &
                            SourceTermsMP, UpdateBoundaryConditionsMP, &
                            UpdateDependentVarsMP, ViscousFluxesMP, &
                            ZeroDummyCellsMP, ZeroResidualsMP
#ifdef RFLU
  USE ModInterfaces, ONLY : RFLU_EquilibriumEulerian, &
                            RFLU_SetVarsContWrapper, &
                            RFLU_SetVarsDiscWrapper, & 
                            RFLU_UpdateBoundaryValues
#ifdef GENX
  USE RFLU_ModRocstarTools, ONLY: RFLU_GENX_InitBFLAG
#endif

#ifdef PLAG
  USE PLAG_RFLU_ModComm
#endif
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iRegLocal, istage

! ... local variables
  INTEGER :: flowModel

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RungeKuttaMP',&
  'RungeKuttaMP.F90' )

! loop over stages and regions ================================================
  

  DO istage=1,global%nrkSteps
#ifdef RFLO
    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor
#endif
#ifdef RFLU
    DO iRegLocal=1,global%nRegionsLocal
      iReg = iRegLocal
#endif

! ----- set pointer and get models --------------------------------------------

        pRegion => regions(iReg)

        flowModel  = regions(iReg)%mixtInput%flowModel
        regions(iReg)%irkStep = istage

#ifdef RFLU
! ----- Set RK time -----------------------------------------------------------

        CALL RFLU_SetTimeRK(pRegion,iStage)

! ----- RFLU fill GENX incoming buffers ---------------------------------------

#ifdef GENX
	CALL RFLU_GENX_InitBFLAG(pRegion)
#endif
        CALL RFLU_UpdateBoundaryValues(regions(iReg),istage)

! ----- Scale grid speeds -----------------------------------------------------

        CALL RFLU_SetGridSpeedScaleFactor(pRegion)
        CALL RFLU_ScaleGridSpeeds(pRegion)
#endif

! ----- store previous solution; set dissipation to zero ----------------------

        CALL RKInitMP( regions(iReg),istage )

! ----- compute cell-gradients for higher-order scheme ------------------------

        CALL CellGradientsMP( regions(iReg) )

! ----- compute numerical dissipation -----------------------------------------

        CALL NumericalDissipationMP( regions(iReg) )

! ----- compute viscous fluxes ------------------------------------------------

        IF ( flowModel == FLOW_NAVST ) THEN
          CALL ViscousFluxesMP( regions(iReg) )
        END IF ! flowModel

! ----- compute convective fluxes; form residual ------------------------------

        CALL ConvectiveFluxesMP( regions(iReg) )

! ----- zero residuals --------------------------------------------------------

        CALL ZeroResidualsMP(regions(iReg))

! ----- add source terms ------------------------------------------------------

        CALL SourceTermsMP( regions(iReg) )

#ifdef RFLU
! ----- add Equilibrium Eulerian corrections ----------------------------------

        CALL RFLU_EquilibriumEulerian( pRegion )
#endif

! ----- zero out residuals in dummy cells -------------------------------------

        CALL ZeroDummyCellsMP( regions(iReg) )
#ifdef RFLO
      ENDIF  ! region on this processor and active
#endif
#ifdef RFLU
! ----- Descale grid speeds -----------------------------------------------------
        CALL RFLU_DescaleGridSpeeds(pRegion)
#endif
    ENDDO    ! iReg


#ifdef RFLU
    IF(global%zoomFactor > 1) THEN
       CALL RFLU_TimeZoomDriver(regions)
    ENDIF ! global%zoomFactor
#endif
    

#ifdef RFLO
    DO iReg=1,global%nRegions
      IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
          regions(iReg)%active==ACTIVE) THEN              ! on my processor
#endif
#ifdef RFLU
    DO iRegLocal=1,global%nRegionsLocal
      iReg = iRegLocal
#endif

! ----- set Region pointer

        pRegion => regions(iReg)

! ----- update solution; sum up residuals -------------------------------------

        CALL RKUpdateMP( regions(iReg),iReg,istage )

#ifdef RFLU
        IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_BXV_ComputeVarsCv(pRegion)
        END IF ! RFLU_NSCBC_DecideHaveNSCBC
#endif

! ----- perform checks and enforce after-update conditions --------------------

        CALL AfterUpdateMP( pRegion,istage )



#ifdef RFLO
! ----- initiate communication kernel -----------------------------------------

        CALL InitCommunicationMP( regions,iReg,istage )
#endif

! ----- update dependent variables --------------------------------------------

#ifdef RFLU
        CALL RFLU_MPI_ISendWrapper(pRegion)
        CALL RFLU_SetVarsContWrapper(pRegion,1,pRegion%grid%nCells)

! ----- update dependent variables on boundary faces --------------------------

        IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_BXV_SetDependentVars(pRegion)
        END IF !

#endif
#ifdef RFLO
        CALL UpdateDependentVarsMP(regions(iReg))
#endif

#ifdef RFLO
      ENDIF  ! region on this processor and active
#endif
    ENDDO    ! iReg

#ifdef RFLO
! - facilitate global non-dummy communication ---------------------------------

    CALL GlobalCommunicationMP( regions )

! - update boundary conditions ------------------------------------------------

    CALL UpdateBoundaryConditionsMP( regions,istage )
#endif

#ifdef RFLU
    CALL RFLU_MPI_CopyWrapper(regions)

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)

      CALL RFLU_MPI_RecvWrapper(pRegion)
      CALL RFLU_SetVarsContWrapper(pRegion,pRegion%grid%nCells+1, & 
                                   pRegion%grid%nCellsTot)
      CALL RFLU_RELP_TransformWrapper(pRegion)                                   
    END DO ! iReg

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)

      CALL RFLU_MPI_ClearRequestWrapper(pRegion)
    END DO ! iReg

#ifdef PLAG
    IF ( global%plagUsed .EQV. .TRUE. ) THEN
      CALL PLAG_RFLU_CommDriver(regions)

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_SetVarsDiscWrapper(pRegion)
      END DO ! iReg
    END IF ! global%plagUsed 
#endif
#endif
  END DO ! istage

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE RungeKuttaMP

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RungeKuttaMP.F90,v $
! Revision 1.20  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.19  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.18  2007/04/20 16:07:48  mtcampbe
! Updating for burnout support function RFLU_GENX_InitBFLAG
!
! Revision 1.17  2007/04/14 14:23:24  mtcampbe
! Updated for TZ
!
! Revision 1.16  2007/03/27 00:39:50  haselbac
! Removed call to PLAG_CalcnPclsTotGlobal, now in RFLU_TimeStepping
!
! Revision 1.15  2007/03/20 22:02:29  fnajjar
! Included call to PLAG_CalcnPclsTotGlobal
!
! Revision 1.14  2006/08/21 16:10:01  haselbac
! Adapted to name change
!
! Revision 1.13  2006/08/19 15:48:25  mparmar
! Added computations of boundary Cv and Dv for NSCBC implementation
!
! Revision 1.12  2006/08/18 21:09:27  fnajjar
! Removed IF around PLAG_RFLU_CommDriver for serial periodic cases
!
! Revision 1.11  2006/03/25 21:40:03  haselbac
! Added call to transforming data on related patches, cosmetics
!
! Revision 1.10  2005/12/03 19:44:54  haselbac
! Apparent bug fix: Separated call to RFLU_MPI_ClearRequestWrapper into separate loop
!
! Revision 1.9  2005/12/01 21:52:14  fnajjar
! Added IF statement around PLAG_RFLU_CommDriver, only active for more than one nRegions
!
! Revision 1.8  2005/11/10 22:21:07  fnajjar
! ACH: Proper fix for updating PLAG dv
!
! Revision 1.7  2005/11/10 16:51:28  fnajjar
! Added plagUsed IF statement around PLAG routines
!
! Revision 1.6  2005/11/02 14:53:24  haselbac
! Fady: Temporary fix so comm particles get non-cv vars updated properly
!
! Revision 1.5  2005/05/18 22:04:41  fnajjar
! Added PLAG communication routines; only initial implementation
!
! Revision 1.4  2005/04/29 00:06:09  haselbac
! Added routines to clear send requests
!
! Revision 1.3  2005/04/15 15:06:06  haselbac
! Converted to MPI
!
! Revision 1.2  2005/03/31 16:31:02  haselbac
! Added call to RFLU_SetTimeRK
!
! Revision 1.1  2004/12/01 16:51:13  haselbac
! Initial revision after changing case
!
! Revision 1.18  2004/11/14 19:36:23  haselbac
! Replaced call to UpdateDependentVarsMP by RFLU_SetVarsWrapper
!
! Revision 1.17  2004/07/30 22:47:33  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.16  2004/04/14 02:07:02  haselbac
! Added grid-speed scaling calls for RFLU
!
! Revision 1.15  2004/03/25 21:14:20  jferry
! changed AfterUpdate to call most subroutines only after final RK stage
!
! Revision 1.14  2004/03/02 21:47:28  jferry
! Added After Update interactions
!
! Revision 1.13  2004/02/26 21:11:58  wasistho
! added globalCommunication
!
! Revision 1.12  2004/02/26 21:01:46  haselbac
! Enclosed updateBoundaryConditionsMP within ifdef RFLO
!
! Revision 1.11  2004/01/29 22:52:47  haselbac
! Added calls to RFLU_EnforceBoundsWrapper and updateDependentVarsMP
!
! Revision 1.10  2003/12/04 03:23:06  haselbac
! Added call to CellGradientsMP and validity check
!
! Revision 1.9  2003/11/25 21:01:45  haselbac
! Added calls to RFLU_UpdateDummyCells and ZeroResidualsMP
!
! Revision 1.8  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/10/03 20:42:07  haselbac
! Added Rocflu calls
!
! Revision 1.4  2003/10/01 23:52:09  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.3  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.2  2003/04/10 01:22:41  jblazek
! Got rid of pRegion in ViscousFluxesMP.
!
! Revision 1.1  2003/03/28 19:42:55  fnajjar
! Initial import for RocfluidMP
!
! ******************************************************************************







