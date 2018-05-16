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
! Purpose: read initial flow solution from file, initialize dependent
!          variables and module specific data.
!
! Description: none.
!
! Input: regions = dimensions, user input.
!
! Output: regions = mixture and module variables.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_GetFlowSolution.F90,v 1.9 2009/10/30 15:26:29 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_GetFlowSolution( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE RFLO_ModInitialization,     ONLY : RFLO_InitNonCvSolution
  USE RFLO_ModBoundaryConditions, ONLY : RFLO_BcondInjectionInit 
  USE ModInterfaces, ONLY : RFLO_ReadSolution, MixtureProperties, &
         RFLO_LimiterReference, RFLO_GetDimensDummy, RFLO_GetCellOffset, &
         RFLO_GetBoundaryValues, RFLO_ReadRandomState

#ifdef INRT
  USE ModInterfacesInteract, ONLY : INRT_TwoDimAverage
#endif
#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY : PLAG_ReadSolution, PLAG_InitSolution
#endif
#ifdef PEUL
  USE ModInterfacesEulerian, ONLY : PEUL_ReadSolution, PEUL_InitSolution
#endif
#ifdef RADI
  USE ModInterfacesRadiation, ONLY : RADI_InitSolution
#endif
#ifdef SPEC
  USE ModInterfacesSpecies, ONLY : SPEC_ReadSolution, SPEC_InitSolution
#endif
#ifdef TURB
  USE ModInterfacesTurbulence, ONLY : TURB_InitSolution, TURB_RFLO_ReadSolution
#endif
#ifdef PERI
  USE ModInterfacesPeriodic, ONLY : PERI_InitSolution
#ifdef RFLO 
  USE PERI_ModHybridDES, ONLY : PERI_RFLO_ReadMean
#endif 
#endif
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iLev, ibc, iec, original_format

! ... local variables
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend, iCOff, ijCOff

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_GetFlowSolution',&
  'RFLO_GetFlowSolution.F90' )

! read/receive solution for all regions on this processor

  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME//'   - mixture'

#ifndef GENX
  CALL RFLO_ReadSolution( regions )
  CALL RFLO_ReadRandomState( regions )
#ifdef PLAG
  CALL PLAG_ReadSolution( regions )
#endif
#ifdef PEUL
  CALL PEUL_ReadSolution( regions )
#endif
#ifdef SPEC
  CALL SPEC_ReadSolution( regions )
#endif
#ifdef TURB
  IF (global%turbActive) THEN
    IF ((global%flowType == FLOW_UNSTEADY .AND. &
         global%timeStamp > 0._RFREAL)     .OR. &
        (global%flowType == FLOW_STEADY   .AND. &
         global%currentIter > 0)) THEN
      CALL TURB_RFLO_ReadSolution( regions )
    ENDIF
#ifdef PERI
    IF (global%flowType == FLOW_UNSTEADY) THEN
      CALL PERI_RFLO_ReadMean( regions )
    ENDIF
#endif
  ENDIF  ! turbActive
#endif
#endif

#ifdef GENX
#ifdef NATIVE_MP_IO
      original_format = global%solutFormat
      global%gridFormat  = FORMAT_ASCII
      global%solutFormat = FORMAT_ASCII
#ifdef PLAG
  CALL PLAG_ReadSolution( regions )
#endif
#ifdef PEUL
  CALL PEUL_ReadSolution( regions )
#endif
#ifdef SPEC
  CALL SPEC_ReadSolution( regions )
#endif
#ifdef TURB
  IF (global%turbActive) THEN
    IF ((global%flowType == FLOW_UNSTEADY .AND. &
         global%timeStamp > 0._RFREAL)     .OR. &
        (global%flowType == FLOW_STEADY   .AND. &
         global%currentIter > 0)) THEN
      CALL TURB_RFLO_ReadSolution( regions )
    ENDIF
#ifdef PERI
    IF (global%flowType == FLOW_UNSTEADY) THEN
      CALL PERI_RFLO_ReadMean( regions )
    ENDIF
#endif
  ENDIF  ! turbActive
#endif
#endif
#endif


! loop over all regions - get variables for physical modules

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

      CALL RFLO_InitNonCvSolution( regions(iReg) )

#ifdef INRT
      IF (global%inrtUsed) THEN
        IF (regions(iReg)%inrtInput%twoDAverage /= 0) THEN
          CALL INRT_TwoDimAverage(regions(iReg))
        END IF
      END IF ! inrtUsed
#endif

#ifdef PEUL
      IF (global%peulUsed) &
        CALL PEUL_InitSolution( iReg,regions(iReg) )
#endif
#ifdef RADI
      IF (regions(iReg)%mixtInput%radiUsed) &
        CALL RADI_InitSolution( iReg,regions(iReg) )
#endif
#ifdef SPEC
      CALL SPEC_InitSolution( iReg,regions(iReg) )
#endif
#ifdef PERI
      IF (regions(iReg)%periInput%flowKind /= OFF) THEN
        CALL PERI_InitSolution( regions,iReg )
      ENDIF
#endif

! --- calculate derived values for the mixture

      DO iLev=1,regions(iReg)%nGridLevels
        CALL RFLO_GetDimensDummy( regions(iReg),iLev,idcbeg,idcend, &
                                  jdcbeg,jdcend,kdcbeg,kdcend )
        CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )
        ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
        iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)
        CALL MixtureProperties( regions(iReg),ibc,iec,.true. )
      ENDDO

! --- convert mdot into density*velocity for injection boundaries

      CALL RFLO_BcondInjectionInit( regions(iReg) )

! --- get values for external BC`s

#ifdef GENX
      IF (global%timeStamp>0._RFREAL .AND. &
          regions(iReg)%mixtInput%externalBc) THEN    ! only if restart
        CALL RFLO_GetBoundaryValues( regions(iReg) )
      ENDIF
#else
      IF (regions(iReg)%mixtInput%externalBc) THEN
        CALL RFLO_GetBoundaryValues( regions(iReg) )
      ENDIF
#endif

! --- initialize Lagrangian particle solution field based on mixture properties

#ifdef PLAG
      IF (global%plagUsed) &
        CALL PLAG_InitSolution( iReg,regions(iReg) )
#endif

! --- initiate turbulence solution based on mixture properties

#ifdef TURB
      IF ((regions(iReg)%mixtInput%flowModel == FLOW_NAVST) .AND. &
          (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE)) THEN
        CALL TURB_InitSolution( regions(iReg) )
      ENDIF
#endif

    ENDIF     ! region on this processor and active
  ENDDO       ! iReg

! calculate reference values for limiter (upwind schene)

  CALL RFLO_LimiterReference( regions )

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_GetFlowSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_GetFlowSolution.F90,v $
! Revision 1.9  2009/10/30 15:26:29  mtcampbe
! #ifdef bugfix
!
! Revision 1.8  2009/10/26 00:19:31  mtcampbe
! Updates for completion of NATIVE_MP_IO
!
! Revision 1.7  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/03/22 03:03:01  wasistho
! added call to RFLO_InitNonCvSolution
!
! Revision 1.4  2006/01/25 01:42:14  wasistho
! call bcondInjectionInit also in Genx
!
! Revision 1.3  2005/03/07 05:04:32  wasistho
! install hybrid DESSA turbulence model
!
! Revision 1.2  2004/12/28 22:48:49  wasistho
! moved RFLO_Bcond* and RFLO_BoundaryCond* routines into RFLO_ModBoundaryConditions
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.37  2004/07/02 22:27:59  jiao
! Changed not to read solutions from native file formats if in GEN3.
!
! Revision 1.36  2004/06/19 03:29:11  wasistho
! removed argument iReg in TURB_InitSolution
!
! Revision 1.35  2004/03/11 03:31:46  wasistho
! changed rocturb nomenclature
!
! Revision 1.34  2004/03/05 22:09:02  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.33  2004/03/02 21:47:29  jferry
! Added After Update interactions
!
! Revision 1.32  2003/11/21 22:35:51  fnajjar
! Update Random Number Generator
!
! Revision 1.31  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.27  2003/09/26 21:44:28  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.26  2003/08/28 20:35:14  wasistho
! excluced ModInterfacesTurbulence,Radiation,Periodic from ModInterfaces
!
! Revision 1.25  2003/08/11 20:16:22  wasistho
! distinguish steady and unsteady in reading turb. solution
!
! Revision 1.24  2003/08/08 02:09:05  wasistho
! fixed turb. restart for steady flow
!
! Revision 1.23  2003/08/01 22:14:51  wasistho
! turbWrite to turbActive
!
! Revision 1.22  2003/07/22 02:57:30  wasistho
! prepare more accurate rocturb restart
!
! Revision 1.21  2003/07/17 01:03:05  wasistho
! initial activation rocrad
!
! Revision 1.20  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.19  2003/04/09 15:05:05  jferry
! added check that particles are used before calling PLAG routines
!
! Revision 1.18  2003/03/29 03:30:04  wasistho
! install ROCPERI
!
! Revision 1.17  2003/03/11 16:04:19  jferry
! Enclosed USE statements for multi-physics routines within ifdefs
!
! Revision 1.16  2003/02/11 22:53:19  jferry
! Initial import of Rocsmoke
!
! Revision 1.15  2002/12/04 15:38:53  f-najjar
! Change calling sequence for PLAG_InitSolution after mixture properties
!
! Revision 1.14  2002/10/16 18:56:15  jblazek
! Within GenX, BC data at t>0 are obtained in GetFlowSolution.
!
! Revision 1.13  2002/10/16 18:30:38  jblazek
! Within GenX, BC data at t=0 are updated in FlowSolver before calling
! the time-stepping routine.
!
! Revision 1.12  2002/10/04 19:33:18  jblazek
! Corrected one more bug in GenX restart.
!
! Revision 1.11  2002/10/03 21:25:39  jblazek
! Changed init. of burnig boundaries for GenX.
!
! Revision 1.10  2002/10/02 22:21:59  jiao
! Debugged GenX restart.
!
! Revision 1.9  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.8  2002/09/17 22:51:23  jferry
! Removed Fast Eulerian particle type
!
! Revision 1.7  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.6  2002/08/24 03:16:37  wasistho
! put safety within #ifdef TURB
!
! Revision 1.5  2002/08/16 21:33:48  jblazek
! Changed interface to MixtureProperties.
!
! Revision 1.4  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.3  2002/06/12 21:56:29  jblazek
! Added read/write solution for physical modules.
!
! Revision 1.2  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.1  2002/02/25 22:36:52  jblazek
! Simplified solver initialization routine.
!
!******************************************************************************







