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
! Purpose: update conserved variable fields for RocfluidMP framework.
!
! Description: none.
!
! Input: region = data of current region,
!        iReg   = index of currect region,
!        istage = RK current stage.
!
! Output: region%levels%*%cv,rhs = new solution after one time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RkUpdateMP.F90,v 1.9 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RKUpdateMP( region,iReg,istage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModMixture,    ONLY : t_mixt
#ifdef RFLO
#ifdef PEUL
  USE ModPartEul,    ONLY : t_peul
#endif
#endif
#ifdef RFLU
  USE ModSpecies,    ONLY : t_spec

  USE RFLU_ModNSCBC, ONLY: RFLU_NSCBC_DecideHaveNSCBC
#endif
  USE ModError
  USE ModParameters
  USE ModBndPatch, ONLY: t_patch

  USE ModInterfaces, ONLY: RkUpdateGeneric

#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset

#include "Indexing.h"
#endif

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY : PLAG_RkUpdateWrapper
#endif

#ifdef PERI
  USE ModInterfacesPeriodic,   ONLY : PERI_SolutionUpdate
#endif

#ifdef TURB
  USE ModInterfacesTurbulence, ONLY : TURB_SolutionUpdate
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region
  INTEGER, INTENT(IN) :: iReg, istage

! ... local variables
#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff
#endif
  INTEGER :: ibc, iec

  LOGICAL :: peulUsed, plagUsed, specUsed, periUsed, turbUsed

  TYPE(t_mixt),   POINTER :: mixt
#ifdef RFLO
#ifdef PEUL
  TYPE(t_peul),   POINTER :: peul
#endif
#endif
#ifdef RFLU
  TYPE(t_spec),   POINTER :: spec
#endif
  TYPE(t_global), POINTER :: global

  INTEGER :: iPatch
  TYPE(t_patch), POINTER :: pPatch
!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RKUpdateMP',&
  'RkUpdateMP.F90' )

! set flags ===================================================================

  peulUsed = global%peulUsed
  plagUsed = global%plagUsed
  specUsed = global%specUsed

#ifdef PERI
  periUsed = (region%periInput%flowKind /= OFF)
#else
  periUsed = .FALSE.
#endif

#ifdef TURB
  turbUsed = (region%mixtInput%flowModel == FLOW_NAVST) .AND. &
             (region%mixtInput%turbModel /= TURB_MODEL_NONE)
#else
  turbUsed = .FALSE.
#endif

  region%irkStep = istage

! get dimensions and pointers -------------------------------------------------

#ifdef RFLO
  iLev = region%currLevel
  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)
  mixt => region%levels(iLev)%mixt
#ifdef PEUL
  peul => region%levels(iLev)%peul
#endif
#endif

#ifdef RFLU
  ibc = 1
  iec = region%grid%nCellsTot
  mixt => region%mixt
  spec => region%spec
#endif

! update solution for flow itself and all multiphysics modules ----------------

  IF ( region%mixtInput%frozenFlag .EQV. .FALSE. ) THEN 
    CALL RkUpdateGeneric(region,VAR_TYPE_CELL,istage,ibc,iec,1,CV_MIXT_NEQS, &
                         mixt%cv,mixt%cvOld,mixt%rhs,mixt%rhsSum)
  END IF ! region%mixtInput%frozenFlag

#ifdef RFLU
! loop over patches in a region and update boundary variables -------------
  DO iPatch = 1,region%grid%nPatches
    pPatch => region%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      CALL RkUpdateGeneric(region,VAR_TYPE_POINT,istage,1,pPatch%nBFaces, &
                           1,CV_MIXT_NEQS,pPatch%mixt%cv,pPatch%mixt%cvOld, &
                           pPatch%mixt%rhs,pPatch%mixt%rhsSum)

    END IF ! pPatch%bcKind
  END DO ! iPatch
#endif

#ifdef RFLO
#ifdef PEUL
  IF ( peulUsed ) THEN
    CALL RkUpdateGeneric(region,VAR_TYPE_CELL,istage,ibc,iec,1,peul%nCv, &
                         peul%cv,peul%cvOld,peul%rhs,peul%rhsSum)
  END IF ! peulUsed
#endif
#endif

#ifdef RFLU
#ifdef SPEC
  IF ( specUsed ) THEN
    CALL RkUpdateGeneric(region,VAR_TYPE_CELL,istage,ibc,iec,1, & 
                         region%specInput%nSpecies,spec%cv,spec%cvOld, &
                         spec%rhs,spec%rhsSum)
  END IF ! specUsed
#endif
#endif

#ifdef PLAG
  IF ( plagUsed ) THEN
    CALL PLAG_rkUpdateWrapper(region,iReg,istage)
  END IF ! plagUsed
#endif

#ifdef PERI
  IF ( periUsed ) THEN
    CALL PERI_SolutionUpdate(region)
  END IF ! periUsed
#endif

#ifdef TURB
  IF ( turbUsed ) THEN
    CALL TURB_SolutionUpdate(region,istage,ibc,iec)
  END IF ! turbUsed
#endif

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE RKUpdateMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RkUpdateMP.F90,v $
! Revision 1.9  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2006/08/19 15:38:33  mparmar
! Added update of Runge-Kutta scheme for boundary arrays
!
! Revision 1.6  2006/02/13 21:00:58  wasistho
! added ifdef PEUL
!
! Revision 1.5  2005/11/10 22:20:38  fnajjar
! ACH: Added IF on frozenFlag
!
! Revision 1.4  2005/04/06 02:16:49  wasistho
! mv call to PERI_CoMeanCorrection to UpdateBoundaryConditionsMP
!
! Revision 1.3  2005/03/11 04:22:57  wasistho
! commented PERI_coMeanCorrection temporarily while in testing
!
! Revision 1.2  2005/03/07 05:05:04  wasistho
! install hybrid DESSA turbulence model
!
! Revision 1.1  2004/12/01 16:51:11  haselbac
! Initial revision after changing case
!
! Revision 1.21  2004/11/17 23:44:27  wasistho
! used generic RK-update for rocturb
!
! Revision 1.20  2004/11/17 16:25:11  haselbac
! Added varType to interface of RkUpdateGeneric
!
! Revision 1.19  2004/07/28 15:29:18  jferry
! created global variable for spec use
!
! Revision 1.18  2004/07/23 22:43:15  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.17  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.16  2004/03/03 23:55:08  jferry
! Made module calls more uniform
!
! Revision 1.15  2004/03/02 21:50:30  jferry
! Changed rkInit and rkUpdate routines to call generic procedures
!
! Revision 1.14  2004/02/26 21:01:45  haselbac
! Removed ifdef RFLO around PLAG_rkUpdateWrapper
!
! Revision 1.13  2004/02/02 22:48:32  haselbac
! Added ifdef RFLO - temporary measure
!
! Revision 1.12  2003/11/25 21:01:44  haselbac
! Added rocspecies support with rkUpdateGeneric routine
!
! Revision 1.11  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.8  2003/09/26 21:45:34  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.7  2003/08/28 20:32:10  wasistho
! excluced ModInterfacesTurbulence,Radiation,Periodic from ModInterfaces
!
! Revision 1.6  2003/08/14 01:46:30  wasistho
! fixed ifdef around TURB_solutionUpdate
!
! Revision 1.5  2003/08/06 15:53:09  wasistho
! added vorticities computation
!
! Revision 1.4  2003/07/08 21:21:36  jblazek
! Modified start up procedure for dual-time stepping.
!
! Revision 1.3  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.2  2003/04/05 01:58:45  wasistho
! install ROCPERI
!
! Revision 1.1  2003/03/28 19:47:19  fnajjar
! Initial import for RocfluidMP
!
!******************************************************************************







