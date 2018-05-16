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
! Purpose: set initial solution field for first RK stage.
!
! Description: none.
!
! Input: region = data of current region,
!        istage  = current RK stage.
!
! Output: region%levels%*%cvOld,diss = initializations for new time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RkInitMP.F90,v 1.6 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RKInitMP( region, istage )

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
#endif
  USE ModError
  USE ModParameters
  USE ModBndPatch, ONLY: t_patch 

  USE ModInterfaces, ONLY : RkInitGeneric
#ifdef RFLU
  USE ModInterfaces, ONLY : RkInitSD

  USE RFLU_ModNSCBC, ONLY: RFLU_NSCBC_DecideHaveNSCBC
#endif
#ifdef RFLO
  USE ModInterfaces, ONLY : ScaleGridSpeeds, RFLO_GetDimensDummy, &
                            RFLO_GetCellOffset

#include "Indexing.h"
#endif

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY : PLAG_RkInit
#endif

#ifdef TURB
  USE ModInterfacesTurbulence, ONLY : TURB_RkInit
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region
  INTEGER,        INTENT(IN)            :: istage

! ... local variables
#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff
#endif
  INTEGER :: ibc, iec

  LOGICAL :: moveGrid, peulUsed, specUsed, plagUsed, turbUsed

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

  CALL RegisterFunction( global,'RKInitMP',&
  'RkInitMP.F90' )

! set flags ===================================================================

  peulUsed = global%peulUsed
  plagUsed = global%plagUsed
  specUsed = global%specUsed

#ifdef TURB
  turbUsed = (region%mixtInput%flowModel == FLOW_NAVST) .AND. &
             (region%mixtInput%turbModel /= TURB_MODEL_NONE)
#else
  turbUsed = .FALSE.
#endif

  moveGrid = region%mixtInput%moveGrid

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

! initialize for flow itself and all multiphysics modules ---------------------

  CALL RkInitGeneric(region,istage,ibc,iec,1,CV_MIXT_NEQS, &
                     mixt%cv,mixt%cvOld,mixt%diss)

#ifdef RFLU
! loop over patches in a region and initialize boundary variables -------------
  DO iPatch = 1,region%grid%nPatches
    pPatch => region%patches(iPatch)

    IF ( pPatch%bcKind == BC_KIND_NSCBC ) THEN
      CALL RkInitGeneric(region,istage,1,pPatch%nBFaces,1,CV_MIXT_NEQS, &
                         pPatch%mixt%cv,pPatch%mixt%cvOld,pPatch%mixt%rhs)
    END IF ! pPatch%bcKind
  END DO ! iPatch
#endif

#ifdef RFLO
#ifdef PEUL
  IF ( peulUsed ) THEN
    CALL RkInitGeneric(region,istage,ibc,iec,1,peul%nCv, &
                       peul%cv,peul%cvOld,peul%diss)
  END IF ! peulUsed
#endif
#endif

#ifdef RFLU
#ifdef SPEC
  IF ( specUsed ) THEN
    CALL RkInitGeneric(region,istage,ibc,iec,1,region%specInput%nSpecies, &
                       spec%cv,spec%cvOld,spec%diss)
  END IF ! specUsed
#endif
  SELECT CASE ( region%mixtInput%indSd )

  CASE (0)
    CALL RkInitSD(region,  0,  1,1,SD_ZMOM-SD_XMOM+1,mixt%sd)

  CASE (1)
    CALL RkInitSD(region,ibc,iec,1,SD_ZMOM-SD_XMOM+1,mixt%sd)

  CASE DEFAULT
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)

  END SELECT ! indSd
#endif

#ifdef PLAG
  IF ( plagUsed ) THEN
    CALL PLAG_RkInit(region,istage)
  END IF ! plagUsed
#endif

#ifdef TURB
  IF ( turbUsed ) THEN
    CALL TURB_RkInit(region,istage)
  ENDIF ! turbUsed
#endif

! scale grid speeds (volume change) -------------------------------------------

#ifdef RFLO
  IF ( moveGrid .EQV. .TRUE. ) THEN
    CALL ScaleGridSpeeds(region)
  END IF ! moveGrid
#endif

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE rkInitMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RkInitMP.F90,v $
! Revision 1.6  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/08/19 15:38:31  mparmar
! Added initialization of Runge-Kutta scheme for boundary arrays
!
! Revision 1.3  2006/02/13 21:00:51  wasistho
! added ifdef PEUL
!
! Revision 1.2  2005/03/31 16:30:04  haselbac
! Replaced nSd by parameters
!
! Revision 1.1  2004/12/01 16:51:04  haselbac
! Initial revision after changing case
!
! Revision 1.19  2004/07/30 22:47:33  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.18  2004/07/28 15:29:18  jferry
! created global variable for spec use
!
! Revision 1.17  2004/07/23 22:43:15  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.16  2004/04/14 02:06:42  haselbac
! No longer call ScaleGridSpeeds for RFLU
!
! Revision 1.15  2004/03/05 22:09:00  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.14  2004/03/03 23:55:08  jferry
! Made module calls more uniform
!
! Revision 1.13  2004/03/02 21:50:30  jferry
! Changed rkInit and rkUpdate routines to call generic procedures
!
! Revision 1.12  2004/02/26 21:13:06  wasistho
! changed TURB_ransRkInit to TURB_rkInit
!
! Revision 1.11  2004/02/26 21:01:43  haselbac
! Removed ifdef RFLO around PLAG_rkInitMP
!
! Revision 1.10  2004/02/02 22:48:21  haselbac
! Added ifdef RFLO - temporary measure
!
! Revision 1.9  2003/11/25 21:01:43  haselbac
! Added rocspecies support with rkInitGeneric routine
!
! Revision 1.8  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/10/03 20:12:53  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.4  2003/10/01 23:52:09  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.3  2003/09/26 21:45:34  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.2  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.1  2003/03/28 19:46:59  fnajjar
! Initial import for RocfluidMP
!
!******************************************************************************







