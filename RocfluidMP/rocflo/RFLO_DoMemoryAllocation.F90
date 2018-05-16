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
! Purpose: allocate memory required for all physical modules.
!
! Description: none.
!
! Input: regions = dimensions and user input.
!
! Output: regions = allocated flow and work variables.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_DoMemoryAllocation.F90,v 1.3 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_DoMemoryAllocation( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_AllocateMemory, AllocateMemoryWork, &
                            RFLO_AllocateDataBuffers
#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY : PLAG_AllocateMemory,PLAG_AllocateDataBuffers
#endif
#ifdef PEUL
  USE ModInterfacesEulerian, ONLY : PEUL_AllocateMemory, PEUL_AllocateDataBuffers
#endif
#ifdef PERI
  USE ModInterfacesPeriodic, ONLY : PERI_AllocateMemory
#endif
#ifdef RADI
  USE ModInterfacesRadiation, ONLY : RADI_AllocateMemory
#endif
#ifdef SPEC
  USE ModInterfacesSpecies, ONLY : SPEC_AllocateMemory
#endif
#ifdef TURB
  USE ModInterfacesTurbulence, ONLY : TURB_AllocateMemory, &
                                      TURB_RFLO_RansAllocDataBuffers
#endif
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  INTEGER :: errorFlag

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_DoMemoryAllocation',&
  'RFLO_DoMemoryAllocation.F90' )

! loop over all regions

  global%nRequests = 0

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

      CALL RFLO_AllocateMemory( regions(iReg) )         ! mixture variables

      CALL AllocateMemoryWork( regions(iReg) )          ! work arrays

! --- allocate buffers for data exchange between regions

      CALL RFLO_AllocateDataBuffers( regions,iReg )

! --- memory required by physical modules

#ifdef PLAG
      IF (global%plagUsed) THEN
        CALL PLAG_AllocateMemory( regions(iReg) )

! --- allocate buffers for data exchange between regions

        CALL PLAG_AllocateDataBuffers( regions, iReg )
      END IF ! plagUsed
#endif

#ifdef PEUL
      IF (global%peulUsed) THEN
        CALL PEUL_AllocateMemory( regions(iReg) )
        CALL PEUL_AllocateDataBuffers( regions, iReg )
      ENDIF ! peulUsed
#endif

#ifdef RADI
      IF (regions(iReg)%mixtInput%radiUsed) THEN
        CALL RADI_AllocateMemory( regions(iReg) )
      ENDIF
#endif

#ifdef SPEC
      CALL SPEC_AllocateMemory( regions(iReg) )
#endif

#ifdef TURB
      IF ((regions(iReg)%mixtInput%flowModel == FLOW_NAVST) .AND. &
          (regions(iReg)%mixtInput%turbModel /= TURB_MODEL_NONE)) THEN
        CALL TURB_AllocateMemory( regions(iReg) )
        CALL TURB_RFLO_RansAllocDataBuffers( regions, iReg )
      ENDIF
#endif

#ifdef PERI
      IF (regions(iReg)%periInput%flowKind /= OFF) THEN
        CALL PERI_AllocateMemory( regions(iReg) )
      ENDIF
#endif

    ENDIF     ! region on this processor and active
  ENDDO       ! iReg

! allocate array for send requests

#ifdef MPI
   ALLOCATE( global%requests(global%nRequests),stat=errorFlag )
   global%error = errorFlag
   IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
#endif

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_DoMemoryAllocation

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_DoMemoryAllocation.F90,v $
! Revision 1.3  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.26  2004/07/27 21:27:14  jferry
! removed rocinteract allocation routines (moved to rocpart)
!
! Revision 1.25  2004/03/11 03:30:54  wasistho
! changed rocturb nomenclature
!
! Revision 1.24  2004/03/05 22:09:02  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.23  2004/03/02 21:49:22  jferry
! Added inrtUsed flag to mixture data structure
!
! Revision 1.22  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.18  2003/10/03 20:18:24  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.17  2003/09/26 21:44:28  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.16  2003/08/28 20:35:06  wasistho
! excluced ModInterfacesTurbulence,Radiation,Periodic from ModInterfaces
!
! Revision 1.15  2003/07/17 01:02:33  wasistho
! initial activation rocrad
!
! Revision 1.14  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.13  2003/04/09 15:05:05  jferry
! added check that particles are used before calling PLAG routines
!
! Revision 1.12  2003/04/09 14:15:10  fnajjar
! Added call to PEUL_allocateDataBuffers
!
! Revision 1.11  2003/03/29 03:29:55  wasistho
! install ROCPERI
!
! Revision 1.10  2003/03/11 16:04:19  jferry
! Enclosed USE statements for multi-physics routines within ifdefs
!
! Revision 1.9  2003/03/04 22:12:34  jferry
! Initial import of Rocinteract
!
! Revision 1.8  2003/02/11 22:53:19  jferry
! Initial import of Rocsmoke
!
! Revision 1.7  2003/01/13 18:58:17  f-najjar
! Added PLAG_allocateDataBuffers
!
! Revision 1.6  2002/10/25 18:36:47  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.5  2002/09/17 22:51:23  jferry
! Removed Fast Eulerian particle type
!
! Revision 1.4  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.3  2002/08/24 03:16:13  wasistho
! put safety within #ifdef TURB
!
! Revision 1.2  2002/03/29 23:15:22  jblazek
! Corrected bug in MPI send.
!
! Revision 1.1  2002/02/25 22:36:52  jblazek
! Simplified solver initialization routine.
!
!******************************************************************************







