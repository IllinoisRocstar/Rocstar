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
! Purpose: read and process user input for the core solver
!          and all physical modules.
!
! Description: none.
!
! Input: regions = dimensions and topology (finest grid).
!
! Output: regions = dimensions, topology and user input on all grid levels.
!
! Notes:
!
!******************************************************************************
!
! $Id: RFLO_GetUserInput.F90,v 1.6 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_GetUserInput( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInterfaces, ONLY : RFLO_UserInput, RFLO_PrintUserInput, &
    RFLO_CopyTopologyLevels, RFLO_CheckBcInput, RFLO_ReadTbcInputFile, & 
    RFLO_CheckDerivedUserInput
#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY : PLAG_UserInput
#endif
#ifdef PEUL
  USE ModInterfacesEulerian, ONLY : PEUL_UserInput
#endif
#ifdef PERI
  USE ModInterfacesPeriodic, ONLY : PERI_UserInput
#endif
#ifdef RADI
  USE ModInterfacesRadiation, ONLY : RADI_UserInput
#endif
#ifdef SPEC
  USE ModInterfacesSpecies, ONLY : SPEC_UserInput
#endif
#ifdef TURB
  USE ModInterfacesTurbulence, ONLY : TURB_UserInput
#endif
#ifdef INRT
  USE ModInterfacesInteract, ONLY : INRT_UserInput
#endif
  USE ModMPI
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... local variables
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_GetUserInput',&
  'RFLO_GetUserInput.F90' )

! input for mixture & base solver

  CALL RFLO_UserInput( regions )

! input for physical modules

#ifdef PLAG
  CALL PLAG_UserInput( regions )
#endif

#ifdef PEUL
  CALL PEUL_UserInput( regions )
#endif

#ifdef PERI
  CALL PERI_UserInput( regions )
#endif

#ifdef RADI
  CALL RADI_UserInput( regions )
#endif

#ifdef SPEC
  CALL SPEC_UserInput( regions )
#endif

#ifdef TURB
  CALL TURB_UserInput( regions )
#endif

#ifdef INRT
  CALL INRT_UserInput( regions ) ! must follow all other *_UserInput routines
#endif


  CALL RFLO_CopyTopologyLevels( regions )


  CALL RFLO_CheckBcInput( regions )


  CALL RFLO_ReadTbcInputFile( regions )



! Check user input which requires knowledge of MP module input



  CALL RFLO_CheckDerivedUserInput(regions)

! write user input

  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_MED) THEN
    CALL RFLO_PrintUserInput( regions )
  ENDIF

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_GetUserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_GetUserInput.F90,v $
! Revision 1.6  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/26 14:03:24  mtcampbe
!
! Port to NCSA Abe cluster. RFLO=1 works as-is. RFLU=1 needs some source mods
! to match Fortran name mangling.  Next checkin will make this automatic.
!
! Revision 1.4  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2008/10/23 18:20:57  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.2  2004/11/30 20:12:15  fnajjar
! Added call to RFLO_CheckDerivedUserInput
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.14  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.10  2003/09/26 21:44:28  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.9  2003/08/28 20:35:21  wasistho
! excluced ModInterfacesTurbulence,Radiation,Periodic from ModInterfaces
!
! Revision 1.8  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.7  2003/03/29 03:30:11  wasistho
! install ROCPERI
!
! Revision 1.6  2003/03/11 16:04:19  jferry
! Enclosed USE statements for multi-physics routines within ifdefs
!
! Revision 1.5  2003/03/04 22:12:34  jferry
! Initial import of Rocinteract
!
! Revision 1.4  2003/02/11 22:49:53  jferry
! Re-worked BC and TBC input routines to add multi-physics capability
!
! Revision 1.1  2003/02/11 22:30:21  jferry
! Re-worked BC and TBC input routines to add multi-physics capability
!
! Revision 1.3  2002/09/17 22:51:23  jferry
! Removed Fast Eulerian particle type
!
! Revision 1.2  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/02/25 22:36:53  jblazek
! Simplified solver initialization routine.
!
!******************************************************************************







