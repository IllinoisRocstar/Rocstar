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
! Purpose: initializes data structures for interactions
!
! Description: specifies data common to all interactions, then calls routines
!              specific to each interaction
!
! Input: region = data of current region.
!
! Output: fills in region%inrtInput
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_Initialize.F90,v 1.3 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_Initialize( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInteract
  USE ModError
  USE INRT_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: iInrt,iPeul

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: nNodes,nPlag,nPeul
  INTEGER :: indMixt,indPlag0,indPeul0,indIntl
  INTEGER :: indPlagJoint,indPlagVapor,iPlagJoint
  INTEGER :: errorFlag

  TYPE(t_inrt_input),    POINTER :: input
  TYPE(t_inrt_interact), POINTER :: inrt
  TYPE(t_inrt_edge),     POINTER :: edge
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_Initialize.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_Initialize',&
  'INRT_Initialize.F90' )

! begin -----------------------------------------------------------------------

  ALLOCATE( region%inrtInput,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  input => region%inrtInput

! Determine the quantities needed to index the Activeness and Permission arrays

  nPlag = 0

#ifdef PLAG
  IF (global%plagUsed) nPlag = region%plagInput%nCont
#endif

! it does not matter which particle index is used to store momentum and energy
! that is common to all particle constituents.  We choose the first index:

  iPlagJoint = 1

  nPeul = 0

#ifdef RFLO
#ifdef PEUL
  IF (global%peulUsed) nPeul = region%peulInput%nPtypes
#endif
#endif

#ifdef RFLU
#ifdef SPEC
  IF (global%specUsed) nPeul = region%specInput%nSpecies
#endif
#endif

! Number of Nodes:
! 1     for Gas
! nPlag for particle constituent masses
! 1     for particle vapor energy
! nPeul for smoke masses
! 1     for Internal Node

  nNodes = 1 + nPlag + 1 + nPeul + 1

! Indices for Nodes

  indMixt  = 1
  indPlag0 = indMixt
  indPeul0 = indPlag0 + nPlag + 1
  indIntl  = indPeul0 + nPeul + 1

  indPlagJoint = indPlag0 + iPlagJoint
  indPlagVapor = indPlag0 + nPlag + 1

! Store all the above information in region%inrtInput

  input%nNodes       = nNodes
  input%nPlag        = nPlag
  input%nPeul        = nPeul
  input%indMixt      = indMixt
  input%indPlag0     = indPlag0
  input%indPeul0     = indPeul0
  input%indIntl      = indIntl
  input%indPlagJoint = indPlagJoint
  input%indPlagVapor = indPlagVapor

! indicate that the maximum number of edges in a diagram is 0 so far

  input%maxConEdges  = 0 ! max over diagrams indexed by cells
  input%maxDisEdges  = 0 ! max over diagrams indexed by particles

! indicate that an INRT_DEFAULT section has not yet been read

  input%defaultRead  = .FALSE.

! active phases are assumed to form a consistent system unless shown otherwise

  input%consistent   = .TRUE.

! compute auxillary quantities by default (allow User override for testing)

  input%computeAux   = .TRUE.

! allocate globActiveness and globPermission

  ALLOCATE( input%globActiveness(nNodes),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( input%globPermission(nNodes),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! give proper default values

  input%globActiveness = INRT_ACT_ACTIVE
  input%globPermission = INRT_PERM_PALL

  input%globActiveness(indIntl) = INRT_ACT_BAD ! undefined for Internal Node

! allocate interactions

  ALLOCATE( input%inrts(INRT_TYPE_TOTAL),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  DO iInrt = 1,INRT_TYPE_TOTAL

    inrt => input%inrts(iInrt)

! - set default values for interactions

    inrt%name        = 'Not defined yet'
    inrt%used        = .FALSE.
    inrt%pclsUsed    = .TRUE.
    inrt%order       = 0
    inrt%nIntl       = 0
    inrt%nInputEdges = 0
    inrt%nSwitches   = 0
    inrt%nData       = 0
    inrt%nEdges      = 0

! - allocate activeness and permission

    ALLOCATE( inrt%activeness(nNodes),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    ALLOCATE( inrt%permission(nNodes),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! - give proper default values

    inrt%activeness = INRT_ACT_BAD  ! require to be overwritten later
    inrt%permission = INRT_PERM_BAD ! require to be overwritten later

! - nullify other pointers initially

    NULLIFY( inrt%switches )
    NULLIFY( inrt%data )
    NULLIFY( inrt%edges )

  END DO ! iInrt

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_Initialize

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_Initialize.F90,v $
! Revision 1.3  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:56:28  fnajjar
! Initial revision after changing case
!
! Revision 1.10  2004/07/28 15:29:20  jferry
! created global variable for spec use
!
! Revision 1.9  2004/07/27 21:30:00  jferry
! integrated maxConEdges and maxDisEdges variables more fully
!
! Revision 1.8  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.7  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2004/03/02 21:47:29  jferry
! Added After Update interactions
!
! Revision 1.5  2003/04/03 21:10:18  jferry
! implemented additional safety checks for rocinteract
!
! Revision 1.4  2003/04/02 22:32:04  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.3  2003/03/24 23:30:52  jferry
! overhauled rocinteract to allow interaction design to use user input
!
! Revision 1.2  2003/03/11 16:09:39  jferry
! Added comments
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************







