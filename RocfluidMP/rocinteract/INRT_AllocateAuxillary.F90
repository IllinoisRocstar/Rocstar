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
! Purpose: Allocates data structures for an individual interaction
!
! Description: none.
!
! Input: global    = global data
!        inrt      = interaction
!        nEdges    = number of Edges in interaction
!        nSwitches = number of integer parameters for interaction
!        nData     = number of real parameters for interaction
!
! Output: sets value of inrt%nEdges,
!         allocates and initializes inrt%edges, inrt%switches, and inrt%data
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_AllocateAuxillary.F90,v 1.3 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_AllocateAuxillary( global,inrt,nEdges,nSwitches,nData )

  USE ModDataTypes
  USE ModGlobal,   ONLY : t_global
  USE ModInteract, ONLY : t_inrt_interact
  USE ModError
  USE INRT_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_global),        POINTER    :: global
  TYPE(t_inrt_interact), POINTER    :: inrt
  INTEGER,               INTENT(IN) :: nEdges,nSwitches,nData

! ... loop variables
  INTEGER :: iEdge

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER           :: errorFlag

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: INRT_AllocateAuxillary.F90,v $ $Revision: 1.3 $'

  CALL RegisterFunction( global,'INRT_AllocateAuxillary',&
  'INRT_AllocateAuxillary.F90' )

! begin -----------------------------------------------------------------------

  inrt%nSwitches = nSwitches
  inrt%nData     = nData
  inrt%nEdges    = nEdges

  IF (nSwitches > 0) THEN

    ALLOCATE( inrt%switches(nSwitches),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    inrt%switches = -1

  END IF ! nSwitches

  IF (nData > 0) THEN

    ALLOCATE( inrt%data(nData),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    inrt%data = -1._RFREAL

  END IF ! nData

  IF (nEdges > 0) THEN

    ALLOCATE( inrt%edges(nEdges),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    DO iEdge = 1,nEdges
      inrt%edges(iEdge)%tEdge = INRT_EDGE_BAD
      inrt%edges(iEdge)%iNode = -1
      inrt%edges(iEdge)%token = INRT_PERM_PALL
    END DO ! iEdge

  END IF ! nEdges

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_AllocateAuxillary

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_AllocateAuxillary.F90,v $
! Revision 1.3  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:56:08  fnajjar
! Initial revision after changing case
!
! Revision 1.4  2003/04/03 21:10:17  jferry
! implemented additional safety checks for rocinteract
!
! Revision 1.3  2003/04/02 22:32:03  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.2  2003/03/11 16:09:39  jferry
! Added comments
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************







