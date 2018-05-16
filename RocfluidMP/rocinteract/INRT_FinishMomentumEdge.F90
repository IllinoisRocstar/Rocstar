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
! Purpose: Finishes definition for the special case of a momentum Edge.
!
! Description: none.
!
! Input:  iXEdge:  index of the momentum Edge (which is its x-component)
!         iEnd:    end of Edge to place an Insulate Token on (iEnd = 1 or 2)
!
! Output: inrt:    pointer to the updated interaction data structure
!
! Notes:
!
!   Defining an Edge requires that two things be specified:
!   the type of Edge (mass, momentum, or energy), and
!   the indices of the Nodes on both ends.
!
!   Momentum Edges require additional set-up, which this routine provides.
!   First, they require that two dummy momentum Edges be created (because
!   routines that store primary source terms can only store scalar quantities,
!   so there needs to be extra room to send the vector momentum).  Second,
!   they require that an Insulate Token be placed at one end as part of the
!   definition of the interaction.
!
!******************************************************************************
!
! $Id: INRT_FinishMomentumEdge.F90,v 1.3 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_FinishMomentumEdge(global,inrt,iXEdge,iEnd)

  USE ModDataTypes
  USE ModGlobal,   ONLY : t_global
  USE ModInteract, ONLY : t_inrt_interact
  USE ModError
  USE INRT_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_global),        POINTER    :: global
  TYPE(t_inrt_interact), POINTER    :: inrt
  INTEGER,               INTENT(IN) :: iXEdge,iEnd

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: INRT_FinishMomentumEdge.F90,v $ $Revision: 1.3 $'

  CALL RegisterFunction( global,'INRT_FinishMomentumEdge',&
  'INRT_FinishMomentumEdge.F90' )

! begin -----------------------------------------------------------------------

! set permission Token (Permit Mass and Momemtum) on specified end

  inrt%edges(iXEdge)%token(iEnd) = INRT_PERM_PMOME

! fill in data for dummy Edges

  inrt%edges(iXEdge+1)%tEdge = INRT_EDGE_MOME_DUM
  inrt%edges(iXEdge+1)%iNode = inrt%edges(iXEdge)%iNode
  inrt%edges(iXEdge+1)%token = INRT_PERM_BLOCK

  inrt%edges(iXEdge+2)%tEdge = INRT_EDGE_MOME_DUM
  inrt%edges(iXEdge+2)%iNode = inrt%edges(iXEdge)%iNode
  inrt%edges(iXEdge+2)%token = INRT_PERM_BLOCK

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_FinishMomentumEdge

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_FinishMomentumEdge.F90,v $
! Revision 1.3  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:56:27  fnajjar
! Initial revision after changing case
!
! Revision 1.3  2003/04/02 22:32:04  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.2  2003/03/11 16:09:39  jferry
! Added comments
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************







