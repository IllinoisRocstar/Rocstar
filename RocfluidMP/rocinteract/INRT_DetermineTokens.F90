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
! Purpose: puts permission Tokens on Edges
!
! Description: sets values of Tokens based in several criteria:
!
!   (a) the Permission level of Nodes
!   (b) the relative Activeness of the Nodes at either end
!   (c) if it is an upwind Node of a mass Edge
!   (d) if the Edge contains an internal Node
!
! Input: region = region data
!        inrt = interaction
!
! Output: modifies inrt
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_DetermineTokens.F90,v 1.3 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_DetermineTokens( region,inrt )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInteract
  USE ModError
  USE INRT_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_region),        INTENT(INOUT) :: region
  TYPE(t_inrt_interact), POINTER       :: inrt

! ... loop variables
  INTEGER :: iEdge

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_inrt_input), POINTER :: input
  TYPE(t_inrt_edge),  POINTER :: edge
  TYPE(t_global),     POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_DetermineTokens.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_DetermineTokens',&
  'INRT_DetermineTokens.F90' )

! begin -----------------------------------------------------------------------

  input => region%inrtInput

  DO iEdge=1,inrt%nEdges

    edge => inrt%edges(iEdge)

! - Permission Tokens already placed on dummy Edges

    IF (edge%tEdge == INRT_EDGE_MOME_DUM) CYCLE

! - For Ghost mass Edge, set downwind permission Token to 0 (Block), and
! - Upwind Token to either 1 (Permit Mass) or 0 (Block).
! - Note that activeness plays no role for Ghost mass Edges.

    IF (edge%tEdge == INRT_EDGE_MASS_GHO) THEN

      edge%token(1) = MIN(INRT_PERM_PMASS,inrt%permission(edge%iNode(1)))
      edge%token(2) = INRT_PERM_BLOCK
      CYCLE

    END IF ! INRT_EDGE_MASS_GHO

! - Decrease permission Tokens to Permission level of corresponding Nodes

    edge%token(1) = MIN(edge%token(1),inrt%permission(edge%iNode(1)))
    edge%token(2) = MIN(edge%token(2),inrt%permission(edge%iNode(2)))

! - If Nodes on either end of Edge differ in Activeness, decrease the
! - permission Token on the more active end to 0 (Block)

    IF (inrt%activeness(edge%iNode(1)) > inrt%activeness(edge%iNode(2))) &
      edge%token(1) = MIN(edge%token(1),INRT_PERM_BLOCK)

    IF (inrt%activeness(edge%iNode(2)) > inrt%activeness(edge%iNode(1))) &
      edge%token(2) = MIN(edge%token(2),INRT_PERM_BLOCK)

! - Decrease permission Token of the upwind end of a mass Edge to 1
! - (Permit Mass only)

    IF (edge%tEdge == INRT_EDGE_MASS) &
      edge%token(1) = MIN(edge%token(1),INRT_PERM_PMASS)

! - Decrease permission Token of the upwind end of an Edge if there is
! - an Internal Node there

    IF (edge%iNode(1) == input%indIntl) &
      edge%token(1) = MIN(edge%token(1),INRT_PERM_BLOCK)

! - Decrease permission Token to 0 (Block) if it is equivalent to Block
! - for its Edge type

    SELECT CASE (edge%tEdge)

    CASE (INRT_EDGE_MOME)
      IF (edge%token(1) < INRT_PERM_PMOME) edge%token(1) = INRT_PERM_BLOCK
      IF (edge%token(2) < INRT_PERM_PMOME) edge%token(2) = INRT_PERM_BLOCK

    CASE (INRT_EDGE_ENER)
      IF (edge%token(1) < INRT_PERM_PALL ) edge%token(1) = INRT_PERM_BLOCK
      IF (edge%token(2) < INRT_PERM_PALL ) edge%token(2) = INRT_PERM_BLOCK

    END SELECT ! edge%tEdge

  END DO ! iEdge

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_DetermineTokens

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_DetermineTokens.F90,v $
! Revision 1.3  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:56:26  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2003/09/19 20:35:26  jferry
! Implemented oxidizer species for burning interaction
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







