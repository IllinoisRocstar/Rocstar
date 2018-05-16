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
! Purpose: define data types related to multi-phase interactions
!
! Description:
!
!  t_inrt_input: information common to all interactions, as well as the
!                individual interactions
!
!  t_inrt_interaction: defines an individual interaction
!
!  t_inrt_edge: an edge of an interaction, representing a transfer between
!               two phases
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModInteract.F90,v 1.11 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

MODULE ModInteract

  USE ModDataTypes
  IMPLICIT NONE

! data types ------------------------------------------------------------------

  TYPE t_inrt_input
    INTEGER :: nNodes       ! Total number of Nodes (for all Node Tag arrays)
    INTEGER :: nPlag        ! Number of Lagrangian particle constituents
    INTEGER :: nPeul        ! Number of smoke types
    INTEGER :: indMixt      ! Index of mixture (in a Node Tag array)
    INTEGER :: indPlag0     ! Zero-index of Lagrangian particle constituents
    INTEGER :: indPeul0     ! Zero-index of smoke types
    INTEGER :: indIntl      ! Index of internal Node (when it exists)
    INTEGER :: indPlagJoint ! Index used for joint momentum and energy storage
    INTEGER :: indPlagVapor ! Index used for Vapor Energy
    INTEGER :: maxConEdges  ! Maximum number of Edges over continuum inrts
    INTEGER :: maxDisEdges  ! Maximum number of Edges over discrete inrts
    LOGICAL :: defaultRead  ! Whether data has been read from input deck
    LOGICAL :: consistent   ! Whether active phases form a consistent system
    LOGICAL :: computeAux   ! Whether to compute any auxillary Edge quantities
    INTEGER :: twoDAverage  ! Whether to average k direction
    INTEGER, POINTER :: globActiveness(:) ! Global Activeness for each Node
    INTEGER, POINTER :: globPermission(:) ! Global Permission for each Node
    TYPE(t_inrt_interact), POINTER :: inrts(:)  ! Individual interactions
  END TYPE t_inrt_input

  TYPE t_inrt_interact
    CHARACTER(CHRLEN) :: name ! name of interaction
    LOGICAL :: used           ! Whether this interaction is used
    LOGICAL :: pclsUsed       ! see below
    INTEGER :: order          ! see below
    INTEGER :: nIntl          ! see below
    INTEGER :: nInputEdges    ! see below
    INTEGER :: nSwitches      ! # of Switches in interaction
    INTEGER :: nData          ! # of Data     in interaction
    INTEGER :: nEdges         ! # of Edges    in interaction
    INTEGER,      POINTER :: switches(:)   ! Integer auxillary information
    REAL(RFREAL), POINTER :: data(:)       ! Real    auxillary information
    INTEGER,      POINTER :: activeness(:) ! Activeness for each Node
    INTEGER,      POINTER :: permission(:) ! Permission for each Node
    TYPE(t_inrt_edge), POINTER :: edges(:) ! List of Edges
  END TYPE t_inrt_interact

! The following parameters must be set for each interaction:
!
! pclsUsed:    default value: .TRUE.
!              .TRUE.  for an interaction involving particles.
!              .FALSE. for an interaction involving only cell-based quantities.
!
! order:       default value: 0
!              0 for a normal interaction.
!              Interactions with a given value are guaranteed to be computed
!              after all those with smaller values.
!              Allowed values are any integers except those near HUGE or -HUGE.
!
! nIntl:       default value: 0
!              0 for an interaction without an Internal Node.
!              1 for an interaction with an Internal Node.
!              0 and 1 are the only allowed values.
!
! nInputEdges: default value: 0
!              Number of edges that precede the internal Node, if it exists.

  TYPE t_inrt_edge
    INTEGER :: tEdge    ! Type of Edge
    INTEGER :: iNode(2) ! Indices of Nodes on either end
    INTEGER :: token(2) ! Permission Tokens on either end of Edge
  END TYPE t_inrt_edge

END MODULE ModInteract

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInteract.F90,v $
! Revision 1.11  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2004/07/28 15:42:12  jferry
! deleted defunct constructs: useDetangle, useSmokeDrag, useSmokeHeatTransfer
!
! Revision 1.8  2004/07/26 17:05:50  fnajjar
! moved allocation of inrtSources into Rocpart
!
! Revision 1.7  2004/03/05 22:09:01  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2004/03/02 21:47:28  jferry
! Added After Update interactions
!
! Revision 1.5  2003/04/03 21:10:17  jferry
! implemented additional safety checks for rocinteract
!
! Revision 1.4  2003/04/02 22:32:03  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.3  2003/03/24 23:30:52  jferry
! overhauled rocinteract to allow interaction design to use user input
!
! Revision 1.2  2003/03/11 16:09:39  jferry
! Added comments
!
! Revision 1.1  2003/03/04 22:12:34  jferry
! Initial import of Rocinteract
!
!******************************************************************************






