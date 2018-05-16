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
! Purpose: Defines the interaction HeatTransferNonBurn
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: modifies region%inrtInput%inrts
!
! Notes:
!
!   Whereas INRT_Initialize sets up everything generic to all interactions,
!   this routine sets up things specific to this interaction.
!
!   In particular, this is where the designer
!
!   (a) gives values for all the interactions parameters (or leaves them
!       as default values given in INRT_Initialize).  These parameters are
!       described in ModInteract.
!
!   (b) defines the Edges of the interaction.  Here the type of an Edge is
!       given (i.e., whether it transports mass, momentum, or energy), as well
!       as the indices of the Nodes at either end.  For each momentum Edge
!       he must call the routine INRT_FinishMomentumEdge.
!
!******************************************************************************
!
! $Id: INRT_DefineHeatTransferNonBurn.F90,v 1.3 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_DefineHeatTransferNonBurn( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInteract
  USE ModError
  USE INRT_ModParameters

  USE INRT_ModInterfaces, ONLY : INRT_AllocateAuxillary
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: nEdges,nPlag01

  TYPE(t_inrt_input),    POINTER :: input
  TYPE(t_inrt_interact), POINTER :: inrt
  TYPE(t_inrt_edge),     POINTER :: edge
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: INRT_DefineHeatTransferNonBurn.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_DefineHeatTransferNonBurn',&
  'INRT_DefineHeatTransferNonBurn.F90' )

! begin -----------------------------------------------------------------------

! allocate memory for edges, switches and data

  input => region%inrtInput
  inrt  => input%inrts(INRT_TYPE_HTRANSNB)

  inrt%name = "Heat_Transfer_Non_Burn"

  nPlag01 = 1
  IF (input%nPlag == 0) nPlag01 = 0

  nEdges = nPlag01 * INRT_HTRANSNB_NEDGES

  CALL INRT_AllocateAuxillary(global,inrt,nEdges, &
    INRT_SWI_HTRANSNB_TOTAL,INRT_DAT_HTRANSNB_TOTAL)

! set parameters for this interaction:  see comment in ModInteract

! ** for this interaction, all parameters keep their default values **

! define Edges

  IF (nEdges > 0) THEN

    edge => inrt%edges(INRT_HTRANSNB_L_ENER_G)

    edge%tEdge    = INRT_EDGE_ENER
    edge%iNode(1) = input%indPlagJoint
    edge%iNode(2) = input%indMixt

  END IF ! nEdges

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_DefineHeatTransferNonBurn

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_DefineHeatTransferNonBurn.F90,v $
! Revision 1.3  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:56:23  fnajjar
! Initial revision after changing case
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
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************







