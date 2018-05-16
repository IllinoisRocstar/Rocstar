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
! Purpose: Defines the interaction Burning
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
! $Id: INRT_DefineBurning.F90,v 1.5 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_DefineBurning( region,matIndIn,matIndOut,matIndOx, &
                               oxUsed,plagOutExists )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInteract
  USE ModError
  USE ModParameters
  USE INRT_ModParameters

  USE INRT_ModInterfaces, ONLY : INRT_AllocateAuxillary
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region
  INTEGER,        INTENT(IN)    :: matIndIn,matIndOut,matIndOx
  LOGICAL,        INTENT(INOUT) :: oxUsed
  LOGICAL,        INTENT(OUT)   :: plagOutExists

! ... loop variables
  INTEGER :: iPlag,iPeul,iPeulOutEdge

! ... local variables
  INTEGER, PARAMETER :: NPEULOUTEDGES_MAX = 10

  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: nPlag,nPeul,nEdges,nPlagInEdges,nPlagOutEdges,nPeulOutEdges
  INTEGER :: nPeulOxEdges,matIndPlag,matIndPeul
  INTEGER :: iPlagIn,iPlagOut,iPeulOx,iPeulOxEdge
  INTEGER :: iPeulArr(NPEULOUTEDGES_MAX)

  TYPE(t_inrt_input),    POINTER :: input
  TYPE(t_inrt_interact), POINTER :: inrt
  TYPE(t_inrt_edge),     POINTER :: edge
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_DefineBurning.F90,v $ $Revision: 1.5 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_DefineBurning',&
  'INRT_DefineBurning.F90' )

! begin -----------------------------------------------------------------------

! allocate memory for edges, switches and data

  input => region%inrtInput
  inrt  => input%inrts(INRT_TYPE_BURNING)

  inrt%name = "Burning"

  nPlag = input%nPlag
  nPeul = input%nPeul

! loop over particle types looking for those with material index
! matIndIn or matIndOut

  nPlagInEdges  = 0
  nPlagOutEdges = 0
  DO iPlag = 1,nPlag

#ifdef PLAG
    matIndPlag = region%plagInput%materialIndex(iPlag)
#endif

    IF (matIndPlag == matIndIn) THEN

      nPlagInEdges = nPlagInEdges + 1
      iPlagIn = iPlag

    ELSE IF (matIndPlag == matIndOut) THEN

      nPlagOutEdges = nPlagOutEdges + 1
      iPlagOut = iPlag

    END IF ! matIndPlag

  END DO ! iPlag

  IF (nPlagInEdges > 1 .OR. nPlagOutEdges > 1) &
    CALL ErrorStop( global,ERR_INRT_MULTPLAGMAT,__LINE__ )

  IF (nPlagInEdges < 1) &
    CALL ErrorStop( global,ERR_INRT_MISSPLAGMAT,__LINE__ )

  plagOutExists = .TRUE.
  IF (nPlagOutEdges < 1) THEN
    plagOutExists = .FALSE. ! flag that no transfer can occur for this Edge
    iPlagOut = iPlagIn      ! therefore the Node index used does not matter
  END IF ! nPlagOutEdges

! loop over smoke types looking for those with material index matIndOut
! or, if oxidizer is used, those with index matIndOx

  iPeulOutEdge = 0
  iPeulOxEdge  = 0
  DO iPeul = 1,nPeul

#ifdef RFLO
#ifdef PEUL
    matIndPeul = region%peulInput%ptypes(iPeul)%material%index
#endif
#endif
#ifdef RFLU
    matIndPeul = region%specInput%specType(iPeul)%pMaterial%index
#endif

    IF (matIndPeul == matIndOut) THEN

      iPeulOutEdge = iPeulOutEdge + 1

      IF (iPeulOutEdge > NPEULOUTEDGES_MAX) &
        CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

      iPeulArr(iPeulOutEdge) = iPeul

    ELSE IF (oxUsed) THEN

      IF (matIndPeul == matIndOx) THEN

        iPeulOxEdge = iPeulOxEdge + 1

        IF (iPeulOxEdge > 1) &
          CALL ErrorStop( global,ERR_INRT_ONLY1,__LINE__ )

        iPeulOx = iPeul

      END IF ! matIndPeul

    END IF ! oxUsed

  END DO ! iPeul

  nPeulOutEdges = iPeulOutEdge
  nPeulOxEdges  = iPeulOxEdge
  nEdges = INRT_BURNING_NEDGES0 + nPeulOutEdges + nPeulOxEdges

  IF (oxUsed .AND. nPeulOxEdges < 1) THEN
    WRITE(STDOUT,'(A)') SOLVER_NAME//'### INRT_WARNING: no oxidizer smoke '// &
      'type: setting OX_USED = NO'
    oxUsed = .FALSE.
  END IF

  CALL INRT_AllocateAuxillary(global,inrt,nEdges, &
    INRT_SWI_BURNING_TOTAL,INRT_DAT_BURNING_TOTAL0 + nPeulOutEdges)

! set parameters for this interaction:  see comment in ModInteract

  inrt%nIntl       = 1
  inrt%nInputEdges = 2 + nPeulOxEdges

! define Edges

  edge => inrt%edges(INRT_BURNING_G_MASS_X)

  edge%tEdge    = INRT_EDGE_MASS
  edge%iNode(1) = input%indMixt
  edge%iNode(2) = input%indIntl

  edge => inrt%edges(INRT_BURNING_L_MASS_X)

  edge%tEdge    = INRT_EDGE_MASS
  edge%iNode(1) = input%indPlag0 + iPlagIn
  edge%iNode(2) = input%indIntl

  IF (nPeulOxEdges > 0) THEN

    edge => inrt%edges(INRT_BURNING_S_MASS_X0 + nPeulOxEdges)

    edge%tEdge    = INRT_EDGE_MASS_GHO
    edge%iNode(1) = input%indPeul0 + iPeulOx
    edge%iNode(2) = input%indIntl

  END IF ! nPeulOxEdges

  edge => inrt%edges(INRT_BURNING_X_ENER_G + nPeulOxEdges)

  edge%tEdge    = INRT_EDGE_ENER
  edge%iNode(1) = input%indIntl
  edge%iNode(2) = input%indMixt

  edge => inrt%edges(INRT_BURNING_X_ENER_LV + nPeulOxEdges)

  edge%tEdge    = INRT_EDGE_ENER
  edge%iNode(1) = input%indIntl
  edge%iNode(2) = input%indPlagVapor

  edge => inrt%edges(INRT_BURNING_X_MASS_G + nPeulOxEdges)

  edge%tEdge    = INRT_EDGE_MASS
  edge%iNode(1) = input%indIntl
  edge%iNode(2) = input%indMixt

  edge => inrt%edges(INRT_BURNING_X_MASS_L + nPeulOxEdges)

  edge%tEdge    = INRT_EDGE_MASS
  edge%iNode(1) = input%indIntl
  edge%iNode(2) = input%indPlag0 + iPlagOut

  DO iPeulOutEdge = 1,nPeulOutEdges

    edge => inrt%edges(INRT_BURNING_X_MASS_S0  + nPeulOxEdges + iPeulOutEdge)

    edge%tEdge    = INRT_EDGE_MASS
    edge%iNode(1) = input%indIntl
    edge%iNode(2) = input%indPeul0 + iPeulArr(iPeulOutEdge)

  END DO ! iPeulOutEdge

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_DefineBurning

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_DefineBurning.F90,v $
! Revision 1.5  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/02/15 20:40:00  wasistho
! put plag within ifdef
!
! Revision 1.2  2006/02/15 20:17:42  wasistho
! put peul within ifdef
!
! Revision 1.1  2004/12/01 21:56:21  fnajjar
! Initial revision after changing case
!
! Revision 1.9  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.8  2004/03/08 21:57:36  jferry
! better error checking for burning without smoke case
!
! Revision 1.7  2004/03/02 21:47:29  jferry
! Added After Update interactions
!
! Revision 1.6  2003/09/25 15:47:58  jferry
! minor change: error message parameter altered
!
! Revision 1.5  2003/09/19 20:35:26  jferry
! Implemented oxidizer species for burning interaction
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







