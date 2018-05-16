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
! Purpose: converts the transfers of primary quantities along Edges for each
!          cell into augmentations of RHS terms for all quantities, for an
!          interaction not involving Lagrangian particles
!
! Description: none.
!
! Input: iInrt = index of interaction
!
! Output: augments region%levels(iLev)%...%rhs structures
!
! Notes:
!
!   The RHS structures use opposite sign as the input source structure
!
!******************************************************************************
!
! $Id: INRT_AugmentConSources.F90,v 1.3 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_AugmentConSources( region,iInrt )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE INRT_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region
  INTEGER,        INTENT(IN)    :: iInrt

! ... local variables
  CHARACTER(CHRLEN)  :: RCSIdentString
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_AugmentConSources.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_AugmentConSources',&
  'INRT_AugmentConSources.F90' )

! begin -----------------------------------------------------------------------

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_AugmentConSources

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_AugmentConSources.F90,v $
! Revision 1.3  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:56:09  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2004/07/28 15:42:12  jferry
! deleted defunct constructs: useDetangle, useSmokeDrag, useSmokeHeatTransfer
!
! Revision 1.4  2004/03/02 21:48:09  jferry
! First phase of replacing Detangle interaction
!
! Revision 1.3  2003/09/25 15:48:43  jferry
! implemented Boiling Regulation interaction
!
! Revision 1.2  2003/03/11 16:09:39  jferry
! Added comments
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************







