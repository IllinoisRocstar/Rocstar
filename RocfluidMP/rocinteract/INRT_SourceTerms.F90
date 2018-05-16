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
! Purpose: compute the source terms for all interactions.
!
! Description:
!   The interactions are computed in the order specified by the integer
!   inrt%order for each interaction inrt.  An interaction is computed by
!   calling a routine specific to the interaction (INRT_Calc*) which computes
!   the primary quantities transferred, followed by one of two generic routines
!   (INRT_Augment*Sources) which then computes secondary quantities and augment
!   source terms as the Tokens on the interaction Edges allow.
!
! Input: region
!
! Output: region with source terms augmented by all the interactions
!
! Notes: the outer loop structure is designed so that arbitrary values of
!   inrt%order may be specified:  this is less restrictive than looping over
!   a given range of allowed values of inrt%order.
!
!   The rather strange value of "big" is used to avoid the possibility that
!   HUGE(1) = +2^31, and thus -HUGE(1) = -2^31 = +2^31 > 0.  In practice,
!   however, HUGE(1) = 2^31-1, so -HUGE(1) is indeed negative, but it seems
!   improper to rely on this.
!
!******************************************************************************
!
! $Id: INRT_SourceTerms.F90,v 1.3 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_SourceTerms( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInteract,   ONLY : t_inrt_interact
  USE ModError
  USE INRT_ModParameters

  USE INRT_ModInterfaces, ONLY : INRT_AugmentDisSources, &
    INRT_AugmentConSources, INRT_CalcDrag, INRT_CalcHeatTransferNonBurn, &
    INRT_CalcBurning, INRT_CalcScouring
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: order, iInrt

! ... local variables
  CHARACTER(CHRLEN)  :: RCSIdentString

  INTEGER :: big,currentOrder,nextOrder

  TYPE(t_inrt_interact), POINTER :: inrt
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_SourceTerms.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_SourceTerms',&
  'INRT_SourceTerms.F90' )

! begin -----------------------------------------------------------------------

  big = HUGE(1) - 3 ! see Notes above for explanation

  currentOrder = -big

  DO
    nextOrder = big

    DO iInrt = 1,INRT_TYPE_TOTAL

      inrt => region%inrtInput%inrts(iInrt)

      IF (inrt%used) THEN                    ! interaction used?

        IF (inrt%order == currentOrder) THEN ! time to do it?

! ------- compute primary source terms for interaction

          SELECT CASE (iInrt)

          CASE (INRT_TYPE_BURNING)
            CALL INRT_CalcBurning(region)

          CASE (INRT_TYPE_DRAG)
            CALL INRT_CalcDrag(region)

          CASE (INRT_TYPE_HTRANSNB)
            CALL INRT_CalcHeatTransferNonBurn(region)

          CASE (INRT_TYPE_SCOURING)
            CALL INRT_CalcScouring(region)

          CASE (INRT_TYPE_BOILRGN)
            CONTINUE ! nothing to calculate at this stage

          CASE DEFAULT
            CALL ErrorStop( global,ERR_INRT_NOINRT,__LINE__,inrt%name )

          END SELECT ! iInrt

! ------- update RHS for primary and secondary quantities

          IF (inrt%pclsUsed) THEN
            CALL INRT_AugmentDisSources(region,iInrt)
          ELSE
            CALL INRT_AugmentConSources(region,iInrt)
          END IF ! inrt%pclsUsed

        ELSE IF (inrt%order > currentOrder) THEN ! future interaction?

          nextOrder = MIN(nextOrder,inrt%order)  ! if so, see if it is next

        END IF ! inrt%order

      END IF ! inrt%used

    END DO ! iInrt

    IF (nextOrder == big) EXIT ! exit loop if no future interactions

    currentOrder = nextOrder   ! set currentOrder to index next interaction

  END DO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_SourceTerms

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_SourceTerms.F90,v $
! Revision 1.3  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:56:46  fnajjar
! Initial revision after changing case
!
! Revision 1.8  2004/07/28 15:42:13  jferry
! deleted defunct constructs: useDetangle, useSmokeDrag, useSmokeHeatTransfer
!
! Revision 1.7  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2003/09/25 15:48:43  jferry
! implemented Boiling Regulation interaction
!
! Revision 1.5  2003/04/03 21:10:18  jferry
! implemented additional safety checks for rocinteract
!
! Revision 1.4  2003/04/03 16:18:28  fnajjar
! Include routines for burning and scouring
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







