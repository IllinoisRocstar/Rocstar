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
! Purpose: Compute maximum number of edges in any particle-related and
!          cell-related interaction
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels(iLev)%inrt%inrtInput%maxDisEdges and maxConEdges
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_ComputeMaxEdges.F90,v 1.3 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_ComputeMaxEdges( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModInteract,   ONLY : t_inrt_interact
  USE INRT_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: iInrt

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: maxConEdges,maxDisEdges

  TYPE(t_inrt_interact), POINTER :: inrt
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: INRT_ComputeMaxEdges.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_ComputeMaxEdges',&
  'INRT_ComputeMaxEdges.F90' )

! begin -----------------------------------------------------------------------

  maxConEdges = 0
  maxDisEdges = 0

! compute maximum number of edges

  DO iInrt = 1,INRT_TYPE_TOTAL

    inrt => region%inrtInput%inrts(iInrt)

    IF (inrt%used) THEN
      IF (inrt%pclsUsed) THEN
        maxDisEdges = MAX(maxDisEdges,inrt%nEdges)
      ELSE
        maxConEdges = MAX(maxConEdges,inrt%nEdges)
      ENDIF
    END IF ! inrt%used

  END DO ! iInrt

  region%inrtInput%maxConEdges = maxConEdges
  region%inrtInput%maxDisEdges = maxDisEdges

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_ComputeMaxEdges

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_ComputeMaxEdges.F90,v $
! Revision 1.3  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:56:19  fnajjar
! Initial revision after changing case
!
! Revision 1.1  2004/07/26 17:06:54  fnajjar
! initial import
!
!******************************************************************************







