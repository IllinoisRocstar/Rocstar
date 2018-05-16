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
! Purpose: initialize interaction data structure and get user input
!
! Description: none
!
! Input: regions
!
! Output: regions with interactions initialized
!
! Notes: regions must already have data about the multiphysics modules,
!        e.g., the number of smoke types and Lagrangian particle constituents
!
!******************************************************************************
!
! $Id: INRT_UserInput.F90,v 1.3 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_UserInput( regions ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE INRT_ModParameters

  USE INRT_ModInterfaces, ONLY : INRT_Initialize,INRT_ReadInputFile, &
                                 INRT_ComputeMaxEdges,INRT_CheckUserInput
  IMPLICIT NONE

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_UserInput.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'INRT_UserInput',&
  'INRT_UserInput.F90' )

! begin -----------------------------------------------------------------------

  global%inrtUsed = (INRT_TYPE_TOTAL > 0) .AND. &
                    (global%peulUsed .OR. global%specUsed .OR. global%plagUsed)

  IF (global%inrtUsed) THEN

! - initialize interactions

    DO iReg=LBOUND(regions,1),UBOUND(regions,1)
      CALL INRT_Initialize( regions(iReg) )
    END DO ! iReg

! - read user input

    CALL INRT_ReadInputFile( regions ) ! can modify value of inrtUsed

  END IF ! inrtUsed

  IF (global%inrtUsed) THEN

! - Compute maximum number of discrete- and continuum-based Edges

    DO iReg=LBOUND(regions,1),UBOUND(regions,1)
      CALL INRT_ComputeMaxEdges( regions(iReg) )
    END DO ! iReg

! - check user input

    DO iReg=LBOUND(regions,1),UBOUND(regions,1)
      CALL INRT_CheckUserInput( regions(iReg) )
    END DO ! iReg

  END IF ! inrtUsed

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_UserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_UserInput.F90,v $
! Revision 1.3  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:56:48  fnajjar
! Initial revision after changing case
!
! Revision 1.9  2004/07/28 15:29:20  jferry
! created global variable for spec use
!
! Revision 1.8  2004/07/26 17:05:51  fnajjar
! moved allocation of inrtSources into Rocpart
!
! Revision 1.7  2004/07/23 22:43:17  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.6  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.5  2004/03/02 21:49:23  jferry
! Added inrtUsed flag to mixture data structure
!
! Revision 1.4  2003/09/26 21:46:54  fnajjar
! Modified ModInterfaces call to ModInterfacesInteract
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







