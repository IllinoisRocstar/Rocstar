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
! Purpose: update solution pertinent to periodic flows rocperi if applicable
!
! Description: Depending on flow kind, different quantities are updated to
!              maintain balance (e.g. momentum balance), while maintaining
!              other important physical requirements (e.g. mass balance).
!
! Input: region = data of current region
!
! Output: PERI flow solution updated, e.g. cpr mean pgrad for cpr case
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PERI_SolutionUpdate.F90,v 1.5 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_SolutionUpdate( region ) ! PUBLIC

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE PERI_ModInterfaces, ONLY : PERI_CoPgradUpdate
  USE ModError
  USE ModParameters
  USE PERI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global
  INTEGER :: iupdate

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_SolutionUpdate.F90,v $ $Revision: 1.5 $'

  global => region%global
  CALL RegisterFunction( global,'PERI_SolutionUpdate',&
  'PERI_SolutionUpdate.F90' )

! update flow solution --------------------------------------------------------

  iupdate = 0

  IF (global%flowType == FLOW_UNSTEADY .AND. &
      global%solverType == SOLV_EXPLICIT) THEN
    IF (region%periInput%flowKind == PERI_FLOW_CPR .AND. &
        region%irkStep == global%nrkSteps) THEN
      iupdate = 1
    ELSEIF (region%periInput%flowKind == PERI_FLOW_CHANNEL) THEN
      iupdate = 1
    ENDIF
  ELSE
    IF (region%periInput%flowKind == PERI_FLOW_CPR .OR. &
        region%periInput%flowKind == PERI_FLOW_CHANNEL) THEN
      iupdate = 1
    ENDIF
  ENDIF
  
  IF (iupdate == 1) THEN
    CALL PERI_CoPgradUpdate( region )
  ENDIF

! finalize -------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PERI_SolutionUpdate

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_SolutionUpdate.F90,v $
! Revision 1.5  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/06/11 21:49:40  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/06/09 01:21:34  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/06/08 23:56:56  wasistho
! changed nomenclature
!
! Revision 1.6  2004/01/21 03:52:57  wasistho
! modify due to newly unsteady explicit and implicit (dual t-s) option
!
! Revision 1.5  2003/04/25 23:17:12  wasistho
! update cpr pgrad per time-step for unsteady
!
! Revision 1.4  2003/04/05 02:03:51  wasistho
! regions to region in PERI_solutionUpdate
!
! Revision 1.3  2003/04/03 20:56:20  wasistho
! replace regions to region in pgradUpdate
!
! Revision 1.2  2003/04/02 03:03:36  wasistho
! removed unused files
!
! Revision 1.1.1.1  2003/03/29 03:36:30  wasistho
! install ROCPERI
!
!
!******************************************************************************







