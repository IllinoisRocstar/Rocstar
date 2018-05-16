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
! Purpose: wrapper for patch update step in PLAG module.
!
! Description: none.
!
! Input: regions = data of all regions,
!        istage  = current RK stage
!
! Output: regions%levels%mixt = new solution after one time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_PatchUpdateWrapper.F90,v 1.3 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_PatchUpdateWrapper( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI
  USE ModParameters

  USE ModInterfacesLagrangian, ONLY : PLAG_NonCvUpdate

  USE PLAG_ModInterfaces, ONLY : PLAG_PatchBufferSendRecv, PLAG_PatchUpdate

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_region), POINTER :: region

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_PatchUpdateWrapper.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PLAG_PatchUpdateWrapper',&
  'PLAG_PatchUpdateWrapper.F90' )

! check if module is active in any region======================================

  IF (.NOT. global%plagUsed) GOTO 999

! loop over stages and regions ================================================

  DO iReg=1,global%nRegions

    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE ) THEN             ! on my processor
      CALL PLAG_PatchUpdate( regions, iReg )
    ENDIF ! regions

  ENDDO ! iReg

  CALL PLAG_PatchBufferSendRecv( regions )

  DO iReg=1,global%nRegions

    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE ) THEN             ! on my processor
      region => regions(iReg)
      CALL PLAG_NonCvUpdate( region )
    ENDIF ! regions

  ENDDO ! iReg

! finalize ====================================================================

999 CONTINUE
  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_PatchUpdateWrapper

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_PatchUpdateWrapper.F90,v $
! Revision 1.3  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:58:00  fnajjar
! Initial revision after changing case
!
! Revision 1.7  2004/11/14 20:30:38  haselbac
! Adapted interface
!
! Revision 1.6  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.5  2004/02/26 21:14:28  haselbac
! Adapted USE <interface> statements and calls
!
! Revision 1.4  2004/02/13 16:23:40  fnajjar
! Include correct Interface call for PLAG_nonCvUpdate
!
! Revision 1.3  2003/06/19 15:47:08  fnajjar
! Fixed comments on Input field
!
! Revision 1.2  2003/05/05 18:11:08  fnajjar
! added check that particles are used in some region
!
! Revision 1.1  2003/03/28 19:53:10  fnajjar
! Initial import of wrapper routines for RocfluidMP
!
!******************************************************************************







