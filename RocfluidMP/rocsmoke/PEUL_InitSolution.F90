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
! Purpose: initialize memory for all variables associated with the Eulerian
!          particles (PEUL) for all active regions on current processor.
!
! Description: none.
!
! Input: iReg = current region number
!        region = current region
!
! Output: region%peul = peul variables
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_InitSolution.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_InitSolution( iReg, region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartEul,    ONLY : t_peul
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  INTEGER,        INTENT(IN)    :: iReg
  TYPE(t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: iCont, iLev

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_peul),   POINTER :: peul
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_InitSolution.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global, 'PEUL_InitSolution',&
  'PEUL_InitSolution.F90' )

! loop over all grid levels

  DO iLev=1,region%nGridLevels

    peul => region%levels(iLev)%peul

! nothing to do yet

  ENDDO ! iLev

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_InitSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_InitSolution.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:37  haselbac
! Initial revision after changing case
!
! Revision 1.2  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.1  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







