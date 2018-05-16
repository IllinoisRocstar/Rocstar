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
! Purpose: define data types related to periodic flows (rocperi).
!
! Description: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModPeriodic.F90,v 1.6 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModPeriodic

  USE ModDataTypes
  IMPLICIT NONE

! input -----------------------------------------------------------------------

  TYPE t_peri_input

! - general PERI input

    INTEGER      :: flowKind        ! 1:cpr, 2:channel 
    INTEGER      :: nVar            ! number of MPI-send/recv variables 
    INTEGER      :: split(3)        ! parallel decomposition directions
    REAL(RFREAL) :: minmax(2)       ! min/max of wall normal coordinat
    REAL(RFREAL) :: meanPgrad       ! mean pressure gradient
    REAL(RFREAL) :: bulkmFlux       ! bulk mass flux 
                                    ! = 2.delta.mRate/cprEpsilon (CPR)
                                    ! = 2.delta.ubulk (CHANNEL)
! - CPR input

    REAL(RFREAL) :: minjRate        ! mass injection rate = rhoinj.vinj
    REAL(RFREAL) :: cprEpsilon      ! ratio between injection/bulk mrate
    REAL(RFREAL) :: headPres        ! head end pressure
    REAL(RFREAL) :: headTemp        ! head end temperature

! - CHANNEL input

    INTEGER      :: pgradType       ! type of pressure gradient calc. method
    REAL(RFREAL) :: cnlRetau        ! Reynolds number based on friction velocity
    REAL(RFREAL) :: cnlCvel         ! mean center velocity by Dean`s relation
    REAL(RFREAL) :: cnlUtau         ! friction velocity

  END TYPE t_peri_input

! data ------------------------------------------------------------------------

  TYPE t_peri

! - CPR data

    REAL(RFREAL), POINTER :: cprVar(:,:), varSend(:,:), varRecv(:,:)
    REAL(RFREAL), POINTER :: cvMean(:,:)

  END TYPE t_peri

END MODULE ModPeriodic

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModPeriodic.F90,v $
! Revision 1.6  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/03/07 05:06:12  wasistho
! install hybrid DESSA turbulence model
!
! Revision 1.3  2004/06/10 22:46:33  wasistho
! removed Rflu periInput
!
! Revision 1.2  2003/09/18 01:56:00  wasistho
! added ijksplit and pgradType in PERI_PgradUpdate
!
! Revision 1.1  2003/03/29 03:29:13  wasistho
! install ROCPERI
!
!
!******************************************************************************






