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
! Purpose: define parameters for Eulerian particles
!
! Description: none
!
! Notes: none
!
!******************************************************************************
!
! $Id: PEUL_ModParameters.F90,v 1.6 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

MODULE PEUL_ModParameters

  IMPLICIT NONE

! Eulerian particles: PEUL ----------------------------------------------------

  INTEGER, PARAMETER :: CV_PEUL_DENS = 1, &    ! index of density (rho)
                        CV_PEUL_NEQS = 1       ! total no. of equations

! boundary conditions ---------------------------------------------------------

  INTEGER, PARAMETER :: BCDAT_PEUL_INFLOW_DENS = 1, & ! index of density in BC
                        BCDAT_PEUL_FARF_DENS   = 1, & ! index of density in BC
                        BCDAT_PEUL_INJECT_FRAC = 1    !     of mass frac in BC

! Methods ---------------------------------------------------------------------

  INTEGER, PARAMETER :: PEUL_METHV_FLUIDVEL  = 0, & ! set smoke = fluid vel
                        PEUL_METHV_EQEUL     = 1    ! Eq Eul method

! Values of switches ----------------------------------------------------------

  INTEGER, PARAMETER :: PEUL_NEG_REPORT_NONE = 0, &
                        PEUL_NEG_REPORT_USED = 1

  INTEGER, PARAMETER :: PEUL_CLIP_MODEL_NONE = 0, &
                        PEUL_CLIP_MODEL_USED = 1

END MODULE PEUL_ModParameters

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_ModParameters.F90,v $
! Revision 1.6  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2004/05/03 15:09:42  jferry
! added equilibrium Eulerian capability for smoke
!
! Revision 1.3  2004/03/05 21:11:45  wasistho
! removed PEUL mpi-tag-shift
!
! Revision 1.2  2004/03/02 21:42:47  jferry
! Added clipping options and corner and edge cell updates
!
! Revision 1.1  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
!******************************************************************************






