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
! Purpose: define various integer parameters for RVAV routines
!
! Description: none
!
! Notes: none
!
!******************************************************************************
!
! $Id: RVAV_ModParameters.F90,v 1.7 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE RVAV_ModParameters

  IMPLICIT NONE

! file types ------------------------------------------------------------------

  INTEGER, PARAMETER :: IF_RVAV_INPUT     = 10    ! this is the fileID
  INTEGER, PARAMETER :: FILE_COMPUTED     = 10, &
                        FILE_EXPERIMENTAL = 20, &
                        FILE_ANALYTICAL   = 30
  INTEGER, PARAMETER :: IF_RVAV_FILE_S2   = 100   ! fileID for Stream 2
  INTEGER, PARAMETER :: IF_RVAV_FILE_BL   = 110   ! fileID for Blasius solution

! output ----------------------------------------------------------------------

  INTEGER, PARAMETER :: COMPUTE_ERRORS_ONLY = 10, & 
                        PLOT_ERRORS_ONLY    = 20, &
                        COMPUTE_AND_PLOT    = 30

! numerical equivalents for the conserved and extracted variables in RVAV -----
  INTEGER, PARAMETER :: RVAV_RHO     =   1, & ! all conserved variables have
                        RVAV_RHOU    =   2, & ! single digit identifiers
                        RVAV_RHOV    =   3, &
                        RVAV_RHOW    =   4, &
                        RVAV_RHOET   =   5, &
! ... the primitive (extracted) variables
                        RVAV_U       =  11, & ! all primitive variable with the
                                              ! exception of density have two 
                                              ! digit identifiers
                        RVAV_V       =  12, &
                        RVAV_W       =  13, &
                        RVAV_P       =  14, &
                        RVAV_T       =  15, &
! ... all derived (extracted) variables
                        RVAV_TSTAG   =  101, & ! all derived variables have
                        RVAV_PSTAG   =  102, & ! three digit identifiers
                        RVAV_S       =  201, &
                        RVAV_OMEGAX  =  301, &
                        RVAV_OMEGAY  =  302, &
                        RVAV_OMEGAZ  =  303

! analytical solution type for similarity analysis ----------------------------

  INTEGER, PARAMETER :: RVAV_CULICK   =   131, & 
                        RVAV_BLASIUS  =   132, & 
                        RVAV_GAMMBUMP =   133

END MODULE RVAV_ModParameters

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_ModParameters.F90,v $
! Revision 1.7  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:18:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2002/07/25 00:51:47  jblazek
! Option for TVD type pressure switch.
!
! Revision 1.4  2002/06/24 15:47:35  f-najjar
! Included File Parameter for Blasius
!
! Revision 1.3  2002/06/19 20:26:22  f-najjar
! Included GAMM Bump Definition
!
! Revision 1.2  2002/06/15 17:53:19  f-najjar
! Parameters for Similarity Solution
!
! Revision 1.1.1.1  2002/06/03 21:41:29  f-najjar
! Initial Import of RocVaV
!
!******************************************************************************






