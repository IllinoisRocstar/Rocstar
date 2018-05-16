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
  MODULE M_Rocburn_ZN_Global_Data 

  IMPLICIT NONE

! INTEGER, PARAMETER :: idp=SELECTED_REAL_KIND(15,100)
!
! ----------------------------------------------------------------------
!
! delcare global variables for Rocburn_1D_ZN
! 
! note that all the global vairables need to be included in
! the derived type G_BUNR_1D
!
!
!  ===================
!  D A T A   T Y P E S
!  ===================
   INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(P=14,R=30)

!
!  ----------------------------------------
!  ROCBURN_ZN global data
!  ----------------------------------------


   TYPE, PUBLIC :: G_BURN_1D

!
!    MPI related
!

     INTEGER  ::  MPI_COMM_ROCBURN, rank

!
!    flags
!

     INTEGER  ::  Model_combustion

!
!    propellant thermophysical properties
!
     REAL(DBL)       :: Ac, Ec, Qc, alfac, rhoc, C, lamc
     REAL(DBL)       :: Bg,     Qg, lamg, MW, R
     REAL(DBL)       :: Ka
     REAL(DBL)       :: Tf_adiabatic, To

!
!    control variables
!

     REAL(DBL)               :: delt_max, xmax, beta, tol_Ts
     INTEGER               :: igrid, itermax

!
!    variables for grid generation
!
     INTEGER :: nxmax, nx
     REAL(DBL), POINTER ::  x(:), z(:), zx(:), zxx(:)
     REAL(DBL) :: delz

!
!     varialbles for Zeldovich-Novozhilov (ZN) approach
!

      REAL(DBL)       :: a_p, n_p
      REAL(DBL)       :: a_T, n_T

!
!     varialbles for ignition modeling
!

      INTEGER         :: ign_flag
      REAL(DBL)       :: To_cold

!
! ----------------------------------------------------------------------

    END TYPE G_BURN_1D
      
  END MODULE M_Rocburn_ZN_Global_Data






