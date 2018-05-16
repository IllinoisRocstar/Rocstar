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
  MODULE M_Rocburn_APN_Global_Data 

  IMPLICIT NONE

! INTEGER, PARAMETER :: idp=SELECTED_REAL_KIND(15,100)
!
! ----------------------------------------------------------------------
!
! delcare global variables for Rocburn_1D_APN
!
!
!  ===================
!  D A T A   T Y P E S
!  ===================
   INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(P=14,R=30)
   INTEGER, PARAMETER :: MATMAX = 4

!
!  ----------------------------------------
!  ROCBURN_ZN global data
!  ----------------------------------------


   TYPE, PUBLIC :: G_BURN_1D

     !
     !   MPI related
     !

     INTEGER  ::  MPI_COMM_ROCBURN, rank

     !
     !   variables for burning rate
     !

     REAL(DBL)       :: a_p(MATMAX), n_p(MATMAX)
     REAL(DBL)       :: Tf_adiabatic(MATMAX), To             ! To has no use, only place holder for Rocburn_2D
     REAL(DBL)       :: xmax(MATMAX)

     !
     !   variables for grid generation (no use, only place holder for Rocburn_2D)
     !

     INTEGER :: nxmax
     INTEGER :: nmat
     INTEGER :: verbosity

   END TYPE G_BURN_1D
!
! ----------------------------------------------------------------------

      
  END MODULE M_Rocburn_APN_Global_Data






