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
  MODULE M_ROCBURN_1D_APN
! ----------------------------------------------------------------
!
! This module provides the following subroutines for ROCBURN_1D_APN  
! to support ROCBURN_2D:
!
!     INITIALIZE_0D
!     INITIALIZE_1D
!     GET_BURNING_RATE_1D
!
!
! The combustion model used in this quasi-steady empirical model
! rb=a*P^n
!
!   Creation Date    :  Sep. 10, 2002
!
!   Modifications    :
!
!    No.     Date         Authors       Description
!
! ----------------------------------------------------------------
!
!

  USE M_Rocburn_APN_Global_Data

  IMPLICIT NONE

  INCLUDE 'mpif.h'


  CONTAINS

!
!     ------------------------------------------------------------------------
!                             INTERNAL PROCEDURES
!     ------------------------------------------------------------------------

    SUBROUTINE CHECK_ALLOC( ierr)
      INTEGER, INTENT(IN) :: ierr
      IF(ierr /= 0) PRINT *, "ROCBRUN_APN ERROR: unable to allocate memory"
      CALL MPI_ABORT( MPI_COMM_WORLD, -1)
    END SUBROUTINE CHECK_ALLOC
! ============================================================================






    SUBROUTINE INITIALIZE_0D( G_APN, comm, Indir, nxmax, To_read)
      TYPE (G_BURN_1D), POINTER :: G_APN
      INTEGER, INTENT(IN)       :: comm
      CHARACTER(*), INTENT(IN)  :: Indir
      INTEGER, INTENT(OUT)      :: nxmax
      REAL(DBL), INTENT (OUT)   :: To_read

! ---------------------------------------------------------------------------
!     local variables
!


      INTEGER :: ierror

      INTERFACE 
        SUBROUTINE APN_input_0d(G_APN, Indir)
          USE M_Rocburn_APN_Global_Data
          TYPE(G_BURN_1D), POINTER :: G_APN
          CHARACTER(*), INTENT(IN) :: Indir
        END SUBROUTINE APN_input_0d

      END INTERFACE

      ALLOCATE(G_APN)

!
!     read Rocburn_APN 0d data
!

      G_APN%MPI_COMM_ROCBURN  = comm
      CALL MPI_COMM_RANK(comm, G_APN%rank, ierror)

      CALL APN_input_0d(G_APN, Indir)

      nxmax   = G_APN%nxmax
      To_read = G_APN%To

      RETURN

    END SUBROUTINE INITIALIZE_0D
! ============================================================================






    SUBROUTINE INITIALIZE_1D( G_APN, bflag, P_mks, To, rhoc_mks, p_coor, rb_mks, Toa, fr, Tn, Tflame)


      IMPLICIT NONE

      TYPE (G_BURN_1D), POINTER :: G_APN
      INTEGER, INTENT(INOUT)     :: bflag
      REAL(DBL), INTENT (IN)  :: P_mks, To, rhoc_mks, p_coor(3)
      REAL(DBL), INTENT (OUT) :: rb_mks, Toa, fr
      REAL(DBL), INTENT (OUT) :: Tn(:)
      REAL(DBL), INTENT (OUT) :: Tflame

! ---------------------------------------------------------------------------
!     local variables

      INTEGER :: mat

! ---------------------------------------------------------------------------

      IF( bflag == 0 ) THEN
!
!        not burning or ignition simulation required
!        set condensed phase temperature profile to initial temperature
!

         rb_mks = 0.
         Tn = G_APN%To
         Tflame = G_APN%To

      ELSE
!
!        burning from onset and no ignition simulation required
!        
!        P    = P_mks *9.869232667E-6           ! Pa to atm

!RAF
!RAF Find out which material applies for this face.  I need the x coord.
!RAF
 
         DO mat=1,G_APN%nmat
           IF (p_coor(1) <= G_APN%xmax(mat)) EXIT
         END DO
         mat = MIN(mat, G_APN%nmat)

         rb_mks = (G_APN%a_p(mat)*(P_mks*9.869232667E-6)**G_APN%n_p(mat))*0.01   ! cm/s to m/s
         Tflame = G_APN%Tf_adiabatic(mat)

         Tn = 700.0                             ! do not acutually need temperature profile
                                                ! and do not use in the code

      END IF

      RETURN

    END SUBROUTINE INITIALIZE_1D
!! ============================================================================






    SUBROUTINE GET_BURNING_RATE_1D( G_APN, delt, P_mks, To, Tn,   &
              qc_mks, qc_old_mks, qr_mks, qr_old_mks, rhoc_mks, &
              Toa, rb_mks, fr, bflag, Tnp1, Tflame, p_coor)
      TYPE (G_BURN_1D), POINTER :: G_APN
      REAL(DBL), INTENT (IN)      :: delt, P_mks, To
      REAL(DBL), INTENT (IN)      :: Tn(:)
      REAL(DBL), INTENT (IN)      :: qc_mks, qc_old_mks, qr_mks, qr_old_mks
      REAL(DBL), INTENT (IN)      :: rhoc_mks
      REAL(DBL), INTENT (INOUT)   :: Toa, rb_mks, fr
      INTEGER,   INTENT (INOUT)   :: bflag
      REAL(DBL), INTENT (OUT)     :: Tnp1(:)
      REAL(DBL), INTENT (OUT)     :: Tflame
      REAL(DBL), INTENT (IN)      :: p_coor(3)


! ---------------------------------------------------------------------------
!     local variables

      INTEGER :: mat

! ---------------------------------------------------------------------------



      IF(bflag == 1) THEN

!
!        propellant burning already, calculate burning rate using 
!
!        P      = P_mks *9.869232667E-6           ! Pa to atm

!RAF
!RAF Find out which material applies for this face.  I need the x coord.
!RAF
 
         DO mat=1,G_APN%nmat
           IF (p_coor(1) <= G_APN%xmax(mat)) EXIT
         END DO
         mat = MIN(mat, G_APN%nmat)

!        rb_mks = (G_APN%a_p(mat)*(P_mks*9.869232667E-6)**G_APN%n_p(mat))*0.01   ! cm/s to m/s
!KJM
!KJM Change the rb_mks calculation to be normalized to 1000psi or 68.046 atm
!KJM
         rb_mks = (G_APN%a_p(mat)*(P_mks*9.869232667E-6/68.046)**G_APN%n_p(mat))*0.01
         Tflame = G_APN%Tf_adiabatic(mat)

      ELSE 
!
!        propellant not burning yet
!
  
         rb_mks = 0.0
         Tflame = G_APN%To

      END IF

      RETURN

    END SUBROUTINE GET_BURNING_RATE_1D

    SUBROUTINE FINALIZE_0D( G_APN)
      TYPE (G_BURN_1D), POINTER :: G_APN

      DEALLOCATE( G_APN)
    END SUBROUTINE FINALIZE_0D


  END MODULE M_ROCBURN_1D_APN






