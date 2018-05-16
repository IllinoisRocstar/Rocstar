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
  MODULE M_ROCBURN_1D_ZN
! ----------------------------------------------------------------
!
! This module provides the following subroutines for ROCBURN_1D_ZN  
! to support ROCBURN_2D:
!
!     INITIALIZE_0D
!     INITIALIZE_1D
!     GET_FILM_COEFF_1D
!     GET_TIME_STEP_1D
!     GET_BURNING_RATE_1D
!
!
! The combustion model used in this module is Zeldovich-Novozhilov 
! model implemented by J. Weber, K.C. Tang, and M. Q. Brewster.
!
!   Creation Date    :  Sep. 2, 2002
!
!   Modifications    :
!
!    No.     Date         Authors       Description
!
! ----------------------------------------------------------------
!
!

  USE M_Rocburn_ZN_Global_Data

  IMPLICIT NONE

  INCLUDE 'mpif.h'


! ----------------------------------------------------------------
! local variables

  CONTAINS

!
!     ------------------------------------------------------------------------
!                             INTERNAL PROCEDURES
!     ------------------------------------------------------------------------

  SUBROUTINE CHECK_ALLOC( ierr)
    INTEGER, INTENT(IN) :: ierr
    IF(ierr /= 0) THEN
       PRINT *, "ROCBRUN_ZN ERROR: unable to allocate memory"
       CALL MPI_ABORT( MPI_COMM_WORLD, -1)
    ENDIF
  END SUBROUTINE CHECK_ALLOC

  SUBROUTINE INITIALIZE_0D(G_ZN, comm, Indir, nxmax, To_read)
      TYPE(G_BURN_1D), POINTER :: G_ZN
      INTEGER, INTENT(IN)      :: comm
      CHARACTER(*), INTENT(IN) :: Indir
      INTEGER, INTENT(OUT)     :: nxmax
      REAL(DBL), INTENT (OUT)  :: To_read

! ---------------------------------------------------------------------------
!     local variables
!


      INTEGER :: ierror

      INTERFACE 
        SUBROUTINE ZN_input_0d(G_ZN, Indir)
          USE M_Rocburn_ZN_Global_Data
          TYPE(G_BURN_1D), POINTER :: G_ZN
          CHARACTER(*), INTENT(IN)  :: Indir
        END SUBROUTINE ZN_input_0d

        SUBROUTINE ZN_gen_grid(G_ZN, gridtype, numx, x, z, zx, zxx)
          USE M_Rocburn_ZN_Global_Data
          TYPE(G_BURN_1D), POINTER :: G_ZN
          INTEGER, INTENT(IN) :: gridtype, numx
          REAL(DBL), INTENT (OUT) :: x(:), z(:), zx(:), zxx(:)
        END SUBROUTINE ZN_gen_grid
      END INTERFACE

      ALLOCATE(G_ZN)

!
!     read Rocburn_ZN 0d data
!

      G_ZN%MPI_COMM_ROCBURN     = comm
      CALL MPI_COMM_RANK(comm, G_ZN%rank, ierror)

      CALL ZN_input_0d(G_ZN, TRIM(Indir))

      nxmax   = G_ZN%nxmax
      To_read = G_ZN%To

!
!     allocate memory for variables related to grid generation
!
      ALLOCATE(G_ZN%x(1:G_ZN%nxmax), STAT=ierror); CALL CHECK_ALLOC (ierror)
      ALLOCATE(G_ZN%z(1:G_ZN%nxmax), STAT=ierror); CALL CHECK_ALLOC (ierror)
      ALLOCATE(G_ZN%zx(1:G_ZN%nxmax), STAT=ierror); CALL CHECK_ALLOC (ierror)
      ALLOCATE(G_ZN%zxx(1:G_ZN%nxmax), STAT=ierror); CALL CHECK_ALLOC (ierror)

!
!     grid generation
!
      CALL ZN_gen_grid(G_ZN, G_ZN%igrid, G_ZN%nx,                  &
                       G_ZN%x, G_ZN%z, G_ZN%zx, G_ZN%zxx)

      RETURN

    END SUBROUTINE INITIALIZE_0D
! ============================================================================






    SUBROUTINE INITIALIZE_1D( G_ZN, bflag, P_mks, To, rhoc_mks, p_coor, rb_mks, Toa, fr, Tn, Tflame)


      IMPLICIT NONE

      TYPE(G_BURN_1D), POINTER :: G_ZN
      INTEGER, INTENT(INOUT)  :: bflag
      REAL(DBL), INTENT (IN)  :: P_mks, To, rhoc_mks, p_coor(3)
      REAL(DBL), INTENT (OUT) :: rb_mks, Toa, fr
      REAL(DBL), INTENT (OUT) :: Tn(:)
      REAL(DBL), INTENT (OUT) :: Tflame
     
      INTERFACE
        SUBROUTINE ZN_ssWSB(G_ZN, P, qr, To, rhoc, rb, Ts, Tf, fr, Tn)
          USE M_Rocburn_ZN_Global_Data
          TYPE(G_BURN_1D), POINTER :: G_ZN
          REAL(DBL), INTENT(IN) :: P, qr, To, rhoc
          REAL(DBL), INTENT(OUT) :: rb, Ts, Tf, fr
          REAL(DBL), INTENT(OUT) :: Tn(:)
        END SUBROUTINE ZN_ssWSB
      END INTERFACE

! ---------------------------------------------------------------------------
!     local variables
      REAL(DBL) :: P, rhoc, rb, Ts, Tf
      REAL(DBL) :: qr

! ---------------------------------------------------------------------------

      IF( (bflag == 0).OR.(G_ZN%ign_flag == 1) ) THEN
!
!        not burning or ignition simulation required
!        set condensed phase temperature profile to initial temperature
!
         Tn = G_ZN%To

      ELSE
!
!        burning from onset and no ignition simulation required, 
!        calculate steady state solution and use as initial condition
!        
         P    = P_mks *9.869232667E-6           ! Pa to atm
         rhoc = rhoc_mks*1.0E-3                 ! Kg/m^3 to g/cm^3
         qr   = 0.0

         CALL ZN_ssWSB(G_ZN, P, qr  , G_ZN%To, rhoc, rb, Ts, Tf, fr, Tn)

         rb_mks = rb*0.01                       ! cm/s to m/s

      END IF

      Tflame = G_ZN%Tf_adiabatic
      Toa    = G_ZN%To

      RETURN

    END SUBROUTINE INITIALIZE_1D
!! ============================================================================






    SUBROUTINE GET_FILM_COEFF_1D(G_ZN, p_coor, Ts, T_euler, P, qc, qcPrime)
       TYPE (G_BURN_1D), POINTER :: G_ZN

       REAL(DBL), INTENT (IN)      :: p_coor(3), Ts, T_euler, P
       REAL(DBL), INTENT (OUT)     :: qc, qcPrime  ! do not change qcPrime
!
!      place holder for calculating convective heat flux qc for
!      Rocburn_1D_ZN
!
!      currently not available; set qc to 0
!
       qc = 0.0
       qcPrime = 0.0

       RETURN

    END SUBROUTINE GET_FILM_COEFF_1D
! ============================================================================






    SUBROUTINE GET_TIME_STEP_1D(G_ZN, rb, Toa, dt_max)
      TYPE (G_BURN_1D), POINTER :: G_ZN

       REAL(DBL), INTENT (IN)      :: rb, Toa
       REAL(DBL), INTENT (OUT)     :: dt_max

! ---------------------------------------------------------------------------
!     local variables

      REAL(DBL) :: dt_c
! ---------------------------------------------------------------------------

      IF( ABS(Toa - G_ZN%To) >= 0.9*G_ZN%To) THEN
         dt_c  = 0.1*G_ZN%alfac/(rb*rb*1.0e4)  ! G_ZN%alfac in cm^2/sec, rb in m/s
         dt_max= MIN(G_ZN%delt_max, dt_c)      ! dt_max in sec
      ELSE
         dt_max= G_ZN%delt_max                 ! dt_max in sec
      END IF

      RETURN

    END SUBROUTINE GET_TIME_STEP_1D
! ============================================================================







    SUBROUTINE GET_BURNING_RATE_1D( G_ZN, delt, P_mks, To, Tn,   &
              qc_mks, qc_old_mks, qr_mks, qr_old_mks, rhoc_mks, &
              Toa, rb_mks, fr, bflag, Tnp1, Tflame)
      TYPE (G_BURN_1D), POINTER :: G_ZN

      REAL(DBL), INTENT (IN)      :: delt, P_mks, To
      REAL(DBL), INTENT (IN)      :: Tn(:)
      REAL(DBL), INTENT (IN)      :: qc_mks, qc_old_mks, qr_mks, qr_old_mks
      REAL(DBL), INTENT (IN)      :: rhoc_mks
      REAL(DBL), INTENT (INOUT)   :: Toa, rb_mks, fr
      INTEGER,   INTENT (INOUT)   :: bflag
      REAL(DBL), INTENT (OUT)     :: Tnp1(:)
      REAL(DBL), INTENT (OUT)     :: Tflame

      INTERFACE
        SUBROUTINE ZN_ssWSB(G_ZN, P, qr, To, rhoc, rb, Ts, Tf, fr, Tn)
          USE M_Rocburn_ZN_Global_Data
          TYPE(G_BURN_1D), POINTER :: G_ZN
          REAL(DBL), INTENT(IN) :: P, qr, To, rhoc
          REAL(DBL), INTENT(OUT) :: rb, Ts, Tf, fr
          REAL(DBL), INTENT(OUT) :: Tn(:)
        END SUBROUTINE ZN_ssWSB

        SUBROUTINE ZN_calc_burning_rate(G_ZN, delt, P, qr, To, rhoc, qr_old, fr_old,  &
                                        Toa, rb, Ts, fr, Tn, Tnp1)
          USE M_Rocburn_ZN_Global_Data
          TYPE(G_BURN_1D), POINTER :: G_ZN
          REAL(DBL), INTENT(IN)  :: delt, P, qr, To, rhoc
          REAL(DBL), INTENT(IN)  :: qr_old, fr_old, Toa
          REAL(DBL), INTENT(OUT) :: rb, Ts, fr
          REAL(DBL), INTENT(IN)  :: Tn(:)
          REAL(DBL), INTENT(OUT) :: Tnp1(:)
        END SUBROUTINE ZN_calc_burning_rate
      END INTERFACE


! ---------------------------------------------------------------------------
!     local variables
      REAL(DBL) :: P, qc, qc_old, qr, qr_old, rhoc, rb
      REAL(DBL) :: Ts, Tf
! ---------------------------------------------------------------------------


      Tflame = G_ZN%Tf_adiabatic

      IF(bflag/=0) THEN

        IF(Tn(1) > G_Zn%To) THEN
!
!         propellant burning already, calculate burning rate using ZN_cal_burning_rate
!
          P      = P_mks *9.869232667E-6           ! Pa to atm
          rb     = rb_mks * 100.0                  ! m/s to cm/s
          qr     = qr_mks * 0.2388459E-4           ! W/m^2 to cal/cm^2/s  1 J = 0.2388459 cal
          qc     = qc_mks *  0.2388459E-4          ! W/m^2 to cal/cm^2/s
          qr_old = qr_old_mks * 0.2388459E-4       ! W/m^2 to cal/cm^2/s
          qc_old = qc_old_mks *  0.2388459E-4      ! W/m^2 to cal/cm^2/s
          rhoc   = rhoc_mks*1.0E-3                 ! Kg/m^3 to g/cm^3
          Ts     = Tn(1)

          CALL ZN_calc_burning_rate(G_ZN, delt, P, qr, To, rhoc, qr_old, fr,     &
                               Toa, rb, Ts, fr, Tn, Tnp1)


          rb_mks = rb*0.01                         ! cm/s to m/s

        ELSE 
!
!         propellant burning for the first time , set initial condition using ZN_ssWSB
!
          P    = P_mks *9.869232667E-6             ! Pa to atm
          rhoc = rhoc_mks*1.0E-3                   ! Kg/m^3 to g/cm^3
          qr     = qr_mks * 0.2388459E-4           ! W/m^2 to cal/cm^2/s  1 J = 0.2388459 cal

          CALL ZN_ssWSB(G_ZN, P, qr, G_ZN%To, rhoc, rb, Ts, Tf, fr, Tnp1)

          rb_mks = rb*0.01                         ! cm/s to m/s
        END IF


      ELSE
!
!       propellant not burning yet, check for ignition simulation requirement
!
        IF(G_ZN%ign_flag == 1) THEN
!
!         ignition simulation required
!         
!         place holder for ignition model
!
          WRITE(*,*) 'ROCBURN_ZN: rank=',G_ZN%rank
          WRITE(*,*) '  Error: igntion model not ready'
          WRITE(*,*) '  job aborted'
          CALL MPI_ABORT( MPI_COMM_WORLD, -1)
          STOP

        ELSE
!
!         ignition simulation not required
!
          bflag = 1
          RETURN

        END IF

      END IF

      RETURN

    END SUBROUTINE GET_BURNING_RATE_1D

!***************************************************************************
    SUBROUTINE FINALIZE_0D(G_ZN)
 
      TYPE (G_BURN_1D), POINTER :: G_ZN

      DEALLOCATE( G_ZN)

    END SUBROUTINE FINALIZE_0D


  END MODULE M_ROCBURN_1D_ZN






