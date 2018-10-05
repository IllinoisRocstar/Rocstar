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
MODULE M_ROCBURN_2D 

!
! ---------------------------------------------------------------------------
!
!  This module contains INITIALIZE, UPDATE, and FINALIZE
!       to be registered with Roccom. 
!
!  Author:            K-C Tang, L. Massa, X. Jiao
!
! ---------------------------------------------------------------------------
!

  USE M_ROCBURN_INTERFACE_DATA
  USE M_CALCDIST

  IMPLICIT NONE
  INCLUDE 'comf90.h'
  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: MODEL_APN = 1, MODEL_PY = 2, MODEL_ZN = 3
  INTEGER, PARAMETER :: NO_TBL = 0

  CHARACTER(*), PARAMETER :: ioWin = "Burn"
  CHARACTER(*), PARAMETER :: intWin = "BurnInt"

! ---------------------------------------------------------------------------  

CONTAINS

  SUBROUTINE CHECK_ALLOC( ierr) 
    INTEGER, INTENT(IN) :: ierr
    IF(ierr /= 0) THEN
       PRINT *, "ROCBRUN ERROR: unable to allocate memory"
       CALL MPI_ABORT( MPI_COMM_WORLD, -1)
    END IF
  END SUBROUTINE CHECK_ALLOC

  SUBROUTINE INIT_WRAPPER(G_b, initial_time, comm, MAN_INIT, inSurf,  &
                          inInt, IN_obt_attr)

    TYPE(list_block), POINTER   :: G_b
    REAL(DBL), INTENT(IN)       :: initial_time
    INTEGER, INTENT(IN)         :: comm, MAN_INIT
    CHARACTER(*), INTENT(IN)    :: inSurf, inInt
    INTEGER, INTENT(IN)         :: IN_obt_attr
    
    G_b%MPI_COMM_ROCBURN  = comm
    G_b%pseudo_time = initial_time

    CALL COM_CALL_FUNCTION(G_b%INIT,6,MAN_INIT,inSurf,inInt, &
         G_b%INIT_0D, G_b%INIT_1D, IN_obt_attr)

  END SUBROUTINE INIT_WRAPPER

!
!     ------------------------------------------------------------------------
!                             INTERNAL PROCEDURES
!     ------------------------------------------------------------------------

  SUBROUTINE INITIALIZE( G_b, MAN_INIT, inSurf, inInt, INIT_0D, INIT_1D,  &
                         IN_obt_attr)

    TYPE(list_block), POINTER   :: G_b
    INTEGER, INTENT(IN)         :: MAN_INIT, IN_obt_attr
    CHARACTER(*), INTENT(IN)    :: inSurf, inInt 
!
!  INIT_0D and INIT_1D  are external subroutine arguments.  
!  Their dummy arguments are described as follow:
!

    INTERFACE
       SUBROUTINE INIT_0D( g_1d, comm, Indir, nxmax, To_read)
         USE M_ROCBURN_INTERFACE_DATA, ONLY : G_BURN_1D, DBL
         TYPE (G_BURN_1D), POINTER :: g_1d
         INTEGER, INTENT(IN)       :: comm
         CHARACTER(*), INTENT(IN)  :: Indir
         INTEGER, INTENT(OUT)      :: nxmax
         REAL(DBL), INTENT (OUT)   :: To_read
       END SUBROUTINE INIT_0D

       SUBROUTINE INIT_1D( g_1d, bflag, P, To, rhoc, p_coor, rb, &
                           Toa, fr, Tn, Tflame)
         USE M_ROCBURN_INTERFACE_DATA, ONLY : G_BURN_1D, DBL
         TYPE (G_BURN_1D), POINTER :: g_1d
         INTEGER, INTENT(INOUT)    :: bflag
         REAL(DBL), INTENT (IN)    :: P, To, rhoc, p_coor(3)
         REAL(DBL), INTENT (OUT)   :: rb, Toa, fr
         REAL(DBL), INTENT (OUT)   :: Tn(:)
         REAL(DBL), INTENT (OUT)   :: Tflame
       END SUBROUTINE INIT_1D
    END INTERFACE


! ---------------------------------------------------------------------------
! local variables
!
    INTEGER                      :: nxmax
    INTEGER                      :: ierror, ib, ic, nblks, bid,n_cell
    INTEGER, POINTER             :: blk_ids(:)
    TYPE(BLOCK), POINTER         :: blk
    LOGICAL                      :: is_APN, comp_filmcoeff
    REAL(DBL)                    :: zero, zerov(3)
!-----------------------------------------------------------------------------

    G_b%burn_iter = 0                    !start with zero also in case of restart
    CALL MPI_COMM_RANK(G_b%MPI_COMM_ROCBURN, G_b%rank, ierror)

!    IF(G_b%rank == 0 .AND. G_b%verbosity .gt. 1) THEN
!       WRITE(*,*)  'Rocburn: received initial_time= ',G_b%pseudo_time
!    END IF
 
!
!   Decode model_code
!
    is_APN                = G_b%burn_model == MODEL_APN
    comp_filmcoeff        = G_b%TBL_flag == NO_TBL

!    IF(G_b%rank == 0 .AND. G_b%verbosity .gt. 1) THEN
!       WRITE(*,*) 'Rocburn: TBL_flag      = ', G_b%TBL_flag
!       WRITE(*,*) 'Rocburn: burn_model    = ', G_b%burn_model
!    END IF

!
!   initialize global vairables (0D level) for inidividual 
!   combustion model
!
    CALL INIT_0D(G_b%g_1d, G_b%MPI_COMM_ROCBURN, TRIM(G_b%mname)//"/", nxmax, G_b%To_read)
    ALLOCATE (G_b%Tn( nxmax), STAT=ierror); CALL CHECK_ALLOC( ierror)

!
!   Create interface data in Roccom and allocate memory for them
!
    CALL COM_new_window( ioWin)
!   Use the subset of fluid or solid mesh.
!   It must use ghost nodes/cells as well in order to visualize.
    CALL COM_clone_dataitem( ioWin//".mesh", inSurf//".mesh")
    CALL COM_clone_dataitem( ioWin//'.bflag', inSurf//'.bflag')

!  
!   Incoming data
!

    CALL COM_new_dataitem( ioWin//".pf_alp", 'e', COM_DOUBLE, 1, "Pa")
    CALL COM_resize_array( ioWin//".pf_alp")

    IF ( .NOT. is_APN) THEN
!RAF       CALL COM_new_dataitem( ioWin//".centers", 'e', COM_DOUBLE, 3, "m")
       CALL COM_new_dataitem( ioWin//".qr_alp", 'e', COM_DOUBLE, 1, "W/m^2")
       CALL COM_new_dataitem( ioWin//".qc_alp", 'e', COM_DOUBLE, 1, "W/m^2")
       CALL COM_new_dataitem( ioWin//".rhos_alp", 'e', COM_DOUBLE, 1, "kg/m^3")
       CALL COM_new_dataitem( ioWin//".Tf_alp", 'e', COM_DOUBLE, 1, "K")
!!!    CALL COM_new_dataitem( ioWin//".To_alp",  'e', COM_DOUBLE, 1, "K")

!RAF       CALL COM_resize_array(ioWin//".centers")
       CALL COM_resize_array(ioWin//".qr_alp")
       CALL COM_resize_array(ioWin//".qc_alp")
       CALL COM_resize_array(ioWin//".rhos_alp")
       CALL COM_resize_array(ioWin//".Tf_alp")
    END IF
    CALL COM_new_dataitem( ioWin//".centers", 'e', COM_DOUBLE, 3, "m")
    CALL COM_resize_array(ioWin//".centers")
!
!   Outgoing data
!
    CALL COM_new_dataitem( ioWin//".rb", 'e', COM_DOUBLE, 1, "m/s")
    CALL COM_new_dataitem( ioWin//".Tflm", 'e', COM_DOUBLE, 1, "K")
    CALL COM_resize_array(ioWin//".rb")
    CALL COM_resize_array(ioWin//".Tflm")

    CALL COM_window_init_done( ioWin)

!
!   Create internal data that need to be saved for predictor-corrector
!   iterations or restart in Roccom and allocate memory for them.
!
    CALL COM_new_window( intWin)
    CALL COM_use_dataitem( intWin//".mesh", ioWin//".mesh")  ! Same mesh

!
!   Data for interpolation if subcycling for individual cells are needed
!
    IF ( .NOT. is_APN) THEN
       CALL COM_clone_dataitem( intWin//".pf_old", ioWin//".pf_alp")
       CALL COM_clone_dataitem( intWin//".qc_old", ioWin//".qc_alp")
       CALL COM_clone_dataitem( intWin//".qr_old", ioWin//".qr_alp")
       CALL COM_clone_dataitem( intWin//".rhos_old", ioWin//".rhos_alp")
       CALL COM_clone_dataitem( intWin//".Tf_old", ioWin//".Tf_alp")
!!!    CALL COM_clone_dataitem( intWin//".To_old",  ioWin//".To_alp")
    END IF

!
!   Profile history
!
    IF ( .NOT. is_APN) THEN
       CALL COM_clone_dataitem( intWin//".Toa", ioWin//".Tflm")
       IF ( comp_filmcoeff) THEN
          CALL COM_new_dataitem( intWin//".dist", 'e', COM_DOUBLE, 1, "m")
          CALL COM_resize_array(intWin//".dist")
       ENDIF
       CALL COM_new_dataitem( intWin//".temp", 'e', COM_DOUBLE, nxmax, "K")
       CALL COM_new_dataitem( intWin//".fr", 'e', COM_DOUBLE, 1, "")
       CALL COM_resize_array(intWin//".temp")
       CALL COM_resize_array(intWin//".fr")
    END IF

    CALL COM_window_init_done( intWin)

!
!   Get size information from Roccom
!
    CALL COM_get_panes( ioWin, nblks, blk_ids)
    ALLOCATE (G_b%blocks( nblks), STAT=ierror); CALL CHECK_ALLOC( ierror)

!
!  Obtain memory address from Roccom and build up the blocks
!
    DO ib = 1, nblks
       blk => G_b%blocks(ib)
       blk%iblock = blk_ids(ib)
       bid = blk_ids(ib)

!
!      incoming data
!
       CALL COM_get_size( ioWin//".pf_alp", bid, blk%nfaces)
       CALL COM_get_array( ioWin//".pf_alp", bid, blk%pres)
       CALL COM_get_array( ioWin//".bflag", bid, blk%burn_flag)

       IF ( .NOT. is_APN) THEN
!RAF          CALL COM_get_array( ioWin//".centers", bid, blk%coor)

          CALL COM_get_array( ioWin//".qr_alp", bid, blk%qr)
          CALL COM_get_array( ioWin//".qc_alp", bid, blk%qc)
          CALL COM_get_array( ioWin//".rhos_alp", bid, blk%rhoc)
          CALL COM_get_array( ioWin//".Tf_alp", bid, blk%Tg)
!!!       CALL COM_get_array( ioWin//".To_alp",   bid, blk%To)
       END IF
       CALL COM_get_array( ioWin//".centers", bid, blk%coor)

!
!      outgoing data
!

       CALL COM_get_array( ioWin//".rb", bid, blk%rb)
       CALL COM_get_array( ioWin//".Tflm", bid, blk%Tf)
       ! The above dataitems need be initialized by INIT_1D
       
!
!      Stored internal data
!
!      Data for interpolation if subcycling for individual cells are needed
!
!
       IF ( .NOT. is_APN) THEN
          CALL COM_get_array( intWin//".pf_old",   bid, blk%pres_old)
          CALL COM_get_array( intWin//".qc_old",   bid, blk%qc_old)
          CALL COM_get_array( intWin//".qr_old",   bid, blk%qr_old)
          CALL COM_get_array( intWin//".rhos_old", bid, blk%rhoc_old)
          CALL COM_get_array( intWin//".Tf_old",   bid, blk%Tg_old)
!!!       CALL COM_get_array( intWin//".To_old",   bid, blk%To_old)
       END IF

!
!      Profile history
!
       IF ( .NOT. is_APN) THEN
          CALL COM_get_array( intWin//".temp", bid, blk%temp)
          CALL COM_get_array( intWin//".Toa",  bid, blk%Toa)
          CALL COM_get_array( intWin//".fr",   bid, blk%fr)
          ! The above dataitems need be initialized by INIT_1D

          IF ( comp_filmcoeff) THEN
             CALL COM_get_array( intWin//".dist", bid, blk%dist)
             ! blk%dist should be initialized by CALCDIST_2D
          END IF
       END IF
    END DO  ! ib

!   Call Rocin to copy dataitems (without the mesh) into the new windows.
    CALL COM_call_function( IN_obt_attr,2, &
         COM_get_dataitem_handle_const(TRIM(inSurf)//".all"), &
         COM_get_dataitem_handle(TRIM(ioWin)//".all" ))
    CALL COM_call_function( IN_obt_attr,2, &
         COM_get_dataitem_handle_const(TRIM(inInt)//".data"), &
         COM_get_dataitem_handle(TRIM(intWin)//".data" ))

!   Call Rocman to prepare for data transfer, predictor-corrector iterations,
!   and restart. If initial_time is nonzero, Rocman will load data buffers
!   from restart files; if initial_time is zero, Rocman will initialize bflag,
!   pressure, qr, and qc.
    CALL COM_call_function( MAN_INIT, 3, ioWin, intWin, G_b%TBL_flag)

!   
!    initialize 1D level vairables for inidividual combustion model
!
!
    n_cell = 0
    IF ( G_b%pseudo_time == 0.0 .or. is_APN) THEN
       zero = 0.0
       zerov = zero 

       DO ib = 1, nblks
          blk => G_b%blocks(ib)

          IF ( .NOT. is_APN) THEN
             DO ic = 1, blk%nfaces
                n_cell = n_cell + 1
                CALL INIT_1D( G_b%g_1d, blk%burn_flag(ic), &
                     blk%pres(ic), G_b%To_read, blk%rhoc(ic), blk%coor(1:3,ic), &
                     blk%rb(ic), blk%Toa(ic), blk%fr(ic), blk%temp(:, ic), blk%Tf(ic))
             END DO    !  ic
             IF ( blk%nfaces>0) THEN
                blk%qr_old   = blk%qr
                blk%qc_old   = blk%qc
                blk%rhoc_old = blk%rhoc
                blk%Tg_old   = blk%Tg
                blk%pres_old = blk%pres            
             END IF
          ELSE
             DO ic = 1, blk%nfaces
                blk%burn_flag(ic) = 1
                CALL INIT_1D( G_b%g_1d, blk%burn_flag(ic), &
                     blk%pres(ic), zero, zero, blk%coor(1:3,ic),       &
                     blk%rb(ic), zero, zero, G_b%Tn(:), blk%Tf(ic))
             END DO    !  ic
          END IF    !APN
       
       END DO  ! ib

    ELSE

       DO ib = 1, nblks
          blk => G_b%blocks(ib)

          n_cell = n_cell + blk%nfaces
       ENDDO

    END IF   !initial_time == 0.0
    
    CALL MPI_ALLREDUCE(n_cell,G_b%total_cell,1,MPI_INTEGER,&
         MPI_SUM, G_b%MPI_COMM_ROCBURN,ierror)


!   Deallocate temporary buffer space
    CALL COM_free_buffer(blk_ids)
 
  END SUBROUTINE INITIALIZE
!*****************************************************************************


  SUBROUTINE FINALIZE( G_b) 
    TYPE(list_block), POINTER   :: G_b

    INTERFACE
       SUBROUTINE FINALIZE_0D( g_1d)
         USE M_ROCBURN_INTERFACE_DATA, ONLY : G_BURN_1D
         TYPE (G_BURN_1D), POINTER :: g_1d
       END SUBROUTINE FINALIZE_0D
    END INTERFACE
    
    CALL COM_delete_window( intWin)   ! Automaticall deallocate all the buffers
    CALL COM_delete_window( ioWin)    ! Automaticall deallocate all the buffers
    
    !   CALL FINALIZE_0D( G_b%g_1d)   ! Disabled because of error on Origin 2K

    !   Deallocate buffer space
    DEALLOCATE( G_b%Tn)
    DEALLOCATE( G_b%blocks)

  END SUBROUTINE FINALIZE

!
! ==========================================================================
!

  SUBROUTINE UPDATE_WRAPPER(G_b, timestamp, dt, MAN_UPDATE)


    TYPE(list_block), POINTER   :: G_b
    REAL(DBL), INTENT(IN)       :: timestamp, dt
    INTEGER, INTENT(IN)         :: MAN_UPDATE
    

    CALL COM_CALL_FUNCTION(G_b%UPDATE,6,timestamp,dt, MAN_UPDATE, &
         G_b%GET_FILM_COEFF, G_b%GET_TIME_STEP, G_b%GET_BURN_RATE)

  END SUBROUTINE UPDATE_WRAPPER

  SUBROUTINE UPDATE( G_b, timestamp, dt, MAN_UPDATE, GET_FILM_COEFF_1D, &
       GET_TIME_STEP_1D, GET_BURNING_RATE_1D)


!!!-------------------------------------------------
    TYPE(list_block), POINTER :: G_b
    REAL(DBL), INTENT (IN)    :: timestamp, dt
    INTEGER, INTENT(IN)       :: MAN_UPDATE

!!!
!!! get_film_coeff_1d, get_time_step1d, and get_burning_rate1d 
!!! are subroutine arguments.  Their dummy arguments are described as follow:

    INTERFACE
       SUBROUTINE GET_FILM_COEFF_1D( g_1d, p_coor, Ts, T_euler, P, Qc, Qcprime)
         USE  M_ROCBURN_INTERFACE_DATA, ONLY : G_BURN_1D, DBL
         TYPE (G_BURN_1D), POINTER   :: g_1d
         REAL(DBL), INTENT (IN)      :: p_coor(3), Ts, T_euler, P
         REAL(DBL), INTENT (OUT)     :: Qc,Qcprime
       END SUBROUTINE GET_FILM_COEFF_1D

       SUBROUTINE GET_TIME_STEP_1D( g_1d, rb, Toa, dt_max)
         USE  M_ROCBURN_INTERFACE_DATA, ONLY : G_BURN_1D, DBL
         TYPE (G_BURN_1D), POINTER   :: g_1d
         REAL(DBL), INTENT (IN)      :: rb, Toa
         REAL(DBL), INTENT (OUT)     :: dt_max
       END SUBROUTINE GET_TIME_STEP_1D

       SUBROUTINE GET_BURNING_RATE_1D ( g_1d, delt, P, To, Tn,   &
            qc, qc_old, qr, qr_old, rhoc, &
            Toa, rb, fr, bflag, Tnp1, Tflame, p_coor)
         USE  M_ROCBURN_INTERFACE_DATA, ONLY : G_BURN_1D, DBL
         TYPE (G_BURN_1D), POINTER   :: g_1d
         REAL(DBL), INTENT (IN)      :: delt, P, To
         REAL(DBL), INTENT (IN)      :: Tn(:)
         REAL(DBL), INTENT (IN)      :: qc, qc_old, qr, qr_old
         REAL(DBL), INTENT (IN)      :: rhoc
         REAL(DBL), INTENT (INOUT)   :: Toa, rb, fr
         INTEGER,   INTENT (INOUT)   :: bflag
         REAL(DBL), INTENT (OUT)     :: Tnp1(:)
         REAL(DBL), INTENT (OUT)     :: Tflame
         REAL(DBL), INTENT (IN)      :: p_coor(3)
       END SUBROUTINE GET_BURNING_RATE_1D
    END INTERFACE

!!! 
!!! local variables 
!!!
    TYPE(BLOCK), POINTER :: blk
    INTEGER     :: ic, ib, one_int, ierror          !Dummy indexes
    LOGICAL     :: is_APN, comp_filmcoeff
    REAL(DBL)   :: zero, one, ten
    REAL(DBL)   :: dt_max, dt_mks
    INTEGER     :: nblks, n_subcycle_0d, i_subcycle_0d
    INTEGER     :: i_pseudo_iter, max_pseudo_iter, n_cell_ignited
    REAL(DBL)   :: inv_n_subcycle_0d, alpha
    REAL(DBL)   :: delta_P, delta_qc, delta_qr, delta_Tg, delta_rhoc
    REAL(DBL)   :: qr_old_mks, qc_old_mks, qr_mks, qc_mks, P_mks, rhoc_mks
    REAL(DBL)   :: Tflame_APN,pre_APN,exp_APN,out_APN
!!!----------------------------------------------------------------------------------------- 
!!!
!!!   point blk to first data patch
!!!

    IF ( ASSOCIATED( G_b%blocks)) THEN
       nblks = UBOUND( G_b%blocks, 1)
    ELSE
       nblks = 0
    END IF

    is_APN                = G_b%burn_model == MODEL_APN
    comp_filmcoeff        = G_b%TBL_flag == NO_TBL
    zero=0.0
    one = 1.0
    ten = 10.0

    alpha=1.0d0
    CALL COM_call_function( MAN_UPDATE, 1, alpha)


    G_b%pseudo_time = timestamp
    G_b%burn_iter = G_b%burn_iter + 1
    n_cell_ignited = 0

    IF (is_APN) THEN   !0D MODEL
!RAF
!RAF For best speed, we could use this optimization if nmat = 1
!RAF
!RAF       one_int = 1
!RAF       CALL GET_BURNING_RATE_1D ( G_b%g_1d, zero, one, &
!RAF            zero, G_b%Tn, zero, zero, zero, zero, zero, zero, &
!RAF            pre_APN, zero, one_int, G_b%Tn, Tflame_APN)                   !get the prexponential term
!RAF
!RAF       CALL GET_BURNING_RATE_1D ( G_b%g_1d, zero, ten, &
!RAF            zero, G_b%Tn, zero, zero, zero, zero, zero, zero, &
!RAF            out_APN, zero, one_int, G_b%Tn, Tflame_APN)                  !A*10^n
!RAF       IF ( out_APN == 0) THEN
!RAF          exp_APN = 0
!RAF       ELSE
!RAF          exp_APN = log10(out_APN/pre_APN)                            !get the exponent
!RAF       END IF

       DO ib=1, nblks
          blk => G_b%blocks(ib)
          DO ic =1, blk%nfaces
            blk%burn_flag(ic) = 1

!RAF             blk%rb(ic) = pre_APN*blk%pres(ic)**exp_APN
!RAF             blk%Tf(ic) = Tflame_APN                                 !This assumes no ignition model for APN

             CALL GET_BURNING_RATE_1D ( G_b%g_1d, zero, blk%pres(ic), &
                  zero, G_b%Tn,   &
                  zero, zero,     &
                  zero, zero,     &
                  zero, zero,     &
                  blk%rb(ic), zero, &
                  blk%burn_flag(ic), G_b%Tn, blk%Tf(ic), blk%coor(1:3,ic))


          ENDDO
       ENDDO

    ELSE         ! 1D MODELS


       DO ib = 1, nblks
          blk => G_b%blocks(ib)

          IF( comp_filmcoeff .AND. blk%nfaces>0) &
               CALL CALCDIST_2D( G_b, blk%coor, blk%dist)

          DO ic = 1, blk%nfaces

             CALL GET_TIME_STEP_1D( G_b%g_1d, blk%rb(ic), blk%Toa(ic), dt_max)
!!!             G_b%To_read = G_b%To(ic)

             n_subcycle_0d = INT(dt/dt_max) + 1
              
             IF ( n_subcycle_0d == 1) THEN   ! no subcycling needed

                G_b%Tn = blk%temp(:,ic)  !condensed phase temperature profile (K)

                IF( comp_filmcoeff) THEN
                   CALL GET_FILM_COEFF_1D( G_b%g_1d, blk%coor(1:3,ic), G_b%Tn(1), &
                        blk%Tg(ic), blk%pres(ic), blk%qc(ic),blk%qr(ic))
                END IF
                CALL GET_BURNING_RATE_1D ( G_b%g_1d, dt, blk%pres(ic), &
                     G_b%To_read, G_b%Tn,   &
                     blk%qc(ic), blk%qc_old(ic), &
                     blk%qr(ic), blk%qr_old(ic), &
                     blk%rhoc(ic), blk%Toa(ic), &
                     blk%rb(ic), blk%fr(ic), &
                     blk%burn_flag(ic), blk%temp(:,ic), blk%Tf(ic),  &
                     blk%coor(1:3,ic))

             ELSE
!!!
!!!          local subcycling needed
!!!
                G_b%Tn  =  blk%temp(:,ic)  ! condensed phase temperature profile

                IF( comp_filmcoeff ) THEN
                   CALL GET_FILM_COEFF_1D( G_b%g_1d, blk%coor(1:3,ic), G_b%Tn(1), &
                        blk%Tg(ic), blk%pres_old(ic), blk%qc(ic),blk%qr(ic))
                END IF

                inv_n_subcycle_0d = 1.0/float(n_subcycle_0d)  
                dt_mks   = dt*inv_n_subcycle_0d
                delta_P    = blk%pres(ic) - blk%pres_old(ic)


!!$                delta_qc   = blk%qc(ic)   - blk%qc_old(ic)
!!$                delta_qr   = blk%qr(ic)   - blk%qr_old(ic)
                delta_Tg   = blk%Tg(ic)   - blk%Tg_old(ic)
                delta_rhoc = blk%rhoc(ic) - blk%rhoc_old(ic)

                qr_old_mks   =  blk%qr_old(ic)
                qc_old_mks   =  blk%qc_old(ic)


                DO i_subcycle_0d = 1, n_subcycle_0d
                   alpha =  float(i_subcycle_0d)*inv_n_subcycle_0d
                   P_mks = blk%pres_old(ic) + delta_P*alpha

                   qr_mks = blk%qr(ic)
                   qc_mks = blk%qc(ic)
                   rhoc_mks = blk%rhoc_old(ic) + delta_rhoc*alpha

                   CALL GET_BURNING_RATE_1D ( G_b%g_1d, dt_mks, P_mks, &
                        G_b%To_read, G_b%Tn,   &
                        qc_mks, qc_old_mks, qr_mks, qr_old_mks, &
                        rhoc_mks, blk%Toa(ic), blk%rb(ic), blk%fr(ic), &
                        blk%burn_flag(ic), blk%temp(:,ic), blk%Tf(ic), &
                        blk%coor(1:3,ic))

!!$                   qr_old_mks  = qr_mks
!!$                   qc_old_mks  = qc_mks
                   G_b%Tn      = blk%temp(:,ic)

                END DO  ! i_subcycle_0d

             END IF  ! if n_subcycle_0d ==1 

             n_cell_ignited = n_cell_ignited +  blk%burn_flag(ic)

          END DO  ! Cells

          IF ( blk%nfaces>0) THEN
             blk%pres_old = blk%pres
             !!     blk%To_old   = blk%To
             blk%qc_old   = blk%qc
             blk%qr_old   = blk%qr
             blk%Tg_old   = blk%Tg
             blk%rhoc_old = blk%rhoc
          END IF
       END DO  ! blocks
!
!      ROCBURNPY STD OUTPUT
!
       CALL MPI_ALLREDUCE(n_cell_ignited,G_b%burn_cell,1,MPI_INTEGER,&
            MPI_SUM, G_b%MPI_COMM_ROCBURN,ierror)
       ! This statement needs access to the verbosity
       IF(G_b%rank == 0)  &
            write(*,*)'ROCBURN iter :: ',G_b%burn_iter,'CELLS IGNITED',G_b%burn_cell,'PERCENT',&
               dble(G_b%burn_cell)/dble(G_b%total_cell)*100.0d0,G_b%total_cell

   END IF     !IF is_APN



  END SUBROUTINE UPDATE


!     -------------------------------------------------------------------
!                       END OF INTERNAL PROCEDURES
!     -------------------------------------------------------------------

END MODULE M_ROCBURN_2D






