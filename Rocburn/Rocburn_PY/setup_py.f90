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
MODULE  SETUP_PY
  USE data_py

CONTAINS


!****************************************************************
  SUBROUTINE get_time_step_1d(bp,rb,Toa,dt_max) 
      IMPLICIT NONE  
!---------------------------------------------------------------------
      TYPE(parameter_structure),POINTER  :: bp
      REAL(DBL), INTENT(IN)  :: rb,Toa
      REAL(DBL), INTENT(OUT) :: dt_max
!------------------------------------------------------------------------

      dt_max = bp%delt_max 

!------------------------------------------------------------------------
      RETURN
    END SUBROUTINE get_time_step_1d
!***********************************************************************
      




!***************************************************************************
    SUBROUTINE burn_init_0d(bp, comm, IN_DIR, nx_read,To_read)
 
      INCLUDE 'mpif.h'

      TYPE(parameter_structure),POINTER  :: bp
      INTEGER, INTENT(IN)    :: comm
      REAL(DBL), INTENT(OUT) :: To_read
      INTEGER , INTENT(OUT)  :: nx_read
      INTEGER                :: ierror, rank
      CHARACTER(*), INTENT(IN)    :: IN_DIR  
!----------------------------------------------------------------------------

      CALL MPI_COMM_RANK(comm, rank, ierror)

      ALLOCATE(bp)
      bp%comm = comm
!
!     read propellant thermophysical properties
!
      CALL read_properties(bp,IN_DIR)
!
!     set return values of this subroutine
!
      nx_read = bp%numx
      To_read = bp%To
!
!     grid generation
!
      CALL grid(bp)
      if(bp%TABUSE == 1) CALL readtable(bp, IN_DIR)


      RETURN 
    END SUBROUTINE burn_init_0d
!********************************************************

!
!==============================================================================
!
    SUBROUTINE read_properties(bp,IN_DIR)

      INCLUDE 'mpif.h'

      TYPE(parameter_structure),POINTER :: bp
      CHARACTER*(*), INTENT(IN)          :: IN_DIR
      INTEGER                           :: ir
      CHARACTER(LEN=90)                 :: filnam
      INTEGER                           :: ioerr,rank,ierror
!---------------------------------------------------------------------

      CALL MPI_COMM_RANK(bp%comm, rank, ierror)        

      filnam = 'RocburnPYControl.txt'
      IF (IN_DIR(LEN_TRIM(IN_DIR):LEN_TRIM(IN_DIR)) == '/')  THEN
         filnam = TRIM(IN_DIR) // TRIM(filnam)
      ELSE
         filnam = TRIM(IN_DIR) // '/' // TRIM(filnam)
      ENDIF

      filnam = TRIM(filnam)
      ir = 10


      OPEN (unit=ir,file=filnam,status='old',IOSTAT=ioerr)
      if(ioerr > 0) THEN
        write(*,*)'LOOKING FOR file ',filnam
        write(*,*)'FILE NOT FOUND'
        CALL MPI_Abort (bp%comm, 1, ierror)
      endif



      READ(ir,*)   bp%a_p
      IF (rank==0) WRITE(*,*) 'ROCBURN:  a_p =', bp%a_p

      READ(ir,*)   bp%n_p
      IF (rank==0) WRITE(*,*) 'ROCBURN:  n_p =', bp%n_p

      READ(ir,*)   bp%Pref
      IF (rank==0) WRITE(*,*) 'ROCBURN:  Pref =', bp%Pref

      READ(ir,*)   bp%Ac
      IF (rank==0) WRITE(*,*) 'ROCBURN:  Ac= ',bp%Ac

      READ(ir,*)   bp%eg_ru
      IF (rank==0) WRITE(*,*) 'ROCBURN:  eg_ru= ',bp%eg_ru

      READ(ir,*)   bp%ec_ru
      IF (rank==0) WRITE(*,*) 'ROCBURN:  ec_ru ',bp%ec_ru

      READ(ir,*)   bp%alfac
      IF (rank==0) WRITE(*,*) 'ROCBURN:  alfac= ',bp%alfac

      READ(ir,*)   bp%C
      IF (rank==0) WRITE(*,*) 'ROCBURN:  C= ',bp%C

      READ(ir,*)   bp%lamg
      IF (rank==0) WRITE(*,*) 'ROCBURN:  lamg= ',bp%lamg

      READ(ir,*)   bp%delt
      IF (rank==0) WRITE(*,*) 'ROCBURN:  delt= ',bp%delt

      READ(ir,*)   bp%igrid
      IF (rank==0) WRITE(*,*) 'ROCBURN:  igrid= ',bp%igrid

      READ(ir,*)   bp%numx
      IF (rank==0) WRITE(*,*) 'ROCBURN:  numx= ',bp%numx

      READ(ir,*)   bp%xmax
      IF (rank==0) WRITE(*,*) 'ROCBURN:  xmax= ',bp%xmax

      READ(ir,*)   bp%beta
      IF (rank==0) WRITE(*,*) 'ROCBURN:  beta= ',bp%beta

      READ(ir,*)   bp%Tstar0
      IF (rank==0) WRITE(*,*) 'ROCBURN:  Tf_adiabatic= ', bp%Tstar0

      READ(ir,*)   bp%To
      IF (rank==0) WRITE(*,*) 'ROCBURN:  To= ',bp%To

      READ(ir,*)   bp%Tignition
      IF (rank==0) WRITE(*,*) 'ROCBURN:  Tignition= ',bp%Tignition

      READ(ir,*)   bp%Tsurf
      IF (rank==0) WRITE(*,*) 'ROCBURN:  Initial Temp = ',bp%Tsurf

      read(ir,*)   bp%film_cons
      IF (rank==0) WRITE(*,*) 'ROCBURN:  film_cons= ',bp%film_cons

      read(ir,*)   bp%ixsymm
      IF (rank==0) WRITE(*,*) 'ROCBURN: ixsymm = ',bp%ixsymm

      read(ir,*)   bp%x_surf_burn
      IF (rank==0) WRITE(*,*) 'ROCBURN: x_surf_burn = ',bp%x_surf_burn

      read(ir,*, ERR = 101,END = 101)   bp%P_range(2)
      IF (rank==0) WRITE(*,*) 'ROCBURN: press_max = ',bp%P_range(2)
      read(ir,*, ERR = 101,END = 101)   bp%P_range(1)
      IF (rank==0) WRITE(*,*) 'ROCBURN: press_min = ',bp%P_range(1)
      if(bp%P_range(2) + bp%P_range(1) < 100.0) goto 101

      read(ir,*, ERR = 101,END = 101)   bp%rb_range(2)
      IF (rank==0) WRITE(*,*) 'ROCBURN: rb_max = ',bp%rb_range(2)
      read(ir,*, ERR = 101,END = 101)   bp%rb_range(1)
      IF (rank==0) WRITE(*,*) 'ROCBURN: rb_min = ',bp%rb_range(1)

      read(ir,*, ERR = 101,END = 101)   bp%Tf_range(2)
      IF (rank==0) WRITE(*,*) 'ROCBURN: Tf_max = ',bp%Tf_range(2)
      read(ir,*, ERR = 101,END = 101)   bp%Tf_range(1)
      IF (rank==0) WRITE(*,*) 'ROCBURN: Tf_min = ',bp%Tf_range(1)

      READ(ir,*, ERR = 102,END = 102)   bp%TABUSE
      IF (rank==0) WRITE(*,*) 'ROCBURN:  TABUSE =', bp%TABUSE

      READ(ir,*, ERR = 102,END = 102)   bp%TABNAM
      IF (rank==0) WRITE(*,*) 'ROCBURN:  TABNAM =', bp%TABNAM

      CLOSE(ir)

      goto 199

! RANGE NOT FOUND
101   CONTINUE
      CLOSE(ir)
      IF (rank==0) WRITE(*,*) 'MISSING variable RANGES, ...&
           & probably you are using an OLD input deck'
      IF (rank==0) WRITE(*,*) 'Hardwiring limits'
      bp%P_range = (/1.d3,1.d8/)
      bp%rb_range = (/-1.d-9,1.d2/)
      bp%Tf_range = (/290.0d0,10000.0d0/)
      goto 199

!NO TABLE
102   CONTINUE
      CLOSE(ir)
      IF (rank==0) WRITE(*,*) 'MISSING TABNAM TABUSE'
      bp%TABUSE = 0
      goto 199


199   CONTINUE

! Check input deck
      IF (rank==0) CALL CHECK_INPUT_RANGE(bp)

!--------------------------------------------------------------------------
        RETURN 
      END SUBROUTINE read_properties
!******************************************************

!=============================================================================
!   S_U_B_R_O_U_T_I_N_E  grid
!   arguments   :                                                             c
!                                                                             c
!      gridtype : = 1  expoential grid                                        c
!                 = 2  boundary layer grid                                    c
!      numx     : maximum number of grid points                               c
!      beta     : stretch parameter                                           c
!      xmax     : the maximum dimension in x-coordinate                       c
!      x        : the x-coordinate array                                      c
!      z        : the z-coordinate arary                                      c
!      zx       : the first derivative of z                                   c
!      zxx      : the second derivative of z                                  c
!      delz     : the z spacing                                               c
!                                                                             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE grid(bp)

        IMPLICIT NONE  
        INCLUDE 'mpif.h'
!----------------------------------------------------------------------
        TYPE(parameter_structure),POINTER  :: bp
        INTEGER            :: gridtype, numx
        REAL(DBL)          :: bpm,bp1,bm1,czx
        REAL(DBL)          :: czxx,term,term2
        REAL(DBL)          :: beta,xmax
        REAL(DBL), POINTER :: x(:), z(:), zx(:), zxx(:)
        INTEGER            :: i,rank,ierror
!----------------------------------------------------------------------
        
        CALL MPI_COMM_RANK(bp%comm, rank, ierror)        

        gridtype = bp%igrid
        numx = bp%numx
        beta = bp%beta
        xmax = bp%xmax
        bp%delz = one/dble(numx-1)


        ALLOCATE(bp%x(numx))
        ALLOCATE(bp%z(numx))
        ALLOCATE(bp%zx(numx))
        ALLOCATE(bp%zxx(numx))

      
!
!     This is not very fancy but will save time. this sub is called one per run
!
        x => bp%x
        z => bp%z
        zx => bp%zx
        zxx => bp%zxx

        IF(rank ==0) THEN
           PRINT *,' ROCBURN: delz=',bp%delz
           PRINT *,' ROCBURN: nxmax=',numx
        END IF

        IF (gridtype.eq.1) THEN
!
!        exponential grid
!
           DO i=1,numx
              z(i) = dble(i-1)*bp%delz
              x(i) = log(one-z(i))/beta
              IF (i.eq.numx) THEN
                 x(i) = -10.
              END IF
              zx(i)  = -beta*exp(beta*x(i))
              zxx(i) = -beta**two*exp(beta*x(i))
           ENDDO

        ELSE
!
!        Anderson, Tanehill and Pletcher boundary layer grid 
!

           bp1=beta+1.0
           bm1=beta-1.0
           bpm=bp1/bm1
           czx=2.0*beta/log(bpm)/xmax
           czxx=-2.0/xmax*czx
           DO  i=1,numx
              z(i)=dble(i-1)*bp%delz
              term=bpm**(one-z(i))
              x(i)=xmax*(bp1-bm1*term)/(term+one)
              term2=beta**2-(one-x(i)/xmax)**2
              zx(i)=czx/term2
              zxx(i)= czxx*(one-x(i)/xmax)/term2**2
           ENDDO
           
        END IF

!  EVALUATE MAXIMUM DT for time integration

        bp%delt_max = half*bp%delz*bp%delz/maxval(abs(zx))**2/bp%alfac
        bp%delz2inv = half/bp%delz
        bp%delzsqinv = one / (bp%delz*bp%delz)
        bp%dx1 = bp%delz/zx(1)

        IF(rank == 0) THEN
           WRITE (*,*)
           WRITE (*,*)' ROCBURN: Smallest delta x = ',x(1)-x(2)
           WRITE (*,*)' ROCBURN: Largest delta x  = ',x(numx-1)-x(numx)
           WRITE (*,*)' ROCBURN: MAXIMUM DELTA t = ',bp%delt_max
           WRITE (*,*)
        END IF

        NULLIFY(x,z,zx,zxx)
!--------------------------------------------------
        RETURN
      END SUBROUTINE grid
!******************************************************


!******************************************************
      SUBROUTINE READTABLE(bp,IN_DIR)

        IMPLICIT NONE  
        INCLUDE 'mpif.h'
!----------------------------------------------------------------------
        TYPE(parameter_structure),POINTER  :: bp
        REAL(DBL) ::  Ts,pread,g,rb,fx2,alfaEFF
        REAL(DBL) ::  Tss,rss,alpTAB
        REAL(DBL) ::  out0,out1,out_1,err
        REAL(DBL),POINTER  :: vtmp(:)
        INTEGER :: i,j,l,inttmp,n1,n2,n3,n4,count, ierror, rank,ioerr
        CHARACTER*(*), INTENT(IN)          :: IN_DIR
        CHARACTER(LEN=90)                 :: filnam
!----------------------------------------------------------------------

        CALL MPI_COMM_RANK(bp%comm, rank, ierror)        

        allocate(Bp%TABLE)

        filnam = trim(Bp%TABNAM)
        IF (IN_DIR(LEN_TRIM(IN_DIR):LEN_TRIM(IN_DIR)) == '/')  THEN
           filnam = TRIM(IN_DIR) // TRIM(filnam)
        ELSE
           filnam = TRIM(IN_DIR) // '/' // TRIM(filnam)
        ENDIF

!--open
        open(32,file=filnam,status='old',IOSTAT=ioerr)
        if(ioerr > 0) THEN
           write(*,*)'LOOKING FOR file ',filnam
           write(*,*)'FILE NOT FOUND'
           write(*,*)'Use table flag:: ',bp%TABUSE
           CALL MPI_Abort (bp%comm, 1, ierror)
        endif


        read(32,*)Bp%TABLE%nx_table,Bp%TABLE%ny_table,Bp%TABLE%nfield,inttmp
        n1 = Bp%TABLE%nx_table; n2 = Bp%TABLE%ny_table; 
        n3 = max(n1,n2); n4 = Bp%TABLE%nfield

!----->allocation table fields
        allocate(Bp%TABLE%tsurf00    (n1))
        allocate(Bp%TABLE%press00    (n2))
        allocate(Bp%TABLE%heatflux00 (n1,n2))
        allocate(Bp%TABLE%rb00       (n1,n2))
        allocate(Bp%TABLE%fxsq00     (n1,n2))
        allocate(Bp%TABLE%Tgas00     (n1,n2))
        allocate(Bp%TABLE%Tstd00     (n2))
        allocate(Bp%TABLE%rstd00     (n2))
        allocate(Bp%TABLE%alph00     (n1,n2))
!----->allocation work fields
        allocate(Bp%TABLE%wrk1       (n3))
        allocate(Bp%TABLE%wrk2       (n3))
        allocate(Bp%TABLE%wrk3       (n3))
        allocate(Bp%TABLE%wrk4       (n3))
        allocate(Bp%TABLE%wrk5       (n3))
        allocate(Bp%TABLE%wrk6       (n3))
        allocate(Bp%TABLE%wrk7       (n3))
        allocate(Bp%TABLE%wrk8       (n3))
        allocate(Bp%TABLE%wrk9       (n3))
        allocate(Bp%TABLE%wrk10      (n3))

        allocate(vtmp                (n4))

! NOTE reacquire a value for alpha from the table

        read(32,*)alpTAB,Bp%TABLE%chi

        if(rank == 0 .and. abs(alpTAB-bp%alfac) > 1d-10 ) &
             write(*,*) 'WARNING changing ALPHA',bp%alfac,alpTAB
        
        bp%alfac = alpTAB
        
        do j=1,Bp%TABLE%ny_table
           do i=1,Bp%TABLE%nx_table
              read(32,*)  (vtmp(l),l=1,n4)

              Bp%TABLE%tsurf00(i) = vtmp(1)
              Bp%TABLE%press00(j) = vtmp(2)
              Bp%TABLE%heatflux00(i,j) = vtmp(3) 
              Bp%TABLE%rb00(i,j) = vtmp(4) 
              Bp%TABLE%fxsq00(i,j) = vtmp(5)
              alfaEFF = bp%alfac*( one + Bp%TABLE%chi + Bp%TABLE%fxsq00(i,j) )
              Bp%TABLE%alph00(i,j) = alfaEFF
              if(n4 > 5) then
                 Bp%TABLE%Tgas00(i,j) = vtmp(6)
              else 
                 Bp%TABLE%Tgas00(i,j) =  bp%Tstar0 - bp%To + Bp%TABLE&
                      %tsurf00(i) - alfaEFF/Bp%TABLE%rb00(i,j)*Bp%TABLE%heatflux00(i,j)
              endif
              if(n4 > 6) then
                 Bp%TABLE%Tstd00(j) = vtmp(7)  !the last value in col
                 Bp%TABLE%rstd00(j) = vtmp(8)
              endif
           enddo
        enddo

        Bp%TABLE%spline = inttmp == 1
        Bp%TABLE%small_1 = (Bp%TABLE%tsurf00(2)-Bp%TABLE%tsurf00(1))/100.d0
        Bp%TABLE%small_2 = (Bp%TABLE%press00(2)-Bp%TABLE%press00(1))&
             /100.d0


        if(n4 <= 6) then
           do j = 1,Bp%TABLE%ny_table

              err = 1d99
              count = 0
              if( j <=2) then
                 Tss = (Bp%TABLE%tsurf00(1) +  Bp%TABLE%tsurf00(Bp%TABLE&
                      %nx_table))/2
              else
                 Tss = Bp%TABLE%Tstd00(j-1) +  &
                      (Bp%TABLE%Tstd00(j-1) - Bp%TABLE%Tstd00(j-2))/&
                      (Bp%TABLE%press00(j-1) - Bp%TABLE%press00(j-2))*&
                      (Bp%TABLE%press00(j) - Bp%TABLE%press00(j-1))
              endif
                 
              do while (count < 20 .and. err >  Bp%TABLE%small_1)
                 count = count +1
                 call STEADYTEMP(bp,Tss,bp%TABLE%press00(j),out0,rss,g)
                 call STEADYTEMP(bp,Tss-Bp%TABLE%small_1/2d0,bp%TABLE&
                      %press00(j),out_1,rss,g)
                 call STEADYTEMP(bp,Tss+Bp%TABLE%small_1/2d0,bp%TABLE&
                      %press00(j),out1,rss,g)
                 Tss = Tss - out0/(out1-out_1)*Bp%TABLE%small_1
                 err = abs(out0)
              enddo
              Bp%TABLE%Tstd00(j) = Tss
              call STEADYTEMP(bp,Tss,bp%TABLE%press00(j),out0,rss,g)
              Bp%TABLE%rstd00(j) = rss
              if(rank == 0) write(*,'(i3,1p3e12.3,0p,A11,1p5e12.3)')count,bp%TABLE%press00(j),Tss,rss, &
                      '<--- steady',out0,g
           enddo

        endif

        close(32)
        deallocate(vtmp)

!--------------------------------------------------
        RETURN
      END SUBROUTINE READTABLE
!******************************************************

!******************************************************
      SUBROUTINE STEADYTEMP(bp,TIN,pIN,out,rb,g)

        IMPLICIT NONE  
        INCLUDE 'mpif.h'
!----------------------------------------------------------------------
        TYPE(parameter_structure),POINTER  :: bp
        INTEGER :: jj,kk
        REAL(DBL) ::  TIN,pIN,out
        REAL(DBL) ::  g,rb,To,dy,alfa,fx2
!----------------------------------------------------------------------

        To = bp%To


        jj = 0;kk = 0;
        call polin2(bp,bp%TABLE%tsurf00,bp%TABLE%press00,bp%TABLE%heatflux00,&
             bp%TABLE%nx_table,bp%TABLE%ny_table,Tin,pIN,g,jj,kk)
        call polin2(bp,bp%TABLE%tsurf00,bp%TABLE%press00,bp%TABLE%rb00,  &
             bp%TABLE%nx_TABLE,bp%TABLE%ny_table,Tin,pIN,rb,jj,kk)
        call polin2(bp,bp%TABLE%tsurf00,bp%TABLE%press00,bp%TABLE%fxsq00, &
             bp%TABLE%nx_TABLE,bp%TABLE%ny_table,Tin,pIN,fx2,jj,kk)
        alfa = bp%alfac*( one + Bp%TABLE%chi + fx2)

        out = g - (Tin - To)*rb/alfa

        RETURN
      END SUBROUTINE STEADYTEMP
!******************************************************

!****************************************************
      SUBROUTINE CHECK_INPUT_RANGE(bp)
        
        IMPLICIT NONE  
!----------------------------------------------------------------------
        TYPE(parameter_structure),POINTER  :: bp
        INTEGER :: ierror
        LOGICAL :: check(20)
!-----------------------------------------------------------------------
        check = .TRUE.
        check(1)  = bp%a_p > zero
        check(2)  = bp%n_p >zero .AND. bp%n_p < 10.0
        check(3)  = bp%Pref > zero .AND. bp%Pref < 200.0
        check(4)  = bp%Ac > 100.0
        check(5)  = bp%eg_ru > 100.0
        check(6)  = bp%ec_ru >100.0
        check(7)  = bp%alfac > zero .AND. bp%alfac < 1.0
        check(8)  = bp%C > zero .AND. bp%C < 1.0
        check(9)  = bp%lamg > zero .AND. bp%lamg < 1.0
        check(10) = bp%To > 100.0  .AND. bp%To < 1000.0
        check(11) = bp%Tignition > bp%To
        check(12) = bp%Tstar0 > bp%Tignition .AND. bp%Tstar0 < 1.d4
        check(13) = bp%ixsymm >= 0 .AND. bp%ixsymm < 3 
        check(14) = bp%Tsurf > 100.0 .AND. bp%Tsurf < bp%Tstar0

        IF(.NOT. ALL(check) ) THEN
           write(*,*) 'ROCBURN CHECK_INPUT_RANGE'
           write(*,*) 'INPUTS OUT OF RANGE'
           write(*,*) 'CHECK', check
           CALL MPI_Abort (bp%comm, 1, ierror)
           STOP
        ELSE
           write(*,*)'ALL VARIABLES IN SPECIFIED range'
        ENDIF

!-----------------------------------------------------------------------
        RETURN
      END SUBROUTINE CHECK_INPUT_RANGE
!===========================================================================
!     


!***************************************************************************
    SUBROUTINE burn_finalize_0d(bp)

      IMPLICIT NONE   
!-----------------------------------------------------------------------
      TYPE(parameter_structure),POINTER  :: bp
!----------------------------------------------------------------------------

      DEALLOCATE( bp)

    END SUBROUTINE burn_finalize_0d

  END MODULE SETUP_PY
