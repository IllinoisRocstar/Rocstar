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
MODULE  UPDATE_PY
  USE data_py

CONTAINS

  SUBROUTINE BURN_GET_BURNING_RATE1D (bp,delt,P,To,Tn,Qc,Qc_old,Qr,Qr_old,  &
                                      rhoc,Toa,rb,fr,bflag,Tnp1,Tflame,coor)
 
    IMPLICIT NONE

    INCLUDE 'mpif.h'
!----------------------------------------------------------------------------

!   Dummy Variables
!   
    TYPE(parameter_structure),POINTER  :: bp

    REAL(DBL),INTENT(IN)      :: delt,P,To,rhoc
    REAL(DBL),INTENT(IN)      :: Qc,Qc_old,Qr,Qr_old
    REAL(DBL), INTENT (INOUT) :: Toa, rb, fr
    INTEGER,INTENT(INOUT)     :: bflag
    REAL(DBL),INTENT(OUT)     :: Tflame
    REAL(DBL),INTENT(IN)      :: Tn(:)
    REAL(DBL),INTENT(OUT)     :: Tnp1(:)
    REAL(DBL),INTENT(IN)      :: coor(3)
!
!   Local Variables
!  
    INTEGER   :: i, nx
    REAL(DBL) :: rb_old
    REAL(DBL) :: first,second
    REAL(DBL) :: add,coe,rhs
    INTEGER   :: ierror
    TYPE(work_structure) :: bw
!------------------------------------------------------------------

!  CHECK the input for divergence

    IF ( .TRUE. .AND. (.NOT. (P > bp%P_range(1) .AND. P < bp%P_range(2)) )) THEN
       write(*,*)"INPUT PRESSURE OUT OF RANGE",P,bp%P_range
       CALL MPI_Abort (bp%comm, 1, ierror)
       STOP
    ENDIF
!
!  convert input variables
!
    bw%P       = P*pa2atm
    bw%rhoc    = rhoc*kgmc2gcc
    bw%Qc      = Qc*j_msq2cal_cmsq
    bw%Qcprime = Qr*j_msq2cal_cmsq
    bw%Ts      = Tn(1)
    bw%To      = To
    bw%lamc    = bw%rhoc * bp%alfac * bp%C
    nx         = bp%numx

    CALL MFUN(bp,bw)

!!    rb_old = bw%rb
    IF (bflag == 0) THEN                                !KJM
       rb_old = merge(bw%rb, zero, bw%Ts >= bp%Tignition)
    ELSE                                                !KJM
       rb_old = bw%rb                                   !KJM
    ENDIF                                               !KJM

    DO  i=2,nx-1

       first  = (Tn(i+1)-Tn(i-1))*bp%delz2inv
       second = (Tn(i+1)-2.0d0*Tn(i)+Tn(i-1))*bp%delzsqinv

       Tnp1(i)=Tn(i)+delt*(                                    &
            (bw%alfa_eff*bp%zxx(i) - rb_old*bp%zx(i)) * first  &
            + bw%alfa_eff*bp%zx(i)**2 * second)       
    ENDDO



    Tnp1(nx)=To

    IF (bflag == 0) THEN                                !KJM
       bw%ignited = (bw%Ts > bp%Tignition) 
    ELSE                                                !KJM
       bw%ignited = .TRUE.                              !KJM
    ENDIF                                               !KJM

    CALL GFUN(bp,bw)

!   update suface B.C. for this iteration
    coe = three + two*bw%fsprime*bp%dx1;
    rhs = two*bp%dx1*(bw%fs - bw%fsprime*bw%Ts);
    add = four*Tnp1(2) - one*Tnp1(3);
    ! third        coe = 11.0d0 + fsprime*6.0*dx1;
    ! third        rhs = 6.0*dx1*(fs - fsprime*Ts)
    ! third        add = 18.0*Tnp1(2)-9.0*Tnp1(3)+2.0*Tnp1(4)
    bw%Ts=(add-rhs)/coe

    IF (bflag == 0) THEN                                !KJM
       bw%ignited = (bw%Ts > bp%Tignition) 
    ELSE                                                !KJM
       bw%ignited = .TRUE.                              !KJM
    ENDIF                                               !KJM

    bflag = merge(1,0,bw%ignited)
    if(bw%ignited) THEN
       CALL TFUN(bp,bw)
    else
!!       bw%Tstar = bw%Ts
       bw%Tstar = bp%Tstar0  !temporary fix, it will have to be changed if NS simulations
    endif

    Tnp1(1)=bw%Ts

! EVALUATE burning rate

    CALL MFUN(bp,bw)
    IF (bflag == 0) THEN                                !KJM
       bw%rb = merge(bw%rb, zero, bw%Ts >= bp%Tignition)
    ELSE
       bw%rb = bw%rb
    ENDIF
!
!    CONVERT the output back to MKS
!
    rb = bw%rb / m2cm
    Tflame = bw%Tstar

!
!  CHECK the output for divergence
!
    IF (.TRUE. .AND. (.NOT. (rb > bp%rb_range(1) .AND. rb < bp%rb_range(2)) )) THEN
       write(*,*)"OUTPUT RB OUT OF RANGE",rb,bp%rb_range
       do i = 1,nx
	  write(*,*)i,bp%x(i),Tnp1(i),Tn(i)
       enddo
       CALL MPI_Abort (bp%comm, 1, ierror)
       STOP
    ENDIF
    IF (.TRUE. .AND. (.NOT. (Tflame > bp%Tf_range(1) .AND. Tflame < bp%Tf_range(2)) )) THEN
       write(*,*)"OUTPUT TFLAME OUT OF RANGE",Tflame,bp%Tf_range,bflag,bw%Ts,P
       do i = 1,nx
	  write(*,*)i,bp%x(i),Tnp1(i),Tn(i)
       enddo
       CALL MPI_Abort (bp%comm, 1, ierror)
       STOP
    ENDIF


!--------------------------------------------------------------------------------
    RETURN   
  END SUBROUTINE BURN_GET_BURNING_RATE1D
!****************************************************************************


!*****************************************************************************
  SUBROUTINE GFUN(bp,bw)
   
    IMPLICIT NONE

!--------------------------------------------------------------------------------
!   Dummy variables:  
    TYPE(parameter_structure),POINTER   :: bp
    TYPE(work_structure), INTENT(INOUT) :: bw
   
!   Local Variables
    REAL(DBL) :: expterm,invtsq,Tstar,tmp1,tmp2,dy  !Tstar is not passed back
    INTEGER   :: jj,kk
!---------------------------------------------------------------------------------

    bw%c1 = -log(bp%a_p/bp%Ac*(bw%P/bp%Pref)**bp%n_p)/bp%ec_ru   !a_p==capk
    bw%c2 = one / bp%Tstar0
    bw%c3 = two * bp%ec_ru / bp%eg_ru
    
    IF(bw%ignited) THEN

       if(bp%TABUSE == 0) then
          bw%c4 = bp%Tstar0 - bw%To
          bw%c5 = bp%Ac / bp%alfac  
          bw%c6 = - bp%ec_ru

          bw%Ts0     = one/bw%c1
          bw%Tslimit = one/( bw%c1 - bw%c2/bw%c3 ) * 0.99  !99% of blow off temp.

          invtsq  = one / bw%Ts**2
          expterm = bw%c5 * exp(bw%c6 / bw%Ts)
          Tstar   = one/(bw%c2-bw%c3*(bw%c1-one/bw%Ts))

          if ( bw%Ts > bw%Ts0 .AND. bw%Ts > bw%Tslimit ) Tstar = bp%Tstar0
          bw%fs      = expterm * ( bw%c4 + bw%Ts-Tstar )
          bw%fsprime = bw%c6 * bw%fs * invtsq &
               + expterm * (one - bw%c3 * Tstar**2 * invtsq)

       else

          jj = 0;kk = 0;
          call polin2(bp,bp%TABLE%tsurf00,bp%TABLE%press00,bp%TABLE%heatflux00,&
               bp%TABLE%nx_table,bp%TABLE%ny_table,bw%Ts,bw%p,bw%fs,jj,kk)
          call polin2(bp,bp%TABLE%tsurf00, bp%TABLE%press00, bp%TABLE%heatflux00,&
               bp%TABLE%nx_table, bp%TABLE%ny_table, bw%Ts+half*bp%TABLE&
               %small_1, bw%p, tmp2,jj,kk)
          call polin2(bp,bp%TABLE%tsurf00, bp%TABLE%press00, bp%TABLE%heatflux00,&
               bp%TABLE%nx_table, bp%TABLE%ny_table, bw%Ts-half*bp%TABLE&
               %small_1, bw%p, tmp1,jj,kk)
          bw%fsprime = (tmp2-tmp1) / bp%TABLE%small_1
       endif
      
    ELSE

       bw%fs = bw%Qc/bw%lamc
       bw%fsprime = bw%Qcprime / bw%lamc

    ENDIF

!------------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE GFUN
!*****************************************************************************

!*****************************************************************************
  SUBROUTINE MFUN(bp,bw)

!--------------------------------------------------------------------------------------
    IMPLICIT NONE
!   Dummy variables:  
    TYPE(parameter_structure),POINTER   :: bp
    TYPE(work_structure), INTENT(INOUT) :: bw
    INTEGER :: jj,kk
    REAL(DBL) :: dy
!--------------------------------------------------------------------------------------


!   Mass Flux over Density

    if(bp%TABUSE == 0) then
!!=====================================================================================
!! Below are the changes made in order to make the model burn according to aP^n
!! instead of according to the temperature, which appears no to work. Also, TFUN
!! was changed to inject only at the adiabatic flame temperature
!!=====================================================================================
!!       bw%rb = bp%Ac * exp ( - bp%ec_ru / bw%Ts)
       bw%rb = bp%a_p * (bw%P/bp%Pref)**bp%n_p
!!======================================================================================
       bw%alfa_eff = bp%alfac
    else
       jj = 0;kk = 0;
       call polin2(bp,bp%TABLE%tsurf00,bp%TABLE%press00,bp%TABLE%rb00,  &
            bp%TABLE%nx_TABLE,bp%TABLE%ny_table,bw%Ts,bw%p,bw%rb,jj,kk)
       call polin2(bp,bp%TABLE%tsurf00,bp%TABLE%press00,bp%TABLE%alph00, &
            bp%TABLE%nx_TABLE,bp%TABLE%ny_table,bw%Ts,bw%p,bw%alfa_eff,jj,kk)
    endif

!---------------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE MFUN
!*****************************************************************************


!*****************************************************************************
    SUBROUTINE TFUN(bp,bw)
      IMPLICIT NONE

!--------------------------------------------------------------------------------------
!   Dummy variables:
    TYPE(parameter_structure),POINTER   :: bp
    TYPE(work_structure), INTENT(INOUT) :: bw
    INTEGER :: jj,kk
!-------------------------------------------------------------------------------------


!   FLAME TEMPERATURE, Temperature of the injected gas

    if (bp%TABUSE == 0) then
!!=====================================================================================
!! Below are the changes made in order to make the model burn according to aP^n
!! instead of according to the temperature, which appears no to work. Also, TFUN
!! was changed to inject only at the adiabatic flame temperature
!!=====================================================================================
!!       bw%Tstar = one/ ( bw%c2 - bw%c3 * (bw%c1-one/bw%Ts) )
       bw%Tstar = bp%Tstar0
!!=====================================================================================
       if ( bw%Ts > bw%Ts0 .AND. bw%Ts > bw%Tslimit ) bw%Tstar = bp%Tstar0
    else
       jj =0;kk = 0;
       call polin2(bp,bp%TABLE%tsurf00,bp%TABLE%press00,bp%TABLE%Tgas00,&
            bp%TABLE%nx_table,bp%TABLE%ny_table,bw%Ts,bw%p,bw%Tstar,jj,kk)
    endif
    

!-----------------------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE TFUN
!*****************************************************************************

END MODULE  UPDATE_PY
