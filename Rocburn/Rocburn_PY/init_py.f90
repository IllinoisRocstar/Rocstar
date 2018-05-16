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
MODULE  INIT_PY
  USE data_py

CONTAINS


!****************************************************************
  SUBROUTINE burn_get_film_coeff_1d(bp,p_coor,Ts,Teuler,P_in,Qc,Qcprime)
   
    IMPLICIT NONE
!---------------------------------------------------------------
!   DUMMY VARIABLES
    TYPE(parameter_structure),POINTER  :: bp
    REAL(DBL), INTENT(IN)  :: p_coor(3),Ts,Teuler,P_in
    REAL(DBL), INTENT(OUT) :: Qc,Qcprime
!   LOCAL VARIABLES
    REAL(DBL)   :: P, Te, try_film, front_dist
    REAL(DBL)   :: cp_MKS, exp_moser, coe_moser, coe1, par
    REAL(DBL)   :: Prndtl, lambda, mu, x_surf
!----------------------------------------------------------------

    P = P_in                         !MKS
!
!   Avoid Euler temperatures larger than 2.0*(adiabatic flame temperature)
!   this is a temporary hardwired fix, it should never affect the solution
!
    Te = min(Teuler,2.0d0*bp%Tstar0)

    IF(bp%ixsymm >= 1)THEN
       x_surf = p_coor(bp%ixsymm)
       Prndtl = 0.72d0
       cp_MKS = bp%C/j_kg2cal_g
       lambda = (2.581d-7*Te + 3.1788d-5)*4186.8/10.0   !Buckmaster
       mu = Prndtl*lambda

       exp_moser = - 1.0d0/5.0d0
       coe_moser = 800.0d0*14.0d0
       coe1 = 800.0d0/0.0287d0

!
!  for x_surf < bp%x_surf_burn the propellan is all burning
!
       front_dist = (x_surf - bp%x_surf_burn)  
       par = max( (coe1*front_dist*mu*cp_MKS*lambda**(-2.0)), 1.d-6 ) 

       try_film = merge(coe_moser*par**exp_moser, 0.0d0, front_dist >= 0.0d0)

    ELSE

       try_film = bp%film_cons        !MKS

    ENDIF

!
!   Avoid negative heat flux
!
    Qc = max(try_film * (Te - Ts), 0.0d0)     !MKS
    Qcprime = - try_film                      !MKS

!------------------------------------------------------------------------
    RETURN
  END SUBROUTINE burn_get_film_coeff_1d
!*************************************************************************



!*************************************************************************
  SUBROUTINE burn_init_1d(bp,bflag,Pin,To,rhoc,p_coor,rb,Toa,fr,Tn,Tflame)

    IMPLICIT NONE
!-------------------------------------------------------------------------
!   DUMMY VARIABLES
    TYPE(parameter_structure),POINTER  :: bp
    INTEGER, INTENT(INOUT) :: bflag
    REAL(DBL), INTENT(IN)  :: Pin,To,rhoc,p_coor(3)
    REAL(DBL), INTENT(OUT) :: rb,Toa,fr
    REAL(DBL), INTENT(OUT) :: Tn(:)
    REAL(DBL), INTENT (OUT)   :: Tflame
!   LOCAL VARIABLES
    REAL(DBL) :: xcond, c1, Ts, dtemp, P, x_surf, dyTAB, alp
    INTEGER   :: i,jj,kk
!-------------------------------------------------------------------------

!
!   SET BFLAG IF AXISYMMETIC BURNING, OTHERWISE BFLAG COMES FROM FLUIDS
!
    IF(bp%ixsymm >= 1)THEN
       x_surf = p_coor(bp%ixsymm)
       IF(x_surf <= bp%x_surf_burn) THEN
          bflag = 1
       ELSE
          bflag = 0
       ENDIF
!RAF ----------------------------------------------------------------------
    ELSE

!RAF HACK: Turn the initial burning off so we can do problems with igniters.
!RAF HACK: For the lab scale rocket, we must set ixsymm >= 1.

       bflag = 0
!RAF ----------------------------------------------------------------------
    ENDIF
!
!   CONVERT  inputs
!
    P = Pin*pa2atm

!
!   calculate Ts, rb given To, P
!
    if(bflag /= 0 ) THEN

       if(bp%TABUSE == 0) then
          rb    = bp%a_p*( P/bp%Pref )**bp%n_p
          c1    = -log( rb/bp%Ac ) / bp%ec_ru
          Ts    = one/c1
          alp = bp%alfac
       else
          call polint(bp,bp%TABLE%press00,bp%TABLE%Tstd00,bp%TABLE%ny_table,P,Ts,dyTAB)
          call polint(bp,bp%TABLE%press00,bp%TABLE%rstd00,bp%TABLE%ny_table,P,rb,dyTAB)
          jj = 0;kk = 0
          call polin2(bp,bp%TABLE%tsurf00,bp%TABLE%press00,bp%TABLE%alph00, &
               bp%TABLE%nx_TABLE,bp%TABLE%ny_table,Ts,P,alp,jj,kk)
       endif

       xcond = rb / alp
       dtemp = Ts - To
       do  i=1,bp%numx
          Tn(i) = Ts - dtemp * (dexp(bp%x(i)*xcond) - one) & 
               / (dexp(bp%xmax*xcond) - one)
       ENDDO
       Tflame = bp%Tstar0
    else
       dtemp = bp%Tsurf - To 
       do  i=1,bp%numx
          xcond = (bp%x(i) - bp%xmax)/bp%xmax  !-1:0
          Tn(i) = To - dtemp * xcond
       enddo
!set the flame temp (there's no flame at this point) to the surf temp
       tflame = bp%Tsurf 
       rb = 0.0d0
    endif

    Toa = bp%Tsurf
    fr = 0.0d0

    rb = rb / m2cm

!------------------------------------------------------------------------
    RETURN
  END SUBROUTINE burn_init_1d
!***********************************************************************

END MODULE INIT_PY






