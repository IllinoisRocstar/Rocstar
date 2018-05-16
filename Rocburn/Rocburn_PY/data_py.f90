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
!**************************************************************************
MODULE data_py

  IMPLICIT NONE

 INTEGER, PARAMETER :: DBL = SELECTED_REAL_KIND(P=14,R=30)
!=========================================================================

 REAL *8, PARAMETER ::  zero    = 0.d0,   &
                        one     = 1.0d0,  &
                        two     = 2.0d0,  &
                        three   = 3.0d0,  &
                        four    = 4.0d0,  &
                        eight   = 8.0d0,  &

                        half    = 0.5d0,        &
                        quarter = 0.25d0,       &
                        eighth  = 0.125d0,      &
                        third   = one/three,         &
                        two_thirds = two/three,      &
                        root_two = 1.41421356237309515d0,   &

                        j2cal = 1.0d0/4.1868d0,         &
                        kg2g  = 1000.0d0,               &
                        m2cm  = 100.0d0,                &
                        ru    = 1.9859d0,               &
                        mpa2atm = 9.86923266716d0,      &
                        pa2atm = mpa2atm*1d-6,          &     
                        kgmc2gcc = kg2g/m2cm**3,        &
                        j_m2cal_cm = j2cal/m2cm,        &
                        j_msq2cal_cmsq = j2cal/m2cm**2, &
                        gcc2kgmc = one/kgmc2gcc,        &
                        kcalmc2atm = 4186.8/101325.d0,  &
                        j_kg2cal_g = j2cal/kg2g
!=======================================================================
! global structure that contains combustion model parameters

  TYPE parameter_structure

    REAL(DBL) :: a_p,n_p,Pref,Ac,eg_ru,ec_ru,alfac,C,lamg
    REAL(DBL) :: Tstar0, To, Tignition, Tsurf, film_cons, lamc
    
    INTEGER   :: comm, igrid, numx, ixsymm
    INTEGER   :: TABUSE
    REAL(DBL) :: delt, xmax, beta, delz, x_surf_burn
    REAL(DBL) :: delt_max, delz2inv, delzsqinv, dx1
    REAL(DBL),POINTER :: x(:), z(:), zx(:), zxx(:)
    REAL(DBL) :: P_range(2), rb_range(2), Tf_range(2) 
    CHARACTER*50 :: TABNAM

    TYPE(TABLETYPE),POINTER :: TABLE
    
  END TYPE parameter_structure

  TYPE work_structure
    
    REAL(DBL) :: P,rhoc,Qc,Qcprime,rb,Ts,To,lamc
    REAL(DBL) :: fs,fsprime,Tstar,Ts0,Tslimit
    REAL(DBL) :: c1,c2,c3,c4,c5,c6,alfa_eff,fx2
    LOGICAL   :: ignited

 END TYPE work_structure

  TYPE, PUBLIC :: TABLETYPE
!A-  table dimensions
     INTEGER  ::  nx_table,ny_table,nfield
     LOGICAL  ::  spline
!B-  table fields
     REAL(DBL),POINTER ::  tsurf00(:),press00(:), heatflux00(:,:)
     REAL(DBL),POINTER ::  rb00(:,:),fxsq00(:,:),Tgas00(:,:)
     REAL(DBL),POINTER ::  Tstd00(:),rstd00(:),alph00(:,:)
!C-  Work Arrays
     REAL(DBL),POINTER ::  wrk1(:),wrk2(:),wrk3(:),wrk4(:),wrk5(:)
     REAL(DBL),POINTER ::  wrk6(:),wrk7(:),wrk8(:),wrk9(:),wrk10(:)
!D-  Table parameters
     REAL(DBL) :: alpha,chi,small_1,small_2
  END TYPE TABLETYPE

CONTAINS

!data interpolation subroutines

!*****************************************************************************
  SUBROUTINE polin2(bp,x1a,x2a,y12a,m,n,x1,x2,y,j,k)
    implicit NONE
    TYPE(parameter_structure),POINTER   :: bp
    INTEGER   :: i,l,m,n,ixtrap,jsrt,jend
    INTEGER,INTENT(INOUT) :: j,k
    REAL(DBL) :: dy,x1,x2,y
    REAL(DBL) :: x1a(m),x2a(n),y12a(m,n)
    REAL(DBL) :: prod,difm,difp,a(4),del12
    REAL(DBL), POINTER :: yntmp(:),ymtmp(:)
!-----------------------------------------------

!    yntmp => Bp%TABLE%wrk1
!    ymtmp => Bp%TABLE%wrk2

    if(j+k > 0) goto 123

    j=0
    ixtrap = 0
    prod = 1.0
    do while(prod .gt. 0)
       j=j+1
       difm=x1-x1a(j)
       difp=x1-x1a(j+1)
       prod = difm*difp
       if(j .ge. m-1 .AND. prod .gt. 0)then
          ixtrap = 1
          prod = -1.0   !break
          if(abs(x1-x1a(m)) .lt. abs(x1-x1a(1)) ) then
             j = m-1
          else
             j = 2
          endif
       endif
    enddo
    k=0
    ixtrap = 0
    prod = 1.0
    do while(prod .gt. 0)
       k=k+1
       difm=x2-x2a(k)
       difp=x2-x2a(k+1)
       prod = difm*difp
       if(k .ge. n-1 .AND. prod .gt. 0)then
          ixtrap = 1
          prod = -1.0   !break
          if(abs(x2-x2a(n)) .lt. abs(x2-x2a(1)) ) then
             k = n-1
          else
             k = 2
          endif
       endif
    enddo

123 CONTINUE

    del12 = (x2a(k+1) - x2a(k))*(x1a(j+1) - x1a(j))

    a(1) = (x2a(k+1) - x2) * (x1a(j+1) - x1)
    a(2) = (x2a(k+1) - x2) * (x1 - x1a(j))
    a(3) = (x2 - x2a(k)) * (x1a(j+1) - x1)
    a(4) = (x2 - x2a(k)) * (x1 - x1a(j))

    a = a/del12

!use bilinear interpolation
    y = a(1)*y12a(j,k) + a(2)*y12a(j+1,k) + a(3) * y12a(j,k+1) + a(4) * y12a(j+1,k+1)

!---------------
    RETURN
  END SUBROUTINE polin2
!*****************************************************************************



!*****************************************************************************
!     Spline interpoaltion routine. 
    SUBROUTINE polint(bp,xa,ya,n,x,y,dx)
      implicit none
      integer i,j,k,l,m,n,ixtrap
      TYPE(parameter_structure),POINTER   :: bp
      REAL(DBL) :: x,y,prod,difm,difp,yo,der,dx
      REAL(DBL) :: xa(n),ya(n)
      REAL(DBL), POINTER :: h(:),alp(:)
      REAL(DBL), POINTER :: c1(:),c2(:),c3(:)
      REAL(DBL), POINTER :: v1(:),v2(:),v3(:)
!-----------------------------------------------


      h   => Bp%TABLE%wrk3

      prod = 1.0
      i=0
      do i =1,n-1
         h(i) = xa(i+1) - xa(i)
      enddo

!SPLINE
      if(bp%TABLE%spline) then
         alp => Bp%TABLE%wrk4
         c1  => Bp%TABLE%wrk5
         c2  => Bp%TABLE%wrk6
         c3  => Bp%TABLE%wrk7
         v1  => Bp%TABLE%wrk8
         v2  => Bp%TABLE%wrk9 
         v3  => Bp%TABLE%wrk10
         do i =2,n-1
            j = i-1
            alp(i) = 3.0/h(i)*(ya(i+1) -ya(i)) &
                 - 3.0/h(j)*(ya(j+1) -ya(j))
         enddo
         v1(1) = 1.d0
         v2(1) = 0.d0
         v3(1) = 0.d0
         do i = 2,n-1
            v1(i) = 2.0d0*(xa(i+1) -xa(i-1)) - h(i-1)*v2(i-1)
            v2(i) = h(i)/v1(i)
            v3(i) = ( alp(i) - h(i-1)*v3(i-1) ) / v1(i)
         enddo

         v1(n) = 1.d0
         v3(n) = 0.d0
         c2(n) = 0.d0
         do i = n-1,1,-1
            c2(i) = v3(i) - v2(i)* c2(i+1)
            c1(i) = (ya(i+1) -ya(i))/h(i) - &
                 h(i) * ( c2(i+1)+2.0*c2(i) ) / 3.0d0
            c3(i) = (c2(i+1)-c2(i))/(3.0d0*h(i))
         enddo
      endif

      i=0
      ixtrap = 0
      prod = 1.0
      do while(prod .gt. 0)
         i=i+1
         difm=x-xa(i)
         difp=x-xa(i+1)
         prod = difm*difp
         if(i .ge. n-1 .AND. prod .gt. 0)then
!            write(*,101)'FAILED SEARCH',x,(xa(m),m=1,n)
            ixtrap = 1
            prod = -1.0
         endif
      enddo
!101   format(2x,a,1p100e10.2)


      if(bp%TABLE%spline) then
!     SPLINE 
         y =  YA(i) + C1(i)*(x-xa(i)) +  C2(i)*(x-xa(i))**2&
              + C3(i)*(x-xa(i))**3
      else
!     USE first order interpolation to avoid problem w\ oscillatory solutions
         y = YA(i) + (ya(i+1)-ya(i))/h(i) * (x-xa(i))
      endif

      if(ixtrap .eq. 1) then
         if(abs(x-xa(n)) .lt. abs(x-xa(1)) ) then
            yo = ya(n)
            der = (ya(n)-ya(n-1))/h(n-1)
            dx = x-xa(n)
         else
            yo = ya(1)
            der = (ya(2)-ya(1))/h(1)
            dx = x-xa(1)    
         endif
         y = yo+der*dx
      endif

      return
    END SUBROUTINE polint
!*************************************************


!=======================================================================

END MODULE data_py
!************************************************************************






