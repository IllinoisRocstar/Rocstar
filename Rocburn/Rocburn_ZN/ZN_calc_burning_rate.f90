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
! ---------------------------------------------------------------------------
!
!   SUBROUTINE  : ZN_calc_burning_rate
!
!   purpose     : calculate burning rate
!
!   Authors          :  J. Weber, K.C. Tang, and M. Q. Brewster
!
!   Creation Date    :  Sep. 04, 2002
!
!   Modifications    :
!
!    No.     Date         Authors       Description
!
!
! ---------------------------------------------------------------------------
!
!   arguments   :
!
!
!      G_ZN     : Global variables for Rocburn_1D_ZN
!      P        : pressure (atm)
!      qr       : radiant flux (cal/cm^2/s)
!      To       : initial temperature (K)
!      rhoc     : condensed phase density (g/cm^3)
!      qr_old   : radiant flux at previous time step (cal/cm^2/s)
!      fr_old   : fr at previous time step
!      Toa      : apparent initial temperature (K)
!      rb       : burning rate (cm/s)
!      Ts       : surface temperature (K)
!      fr       : fraction of radiation absorbed below surface reaction
!                 zone
!      Tn       : temperature profile at previous time step
!      Tnp1     : temperature profile
!
! ---------------------------------------------------------------------------
!
  SUBROUTINE ZN_calc_burning_rate (G_ZN, delt, P, qr, To, rhoc, qr_old, fr_old,  &
                               Toa, rb, Ts, fr, Tn, Tnp1)

    USE M_Rocburn_ZN_Global_Data
  
    IMPLICIT NONE

    INCLUDE 'mpif.h'

    TYPE(G_BURN_1D), POINTER :: G_ZN
    REAL(DBL), INTENT(IN)    :: delt, P, qr, To, rhoc
    REAL(DBL), INTENT(IN)    :: qr_old, fr_old
    REAL(DBL), INTENT(INOUT) :: Toa
    REAL(DBL), INTENT(INOUT) :: rb, Ts, fr
    REAL(DBL), INTENT(IN)    :: Tn(:)
    REAL(DBL), INTENT(OUT)   :: Tnp1(:)

!
! ============================================================================
!
!     delcare local variables
!
    INTEGER   :: iter, itrmx, i, MPI_COMM_ROCBURN
    REAL(DBL) :: rblast, Tslast
    REAL(DBL) :: fs, rb_old
    REAL(DBL) :: first(1:G_ZN%nxmax), second(1:G_ZN%nxmax)
    REAL(DBL) :: da(1:G_ZN%nxmax), db(1:G_ZN%nxmax), dc(1:G_ZN%nxmax)
    INTEGER   :: ierror

!
! ============================================================================
!

    IF(G_ZN%Model_combustion == 3) THEN
!
!        Model_combustion = 3
!
!        rb=a*P**n
!

         rb=G_ZN%a_p*P**G_ZN%n_p
         Tnp1(1)=1000.0
   
         RETURN

    END IF

    rb_old = rb

    fs=rb*(Ts-Toa)/G_ZN%alfac-fr_old*qr_old/G_ZN%lamc

    DO iter=0,G_ZN%itermax

!
!      get 1st and 2nd derivatives at last time step (n)
!      with this iterations surface temperature B.C.

       CALL cmpct1(G_ZN%nx,G_ZN%delz,Tn,first, da, db, dc)
       CALL cmpct2(G_ZN%nx,G_ZN%delz,Tn,second, da, db, dc)

!
!      determine T profile at this time step (n+1) based on
!      properties a last time step (n) (explicit differencing)
!
       DO i=2,G_ZN%nx-1
          Tnp1(i)=Tn(i)+delt*(                                 &
                 (G_ZN%alfac*G_ZN%zxx(i)-rb_old*G_ZN%zx(i))*first(i) +        &
                  G_ZN%alfac*G_ZN%zx(i)**2*second(i) +                   &
                  fr_old*qr_old/rhoc/G_ZN%C*G_ZN%Ka*exp(G_ZN%Ka*G_ZN%x(i)))  
       END DO ! i

       Tnp1(G_ZN%nx)=G_ZN%To
  

!        update suface B.C. for this iteration
!        Tnp1(1)=(48.0*Tnp1(2)-36.0*Tnp1(3)+16.0*Tnp1(4)                 ! 4th
!    +            -3.0*Tnp1(5)-12.0*G_ZN%delz*fs/G_ZN%zx(1))/25.0        ! 4th
         Tnp1(1)=(18.0*Tnp1(2)-9.0*Tnp1(3)+2.0*Tnp1(4)                   &
                 -6.0*G_ZN%delz*fs/G_ZN%zx(1))/11.0                      ! 3rd
!        Tnp1(1)=(Tnp1(2)-G_ZN%delz*fs/G_ZN%zx(1))                       ! 1st

   

         Ts=Tnp1(1)
!        Tn(1)=Tnp1(1)

!
!        calculate rb, Toa, and fr for given  P, qr, and Ts
!
         CALL ZN(P, qr, To, Ts, Toa, rb, fr)

!
!        calculate fs for this iteration of (n+1)
!
         fs=rb*(Ts-Toa)/G_ZN%alfac-fr*qr/G_ZN%lamc

         IF(iter.ge.1) THEN
            IF(abs(Ts-Tslast).lt.G_ZN%tol_Ts.and.abs(rb-rblast).lt.G_ZN%tol_Ts)  THEN
               RETURN
            END IF
         END IF

         Tslast=Ts
         rblast=rb

    END DO ! iter
        
    PRINT *,' ROCBURN_ZN: rank=',G_ZN%rank
    PRINT *,' ROCBURN_ZN: Error - burning rate: Maximum # of iterations '
    PRINT *,' ROCBURN_ZN: itrmx=', G_ZN%itermax, 'reached'
    PRINT *,' ROCBURN_ZN: tol_Ts=',G_ZN%tol_Ts
    PRINT *,' ROCBURN_ZN: Tn=',Tn
    PRINT *,' ROCBURN_ZN: Tnp1=',Tnp1
    STOP


    CONTAINS

!     ------------------------------------------------------------------------
!                             INTERNAL PROCEDURES
!     ------------------------------------------------------------------------
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                             c
!   subroutine  : ZN                                                          c
!                                                                             c
!   purpose     :  alculate rb,Toa given Ts,P,qr                              c
!                                                                             c
! --------------------------------------------------------------------------- c
!                                                                             c
!                            written by                                       c
!                      K.C. Tang and M.Q. Brewster                            c
!                                                                             c
!                       Last modified : 09/14/2000                            c
! --------------------------------------------------------------------------- c
!                                                                             c
!   input       :                                                             c
!                                                                             c
!      P        : pressure (atm)                                              c
!      qr       : radiative heat flux (al/cm^2-s)                             c
!      To       : propellant temperatuer at deep into the propellant (K)      c
!      Ts       : surface temperature (K)                                     c
!      Toa      : ZN temperature (K) for initial guess                        c
!      rb       : burning rate (cm/s) for initial guess                       c
!                                                                             c
!                                                                             c
!   output      :                                                             c
!                                                                             c
!      Toa      : ZN temperature (K)                                          c
!      rb       : burning rate (cm/s)                                         c
!      fr       : transmissivity of surface reaction zone                     c
!                                                                             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!

!
      SUBROUTINE ZN(P, qr, To, Ts, Toa, rb, fr)

      IMPLICIT NONE

      REAL(DBL) :: P, qr, To
      REAL(DBL) :: Ts, rb, Toa, fr
      REAL(DBL) :: Tf
      REAL(DBL) :: delta_T2, Ag


      REAL(DBL), parameter :: tol=1.d-5

! ============================================================================
!
!     delcare local variables
!

      REAL(DBL) :: J(2,2),xx(2),e(2)
      REAL(DBL) :: xreac
      REAL(DBL) :: term0, term1, term2, term3, term4, term5
      REAL(DBL) :: xcd, Da, xg

      integer iter


      term1=G_ZN%Ac*G_ZN%R*(Ts*rhoc)**2*G_ZN%alfac*G_ZN%C* &
            exp(-G_ZN%Ec/G_ZN%R/Ts)/G_ZN%Ec


!
!     initial guess
!

      xx(1)=rb*rhoc
      xx(2)=Toa

      IF(G_ZN% Model_combustion == 1) THEN

!        Model_combustion = 1
!
!        WSB homogeneouse propellant combustion model
!

         Ag=G_ZN%Bg/(G_ZN%R*G_ZN%R)*G_ZN%MW*G_ZN%MW


         DO iter=0,G_ZN%itermax

!
!           set up residual matrix and Jacobian derivative matrix for
!           Ibiricu and Williams condensed phase expression

            xreac=2.0*rhoc*G_ZN%alfac*G_ZN%R*Ts/xx(1)/G_ZN%Ec
            term0=1.0-G_ZN%Ka*xreac

            fr=dexp(-xreac*G_ZN%Ka)
            term2=G_ZN%C*(Ts-xx(2))-G_ZN%Qc/2.0-fr*qr/xx(1)
            term3=term1/term2

            e(1)=xx(1)**2-term3
            J(1,1)=2.0*xx(1)+term3/term2*fr*qr/xx(1)**2*term0
            J(1,2)=-term3*G_ZN%C/term2

!
!           set up residual matrix and Jacobian derivative matrix for
!           Ward-Son-Brewster gas phase expression


            xcd=G_ZN%lamg/G_ZN%C/xx(1)
            Da=4.0*G_ZN%Bg*G_ZN%C/G_ZN%lamg*(P*G_ZN%MW*xcd/G_ZN%R)**2
            term4=(1.0+Da)**0.5
            xg=2.0*xcd/(term4-1.0)
            term5=G_ZN%Qc-G_ZN%C*(Ts-xx(2))
            Tf=(qr/xx(1)+G_ZN%Qc+G_ZN%Qg)/G_ZN%C+xx(2)

            e(2)=(Tf-Ts)/xg+(qr+xx(1)*term5)/G_ZN%lamg
            J(2,1)=-qr/xg/xx(1)**2/G_ZN%C-(Tf-Ts)/xg/xx(1)/term4+term5/G_ZN%lamg
            J(2,2)=1/xg+1/xcd

!
!           update guesses for m and Toa using Newton-Rhapson
!           & guard against a negative (erroneus) rb solution
!     
            call Newton(xx,J,e)
            if(xx(1).lt.0.0) then
               PRINT *,' ROCBURN_ZN: rank=',G_ZN%rank,   &
                       ' rb negative xx(1)=',xx(1),' Ts=',Ts
               xx(1)=abs(xx(1))
            end if

            if(abs(e(1)).lt.tol .and. abs(e(2)).lt.tol) then

!
!              reset final converged values
!
               rb=xx(1)/rhoc
               Toa=xx(2)
               fr=dexp(-xreac*G_ZN%Ka)
               return

            end if
         END DO ! iter

         WRITE(*,*) ' ROCBURN_ZN: in ZN rank=',G_ZN%rank,  &
                    ' Error-WSB: Maximum # of iterations reached'
         STOP

      end if

      if(G_ZN%Model_combustion.eq.2) then

!        Model_combustion = 2
!
!        Zeldovich-Novozhilov phenomenological model for composite
!        propellant
!

         DO iter=0,G_ZN%itermax

!
!           set up residual matrix and Jacobian derivative matrix for
!           Ibiricu and Williams condensed phase expression

            xreac=2.0*rhoc*G_ZN%alfac*G_ZN%R*Ts/xx(1)/G_ZN%Ec
            term0=1.0-G_ZN%Ka*xreac

            fr=dexp(-xreac*G_ZN%Ka)
            term2=G_ZN%C*(Ts-xx(2))-G_ZN%Qc/2.0-fr*qr/xx(1)
            term3=term1/term2

            e(1)=xx(1)**2-term3
            J(1,1)=2.0*xx(1)+term3/term2*fr*qr/xx(1)**2*term0
            J(1,2)=-term3*G_ZN%C/term2

!
!           set up residual matrix and Jacobian derivative matrix for
!           Ward-Son-Brewster gas phase expression
!


            delta_T2 = Ts - 300.0
            e(2)=G_ZN%a_T*((Ts-xx(2)-qr/(xx(1)*G_ZN%C))/delta_T2)**G_ZN%n_T -    &
                 rhoc*G_ZN%a_p*p**G_ZN%n_p/xx(1)
            J(2,1)=rhoc*G_ZN%a_p*p**G_ZN%n_p/(xx(1)*xx(1))
            J(2,2)=-G_ZN%n_T*G_ZN%a_T*(Ts-xx(2)-qr/(xx(1)*G_ZN%C))**(G_ZN%n_T-1.0)/   &
                    (delta_T2**G_ZN%n_T)
         

!
!           update guesses for m and Toa using Newton-Rhapson
!           & guard against a negative (erroneus) rb solution
!    

            call Newton(xx,J,e)

            if(xx(1).lt.0.0) then
               PRINT *,' ROCBURN_ZN: inside ZN rank=',G_ZN%rank,  &
                       ' rb negative xx(1)=',xx(1),' Ts=',Ts
               xx(1)=abs(xx(1))
            end if

            if(abs(e(1)).lt.tol .and. abs(e(2)).lt.tol) then

!
!              reset final converged values
!
               rb=xx(1)/rhoc
               Toa=xx(2)
               fr=dexp(-xreac*G_ZN%Ka)
               RETURN

            END IF
         END DO ! iter

         WRITE(*,*) ' ROCBURN_ZN: in ZN rank=',G_ZN%rank,    &
                    ' Error-Composite: Maximum # of iterations reached'
         STOP

      end if

      IF(G_ZN%Model_combustion == 3) THEN

!
!
!        Model_combustion = 3
!
!        rb=a*P**n
!

         rb=G_ZN%a_p*P**G_ZN%n_p
         Tnp1(1)=1000.0
   
         RETURN

      END IF

      WRITE(*,*) ' ROCBURN_ZN: in ZN rank=',G_ZN%rank,  &
                 ' Error-ZN: No appropriate combustion model'
      WRITE(*,*) ' ROCBURN_ZN: Combustion Model= ', G_ZN%Model_combustion
      STOP

      END SUBROUTINE ZN



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                             c
!   subroutine  : Newton                                                      c
!                                                                             c
!                                                                             c
!   Author:                                                                   c
!                                                                             c
!      Paul Loner                                                             c
!                                                                             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                             c
      SUBROUTINE Newton(x,J,e)

      implicit none

      REAL(DBL) :: x(2),J(2,2),e(2),detJ,delx(2),dum
      integer i

      detJ=J(1,1)*J(2,2)-J(2,1)*J(1,2)

      IF(detJ.eq.0) then
        PRINT *,' ROCBURN_ZN: in Newton rank=',G_ZN%rank,   &
                ' Error: singular matrix encountered'
        PRINT *,' ROCBURN_ZN: P=', P
        PRINT *,' ROCBURN_ZN: J=', J
        STOP
      ENDIF
!
!     invert the 2x2 Jacobian matrix
!
      dum=J(1,1)
      J(1,1)=J(2,2)
      J(2,2)=dum
      J(1,1)=J(1,1)/detJ
      J(1,2)=-J(1,2)/detJ
      J(2,1)=-J(2,1)/detJ
      J(2,2)=J(2,2)/detJ
!
!     multiply J^-1()*e() and update guesses
!
      do i=1,2
        delx(i)=J(i,1)*e(1)+J(i,2)*e(2)
        x(i)=x(i)-delx(i)
      end do ! i

      RETURN
      END SUBROUTINE Newton

!***********************************************************************
!
      SUBROUTINE cmpct1(imax,h,u,f, a, b, cc)
!
!***********************************************************************
!  Modifications:
!
!    No.     Date         Programmer    Description
!    001  May  03, 2001   K. C. Tang    change memory allocation method for
!                                       temporary arrarys inside a 
!                                       subroutine

      IMPLICIT NONE

      integer i,imax
      REAL(DBL) :: h,u(imax),f(imax)

      REAL(DBL) :: a(:),b(:),cc(:)


!     set up equations for first derivatives
!      print *,' in cmpct1 imax=',imax,' h=',h,' u=',u,' f=',f

!     ALLOCATE(a(imax))
!     ALLOCATE(b(imax))
!     ALLOCATE(cc(imax))

      a(1)=0.0
      b(1)=1.0
      cc(1)=0.0
      f(1)=(-25.0*u(1)+48.0*u(2)-36.0*u(3)+16.0*u(4)-3.0*u(5))/(12.0*h) !4th

      do i=2,imax-1
        a(i)=1.0
        b(i)=4.0
        cc(i)=1.0
        f(i)=3.0*(u(i+1)-u(i-1))/h
      end do ! i

      a(imax)=0.0
      b(imax)=1.0
      cc(imax)=0.0
      f(imax)=0.0

!     solve for first derivatives

      call tridg (imax,a,b,cc,f)

!     DEALLOCATE(a)
!     DEALLOCATE(b)
!     DEALLOCATE(cc)
      RETURN
      END SUBROUTINE cmpct1


!***********************************************************************
!
      SUBROUTINE cmpct2(imax,h,u,s, a, b, cc )
!
!***********************************************************************
!  Modifications:
!
!    No.     Date         Programmer    Description
!    001  May  03, 2001   K. C. Tang    change memory allocation method for
!                                       temporary arrarys inside a 
!                                       subroutine

      implicit none

      integer i,imax
      REAL(DBL) :: h,h2,u(imax),s(imax)
      REAL(DBL) :: a(:),b(:),cc(:)
!     REAL(DBL), ALLOCATABLE :: a(:),b(:),cc(:)


      h2=h*h

!     ALLOCATE(a(imax))
!     ALLOCATE(b(imax))
!     ALLOCATE(cc(imax))

!     set up equations for second derivatives

      a(1)=0.0
      b(1)=1.0
      cc(1)=0.0
      s(1)=(45.0*u(1)-154.0*u(2)+214.0*u(3)-156.0*u(4)+61.0*u(5)-     &
           10.0*u(6))/12.0/h2                                        ! 4th

      DO i=2,imax-1
        a(i)=1.0
        b(i)=10.0
        cc(i)=1.0
        s(i)=12.0*(u(i+1)-2.0*u(i)+u(i-1))/h2
      END DO ! i

      a(imax)=0.0
      b(imax)=1.0
      cc(imax)=0.0
      s(imax)=0.0

!     solve for second derivatives

      call tridg (imax,a,b,cc,s)

!     DEALLOCATE(a)
!     DEALLOCATE(b)
!     DEALLOCATE(cc)

      RETURN
      END SUBROUTINE cmpct2


!***********************************************************************
!
      SUBROUTINE tridg (imax,a,b,c,d)
!
!     solves a tridiagonal matrix equation of the form:
!
!       a(i)*u(i-1) + b(i)*u(i) + c(i)*u(i+1) = d(i)
!
!     using the Thomas algorithm.
!
!***********************************************************************

      IMPLICIT NONE

      INTEGER :: imax
      REAL(DBL)  :: a(:), b(:), c(:), d(:)
!     REAL(DBL), DIMENSION(imax) :: a, b, c, d
 
      INTEGER ::  i

      DO i=2,imax
        a(i)=a(i)/b(i-1)
        b(i)=b(i)-a(i)*c(i-1)
        d(i)=d(i)-a(i)*d(i-1)
      END DO ! i

      d(imax)=d(imax)/b(imax)

      DO     i=imax-1,1,-1
        d(i)=(d(i)-c(i)*d(i+1))/b(i)
      END DO ! i
      RETURN
      END SUBROUTINE tridg


!     -------------------------------------------------------------------
!                       END OF INTERNAL PROCEDURES
!     -------------------------------------------------------------------


  END SUBROUTINE ZN_calc_burning_rate






