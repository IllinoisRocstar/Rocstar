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
!   SUBROUTINE  : ZN_ssWSB
!
!   purpose     : calculate steady state solution
!
!   Authors          :  K.C. Tang and M. Q. Brewster
!
!   Creation Date    :  Sep. 15, 2000
!
!   Modifications    :
!
!    No.     Date         Authors       Description
!     1    09/03/02      K.C. Tang      Global variables G_ZN
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
!      rb       : burning rate (cm/s)                   
!      Ts       : surface temperature (K)
!      Tf       : flame temperature      
!      fr       : fraction of radiation absorbed below surface reaction
!                 zone
!
! ---------------------------------------------------------------------------  
!
  SUBROUTINE ZN_ssWSB(G_ZN, P, qr, To, rhoc, rb, Ts, Tf, fr, Tn)

    USE M_Rocburn_ZN_Global_Data

    IMPLICIT NONE
    INCLUDE 'mpif.h'


    TYPE(G_BURN_1D), POINTER :: G_ZN
    REAL(DBL), INTENT(IN) :: P, qr, To, rhoc
    REAL(DBL), INTENT(OUT) :: rb, Ts, Tf, fr
    REAL(DBL), INTENT(OUT) :: Tn(:)


    REAL(DBL), parameter :: tol=1.0d-5


!
! ============================================================================
!
!   delcare local variables
!
    REAL(DBL) :: J(2,2),xx(2),e(2)
    REAL(DBL) :: xreac
    REAL(DBL) :: term0, term1, term2, term3, term4, term5
    REAL(DBL) :: xcd, Da, xg, xcond, frJ

    INTEGER   :: iter, i
!
! ============================================================================
!

!
!   Generate initial conditions
!
!
!

    IF(G_ZN%Model_combustion == 1) THEN

!     Model_combustion = 1
!
!     steady state WSB homogeneouse propellant combustion model
!
!     guess initial values for doublebase propellant 
!     rb and Ts 
!
       
      rb=G_ZN%a_p*(P**G_ZN%n_p)                          ! cm/s
      Ts=550.0                                           ! K
      xx(1)=rb*rhoc                                      ! g/(cm^2*s)
      xx(2)=Ts                                           ! K

      DO iter=0,G_ZN%itermax
!
!       set up residual matrix and Jacobian derivative matrix for
!       Ibiricu and Williams condensed phase expression
!   

        xreac=2.0*rhoc*G_ZN%alfac*G_ZN%R*xx(2)/xx(1)/G_ZN%Ec           ! cm
        term0=1.0-G_ZN%Ka*xreac                                        ! -

        term1=G_ZN%Ac*G_ZN%R*(xx(2)*rhoc)**2*G_ZN%alfac*G_ZN%C*  &
                     dexp(-G_ZN%Ec/G_ZN%R/xx(2))/G_ZN%Ec               ! cal*g/(cm^4*s^2)
 
        fr=dexp(-xreac*G_ZN%Ka)
        term2=G_ZN%C*(xx(2)-To)-G_ZN%Qc/2.0-fr*qr/xx(1)                ! cal/g
        term3=term1/term2                                              ! g^2/(cm^4*s^2)

        e(1)=xx(1)**2-term3                                            ! g^2/(cm^4*s^2)
        J(1,1)=2.0*xx(1)+term3/term2*fr*qr/xx(1)**2                   &
                         *term0                                        ! g/(cm^2*s)
        J(1,2)=-xx(1)**2*((2.0d0+G_ZN%Ec/G_ZN%R/xx(2))/xx(2) -        &
                (G_ZN%C+G_ZN%Ka*xreac*qr*fr/xx(1)/xx(2))/term2)        ! g^2/(cm^4*s^2*K)

!
!       set up residual matrix and Jacobian derivative matrix for
!       Ward-Son-Brewster gas phase expression

        xcd=G_ZN%lamg/G_ZN%C/xx(1)                                     ! cm
        Da=4.0*G_ZN%Bg*G_ZN%C/G_ZN%lamg*(P*G_ZN%MW*xcd/G_ZN%R)**2      ! -
        term4=(1.0+Da)**0.5                                            ! -
        xg=2.0*xcd/(term4-1.0)                                         ! cm
        term5=G_ZN%Qc-G_ZN%C*(xx(2)-To)                                ! cal/g
        Tf=(qr/xx(1)+G_ZN%Qc+G_ZN%Qg)/G_ZN%C+To                        ! K

        e(2)=(Tf-xx(2))/xg+(qr+xx(1)*term5)/G_ZN%lamg                  ! K/cm
        J(2,1)=-qr/xg/xx(1)**2/G_ZN%C-                          &
                (Tf-xx(2))/xg/xx(1)/term4+term5/G_ZN%lamg              ! (K*cm*s)/g
        J(2,2)=-1.0/xg-1.0/xcd                                         ! 1/cm

!
!       update guesses for m and Ts using Newton-Rhapson
!       & guard against a negative (erroneous) rb solution
!
        CALL Newton(xx,J,e)
        IF(xx(1).LT.0.0) THEN

            PRINT *,' ROCBURN_ZN: rank=',G_ZN%rank
            PRINT *,' rb negative xx(1)=',xx(1)
            xx(1)=0.0d0

        END IF
        

        IF(ABS(e(1)).LT.tol .AND. ABS(e(2)).LT.tol) THEN

!
!           reset final converged values
!
            rb=xx(1)/rhoc                                              ! cm/s
            Ts=xx(2)                                                   ! K
            xreac=2.0*G_ZN%alfac*G_ZN%R*Ts/rb/G_ZN%Ec                  ! 1/cm
            fr=dexp(-xreac*G_ZN%Ka)                                    ! -
            xcd=G_ZN%lamg/G_ZN%C/xx(1)                                 ! cm
            Da=4.0*G_ZN%Bg*G_ZN%C/G_ZN%lamg*(P*G_ZN%MW*xcd/G_ZN%R)**2  ! -
            term4=(1.0+Da)**0.5                                        ! -
            xg=2.0*xcd/(term4-1.0)                                     ! cm
            Tf=(qr/xx(1)+G_ZN%Qc+G_ZN%Qg)/G_ZN%C+To                    ! K

!
!           Sets ICs using I&W S.S. T profile
!

            xcond=G_ZN%alfac/rb
            xreac=xcond*G_ZN%R*Ts*2.0/G_ZN%Ec
            fr=dexp(-xreac*G_ZN%Ka)
            frJ=fr*qr/rhoc/rb/G_ZN%C/(Ts-To)/(1.-G_ZN%Ka*xcond)
            DO i=1,G_ZN%nx
               Tn(i)=To+(Ts-To)*((1.0-frJ)*    &
                  dexp(G_ZN%x(i)/xcond)+ frJ*exp(G_ZN%x(i)*G_ZN%Ka))
            END DO ! i

            IF(ABS(Tn(G_ZN%nx)-To)/To .GT. 0.001) THEN
   
               WRITE(*,*) 'ROCBURN_ZN: rank=',G_ZN%rank
               WRITE(*,*) '  Error: xmax not large enough'
               WRITE(*,*) '  To boundary condition not satisfied'
               CALL MPI_ABORT( MPI_COMM_WORLD, -1)
               STOP

            ELSE
   
               Tn(G_ZN%nx)  =To
   
            END IF

            RETURN 
         END IF
      END DO ! iter
         
      PRINT *,'  ROCBURN_ZN: rank=',G_ZN%rank
      PRINT *,'  ROCBURN_ZN: Error-steady state WSB'    &
             ,'  in ignition_combustion'
      PRINT *,'  ROCBURN_ZN: Maximum # of ',G_ZN%itermax     &
             , ' iterations exceeded'
      CALL MPI_ABORT( MPI_COMM_WORLD, -1)
      STOP
!
!     end of steady state WSB homogeneouse propellant combustion model
!
    END IF


    IF(G_ZN%Model_combustion == 2) THEN

!     Model_combustion = 2
!
!     ZN approach for composite propellant with empirical formulatin
!     for the gas phase reaction law
!
!     initial guess of rb and Ts
!

      rb=G_ZN%a_p*P**G_ZN%n_p                                           ! cm/s
      Ts=727.002                                                        ! K
      xx(1)=rb*rhoc                                                     ! g/(cm^2*s)
      xx(2)=Ts                                                          ! K

      DO iter=0,G_ZN%itermax

!
!        set up residual matrix and Jacobian derivative matrix for
!        Ibiricu and Williams condensed phase expression
!
         xreac=2.0*rhoc*G_ZN%alfac*G_ZN%R*xx(2)/xx(1)/G_ZN%Ec           ! cm
         term0=1.0-G_ZN%Ka*xreac                                        ! -
         term1=G_ZN%Ac*G_ZN%R*(xx(2)*rhoc)**2*G_ZN%alfac*G_ZN%C*                   &
                      dexp(-G_ZN%Ec/G_ZN%R/xx(2))/G_ZN%Ec               ! cal*g/(cm^4*s^2)

         fr=dexp(-xreac*G_ZN%Ka)
         term2=G_ZN%C*(xx(2)-To)-G_ZN%Qc/2.0-fr*qr/xx(1)                ! cal/g
         term3=term1/term2                                              ! g^2/(cm^4*s^2)

         e(1)=xx(1)**2-term3                                            ! g^2/(cm^4*s^2)
         J(1,1)=2.0*xx(1)+term3/term2*fr*qr/xx(1)**2           &
                          *term0                                        ! g/(cm^2*s)
         J(1,2)=-xx(1)**2*((2.0d0+G_ZN%Ec/G_ZN%R/xx(2))/xx(2) -          &
                 (G_ZN%C+G_ZN%Ka*xreac*qr*fr/xx(1)/xx(2))/term2)        ! g^2/(cm^4*s^2*K)
         e(2)=G_ZN%a_T*((xx(2) - To - qr/(xx(1)*G_ZN%C) )/               &
              (xx(2)-300.0))**G_ZN%n_T - rhoc*G_ZN%a_p*p**G_ZN%n_p/xx(1)
         J(2,1)=rhoc*G_ZN%a_p*p**G_ZN%n_p/(xx(1)*xx(1))
         J(2,2)=G_ZN%n_T*G_ZN%a_T*((xx(2)-To)/(xx(2)-300.0))**(G_ZN%n_T-1.0)  &
                *(To-300.0)/(xx(2)-300.0)**2
!
!        update guesses for m and Ts using Newton-Rhapson
!        & guard against a negative (erroneous) rb solution
!
         CALL Newton(xx,J,e)
         IF(xx(1).lt.0.0) THEN

             PRINT *,' ROCBURN_ZN: rank=',G_ZN%rank
             PRINT *,' rb negative xx(1)=',xx(1)
             xx(1)=0.0d0

         END IF

         IF(ABS(e(1)).LT.tol .AND. ABS(e(2)).LT.tol) THEN

!
!           reset final converged values
!
            rb=xx(1)/rhoc                                              ! cm/sec
            Ts=xx(2)                                                   ! K
            xreac=2.0*G_ZN%alfac*G_ZN%R*Ts/rb/G_ZN%Ec                  ! 1/cm
            fr=dexp(-xreac*G_ZN%Ka)                                    ! -
            Tf=(qr/xx(1)+G_ZN%Qc+G_ZN%Qg)/G_ZN%C+To                    ! K

!
!           Sets ICs using I&W S.S. T profile
!

            xcond=G_ZN%alfac/rb
            xreac=xcond*G_ZN%R*Ts*2.0/G_ZN%Ec
            fr=dexp(-xreac*G_ZN%Ka)
            frJ=fr*qr/rhoc/rb/G_ZN%C/(Ts-To)/(1.-G_ZN%Ka*xcond)
            DO i=1,G_ZN%nx
               Tn(i)=To+(Ts-To)*((1.0-frJ)*         &
               dexp(G_ZN%x(i)/xcond)+frJ*exp(G_ZN%x(i)*G_ZN%Ka))
            END DO ! i

            IF(abs(Tn(G_ZN%nx)-To)/To .gt. 0.001) THEN
  
               WRITE(*,*) ' ROCBURN: rank=',G_ZN%rank
               WRITE(*,*) ' ROCBURN: Error: xmax'         &
                         ,' not large enough'
               WRITE(*,*) ' ROCBURN: To boundary'         &
                         ,' condition not satisfied'
               CALL MPI_ABORT( MPI_COMM_WORLD, -1)
               STOP

            ELSE
  
               Tn(G_ZN%nx)  =To
  
            END IF

            RETURN

         END IF

      END DO ! iter

      PRINT *,'  ROCBURN_ZN: rank=',G_ZN%rank
      PRINT *,'  ROCBURN_ZN: Error-steady state ZN'      &
             ,'  in ignition_combustion'
      PRINT *,'  ROCBURN_ZN: Maximum # of ',G_ZN%itermax &
             , ' iterations exceeded'
      CALL MPI_ABORT( MPI_COMM_WORLD, -1)
      STOP

!
!        end of steady state ZN phenomenological combustion model
!
    END IF

    IF(G_ZN%Model_combustion == 3) THEN

!
!     Model_combustion = 3
!
!     rb=a*P**n
!
!        initial guess of rb and Ts
!

      rb=G_ZN%a_p*P**G_ZN%n_p                                      ! cm/s
      Tn   =1000.0d0
      Ts   =600.00d0
         
      RETURN

    END IF

    WRITE(*,*) ' ROCBURN_ZN: rank=',G_ZN%rank
    WRITE(*,*) ' ROCBURN_ZN: Error-ssWSB:'                   
    WRITE(*,*) ' No appropriate combustion model'
    WRITE(*,*) ' ROCBURN_ZN: Model_combustion=', G_ZN%Model_combustion
    CALL MPI_ABORT( MPI_COMM_WORLD, -1)
    STOP

    RETURN

    CONTAINS

!     ------------------------------------------------------------------------
!                             INTERNAL PROCEDURES
!     ------------------------------------------------------------------------


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

      IMPLICIT NONE

      REAL(DBL) :: x(2),J(2,2),e(2),detJ,delx(2),dum
      integer i

      detJ=J(1,1)*J(2,2)-J(2,1)*J(1,2)

      if(detJ.eq.0) then
        write(*,*) 'ROCBUNR_ZN Error: singular matrix encountered-could not invert'
        stop
      endif
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
      do 10 i=1,2
        delx(i)=J(i,1)*e(1)+J(i,2)*e(2)
        x(i)=x(i)-delx(i)
   10 continue

      RETURN
      END SUBROUTINE Newton

!     -------------------------------------------------------------------
!                       END OF INTERNAL PROCEDURES
!     -------------------------------------------------------------------

  END SUBROUTINE ZN_ssWSB






