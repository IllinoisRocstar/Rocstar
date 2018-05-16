c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c   subroutine  : burning_rate                                                c
c                                                                             c
c   purpose     : calculate the burning rate                                  c
c                                                                             c
c --------------------------------------------------------------------------- c
c                                                                             c
c   arguments   :                                                             c
c
c --------------------------------------------------------------------------- c
c                                                                             c
c                            written by                                       c
c                      K.C. Tang and M.Q. Brewster                            c
c                                                                             c
c                       Last modified : 05/03/2001                            c
c                                                                             c
!  Modifications:
!
!    No.     Date         Programmer    Description
!    001  May  03, 2001   K. C. Tang    change memory allocation method for
!                                       temporary arrarys inside a 
!                                       subroutine
c                                                                             c
c --------------------------------------------------------------------------  c
c                                                                             c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE ZN_burning_rate (P, qr, To, qr_old, fr_old, 
     &                         Toa, rb, Ts, fr)
  
!     USE burn_global_data
!     USE burn_rate_table_data
      IMPLICIT NONE

      INCLUDE 'mpif.h'

      REAL(DBL) :: rb, Ts, fr, Toa
      REAL(DBL) :: P, qr, To
      REAL(DBL) :: qr_old, rb_old, fr_old
 

c
c ============================================================================
c
c     delcare local variables
c
      INTEGER   :: iter, itrmx, i
      REAL(DBL) :: rblast, Tslast
      REAL(DBL) :: fs
      REAL(DBL), DIMENSION (nxmax) ::  first, second
      INTEGER   :: ierror

c
c ============================================================================
c



      CALL MPI_COMM_RANK(MPI_COMM_ROCBURN, rank, ierror)

! 04.18.01 KCT
! Move Combustion model 3 (rb=a*P**n) from subroutine ZN
! to burn_rate 
!  

      IF(Model_combustion.eq.3) then

c
c
c        Model_combustion = 3
c
c        rb=a*P**n
c

         rb=a_p*P**n_p
         Tnp1(1)=1000.0
   
         RETURN

      END IF

! 04.18.01 KCT
! Move Combustion model 3 (rb=a*P**n) from subroutine ZN
! to subroutine burn_rate 
!  

      itrmx = 100  

      rb_old = rb

      fs=rb*(Ts-Toa)/alfac-fr_old*qr_old/lamc

      DO i=1,nxmax
        IF( (Tn(i) > 10000.0).OR.(Tn(i) < 0.001)) THEN
           PRINT *,' ROCBURN: rank=',rank,
     +             ' inside ZN  Tn(',i,')=',Tn(i)
        END IF
      END DO
!     PRINT *,' ROCBURN: inside burn_rate'
!     PRINT *,' ROCBURN: delt=',delt,' alfac=',alfac,' rb_old=',rb_old
!     PRINT *,' ROCBURN: Ts=',Ts,'delz=',delz, 'P=',P,' qr=',qr,' To=',To
!     PRINT *,' ROCBURN: rb=',rb,' lamc=',lamc
!     PRINT *,' ROCBURN: rhoc=',rhoc,' Ka    ',Ka   ,' fr_old=',fr_old, 'Toa=',Toa
!     PRINT *,' ROCBURN: C   =',C   ,' qr_old=',qr_old


      DO 1000 iter=0,itrmx

c
c        get 1st and 2nd derivatives at last time step (n)
c        with this iterations surface temperature B.C.

         CALL cmpct1(nx,delz,Tn,first)
         DO i=1,nxmax
           IF( (Tn(i) > 10000.0).OR.(Tn(i) < 0.001)) THEN
              PRINT *,' ROCBURN: rank=',rank
              PRINT *,' ROCBURN: before ZN  after ct1 Tn(',i,')=',Tn(i)
           END IF
         END DO
         CALL cmpct2(nx,delz,Tn,second)
         DO i=1,nxmax
           IF( (Tn(i) > 10000.0).OR.(Tn(i) < 0.001)) THEN
              PRINT *,' ROCBURN: rank=',rank
              PRINT *,' ROCBURN: before ZN  after ct2 Tn(',i,')=',Tn(i)
           END IF
         END DO

c
c        determine T profile at this time step (n+1) based on
c        properties a last time step (n) (explicit differencing)
c
         DO 2000 i=2,nx-1
            Tnp1(i)=Tn(i)+delt*(
     +              (alfac*zxx(i)-rb_old*zx(i))*first(i) +
     +              alfac*zx(i)**2*second(i) +
     +              fr_old*qr_old/rhoc/C*Ka*exp(Ka*x(i)))
 2000    CONTINUE


         DO i=1,nxmax
           IF( (Tn(i) > 10000.0).OR.(Tn(i) < 0.001)) THEN
              PRINT *,' ROCBURN: rank=',rank,
     +                ' before ZN  after Tnp1 Tn(',i,')=',Tn(i)
           END IF
         END DO

         Tnp1(nx)=To
  

c        update suface B.C. for this iteration
c        Tnp1(1)=(48.0*Tnp1(2)-36.0*Tnp1(3)+16.0*Tnp1(4)       ! 4th
c    +            -3.0*Tnp1(5)-12.0*delz*fs/zx(1))/25.0        ! 4th
         Tnp1(1)=(18.0*Tnp1(2)-9.0*Tnp1(3)+2.0*Tnp1(4)         ! 3rd
     +           -6.0*delz*fs/zx(1))/11.0                      ! 3rd
c        Tnp1(1)=(Tnp1(2)-delz*fs/zx(1))                       ! 1st

   

         Ts=Tnp1(1)
         Tn(1)=Tnp1(1)

!        DO i=1,nxmax
!          IF( (Tn(i) > 10000.0).OR.(Tn(i) < 0.001)) THEN
!             PRINT *,' ROCBURN: rank=',rank,
!    +                ' before ZN  after Tnp1 Tn(',i,')=',Tn(i)
!          END IF
!        END DO

!        do i=1,nxmax
!          IF( (Tnp1(i) > 10000.0).OR.(Tnp1(i) < 0.001)) THEN
!             PRINT *,' ROCBURN: rank=',rank,
!    +                ' before ZN  after Tnp1 Tnp1(',i,')=',Tnp1(i)
!          END IF
!        end do

!        IF( (Tnp1(1) > 10000.0).OR.(Tnp1(1) < 0.001)) THEN
!           PRINT *,' ROCBURN: rank=',rank
!           PRINT *,' ROCBURN: before ZN  ',
!    +              ' after Tnp1 Tnp1(',i,')=',Tnp1(i)
!           PRINT *,' ROCBURN: delt=',delt,
!    +              ' alfac=',alfac,' rb_old=',rb_old
!           PRINT *,' ROCBURN: rhoc=',rhoc,' Ka    ',Ka,
!    +              ' fr_old=',fr_old
!           PRINT *,' ROCBURN: C   =',C   ,' qr_old=',qr_old
!           PRINT *,'  '
!           PRINT *,' ROCBURN: x=',x
!           PRINT *,'  '
!           PRINT *,' ROCBURN: zx=',zx
!           PRINT *,'  '
!           PRINT *,' ROCBURN: zxx=',zxx
!        END IF
c
c        calculate rb, Toa, and fr for given  P, qr, and Ts
c
c        print *,' before ZN Ts=',Ts, 'zx(1)=',zx(1),' fs=',fs,
c    +           ' delz=',delz
c        print *,'Tn=',Tn
c        print *,'Tnp1=',Tnp1
 
         CALL ZN(P, qr, To, Ts, Toa, rb, fr)

c
c        calculate fs for this iteration of (n+1)
c
         fs=rb*(Ts-Toa)/alfac-fr*qr/lamc

!        PRINT *,'iter=',iter,' fs=',fs,' rb=',rb,' Ts=',Ts
         IF(iter.ge.1) THEN
            IF(abs(Ts-Tslast).lt.tol_Ts.and.abs(rb-rblast).lt.tol_Ts) 
     &                                                          THEN
!        PRINT *,'inside burn_rate rb=',rb,' iter=',iter,'itrmx=',itrmx
!        PRINT *,'tol_Ts=',tol_Ts
               RETURN
            END IF
         END IF

         Tslast=Ts
         rblast=rb

 1000 CONTINUE
    
      PRINT *,' ROCBURN: rank=',rank
      PRINT *,' ROCBURN: Error - burning rate: Maximum # of iterations '
      PRINT *,' ROCBURN: itrmx=', itrmx, 'reached'
      PRINT *,' ROCBURN: tol_Ts=',tol_Ts
      PRINT *,' ROCBURN: Tn=',Tn
      PRINT *,' ROCBURN: Tnp1=',Tnp1
      STOP


      CONTAINS

!     ------------------------------------------------------------------------
!                             INTERNAL PROCEDURES
!     ------------------------------------------------------------------------
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c   subroutine  : ZN                                                          c
c                                                                             c
c   purpose     :  alculate rb,Toa given Ts,P,qr                              c
c                                                                             c
c --------------------------------------------------------------------------- c
c                                                                             c
c                            written by                                       c
c                      K.C. Tang and M.Q. Brewster                            c
c                                                                             c
c                       Last modified : 09/14/2000                            c
c --------------------------------------------------------------------------- c
c                                                                             c
c   input       :                                                             c
c                                                                             c
c      P        : pressure (atm)                                              c
c      qr       : radiative heat flux (al/cm^2-s)                             c
c      To       : propellant temperatuer at deep into the propellant (K)      c
c      Ts       : surface temperature (K)                                     c
c      Toa      : ZN temperature (K) for initial guess                        c
c      rb       : burning rate (cm/s) for initial guess                       c
c                                                                             c
c                                                                             c
c   output      :                                                             c
c                                                                             c
c      Toa      : ZN temperature (K)                                          c
c      rb       : burning rate (cm/s)                                         c
c      fr       : transmissivity of surface reaction zone                     c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c
      SUBROUTINE ZN(P, qr, To, Ts, Toa, rb, fr)

      IMPLICIT NONE

      REAL(DBL) :: P, qr, To
      REAL(DBL) :: Ts, rb, Toa, fr
      REAL(DBL) :: Tf
      REAL(DBL) :: delta_T2, Ag


      REAL(DBL), parameter :: tol=1.d-5

c ============================================================================
c
c     delcare local variables
c

      REAL(DBL) :: J(2,2),xx(2),e(2)
      REAL(DBL) :: xreac
      REAL(DBL) :: term0, term1, term2, term3, term4, term5
      REAL(DBL) :: xcd, Da, xg

      integer iter


      term1=Ac*R*(Ts*rhoc)**2*alfac*C*exp(-Ec/R/Ts)/Ec


c
c     initial guess
c

      xx(1)=rb*rhoc
      xx(2)=Toa

      IF(Model_combustion == 1) THEN

c
c
c        Model_combustion = 1
c
c
c
c        WSB homogeneouse propellant combustion model
c
c

         Ag=Bg/(R*R)*MW*MW


         do 1000 iter=0,itermax

c
c           set up residual matrix and Jacobian derivative matrix for
c           Ibiricu and Williams condensed phase expression

            xreac=2.0*rhoc*alfac*R*Ts/x(1)/Ec
            term0=1.0-Ka*xreac

            fr=dexp(-xreac*Ka)
            term2=C*(Ts-xx(2))-Qc/2.0-fr*qr/xx(1)
            term3=term1/term2

            e(1)=xx(1)**2-term3
            J(1,1)=2.0*xx(1)+term3/term2*fr*qr/xx(1)**2*term0
            J(1,2)=-term3*C/term2

c
c           set up residual matrix and Jacobian derivative matrix for
c           Ward-Son-Brewster gas phase expression


            xcd=lamg/C/xx(1)
            Da=4.0*Bg*C/lamg*(P*MW*xcd/R)**2
            term4=(1.0+Da)**0.5
            xg=2.0*xcd/(term4-1.0)
            term5=Qc-C*(Ts-xx(2))
            Tf=(qr/xx(1)+Qc+Qg)/C+xx(2)

            e(2)=(Tf-Ts)/xg+(qr+xx(1)*term5)/lamg
            J(2,1)=-qr/xg/xx(1)**2/C-(Tf-Ts)/xg/xx(1)/term4+term5/lamg
            J(2,2)=1/xg+1/xcd

c
c           update guesses for m and Toa using Newton-Rhapson
c           & guard against a negative (erroneus) rb solution
c     
            call Newton(xx,J,e)
            if(xx(1).lt.0.0) then
               PRINT *,' ROCBURN: rank=',rank,
     +                 ' rb negative xx(1)=',xx(1),' Ts=',Ts
               xx(1)=abs(xx(1))
            end if

            if(abs(e(1)).lt.tol .and. abs(e(2)).lt.tol) then

c
c              reset final converged values
c
               rb=xx(1)/rhoc
               Toa=xx(2)
               fr=dexp(-xreac*Ka)
               return

            end if
 1000    continue

         WRITE(*,*) ' ROCBURN: in ZN rank=',rank,
     +              ' Error-WSB: Maximum # of iterations reached'
         STOP

      end if

      if(Model_combustion.eq.2) then

c
c
c        Model_combustion = 2
c
c        Zeldovich-Novozhilov phenomenological model for composite
c        propellant
c

         do 2000 iter=0,itermax

c
c           set up residual matrix and Jacobian derivative matrix for
c           Ibiricu and Williams condensed phase expression

            xreac=2.0*rhoc*alfac*R*Ts/xx(1)/Ec
            term0=1.0-Ka*xreac

            fr=dexp(-xreac*Ka)
            term2=C*(Ts-xx(2))-Qc/2.0-fr*qr/xx(1)
            term3=term1/term2

            e(1)=xx(1)**2-term3
            J(1,1)=2.0*xx(1)+term3/term2*fr*qr/xx(1)**2*term0
            J(1,2)=-term3*C/term2

c
c           set up residual matrix and Jacobian derivative matrix for
c           Ward-Son-Brewster gas phase expression
c


            delta_T2 = Ts - 300.0
            e(2)=a_T*((Ts-xx(2)-qr/(xx(1)*C))/delta_T2)**n_T -
     +           rhoc*a_p*p**n_p/xx(1)
            J(2,1)=rhoc*a_p*p**n_p/(xx(1)*xx(1))
            J(2,2)=-n_T*a_T*(Ts-xx(2)-qr/(xx(1)*C))**(n_T-1.0)/
     +              (delta_T2**n_T)
         

c
c           update guesses for m and Toa using Newton-Rhapson
c           & guard against a negative (erroneus) rb solution
c    

            call Newton(xx,J,e)

            if(xx(1).lt.0.0) then
               PRINT *,' ROCBURN: inside ZN rank=',rank,
     +                 ' rb negative xx(1)=',xx(1),' Ts=',Ts
               xx(1)=abs(xx(1))
            end if

            if(abs(e(1)).lt.tol .and. abs(e(2)).lt.tol) then

c
c              reset final converged values
c
               rb=xx(1)/rhoc
               Toa=xx(2)
               fr=dexp(-xreac*Ka)
               RETURN

            END IF
 2000    CONTINUE

         WRITE(*,*) ' ROCBURN: in ZN rank=',rank,
     +              ' Error-Composite: Maximum # of iterations reached'
         STOP

      end if

      IF(Model_combustion.eq.3) THEN

c
c
c        Model_combustion = 3
c
c        rb=a*P**n
c

         rb=a_p*P**n_p
         Tnp1(1)=1000.0
   
         RETURN

      END IF

      WRITE(*,*) ' ROCBURN: in ZN rank=',rank,
     +           ' Error-ZN: No appropriate combustion model'
      WRITE(*,*) ' ROCBURN: Combustion Model= ', Model_combustion
      STOP

      END SUBROUTINE ZN



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c   subroutine  : Newton                                                      c
c                                                                             c
c                                                                             c
c   Author:                                                                   c
c                                                                             c
c      Paul Loner                                                             c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
      SUBROUTINE Newton(x,J,e)

      implicit none

      REAL(DBL) :: x(2),J(2,2),e(2),detJ,delx(2),dum
      integer i

      detJ=J(1,1)*J(2,2)-J(2,1)*J(1,2)

      IF(detJ.eq.0) then
        PRINT *,' ROCBURN: in Newton rank=',rank,
     +          ' Error: singular matrix encountered'
        PRINT *,' ROCBURN: P=', P
        PRINT *,' ROCBURN: J=', J
        STOP
      ENDIF
c
c     invert the 2x2 Jacobian matrix
c
      dum=J(1,1)
      J(1,1)=J(2,2)
      J(2,2)=dum
      J(1,1)=J(1,1)/detJ
      J(1,2)=-J(1,2)/detJ
      J(2,1)=-J(2,1)/detJ
      J(2,2)=J(2,2)/detJ
c
c     multiply J^-1()*e() and update guesses
c
      do 10 i=1,2
        delx(i)=J(i,1)*e(1)+J(i,2)*e(2)
        x(i)=x(i)-delx(i)
   10 continue

      RETURN
      END SUBROUTINE Newton

c***********************************************************************
c
      SUBROUTINE cmpct1(imax,h,u,f)
c
c***********************************************************************
!  Modifications:
!
!    No.     Date         Programmer    Description
!    001  May  03, 2001   K. C. Tang    change memory allocation method for
!                                       temporary arrarys inside a 
!                                       subroutine

      IMPLICIT NONE

      integer i,imax
      REAL(DBL) :: h,u(imax),f(imax)

      REAL(DBL), ALLOCATABLE :: a(:),b(:),cc(:)


c     set up equations for first derivatives
c      print *,' in cmpct1 imax=',imax,' h=',h,' u=',u,' f=',f

      ALLOCATE(a(imax))
      ALLOCATE(b(imax))
      ALLOCATE(cc(imax))

      a(1)=0.0
      b(1)=1.0
      cc(1)=0.0
      f(1)=(-25.0*u(1)+48.0*u(2)-36.0*u(3)+16.0*u(4)-3.0*u(5))/(12.0*h) !4th

      do 300 i=2,imax-1
        a(i)=1.0
        b(i)=4.0
        cc(i)=1.0
        f(i)=3.0*(u(i+1)-u(i-1))/h
  300 continue

      a(imax)=0.0
      b(imax)=1.0
      cc(imax)=0.0
      f(imax)=0.0

c     solve for first derivatives

      call tridg (imax,a,b,cc,f)

      DEALLOCATE(a)
      DEALLOCATE(b)
      DEALLOCATE(cc)
      RETURN
      END SUBROUTINE cmpct1


c***********************************************************************
c
      SUBROUTINE cmpct2(imax,h,u,s)
c
c***********************************************************************
!  Modifications:
!
!    No.     Date         Programmer    Description
!    001  May  03, 2001   K. C. Tang    change memory allocation method for
!                                       temporary arrarys inside a 
!                                       subroutine

      implicit none

      integer i,imax
      REAL(DBL) :: h,h2,u(imax),s(imax)
      REAL(DBL), ALLOCATABLE :: a(:),b(:),cc(:)


      h2=h*h

      ALLOCATE(a(imax))
      ALLOCATE(b(imax))
      ALLOCATE(cc(imax))

c     set up equations for second derivatives

      a(1)=0.0
      b(1)=1.0
      cc(1)=0.0
      s(1)=(45.0*u(1)-154.0*u(2)+214.0*u(3)-156.0*u(4)+61.0*u(5)-    ! 4th
     +     10.0*u(6))/12.0/h2                                        ! 4th

      do 400 i=2,imax-1
        a(i)=1.0
        b(i)=10.0
        cc(i)=1.0
        s(i)=12.0*(u(i+1)-2.0*u(i)+u(i-1))/h2
  400 continue

      a(imax)=0.0
      b(imax)=1.0
      cc(imax)=0.0
      s(imax)=0.0

c     solve for second derivatives

      call tridg (imax,a,b,cc,s)

      DEALLOCATE(a)
      DEALLOCATE(b)
      DEALLOCATE(cc)

      RETURN
      END SUBROUTINE cmpct2


c***********************************************************************
c
      SUBROUTINE tridg (imax,a,b,c,d)
c
c     solves a tridiagonal matrix equation of the form:
c
c       a(i)*u(i-1) + b(i)*u(i) + c(i)*u(i+1) = d(i)
c
c     using the Thomas algorithm.
c
c***********************************************************************

      IMPLICIT NONE

      INTEGER :: imax
      REAL(DBL), DIMENSION(imax) :: a, b, c, d
 
      INTEGER ::  i
c     integer i,imax
c     real a(1),b(1),c(1),d(1)

      do 100 i=2,imax
        a(i)=a(i)/b(i-1)
        b(i)=b(i)-a(i)*c(i-1)
        d(i)=d(i)-a(i)*d(i-1)
  100 continue

      d(imax)=d(imax)/b(imax)

      do 200 i=imax-1,1,-1
        d(i)=(d(i)-c(i)*d(i+1))/b(i)
  200 continue
      RETURN
      END SUBROUTINE tridg


!     -------------------------------------------------------------------
!                       END OF INTERNAL PROCEDURES
!     -------------------------------------------------------------------


      END SUBROUTINE burning_rate
