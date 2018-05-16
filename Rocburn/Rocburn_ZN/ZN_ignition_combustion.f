c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c   subroutine  : ignition_comb                                               c
c                                                                             c
c   purpose     : provide initial condition for combustion model after the    c
c                 burning surface element being designated as ignited         c
c                                                                             c
c --------------------------------------------------------------------------- c
c                                                                             c
c   arguments   :                                                             c
c                                                                             c
c                                                                             c
c --------------------------------------------------------------------------- c
c                                                                             c
c                            written  by                                      c
c                     K.C Tang and M.Q Brewster                               c 
c                                                                             c
c                     Last modified :  09/15/2000                             c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE ignition_combustion(P, qr, To, rb, Ts, Tf, fr)

      USE M_Rocburn_ZN_Global_Data
      IMPLICIT NONE

      REAL(DBL), INTENT(IN) :: P, qr, To
      REAL(DBL), INTENT(OUT) :: rb, Ts, Tf, fr

  
      REAL(DBL), parameter :: tol=1.0d-5


c
c ============================================================================
c
c     delcare local variables
c
      REAL(DBL) :: J(2,2),xx(2),e(2)
      REAL(DBL) :: xreac
      REAL(DBL) :: term0, term1, term2, term3, term4, term5
      REAL(DBL) :: xcd, Da, xg, xcond, frJ

      INTEGER   :: iter, i

c
c ============================================================================
c


c
c     initial condition
c
c     calculate Ts, rb for given To, P, qr
c
c


c
c     Generate initial conditions
c
c
c

      if(Model_combustion.eq.1) then

c
c
c        Model_combustion = 1
c
c
c
c        steady state WSB homogeneouse propellant combustion model
c
c        guess initial values for doublebase propellant 
c        rb and Ts 
c


       
         rb=a_p*(P**n_p)                                    ! cm/s
         Ts=550.0                                           ! K
         xx(1)=rb*rhoc                                      ! g/(cm^2*s)
         xx(2)=Ts                                           ! K
         do 1000 iter=0,itermax

c
c           set up residual matrix and Jacobian derivative matrix for
c           Ibiricu and Williams condensed phase expression
c   

            xreac=2.0*rhoc*alfac*R*xx(2)/xx(1)/Ec           ! cm
            term0=1.0-Ka*xreac                              ! -

            term1=Ac*R*(xx(2)*rhoc)**2*alfac*C*
     +                   dexp(-Ec/R/xx(2))/Ec               ! cal*g/(cm^4*s^2)
 
            fr=dexp(-xreac*Ka)
            term2=C*(xx(2)-To)-Qc/2.0-fr*qr/xx(1)           ! cal/g
            term3=term1/term2                               ! g^2/(cm^4*s^2)

            e(1)=xx(1)**2-term3                             ! g^2/(cm^4*s^2)
            J(1,1)=2.0*xx(1)+term3/term2*fr*qr/xx(1)**2
     +                       *term0                         ! g/(cm^2*s)
            J(1,2)=-xx(1)**2*((2.0d0+Ec/R/xx(2))/xx(2) - 
     +              (C+Ka*xreac*qr*fr/xx(1)/xx(2))/term2)   ! g^2/(cm^4*s^2*K)

c
c           set up residual matrix and Jacobian derivative matrix for
c           Ward-Son-Brewster gas phase expression

            xcd=lamg/C/xx(1)                                ! cm
            Da=4.0*Bg*C/lamg*(P*MW*xcd/R)**2                ! -
            term4=(1.0+Da)**0.5                             ! -
            xg=2.0*xcd/(term4-1.0)                          ! cm
            term5=Qc-C*(xx(2)-To)                           ! cal/g
            Tf=(qr/xx(1)+Qc+Qg)/C+To                        ! K

            e(2)=(Tf-xx(2))/xg+(qr+xx(1)*term5)/lamg        ! K/cm
            J(2,1)=-qr/xg/xx(1)**2/C-
     +              (Tf-xx(2))/xg/xx(1)/term4+term5/lamg    ! (K*cm*s)/g
            J(2,2)=-1.0/xg-1.0/xcd                          ! 1/cm

c
c           update guesses for m and Ts using Newton-Rhapson
c           & guard against a negative (erroneous) rb solution
c
            call Newton(xx,J,e)
            if(xx(1).lt.0.0) then

                print *,' rank=',rank
                print *,' rb negative xx(1)=',xx(1)
                xx(1)=0.0d0

            endif
        

            if(abs(e(1)).lt.tol .and. abs(e(2)).lt.tol) then

c
c               reset final converged values
c
                rb=xx(1)/rhoc                                 ! cm/s
                Ts=xx(2)                                      ! K
                xreac=2.0*alfac*R*Ts/rb/Ec                    ! 1/cm
                fr=dexp(-xreac*Ka)                            ! -
                xcd=lamg/C/xx(1)                              ! cm
                Da=4.0*Bg*C/lamg*(P*MW*xcd/R)**2              ! -
                term4=(1.0+Da)**0.5                           ! -
                xg=2.0*xcd/(term4-1.0)                        ! cm
                Tf=(qr/xx(1)+Qc+Qg)/C+To                      ! K

c
c               set up initial temperature profile
c

                if (Tintype .eq. 1) then
c
c
c                  Sets ICs using I&W S.S. T profile
c

                   xcond=alfac/rb
                   xreac=xcond*R*Ts*2.0/Ec
                   fr=dexp(-xreac*Ka)
                   frJ=fr*qr/rhoc/rb/C/(Ts-To)/(1.-Ka*xcond)
                   do 3000 i=1,nx
                      Tn(i)=To+(Ts-To)*((1.0-frJ)*
     +                   dexp(x(i)/xcond)+ frJ*exp(x(i)*Ka))
                      Tnp1(i)=Tn(i)
 3000              continue

                   if(abs(Tn(nx)-To)/To .gt. 0.001) then
   
                      write(*,*) 'rank=',rank
                      write(*,*) 'Error: xmax not large enough'
                      write(*,*) 'To boundary condition not satisfied'
                      stop

                   else
   
                      Tn(nx)  =To
                      Tnp1(nx)=To
   
                   end if

                else

c
c                  Set ICs using data in Tinit.dat
c

                   open(unit=16,file='Tinit.dat',status='unknown')
                   read(16,*)
                   read(16,*)
                   read(16,*)
                   read(16,*)
                   read(16,*)
                   read(16,*)
                   do 3010 i = 1,nx
                      read(16,*) x(i), Tn(i)
                      Tnp1(i)=Tn(i)
 3010              continue

                endif

                return 

             endif

 1000    continue
         
         PRINT *,'  ROCBURN: rank=',rank
         PRINT *,'  ROCBURN: Error-steady state WSB'
     +          ,'  in ignition_combustion'
         PRINT *,'  ROCBURN: Maximum # of ',itermax
     +          , ' iterations exceeded'
         STOP
c
c        end of steady state WSB homogeneouse propellant combustion model
c
      end if
c
c
c
c
c
c


      if(Model_combustion.eq.2) then

c
c
c
c        Model_combustion = 2
c
c
c
c        ZN approach for composite propellant
c
c
c
c           initial guess of rb and Ts
c

         rb=a_p*P**n_p                                      ! cm/s
         Ts=727.002                                         ! K
         xx(1)=rb*rhoc                                      ! g/(cm^2*s)
         xx(2)=Ts                                           ! K

         do 2000 iter=0,itermax

c
c           set up residual matrix and Jacobian derivative matrix for
c           Ibiricu and Williams condensed phase expression
c
            xreac=2.0*rhoc*alfac*R*xx(2)/xx(1)/Ec           ! cm
            term0=1.0-Ka*xreac                              ! -
            term1=Ac*R*(xx(2)*rhoc)**2*alfac*C*
     +                   dexp(-Ec/R/xx(2))/Ec               ! cal*g/(cm^4*s^2)

            fr=dexp(-xreac*Ka)
            term2=C*(xx(2)-To)-Qc/2.0-fr*qr/xx(1)           ! cal/g
            term3=term1/term2                               ! g^2/(cm^4*s^2)

            e(1)=xx(1)**2-term3                             ! g^2/(cm^4*s^2)
            J(1,1)=2.0*xx(1)+term3/term2*fr*qr/xx(1)**2
     +                       *term0                         ! g/(cm^2*s)
            J(1,2)=-xx(1)**2*((2.0d0+Ec/R/xx(2))/xx(2) -
     +              (C+Ka*xreac*qr*fr/xx(1)/xx(2))/term2)   ! g^2/(cm^4*s^2*K)
            e(2)=a_T*((xx(2) - To - qr/(xx(1)*C) )/
     +           (xx(2)-300.0))**n_T - rhoc*a_p*p**n_p/xx(1)
            J(2,1)=rhoc*a_p*p**n_p/(xx(1)*xx(1))
            J(2,2)=n_T*a_T*((xx(2)-To)/(xx(2)-300.0))**(n_T-1.0)
     +             *(To-300.0)/(xx(2)-300.0)**2
c
c           update guesses for m and Ts using Newton-Rhapson
c           & guard against a negative (erroneous) rb solution
c
            call Newton(xx,J,e)
            if(xx(1).lt.0.0) then

                print *,' rank=',rank
                print *,' rb negative xx(1)=',xx(1)
                xx(1)=0.0d0

            end if

            if(abs(e(1)).lt.tol .and. abs(e(2)).lt.tol) then

c
c              reset final converged values
c
               rb=xx(1)/rhoc                                 ! cm/sec
               Ts=xx(2)                                      ! K
               xreac=2.0*alfac*R*Ts/rb/Ec                    ! 1/cm
               fr=dexp(-xreac*Ka)                            ! -
               Tf=(qr/xx(1)+Qc+Qg)/C+To                      ! K

c
c               set up initial temperature profile
c

                if (Tintype .eq. 1) then
c
c
c                  Sets ICs using I&W S.S. T profile
c

                   xcond=alfac/rb
                   xreac=xcond*R*Ts*2.0/Ec
                   fr=dexp(-xreac*Ka)
                   frJ=fr*qr/rhoc/rb/C/(Ts-To)/(1.-Ka*xcond)
                   do 3300 i=1,nx
                      Tn(i)=To+(Ts-To)*((1.0-frJ)*
     +                   dexp(x(i)/xcond)+frJ*exp(x(i)*Ka))
                      Tnp1(i)=Tn(i)
 3300              continue

                   if(abs(Tn(nx)-To)/To .gt. 0.001) then
  
                      WRITE(*,*) ' ROCBURN: rank=',rank
                      WRITE(*,*) ' ROCBURN: Error: xmax'
     +                          ,' not large enough'
                      WRITE(*,*) ' ROCBURN: To boundary'
     +                          ,' condition not satisfied'
                      STOP

                   else
  
                      Tn(nx)  =To
                      Tnp1(nx)=To
  
                   end if

                else

c
c                  Set ICs using data in Tinit.dat
c

                   open(unit=16,file='Tinit.dat',status='unknown')
                   read(16,*)
                   read(16,*)
                   read(16,*)
                   read(16,*)
                   read(16,*)
                   read(16,*)
                   do 3310 i = 1,nx
                      read(16,*) x(i), Tn(i)
                      Tnp1(i)=Tn(i)
 3310              continue

                endif

                return

            end if

 2000    continue

         PRINT *,'  ROCBURN: rank=',rank
         PRINT *,'  ROCBURN: Error-steady state ZN'
     +          ,'  in ignition_combustion'
         PRINT *,'  ROCBURN: Maximum # of ',itermax
     +          , ' iterations exceeded'
         STOP

c
c        end of steady state ZN phenomenological combustion model
c
      end if

      if(Model_combustion.eq.3) then

c
c        Model_combustion = 3
c
c        rb=a*P**n
c
c           initial guess of rb and Ts
c

         rb=a_p*P**n_p                                      ! cm/s
         Tn   =1000.0d0
         Tnp1 =1000.0d0
         Ts   =600.00d0

         
         return

      end if

      WRITE(*,*) ' ROCBURN: rank=',rank
      WRITE(*,*) ' ROCBURN: Error-ignition:'
     +          ,' No appropriate combustion model'
      STOP


      RETURN

      CONTAINS

!     ------------------------------------------------------------------------
!                             INTERNAL PROCEDURES
!     ------------------------------------------------------------------------


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

      IMPLICIT NONE

      REAL(DBL) :: x(2),J(2,2),e(2),detJ,delx(2),dum
      integer i

      detJ=J(1,1)*J(2,2)-J(2,1)*J(1,2)

      if(detJ.eq.0) then
        write(*,*) 'Error: singular matrix encountered-could not invert'
        stop
      endif
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

!     -------------------------------------------------------------------
!                       END OF INTERNAL PROCEDURES
!     -------------------------------------------------------------------

      END SUBROUTINE ignition_combustion
