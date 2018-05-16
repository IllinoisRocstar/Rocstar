c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c   subroutine  : input_combustion                                            c
c                                                                             c
c   purpose     : read inputs for burning rate model                          c
c                                                                             c
c --------------------------------------------------------------------------- c
c                                                                             c
c                            written by                                       c
c                      K.C. Tang and M.Q. Brewster                            c
c                                                                             c
c                       Last modified : 11/27/2001 by                         c
c                      Luca Massa and Robert Fiedler                          c
c --------------------------------------------------------------------------  c
c                                                                             c
c                                                                             c
c   arguments   :                                                             c
c                                                                             c
c      qr       : raditiave heat flux (cal/cm^2-s)                            c
c      To       : propellant temperature deep into the propellant (K)         c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE input_combustion(To,Ts,rb)

      USE burn_global_data
      IMPLICIT NONE

      REAL(DBL), INTENT(OUT) :: To,Ts,rb
      INTEGER                :: ierror

!    
! ??     Model_combustion = 1  :: WSB homogeneouse propellant 
!                                 combustion model, no ignition
!        Model_combustion = 2  :: ZN phenomenological combustion model 
!                                 for composite and homogeneous 
!                                 propellant, no ignition
!        Model_combustion = 3  :: rb=a*P**n, no ignition
!    
! ??     Model_combustion = 4  :: WSB homogeneouse propellant 
!                                 combustion model, with ignition
!        Model_combustion = 5  :: ZN phenomenological combustion model 
!                                 for composite and homogeneous 
!                                 propellant, with ignition
!        Model_combustion = 6  :: rb=a*P**n, with ignition

c
c     read propellant thermophysical properties
c
      CALL read_properties

c
c     input propellant initial temperature (to be supplied by ROCSOLID)
c

!RAF 11/20/01
!RAF      To       = 300.0                                 ! K
!RAF
      a_T      = 1.0
      n_T      = 1.0
!RAF 11/27/01
      Tn = To 
      Tnp1 = Tn
      Ts = To 
      rb = Ac * exp ( - ec_ru/ Ts)
!RAF
c
c     input radiation heat flux (to be supplied by ROCFLOW)            
c

!RAF      qr = 0.000                                       ! (cal/(cm^2-s))


c
c     grid generation
c

!RAF 11/27/01
!RAF Use a parameter to set the number of grid points.
!RAF      nx=101
      nx = nxmax
!RAF

      CALL grid(igrid,nx)
!      delt_max=2.0d-6
      IF(rank == 0) THEN
         PRINT *,' ROCBURN: rank=',rank,' delt_max=',delt_max
         PRINT *,' ROCBURN: Model_combustion=',Model_combustion
         SELECT CASE(Model_combustion)
           CASE(1)
              PRINT *, ' ROCBURN: rank=',rank,', Model_cobmustion = ',
     +                  Model_combustion,    
     +                 '; burn_update uses WSB dynamic model'
           CASE(2)
              PRINT *, ' ROCBURN: rank=',rank,', Model_cobmustion = ',
     +                  Model_combustion,  
     +                 '; burn_update uses ZN  dynamic model'
           CASE(3)
              PRINT *, ' ROCBURN: rank=',rank,', Model_cobmustion = ',
     +                  Model_combustion,  
     +                 '; burn_update uses rb=a*P**n   model'
           CASE(4)
              PRINT *, ' ROCBURN: rank=',rank,', Model_cobmustion = ',
     +                  Model_combustion,    
     +                 '; burn_update uses WSB dynamic model and ',
     +                 'ignition'
           CASE(5)
              PRINT *, ' ROCBURN: rank=',rank,', Model_cobmustion = ',
     +                  Model_combustion,  
     +                 '; burn_update uses ZN  dynamic model and ',
     +                 'ignition'
           CASE(6)
              PRINT *, ' ROCBURN: rank=',rank,', Model_cobmustion = ',
     +                  Model_combustion,  
     +                 '; burn_update uses rb=a*P**n   model and ',
     +                 'ignition'
           CASE DEFAULT
              PRINT *, ' ROCBURN: rank=',rank,', Model_cobmustion = ',
     +                  Model_combustion,  
     +                 '; No combustion model available'
              PRINT *, ' ROCBURN: rank=',rank,' Job Aborted'
              CALL MPI_Abort(MPI_COMM_ROCBURN,1,ierror)
              STOP
         END SELECT
      END IF

      RETURN 

      CONTAINS

!     ------------------------------------------------------------------------
!                             INTERNAL PROCEDURES
!     ------------------------------------------------------------------------

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c   subroutine  : read_properties                                             c
c                                                                             c
c   purpose     : input the propellant thermophysical properties              c
c                                                                             c
c --------------------------------------------------------------------------- c
c                                                                             c
c   arguments   :                                                             c
c                                                                             c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE read_properties
 
      USE burn_global_data

      IMPLICIT NONE

      integer ir
      REAL(DBL) :: inpmin,inpmax

!RAF 5/1/01
      integer ioerr
!RAF 5/1/01


      ir = 10

      OPEN (unit=ir,file='Rocburn/flame_properties.dat',status='old')

      READ(ir,*)   Model_combustion
      IF (rank==0) PRINT *,' ROCBURN: rank=',rank,
     +        ' ; input propellant thermophysical properties',
     +  Model_combustion
      READ(ir,*)   a_p
      IF (rank==0) PRINT *, ' ROCBURN: Model_combustion =',
     +  Model_combustion,a_p
      IF (rank==0) WRITE(*,*) 'ROCBURN:  a_p =', a_p
      READ(ir,*)   n_p
      IF (rank==0) WRITE(*,*) 'ROCBURN:  n_p =', n_p
      READ(ir,*)   Pref
      IF (rank==0) WRITE(*,*) 'ROCBURN:  Pref =', Pref
      READ(ir,*)   Ac
      IF (rank==0) WRITE(*,*) 'ROCBURN:  Ac= ',Ac
      READ(ir,*)   eg_ru
      IF (rank==0) WRITE(*,*) 'ROCBURN:  eg_ru= ',eg_ru
      READ(ir,*)   ec_ru
      IF (rank==0) WRITE(*,*) 'ROCBURN:  ec_ru ',ec_ru
      READ(ir,*)   Qc
      IF (rank==0) WRITE(*,*) 'ROCBURN:  Qc= ',Qc
      READ(ir,*)   Qg
      IF (rank==0) WRITE(*,*) 'ROCBURN:  Qg= ',Qg
      READ(ir,*)   alfac
      IF (rank==0) WRITE(*,*) 'ROCBURN:  alfac= ',alfac
      READ(ir,*)   C
      IF (rank==0) WRITE(*,*) 'ROCBURN:  C= ',C
      READ(ir,*)   rhoc
      IF (rank==0) WRITE(*,*) 'ROCBURN:  rhoc= ',rhoc
      rhoc_mks = rhoc*1000.0D0      ! from gm/(cm^3) to Kg/(m^3)
      READ(ir,*)   lamg
      IF (rank==0) WRITE(*,*) 'ROCBURN:  lamg= ',lamg
      READ(ir,*)   MW
      IF (rank==0) WRITE(*,*) 'ROCBURN:  MW= ',MW
      READ(ir,*)   Ka
      IF (rank==0) WRITE(*,*) 'ROCBURN:  Ka= ',Ka
      READ(ir,*)   delt
      IF (rank==0) WRITE(*,*) 'ROCBURN:  delt= ',delt
      READ(ir,*)   igrid
      IF (rank==0) WRITE(*,*) 'ROCBURN:  igrid= ',igrid
      READ(ir,*)   xmax
      IF (rank==0) WRITE(*,*) 'ROCBURN:  xmax= ',xmax
      READ(ir,*)   beta
      IF (rank==0) WRITE(*,*) 'ROCBURN:  beta= ',beta
      READ(ir,*)   tol_Ts
      IF (rank==0) WRITE(*,*) 'ROCBURN:  tol_Ts= ',tol_Ts
      READ(ir,*)   itermax
      IF (rank==0) WRITE(*,*) 'ROCBURN:  itermax= ',itermax
      READ(ir,*)   Tintype
      IF (rank==0) WRITE(*,*) 'ROCBURN:  Tintype= ',Tintype

!RAF 5/1/01 Read in a directory name for the restart files
!RAF        Add a trailing / if one is not present.  The
!RAF        default is the current directory.

      read(ir,'(A)',IOSTAT=ioerr) Outdir
      if (ioerr /= 0) then
        Outdir = './'
      else
        if (LEN_TRIM(Outdir) == 0) then
           Outdir = './'
        endif
        if (Outdir(LEN_TRIM(Outdir):LEN_TRIM(Outdir)) /= '/') then
           Outdir = TRIM(Outdir) // '/'
        endif
      endif

!RAF 5/1/01

      READ(ir,*)   Tstar0
      IF (rank==0) WRITE(*,*) 'ROCBURN:  Tf_adiabatic= ', Tstar0

!RAF 11/20/01
      READ(ir,*)   To
      IF (rank==0) WRITE(*,*) 'ROCBURN:  To= ',To
      READ(ir,*)   Tignition
      IF (rank==0) WRITE(*,*) 'ROCBURN:  Tignition= ',Tignition
!RAF 11/20/01
      read(ir,*,IOSTAT=ioerr) film_MKS
      if (ioerr /= 0) then
        film_MKS = 2000.0842d0  ! Default value
      endif
      IF (rank==0) WRITE(*,*) 'ROCBURN:  film_MKS= ',film_MKS
      film_coe = film_MKS / 41868.0d0

      CLOSE(ir)

      lamc=alfac*rhoc*C
   
      R=1.9872
    
!RAF      Pref = 100.0    !(atm)

!..............................................................................
! Check the input deck
!..............................................................................

      inpmin = min(a_p,Ac,Pref,n_p,ec_ru,Tstar0,eg_ru,alfac)
      inpmax = max(a_p,Ac,Pref,n_p,ec_ru,Tstar0,eg_ru,alfac)

      if( inpmin < alfac - 1.d-8 .OR. inpmin <= 1.d-4 .OR.
     &    inpmax > Ac + 1.d-8 ) THEN
         PRINT *,'ROCBURN: Rank ',rank,' ERROR IN THE INPUT DECK',inpmin,inpmax
         CALL MPI_Abort(MPI_COMM_ROCBURN,1,ioerr)
         STOP
      
      else
      
      IF (rank==0) PRINT *,'ROCBURN: rank=',rank, 
     +        ' done input thermophysical properties'

      endif

      RETURN 
      END SUBROUTINE read_properties

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c   subroutine  : grid                                                        c
c                                                                             c
c   purpose     : grid generation                                             c
c                                                                             c
c --------------------------------------------------------------------------- c
c                                                                             c
c   arguments   :                                                             c
c                                                                             c
c      gridtype : = 1  expoential grid                                        c
c                 = 2  boundary layer grid                                    c
c      numx     : maximum number of grid points                               c
c      beta     : stretch parameter                                           c
c      xmax     : the maximum dimension in x-coordinate                       c
c      x        : the x-coordinate array                                      c
c      z        : the z-coordinate arary                                      c
c      zx       : the first derivative of z                                   c
c      zxx      : the second derivative of z                                  c
c      delz     : the z spacing                                               c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE grid(gridtype,numx)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: gridtype, numx
!     REAL(DBL) :: x(:), z(:), zx(:), zxx(:)
!     REAL(DBL) :: beta, xmax, delz

      REAL(DBL) :: bpm,bp1,bm1,czx,czxx,term,term2
      INTEGER :: i



      IF(xmax.ge.0.0) then
         PRINT *,' ROCBURN: rank=',rank,' xmax >= 0'
         PRINT *,' ROCBURN: job aborted!'
         STOP
      END IF

      if(numx.gt.(nxmax)) then
         PRINT *,' ROCBURN: rankk=',rank,' numx=',numx,' nxmax=',nxmax
         PRINT *,' ROCBURN: nx > maximum number of',
     +           ' allowed spatial nodes'
         PRINT *,' ROCBURN: job aborted!'
         STOP
      end if

      delz = 1.0/dble(numx-1)

      IF(rank ==0) THEN
         PRINT *,' ROCBURN: delz=',delz
         PRINT *,' ROCBURN: nxmax=',nxmax
      END IF

      IF (gridtype.eq.1) THEN
c
c        exponential grid
c

         DO 10 i=1,numx
            z(i) = dble(i-1)*delz
            x(i) = log(1.-z(i))/beta
            IF (i.eq.numx) THEN
                    x(i) = -10.
            END IF
            zx(i)  = -beta*exp(beta*x(i))
            zxx(i) = -beta**2*exp(beta*x(i))
   10    CONTINUE

      ELSE
c
c        Anderson, Tanehill and Pletcher boundary layer grid control
c

         bp1=beta+1.0
         bm1=beta-1.0
         bpm=bp1/bm1
         czx=2.0*beta/log(bpm)/xmax
         czxx=-2.0/xmax*czx
         DO 30 i=1,numx
            z(i)=dble(i-1)*delz
            term=bpm**(1.0-z(i))
            x(i)=xmax*(bp1-bm1*term)/(term+1.0)
            term2=beta**2-(1.0-x(i)/xmax)**2
            zx(i)=czx/term2
            zxx(i)= czxx*(1.0-x(i)/xmax)/term2**2
   30    CONTINUE

      END IF

!..............................................................................................
!  EVALUATE MAXIMUM DT for the diffusion
!.............................................................................................
      delt_max = 0.5d0*delz*delz/maxval(abs(zx))**2/alfac
      delz2inv = 0.5d0/delz
      delzsqinv = 1.0d0 / (delz*delz)
      dx1 = delz/zx(1)

      IF(rank == 0) THEN
         WRITE (*,*)
         WRITE (*,*)' ROCBURN: Smallest delta x = ',x(1)-x(2)
         WRITE (*,*)' ROCBURN: Largest delta x  = ',x(numx-1)-x(numx)
         WRITE (*,*)
      END IF
      RETURN
      END SUBROUTINE grid


!     -------------------------------------------------------------------
!                       END OF INTERNAL PROCEDURES
!     -------------------------------------------------------------------

      END SUBROUTINE input_combustion
