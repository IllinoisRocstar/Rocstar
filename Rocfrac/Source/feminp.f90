!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
SUBROUTINE feminp(glb, myid)

  USE ROCSTAR_RocFrac

!!****f* Rocfrac/Source/feminp.f90
!!
!!  NAME
!!     feminp
!!
!!  FUNCTION
!!
!!     READ INPUT INFORMATION (i.e. Analysis Deck File)
!!
!!  INPUTS
!!     glb -- global array
!!     myid -- processor id (starting at 0)
!!
!!****

  IMPLICIT NONE
  INCLUDE 'mpif.h'
  TYPE(ROCFRAC_GLOBAL) :: glb

! -- local
  CHARACTER*200 :: keywd    ! keyword parameter
  
  INTEGER :: ios            ! io error
  INTEGER :: ierr           ! mpi error
  
  INTEGER :: i              ! loop counter
! -- global
  INTEGER :: myid           ! processor id
  INTEGER :: lll            ! length of keyword, with no trailing blanks

  REAL*8, ALLOCATABLE, DIMENSION(:) :: tmp_E, tmp_xnu, tmp_rho, tmp_alpha
  INTEGER, ALLOCATABLE,DIMENSION(:) :: tmp_iSolnType
  REAL*8, ALLOCATABLE, DIMENSION(:) :: tmp_Cp, tmp_KappaHT
  INTEGER, ALLOCATABLE,DIMENSION(:) :: tmp_iSolnTypeHT

  ALLOCATE(tmp_E(1:10), tmp_xnu(1:10), tmp_rho(1:10), tmp_alpha(1:10), tmp_iSolnType(1:10) )
  ALLOCATE(tmp_Cp(1:10), tmp_KappaHT(1:10),tmp_iSolnTypeHT(1:10))
!
! -- Set Default values
  
!  IF(myid.EQ.0) PRINT*,'ROCFRAC:: INSIDE feminp.f90'
  
  glb%DummyTractVal = 0.d0
  glb%DummyBurnRate = 0.d0
  
  glb%NumNodeIO = 0
  glb%ALEenabled = .FALSE.       ! ALE disabled by default
  glb%ipstatic = .FALSE.     ! don't subract the intial pressure
  glb%ReStart = .FALSE.
  glb%iElType = 4 ! default to 4 node tets
  glb%IONEWER = .FALSE.
  glb%DampEnabled = .FALSE.
  glb%iElIntgratn = 0
  glb%NdMassLump = 0
  glb%EnforceTractionS = .FALSE.! To enforce traction where no fluids pressure
  glb%NumEntries = 0
  glb%NumMatVol = 0
  glb%NumMatVolHT = 0
  glb%NumMatCoh = 0
  glb%cd_fastest = 0.d0
  glb%NdBasedEl = .FALSE.
  glb%UnDefConfig = .FALSE.
  glb%HeatTransSoln = .FALSE.
  glb%Temperature0 = 0.
  glb%ArtificialDamping = .FALSE.
  glb%EnforceTractionS = .FALSE.
  glb%EnforceTractionSF = .FALSE.
  glb%DebondPart = .FALSE.
  glb%DebondPart_MATOUS = .FALSE.
  glb%ThermalExpansion = .FALSE.
  glb%AmplitudeTable = .FALSE.
  glb%debug_state = .FALSE.
  glb%NumProbesEl = 0
  glb%NumProbesNd = 0
  glb%OverlayExist = .false.
  glb%DummyFlux = 0.d0

  glb%NumElOverlay = 0
  glb%NumNpOverlay = 0
  glb%Verb = 1
!
! -- Open Analysis Deck File

  OPEN(glb%io_input,FILE='./Rocfrac/RocfracControl.txt',STATUS='old',IOSTAT=ios)
  IF(ios.NE.0)THEN
     IF(myid.EQ.0) PRINT*, 'ROCFRAC:: Unable to find RocfracControl.txt - STOPPING'
     CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)
     STOP
  ENDIF
!
! Input deck summary file

  IF(myid.EQ.0)THEN
     OPEN(glb%io_sum,FILE='Rocfrac/Modin/InputSummary.res',STATUS='unknown',IOSTAT=ios)
     IF(ios.NE.0)THEN ! Try without Modin
        OPEN(glb%io_sum,FILE='Rocfrac/InputSummary.res',STATUS='unknown',IOSTAT=ios)
     END IF

     IF(ios.NE.0)THEN
        IF(myid.EQ.0) PRINT*, 'ROCFRAC:: Unable to find InputSummary.res under Rocfrac/Modin/ or Rocfrac/ - STOPPING'
        CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)
        STOP
     ENDIF
  ENDIF

!
! -- Read Analysis Deck File

  REWIND glb%io_input
1 READ(glb%io_input,'(A)',IOSTAT=ios) keywd
!  if(myid.EQ.0) print*,keywd
  IF(ios.LT.0) THEN ! Negative ios means end-of-file
     PRINT*,' *END parameter not found - STOPPING'
     CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)
     STOP
  ENDIF

2 CONTINUE

!
! Comment field

  IF(keywd(1:1).NE.'*')THEN
     GOTO 1
  ENDIF

  IF(keywd(1:2).EQ.'**')THEN
     GOTO 1
  ENDIF

  lll = LEN_TRIM(keywd)

  IF(myid.EQ.0) WRITE(glb%io_sum,'(2X,A,A)') keywd(1:lll)
  
  IF(keywd(1:4).EQ.'*END') THEN
     GOTO 3

! Amplitude

  ELSE IF(keywd(1:10).EQ.'*AMPLITUDE') THEN
     CALL AMPLITUDE_SUB(glb,keywd)
     GOTO 1

  ELSE IF(keywd(1:13).EQ.'*HYPERELASTIC')THEN
     CALL MATMODEL_HYPERELASTIC(glb, keywd, &
          tmp_E, tmp_xnu, tmp_rho, tmp_alpha, tmp_iSolnType)
     GOTO 1

  ELSE IF(keywd(1:14).EQ.'*HEAT TRANSFER')THEN
     CALL HEAT_TRANSFER_SUB(glb, keywd, tmp_KappaHT, tmp_Cp, tmp_iSolnTypeHT)
     GOTO 1

  ELSE IF(keywd(1:8).EQ.'*ELASTIC')THEN
     CALL MATMODEL_ELASTIC(glb, keywd, &
          tmp_E, tmp_xnu, tmp_rho, tmp_alpha, tmp_iSolnType)
     GOTO 1

  ELSE IF(keywd(1:8).EQ.'*ELEMENT')THEN
     CALL ELEMENT_SUB(glb,keywd)
     GOTO 1
  ELSE IF(keywd(1:8).EQ.'*DAMPING') THEN
     READ(glb%io_input,*) glb%KappaDamp
     glb%DampEnabled = .TRUE.
     GOTO 1
!
!  Prefix for input/output files
!
  ELSE IF(keywd(1:7).EQ.'*PREFIX') THEN
     CALL PREFIX_SUB(glb)
     GOTO 1
  ELSE IF(keywd(1:4).EQ.'*ALE') THEN
     CALL ALE_SUB(glb)
     GOTO 1
  ELSE IF(keywd(1:5).EQ.'*NRUN') THEN
     CALL NRUN_SUB(glb)
     GOTO 1
  ELSE IF(keywd(1:8).EQ.'*DYNAMIC') THEN
     CALL DYNAMIC_SUB(glb,keywd)
     GOTO 1
!*** OBSOLETE REMOVE
  ELSE IF(keywd(1:7).EQ.'*MATVOL') THEN ! Read volumetric material props.
     CALL MATVOL_SUB(glb, &
          tmp_E, tmp_xnu, tmp_rho, tmp_alpha, tmp_iSolnType)
     GOTO 1
  ELSE IF(keywd(1:7).EQ.'*MATCOH') THEN ! Read cohesive material props
     CALL MATCOH_SUB(glb)
     GOTO 1
  ELSE IF(keywd(1:7).EQ.'*PLOAD1') THEN
     CALL PLOAD1_SUB(glb)
     GOTO 1
  ELSE IF(keywd(1:5).EQ.'*NODE') THEN
     CALL NODEIO_SUB(myid)
     GOTO 1
  ELSE IF(keywd(1:11).EQ.'*DUMMYTRACT') THEN
     CALL DUMMYTRACT_SUB(glb,keywd)
     GOTO 1
  ELSE IF(keywd(1:10).EQ.'*DUMMYBURN') THEN
     READ(glb%io_input,*) glb%DummyBurnRate
     GOTO 1
  ELSE IF(keywd(1:10).EQ.'*DUMMYFLUX') THEN
     READ(glb%io_input,*) glb%DummyFlux
     GOTO 1
  ELSE IF(keywd(1:11).EQ.'*BOUNDARYMM') THEN
     CALL BOUNDARYMM_SUB(glb)
     GOTO 1
  ELSE IF(keywd(1:11).EQ.'*BOUNDARYHT') THEN
     CALL BOUNDARYHT_SUB(glb)
     GOTO 1
  ELSE IF(keywd(1:9).EQ.'*BOUNDARY') THEN
     CALL BOUNDARY_SUB(glb)
     GOTO 1
  ELSE IF(keywd(1:9).EQ.'*IPSTATIC')THEN ! to subract out the initial pressure
     glb%ipstatic = .TRUE.
     GOTO 1
  ELSE IF(keywd(1:8).EQ.'*IONEWER')THEN ! new gen 2.5 input format
     glb%IONEWER = .TRUE.
     GOTO 1
  ELSE IF(keywd(1:10).EQ.'*DEFCONFIG')THEN
     glb%UnDefConfig = .TRUE.
     GOTO 1

  ELSE IF(keywd(1:16).EQ.'*MICROMECHANICAL') THEN ! micromechanical model
     CALL MICROMECHANICAL_SUB(glb, keywd, &
          tmp_E, tmp_xnu, tmp_rho, tmp_alpha, tmp_iSolnType)
     GOTO 1

  ELSE IF(keywd(1:11).EQ.'*ARTDAMPING')THEN
     glb%ArtificialDamping = .TRUE.
     GOTO 1
  ELSE IF(keywd(1:6).EQ.'*PROBE')THEN
     CALL PROBE_SUB(glb)
     GOTO 1
  ELSE IF(keywd(1:18).EQ.'*INITIAL CONDITION')THEN
     CALL INITIALCONDITION_SUB(glb,keywd)
     GOTO 1
  ELSE IF(keywd(1: ).EQ.'*VERBOSITY')THEN
     READ(glb%io_input,*)  glb%Verb
     GOTO 1
  ELSE IF(keywd(1:6).EQ.'*DEBUG')THEN
     glb%debug_state = .TRUE.
     GOTO 1

  ELSE IF(keywd(1:1).EQ.'*')THEN
     PRINT*,'ROCFRAC: ERROR'
     PRINT*,'ROCFRAC: CONTROL DECK OPTION ',TRIM(keywd),' NOT SUPPORTED'
     PRINT*,'ROCFRAC: STOPPING'
     CALL MPI_BARRIER(glb%MPI_COMM_ROCFRAC,ierr)
     CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)
     STOP
  ELSE
     GOTO 1
  ENDIF
3 CONTINUE

  CLOSE(glb%io_input)
!
!---- Write analysis summary file

  ALLOCATE(glb%E       (1:glb%NumMatVol) )
  ALLOCATE(glb%xnu     (1:glb%NumMatVol) )
  ALLOCATE(glb%rho     (1:glb%NumMatVol) )
  ALLOCATE(glb%alpha   (1:glb%NumMatVol) )
  ALLOCATE(glb%xmu     (1:glb%NumMatVol) )
  ALLOCATE(glb%xlambda (1:glb%NumMatVol) )
  ALLOCATE(glb%xkappa  (1:glb%NumMatVol) )

  ALLOCATE(glb%iSolnType(1:glb%NumMatVol) )


  IF(glb%DebondPart_Matous) glb%rho(1) = tmp_rho(1)

  DO i = 1, glb%NumMatVol
     glb%E(i) = tmp_E(i)
     glb%xnu(i) = tmp_xnu(i)
     glb%rho(i) = tmp_rho(i)
     glb%alpha(i) = tmp_alpha(i)
     glb%iSolnType(i) = tmp_iSolnType(i)

     glb%xmu(i) = glb%E(i)/(2.d0*(1.d0+glb%xnu(i)))
    
     glb%xlambda(i) = 2.d0*glb%xmu(i)*glb%xnu(i)/(1.d0-2.d0*glb%xnu(i))

! bulk modulus
     glb%xkappa(i) = glb%xlambda(i) + 2.d0/3.d0*glb%xmu(i)

  ENDDO
    
  IF(glb%HeatTransSoln)THEN
     glb%ThermalDiffusivity = -1.

     ALLOCATE(glb%iSolnTypeHT(1:glb%NumMatVolHT) )
     ALLOCATE(glb%KappaHT (1:glb%NumMatVolHT) )
     ALLOCATE(glb%Cp      (1:glb%NumMatVolHT) )

     DO i = 1, glb%NumMatVolHT
        glb%KappaHT(i) = tmp_KappaHT(i)
        glb%Cp(i) = tmp_Cp(i)
        glb%iSolnTypeHT(i) = tmp_iSolnTypeHT(i)

        

        glb%ThermalDiffusivity = MAX( glb%ThermalDiffusivity, glb%KappaHT(i)/(glb%rho(i)*glb%Cp(i)) )


     ENDDO
  ENDIF

  DEALLOCATE(tmp_KappaHT, tmp_Cp, tmp_iSolnTypeHT)
  DEALLOCATE(tmp_E,tmp_xnu,tmp_rho,tmp_alpha,tmp_iSolnType)

  IF(myid.EQ.0)THEN
     WRITE(glb%io_sum,50) glb%prefx
     WRITE(glb%io_sum,60) 
     IF (glb%restart) WRITE(glb%io_sum,65)
     WRITE(glb%io_sum,80) glb%NumMatVol,glb%NumMatCoh
     WRITE(glb%io_sum,85) glb%CourantRatio
     WRITE(glb%io_sum,86)
     DO i = 1,glb%NumMatCoh
        WRITE(glb%io_sum,88) glb%deltan(i),glb%deltat(i),glb%SigmaMax(i),glb%TauMax(i),glb%Sinit(i)


     ENDDO
     WRITE(glb%io_sum,90)
     WRITE(glb%io_sum,91)
     DO i = 1,glb%NumMatVol
        WRITE(glb%io_sum,92) i,glb%E(i),glb%xnu(i),glb%rho(i),glb%alpha(i)
        IF(glb%iSolnType(i).EQ.0)THEN
           WRITE(glb%io_sum,*) '          Material Model = Arruda-Boyce '
        ELSE IF(glb%iSolnType(i).EQ.-1)THEN
           WRITE(glb%io_sum,*) '          Material Model = NeoHookean Incompressible '
        ELSE IF(glb%iSolnType(i).EQ.1)THEN
           WRITE(glb%io_sum,*) '          Material Model = Elastic, Large Deformation '
        ELSE IF(glb%iSolnType(i).EQ.2)THEN
           WRITE(glb%io_sum,*) '          Material Model = Elastic, Small Deformation '
        ELSE
           WRITE(glb%io_sum,*) ' Error, Not a valid Material model, = ', glb%iSolnType(i)
        ENDIF
           
     ENDDO

     IF(glb%HeatTransSoln) WRITE(glb%io_sum,*) 'Heat Transfer Solution'
     CLOSE(glb%io_sum)
  ENDIF
  
  RETURN
   
!--------------------------------FORMATS--------------------------------
!

50 FORMAT(//,'DYNAMIC 3D LINEAR ELASTIC FEA',///,'Job id: ',a20,//)
60 FORMAT('*** ISOTROPIC ANALYSIS ***',/)
65 FORMAT(' ***    RESTART   ***',/)
80 FORMAT(1x,'Number of material types (NUMAT_VOL)        =',i12 &
        /,1x,'Number of material types (NUMAT_COH)        =',i12)
85 FORMAT(/,1x,'Steps per characteristic length (STEPS)     =',e12.4)
86 FORMAT(///,1x,'COHESIVE ELEMENT DATA')
88 FORMAT(/,1x,'Characteristic lengths for the cohesive law:', &
        /,1x,'     normal (DELTAN)                        =',e12.4, &
        /,1x,'     tangential (DELTAT)                    =',e12.4, &
        /,1x,'Maximum normal stress (GLB%SIGMAMAX)              =',e12.4, &
        /,1x,'Maximum shearing stress (TAUMAX)            =',e12.4, &
        /,1x,'The initial Sthreshold S(init)              =',e12.4)
90 FORMAT(///,1x,'VOLUMETRIC ELEMENT DATA') 
91 FORMAT(/,4x,'MATERIAL SETS')
92 FORMAT(/,9x,'SET',4x,i4 &
        /,12x,'E             =',10x,e13.5, &
        /,12x,'Nu            =',10x,e13.5, &
        /,12x,'rho           =',10x,e13.5, &
        /,12x,'alpha         =',10x,e13.5)
  
END SUBROUTINE feminp


         
SUBROUTINE PREFIX_SUB(glb)

!!****f* Rocfrac/Source/feminp/PREFIX_SUB
!!
!!  NAME
!!     PREFIX_SUB
!!
!!  FUNCTION
!!
!!     Reads prefix keyword (i.e. Analysis Deck File)
!!
!!  INPUTS
!!     glb -- global array
!!
!!****
  
  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  
  TYPE(ROCFRAC_GLOBAL) :: glb
  
  READ(glb%io_input,'(A)') glb%prefx
  glb%prefx_lngth = LEN_TRIM(glb%prefx)
  
  RETURN
END SUBROUTINE PREFIX_SUB

SUBROUTINE ALE_SUB(glb)

!!****f* Rocfrac/Source/feminp/ALE_SUB
!!
!!  NAME
!!     ALE_SUB
!!
!!  FUNCTION
!!
!!     ALE keyword turns on ALE routines
!!
!!  INPUTS
!!     glb -- global array
!!
!!****
  
  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  
  TYPE(ROCFRAC_GLOBAL) :: glb

! kappa parameter for mesh motion

  READ(glb%io_input,*) glb%kappa
  glb%ALEenabled = .TRUE.
  
  !print*,'Running with ALE'
  
  RETURN
END SUBROUTINE ALE_SUB

SUBROUTINE NRUN_SUB(glb)
  
  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  
  TYPE(ROCFRAC_GLOBAL) :: glb
  
  READ(glb%io_input,*) glb%CourantRatio

!RAF  Make this backward compatible with old input decks --
!RAF  You would not want a ratio > 1, right?

  IF (glb%CourantRatio > 1.0d0) THEN
    glb%CourantRatio = 1.0d0 / glb%CourantRatio
  ENDIF
  
  RETURN
END SUBROUTINE NRUN_SUB

SUBROUTINE DYNAMIC_SUB(glb,keywd)

!!****f* Rocfrac/Source/feminp/DYNAMIC_SUB
!!
!!  NAME
!!     DYNAMIC_SUB
!!
!!  FUNCTION
!!
!!     Courant limit multiplier
!!
!!  INPUTS
!!     glb -- global array
!!     keywd -- keywd for control deck
!!
!!****
  
  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE

  CHARACTER(len=200) :: keywd
  CHARACTER(len=16) :: ScaleFactor
  INTEGER :: k1, k2
  
  
  TYPE(ROCFRAC_GLOBAL) :: glb

  CALL locchr(keywd,'SCALE FACTOR                ',12,8,k1,k2)

  ScaleFactor = keywd(k1:k2)

  CALL rchar(ScaleFactor,glb%CourantRatio)
  
  RETURN
END SUBROUTINE DYNAMIC_SUB

SUBROUTINE DUMMYTRACT_SUB(glb,keywd)

!!****f* Rocfrac/Source/feminp/DUMMYTRACT_SUB
!!
!!  NAME
!!    DUMMYTRACT_SUB
!!
!!  FUNCTION
!!    Marks if applying the dummy traction to what surface 
!!     
!!
!!  INPUTS
!!     glb -- global array
!!     keywd -- keywd for control deck
!!
!!****
  
  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE

  CHARACTER(len=200) :: keywd
  INTEGER :: k1, k2
  CHARACTER :: Option*16
  
  TYPE(ROCFRAC_GLOBAL) :: glb

  CALL locchr(keywd,'INTERFACE                ',9,8,k1,k2)

  Option = keywd(k1:k2)

  IF(Option.EQ.'S') glb%EnforceTractionS = .true.
  IF(Option.EQ.'SF') glb%EnforceTractionSF = .true.
  IF(Option.EQ.'ALL')THEN
     glb%EnforceTractionS = .true.
     glb%EnforceTractionSF = .true.
  ENDIF

  READ(glb%io_input,*) glb%DummyTractVal
  
  RETURN
END SUBROUTINE DUMMYTRACT_SUB

SUBROUTINE MATMODEL_HYPERELASTIC(glb,keywd, &
     tmp_E, tmp_xnu, tmp_rho, tmp_alpha, tmp_iSolnType)

!!****f* Rocfrac/Source/feminp/MATMODEL_HYPERELASTIC
!!
!!  NAME
!!     MATMODEL_HYPERELASTIC
!!
!!  FUNCTION
!!     reads material model type, and the 
!!     hyperelastic material parameters
!!     
!!
!!  INPUTS
!!     glb -- global array
!!     keywd -- keywd for control deck
!!  
!!  OUTPUT
!!     tmp_E -- Young's Modulus
!!     tmp_xnu -- Possion's ratio
!!     tmp_rho -- Density
!!     tmp_alpha -- thermal coefficient of expansion
!!     tmp_iSolnType -- material type model 
!!
!!****

  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE

  TYPE(ROCFRAC_GLOBAL) :: glb
  
  CHARACTER :: keywd*200
  INTEGER :: key ! 0 = no key word found, 1 = key word found)
  CHARACTER*26 :: MatType
  INTEGER :: NumMatType
  INTEGER :: i,ii,ierr
  REAL*8, DIMENSION(1:10) :: tmp_E, tmp_xnu, tmp_rho, tmp_alpha
  INTEGER, DIMENSION(1:10) :: tmp_iSolnType

  CALL CONCHR(keywd,'NEOHOOKINC                ',10,13,key)
  IF(key.EQ.1) MatType = 'NEOHOOKINC'

  CALL CONCHR(keywd,'ARRUDA-BOYCE              ',12,13,key)
  IF(key.EQ.1) MatType = 'ARRUDA-BOYCE'

  SELECT CASE (TRIM(MatType))

  CASE ('ARRUDA-BOYCE')

     READ(glb%io_input,*) NumMatType

     DO i = 1, NumMatType
        glb%NumMatVol = glb%NumMatVol + 1

        ii = glb%NumMatVol

        IF(ii.GT.10)THEN
           PRINT*,'ROCFRAC :: ERROR'
           PRINT*,'Number of materials GREATER then 10'
           CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)
        ENDIF

        READ(glb%io_input,*) tmp_E(ii), tmp_xnu(ii), tmp_rho(ii), tmp_alpha(ii)

        glb%cd_fastest = MAX( glb%cd_fastest, &
             SQRT(tmp_E(ii)*(1.d0-tmp_xnu(ii))/tmp_rho(ii)/(1.d0+tmp_xnu(ii))/(1.d0-2.d0*tmp_xnu(ii) )) )
        tmp_iSolnType(ii) = 0

     ENDDO
     
  CASE ('NEOHOOKINC')

     READ(glb%io_input,*) NumMatType

     DO i = 1, NumMatType
        glb%NumMatVol = glb%NumMatVol + 1

        ii = glb%NumMatVol

        IF(ii.GT.10)THEN
           PRINT*,'ROCFRAC :: ERROR'
           PRINT*,'Number of materials GREATER then 10'
           CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)
        ENDIF

        READ(glb%io_input,*) tmp_E(ii), tmp_xnu(ii), tmp_rho(ii), tmp_alpha(ii)

        glb%cd_fastest = MAX( glb%cd_fastest, &
             SQRT(tmp_E(ii)*(1.d0-tmp_xnu(ii))/tmp_rho(ii)/(1.d0+tmp_xnu(ii))/(1.d0-2.d0*tmp_xnu(ii) )) )
        tmp_iSolnType(ii) = -1

     ENDDO

  CASE default

     PRINT*,'ROCFRAC :: ERROR'
     PRINT*,'ROCFRAC :: *HYPERELASTIC KEYWORD ',TRIM(MatType), ' NOT FOUND'
     PRINT*,'ROCFRAC :: STOPPING'
     CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)
     
     
  END SELECT


END SUBROUTINE MATMODEL_HYPERELASTIC

SUBROUTINE MATMODEL_ELASTIC(glb,keywd, &
     tmp_E, tmp_xnu, tmp_rho, tmp_alpha, tmp_iSolnType)

!!****f* Rocfrac/Source/feminp/MATMODEL_ELASTIC
!!
!!  NAME
!!     MATMODEL_ELASTIC
!!
!!  FUNCTION
!!     reads material model type, and the 
!!     elastic material parameters
!!     
!!
!!  INPUTS
!!     glb -- global array
!!     keywd -- keywd for control deck
!!  
!!  OUTPUT
!!     tmp_E -- Young's Modulus
!!     tmp_xnu -- Possion's ratio
!!     tmp_rho -- Density
!!     tmp_alpha -- thermal coefficient of expansion
!!     tmp_iSolnType -- material type model 
!!
!!****

  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE

  TYPE(ROCFRAC_GLOBAL) :: glb
  
  CHARACTER :: keywd*200
  INTEGER :: key ! 0 = no key word found, 1 = key word found)
  CHARACTER*26 :: MatType
  INTEGER :: NumMatType
  INTEGER :: i,ii,ierr
  REAL*8, DIMENSION(1:10) :: tmp_E, tmp_xnu, tmp_rho, tmp_alpha
  INTEGER, DIMENSION(1:10) :: tmp_iSolnType

  INTEGER :: k1, k2
  CHARACTER :: NLGeom*16
  INTEGER :: NLGeomType

  CALL locchr(keywd,'NLGEOM                    ',6,8,k1,k2)

  NLGeom = keywd(k1:k2)

  NLGeomType = 1            ! defaults to nl 
  IF(NLGeom.EQ.'NO')THEN
     NLGeomType = 2
  ENDIF


  READ(glb%io_input,*) NumMatType

  DO i = 1, NumMatType

     glb%NumMatVol = glb%NumMatVol + 1

     ii = glb%NumMatVol
     
     IF(ii.GT.10)THEN
        PRINT*,'ROCFRAC :: ERROR'
        PRINT*,'Number of ELASTIC materials GREATER then 10'
        CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)
     ENDIF
     
     READ(glb%io_input,*) tmp_E(ii), tmp_xnu(ii), tmp_rho(ii), tmp_alpha(ii)

     IF(tmp_alpha(ii).NE.0.d0) glb%ThermalExpansion = .TRUE.

     glb%cd_fastest = MAX( glb%cd_fastest, &
          SQRT(tmp_E(ii)*(1.d0-tmp_xnu(ii))/tmp_rho(ii)/(1.d0+tmp_xnu(ii))/(1.d0-2.d0*tmp_xnu(ii) )) )
     tmp_iSolnType(ii) = NLGeomType
     
  ENDDO

END SUBROUTINE MATMODEL_ELASTIC

!*** OBSOLETE REMOVE 

SUBROUTINE MATVOL_SUB(glb, &
     tmp_E, tmp_xnu, tmp_rho, tmp_alpha, tmp_iSolnType)
  
  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  
  TYPE(ROCFRAC_GLOBAL) :: glb
  
  INTEGER :: i,ii,ierr              ! loop counter
  integer :: NumMatType
  REAL*8, DIMENSION(1:10) :: tmp_E, tmp_xnu, tmp_rho, tmp_alpha
  INTEGER, DIMENSION(1:10) :: tmp_iSolnType

  READ(glb%io_input,*) NumMatType

  DO i = 1, NumMatType
     glb%NumMatVol = glb%NumMatVol + 1
     
     ii = glb%NumMatVol
     
     IF(ii.GT.10)THEN
        PRINT*,'ROCFRAC :: ERROR'
        PRINT*,'Number of materials GREATER then 10'
        CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)
        STOP
     ENDIF
     
     READ(glb%io_input,*) tmp_E(ii), tmp_xnu(ii), tmp_rho(ii), tmp_alpha(ii),tmp_iSolnType(ii)
     

     glb%cd_fastest = MAX( glb%cd_fastest, &
          SQRT(tmp_E(ii)*(1.d0-tmp_xnu(ii))/tmp_rho(ii)/(1.d0+tmp_xnu(ii))/(1.d0-2.d0*tmp_xnu(ii) )) )
  ENDDO
  RETURN
END SUBROUTINE MATVOL_SUB

SUBROUTINE MATCOH_SUB(glb)
  
  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  
  TYPE(ROCFRAC_GLOBAL) :: glb
  
  INTEGER :: i              ! loop counter
  
  READ(glb%io_input,*) glb%NumMatCoh
  
  ALLOCATE(glb%deltan(1:glb%NumMatCoh))
  ALLOCATE(glb%deltat(1:glb%NumMatCoh))
  ALLOCATE(glb%SigmaMax(1:glb%NumMatCoh))
  ALLOCATE(glb%TauMax(1:glb%NumMatCoh))
  ALLOCATE(glb%Sinit (1:glb%NumMatCoh))
  
  DO i = 1, glb%NumMatCoh
     READ(glb%io_input,*) glb%deltan(i),glb%deltat(i),glb%SigmaMax(i),glb%TauMax(i), glb%Sinit(i)
  ENDDO

  RETURN
END SUBROUTINE MATCOH_SUB

SUBROUTINE PLOAD1_SUB(glb)
  
  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  
  TYPE(ROCFRAC_GLOBAL):: glb
  
  INTEGER :: pload1_input ! fix not donei              ! loop counter
  
  READ(glb%io_input,*) pload1_input
  
  RETURN
END SUBROUTINE PLOAD1_SUB

SUBROUTINE NODEIO_SUB(myid)
  
  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  
  TYPE(ROCFRAC_GLOBAL) :: glb
  
  INTEGER :: i              ! loop counter
  INTEGER :: myid
  
  READ(glb%io_input,*) glb%NumNodeIO, glb%NumNodeIOpid
  
  IF(glb%NumNodeIOpid.EQ.myid)THEN
     
     ALLOCATE(glb%NodeIO(1:glb%NumNodeIO))
     
     DO i = 1, glb%NumNodeIO
        READ(glb%io_input,*) glb%NodeIO
     ENDDO
     
  ELSE
     
     glb%NumNodeIO = 0
     
  ENDIF
  
  RETURN
END SUBROUTINE NODEIO_SUB


SUBROUTINE BOUNDARY_SUB(glb)

!!****f* Rocfrac/Source/feminp/BOUNDARY_SUB
!!
!!  NAME
!!     BOUNDARY_SUB
!!
!!  FUNCTION
!!     This option is used to prescibe boundary conditions 
!!     at nodes.
!!     
!!
!!  INPUTS
!!     glb -- global array
!!
!!****

  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  
  TYPE(ROCFRAC_GLOBAL) :: glb
  
  INTEGER :: i,iaux              ! loop counter
  INTEGER :: NumBcFlags
  
  READ(glb%io_input,*) NumBcFlags
  
  ALLOCATE(glb%bcCond(1:NumBcFlags))
  
  DO i = 1, NumBcFlags
     READ(glb%io_input,*)  iaux,glb%bcCond(i)%BCtypeX,glb%bcCond(i)%BCtypeY,glb%bcCond(i)%BCtypeZ, &
          glb%bcCond(i)%BCvalueX,glb%bcCond(i)%BCvalueY,glb%bcCond(i)%BCvalueZ
  ENDDO
  
  RETURN
END SUBROUTINE BOUNDARY_SUB

SUBROUTINE BOUNDARYHT_SUB(glb)

!!****f* Rocfrac/Source/feminp/BOUNDARYHT_SUB
!!
!!  NAME
!!     BOUNDARYHT_SUB
!!
!!  FUNCTION
!!     This option is used to prescibe heat transfer
!!     boundary conditions at nodes.
!!     
!!
!!  INPUTS
!!     glb -- global array
!!
!!****

  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  
  TYPE(ROCFRAC_GLOBAL) :: glb
  
  INTEGER :: i,iaux              ! loop counter
  INTEGER :: NumBcFlagsHT
  
  READ(glb%io_input,*) NumBcFlagsHT

  
  ALLOCATE(glb%bcCondHT(1:NumBcFlagsHT))
  
  DO i = 1, NumBcFlagsHT
     READ(glb%io_input,*)  iaux,glb%bcCondHT(i)%BCtypeX,glb%bcCondHT(i)%BCtypeY,glb%bcCondHT(i)%BCtypeZ, &
          glb%bcCondHT(i)%BCvalueX,glb%bcCondHT(i)%BCvalueY,glb%bcCondHT(i)%BCvalueZ
  ENDDO
  
  RETURN
END SUBROUTINE BOUNDARYHT_SUB

SUBROUTINE BOUNDARYMM_SUB(glb)

!!****f* Rocfrac/Source/feminp/BOUNDARYMM_SUB
!!
!!  NAME
!!     BOUNDARYMM_SUB
!!
!!  FUNCTION
!!     This option is used to prescibe mesh motion boundary 
!!     conditions at nodes.
!!     
!!  INPUTS
!!     glb -- global array
!!
!!****

  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  
  TYPE(ROCFRAC_GLOBAL) :: glb
  
  INTEGER :: i,iaux              ! loop counter
  INTEGER :: NumBcFlagsmm
  
  READ(glb%io_input,*) NumBcFlagsmm
  
  ALLOCATE(glb%bcCondmm(1:NumBcFlagsmm))
  
  DO i = 1, NumBcFlagsmm
     READ(glb%io_input,*)  iaux,glb%bcCondmm(i)%BCtypeX,glb%bcCondmm(i)%BCtypeY,glb%bcCondmm(i)%BCtypeZ, &
          glb%bcCondmm(i)%BCvalueX,glb%bcCondmm(i)%BCvalueY,glb%bcCondmm(i)%BCvalueZ
  ENDDO
  
  RETURN
END SUBROUTINE BOUNDARYMM_SUB

SUBROUTINE AMPLITUDE_SUB(glb,keywd)

  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  
  TYPE(ROCFRAC_GLOBAL) :: glb

  CHARACTER :: keywd*200
  INTEGER :: k1, k2
  
  INTEGER :: i              ! loop counter

  CHARACTER :: AmpType*16
  REAL*8 :: der1
  REAL*8 :: slp, intcpt, rise, run
  REAL*8, POINTER, DIMENSION(:,:) :: TableVal
  REAL*8 :: x1, y1, x2, y2

  CALL locchr(keywd,'DEFINITION                ',10,10,k1,k2)

  AmpType = keywd(k1:k2)

  IF(AmpType.EQ.'TABULAR')THEN

     glb%AmplitudeTable = .true.

     READ(glb%io_input,*) glb%NumEntries

! Time
! 1st Derivative
! 2nd Derivative

     ALLOCATE(TableVal(1:2,glb%NumEntries))
     DO i = 1, glb%NumEntries
        READ(glb%io_input,*) TableVal(1:2,i)
     ENDDO

     glb%NumEntries = glb%NumEntries - 1

     ALLOCATE( glb%AmpTable(1:3,glb%NumEntries))

     DO i = 1, glb%NumEntries
        
        x1 = TableVal(1,i)
        y1 = TableVal(2,i)
        x2 = TableVal(1,i+1)
        y2 = TableVal(2,i+1)

        rise = y2-y1
        run = x2-x1

        slp = (y2-y1)/(x2-x1)

        IF(ABS(slp*x1).LT.1.0e-6)THEN
           intcpt = y1
        ELSE
           intcpt = y1/(slp*x1)
        ENDIF

        glb%AmpTable(1,i) = TableVal(1,i)
        glb%AmpTable(2,i) = slp
        glb%AmpTable(3,i) = intcpt


     ENDDO

     DEALLOCATE(TableVal)

  ENDIF

  RETURN
END SUBROUTINE AMPLITUDE_SUB

SUBROUTINE ELEMENT_SUB(glb,keywd)

!!****f* Rocfrac/Source/feminp/ELEMENT_SUB
!!
!!  NAME
!!     ELEMENT_SUB
!!
!!  FUNCTION
!!     Specifies the element type
!!     
!!  INPUTS
!!     glb -- global array
!!
!!****

  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  
  TYPE(ROCFRAC_GLOBAL) :: glb

  CHARACTER(len=200) :: keywd
  INTEGER :: k1, k2

  CHARACTER(len=16) :: ElType

  CALL locchr(keywd,'TYPE                      ',4,8,k1,k2)

  ElType = keywd(k1:k2)

  SELECT CASE (TRIM(ElType))
  CASE ('V3D4')
     glb%iElType = 4
  CASE ('V3D4NCC')
     glb%iElType = 4
     glb%NdMassLump = 1
     glb%NdBasedEl = .TRUE.
  CASE ('V3D4N')
     glb%iElType = 4
     glb%NdBasedEl = .TRUE.
  CASE ('V3D10R')
     glb%iElType = 10
     glb%iElIntgratn = 1
  CASE ('V3D10')
     glb%iElType = 10
  CASE ('V3D10BBAR')
     glb%iElType = 10
     glb%iElIntgratn = 1
  CASE ('V3D8ME')
     glb%iElType = 8
  CASE default
     PRINT*,' ERROR:'
     PRINT*,'*ELEMENT TYPE NOT FOUND'
     STOP
  END SELECT

  RETURN
END SUBROUTINE ELEMENT_SUB
        
SUBROUTINE HEAT_TRANSFER_SUB(glb, keywd, tmp_KappaHT, tmp_Cp, tmp_iSolnTypeHT)

  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE

  TYPE(ROCFRAC_GLOBAL) :: glb
  
  CHARACTER :: keywd*200
  INTEGER :: key ! 0 = no key word found, 1 = key word found)
  CHARACTER*26 :: MatType
  INTEGER :: NumMatTypeHT
  INTEGER :: i,ii,ierr, j
  REAL*8, DIMENSION(1:10) :: tmp_KappaHT, tmp_Cp
  INTEGER, DIMENSION(1:10) :: tmp_iSolnTypeHT

  READ(glb%io_input,*) NumMatTypeHT

  DO i = 1, NumMatTypeHT
     glb%NumMatVolHT = glb%NumMatVolHT + 1
     
     ii = glb%NumMatVolHT
     
     IF(ii.GT.10)THEN
        PRINT*,'ROCFRAC :: ERROR'
        PRINT*,'Number of materials GREATER then 10'
        CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)
     ENDIF
     
     READ(glb%io_input,*) j, tmp_KappaHT(ii),tmp_Cp(ii)
   
     tmp_iSolnTypeHT(ii) = j   ! > 0 heat transfer solution wanted, 0 = no heat transfer solution

!!$! courant limit for temperature problem
!!$!
!!$!          2                          
!!$! t <=   Dx
!!$!      ------      (alpha = Kappa / (rho*cp) )
!!$!      3 alpha
!!$!
!!$     AlphaHT = 3.* (KappaHT/RhoCp)

     
  ENDDO

  glb%HeatTransSoln = .TRUE.    

END SUBROUTINE HEAT_TRANSFER_SUB
     


SUBROUTINE PROBE_SUB(glb)

!!****f* Rocfrac/Source/feminp/PROBE_SUB
!!
!!  NAME
!!     PROBE_SUB
!!
!!  FUNCTION
!!     This option is used to mark a node near a
!!     specified coordinate
!!     
!!
!!  INPUTS
!!     glb -- global array
!!
!!****

  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  
  TYPE(ROCFRAC_GLOBAL) :: glb
  
  INTEGER :: i             ! loop counter
  
  READ(glb%io_input,*)  glb%NumProbesNd

  ALLOCATE(glb%ProbeCoorNd(1:3,1:glb%NumProbesNd))
  ALLOCATE(glb%ProbeNd(1:glb%NumProbesNd))

  DO i = 1, glb%NumProbesNd
     READ(glb%io_input,*) glb%ProbeCoorNd(1:3,i)
  ENDDO
  
  RETURN
END SUBROUTINE PROBE_SUB


SUBROUTINE MICROMECHANICAL_SUB(glb, keywd, &
     tmp_E, tmp_xnu, tmp_rho, tmp_alpha, tmp_iSolnType)

!!****f* Rocfrac/Source/feminp/MICROMECHANICAL
!!
!!  NAME
!!     MICROMECHANICAL
!!
!!  FUNCTION
!!     micro mechanical material model
!!     
!!  INPUTS
!!     glb -- global array
!!     keywd -- keywd for control deck
!!  OUTPUT
!!     tmp_E -- Young's Modulus
!!     tmp_xnu -- Possion's ratio
!!     tmp_rho -- Density
!!     tmp_alpha -- thermal coefficient of expansion
!!     tmp_iSolnType -- material type model 
!!
!!****

  USE ROCSTAR_RocFrac

  IMPLICIT NONE

  TYPE(ROCFRAC_GLOBAL) :: glb

  CHARACTER :: keywd*200
  INTEGER :: key ! 0 = no key word found, 1 = key word found)
  INTEGER :: k1, k2, i
  CHARACTER :: ModelType*16
  REAL*8, DIMENSION(1:10) :: tmp_E, tmp_xnu, tmp_rho, tmp_alpha
  INTEGER, DIMENSION(1:10) :: tmp_iSolnType

  INTEGER :: NumMatVol_loc

  integer :: NLGeomType

  NLGeomType = 1            ! defaults to nl 

  CALL locchr(keywd,'MODEL                     ',5,8,k1,k2)

  IF(k1.GT.0.and.k2.GT.k1)THEN

     ModelType = keywd(k1:k2)

     IF(ModelType.EQ.'HUANG')THEN
       ! currently not used
        glb%NSTATEV = 1


        glb%NMATRIX = 3
        ALLOCATE(glb%MATRIX(1:glb%NMATRIX))

        READ(glb%io_input,*) glb%MATRIX(1:glb%NMATRIX), tmp_rho(1)

        READ(glb%io_input,*)  glb%NPARTICLETYPE

        glb%NPARTICLE = 4

        allocate(glb%PARTICLE(1:glb%NPARTICLE,1:glb%NPARTICLETYPE))
        DO i = 1, glb%NPARTICLETYPE
           READ(glb%io_input,*) glb%PARTICLE(1:4,i)
        enddo

!-----  two particle sizes: 
!     the code in the following is limited to two particle sizes, with the
!     same elastic moduli.

        IF (glb%particle(1,1)-glb%particle(1,2).NE.0.d0 .OR. &
             glb%particle(2,1)-glb%particle(2,2).NE.0.d0) THEN
           PRINT*,'ROCFRAC: Error: the particles do not have the same moduli.'
           STOP
        END IF

        IF (ABS( glb%particle(3,1)+glb%particle(3,2) +glb%matrix(3)-1.d0).GT.0.001) THEN
           PRINT*,'ROCFRAC: Error: the volume fractions of particles and matrix'
           PRINT*,'do not add up to 1'
           STOP
        END IF

  
        glb%NINTERFAC = 3
        allocate(glb%INTERFAC(1:glb%NINTERFAC))
        READ(glb%io_input,*) glb%INTERFAC(1:3)
 
        glb%NumMatVol = 1
        glb%DebondPart = .true.


        CALL HomogenizedMat( tmp_rho(1), glb%Matrix(1), glb%Matrix(2), &
             glb%PARTICLE(1,1), glb%PARTICLE(2,1), glb%MATRIX(3), glb%cd_fastest )
  
!!$        glb%cd_fastest = 0.d0
!!$        DO i = 1, 1
!!$           glb%cd_fastest = MAX( glb%cd_fastest, &
!!$                SQRT(glb%MATRIX(1)*(1.d0-glb%MATRIX(2))/tmp_rho(1)/(1.d0+glb%MATRIX(2))/(1.d0-2.d0*glb%MATRIX(2))) )
!!$        ENDDO
!!$        DO i = 1, glb%NPARTICLETYPE
!!$           glb%cd_fastest = MAX( glb%cd_fastest, &
!!$                SQRT(glb%PARTICLE(1,i)*(1.d0-glb%PARTICLE(2,i))/tmp_rho(1)/(1.d0+glb%PARTICLE(2,i))/(1.d0-2.d0*glb%PARTICLE(2,i))) )
!!$        ENDDO


     ELSE IF(ModelType.EQ.'MATOUS')THEN 

           glb%DebondPart_Matous = .TRUE.
           glb%NumMatVol = glb%NumMatVol + 1

           READ(glb%io_input,*) glb%NumMatVol_Part

           ALLOCATE(glb%E1(1:glb%NumMatVol_Part) )
           ALLOCATE(glb%E2(1:glb%NumMatVol_Part) )
           ALLOCATE(glb%E3(1:glb%NumMatVol_Part) )
           ALLOCATE(glb%nu12(1:glb%NumMatVol_Part) )
           ALLOCATE(glb%nu13(1:glb%NumMatVol_Part) )
           ALLOCATE(glb%nu23(1:glb%NumMatVol_Part) )
           ALLOCATE(glb%G12(1:glb%NumMatVol_Part) )
           ALLOCATE(glb%G13(1:glb%NumMatVol_Part) )
           ALLOCATE(glb%G23(1:glb%NumMatVol_Part) )

           DO i = 1, glb%NumMatVol_Part
              READ(glb%io_input,*) glb%E1(i), glb%E2(i), glb%E3(i),  &
                   glb%nu12(i), glb%nu13(i), glb%nu23(i),  &
                   glb%G12(i), glb%G13(i), glb%G23(i), tmp_rho(i)
     
           ENDDO

           READ(glb%io_input,*) glb%alpha1, glb%alpha2, glb%c2, glb%p1, glb%p2, glb%Yin, glb%a_eta, glb%a_zeta

           glb%cm = 1.d0 - glb%c2
           glb%cb = glb%c2
  
     ENDIF

  ELSE
     PRINT*,'ERROR in *MICROMECHANICAL keywd'
     PRINT*,'MODEL type not found'

  ENDIF


  tmp_iSolnType(1) = NLGeomType ! fix should not be 1


END SUBROUTINE MICROMECHANICAL_SUB

SUBROUTINE INITIALCONDITION_SUB(glb,keywd)

!!****f* Rocfrac/Source/feminp/INITIALCONDITION_SUB
!!
!!  NAME
!!     INITIALCONDITION_SUB
!!
!!  FUNCTION
!!     Specifies reference conditions
!!     
!!  INPUTS
!!     glb -- global array
!!
!!****

  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  
  TYPE(ROCFRAC_GLOBAL) :: glb

  CHARACTER(len=200) :: keywd
  INTEGER :: k1, k2

  CALL locchr(keywd,'TYPE                      ',4,16,k1,k2)

  
  IF(keywd(k1:k2).EQ.'TEMPERATURE')THEN
     READ(glb%io_input,*) glb%Temperature0
  ENDIF

  RETURN
END SUBROUTINE INITIALCONDITION_SUB

