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
!******************************************************************************
!
! Purpose: register grid and flow variables with GenX.
!
! Description: none.
!
! Input: regions        = dimensions of boundary patches, types of BC`s
!        handle, solver = GenX stuff.
!
! Output: to Roccom.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_InitGenxInterface.F90,v 1.31 2010/02/18 21:47:38 juzhang Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_InitGenxInterface( regions,handle,solver,inSurf,inVolPlag, &
                                   obtain_dataitem )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes
#ifdef PEUL
  USE ModInterfaces, ONLY : PEUL_InitGenxInterface
#endif
#ifdef PLAG
  USE ModInterfaces, ONLY : PLAG_InitGenxInterface
#endif
#ifdef RADI
  USE ModInterfaces, ONLY : RADI_InitGenxInterface
#endif
#ifdef TURB
  USE ModInterfaces, ONLY : TURB_InitGenxInterface
#endif
  USE ModInterfaces, ONLY : randInitGenxInterface
  USE ModError
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE
  INCLUDE 'comf90.h'

! ... parameters
  CHARACTER(*) :: inSurf, inVolPlag
  INTEGER, INTENT(IN) :: handle, solver, obtain_dataitem

  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  CHARACTER(CHRLEN) :: wins, winv, inVol
#ifdef PLAG
  CHARACTER(CHRLEN) :: inPlag
#endif
#ifdef STATS
  CHARACTER(CHRLEN), POINTER :: statNm(:,:,:)
  INTEGER :: iStat
#endif
  INTEGER :: iLev, bcType, pid, icount, errorFlag, ilb, mpierr
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend, iProc
  INTEGER, POINTER :: dims(:), iConstrType
  LOGICAL :: fileExist

  TYPE(t_global), POINTER :: global
  TYPE(t_patch) , POINTER :: patch

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_InitGenxInterface',&
  'RFLO_InitRocstarInterface.F90' )

! open data windows and register variables ------------------------------------

  wins = TRIM(global%winName)//'_surf'
  winv = TRIM(global%winName)//'_vol'

#ifdef PLAG
! obtain names of inVol and inPlag
  READ (inVolPlag, *) inVol, inPlag
#else
  READ (inVolPlag, *) inVol
#endif

  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  IF ( regions(1)%global%myProcid == MASTERPROC .AND. &
       regions(1)%global%verbLevel>= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating Rocstar surface windows'
 ENDIF

  CALL COM_new_window( TRIM(wins) )

  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  IF ( regions(1)%global%myProcid == MASTERPROC .AND. &
       regions(1)%global%verbLevel>= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating Rocstar surface input dataitems'
 ENDIF
! input data

  CALL COM_new_dataitem( TRIM(wins)//'.du_alp'    ,'n',COM_DOUBLE,3, &
                          'm' )
  CALL COM_new_dataitem( TRIM(wins)//'.mdot_alp'  ,'e',COM_DOUBLE,1, &
                          'kg/(m^2s)' )
  CALL COM_new_dataitem( TRIM(wins)//'.rhofvf_alp','e',COM_DOUBLE,3, &
                          'kg/(m^2s)' )
  CALL COM_new_dataitem( TRIM(wins)//'.Tflm_alp'  ,'e',COM_DOUBLE,1, &
                          'K' )

  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  IF ( regions(1)%global%myProcid == MASTERPROC .AND. &
       regions(1)%global%verbLevel>= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating Rocstar surface output dataitems'
 ENDIF
! output data

  CALL COM_new_dataitem( TRIM(wins)//'.pf'      ,'e',COM_DOUBLE,1, &
                          'Pa' )
  CALL COM_new_dataitem( TRIM(wins)//'.qc'      ,'e',COM_DOUBLE,1, &
                          'kgK/(m^2s)' )
  CALL COM_new_dataitem( TRIM(wins)//'.qr'      ,'e',COM_DOUBLE,1, &
                          'kgK/(m^2s)' )
  CALL COM_new_dataitem( TRIM(wins)//'.rhof_alp','e',COM_DOUBLE,1, &
                          'kg/m^3' )
  CALL COM_new_dataitem( TRIM(wins)//'.nf_alp'  ,'e',COM_DOUBLE,3, &
                          '' )
  CALL COM_new_dataitem( TRIM(wins)//'.tf'      ,'e',COM_DOUBLE,3, &
                          'Pa' )
  CALL COM_new_dataitem( TRIM(wins)//'.Tf'      ,'e',COM_DOUBLE,1, &
                          'K' )
  CALL COM_new_dataitem( TRIM(wins)//'.Tv'      ,'e',COM_DOUBLE,1, &
                          'K' )
  CALL COM_new_dataitem( TRIM(wins)//'.dn'      ,'e',COM_DOUBLE,1, &
                          'm' )
  CALL COM_new_dataitem( TRIM(wins)//'.bflag'   ,'e',COM_INTEGER,1,'' )

  CALL COM_new_dataitem( TRIM(wins)//'.bcflag'  ,'p',COM_INTEGER,1,'' )

  CALL COM_new_dataitem( TRIM(wins)//'.cnstr_type','p',COM_INTEGER,1,'' )

! restart data (si/j/kvel, cv) and additional plot data (dv=p, T, c)

  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  IF ( regions(1)%global%myProcid == MASTERPROC .AND. &
       regions(1)%global%verbLevel>= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Creating Rocstar volume windows'
 ENDIF

  CALL COM_new_window( TRIM(winv) )

  CALL COM_new_dataitem( TRIM(winv)//'.siVel','n' ,COM_DOUBLE,1,'m/s' )
  CALL COM_new_dataitem( TRIM(winv)//'.sjVel','n' ,COM_DOUBLE,1,'m/s' )
  CALL COM_new_dataitem( TRIM(winv)//'.skVel','n' ,COM_DOUBLE,1,'m/s' )

  CALL COM_new_dataitem( TRIM(winv)//'.dtf','e',COM_DOUBLE,1,'s' )

  CALL COM_new_dataitem( TRIM(winv)//'.rhof' ,'e',COM_DOUBLE,1,&
                          'kg/(m^3)')
  CALL COM_new_dataitem( TRIM(winv)//'.rhovf','e',COM_DOUBLE,3,&
                          'kg/(m^2 s)')
  CALL COM_new_dataitem( TRIM(winv)//'.rhoEf','e',COM_DOUBLE,1,&
                          '(J/kg)')

  CALL COM_new_dataitem( TRIM(winv)//'.vf','e',COM_DOUBLE,3,'m/s' )
  CALL COM_new_dataitem( TRIM(winv)//'.Tf','e',COM_DOUBLE,1,'K' )
  CALL COM_new_dataitem( TRIM(winv)//'.pf','e',COM_DOUBLE,1,'Pa' )
  CALL COM_new_dataitem( TRIM(winv)//'.Tv','e',COM_DOUBLE,1,'K' )
  CALL COM_new_dataitem( TRIM(winv)//'.dn','e',COM_DOUBLE,1,'m' )
  CALL COM_new_dataitem( TRIM(winv)//'.af','e',COM_DOUBLE,1,'m/s' )

  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  IF ( regions(1)%global%myProcid == MASTERPROC .AND. &
       regions(1)%global%verbLevel>= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Main windows created.'
 ENDIF

! statistics

#ifdef STATS
  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  IF ( regions(1)%global%myProcid == MASTERPROC .AND. &
       regions(1)%global%verbLevel>= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Adding statistics dataitems'
 ENDIF

  IF ((global%flowType == FLOW_UNSTEADY) .AND. (global%doStat == ACTIVE)) THEN

! - stats global data
    CALL COM_new_dataitem( TRIM(winv)//'.tStat','w',COM_DOUBLE,1,'s' )

    IF (global%mixtNStat > 0) THEN
      statNm => global%mixtStatNm
      DO iStat=1,global%mixtNStat
        CALL COM_new_dataitem( TRIM(winv)//'.'//TRIM(statNm(1,1,iStat)),'e', &
                                COM_DOUBLE,1,TRIM(statNm(1,2,iStat)) )
      ENDDO
    ENDIF
    CALL COM_set_array( TRIM(winv)//'.tStat',0, global%integrTime )

  ENDIF  ! unsteady and dostat
#endif

! store pointers to variables -------------------------------------------------

  ALLOCATE( dims(3),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! loop over all regions

  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  IF ( regions(1)%global%myProcid == MASTERPROC .AND. &
       regions(1)%global%verbLevel>= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Processing Rocstar surface windows'
 ENDIF

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

      iLev   = regions(iReg)%currLevel
      icount = 0

! --- surface data

      DO iPatch=1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(iLev)%patches(iPatch)
        bcType =  patch%bcType
!        WRITE(*,*) 'doing patch with bc ',bcType
!        IF(bcType .NE. BC_SYMMETRY) THEN
        IF ( patch%bcCoupled == BC_EXTERNAL .OR. &   ! data from outside
            (patch%bcCoupled == BC_INTERNAL .AND. &  ! data from internal APN
             bcType == BC_INJECTION_APN)) THEN        
          icount  = icount + 1
          pid     = iReg*REGOFF + icount
!          WRITE(*,*) ' external bc on patch '
! ------- burning pane?

          ALLOCATE( patch%bcFlag(1),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

          CALL COM_set_size(TRIM(wins)//'.bcflag',pid,1)
          CALL COM_set_array( TRIM(wins)//'.bcflag',pid,patch%bcFlag )

          IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) THEN
            patch%bcFlag(1) = 1   ! ignitable patch
          ELSE
            patch%bcFlag(1) = 0   ! non-ignitable patch
          ENDIF

! ------- surface grid

          dims(1) = ABS(patch%l1end-patch%l1beg) + 2    ! nodal values
          dims(2) = ABS(patch%l2end-patch%l2beg) + 2

          ALLOCATE( patch%surfCoord(3,dims(1),dims(2)),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

          CALL COM_set_array_const( TRIM(wins)//'.:st2:',pid,dims )
          CALL COM_set_array( TRIM(wins)//'.nc',pid,patch%surfCoord )

! ------- input data

          dims(1) = ABS(patch%l1end-patch%l1beg) + 2    ! nodal values
          dims(2) = ABS(patch%l2end-patch%l2beg) + 2
          ALLOCATE( patch%duAlp(3,dims(1),dims(2)),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

          CALL COM_set_array( TRIM(wins)//'.du_alp',pid,patch%duAlp )

          dims(1) = ABS(patch%l1end-patch%l1beg) + 1    ! cell values
          dims(2) = ABS(patch%l2end-patch%l2beg) + 1
          ALLOCATE( patch%mdotAlp  (  dims(1),dims(2)),stat=errorFlag )
          ALLOCATE( patch%rhofvfAlp(3,dims(1),dims(2)),stat=errorFlag )
          ALLOCATE( patch%tflmAlp  (  dims(1),dims(2)),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

          CALL COM_set_array( TRIM(wins)//'.mdot_alp'  ,pid,patch%mdotAlp )
          CALL COM_set_array( TRIM(wins)//'.rhofvf_alp',pid,patch%rhofvfAlp )
          CALL COM_set_array( TRIM(wins)//'.Tflm_alp'  ,pid,patch%tflmAlp )

! ------- output data

          dims(1) = ABS(patch%l1end-patch%l1beg) + 1    ! cell values
          dims(2) = ABS(patch%l2end-patch%l2beg) + 1
          ALLOCATE( patch%pf     (  dims(1),dims(2)),stat=errorFlag )
          ALLOCATE( patch%qc     (  dims(1),dims(2)),stat=errorFlag )
          ALLOCATE( patch%qr     (  dims(1),dims(2)),stat=errorFlag )
          ALLOCATE( patch%rhofAlp(  dims(1),dims(2)),stat=errorFlag )
          ALLOCATE( patch%tempf  (  dims(1),dims(2)),stat=errorFlag )
          ALLOCATE( patch%tempv  (  dims(1),dims(2)),stat=errorFlag )
          ALLOCATE( patch%dnml   (  dims(1),dims(2)),stat=errorFlag )
          ALLOCATE( patch%nfAlp  (3,dims(1),dims(2)),stat=errorFlag )
          ALLOCATE( patch%tracf  (3,dims(1),dims(2)),stat=errorFlag )
          ALLOCATE( patch%bFlag  (  dims(1),dims(2)),stat=errorFlag )
          global%error = errorFlag
          IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

          CALL COM_set_array( TRIM(wins)//'.pf'      ,pid,patch%pf      )
          CALL COM_set_array( TRIM(wins)//'.qc'      ,pid,patch%qc      )
          CALL COM_set_array( TRIM(wins)//'.qr'      ,pid,patch%qr      )
          CALL COM_set_array( TRIM(wins)//'.rhof_alp',pid,patch%rhofAlp )
          CALL COM_set_array( TRIM(wins)//'.nf_alp'  ,pid,patch%nfAlp   )
          CALL COM_set_array( TRIM(wins)//'.tf'      ,pid,patch%tracf   )
          CALL COM_set_array( TRIM(wins)//'.Tf'      ,pid,patch%tempf   )
          CALL COM_set_array( TRIM(wins)//'.Tv'      ,pid,patch%tempv   )
          CALL COM_set_array( TRIM(wins)//'.dn'      ,pid,patch%dnml   )
          CALL COM_set_array( TRIM(wins)//'.bflag'   ,pid,patch%bFlag   )

! ------- constraint type

          CALL COM_set_size(TRIM(wins)//'.cnstr_type',pid,1)
          CALL COM_allocate_array( TRIM(wins)//'.cnstr_type',pid,iConstrType )

          iConstrType = 0
          IF (bcType==BC_SYMMETRY_FREE) THEN
            iConstrType = 0
          ELSEIF (bcType==BC_SYMMETRY_FIXED) THEN
            iConstrType = 2
          ELSEIF (bcType==BC_SYMMETRY_XSLIDE) THEN
            iConstrType = 120
          ELSEIF (bcType==BC_SYMMETRY_YSLIDE) THEN
            iConstrType = 121
          ELSEIF (bcType==BC_SYMMETRY_ZSLIDE) THEN
            iConstrType = 122
          ELSEIF (bcType==BC_SYMMETRY_XYSLIDE) THEN
            iConstrType = -122
          ELSEIF (bcType==BC_SYMMETRY_XZSLIDE) THEN
            iConstrType = -121
          ELSEIF (bcType==BC_SYMMETRY_YZSLIDE) THEN
            iConstrType = -120
          ENDIF
          IF (bcType==BC_SLIPWALL_FREE) THEN
            iConstrType = 0
          ELSEIF (bcType==BC_SLIPWALL_FIXED) THEN
            iConstrType = 2
          ELSEIF (bcType==BC_SLIPWALL_XSLIDE) THEN
            iConstrType = 120
          ELSEIF (bcType==BC_SLIPWALL_YSLIDE) THEN
            iConstrType = 121
          ELSEIF (bcType==BC_SLIPWALL_ZSLIDE) THEN
            iConstrType = 122
          ELSEIF (bcType==BC_SLIPWALL_XYSLIDE) THEN
            iConstrType = -122
          ELSEIF (bcType==BC_SLIPWALL_XZSLIDE) THEN
            iConstrType = -121
          ELSEIF (bcType==BC_SLIPWALL_YZSLIDE) THEN
            iConstrType = -120
          ENDIF
          IF (bcType==BC_NOSLIPWALL_FREE) THEN
            iConstrType = 0
          ELSEIF (bcType==BC_NOSLIPWALL_FIXED) THEN
            iConstrType = 2
          ELSEIF (bcType==BC_NOSLIPWALL_XSLIDE) THEN
            iConstrType = 120
          ELSEIF (bcType==BC_NOSLIPWALL_YSLIDE) THEN
            iConstrType = 121
          ELSEIF (bcType==BC_NOSLIPWALL_ZSLIDE) THEN
            iConstrType = 122
          ELSEIF (bcType==BC_NOSLIPWALL_XYSLIDE) THEN
            iConstrType = -122
          ELSEIF (bcType==BC_NOSLIPWALL_XZSLIDE) THEN
            iConstrType = -121
          ELSEIF (bcType==BC_NOSLIPWALL_YZSLIDE) THEN
            iConstrType = -120
          ENDIF
          IF (bcType==BC_OUTFLOW_FREE) THEN
            iConstrType = 0
          ELSEIF (bcType==BC_OUTFLOW_FIXED) THEN
            iConstrType = 2
          ELSEIF (bcType==BC_OUTFLOW_XSLIDE) THEN
            iConstrType = 120
          ELSEIF (bcType==BC_OUTFLOW_YSLIDE) THEN
            iConstrType = 121
          ELSEIF (bcType==BC_OUTFLOW_ZSLIDE) THEN
            iConstrType = 122
          ELSEIF (bcType==BC_OUTFLOW_XYSLIDE) THEN
            iConstrType = -122
          ELSEIF (bcType==BC_OUTFLOW_XZSLIDE) THEN
            iConstrType = -121
          ELSEIF (bcType==BC_OUTFLOW_YZSLIDE) THEN
            iConstrType = -120
          ENDIF
!          PRINT *,'RFLO: bcType = ',bcType,' iConstrType = ',iConstrType

! ------- zero out radiation flux (set by Rocrad if active)

          patch%qr(:,:) = 0._RFREAL

        ELSE     ! internal BC
!          WRITE(*,*) ' internal bc on patch '

#ifndef PRE_RFLOPREP_V2300
          IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
              (bcType>=BC_REGIONINT .AND. bcType<=BC_REGIONINT+BC_RANGE) .OR. &
              (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
              (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE) .OR. &
              (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE)) THEN

          ELSE
            icount  = icount + 1
            pid     = iReg*REGOFF + icount

            ALLOCATE( patch%bcFlag(1),stat=errorFlag )
            global%error = errorFlag
            IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

            CALL COM_set_size(TRIM(wins)//'.bcflag',pid,1)
            CALL COM_set_array( TRIM(wins)//'.bcflag',pid,patch%bcFlag )

            patch%bcFlag(1) = 2   ! non-interacting patch

! --------- surface grid

            dims(1) = ABS(patch%l1end-patch%l1beg) + 2    ! nodal values
            dims(2) = ABS(patch%l2end-patch%l2beg) + 2

            ALLOCATE( patch%surfCoord(3,dims(1),dims(2)),stat=errorFlag )
            global%error = errorFlag
            IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

            CALL COM_set_array_const( TRIM(wins)//'.:st2:',pid,dims )
            CALL COM_set_array( TRIM(wins)//'.nc',pid,patch%surfCoord )

            ALLOCATE( patch%duAlp(3,dims(1),dims(2)),stat=errorFlag )
            global%error = errorFlag
            IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

            CALL COM_set_array( TRIM(wins)//'.du_alp',pid,patch%duAlp )

! --------- constraint type

            CALL COM_set_size(TRIM(wins)//'.cnstr_type',pid,1)
            CALL COM_allocate_array( TRIM(wins)//'.cnstr_type',pid,iConstrType )

            iConstrType = 2
            IF (bcType==BC_SYMMETRY_FREE) THEN
               iConstrType = 0
            ELSEIF (bcType==BC_SYMMETRY_FIXED) THEN
               iConstrType = 2
            ELSEIF (bcType==BC_SYMMETRY_XSLIDE) THEN
               iConstrType = 120
            ELSEIF (bcType==BC_SYMMETRY_YSLIDE) THEN
               iConstrType = 121
            ELSEIF (bcType==BC_SYMMETRY_ZSLIDE) THEN
               iConstrType = 122
            ELSEIF (bcType==BC_SYMMETRY_XYSLIDE) THEN
               iConstrType = -122
            ELSEIF (bcType==BC_SYMMETRY_XZSLIDE) THEN
               iConstrType = -121
            ELSEIF (bcType==BC_SYMMETRY_YZSLIDE) THEN
               iConstrType = -120
            ENDIF
            IF (bcType==BC_SLIPWALL_FREE) THEN
              iConstrType = 0
            ELSEIF (bcType==BC_SLIPWALL_FIXED) THEN
              iConstrType = 2
            ELSEIF (bcType==BC_SLIPWALL_XSLIDE) THEN
              iConstrType = 120
            ELSEIF (bcType==BC_SLIPWALL_YSLIDE) THEN
              iConstrType = 121
            ELSEIF (bcType==BC_SLIPWALL_ZSLIDE) THEN
              iConstrType = 122
            ELSEIF (bcType==BC_SLIPWALL_XYSLIDE) THEN
              iConstrType = -122
            ELSEIF (bcType==BC_SLIPWALL_XZSLIDE) THEN
              iConstrType = -121
            ELSEIF (bcType==BC_SLIPWALL_YZSLIDE) THEN
!               write(*,*) 'RFLOn: BC_*YZSLIDE'
              iConstrType = -120
            ENDIF
            IF (bcType==BC_NOSLIPWALL_FREE) THEN
              iConstrType = 0
            ELSEIF (bcType==BC_NOSLIPWALL_FIXED) THEN
              iConstrType = 2
            ELSEIF (bcType==BC_NOSLIPWALL_XSLIDE) THEN
              iConstrType = 120
            ELSEIF (bcType==BC_NOSLIPWALL_YSLIDE) THEN
              iConstrType = 121
            ELSEIF (bcType==BC_NOSLIPWALL_ZSLIDE) THEN
              iConstrType = 122
            ELSEIF (bcType==BC_NOSLIPWALL_XYSLIDE) THEN
              iConstrType = -122
            ELSEIF (bcType==BC_NOSLIPWALL_XZSLIDE) THEN
              iConstrType = -121
            ELSEIF (bcType==BC_NOSLIPWALL_YZSLIDE) THEN
 !              write(*,*) 'RFLOn: BC_*YZSLIDE'
              iConstrType = -120
            ENDIF
            IF (bcType==BC_OUTFLOW_FREE) THEN
              iConstrType = 0
            ELSEIF (bcType==BC_OUTFLOW_FIXED) THEN
              iConstrType = 2
            ELSEIF (bcType==BC_OUTFLOW_XSLIDE) THEN
              iConstrType = 120
            ELSEIF (bcType==BC_OUTFLOW_YSLIDE) THEN
              iConstrType = 121
            ELSEIF (bcType==BC_OUTFLOW_ZSLIDE) THEN
              iConstrType = 122
            ELSEIF (bcType==BC_OUTFLOW_XYSLIDE) THEN
              iConstrType = -122
            ELSEIF (bcType==BC_OUTFLOW_XZSLIDE) THEN
              iConstrType = -121
            ELSEIF (bcType==BC_OUTFLOW_YZSLIDE) THEN
!               write(*,*) 'RFLOn: BC_*YZSLIDE'
              iConstrType = -120
            ENDIF
!            PRINT *,'RFLOn: bcType = ',bcType ,' iConstrType = ',iConstrType

          ENDIF
#endif
        ENDIF    ! external/internal BC
!     ENDIF
!        WRITE(*,*) ' done with patch having bc',bcType
     ENDDO      ! iPatch

! --- volume data
!      WRITE(*,*) 'now doing volume data'
  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  IF ( regions(1)%global%myProcid == MASTERPROC .AND. &
       regions(1)%global%verbLevel>= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Processing Rocstar volume windows.'
 ENDIF

      pid = iReg*REGOFF

      CALL RFLO_GetDimensDummyNodes( regions(iReg),iLev,idnbeg,idnend, &
                                     jdnbeg,jdnend,kdnbeg,kdnend )
      dims(1) = idnend - idnbeg + 1
      dims(2) = jdnend - jdnbeg + 1
      dims(3) = kdnend - kdnbeg + 1

      CALL COM_set_size( TRIM(winv)//".:st3:",pid,3, regions(iReg)%nDumCells)
      CALL COM_set_array_const( TRIM(winv)//".:st3:",pid,dims )
      CALL COM_set_array( TRIM(winv)//'.nc',pid, &
              regions(iReg)%levels(iLev)%grid%xyz )

!      WRITE(*,*) 'setting arrays'
      IF (regions(iReg)%mixtInput%moveGrid) THEN
        CALL COM_set_array( TRIM(winv)//'.siVel',pid, &
             regions(iReg)%levels(iLev)%grid%siVel )
        CALL COM_set_array( TRIM(winv)//'.sjVel',pid, &
             regions(iReg)%levels(iLev)%grid%sjVel )
        CALL COM_set_array( TRIM(winv)//'.skVel',pid, &
             regions(iReg)%levels(iLev)%grid%skVel )
      ENDIF

      ilb = LBOUND(regions(iReg)%levels(iLev)%mixt%cv,2)

      CALL COM_set_array( TRIM(winv)//'.rhof',pid, &
           regions(iReg)%levels(iLev)%mixt%cv(1,ilb),5)
      CALL COM_set_array( TRIM(winv)//'.1-rhovf',pid, &
           regions(iReg)%levels(iLev)%mixt%cv(2,ilb),5)
      CALL COM_set_array( TRIM(winv)//'.2-rhovf',pid, &
           regions(iReg)%levels(iLev)%mixt%cv(3,ilb),5)
      CALL COM_set_array( TRIM(winv)//'.3-rhovf',pid, &
           regions(iReg)%levels(iLev)%mixt%cv(4,ilb),5)
      CALL COM_set_array( TRIM(winv)//'.rhoEf',pid, &
           regions(iReg)%levels(iLev)%mixt%cv(5,ilb),5)

      CALL COM_set_array( TRIM(winv)//'.1-vf',pid, &
           regions(iReg)%levels(iLev)%mixt%dv(1,ilb),6)
      CALL COM_set_array( TRIM(winv)//'.2-vf',pid, &
           regions(iReg)%levels(iLev)%mixt%dv(2,ilb),6)
      CALL COM_set_array( TRIM(winv)//'.3-vf',pid, &
           regions(iReg)%levels(iLev)%mixt%dv(3,ilb),6)
      CALL COM_set_array( TRIM(winv)//'.Tf',pid, &
           regions(iReg)%levels(iLev)%mixt%dv(4,ilb),6)
      CALL COM_set_array( TRIM(winv)//'.pf',pid, &
           regions(iReg)%levels(iLev)%mixt%dv(5,ilb),6)
      CALL COM_set_array( TRIM(winv)//'.af',pid, &
           regions(iReg)%levels(iLev)%mixt%dv(6,ilb),6)

      CALL COM_set_array( TRIM(winv)//'.dtf',pid, &
           regions(iReg)%levels(iLev)%dt(ilb),1)

! --- statistics

#ifdef STATS
      IF ((global%flowType==FLOW_UNSTEADY) .AND. (global%doStat==ACTIVE)) THEN
        IF (global%mixtNStat > 0) THEN
          DO iStat=1,global%mixtNStat
            CALL COM_set_array( TRIM(winv)//'.'//TRIM(statNm(1,1,iStat)), pid,&
               regions(iReg)%levels(iLev)%mixt%tav(iStat,ilb), global%mixtNStat)
          ENDDO
        ENDIF
      ENDIF
#endif
    ENDIF    ! region on this processor and active
  ENDDO      ! iReg

  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  IF ( regions(1)%global%myProcid == MASTERPROC .AND. &
       regions(1)%global%verbLevel>= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initializing Rocstar interfaces for MP'
 ENDIF

! Genx initialization of physical modules -------------------------------------

  CALL randInitGenxInterface( regions,wins,winv )

#ifdef PEUL
  IF (global%peulUsed)   CALL PEUL_initGenxInterface( regions,wins,winv )
#endif
#ifdef PLAG
!#ifndef NATIVE_MP_IO
  IF (global%plagUsed)   CALL PLAG_initGenxInterface( regions,wins, &
                                                      inPlag,obtain_dataitem )
!#endif
#endif
#ifdef TURB
  IF (global%turbActive) CALL TURB_initGenxInterface( regions,wins,winv )
#endif
#ifdef RADI
  IF (global%radiActive) CALL RADI_initGenxInterface( regions,wins,winv )
#endif

! finalize --------------------------------------------------------------------

  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel>= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Window setup done.'
 ENDIF

  CALL COM_window_init_done( TRIM(wins) )
  CALL COM_window_init_done( TRIM(winv) )
!  WRITE(*,*) 'getting to dataitem grabbing'
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel>= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Populating volume windows.'
 ENDIF

 DO iProc = 0, global%nprocalloc
!    IF(global%myProcid == MASTERPROC) THEN
!       WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting volume data.',iProc
!    ENDIF
    IF(global%myProcid == iProc) THEN
  CALL COM_call_function( obtain_dataitem,2, &
                          COM_get_dataitem_handle_const(TRIM(inVol)//".all"), &
                          COM_get_dataitem_handle(TRIM(winv)//".all") )


    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
 ENDDo

 
  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  IF ( regions(1)%global%myProcid == MASTERPROC .AND. &
       regions(1)%global%verbLevel>= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Populating surface windows.'
 ENDIF

 DO iProc = 0, global%nprocalloc
!    IF(global%myProcid == MASTERPROC) THEN
!       WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Getting surface data.',iProc
!    ENDIF
    IF(global%myProcid == iProc) THEN
       CALL COM_call_function( obtain_dataitem,2, &
            COM_get_dataitem_handle_const(TRIM(inSurf)//".all"), &
            COM_get_dataitem_handle(TRIM(wins)//".all") )
       
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
 ENDDo

  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel>= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Done Populating all windows.'
 ENDIF
!  WRITE(*,*) 'past dataitem grabbing'
#ifdef PLAG
  CALL COM_call_function( handle,3,TRIM(wins),&
                          TRIM(winv)//' '//TRIM(global%winp),solver )
#else
  CALL COM_call_function( handle,3,TRIM(wins),TRIM(winv),solver )
#endif

! set tav from actual time averaged to accumulated values --------------------

#ifdef STATS
  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
      iLev = regions(iReg)%currLevel

      IF ((global%flowType==FLOW_UNSTEADY) .AND. (global%doStat==ACTIVE)) THEN
        IF (global%mixtNStat > 0) THEN
          DO iStat=1,global%mixtNStat
            regions(iReg)%levels(iLev)%mixt%tav(iStat,:) = &
            regions(iReg)%levels(iLev)%mixt%tav(iStat,:)*global%integrTime
          ENDDO
        ENDIF ! mixtNstat
#ifdef TURB
        IF ((global%turbActive .EQV. .true.) .AND. (global%turbNStat > 0)) THEN
          DO iStat=1,global%turbNStat
            regions(iReg)%levels(iLev)%turb%tav(iStat,:) = &
            regions(iReg)%levels(iLev)%turb%tav(iStat,:)*global%integrTime
          ENDDO
        ENDIF ! turbNstat
#endif
      ENDIF  ! unsteady and dostat
    ENDIF    ! region on this processor and active
  ENDDO      ! iReg
#endif

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_InitGenxInterface

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_InitGenxInterface.F90,v $
! Revision 1.31  2010/02/18 21:47:38  juzhang
! Heat transfer bc for non-propellant surface documented in Rocburn_PY_HT.pdf in Rocburn_PY directory is implemented within Rocburn_PY. Major changes were made to Rocburn, Rocman3, RocfluidMP/genx, RocfluidMP/modflo directories.
!
! Revision 1.30  2009/08/27 14:04:49  mtcampbe
! Updated to enable burning motion with symmetry boundaries and enhanced
! burnout code.
!
! Revision 1.29  2009/08/12 04:15:57  mtcampbe
! Major update, bugfix from Abe development, more propagation compatibility,
! some Rocstar IO changes, Ju's temporary clipping fix for turbulence. A bug
! fix for initialization IO.
!
! Revision 1.27  2008/12/06 08:44:00  mtcampbe
! Updated license.
!
! Revision 1.26  2008/11/19 22:17:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.25  2006/08/28 11:42:11  rfiedler
! Add grid motion constraint types for outflow BC.
!
! Revision 1.24  2006/08/24 14:58:46  rfiedler
! Use numbers for constraints instead of ICHAR; print constraints.
!
! Revision 1.22  2006/05/08 22:30:21  wasistho
! added prop-NS capability
!
! Revision 1.21  2006/03/03 07:06:06  wasistho
! backed out from zeroing out incoming alpha variables
!
! Revision 1.20  2006/03/03 06:07:31  wasistho
! initialized incoming alpha variables
!
! Revision 1.19  2006/01/23 22:44:03  wasistho
! added condition for internal injectionAPN
!
! Revision 1.18  2005/12/08 07:03:28  wasistho
! added ifdef TURB
!
! Revision 1.17  2005/12/08 02:34:21  wasistho
! bug fixed move registration of global%integrTime outside regions loop
!
! Revision 1.16  2005/12/08 00:18:46  wasistho
! stored actual time averaged vars in hdf
!
! Revision 1.15  2005/12/07 20:03:34  wasistho
! removed attemp to store actual tav i.o. accumulated
!
! Revision 1.14  2005/12/07 04:43:57  wasistho
! modified statistics treatment
!
! Revision 1.13  2005/12/07 02:23:36  wasistho
! added integrTime with eps
!
! Revision 1.12  2005/12/06 21:52:58  wasistho
! devided and multiply tav with integrTime
!
! Revision 1.11  2005/12/04 09:15:01  wasistho
! added statistics integration time
!
! Revision 1.10  2005/06/19 05:33:21  wasistho
! shift index rocprop slipwalls and change default iConstrType
!
! Revision 1.9  2005/06/17 03:09:26  wasistho
! relocated cnstr_type kernel inside both external internal loops
!
! Revision 1.8  2005/06/17 02:23:46  jiao
! Fixed bug in registering cnstr_type.
!
! Revision 1.7  2005/06/16 22:33:31  wasistho
! added cnstr_type
!
! Revision 1.6  2005/05/11 19:44:46  wasistho
! changed REG_NONINTERACT to PRE_RFLOPREP_V2300
!
! Revision 1.5  2005/05/10 15:01:40  wasistho
! exclude block interfaces in registration of internal surfaces
!
! Revision 1.4  2005/04/20 02:50:25  wasistho
! added error msg for incorrect compile option
!
! Revision 1.3  2005/04/18 20:34:55  wasistho
! added ifdef REG_NONINTERACT
!
! Revision 1.2  2005/04/18 18:11:44  wasistho
! registered non-interacting patches
!
! Revision 1.1  2004/12/01 21:23:52  haselbac
! Initial revision after changing case
!
! Revision 1.30  2004/07/02 22:48:40  jiao
! Fixed function call to Rocman when PLAG is defined.
!
! Revision 1.29  2004/07/02 22:06:41  fnajjar
! Modified PLAG call for Roccom3 import
!
! Revision 1.28  2004/06/30 04:05:56  wasistho
! moved Genx related parameter REGOFF to ModParameters
!
! Revision 1.27  2004/06/29 23:52:10  wasistho
! migrated to Roccom-3
!
! Revision 1.26  2004/06/07 23:05:21  wasistho
! provide Genx statistics names, units, and anytime-activation
!
! Revision 1.25  2004/03/05 22:08:58  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.24  2003/12/05 02:05:53  rfiedler
! Make the data in the fluids time step begin at index ilb like
! the other variables.  RAF
!
! Revision 1.23  2003/12/03 03:02:37  jiao
! Removed all calls involving COM_NULL.
!
! Revision 1.22  2003/12/02 21:20:37  fnajjar
! Included timestep size in output file
!
! Revision 1.21  2003/11/21 22:17:36  fnajjar
! Added PLAG, PEUL, and rand interfaces.
!
! Revision 1.20  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.17  2003/08/14 20:06:58  jblazek
! Corrected bug associated with radiation flux qr.
!
! Revision 1.16  2003/08/09 02:07:47  wasistho
! added TURB and RADI_initGenxInterface
!
! Revision 1.15  2003/05/15 02:57:00  jblazek
! Inlined index function.
!
! Revision 1.14  2003/05/09 17:01:03  jiao
! Renamed the COM_call_function_handlers to COM_call_function.
!
! Revision 1.13  2003/04/25 19:29:24  haselbac
! Jiao: Added support for ghost cells/nodes in unstructured meshes.
!
! Revision 1.12  2002/10/30 22:10:20  jiao
! Split volume data into more descriptive dataitems.
!
! Revision 1.11  2002/10/19 00:40:30  jblazek
! Added utility (rflosurf) to write out surface grids for GenX.
!
! Revision 1.10  2002/10/18 16:49:19  jblazek
! Changed parameter lists to some GenX routines.
!
! Revision 1.9  2002/10/15 23:22:59  jblazek
! dded Rocturb to GenX compilation path.
!
! Revision 1.8  2002/10/15 21:09:25  jiao
! Back to number of dummy cells again ...
!
! Revision 1.7  2002/10/03 21:33:48  jblazek
! Init. of bcflag moved from SendBoundaryValues to InitGenxInterface.
!
! Revision 1.6  2002/10/01 00:05:30  jblazek
! Removed st2 again.
!
! Revision 1.5  2002/09/27 21:28:47  jblazek
! Some more modifications regarding the interface to GenX.
!
! Revision 1.4  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.3  2002/09/25 17:44:56  jblazek
! Added dependent variables to volume window.
!
! Revision 1.2  2002/09/24 23:18:01  jblazek
! Changed bcflag to a pointer.
!
! Revision 1.1  2002/09/20 22:22:34  jblazek
! Finalized integration into GenX.
!
!******************************************************************************







