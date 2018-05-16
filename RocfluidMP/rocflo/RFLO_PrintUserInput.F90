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
! Purpose: write out user input for checking purposes.
!
! Description: none.
!
! Input: regions = user input.
!
! Output: to standard output.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_PrintUserInput.F90,v 1.16 2009/08/12 04:15:58 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_PrintUserInput( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY : PLAG_PrintUserInput
#endif
#ifdef PEUL
  USE ModInterfacesEulerian, ONLY : PEUL_PrintUserInput
#endif
#ifdef PERI
  USE ModInterfacesPeriodic, ONLY : PERI_PrintUserInput
#endif
#ifdef RADI
  USE ModInterfacesRadiation, ONLY : RADI_PrintUserInput
#endif
#ifdef SPEC
  USE ModInterfacesSpecies, ONLY : SPEC_PrintUserInput
#endif
#ifdef TURB
  USE ModInterfacesTurbulence, ONLY : TURB_PrintUserInput
#endif
#ifdef INRT
  USE ModInterfacesInteract, ONLY : INRT_PrintUserInput, INRT_PrintMaterialInput
#endif
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iProc

! ... local variables
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_PrintUserInput',&
  'RFLO_PrintUserInput.F90' )

! start -----------------------------------------------------------------------

  WRITE(STDOUT,1000) SOLVER_NAME

! region mapping

  WRITE(STDOUT,1005) SOLVER_NAME//' Region mapping:'
  IF(global%nRegionsProc > 0) THEN
     WRITE(STDOUT,1010) SOLVER_NAME//' nRegions/Processor:',global%nRegionsProc
  ENDIF
  IF (global%verbLevel.gt.VERBOSE_MED) then
     DO iProc=0,global%nProcAlloc-1
        WRITE(STDOUT,'(/,A,I5,A)',advance='no') SOLVER_NAME//' proc. ',iProc,': '
        DO iReg=1,global%nRegions
           IF (regions(iReg)%procid == iProc) THEN
              WRITE(STDOUT,'(1X,I5)',advance='no') iReg
           ENDIF
        ENDDO
     ENDDO
     WRITE(STDOUT,'(/)')
  ENDIF
! reference values

  WRITE(STDOUT,1005) SOLVER_NAME//' Reference values:'
  WRITE(STDOUT,1020) SOLVER_NAME//'   absvel ',global%refVelocity
  WRITE(STDOUT,1020) SOLVER_NAME//'   press  ',global%refPressure
  WRITE(STDOUT,1020) SOLVER_NAME//'   dens   ',global%refDensity
  WRITE(STDOUT,1020) SOLVER_NAME//'   cp     ',global%refCp
  WRITE(STDOUT,1020) SOLVER_NAME//'   gamma  ',global%refGamma
  WRITE(STDOUT,1020) SOLVER_NAME//'   length ',global%refLength
  WRITE(STDOUT,1020) SOLVER_NAME//'   renum  ',global%refREnum
  WRITE(STDOUT,1020) SOLVER_NAME//'   visc   ',global%refVisc
  WRITE(STDOUT,1020) SOLVER_NAME//'   prlam  ',global%prLam
  WRITE(STDOUT,1020) SOLVER_NAME//'   prturb ',global%prTurb
  WRITE(STDOUT,1020) SOLVER_NAME//'   scnlam ',global%scnLam
  WRITE(STDOUT,1020) SOLVER_NAME//'   scnturb',global%scnTurb

! acceleration

  WRITE(STDOUT,1005) SOLVER_NAME//' Acceleration terms:'
  IF (global%accelOn) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   terms on'
    WRITE(STDOUT,1020) SOLVER_NAME//'   accel-x ',global%accelX
    WRITE(STDOUT,1020) SOLVER_NAME//'   accel-y ',global%accelY
    WRITE(STDOUT,1020) SOLVER_NAME//'   accel-z ',global%accelZ
  ELSE
    WRITE(STDOUT,1030) SOLVER_NAME//'   terms off'
  ENDIF

! materials

#ifdef INRT
  CALL INRT_PrintMaterialInput(global)
#endif

! time stepping

  WRITE(STDOUT,1005) SOLVER_NAME//' Time stepping:'
  IF (global%solverType == SOLV_EXPLICIT) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   explicit scheme'
  ELSE
    WRITE(STDOUT,1030) SOLVER_NAME//'   implicit scheme'
  ENDIF
  IF (global%flowType == FLOW_UNSTEADY) THEN
    WRITE(STDOUT,1020) SOLVER_NAME//'   timestep',global%dtImposed
    WRITE(STDOUT,1020) SOLVER_NAME//'   time    ',global%timeStamp
    WRITE(STDOUT,1020) SOLVER_NAME//'   maxtime ',global%maxTime
    WRITE(STDOUT,1020) SOLVER_NAME//'   writime ',global%writeTime
    WRITE(STDOUT,1020) SOLVER_NAME//'   prntime ',global%printTime
    IF (global%solverType == SOLV_IMPLICIT) THEN
      WRITE(STDOUT,1015) SOLVER_NAME//'   order   ',global%tstepOrder
      WRITE(STDOUT,1015) SOLVER_NAME//'   subiter ',global%maxSubIter
      WRITE(STDOUT,1020) SOLVER_NAME//'   tolsub  ',global%tolSubIter
      IF (global%predictSol) THEN
        WRITE(STDOUT,1030) SOLVER_NAME//'   predict  = yes'
      ELSE
        WRITE(STDOUT,1030) SOLVER_NAME//'   predict  = no'
      ENDIF
      IF (global%dtFixed) THEN
        WRITE(STDOUT,1030) SOLVER_NAME//'   fixed dt = yes'
      ELSE
        WRITE(STDOUT,1030) SOLVER_NAME//'   fixed dt = no'
      ENDIF
    ENDIF
    WRITE(STDOUT,1015) SOLVER_NAME//'   rkScheme',global%rkScheme
  ELSE
    WRITE(STDOUT,1015) SOLVER_NAME//'   iter   ',global%currentIter
    WRITE(STDOUT,1015) SOLVER_NAME//'   maxiter',global%maxIter
    WRITE(STDOUT,1020) SOLVER_NAME//'   restol ',global%resTol
    WRITE(STDOUT,1015) SOLVER_NAME//'   wriiter',global%writeIter
    WRITE(STDOUT,1015) SOLVER_NAME//'   prniter',global%printIter
  ENDIF

! grid motion

  WRITE(STDOUT,1005) SOLVER_NAME//' Grid motion:'
  IF (global%moveGridScheme == MOVEGRID_BLOCKS) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   block-TFI scheme'
  ELSEIF (global%moveGridScheme == MOVEGRID_GLOBAL) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   block-wght-Laplacian scheme'
    WRITE(STDOUT,1015) SOLVER_NAME//'   niter ',global%moveGridNiter
    WRITE(STDOUT,1020) SOLVER_NAME//'   power ',global%moveGridPower
  ELSEIF (global%moveGridScheme == MOVEGRID_FRAME) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   global-wght-Laplacian scheme'
    WRITE(STDOUT,1015) SOLVER_NAME//'   niter     ',global%moveGridNiter
    WRITE(STDOUT,1020) SOLVER_NAME//'   amplifx   ',global%moveGridAmplifX
    WRITE(STDOUT,1020) SOLVER_NAME//'   amplify   ',global%moveGridAmplifY
    WRITE(STDOUT,1020) SOLVER_NAME//'   amplifz   ',global%moveGridAmplifZ
    WRITE(STDOUT,1020) SOLVER_NAME//'   power     ',global%moveGridPower
    WRITE(STDOUT,1015) SOLVER_NAME//'   neighbors ',global%moveGridNbour
    WRITE(STDOUT,1015) SOLVER_NAME//'   orthodir  ',global%moveGridOrthDir
    WRITE(STDOUT,1020) SOLVER_NAME//'   orthowghtx',global%moveGridOrthWghtX
    WRITE(STDOUT,1020) SOLVER_NAME//'   orthowghtY',global%moveGridOrthWghtY
    WRITE(STDOUT,1020) SOLVER_NAME//'   orthowghtZ',global%moveGridOrthWghtZ
    WRITE(STDOUT,1015) SOLVER_NAME//'   nsurfmatch',global%moveGridNsmatch
  ELSEIF (global%moveGridScheme == MOVEGRID_FOMS) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   global-orthogonal-Laplacian scheme'
    WRITE(STDOUT,1015) SOLVER_NAME//'   niter     ',global%moveGridNiter
    WRITE(STDOUT,1020) SOLVER_NAME//'   weight    ',global%moveGridWeight
    WRITE(STDOUT,1020) SOLVER_NAME//'   amplifx   ',global%moveGridAmplifX
    WRITE(STDOUT,1020) SOLVER_NAME//'   amplify   ',global%moveGridAmplifY
    WRITE(STDOUT,1020) SOLVER_NAME//'   amplifz   ',global%moveGridAmplifZ
    WRITE(STDOUT,1020) SOLVER_NAME//'   power     ',global%moveGridPower
    WRITE(STDOUT,1015) SOLVER_NAME//'   neighbors ',global%moveGridNbour
    WRITE(STDOUT,1015) SOLVER_NAME//'   orthodir  ',global%moveGridOrthDir
    WRITE(STDOUT,1020) SOLVER_NAME//'   orthowghtx',global%moveGridOrthWghtX
    WRITE(STDOUT,1020) SOLVER_NAME//'   orthowghtY',global%moveGridOrthWghtY
    WRITE(STDOUT,1020) SOLVER_NAME//'   orthowghtZ',global%moveGridOrthWghtZ
    WRITE(STDOUT,1020) SOLVER_NAME//'   weight    ',global%moveGridWeight
    WRITE(STDOUT,1020) SOLVER_NAME//'   orthocell ',global%moveGridOrthCell
    WRITE(STDOUT,1015) SOLVER_NAME//'   nsurfmatch',global%moveGridNsmatch
  ELSEIF (global%moveGridScheme == MOVEGRID_ELGLOBAL) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   block-elliptic-PDE scheme'
    WRITE(STDOUT,1015) SOLVER_NAME//'   niter     ',global%moveGridNiter
    WRITE(STDOUT,1015) SOLVER_NAME//'   viter     ',global%moveGridViter
    WRITE(STDOUT,1015) SOLVER_NAME//'   siter     ',global%moveGridSiter
    WRITE(STDOUT,1020) SOLVER_NAME//'   power     ',global%moveGridPower
  ELSEIF (global%moveGridScheme == MOVEGRID_ELFRAME) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   global-elliptic-PDE scheme'
    WRITE(STDOUT,1015) SOLVER_NAME//'   niter     ',global%moveGridNiter
    WRITE(STDOUT,1015) SOLVER_NAME//'   viter     ',global%moveGridViter
    WRITE(STDOUT,1015) SOLVER_NAME//'   siter     ',global%moveGridSiter
    WRITE(STDOUT,1020) SOLVER_NAME//'   amplifx   ',global%moveGridAmplifX
    WRITE(STDOUT,1020) SOLVER_NAME//'   amplify   ',global%moveGridAmplifY
    WRITE(STDOUT,1020) SOLVER_NAME//'   amplifz   ',global%moveGridAmplifZ
    WRITE(STDOUT,1020) SOLVER_NAME//'   power     ',global%moveGridPower
    WRITE(STDOUT,1015) SOLVER_NAME//'   neighbors ',global%moveGridNbour
    WRITE(STDOUT,1015) SOLVER_NAME//'   orthodir  ',global%moveGridOrthDir
    WRITE(STDOUT,1020) SOLVER_NAME//'   orthowghtx',global%moveGridOrthWghtX
    WRITE(STDOUT,1020) SOLVER_NAME//'   orthowghtY',global%moveGridOrthWghtY
    WRITE(STDOUT,1020) SOLVER_NAME//'   orthowghtZ',global%moveGridOrthWghtZ
    WRITE(STDOUT,1015) SOLVER_NAME//'   nsurfmatch',global%moveGridNsmatch
  ELSEIF (global%moveGridScheme == MOVEGRID_VMS) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//'   block-VMS scheme'
    WRITE(STDOUT,1015) SOLVER_NAME//'   niter ',global%moveGridNiter
    WRITE(STDOUT,1020) SOLVER_NAME//'   weight',global%moveGridWeight
  ENDIF

#ifdef STATS
! statistics

  IF (global%doStat == ACTIVE) THEN
    IF (global%reStat == ACTIVE) THEN
      WRITE(STDOUT,1005) SOLVER_NAME//' Continued statistics:'
    ELSE
      WRITE(STDOUT,1005) SOLVER_NAME//' New statistics:'
    ENDIF
    IF (global%mixtNStat > 0) THEN
      WRITE(STDOUT,*) SOLVER_NAME//'  mixture NStat  =',global%mixtNStat
      WRITE(STDOUT,*) SOLVER_NAME//'  mixture StatID =', &
      global%mixtStatID(1,:)*10+global%mixtStatID(2,:)
    ENDIF
    IF (global%turbNStat > 0) THEN
      WRITE(STDOUT,*) SOLVER_NAME//'  TURB NStat     =',global%turbNStat
      WRITE(STDOUT,*) SOLVER_NAME//'  TURB StatID    =', &
      global%turbStatID(1,:)*10+global%turbStatID(2,:)
    ENDIF
  ENDIF
#endif

#ifdef PERI
  IF (regions(1)%periInput%flowKind /= OFF) THEN
    CALL PERI_PrintUserInput( regions(1) )
  ENDIF
#endif

! repeat for all regions ------------------------------------------------------


  IF (global%verbLevel.gt.VERBOSE_MED) then
     DO iReg=1,global%nRegions
        WRITE(STDOUT,1025) SOLVER_NAME,iReg
        
        ! - dimensions
        
        WRITE(STDOUT,1005) SOLVER_NAME//'   Dimensions:'
        WRITE(STDOUT,1015) SOLVER_NAME//'     icells  ',regions(iReg)%levels(1)%grid%ipc
        WRITE(STDOUT,1015) SOLVER_NAME//'     jcells  ',regions(iReg)%levels(1)%grid%jpc
        WRITE(STDOUT,1015) SOLVER_NAME//'     kcells  ',regions(iReg)%levels(1)%grid%kpc
        WRITE(STDOUT,1010) SOLVER_NAME//'     dumcells',regions(iReg)%nDumCells
        
        ! - flow model
        
        WRITE(STDOUT,1005) SOLVER_NAME//'   Flow type, turbulence:'
        IF (regions(iReg)%mixtInput%flowModel == FLOW_EULER) THEN
           WRITE(STDOUT,1030) SOLVER_NAME//'     flow model  = Euler'
        ELSE
           WRITE(STDOUT,1030) SOLVER_NAME//'     flow model  = Navier-Stokes'
        ENDIF
        IF (regions(iReg)%mixtInput%moveGrid) THEN
           WRITE(STDOUT,1030) SOLVER_NAME//'     moving grid = yes'
        ELSE
           WRITE(STDOUT,1030) SOLVER_NAME//'     moving grid = no'
        ENDIF
        
        IF (regions(iReg)%mixtInput%flowModel == FLOW_NAVST) THEN
           IF (regions(iReg)%mixtInput%turbModel == TURB_MODEL_NONE) THEN
              WRITE(STDOUT,1030) SOLVER_NAME//'     turb. model = none'
           ELSE
#ifdef TURB
              CALL TURB_PrintUserInput( regions(iReg) )
#endif
           ENDIF
        ENDIF
        
        ! - viscosity model
        
        IF (regions(iReg)%mixtInput%computeTv) THEN
           WRITE(STDOUT,1005) SOLVER_NAME//'   Viscosity model:'
           
           IF (regions(iReg)%mixtInput%viscModel == VISC_SUTHR) THEN
              WRITE(STDOUT,1030) SOLVER_NAME//'     viscosity model        = Sutherland'
              WRITE(STDOUT,1020) SOLVER_NAME//'     reference viscosity   ',&
                   regions(iReg)%mixtInput%refVisc
              WRITE(STDOUT,1020) SOLVER_NAME//'     reference temperature ',&
                   regions(iReg)%mixtInput%refTemp
              WRITE(STDOUT,1020) SOLVER_NAME//'     Sutherland coefficient',&
                   regions(iReg)%mixtInput%suthCoef
              
           ELSEIF (regions(iReg)%mixtInput%viscModel == VISC_FIXED) THEN
              WRITE(STDOUT,1030) SOLVER_NAME//'     viscosity model  = Fixed'
              WRITE(STDOUT,1020) SOLVER_NAME//'     viscosity value ',&
                   regions(iReg)%mixtInput%refVisc
              
           ELSEIF (regions(iReg)%mixtInput%viscModel == VISC_ANTIB) THEN
              WRITE(STDOUT,1030) SOLVER_NAME//'     viscosity model  = Antibes'
              WRITE(STDOUT,1020) SOLVER_NAME//'     reference viscosity',&
                   regions(iReg)%mixtInput%refVisc
              
           ENDIF ! viscModel
        ENDIF   ! computeTv
        
        ! - multi-physics modules
        
        WRITE(STDOUT,1005) SOLVER_NAME//'   Multi-physics modules:'
        
        IF (regions(iReg)%mixtInput%gasModel == GAS_MODEL_TCPERF) THEN
           WRITE(STDOUT,1030) SOLVER_NAME//'     species model  = none'
        ELSE
#ifdef SPEC
           CALL SPEC_PrintUserInput( regions(iReg) )
#endif
        ENDIF
        
        IF (global%peulUsed) THEN
#ifdef PEUL
           WRITE(STDOUT,1030) SOLVER_NAME//'     con.part. used = yes'
           CALL PEUL_PrintUserInput( regions(iReg) )
#endif
        ELSE
           WRITE(STDOUT,1030) SOLVER_NAME//'     con.part. used = no'
        ENDIF
        
        IF (global%plagUsed) THEN
#ifdef PLAG
           WRITE(STDOUT,1030) SOLVER_NAME//'     dis.part. used = yes'
           CALL PLAG_PrintUserInput( regions(iReg) )
#endif
        ELSE
           WRITE(STDOUT,1030) SOLVER_NAME//'     dis.part. used = no'
        ENDIF
        
        IF (regions(iReg)%mixtInput%radiUsed) THEN
#ifdef RADI
           CALL RADI_PrintUserInput( regions(iReg) )
#endif
        ELSE
           WRITE(STDOUT,1030) SOLVER_NAME//'     radiation used = no'
        ENDIF
        
#ifdef INRT
        IF (global%inrtUsed) THEN
           CALL INRT_PrintUserInput( regions(iReg) )
        ENDIF
#endif
        
        ! - numerics
        
        WRITE(STDOUT,1005) SOLVER_NAME//'   Numerics:'
        WRITE(STDOUT,1010) SOLVER_NAME//'     levels  ',regions(iReg)%nGridLevels
        IF (regions(iReg)%mixtInput%timeScheme == TST_HYB5RK) THEN
           WRITE(STDOUT,1030) SOLVER_NAME//'     timedis  = explicit multistage'
        ELSE IF (regions(iReg)%mixtInput%timeScheme == TST_STD4RK) THEN
           WRITE(STDOUT,1030) SOLVER_NAME//'     timedis  = explicit classical Runge-Kutta'
        ENDIF
        WRITE(STDOUT,1020) SOLVER_NAME//'     CFL     ',regions(iReg)%mixtInput%cfl
        WRITE(STDOUT,1020) SOLVER_NAME//'     smoocf  ',regions(iReg)%mixtInput%smoocf
        WRITE(STDOUT,1010) SOLVER_NAME//'     discr   ',regions(iReg)%mixtInput%spaceDiscr
        IF (regions(iReg)%mixtInput%spaceDiscr == DISCR_CEN_SCAL) THEN
           WRITE(STDOUT,1020) SOLVER_NAME//'     k2      ',regions(iReg)%mixtInput%vis2
           WRITE(STDOUT,1020) SOLVER_NAME//'     1/k4    ',1./regions(iReg)%mixtInput%vis4
           IF (regions(iReg)%mixtInput%pSwitchType == PSWITCH_STD) THEN
              WRITE(STDOUT,1030) SOLVER_NAME//'     pswitch  = standard'
           ELSE
              WRITE(STDOUT,1030) SOLVER_NAME//'     pswitch  = mixed TVD and standard'
              WRITE(STDOUT,1020) SOLVER_NAME//'     pswOmega',regions(iReg)%mixtInput%pSwitchOmega
           ENDIF
        ELSE
           WRITE(STDOUT,1010) SOLVER_NAME//'     order   ',regions(iReg)%mixtInput%spaceOrder
           WRITE(STDOUT,1020) SOLVER_NAME//'     limfac  ',regions(iReg)%mixtInput%limfac
           WRITE(STDOUT,1020) SOLVER_NAME//'     entropy ',regions(iReg)%mixtInput%epsentr
        ENDIF
        IF (regions(iReg)%mixtInput%faceEdgeAvg == FE_AVG_UNIFORM) THEN
           WRITE(STDOUT,1030) SOLVER_NAME//'     f/e avg  = uniform'
        ELSE
           WRITE(STDOUT,1030) SOLVER_NAME//'     f/e avg  = grid-dependent linear'
        ENDIF
        
     ENDDO   ! iReg
  ENDIF
! finish ----------------------------------------------------------------------

  WRITE(STDOUT,1035) SOLVER_NAME

  CALL DeregisterFunction( global )

1000 FORMAT(/,A,1X,80('-'))
1005 FORMAT(/,A)
1010 FORMAT(A,' = ',I2)
1015 FORMAT(A,' = ',I8)
1020 FORMAT(A,' = ',E12.5)
1025 FORMAT(/,A,' Region ',I6,':')
1030 FORMAT(A)
1035 FORMAT(/,A,1X,80('-'),/)

END SUBROUTINE RFLO_PrintUserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_PrintUserInput.F90,v $
! Revision 1.16  2009/08/12 04:15:58  mtcampbe
! Major update, bugfix from Abe development, more propagation compatibility,
! some Rocstar IO changes, Ju's temporary clipping fix for turbulence. A bug
! fix for initialization IO.
!
! Revision 1.15  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.14  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.13  2006/05/09 23:33:44  wasistho
! added DTFIXED in timestep section
!
! Revision 1.12  2006/04/27 01:55:14  wasistho
! added screenprint neigbors
!
! Revision 1.11  2006/03/18 13:26:45  wasistho
! added orthDir and orthWghtX,Y,Z
!
! Revision 1.10  2006/03/08 06:34:38  wasistho
! added movegrid_elglobal and elframe
!
! Revision 1.9  2005/11/18 07:21:01  wasistho
! rearranged gridmotion screen print
!
! Revision 1.8  2005/10/31 21:09:36  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.7  2005/10/28 22:47:24  wasistho
! print FOMS orthocell
!
! Revision 1.6  2005/10/28 05:42:44  wasistho
! printout FOMS gridmotion
!
! Revision 1.5  2005/08/28 23:49:05  wasistho
! added orthoWght for block orthogonality of RFLO global-gridmotion
!
! Revision 1.4  2005/08/18 19:49:19  wasistho
! added print NSURFMATCH
!
! Revision 1.3  2005/06/04 01:02:20  wasistho
! distinguished to AMPLIFX,Y,Z
!
! Revision 1.2  2005/06/02 22:59:54  wasistho
! added moveGridAmplif and moveGridPower
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.34  2004/11/17 16:30:44  haselbac
! Added printing of rkScheme
!
! Revision 1.33  2004/09/02 02:59:42  wasistho
! screen output face-edge averaging option
!
! Revision 1.32  2004/03/05 22:09:02  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.31  2004/03/02 21:49:22  jferry
! Added inrtUsed flag to mixture data structure
!
! Revision 1.30  2003/11/20 16:40:40  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.26  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.25  2003/09/26 21:44:28  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.24  2003/08/28 20:37:43  wasistho
! excluced ModInterfacesTurbulence,Radiation,Periodic from ModInterfaces
!
! Revision 1.23  2003/08/28 20:05:39  jblazek
! Added acceleration terms.
!
! Revision 1.22  2003/08/11 21:51:18  jblazek
! Added basic global grid smoothing scheme.
!
! Revision 1.21  2003/07/03 21:48:45  jblazek
! Implemented dual-time stepping.
!
! Revision 1.20  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.19  2003/04/10 23:31:20  fnajjar
! Added printouts for viscosity models
!
! Revision 1.18  2003/03/29 03:30:18  wasistho
! install ROCPERI
!
! Revision 1.17  2003/03/24 23:27:01  jferry
! converted PrintMaterialInput to INRT_PrintMaterialInput
!
! Revision 1.16  2003/03/11 16:04:19  jferry
! Enclosed USE statements for multi-physics routines within ifdefs
!
! Revision 1.15  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
! Revision 1.14  2002/11/04 22:30:37  wasistho
! Put stats print userinput within ifdef STATS
!
! Revision 1.13  2002/11/02 01:57:21  wasistho
! Added TURB statistics
!
! Revision 1.12  2002/09/17 22:51:23  jferry
! Removed Fast Eulerian particle type
!
! Revision 1.11  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.10  2002/07/25 00:40:25  jblazek
! Option for TVD type pressure switch.
!
! Revision 1.9  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.8  2002/06/22 01:13:38  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.7  2002/02/27 18:38:20  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
! Revision 1.6  2002/02/25 22:36:53  jblazek
! Simplified solver initialization routine.
!
! Revision 1.5  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.3  2002/01/08 22:09:17  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.2  2001/12/19 23:09:22  jblazek
! Added routines to read grid and solution.
!
! Revision 1.1  2001/12/11 21:59:29  jblazek
! memory allocation added.
!
!******************************************************************************







