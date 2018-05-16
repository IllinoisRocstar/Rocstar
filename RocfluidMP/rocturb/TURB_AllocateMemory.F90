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
! Purpose: Allocate memory for all variables associated with turbulence
!          for all active regions on current processor.
!
! Description: none.
!
! Input: region = info of current region data
!
! Output: region%turb = info of current turbulence data
!
! Notes: The total size of memory depends on the input selected by user.
!
!******************************************************************************
!
! $Id: TURB_AllocateMemory.F90,v 1.17 2009/08/31 05:17:35 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_AllocateMemory( region ) ! PUBLIC

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModTurbulence, ONLY : t_turb
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset, &
                            RFLO_GetDimensDummy, RFLO_GetCellOffset
#include "Indexing.h"
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), TARGET :: region

! ... loop variables
  INTEGER :: iPatch, ijBeg, ijEnd, ibc, iec, ibn, ien
#ifdef RFLO
  INTEGER :: iLev
#endif

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_turb), POINTER   :: turb
  TYPE(t_patch), POINTER  :: patch1, patch

  INTEGER :: turbModel, modelClass, filterType, errFl
  INTEGER :: nCv, nDv, nSv, nSt, nGrad, nZof, filterWidth(DIRI:DIRK)
#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: n1, n2, iOff, iCOff, ijCOff, iNOff, ijNOff
  INTEGER :: homDir(DIRI:DIRK)
#endif
#ifdef RFLU
  INTEGER :: nPatches, nCellsTot, nFaces, nFacesTot, nBFaces, nBFacesTot
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_AllocateMemory.F90,v $ $Revision: 1.17 $'

  global => region%global
  CALL RegisterFunction( global,'TURB_AllocateMemory',&
  'TURB_AllocateMemory.F90' )

! get variables ---------------------------------------------------------------

  turbModel      = region%mixtInput%turbModel
  modelClass     = region%turbInput%modelClass
  nCv            = region%turbInput%nCv
  nDv            = region%turbInput%nDv
  nSv            = region%turbInput%nSv
  nSt            = region%turbInput%nSt
  nGrad          = region%turbInput%nGrad
  nZof           = region%turbInput%nZof
  filterWidth(:) = region%turbInput%filterWidth(:)  ! RFLU only uses 1st comp.

#ifdef RFLO
  filterType     = region%turbInput%filterType
  homDir(:)      = region%turbInput%homDir(:)

! loop over all grid levels
  
  DO iLev=1,region%nGridLevels

    turb => region%levels(iLev)%turb

! - get cell and node dimensions ----------------------------------------------

    CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
    ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
    iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

    CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
    ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
    ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)
#endif
#ifdef RFLU
! - get variables and pointers ------------------------------------------------
    turb => region%turb

! - get cell and node dimensions
    nCellsTot = region%grid%nCellsTot
    nFaces    = region%grid%nFaces
    nFacesTot = region%grid%nFacesTot
    ibc       = 1
    iec       = nCellsTot
    ibn       = 1
    ien       = nFaces
    nPatches  = region%grid%nPatches

    nBFaces    = 0
    nBFacesTot = 0

    DO iPatch = 1,nPatches
      patch => region%patches(iPatch)

      nBFaces    = nBFaces    + patch%nBTris    + patch%nBQuads
      nBFacesTot = nBFacesTot + patch%nBTrisTot + patch%nBQuadsTot
    END DO ! iPatch
#endif

! - RANS/DES ----------------------------------------------------------------- 
! - RANS transport equations variables: cv, cvOld, rhs, rhsSum, diss, lens,
!                                       srad, epsIrs (only RFLO)

    IF ((modelClass==MODEL_RANS) .AND. (nCv > 0)) THEN
      ALLOCATE( turb%cv(    nCv,ibc:iec),stat=errFl )
      IF (errFl>0) GOTO 88
      ALLOCATE( turb%cvOld( nCv,ibc:iec),stat=errFl )
      IF (errFl>0) GOTO 88
      ALLOCATE( turb%rhs(   nCv,ibc:iec),stat=errFl )
      IF (errFl>0) GOTO 88
      ALLOCATE( turb%rhsSum(nCv,ibc:iec),stat=errFl )
      IF (errFl>0) GOTO 88
      ALLOCATE( turb%diss(  nCv,ibc:iec),stat=errFl )
      IF (errFl>0) GOTO 88
      ALLOCATE( turb%dsterm(nCv,ibc:iec),stat=errFl )
      IF (errFl>0) GOTO 88
      ALLOCATE( turb%lens(      ibc:iec),stat=errFl )
      IF (errFl>0) GOTO 88
      IF (global%solverType==SOLV_IMPLICIT) THEN  ! Dual Tst
         ALLOCATE( turb%cvn(  nCv,ibc:iec),stat=errFl )
         IF (errFl>0) GOTO 88
         ALLOCATE( turb%cvn1( nCv,ibc:iec),stat=errFl )
         IF (errFl>0) GOTO 88
         ALLOCATE( turb%cvn2( nCv,ibc:iec),stat=errFl )
         IF (errFl>0) GOTO 88
         ALLOCATE( turb%sDual(nCv,ibc:iec),stat=errFl ); 
         IF (errFl>0) GOTO 88
      ELSE
         NULLIFY( turb%cvn,turb%cvn1,turb%cvn2,turb%sDual )
!         NULLIFY( turb%cv,turb%cvOld,turb%rhs,turb%rhsSum,turb%diss)
      ENDIF  ! implicit
   ELSE
      NULLIFY( turb%cvn,turb%cvn1,turb%cvn2,turb%sDual )
      NULLIFY( turb%cv,turb%cvOld,turb%rhs,turb%rhsSum,turb%diss)
   ENDIF  ! rans
   
#ifdef RFLO
   IF ((modelClass==MODEL_RANS) .AND. (nCv > 0)) THEN
       ALLOCATE( turb%srad(DIRI:DIRK,ibc:iec),stat=errFl )
       IF (errFl>0) GOTO 88
       IF ((global%flowType==FLOW_STEADY) .AND. &
            (region%turbInput%smoocf > 0._RFREAL)) THEN
          ALLOCATE( turb%epsIrs(DIRI:DIRK,ibc:iec), stat=errFl )
          IF (errFl>0) GOTO 88
       ELSE
          NULLIFY( turb%epsIrs ) 
       ENDIF
    ELSE
       NULLIFY( turb%srad,turb%epsIrs )
    ENDIF
#endif

! - RANS gradient variables
! - note: put within MODEL_RANS for not allocated twice by LES as they may be
!         allocated somewhere as LES work variables to save memory

    IF (modelClass==MODEL_RANS) THEN 
       IF (nGrad > 0) THEN
#ifdef RFLO
          ALLOCATE( turb%gradi(nGrad,ibn:ien),stat=errFl )
          IF (errFl>0) GOTO 88  
          ALLOCATE( turb%gradj(nGrad,ibn:ien),stat=errFl )
          IF (errFl>0) GOTO 88  
          ALLOCATE( turb%gradk(nGrad,ibn:ien),stat=errFl )
          IF (errFl>0) GOTO 88
       ELSE
          NULLIFY( turb%gradi,turb%gradj,turb%gradk )
#endif
#ifdef RFLU
          ALLOCATE( turb%gradi( XCOORD:ZCOORD,nGrad, nFaces),stat=errFl )
          IF (errFl>0) GOTO 88
          ALLOCATE( turb%bGradi(XCOORD:ZCOORD,nGrad,nBFaces),stat=errFl )
          IF (errFl>0) GOTO 88
       ELSE
          NULLIFY( turb%gradi, turb%bGradi )
#endif
       ENDIF
    ENDIF
    
! - total number of no-slip faces for wall dist., only on finest grid
    
#ifdef RFLO
    IF (iLev == 1) THEN
       DO iPatch=1,region%nPatches
          patch  => region%levels(iLev)%patches(iPatch)
          
          IF ((patch%bcType>=BC_NOSLIPWALL .AND. &
               patch%bcType<=BC_NOSLIPWALL+BC_RANGE) .OR. &
               (patch%bcType>=BC_INJECTION .AND. &
               patch%bcType<=BC_INJECTION+BC_RANGE)) THEN
             
             n1    = ABS(patch%l1end-patch%l1beg)
             n2    = ABS(patch%l2end-patch%l2beg)
             iOff  = n1 + 1
             ijBeg = IndIJ( 0, 0,iOff)
             ijEnd = IndIJ(n1,n2,iOff)
             global%turbWallDim = global%turbWallDim + ijEnd-ijBeg+1
          ENDIF    ! bcType 
       ENDDO      ! iPatch
    ENDIF        ! iLev
#endif
#ifdef RFLU
    DO iPatch=1,region%grid%nPatches
       patch  => region%patches(iPatch)
       
       IF ((patch%bcType>=BC_NOSLIPWALL .AND. &
            patch%bcType<=BC_NOSLIPWALL+BC_RANGE) .OR. &
            (patch%bcType>=BC_INJECTION .AND. &
            patch%bcType<=BC_INJECTION+BC_RANGE)) THEN
          
          global%turbWallDim = global%turbWallDim + patch%nBFaces
       ENDIF    ! bcType 
    ENDDO      ! iPatch
#endif

! - GENERAL -------------------------------------------------------------------
! - TURB derived variables
    
    IF (nDv > 0) THEN
       ALLOCATE( turb%dv(nDv,ibc:iec),stat=errFl )
       IF (errFl>0) GOTO 88
    ELSE
       NULLIFY( turb%dv )
    ENDIF
    
! - turbulent stress components
    
    IF (nSv > 0) THEN
       ALLOCATE( turb%sv(nSv,ibc:iec),stat=errFl )
       IF (errFl>0) GOTO 88
    ELSE
       NULLIFY( turb%sv )
    ENDIF
    
! - vorticity components
    
    IF ((region%turbInput%calcVort /= CALCVORT_NO) .OR. &
         (modelClass == MODEL_RANS)) THEN
       ALLOCATE( turb%vort(XCOORD:ZCOORD,ibc:iec),stat=errFl )
       IF (errFl>0) GOTO 88
    ELSE
       NULLIFY( turb%vort )
    ENDIF
    
! - zero-one switch field
    
    IF (nZof > 0) THEN
#ifdef RFLO
       ALLOCATE( turb%zofi(XCOORD:ZCOORD,nZof,ibn:ien),stat=errFl )
       IF (errFl>0) GOTO 88
       ALLOCATE( turb%zofj(XCOORD:ZCOORD,nZof,ibn:ien),stat=errFl )
       IF (errFl>0) GOTO 88
       ALLOCATE( turb%zofk(XCOORD:ZCOORD,nZof,ibn:ien),stat=errFl )
       IF (errFl>0) GOTO 88
#endif
#ifdef RFLU
       ALLOCATE( turb%zofi(XCOORD:ZCOORD,nZof,nFaces),stat=errFl )
       IF (errFl>0) GOTO 88
       ALLOCATE( turb%bZofi(XCOORD:ZCOORD,nZof,nBFaces),stat=errFl )
       IF (errFl>0) GOTO 88
#endif
    ENDIF
    
! - statistics of TURB pertinent variables (s.u. dynamic model coefficient)
    
#ifdef STATS
    IF ((global%flowType == FLOW_UNSTEADY) .AND. &
         (global%doStat == ACTIVE) .AND. &
         (global%turbNStat > 0)) THEN
       ALLOCATE( turb%tav(global%turbNStat,ibc:iec),stat=errFl )
       IF (errFl>0) GOTO 88
    ELSE
       NULLIFY( turb%tav )
    ENDIF
    
    IF (nSt > 0) THEN
       ALLOCATE( turb%st(nSt,ibc:iec),stat=errFl )
       IF (errFl>0) GOTO 88
    ELSE
       NULLIFY( turb%st )
    ENDIF
#endif
    
! - LES -----------------------------------------------------------------------
! - LES filter coefficients, fix arrays
    
#ifdef RFLO
    IF ((turbModel==TURB_MODEL_FIXSMAG) .OR. &  ! involve eddy vis.type
         (turbModel==TURB_MODEL_DYNSMAG) .OR. &
         (turbModel==TURB_MODEL_DYNMIXD)) THEN
       ALLOCATE( turb%fvolI(ibn:ien),stat=errFl )
       IF (errFl>0) GOTO 88
       ALLOCATE( turb%fvolJ(ibn:ien),stat=errFl )
       IF (errFl>0) GOTO 88
       ALLOCATE( turb%fvolK(ibn:ien),stat=errFl )
       IF (errFl>0) GOTO 88
    ELSE
       NULLIFY( turb%fvolI, turb%fvolJ, turb%fvolK )
    ENDIF
    
    IF (((turbModel==TURB_MODEL_SCALSIM)  .OR. & ! involve filtering
         (turbModel==TURB_MODEL_DYNSMAG)  .OR. &
         (turbModel==TURB_MODEL_DYNMIXD)) .AND. &
         (filterType == FILTYPE_NONUNIF)) THEN
       
       IF (homDir(DIRI) == OFF) THEN
          IF (filterWidth(DIRI) == FILWIDTH_ONE) THEN
             ALLOCATE( turb%ccCofi1( 3,ibc:iec),stat=errFl )
             IF (errFl>0) GOTO 88
             ALLOCATE( turb%ccCofi2( 3,ibc:iec),stat=errFl )
             IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofi1I(3,ibn:ien),stat=errFl )
             IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofi1J(3,ibn:ien),stat=errFl )
             IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofi1K(3,ibn:ien),stat=errFl )
             IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofi2I(3,ibn:ien),stat=errFl )
             IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofi2J(3,ibn:ien),stat=errFl )
             IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofi2K(3,ibn:ien),stat=errFl )
             IF (errFl>0) GOTO 88
          ELSEIF ((filterWidth(DIRI) == FILWIDTH_TWO) .OR. &
               (filterWidth(DIRI) == FILWIDTH_ZERO)) THEN
             ALLOCATE( turb%ccCofi2( 3,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ccCofi4( 5,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofi2I(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofi2J(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofi2K(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofi4I(5,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofi4J(5,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofi4K(5,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
          ENDIF
       ENDIF
       IF (homDir(DIRJ) == OFF) THEN
          IF (filterWidth(DIRJ) == FILWIDTH_ONE) THEN
             ALLOCATE( turb%ccCofj1( 3,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ccCofj2( 3,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofj1I(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofj1J(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofj1K(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofj2I(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofj2J(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofj2K(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
          ELSEIF ((filterWidth(DIRJ) == FILWIDTH_TWO) .OR. &
               (filterWidth(DIRJ) == FILWIDTH_ZERO)) THEN
             ALLOCATE( turb%ccCofj2( 3,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ccCofj4( 5,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofj2I(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofj2J(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofj2K(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofj4I(5,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofj4J(5,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofj4K(5,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
          ENDIF
       ENDIF
       IF (homDir(DIRK) == OFF) THEN
          IF (filterWidth(DIRK) == FILWIDTH_ONE) THEN
             ALLOCATE( turb%ccCofk1( 3,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ccCofk2( 3,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofk1I(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofk1J(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofk1K(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofk2I(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofk2J(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofk2K(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
          ELSEIF ((filterWidth(DIRK) == FILWIDTH_TWO) .OR. &
               (filterWidth(DIRK) == FILWIDTH_ZERO)) THEN
             ALLOCATE( turb%ccCofk2( 3,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ccCofk4( 5,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofk2I(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofk2J(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofk2K(3,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofk4I(5,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofk4J(5,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
             ALLOCATE( turb%ffCofk4K(5,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
          ENDIF
       ENDIF
    ELSE
       NULLIFY( turb%ccCofi1 ,turb%ccCofi2 ,turb%ccCofi4  )
       NULLIFY( turb%ccCofj1 ,turb%ccCofj2 ,turb%ccCofj4  )
       NULLIFY( turb%ccCofk1 ,turb%ccCofk2 ,turb%ccCofk4  )
       NULLIFY( turb%ffCofi1I,turb%ffCofi2I,turb%ffCofi4I )
       NULLIFY( turb%ffCofi1J,turb%ffCofi2J,turb%ffCofi4J )
       NULLIFY( turb%ffCofi1K,turb%ffCofi2K,turb%ffCofi4K )
       NULLIFY( turb%ffCofj1I,turb%ffCofj2I,turb%ffCofj4I )
       NULLIFY( turb%ffCofj1J,turb%ffCofj2J,turb%ffCofj4J )
       NULLIFY( turb%ffCofj1K,turb%ffCofj2K,turb%ffCofj4K )
       NULLIFY( turb%ffCofk1I,turb%ffCofk2I,turb%ffCofk4I )
       NULLIFY( turb%ffCofk1J,turb%ffCofk2J,turb%ffCofk4J )
       NULLIFY( turb%ffCofk1K,turb%ffCofk2K,turb%ffCofk4K )
    ENDIF
#endif
#ifdef RFLU
    IF ((turbModel==TURB_MODEL_FIXSMAG) .OR. &  ! involve eddy vis.type
         (turbModel==TURB_MODEL_DYNSMAG) .OR. &
         (turbModel==TURB_MODEL_DYNMIXD)) THEN
       ALLOCATE( turb%fvolI(  nFaces),stat=errFl ); IF (errFl>0) GOTO 88
       ALLOCATE( turb%bfVolI(nBFaces),stat=errFl ); IF (errFl>0) GOTO 88
    ELSE
       NULLIFY( turb%fvolI, turb%bfVolI )
    ENDIF
    
    IF (modelClass == MODEL_LES) THEN
       ALLOCATE( turb%avgCoI( 2, nFaces),stat=errFl ); IF (errFl>0) GOTO 88
       ALLOCATE( turb%bAvgCoI(2,nBFaces),stat=errFl ); IF (errFl>0) GOTO 88
    ELSE
       NULLIFY( turb%avgCoI, turb%bAvgCoI )
    ENDIF
#endif
    
! - nullify unused arrays
    
    IF (modelClass == MODEL_RANS) THEN
       NULLIFY( turb%trace, turb%lij, turb%mij )
       NULLIFY( turb%fVar, turb%ffVar, turb%ccVar, turb%coef, turb%mueT )
#ifdef RFLO
       NULLIFY( turb%fISij, turb%fJSij, turb%fKSij )
#endif
#ifdef RFLU
       NULLIFY( turb%bLij, turb%bMij )
       NULLIFY( turb%bfVar, turb%bffVar, turb%bCoef, turb%bMueT )
       NULLIFY( turb%fISij, turb%bfISij )
#endif
    ENDIF
    
! - WLM (wall layer model)-----------------------------------------------------
! - allocate arrays pertinent to WLM on coarser grid level (RFLO only)
    
#ifdef RFLO
    IF (iLev > 1) THEN
       DO iPatch=1,region%nPatches
          patch1 => region%levels(1)%patches(iPatch)
          patch  => region%levels(iLev)%patches(iPatch)
          
          IF (patch%bcType>=BC_NOSLIPWALL .AND. &
               patch%bcType<=BC_NOSLIPWALL+BC_RANGE .AND. &
               patch1%valBola%switches(WLM_INPUT_MODEL) /= WLM_MODEL_NOMODEL) THEN
             patch%valBola%nData     = patch1%valBola%nData
             patch%valBola%nSwitches = patch1%valBola%nSwitches
             patch%valBola%distrib   = patch1%valBola%distrib
             
             ALLOCATE( patch%valBola%switches(patch%valBola%nSwitches), &
                  stat=errFl ); IF (errFl>0) GOTO 88
             patch%valBola%switches  = patch1%valBola%switches
             
             IF (patch%valBola%nData > 0) THEN  ! vals at level 1 are defined
                n1    = ABS(patch%l1end-patch%l1beg)
                n2    = ABS(patch%l2end-patch%l2beg)
                iOff  = n1 + 1
                ijBeg = IndIJ( 0, 0,iOff)
                ijEnd = IndIJ(n1,n2,iOff)
                
                ALLOCATE( patch%valBola%vals(ijBeg:ijEnd,patch%valBola%nData), &
                     stat=errFl ); IF (errFl>0) GOTO 88
                
             ENDIF  ! nData
          ENDIF    ! bcType 
       ENDDO      ! iPatch
    ENDIF        ! iLev
    
 ENDDO  !iLev
#endif
 
 GOTO 999
 
! finalize ----------------------------------------------------------
 
88 CONTINUE
 
 global%error = errFl
 CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
 
999 CONTINUE
 
 CALL DeregisterFunction( global )
 
END SUBROUTINE TURB_AllocateMemory

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_AllocateMemory.F90,v $
! Revision 1.17  2009/08/31 05:17:35  mtcampbe
! Fix DES allocation bug
!
! Revision 1.16  2009/04/07 18:51:50  mtcampbe
! Fixed intel compiler complaint about line continuations
!
! Revision 1.15  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.14  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.13  2006/01/12 09:49:44  wasistho
! enabled tripping fixed Smagorinsky
!
! Revision 1.12  2004/12/09 22:18:05  wasistho
! allocated added data structures
!
! Revision 1.11  2004/10/25 02:11:45  wasistho
! shift allocation of turb%st within ifdef STATS
!
! Revision 1.10  2004/09/23 22:22:58  wasistho
! added nullify turb%dsterm
!
! Revision 1.9  2004/08/04 02:49:06  wasistho
! removed turb%avgCoI,J,K as it is defined as grid%c2fCoI,J,K
!
! Revision 1.8  2004/06/03 02:10:25  wasistho
! enabled non-uniform fix-Smagorinsky
!
! Revision 1.7  2004/05/28 01:57:21  wasistho
! update unstructured grid LES
!
! Revision 1.6  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.5  2004/03/25 04:40:41  wasistho
! prepared for RFLU
!
! Revision 1.4  2004/03/20 03:28:29  wasistho
! prepared for RFLU
!
! Revision 1.3  2004/03/19 02:42:19  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/13 03:10:45  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/03/05 04:36:59  wasistho
! changed nomenclature
!
! Revision 1.15  2004/02/04 22:30:53  wasistho
! move MPIsum of global%turbWallDim from allocateMemory to calcMetrics
!
! Revision 1.14  2004/01/21 03:51:08  wasistho
! add injection wall condition for counting number of wall faces
!
! Revision 1.13  2003/10/21 20:31:23  wasistho
! added dt relaxation in steady flow due to RANS source term
!
! Revision 1.12  2003/10/07 20:31:24  wasistho
! allocate vorticities also for RaNS/DES
!
! Revision 1.10  2003/08/29 01:40:25  wasistho
! Added TARGET attribute to region variable, since pointers are cached into it
!
! Revision 1.9  2003/08/06 15:55:58  wasistho
! added vorticities computation
!
! Revision 1.8  2003/06/09 23:19:09  wasistho
! extra condition for wlm data allocation
!
! Revision 1.7  2003/05/31 01:46:14  wasistho
! installed turb. wall layer model
!
! Revision 1.6  2003/05/24 02:08:00  wasistho
! turbulence statistics expanded
!
! Revision 1.5  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.4  2002/11/02 02:05:37  wasistho
! Added TURB statistics
!
! Revision 1.3  2002/10/16 07:48:53  wasistho
! Enable Fix Smagorinsky
!
! Revision 1.2  2002/10/16 01:58:44  wasistho
! Changed global%error flag
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************












