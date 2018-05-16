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
! Purpose: loads data buffer size for patches
!          and shrinks particle datastructure. 
!
! Description: none.
!
! Input: region = current region,
!        iReg   = current region number.
!
! Output: regions(iReg)%levels%patch%buffPlag%aiv,arv,cv,dv,tv = buffer data.
!         regions(iReg)%levels%plag%aiv,arv,cv,dv,tv           = Plag data.
!
! Notes:
!
!   The cell index for the particle has been updated to its position
!
!******************************************************************************
!
! $Id: PLAG_PatchLoadDataBuffers.F90,v 1.3 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_PatchLoadDataBuffers( region, iReg )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag, t_plag_input, t_buffer_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global  
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetDimensPhys
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  
  INTEGER        :: iReg

! ... loop variables
  INTEGER :: iPatch, iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: bcType, iLev, nPatches, nPcls
  INTEGER :: iPclsRegionIn, iPclsRegionBuf, nPclsPrev
  INTEGER :: iCPlag, ipcbeg, ipcend, ibeg, iend, idir
  INTEGER :: jCPlag, jpcbeg, jpcend, jbeg, jend, jdir
  INTEGER :: kCPlag, kpcbeg, kpcend, kbeg, kend, kdir
  INTEGER :: nPclsBeg, nPclsEnd
  INTEGER :: errorFlag
  INTEGER :: iPatchCounterSum

  INTEGER, POINTER, DIMENSION(:,:) :: pAiv, pAivOld
  INTEGER, ALLOCATABLE, DIMENSION(:) ::patchCounter

  LOGICAL :: lboundSkip(6)

  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv, pArvOld, pCv, pCvOld, &
                                           pDv, pTv, pRhs, pRhsSum

  TYPE(t_patch),       POINTER :: pPatch
  TYPE(t_plag),        POINTER :: pPlag 
  TYPE(t_buffer_plag), POINTER :: pBuffPlag
  TYPE(t_global),      POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PLAG_PatchLoadDataBuffers.F90,v $ $Revision: 1.3 $'

  global => region%global
    
  CALL RegisterFunction( global, 'PLAG_PatchLoadDataBuffers',&
  'PLAG_PatchLoadDataBuffers.F90' )

! Get dimensions --------------------------------------------------------------

  iLev  = region%currLevel
  nPcls = region%levels(iLev)%plag%nPcls 

  IF (nPcls == 0) GOTO 999 ! exit if no particles

  nPatches = region%nPatches

! Set pointers ----------------------------------------------------------------
   
  pPlag   => region%levels(iLev)%plag   
  pAiv    => pPlag%aiv
  pArv    => pPlag%arv
  pCv     => pPlag%cv
  pDv     => pPlag%dv
  pTv     => pPlag%tv
  pRhs    => pPlag%rhs
  pRhsSum => pPlag%rhsSum

  pAivOld => pPlag%aivOld
  pArvOld => pPlag%arvOld
  pCvOld  => pPlag%cvOld
  
! Get grid dimensions ---------------------------------------------------------
  
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )

! Allocate patch buffer counter -----------------------------------------------

  ALLOCATE( patchCounter(nPatches), STAT=errorFlag )
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) &
    CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! Initialize counters for particles in inside and buffer regions --------------
  
  iPclsRegionIn  = 0
  iPclsRegionBuf = 0
  patchCounter(1:nPatches) = 0
  
! Loop over particles ---------------------------------------------------------

  DO iPcls = 1, nPcls

    iCPlag = pAiv(AIV_PLAG_INDEXI,iPcls)
    jCPlag = pAiv(AIV_PLAG_INDEXJ,iPcls)
    kCPlag = pAiv(AIV_PLAG_INDEXK,iPcls)

! - Check if particle cell is within region -----------------------------------

    IF ( ipcbeg <= iCPlag .AND. iCPlag <= ipcend .AND. &
         jpcbeg <= jCPlag .AND. jCPlag <= jpcend .AND. &
         kpcbeg <= kCPlag .AND. kCPlag <= kpcend ) THEN

! --- Shift particle datastructure only if particle is not in its proper spot -
 
      iPclsRegionIn = iPclsRegionIn + 1 

      IF ( iPclsRegionIn /= iPcls ) THEN
        pAiv(   :,iPclsRegionIn) = pAiv(   :,iPcls)
        pArv(   :,iPclsRegionIn) = pArv(   :,iPcls)
        pCv(    :,iPclsRegionIn) = pCv(    :,iPcls)
        pDv(    :,iPclsRegionIn) = pDv(    :,iPcls)
        pTv(    :,iPclsRegionIn) = pTv(    :,iPcls)
        pRhs(   :,iPclsRegionIn) = pRhs(   :,iPcls)
        pRhsSum(:,iPclsRegionIn) = pRhsSum(:,iPcls)

        pAivOld(:,iPclsRegionIn) = pAivOld(:,iPcls)
        pArvOld(:,iPclsRegionIn) = pArvOld(:,iPcls)
        pCvOld( :,iPclsRegionIn) = pCvOld( :,iPcls)
      ENDIF ! iPclsRegionIn
 
    ELSE 

! --- Set lboundSkip(:) to .TRUE. for values of lbound adjacent to a boundary -

      lboundSkip(1:6) = .TRUE.

      IF (     iCPlag < ipcbeg) THEN
        lboundSkip(1) = .FALSE.
      ELSE IF (iCPlag > ipcend) THEN
        lboundSkip(2) = .FALSE.
      ENDIF ! iCPlag

      IF (     jCPlag < jpcbeg) THEN
        lboundSkip(3) = .FALSE.
      ELSE IF (jCPlag > jpcend) THEN
        lboundSkip(4) = .FALSE.
      ENDIF ! jCPlag

      IF (     kCPlag < kpcbeg) THEN
        lboundSkip(5) = .FALSE.
      ELSE IF (kCPlag > kpcend) THEN
        lboundSkip(6) = .FALSE.
      ENDIF ! kCPlag

! --- Remove particle from active datastructure -------------------------------

! --- Loop over patches -------------------------------------------------------

      DO iPatch=1,nPatches

        pPatch => region%levels(iLev)%patches(iPatch)
        IF (lboundSkip(pPatch%lbound)) CYCLE

        bcType = pPatch%bcType

! ----- Select communication boundary conditions ------------------------------
        
        IF ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
            (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
            (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE)) THEN 

! ------- Check if particle cell is within (first dummy cell layer of) patch --

          CALL RFLO_GetPatchIndices(   region,pPatch,iLev, &
                                       ibeg,iend,jbeg,jend,kbeg,kend )
          CALL RFLO_GetPatchDirection( pPatch,idir,jdir,kdir )

          IF ( ibeg - idir <= iCPlag .AND. iCPlag <= iend - idir .AND. &
               jbeg - jdir <= jCPlag .AND. jCPlag <= jend - jdir .AND. &
               kbeg - kdir <= kCPlag .AND. kCPlag <= kend - kdir ) THEN

! --------- Load communication buffer arrays ----------------------------------

! --------- Every variable packed down above is sent --------------------------

            patchCounter(iPatch) = patchCounter(iPatch)+1
            iPclsRegionBuf       = patchCounter(iPatch)

            pBuffPlag => pPatch%bufferPlag 

            pBuffPlag%aiv(   :,iPclsRegionBuf) = pAiv(   :,iPcls)
            pBuffPlag%arv(   :,iPclsRegionBuf) = pArv(   :,iPcls)
            pBuffPlag%cv(    :,iPclsRegionBuf) = pCv(    :,iPcls)
            pBuffPlag%dv(    :,iPclsRegionBuf) = pDv(    :,iPcls)
            pBuffPlag%tv(    :,iPclsRegionBuf) = pTv(    :,iPcls)
            pBuffPlag%rhs(   :,iPclsRegionBuf) = pRhs(   :,iPcls)
            pBuffPlag%rhsSum(:,iPclsRegionBuf) = pRhsSum(:,iPcls)
  
            pBuffPlag%aivOld(:,iPclsRegionBuf) = pAivOld(:,iPcls)
            pBuffPlag%arvOld(:,iPclsRegionBuf) = pArvOld(:,iPcls)
            pBuffPlag%cvOld( :,iPclsRegionBuf) = pCvOld( :,iPcls)

#ifdef PLAG_DEBUG    
            IF( pBuffPlag%nBuffSize /= 0 )&
            WRITE(STDOUT,'(A,4I4)') &
            '      PLAG_PatchLoadDataBuffers-iReg: nBuffSize,iPatch,patchCounter', &
            iReg,pBuffPlag%nBuffSize,iPatch,patchCounter(iPatch)
#endif    

          ENDIF ! iCPlag, jCPlag, kCPlag
        ENDIF   ! bcType  
      ENDDO     ! iPatch
    ENDIF       ! iCPlag, jCPlag, kCPlag
  ENDDO         ! iPcls

! Get new particle datasize and buffer ----------------------------------------

  nPclsPrev   = pPlag%nPcls    
  pPlag%nPcls = iPclsRegionIn
  iPatchCounterSum = SUM(patchCounter)

#ifdef PLAG_DEBUG    
  IF(nPclsPrev /= pPlag%nPcls) &
  WRITE(STDOUT,'(A,I4,2I8,I4)') &
  '      PLAG_PatchLoadDataBuffers: iReg,nPclsPrev,nPclsCurr,iPatchCounterSum ',&
        iReg,nPclsPrev,pPlag%nPcls,iPatchCounterSum
#endif 

! Reinitialize reshuffled particle datastructure to account for ---------------
!   region with null size particle --------------------------------------------
!   perform if data reshuffled and particle size null -------------------------

  IF ( iPatchCounterSum /=0 .AND. pPlag%nPcls == 0) THEN 
    nPclsBeg = MAX(1,pPlag%nPcls+1)
    nPclsEnd = nPclsPrev
  
    pPlag%aiv(:   ,nPclsBeg:nPclsEnd) = 0
    pPlag%aivOld(:,nPclsBeg:nPclsEnd) = 0
    pPlag%arv(:   ,nPclsBeg:nPclsEnd) = 0.0_RFREAL
    pPlag%arvOld(:,nPclsBeg:nPclsEnd) = 0.0_RFREAL
    pPlag%cv(:    ,nPclsBeg:nPclsEnd) = 0.0_RFREAL
    pPlag%cvOld(: ,nPclsBeg:nPclsEnd) = 0.0_RFREAL
    pPlag%rhs(:   ,nPclsBeg:nPclsEnd) = 0.0_RFREAL
    pPlag%rhsSum(:,nPclsBeg:nPclsEnd) = 0.0_RFREAL

#ifdef PLAG_DEBUG    
    WRITE(STDOUT,'(A,I4,2I6)') &
    '      PLAG_PatchLoadDataBuffers: iReg,nPclsBeg,nPclsEnd = ',&
           iReg,nPclsBeg,nPclsEnd
#endif
  ENDIF ! iPclsRegionBuf

! Deallocate patch buffer counter ----------------------------------------------

  DEALLOCATE( patchCounter, STAT=errorFlag )
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) &
    CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

! finalize --------------------------------------------------------------------

999  CONTINUE
  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_PatchLoadDataBuffers

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_PatchLoadDataBuffers.F90,v $
! Revision 1.3  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:57  fnajjar
! Initial revision after changing case
!
! Revision 1.12  2004/04/09 23:14:57  fnajjar
! Added PatchCounter for multiple patches connecting to different regions
!
! Revision 1.11  2003/05/27 19:19:57  fnajjar
! Removed distPartBurning and all pertinent LOGICAL datastructure
!
! Revision 1.10  2003/05/06 23:56:40  fnajjar
! Bug fix to reinclude iPclsRegionIn
!
! Revision 1.9  2003/05/06 19:35:19  fnajjar
! Bug fix that removes nBuffSize=0 since it is computed in PLAG_patchGetBufferSize
!
! Revision 1.8  2003/05/01 22:58:30  jferry
! overhauled structure in order to optimize performance
!
! Revision 1.7  2003/04/17 00:00:27  fnajjar
! Commented out I/O to STDOUT
!
! Revision 1.6  2003/01/17 20:34:25  f-najjar
! Fixed FORMAT statements with iReg
!
! Revision 1.5  2003/01/17 19:34:08  f-najjar
! Added iReg in WRITE statements
!
! Revision 1.4  2003/01/16 15:39:51  f-najjar
! Included iReg in calling sequence since iRegionGlobal non-existent for Rocflo
!
! Revision 1.3  2003/01/16 15:34:24  f-najjar
! Move IF statement for nPcls check after extracting value
!
! Revision 1.2  2003/01/13 19:06:27  f-najjar
! Added aivOld, arvOld and cvOld to loading sequence
!
! Revision 1.1  2002/10/25 14:18:35  f-najjar
! Initial Import of Rocpart
!
!******************************************************************************







