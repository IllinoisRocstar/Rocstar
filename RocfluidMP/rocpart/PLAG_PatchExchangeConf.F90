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
! Purpose: exchange values between buffers of the corresponding 
!          patches of the adjacent regions (both regions
!          are on the same processor).
!
! Description: none.
!
! Input: region    = current (source) region
!        regionDes = destination region
!        patch     = current patch of region
!        patchDes  = destination patch of regionDes
!        iReg      = index of current region
!        iRegDes   = index of destination region.
!
! Output: regions%levels%plagBuff = variables in buffers.
!
! Notes: intended for conforming grid boundaries only.
!        Patch of current (source) region has data residing in dummy cells and 
!        send them to patch of destination region.
!        RFLO srcLbound corresponds to the destination (Des) region for PLAG
!        since the communication pattern is opposite for PLAG.
!
!******************************************************************************
!
! $Id: PLAG_PatchExchangeConf.F90,v 1.3 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_PatchExchangeConf( region,regionDes,patch,patchDes,iReg,iRegDes )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_buffer_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModIndexing, ONLY   : GetIJK, IndIJKMap
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
        RFLO_GetCellOffset, RFLO_GetPatchMapping
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region, regionDes
  TYPE(t_patch)  :: patch, patchDes

  INTEGER :: iReg, iRegDes

! ... loop variables
  INTEGER :: idum, i, j, k, ii, jj, kk

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir, &
             iCOff, ijCOff, ijkD, iLev
  INTEGER :: ibegDes, iendDes, jbegDes, jendDes, kbegDes, kendDes, &
             idirDes, jdirDes, kdirDes, iCOffDes, ijCOffDes,       &
             ijkCDes, iLevDes
  INTEGER :: lb, lbDes, l1DesDir, l2DesDir, mapMat(3,4),           &
             nDumCells, nDumCellsDes

  INTEGER :: iBuff, iBuffDes, iBuffSrc, jBuffDes, jBuffSrc, &
             kBuffDes, kBuffSrc, ijkBuffSrc, nBuffSizeDes, nBuffSizeSrc
             
  LOGICAL :: align

  INTEGER, POINTER, DIMENSION(:,:) :: pAiv, pAivDes, pAivOld, pAivOldDes

  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv, pArvDes,              &
                                           pArvOld, pArvOldDes,        &
                                           pCv, pCvDes,                &
                                           pCvOld, pCvOldDes,          &
                                           pDv, pDvDes, pRhs, pRhsDes, &
                                           pRhsSum, pRhsSumDes, pTv, pTvDes

  TYPE(t_global),      POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_PatchExchangeConf.F90,v $ $Revision: 1.3 $'

  global => region%global
  
  CALL RegisterFunction( global,'PLAG_PatchExchangeConf',&
  'PLAG_PatchExchangeConf.F90' )

! check if the source region is active ----------------------------------------

  IF (regionDes%active == OFF) THEN
    CALL ErrorStop( global,ERR_SRCREGION_OFF,__LINE__ )
  ENDIF

! Get dimensions --------------------------------------------------------------

  iLev      = region%currLevel
  nDumCells = region%nDumCells

  iLevDes = regionDes%currLevel
  nDumCellsDes = regionDes%nDumCells
  
  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend, &
                             jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchIndices( regionDes,patchDes,iLevDes,ibegDes,iendDes, &
                             jbegDes,jendDes,kbegDes,kendDes )
  CALL RFLO_GetPatchDirection( patch   ,idir   ,jdir,   kdir    )
  CALL RFLO_GetPatchDirection( patchDes,idirDes,jdirDes,kdirDes )
  CALL RFLO_GetCellOffset( region   ,iLev,iCOff   ,ijCOff    )
  CALL RFLO_GetCellOffset( regionDes,iLevDes,iCOffDes,ijCOffDes )

! Set pointers ----------------------------------------------------------------
 
  pAiv  => patch%bufferPlag%aiv
  pArv  => patch%bufferPlag%arv
  pCv   => patch%bufferPlag%cv
  pDv   => patch%bufferPlag%dv
  pTv   => patch%bufferPlag%tv
  pRhs  => patch%bufferPlag%rhs
  pRhsSum  => patch%bufferPlag%rhsSum

  pAivOld  => patch%bufferPlag%aivOld
  pArvOld  => patch%bufferPlag%arvOld
  pCvOld   => patch%bufferPlag%cvOld

  pAivDes  => patchDes%bufferPlag%aiv
  pArvDes  => patchDes%bufferPlag%arv
  pCvDes   => patchDes%bufferPlag%cv
  pDvDes   => patchDes%bufferPlag%dv
  pTvDes   => patchDes%bufferPlag%tv
  pRhsDes  => patchDes%bufferPlag%rhs
  pRhsSumDes  => patchDes%bufferPlag%rhsSum

  pAivOldDes  => patchDes%bufferPlag%aivOld
  pArvOldDes  => patchDes%bufferPlag%arvOld
  pCvOldDes   => patchDes%bufferPlag%cvOld

  nBuffSizeSrc = patch%bufferPlag%nBuffSize

  nBuffSizeDes = nBuffSizeSrc
  patchDes%bufferPlag%nBuffSizeDes = nBuffSizeDes
  
! exit when Buffer Size is Zero 

  IF ( nBuffSizeSrc == 0 )  GOTO 999

! mapping between patches -----------------------------------------------------

                                       l1DesDir =  1
  IF (patch%srcL1beg > patch%srcL1end) l1DesDir = -1
                                       l2DesDir =  1
  IF (patch%srcL2beg > patch%srcL2end) l2DesDir = -1

  lb    = patch%lbound
  lbDes = patch%srcLbound
  align = patch%align

  CALL RFLO_GetPatchMapping( lb,lbDes,l1DesDir,l2DesDir,align,                &
                             idir,jdir,kdir,idirDes,jdirDes,kdirDes,          &
                             ibeg,iend,jbeg,jend,kbeg,kend,                   &
                             ibegDes,iendDes,jbegDes,jendDes,kbegDes,kendDes, &
                             mapMat )

! Loop over buffer particles of current patch ---------------------------------

  DO iBuff=1, nBuffSizeSrc

    iBuffSrc   = pAiv(AIV_PLAG_INDEXI,iBuff)
    jBuffSrc   = pAiv(AIV_PLAG_INDEXJ,iBuff)
    kBuffSrc   = pAiv(AIV_PLAG_INDEXK,iBuff)
    ijkBuffSrc = pAiv(AIV_PLAG_ICELLS,iBuff)
    
    ijkCDes = IndIJKMap(iBuffSrc,jBuffSrc,kbuffSrc,mapMat,iCOffDes,ijCOffDes)
    CALL GetIJK(ijkCDes,iCOffDes,ijCOffDes,nDumCellsDes,iBuffDes,jBuffDes,kBuffDes)

! - Update aiv field ----------------------------------------------------------
    
    pAivDes(AIV_PLAG_ICELLS,iBuff) = ijkCDes    
    pAivDes(AIV_PLAG_INDEXI,iBuff) = iBuffDes
    pAivDes(AIV_PLAG_INDEXJ,iBuff) = jBuffDes
    pAivDes(AIV_PLAG_INDEXK,iBuff) = kBuffDes

    pAivDes(AIV_PLAG_PIDINI,iBuff) = pAiv(AIV_PLAG_PIDINI,iBuff)
    pAivDes(AIV_PLAG_REGINI,iBuff) = pAiv(AIV_PLAG_REGINI,iBuff)
    pAivDes(AIV_PLAG_REGCRT,iBuff) = iRegDes
    pAivDes(AIV_PLAG_BURNSTAT,iBuff) = pAiv(AIV_PLAG_BURNSTAT,iBuff)
    
    pAivOldDes(:,iBuff) = pAivDes(:,iBuff)

! - Update arv, cv, dv, rhs, rhsSum, tv field ---------------------------------

    pArvDes(:,iBuff)    = pArv(:,iBuff)
    pCvDes(:,iBuff)     = pCv(:,iBuff)
    pDvDes(:,iBuff)     = pDv(:,iBuff)
    pRhsDes(:,iBuff)    = pRhs(:,iBuff)
    pRhsSumDes(:,iBuff) = pRhsSum(:,iBuff)
    pTvDes(:,iBuff)     = pTv(:,iBuff)

    pArvOldDes(:,iBuff)    = pArvOld(:,iBuff)
    pCvOldDes(:,iBuff)     = pCvOld(:,iBuff)

  ENDDO ! iBuff

! reset Buffer size on Source Region to Zero to trap null-particle regions ----

  patch%bufferPlag%nBuffSize = 0

! finalize

999 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_PatchExchangeConf

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_PatchExchangeConf.F90,v $
! Revision 1.3  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:55  fnajjar
! Initial revision after changing case
!
! Revision 1.6  2004/02/13 23:22:07  fnajjar
! Included new cv and aiv definitions for particle burning module
!
! Revision 1.5  2003/09/21 00:19:58  fnajjar
! Reset buffer size on source region to zero
!
! Revision 1.4  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.3  2003/02/04 22:32:34  f-najjar
! Bug fix in GetIJK call by loading nDumCellsDes
!
! Revision 1.2  2003/01/17 19:45:21  f-najjar
! Compute nBuffSizeDes for buffer size of destination buffers
!
! Revision 1.1  2003/01/13 19:14:50  f-najjar
! Initial Import of Multiblock capability
!
!******************************************************************************







