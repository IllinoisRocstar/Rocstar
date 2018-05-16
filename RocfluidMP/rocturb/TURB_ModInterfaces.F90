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
! Purpose: set explicit interfaces of TURB subroutines and functions.
!
! Description: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_ModInterfaces.F90,v 1.28 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************
  
MODULE TURB_ModInterfaces

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! Interfaces to external code
! =============================================================================

#ifdef RFLO
  SUBROUTINE TURB_CalcStrainRate( region,ibn,ien,grIndx,gradi,gradj,gradk, &
                                  sRateI,sRateJ,sRateK )
#endif
#ifdef RFLU
  SUBROUTINE TURB_CalcStrainRate( region,ibn,ien,grIndx,gradf,sRateI )
#endif
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: ibn, ien
#ifdef RFLO
    INTEGER        :: grIndx(9)
    REAL(RFREAL), POINTER :: gradi(:,:),gradj(:,:),gradk(:,:)
    REAL(RFREAL), POINTER :: sRateI(:,:),sRateJ(:,:),sRateK(:,:)
#endif
#ifdef RFLU
    INTEGER        :: grIndx(3)
    REAL(RFREAL), POINTER :: gradf(:,:,:)
    REAL(RFREAL), POINTER :: sRateI(:,:)
#endif
  END SUBROUTINE TURB_CalcStrainRate

  SUBROUTINE TURB_CalcVortic( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_CalcVortic

  SUBROUTINE TURB_CheckParamInput( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_CheckParamInput

  SUBROUTINE TURB_DerivedInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_DerivedInputValues

  SUBROUTINE TURB_GetModelStressCell( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_GetModelStressCell

  SUBROUTINE TURB_GetTvCell( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_GetTvCell

  SUBROUTINE TURB_InitInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_InitInputValues

  SUBROUTINE TURB_LesCalcEddyVis( region,ibn,ien,ijk )
    USE ModDataStruct, ONLY : t_region
#ifdef RFLO
    TYPE(t_region), TARGET  :: region
#endif
#ifdef RFLU
    TYPE(t_region), POINTER :: region
#endif
    INTEGER                 :: ibn,ien,ijk
  END SUBROUTINE TURB_LesCalcEddyVis

  SUBROUTINE TURB_LesCoefDynMixd( region,ibn,ien,ijk )
    USE ModDataStruct, ONLY : t_region
#ifdef RFLO
    TYPE(t_region), TARGET  :: region
#endif
#ifdef RFLU
    TYPE(t_region), POINTER :: region
#endif
    INTEGER                 :: ibn,ien,ijk
  END SUBROUTINE TURB_LesCoefDynMixd

  SUBROUTINE TURB_LesCoefDynSmag( region,ibn,ien,ijk )
    USE ModDataStruct, ONLY : t_region
#ifdef RFLO
    TYPE(t_region), TARGET  :: region
#endif
#ifdef RFLU
    TYPE(t_region), POINTER :: region
#endif
    INTEGER                 :: ibn,ien,ijk
  END SUBROUTINE TURB_LesCoefDynSmag

  SUBROUTINE TURB_LesContract( region,ijk )
    USE ModDataStruct, ONLY : t_region
#ifdef RFLO
    TYPE(t_region)          :: region
#endif
#ifdef RFLU
    TYPE(t_region), POINTER :: region
#endif
    INTEGER                 :: ijk
  END SUBROUTINE TURB_LesContract

  SUBROUTINE TURB_LesEsgModel4( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)  :: region
  END SUBROUTINE TURB_LesEsgModel4

  SUBROUTINE TURB_LesFluxFixSmag( region,ibn,ien )
    USE ModDataStruct, ONLY : t_region
#ifdef RFLO
    TYPE(t_region)          :: region
#endif
#ifdef RFLU
    TYPE(t_region), POINTER :: region
#endif
    INTEGER                 :: ibn,ien
  END SUBROUTINE TURB_LesFluxFixSmag

  SUBROUTINE TURB_LesFluxScalSim( region,ibn,ien )
    USE ModDataStruct, ONLY : t_region
#ifdef RFLO
    TYPE(t_region), TARGET  :: region
#endif
#ifdef RFLU
    TYPE(t_region), POINTER :: region
#endif
    INTEGER                 :: ibn,ien
  END SUBROUTINE TURB_LesFluxScalSim

  SUBROUTINE TURB_LesGetEddyVis( region,ibc,iec,ibn,ien )
    USE ModDataStruct, ONLY : t_region
#ifdef RFLO
    TYPE(t_region), TARGET  :: region
#endif
#ifdef RFLU
    TYPE(t_region), POINTER :: region
#endif
    INTEGER        :: ibc,iec,ibn,ien
  END SUBROUTINE TURB_LesGetEddyVis

  SUBROUTINE TURB_LesHij( region,ijk )
    USE ModDataStruct, ONLY : t_region
#ifdef RFLO
    TYPE(t_region)          :: region
#endif
#ifdef RFLU
    TYPE(t_region), POINTER :: region
#endif
    INTEGER                 :: ijk
  END SUBROUTINE TURB_LesHij

  SUBROUTINE TURB_LesLij( region,ijk,nDel,lij )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE TURB_ModParameters
#ifdef RFLO
    TYPE(t_region)          :: region
#endif
#ifdef RFLU
    TYPE(t_region), POINTER :: region
#endif
    INTEGER                 :: ijk,nDel(DIRI:DIRK)
    REAL(RFREAL), POINTER   :: lij(:,:)
  END SUBROUTINE TURB_LesLij

  SUBROUTINE TURB_LesMij( region,ijk )
    USE ModDataStruct, ONLY : t_region
#ifdef RFLO
    TYPE(t_region)          :: region
#endif
#ifdef RFLU
    TYPE(t_region), POINTER :: region
#endif
    INTEGER                 :: ijk
  END SUBROUTINE TURB_LesMij

  SUBROUTINE TURB_LesRkInit( region,iStage )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER :: iStage
  END SUBROUTINE TURB_LesRkInit

  SUBROUTINE TURB_LesTestRhoV( region,ibc,iec )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: ibc,iec
  END SUBROUTINE TURB_LesTestRhoV

  SUBROUTINE TURB_RansEmsInit( region,iStage )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER :: iStage
  END SUBROUTINE TURB_RansEmsInit

  SUBROUTINE TURB_RansEmsUpdate( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_RansEmsUpdate

  SUBROUTINE TURB_RansRkInit( region, iStage )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER :: iStage
  END SUBROUTINE TURB_RansRkInit

  SUBROUTINE TURB_RansSAGetEddyVis( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_RansSAGetEddyVis

  SUBROUTINE TURB_RansSASourceTerms( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_RansSASourceTerms

  SUBROUTINE TURB_RansSAVisFlux( region )
    USE ModDataStruct, ONLY : t_region
#ifdef RFLO
    TYPE(t_region) :: region
#endif
#ifdef RFLU
    TYPE(t_region), POINTER :: region
#endif
  END SUBROUTINE TURB_RansSAVisFlux

  SUBROUTINE TURB_RansSAVisFluxPatch( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_RansSAVisFluxPatch

  SUBROUTINE TURB_RansTotalTv( region,indxMu,indxTCo,tvt )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: indxMu, indxTCo
    REAL(RFREAL), POINTER :: tvt(:,:)
  END SUBROUTINE TURB_RansTotalTv

  SUBROUTINE TURB_RansWallDistOVPatch( region,patch )  
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_RansWallDistOVPatch

  SUBROUTINE TURB_ReadBcInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_ReadBcInputFile

  SUBROUTINE TURB_ReadInputFile( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_ReadInputFile

  SUBROUTINE TURB_ReadTurbSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_ReadTurbSection

  SUBROUTINE TURB_StatCCollector( region,iBegSt,iEndSt,colVar )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: iBegSt, iEndSt
    REAL(RFREAL), POINTER :: colVar(:,:)
  END SUBROUTINE TURB_StatCCollector

#ifdef RFLO
  SUBROUTINE TURB_StatFCollector( region,ijk,iBegSt,iEndSt,colVar )
#endif
#ifdef RFLU
  SUBROUTINE TURB_StatFCollector( region,ijk,iBegSt,iEndSt,colVar,colBVar )
#endif
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: ijk, iBegSt, iEndSt
    REAL(RFREAL), POINTER :: colVar(:,:)
#ifdef RFLU
    REAL(RFREAL), POINTER :: colBVar(:,:)
#endif
  END SUBROUTINE TURB_StatFCollector

  SUBROUTINE TURB_VFluxHybrid( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), TARGET :: region
  END SUBROUTINE TURB_VFluxHybrid

  SUBROUTINE TURB_VFluxHybridPatch( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), TARGET :: region
    TYPE(t_patch)          :: patch
  END SUBROUTINE TURB_VFluxHybridPatch

  SUBROUTINE TURB_VisFluxEddy( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), TARGET :: region
  END SUBROUTINE TURB_VisFluxEddy

  SUBROUTINE TURB_VisFluxEddyPatch( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), TARGET :: region
    TYPE(t_patch)          :: patch
  END SUBROUTINE TURB_VisFluxEddyPatch

  SUBROUTINE TURB_WlmFluxPatch( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_WlmFluxPatch

  SUBROUTINE TURB_WlmInitia( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_WlmInitia

  SUBROUTINE TURB_WlmReyAnalogy( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_WlmReyAnalogy

  SUBROUTINE TURB_WlmTauWallMapping( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_WlmTauWallMapping

  SUBROUTINE TURB_WlmUpdate( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_WlmUpdate

  SUBROUTINE TURB_WlmUpdateBndlay( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_WlmUpdateBndlay

  SUBROUTINE TURB_CoRansWallDistOV( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_CoRansWallDistOV

  SUBROUTINE TURB_CoWlmReadBcSection( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE TURB_CoWlmReadBcSection

#ifdef RFLO
! =============================================================================
! Rocflo-specific routines
! =============================================================================

  SUBROUTINE TURB_FloExtrapIntCellScal( region,fVec )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)          :: region
    REAL(RFREAL), POINTER   :: fVec(:)
  END SUBROUTINE TURB_FloExtrapIntCellScal

  SUBROUTINE TURB_FloExtrapIntCellVec( region,idBeg,idEnd,fVec )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)          :: region
    INTEGER                 :: idBeg,idEnd
    REAL(RFREAL), POINTER   :: fVec(:,:)
  END SUBROUTINE TURB_FloExtrapIntCellVec

  SUBROUTINE TURB_FloExtrapolCellVec( region,idBeg,idEnd,fVec )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)          :: region
    INTEGER                 :: idBeg,idEnd
    REAL(RFREAL), POINTER   :: fVec(:,:)
  END SUBROUTINE TURB_FloExtrapolCellVec

  SUBROUTINE TURB_FloExtrapolFaceVec( region,intDIR,idBeg,idEnd,fVec )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)          :: region
    INTEGER                 :: intDIR,idBeg,idEnd
    REAL(RFREAL), POINTER   :: fVec(:,:)
  END SUBROUTINE TURB_FloExtrapolFaceVec

  SUBROUTINE TURB_FloFaceVolume( region,ijk )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)  :: region
    INTEGER         :: ijk
  END SUBROUTINE TURB_FloFaceVolume

  SUBROUTINE TURB_FloFaceWidth( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_FloFaceWidth

  SUBROUTINE TURB_FloFaceWidthDummy( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_FloFaceWidthDummy

  SUBROUTINE TURB_FloFaceWidthDummyConn( region,lbound,idir,jdir,kdir, &
                                  indBeg,indEnd,jndBeg,jndEnd,kndBeg,kndEnd )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: lbound,idir,jdir,kdir
    INTEGER        :: indBeg, indEnd, jndBeg, jndEnd, kndBeg, kndEnd
  END SUBROUTINE TURB_FloFaceWidthDummyConn

  SUBROUTINE TURB_FloFaceWidthDummyPhys( region,lbound,idir,jdir,kdir, &
                                  indBeg,indEnd,jndBeg,jndEnd,kndBeg,kndEnd )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: lbound,idir,jdir,kdir
    INTEGER        :: indBeg, indEnd, jndBeg, jndEnd, kndBeg, kndEnd
  END SUBROUTINE TURB_FloFaceWidthDummyPhys

  SUBROUTINE TURB_FloLesAverageFace( region,ijk,faceVar,avgFaceVar )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)        :: region
    INTEGER               :: ijk
    REAL(RFREAL), POINTER :: faceVar(:),avgFaceVar(:)
  END SUBROUTINE TURB_FloLesAverageFace

  SUBROUTINE TURB_FloLesGenC2F( region,ijk )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)  :: region
    INTEGER         :: ijk
  END SUBROUTINE TURB_FloLesGenC2F

  SUBROUTINE TURB_FloLesUniFiltCC( region,nDel,idBeg,idEnd,fVar,fbVar )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE TURB_ModParameters
    TYPE(t_region)        :: region
    INTEGER               :: ijk,nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL), POINTER :: fVar(:,:),fbVar(:,:)
  END SUBROUTINE TURB_FloLesUniFiltCC

  SUBROUTINE TURB_FloLesUniFiltCCI( global,nDumi,ibeg,iend,jbeg,jend,kbeg,kend, &
                         iCOff,ijCOff,nDel,idBeg,idEnd,fact1,fact2,fVar, &
                         filtVar )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE TURB_ModParameters
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
    INTEGER               :: nDumi,ibeg,iend,jbeg,jend,kbeg,kend,iCOff,ijCOff
    INTEGER               :: nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL)          :: fact1(FILWIDTH_FOUR),fact2(FILWIDTH_FOUR)
    REAL(RFREAL), POINTER :: fVar(:,:),filtVar(:,:)
  END SUBROUTINE TURB_FloLesUniFiltCCI

  SUBROUTINE TURB_FloLesUniFiltCCJ( global,nDumi,ibeg,iend,jbeg,jend,kbeg,kend, &
                         iCOff,ijCOff,nDel,idBeg,idEnd,fact1,fact2,fVar, &
                         filtVar )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE TURB_ModParameters
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
    INTEGER               :: nDumi,ibeg,iend,jbeg,jend,kbeg,kend,iCOff,ijCOff
    INTEGER               :: nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL)          :: fact1(FILWIDTH_FOUR),fact2(FILWIDTH_FOUR)
    REAL(RFREAL), POINTER :: fVar(:,:),filtVar(:,:)
  END SUBROUTINE TURB_FloLesUniFiltCCJ

  SUBROUTINE TURB_FloLesUniFiltCCK( global,nDumi,ibeg,iend,jbeg,jend,kbeg,kend, &
                         iCOff,ijCOff,nDel,idBeg,idEnd,fact1,fact2,fVar, &
                         filtVar )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE TURB_ModParameters
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
    INTEGER               :: nDumi,ibeg,iend,jbeg,jend,kbeg,kend,iCOff,ijCOff
    INTEGER               :: nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL)          :: fact1(FILWIDTH_FOUR),fact2(FILWIDTH_FOUR)
    REAL(RFREAL), POINTER :: fVar(:,:),filtVar(:,:)
  END SUBROUTINE TURB_FloLesUniFiltCCK

  SUBROUTINE TURB_FloLesUniFiltFF( region,ijk,nDel,idBeg,idEnd,fVar,fbVar )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE TURB_ModParameters
    TYPE(t_region)        :: region
    INTEGER               :: ijk,nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL), POINTER :: fVar(:,:),fbVar(:,:)
  END SUBROUTINE TURB_FloLesUniFiltFF

  SUBROUTINE TURB_FloLesUniFiltFFI( global,ibeg,iend,jbeg,jend,kbeg,kend,iNOff, &
                          ijNOff,nDel,idBeg,idEnd,fact1,fact2,fVar,filtVar )
    USE ModDataTypes
    USE TURB_ModParameters
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
    INTEGER               :: ibeg,iend,jbeg,jend,kbeg,kend,iNOff,ijNOff
    INTEGER               :: nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL)          :: fact1(FILWIDTH_FOUR),fact2(FILWIDTH_FOUR)
    REAL(RFREAL), POINTER :: fVar(:,:),filtVar(:,:)
  END SUBROUTINE TURB_FloLesUniFiltFFI

  SUBROUTINE TURB_FloLesUniFiltFFJ( global,ibeg,iend,jbeg,jend,kbeg,kend,iNOff, &
                          ijNOff,nDel,idBeg,idEnd,fact1,fact2,fVar,filtVar )
    USE ModDataTypes
    USE TURB_ModParameters
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
    INTEGER               :: ibeg,iend,jbeg,jend,kbeg,kend,iNOff,ijNOff
    INTEGER               :: nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL)          :: fact1(FILWIDTH_FOUR),fact2(FILWIDTH_FOUR)
    REAL(RFREAL), POINTER :: fVar(:,:),filtVar(:,:)
  END SUBROUTINE TURB_FloLesUniFiltFFJ

  SUBROUTINE TURB_FloLesUniFiltFFK( global,ibeg,iend,jbeg,jend,kbeg,kend,iNOff, &
                        ijNOff,nDel,idBeg,idEnd,fact1,fact2,fVar,filtVar )
    USE ModDataTypes
    USE TURB_ModParameters
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
    INTEGER               :: ibeg,iend,jbeg,jend,kbeg,kend,iNOff,ijNOff
    INTEGER               :: nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL)          :: fact1(FILWIDTH_FOUR),fact2(FILWIDTH_FOUR)
    REAL(RFREAL), POINTER :: fVar(:,:),filtVar(:,:)
  END SUBROUTINE TURB_FloLesUniFiltFFK

  SUBROUTINE TURB_FloLesGenFiltCC( region,nDel,idBeg,idEnd,fVar,fbVar )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE TURB_ModParameters
    TYPE(t_region)        :: region
    INTEGER               :: nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL), POINTER :: fVar(:,:),fbVar(:,:)
  END SUBROUTINE TURB_FloLesGenFiltCC

  SUBROUTINE TURB_FloLesGenFiltCCI( region,nDumi,ibeg,iend,jbeg,jend,kbeg,kend, &
                               iCOff,ijCOff,nDel,idBeg,idEnd,fVar,filtVar )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE TURB_ModParameters
    TYPE(t_region), TARGET :: region
    INTEGER                :: nDumi,ibeg,iend,jbeg,jend,kbeg,kend,iCOff,ijCOff
    INTEGER                :: nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL), POINTER  :: fVar(:,:),filtVar(:,:)
  END SUBROUTINE TURB_FloLesGenFiltCCI

  SUBROUTINE TURB_FloLesGenFiltCCJ( region,nDumi,ibeg,iend,jbeg,jend,kbeg,kend, &
                               iCOff,ijCOff,nDel,idBeg,idEnd,fVar,filtVar )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE TURB_ModParameters
    TYPE(t_region), TARGET :: region
    INTEGER                :: nDumi,ibeg,iend,jbeg,jend,kbeg,kend,iCOff,ijCOff
    INTEGER                :: nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL), POINTER  :: fVar(:,:),filtVar(:,:)
  END SUBROUTINE TURB_FloLesGenFiltCCJ

  SUBROUTINE TURB_FloLesGenFiltCCK( region,nDumi,ibeg,iend,jbeg,jend,kbeg,kend, &
                               iCOff,ijCOff,nDel,idBeg,idEnd,fVar,filtVar )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE TURB_ModParameters
    TYPE(t_region), TARGET :: region
    INTEGER                :: nDumi,ibeg,iend,jbeg,jend,kbeg,kend,iCOff,ijCOff
    INTEGER                :: nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL), POINTER  :: fVar(:,:),filtVar(:,:)
  END SUBROUTINE TURB_FloLesGenFiltCCK

  SUBROUTINE TURB_FloLesGenFiltFF( region,ijk,nDel,idBeg,idEnd,fVar,fbVar )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE TURB_ModParameters
    TYPE(t_region)        :: region
    INTEGER               :: ijk,nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL), POINTER :: fVar(:,:),fbVar(:,:)
  END SUBROUTINE TURB_FloLesGenFiltFF

  SUBROUTINE TURB_FloLesGenFiltFFI( region,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
                                    iNOff,ijNOff,nDel,idBeg,idEnd,fVar,filtVar )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE TURB_ModParameters
    TYPE(t_region), TARGET :: region
    INTEGER                :: ijk,ibeg,iend,jbeg,jend,kbeg,kend,iNOff,ijNOff
    INTEGER                :: nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL), POINTER  :: fVar(:,:),filtVar(:,:)
  END SUBROUTINE TURB_FloLesGenFiltFFI

  SUBROUTINE TURB_FloLesGenFiltFFJ( region,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
                                    iNOff,ijNOff,nDel,idBeg,idEnd,fVar,filtVar )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE TURB_ModParameters
    TYPE(t_region), TARGET :: region
    INTEGER                :: ijk,ibeg,iend,jbeg,jend,kbeg,kend,iNOff,ijNOff
    INTEGER                :: nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL), POINTER  :: fVar(:,:),filtVar(:,:)
  END SUBROUTINE TURB_FloLesGenFiltFFJ

  SUBROUTINE TURB_FloLesGenFiltFFK( region,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
                                    iNOff,ijNOff,nDel,idBeg,idEnd,fVar,filtVar )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE TURB_ModParameters
    TYPE(t_region), TARGET :: region
    INTEGER                :: ijk,ibeg,iend,jbeg,jend,kbeg,kend,iNOff,ijNOff
    INTEGER                :: nDel(DIRI:DIRK),idBeg,idEnd
    REAL(RFREAL), POINTER  :: fVar(:,:),filtVar(:,:)
  END SUBROUTINE TURB_FloLesGenFiltFFK

  SUBROUTINE TURB_FloLesGenCoCC( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), TARGET :: region
  END SUBROUTINE TURB_FloLesGenCoCC

  SUBROUTINE TURB_FloLesGenCoCCHi( global,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
                                minIdx,maxIdx,segId,iCOff,ijCOff,iNOff,ijNOff, &
                                ds,segm,ccCofA,ccCofB )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
    INTEGER               :: ijk,ibeg,iend,jbeg,jend,kbeg,kend,minIdx,maxIdx
    INTEGER               :: segId,iCOff,ijCOff,iNOff,ijNOff
    REAL(RFREAL)          :: ds(minIdx:maxIdx)
    REAL(RFREAL), POINTER :: segm(:,:),ccCofA(:,:),ccCofB(:,:)
  END SUBROUTINE TURB_FloLesGenCoCCHi

  SUBROUTINE TURB_FloLesGenCoCCLo( global,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
                                minIdx,maxIdx,segId,iCOff,ijCOff,iNOff,ijNOff, &
                                ds,segm,ccCofA,ccCofB )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
    INTEGER               :: ijk,ibeg,iend,jbeg,jend,kbeg,kend,minIdx,maxIdx
    INTEGER               :: segId,iCOff,ijCOff,iNOff,ijNOff
    REAL(RFREAL)          :: ds(minIdx:maxIdx)
    REAL(RFREAL), POINTER :: segm(:,:),ccCofA(:,:),ccCofB(:,:)
  END SUBROUTINE TURB_FloLesGenCoCCLo

  SUBROUTINE TURB_FloLesGenCoFF( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), TARGET :: region
  END SUBROUTINE TURB_FloLesGenCoFF

  SUBROUTINE TURB_FloLesGenCoFFHi( global,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
                                minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm, &
                                ffCofA,ffCofB )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
    INTEGER               :: ijk,ibeg,iend,jbeg,jend,kbeg,kend,minIdx,maxIdx
    INTEGER               :: segId,iNOff,ijNOff
    REAL(RFREAL)          :: ds(minIdx:maxIdx)
    REAL(RFREAL), POINTER :: segm(:,:),ffCofA(:,:),ffCofB(:,:)
  END SUBROUTINE TURB_FloLesGenCoFFHi

  SUBROUTINE TURB_FloLesGenCoFFLo( global,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
                                minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm, &
                                ffCofA,ffCofB )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
    INTEGER               :: ijk,ibeg,iend,jbeg,jend,kbeg,kend,minIdx,maxIdx
    INTEGER               :: segId,iNOff,ijNOff
    REAL(RFREAL)          :: ds(minIdx:maxIdx)
    REAL(RFREAL), POINTER :: segm(:,:),ffCofA(:,:),ffCofB(:,:)
  END SUBROUTINE TURB_FloLesGenCoFFLo

  SUBROUTINE TURB_FloLesGenCoFCHi( global,filtDir,ibeg,iend,jbeg,jend,kbeg, &
                                kend,minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm, &
                                ffCofA,ffCofB )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
    INTEGER               :: filtDir,ibeg,iend,jbeg,jend,kbeg,kend
    INTEGER               :: minIdx,maxIdx,segId,iNOff,ijNOff
    REAL(RFREAL)          :: ds(minIdx:maxIdx)
    REAL(RFREAL), POINTER :: segm(:,:),ffCofA(:,:),ffCofB(:,:)
  END SUBROUTINE TURB_FloLesGenCoFCHi

  SUBROUTINE TURB_FloLesGenCoFCLo( global,filtDir,ibeg,iend,jbeg,jend,kbeg, &
                                kend,minIdx,maxIdx,segId,iNOff,ijNOff,ds,segm, &
                                ffCofA,ffCofB )
    USE ModDataTypes
    USE ModGlobal, ONLY : t_global
    TYPE(t_global), POINTER :: global
    INTEGER               :: filtDir,ibeg,iend,jbeg,jend,kbeg,kend
    INTEGER               :: minIdx,maxIdx,segId,iNOff,ijNOff
    REAL(RFREAL)          :: ds(minIdx:maxIdx)
    REAL(RFREAL), POINTER :: segm(:,:),ffCofA(:,:),ffCofB(:,:)
  END SUBROUTINE TURB_FloLesGenCoFCLo

  SUBROUTINE TURB_FloRansBcondInflow( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_FloRansBcondInflow

  SUBROUTINE TURB_FloRansBcondInjection( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_FloRansBcondInjection

  SUBROUTINE TURB_FloRansBcondNoslipWall( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_FloRansBcondNoslipWall

  SUBROUTINE TURB_FloRansBcondSymmetry( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_FloRansBcondSymmetry 

  SUBROUTINE TURB_FloRansBcondZeroGrad( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_FloRansBcondZeroGrad

  SUBROUTINE TURB_FloRansCentralDissipation( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_FloRansCentralDissipation

  SUBROUTINE TURB_FloRansCorrCornEdgeCells( region,patch,bcType )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
    INTEGER        :: bcType
  END SUBROUTINE TURB_FloRansCorrCornEdgeCells

  SUBROUTINE TURB_FloRansExchCornEdgeCells( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE TURB_FloRansExchCornEdgeCells

  SUBROUTINE TURB_FloRansExchangeDummyConf( region,regionSrc,patch,patchSrc )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch, patchSrc
  END SUBROUTINE TURB_FloRansExchangeDummyConf

  SUBROUTINE TURB_FloRansRecvCornEdgeCells( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE TURB_FloRansRecvCornEdgeCells

  SUBROUTINE TURB_FloRansRecvDummyVals( region,regionSrc,patch,patchSrc )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch, patchSrc
  END SUBROUTINE TURB_FloRansRecvDummyVals

  SUBROUTINE TURB_FloRansSACentFlux( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_FloRansSACentFlux

  SUBROUTINE TURB_FloRansSACentFluxPatch( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_FloRansSACentFluxPatch

  SUBROUTINE TURB_FloRansSARoe1stFlux( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_FloRansSARoe1stFlux

  SUBROUTINE TURB_FloRansSARoe2ndFlux( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_FloRansSARoe2ndFlux

  SUBROUTINE TURB_FloRansSARoeFluxPatch( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_FloRansSARoeFluxPatch

  SUBROUTINE TURB_FloRansSendCornEdgeCells( regions,iReg )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
    INTEGER :: iReg
  END SUBROUTINE TURB_FloRansSendCornEdgeCells

  SUBROUTINE TURB_FloRansSendDummyConf( region,regionSrc,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_FloRansSendDummyConf

  SUBROUTINE TURB_FloRansSetCornEdgeCells( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE TURB_FloRansSetCornEdgeCells

  SUBROUTINE TURB_FloWlmMetric( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_FloWlmMetric

  SUBROUTINE TURB_FloWlmUpdateLoglay( region,patch )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE TURB_FloWlmUpdateLoglay
#endif

#ifdef RFLU
! =============================================================================
! Rocflu-specific routines
! =============================================================================

  SUBROUTINE TURB_FluCv2Cons( region,cvStateFuture )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), TARGET :: region
    INTEGER, INTENT(IN) :: cvStateFuture
  END SUBROUTINE TURB_FluCv2Cons

  SUBROUTINE TURB_FluCv2Prim( region,cvStateFuture )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), TARGET :: region
    INTEGER, INTENT(IN) :: cvStateFuture
  END SUBROUTINE TURB_FluCv2Prim

  SUBROUTINE TURB_FluFaceVolume( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)  :: region
  END SUBROUTINE TURB_FluFaceVolume

  SUBROUTINE TURB_FluLesBLij( region,nDel,lij )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE TURB_ModParameters
    TYPE(t_region), POINTER :: region
    INTEGER                 :: nDel(DIRI:DIRK)
    REAL(RFREAL), POINTER   :: lij(:,:)
  END SUBROUTINE TURB_FluLesBLij

  SUBROUTINE TURB_FluLesBMij( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: region
  END SUBROUTINE TURB_FluLesBMij

  SUBROUTINE TURB_FluLesC2F( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)  :: region
  END SUBROUTINE TURB_FluLesC2F
#endif

  END INTERFACE

END MODULE TURB_ModInterfaces

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_ModInterfaces.F90,v $
! Revision 1.28  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.27  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.26  2004/11/17 23:45:21  wasistho
! used generic RK-update for rocturb
!
! Revision 1.25  2004/10/22 23:17:03  wasistho
! added TURB_StatCCollector
!
! Revision 1.24  2004/08/04 02:45:37  wasistho
! removed turb%avgCoI,J,K as it is defined as grid%c2fCoI,J,K
!
! Revision 1.23  2004/07/30 22:35:25  wasistho
! replaced floLesUniC2F by floLesGenC2F
!
! Revision 1.22  2004/05/28 01:57:12  wasistho
! update unstructured grid LES
!
! Revision 1.21  2004/03/29 21:09:30  wasistho
! add flu routines
!
! Revision 1.20  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.19  2004/03/23 03:35:00  wasistho
! prepared for RFLU
!
! Revision 1.18  2004/03/19 22:33:47  wasistho
! confined flo interfaces within ifdef RFLO
!
! Revision 1.17  2004/03/13 03:15:03  wasistho
! get rid of flo/flu identifier in TURB_Co.. routines
!
! Revision 1.16  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.15  2004/03/08 23:29:09  wasistho
! changed turb nomenclature
!
! Revision 1.14  2004/02/26 21:20:26  wasistho
! added TURB_ransRkInit, TURB_lesRkInit, TURB_ransEmsInit
!
! Revision 1.13  2004/02/12 03:46:14  wasistho
! filled in RaNS lengthscale in dummy cells
!
! Revision 1.12  2004/01/23 00:35:43  wasistho
! added new RaNS edge/corners routines
!
! Revision 1.11  2003/10/27 04:52:26  wasistho
! added RaNS upwind schemes
!
! Revision 1.10  2003/10/16 20:19:00  wasistho
! installed RaNS in steady state flow (Exp.Mult.Stg)
!
! Revision 1.9  2003/10/09 23:06:53  wasistho
! renamed CalcEddyVis to LesCalcEddyVis
!
! Revision 1.8  2003/10/07 02:04:15  wasistho
! initial installation of RaNS-SA and DES
!
! Revision 1.7  2003/08/06 15:57:28  wasistho
! added CalcVortic and SolutionUpdate for vorticities
!
! Revision 1.6  2003/06/05 19:18:37  wasistho
! implemented heat transfer model
!
! Revision 1.5  2003/05/31 01:45:32  wasistho
! installed turb. wall layer model
!
! Revision 1.4  2003/05/24 02:09:48  wasistho
! turbulence statistics expanded
!
! Revision 1.3  2003/05/16 05:45:41  wasistho
! modified array range of CC-filtered
!
! Revision 1.2  2002/10/14 23:53:46  wasistho
! Install Rocturb
!
!
!******************************************************************************






