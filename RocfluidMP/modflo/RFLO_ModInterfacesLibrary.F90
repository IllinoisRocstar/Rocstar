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
! Purpose: set explicit interfaces to subroutines and functions.
!
! Description: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ModInterfacesLibrary.F90,v 1.23 2008/12/06 08:44:16 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************
  
MODULE RFLO_ModInterfacesLibrary

  IMPLICIT NONE

  INTERFACE

! =============================================================================
! ROCFLO specific library
! =============================================================================

  SUBROUTINE RFLO_ArcLengthBounds( region,xyz,arcLen12,arcLen34,arcLen56 )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    REAL(RFREAL), POINTER :: arcLen12(:,:), arcLen34(:,:), arcLen56(:,:)
    REAL(RFREAL), POINTER :: xyz(:,:)
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_ArcLengthBounds

  SUBROUTINE RFLO_BoundaryDeformation( region,boundMoved,edgeMoved, &
                                       arcLen12,arcLen34,arcLen56,  &
                                       xyzOld,dNode )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    LOGICAL :: boundMoved(6), edgeMoved(12)
    REAL(RFREAL), POINTER :: arcLen12(:,:), arcLen34(:,:), arcLen56(:,:)
    REAL(RFREAL), POINTER :: dNode(:,:), xyzOld(:,:)
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_BoundaryDeformation

  SUBROUTINE RFLO_CalcCellCentroids( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CalcCellCentroids

  SUBROUTINE RFLO_CalcFaceCentroids( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_CalcFaceCentroids

  SUBROUTINE RFLO_CalcGradVector( region,iBegV,iEndV,iBegG,iEndG, &
                                  var,gradi,gradj,gradk )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)        :: region
    INTEGER               :: iBegV, iEndV, iBegG, iEndG
    REAL(RFREAL), POINTER :: var(:,:), gradi(:,:), gradj(:,:), gradk(:,:)
  END SUBROUTINE RFLO_CalcGradVector

  SUBROUTINE RFLO_CalcGradFaces( region,ilev,iBegV,iEndV,iBegG,iEndG, &
                                 var,gradi,gradj,gradk )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)        :: region
    INTEGER               :: ilev, iBegV, iEndV, iBegG, iEndG
    REAL(RFREAL), POINTER :: var(:,:), gradi(:,:), gradj(:,:), gradk(:,:)
  END SUBROUTINE RFLO_CalcGradFaces

  SUBROUTINE RFLO_CalcGradConnBc( region,patch,iConBc,iBegV,iEndV,iBegG,iEndG, &
                                  var,gradi,gradj,gradk )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch, ONLY   : t_patch
    TYPE(t_region)        :: region
    TYPE(t_patch)         :: patch
    INTEGER               :: iConBc(6), iBegV, iEndV, iBegG, iEndG
    REAL(RFREAL), POINTER :: var(:,:), gradi(:,:), gradj(:,:), gradk(:,:)
  END SUBROUTINE RFLO_CalcGradConnBc

  SUBROUTINE RFLO_CalcGradPhysBc( region,patch,iBegV,iEndV,iBegG,iEndG, &
                                  var,gradi,gradj,gradk )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch, ONLY   : t_patch
    TYPE(t_region)        :: region
    TYPE(t_patch)         :: patch
    INTEGER               :: iBegV, iEndV, iBegG, iEndG
    REAL(RFREAL), POINTER :: var(:,:), gradi(:,:), gradj(:,:), gradk(:,:)
  END SUBROUTINE RFLO_CalcGradPhysBc

  SUBROUTINE RFLO_CalcGradDummy( region,patch,iBegV,iEndV,iBegG,iEndG, &
                                 gradi,gradj,gradk )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    USE ModBndPatch, ONLY   : t_patch
    TYPE(t_region)        :: region
    TYPE(t_patch)         :: patch
    INTEGER               :: iBegV, iEndV, iBegG, iEndG
    REAL(RFREAL), POINTER :: gradi(:,:), gradj(:,:), gradk(:,:)
  END SUBROUTINE RFLO_CalcGradDummy

  SUBROUTINE RFLO_CalcGradDummyPhys( region,lbound,idir,jdir,kdir, &
                                     indBeg,indEnd,jndBeg,jndEnd,kndBeg,kndEnd, &
                                     iBegV,iEndV,iBegG,iEndG, &
                                     gradi,gradj,gradk )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)        :: region
    INTEGER               :: lbound,idir,jdir,kdir
    INTEGER               :: indBeg, indEnd, jndBeg, jndEnd, kndBeg, kndEnd
    INTEGER               :: iBegV, iEndV, iBegG, iEndG
    REAL(RFREAL), POINTER :: gradi(:,:), gradj(:,:), gradk(:,:)
  END SUBROUTINE RFLO_CalcGradDummyPhys

  SUBROUTINE RFLO_CalcGradDummySymm( region,lbound,idir,jdir,kdir, &
                                     inode,jnode,knode,ibeg,iend,jbeg,jend, &
                                     kbeg,kend,indBeg,indEnd,jndBeg,jndEnd, &
                                     kndBeg,kndEnd,iBegV,iEndV,iBegG,iEndG, &
                                     gradi,gradj,gradk )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)        :: region
    INTEGER               :: lbound,idir,jdir,kdir,inode,jnode,knode
    INTEGER               :: ibeg, iend, jbeg, jend, kbeg, kend
    INTEGER               :: indBeg, indEnd, jndBeg, jndEnd, kndBeg, kndEnd
    INTEGER               :: iBegV, iEndV, iBegG, iEndG
    REAL(RFREAL), POINTER :: gradi(:,:), gradj(:,:), gradk(:,:)
  END SUBROUTINE RFLO_CalcGradDummySymm

  SUBROUTINE RFLO_CalcGradDummyConn( region,lbound,idir,jdir,kdir, &
                                     indBeg,indEnd,jndBeg,jndEnd,kndBeg,kndEnd, &
                                     iBegV,iEndV,iBegG,iEndG,gradi,gradj,gradk )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)        :: region
    INTEGER               :: lbound,idir,jdir,kdir
    INTEGER               :: indBeg, indEnd, jndBeg, jndEnd, kndBeg, kndEnd
    INTEGER               :: iBegV, iEndV, iBegG, iEndG
    REAL(RFREAL), POINTER :: gradi(:,:), gradj(:,:), gradk(:,:)
  END SUBROUTINE RFLO_CalcGradDummyConn

  SUBROUTINE RFLO_CopyEdgeFaceNorm( region,iFBeg,iFEnd,fvari,fvarj,fvark )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: iFBeg, iFEnd
    REAL(RFREAL), POINTER :: fvari(:,:), fvarj(:,:), fvark(:,:)
  END SUBROUTINE RFLO_CopyEdgeFaceNorm

  SUBROUTINE RFLO_CopyEdgeFaceParal( region,iFBeg,iFEnd,fvari,fvarj,fvark )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
    INTEGER        :: iFBeg, iFEnd
    REAL(RFREAL), POINTER :: fvari(:,:), fvarj(:,:), fvark(:,:)
  END SUBROUTINE RFLO_CopyEdgeFaceParal

  SUBROUTINE RFLO_ChangeInteriorGrid( region,boundMoved,edgeMoved, &
                                      arcLen12,arcLen34,arcLen56, &
                                      xyzOld,xyz )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    LOGICAL :: boundMoved(6), edgeMoved(12)
    REAL(RFREAL), POINTER :: arcLen12(:,:), arcLen34(:,:), arcLen56(:,:)
    REAL(RFREAL), POINTER :: xyz(:,:), xyzOld(:,:)
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_ChangeInteriorGrid

  SUBROUTINE RFLO_CheckValidity( region )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), TARGET :: region
  END SUBROUTINE RFLO_CheckValidity

  SUBROUTINE RFLO_DerivedInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_DerivedInputValues

  SUBROUTINE RFLO_EdgeDeformation( region,boundMoved,edgeMoved, &
                                   arcLen12,arcLen34,arcLen56,xyzOld,dNode )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    LOGICAL :: boundMoved(6), edgeMoved(12)
    REAL(RFREAL), POINTER :: arcLen12(:,:), arcLen34(:,:), arcLen56(:,:)
    REAL(RFREAL), POINTER :: dNode(:,:), xyzOld(:,:)
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_EdgeDeformation

  SUBROUTINE RFLO_EdgeDeformationStraight( region,boundMoved,boundFlat, &
                                    edgeMoved,arcLen12,arcLen34,arcLen56, &
                                    xyzOrig,xyzOld,dNode )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    LOGICAL :: boundMoved(6), boundFlat(6), edgeMoved(12)
    REAL(RFREAL), POINTER :: arcLen12(:,:), arcLen34(:,:), arcLen56(:,:)
    REAL(RFREAL), POINTER :: dNode(:,:), xyzOld(:,:), xyzOrig(:,:)
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_EdgeDeformationStraight

  SUBROUTINE RFLO_ExchangeDnodeCopy( region,regionSrc,patch,patchSrc, &
                                     average,dNode,dNodeSrc )
    USE ModDataTypes
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    LOGICAL :: average
    REAL(RFREAL), POINTER :: dNode(:,:), dNodeSrc(:,:)
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch, patchSrc
  END SUBROUTINE RFLO_ExchangeDnodeCopy

  SUBROUTINE RFLO_ExchangeDnodeRecv( region,regionSrc,patch,patchSrc, &
                                     average,dNode )
    USE ModDataTypes
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    LOGICAL :: average
    REAL(RFREAL), POINTER :: dNode(:,:)
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch, patchSrc
  END SUBROUTINE RFLO_ExchangeDnodeRecv

  SUBROUTINE RFLO_ExchangeDnodeSend( region,regionSrc,patch,dNode )
    USE ModDataTypes
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    REAL(RFREAL), POINTER :: dNode(:,:)
    TYPE(t_region) :: region, regionSrc
    TYPE(t_patch)  :: patch
  END SUBROUTINE RFLO_ExchangeDnodeSend

  SUBROUTINE RFLO_ExtrapIntCellScal( region,fVec )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)          :: region
    REAL(RFREAL), POINTER   :: fVec(:)
  END SUBROUTINE RFLO_ExtrapIntCellScal

  SUBROUTINE RFLO_ExtrapIntCellVec( region,idBeg,idEnd,fVec )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)          :: region
    INTEGER                 :: idBeg,idEnd
    REAL(RFREAL), POINTER   :: fVec(:,:)
  END SUBROUTINE RFLO_ExtrapIntCellVec

  SUBROUTINE RFLO_GenerateCoarseGrids( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GenerateCoarseGrids

  SUBROUTINE RFLO_GetDimensDummy( region,iLev,idcbeg,idcend,jdcbeg,jdcend, &
                                  kdcbeg,kdcend )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetDimensDummy

  SUBROUTINE RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend,jpcbeg,jpcend, &
                                 kpcbeg,kpcend )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetDimensPhys

  SUBROUTINE RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend,jdnbeg,jdnend,&
                                       kdnbeg,kdnend )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetDimensDummyNodes

  SUBROUTINE RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend,jpnbeg,jpnend, &
                                      kpnbeg,kpnend )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetDimensPhysNodes

  SUBROUTINE RFLO_GetCellOffset( region,iLev,iCellOffset,ijCellOffset )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, iCellOffset, ijCellOffset
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetCellOffset

  SUBROUTINE RFLO_GetNodeOffset( region,iLev,iNodeOffset,ijNodeOffset )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, iNodeOffset, ijNodeOffset
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetNodeOffset

  SUBROUTINE RFLO_GetCornerCellsIndices( region,iLev,icorner, &
                                         icbeg,icend,jcbeg,jcend,kcbeg,kcend )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, icorner, icbeg, icend, jcbeg, jcend, kcbeg, kcend
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetCornerCellsIndices

  SUBROUTINE RFLO_GetEdgeCellsIndices( region,iLev,iedge, &
                                       iebeg,ieend,jebeg,jeend,kebeg,keend )
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, iedge
    INTEGER        :: iebeg, ieend, jebeg, jeend, kebeg, keend
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GetEdgeCellsIndices

  SUBROUTINE RFLO_GetPatchDirection( patch,idir,jdir,kdir )
    USE ModBndPatch, ONLY : t_patch
    INTEGER       :: idir, jdir, kdir
    TYPE(t_patch) :: patch
  END SUBROUTINE RFLO_GetPatchDirection

  SUBROUTINE RFLO_GetPatchIndices( region,patch,iLev, &
                                   ibeg,iend,jbeg,jend,kbeg,kend )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, ibeg, iend, jbeg, jend, kbeg, kend
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RFLO_GetPatchIndices

  SUBROUTINE RFLO_GetPatchIndicesNodes( region,patch,iLev,ibeg,iend, &
                                        jbeg,jend,kbeg,kend )
    USE ModBndPatch, ONLY   : t_patch
    USE ModDataStruct, ONLY : t_region
    INTEGER        :: iLev, ibeg, iend, jbeg, jend, kbeg, kend
    TYPE(t_region) :: region
    TYPE(t_patch)  :: patch
  END SUBROUTINE RFLO_GetPatchIndicesNodes

  SUBROUTINE RFLO_GetPatchMapping( lb,lbs,l1SrcDir,l2SrcDir,align, &
                                   idir,jdir,kdir,idirSrc,jdirSrc,kdirSrc, &
                                   ibeg,iend,jbeg,jend,kbeg,kend, &
                                   ibegSrc,iendSrc,jbegSrc,jendSrc, &
                                   kbegSrc,kendSrc,mapMat )
    INTEGER :: lb, lbs, l1SrcDir, l2SrcDir
    INTEGER :: idir, jdir, kdir, idirSrc, jdirSrc, kdirSrc
    INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
    INTEGER :: ibegSrc, iendSrc, jbegSrc, jendSrc, kbegSrc, kendSrc
    INTEGER :: mapMat(3,4)
    LOGICAL :: align
  END SUBROUTINE RFLO_GetPatchMapping

  SUBROUTINE RFLO_GridRemesh( region )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region) :: region
  END SUBROUTINE RFLO_GridRemesh

  SUBROUTINE RFLO_InitInputValues( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_InitInputValues

  SUBROUTINE RFLO_InterpolDistrib( n1f,n2f,n1c,n2c,nData,valf,valc )
    USE ModDataTypes
    INTEGER :: n1f, n2f, n1c, n2c, nData
    REAL(RFREAL), POINTER :: valf(:,:), valc(:,:)
  END SUBROUTINE RFLO_InterpolDistrib

  SUBROUTINE RFLO_RandomInit( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_RandomInit

  SUBROUTINE RFLO_ReadRegionTopology( global,regions )
    USE ModDataStruct, ONLY : t_region
    USE ModGlobal, ONLY     : t_global
    TYPE(t_global), POINTER :: global
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadRegionTopology

  SUBROUTINE RFLO_ReadDataFileInt( global,fileId,form,nDim1,nDim2,ivar )
    USE ModGlobal, ONLY : t_global
    INTEGER :: fileId, form, nDim1, nDim2
    INTEGER :: ivar(:,:)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLO_ReadDataFileInt

  SUBROUTINE RFLO_ReadDataFileReal( global,fileId,form,nDim1,nDim2,var )
    USE ModGlobal, ONLY : t_global
    USE ModDataTypes
    INTEGER :: fileId, form, nDim1, nDim2
    REAL(RFREAL) :: var(:,:)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLO_ReadDataFileReal

  SUBROUTINE RFLO_ReadGrid( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadGrid

  SUBROUTINE RFLO_ReadGridRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadGridRegion

  SUBROUTINE RFLO_ReadRandomState( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadRandomState

  SUBROUTINE RFLO_ReadSolution( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadSolution

  SUBROUTINE RFLO_ReadSolutionRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadSolutionRegion

  SUBROUTINE RFLO_SetMstageCoeffs( global,input,nrkStages )
    USE ModGlobal, ONLY: t_global
    USE ModMixture, ONLY : t_mixt_input
    TYPE(t_mixt_input)  :: input
    INTEGER :: nrkSteps
    TYPE(t_global), POINTER :: global    
  END SUBROUTINE RFLO_SetMstageCoeffs

  SUBROUTINE RFLO_Tfint1d( s,p1,p2,xyz )
    USE ModDataTypes
    REAL(RFREAL) :: s, p1(3), p2(3), xyz(3)
  END SUBROUTINE RFLO_Tfint1d

  SUBROUTINE RFLO_Tfint2d( s1,s2,s3,s4,e1,e2,e3,e4,p1,p2,p3,p4,xyz )
    USE ModDataTypes
    REAL(RFREAL) :: s1, s2, s3, s4
    REAL(RFREAL) :: e1(3), e2(3), e3(3), e4(3)
    REAL(RFREAL) :: p1(3), p2(3), p3(3), p4(3)
    REAL(RFREAL) :: xyz(3)
  END SUBROUTINE RFLO_Tfint2d

  SUBROUTINE RFLO_WriteDataFileInt( global,fileId,form,nDim1,nDim2,ivar )
    USE ModGlobal, ONLY : t_global
    INTEGER :: fileId, form, nDim1, nDim2
    INTEGER :: ivar(:,:)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLO_WriteDataFileInt

  SUBROUTINE RFLO_WriteDataFileReal( global,fileId,form,nDim1,nDim2,var )
    USE ModGlobal, ONLY : t_global
    USE ModDataTypes
    INTEGER :: fileId, form, nDim1, nDim2
    REAL(RFREAL) :: var(:,:)
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLO_WriteDataFileReal

  SUBROUTINE RFLO_WriteGrid( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_WriteGrid

  SUBROUTINE RFLO_WriteGridRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_WriteGridRegion

  SUBROUTINE RFLO_WriteRandomState( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_WriteRandomState

  SUBROUTINE RFLO_WriteRegionTopology( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_WriteRegionTopology

  SUBROUTINE RFLO_WriteSolution( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_WriteSolution

  SUBROUTINE RFLO_WriteSolutionRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_WriteSolutionRegion

#ifdef STATS
  SUBROUTINE RFLO_ReadStat( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadStat

  SUBROUTINE RFLO_ReadStatRegion( iReg,regions )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: iReg
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_ReadStatRegion

  SUBROUTINE RFLO_WriteStat( regions )
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region), POINTER :: regions(:)
  END SUBROUTINE RFLO_WriteStat
#endif

  SUBROUTINE RFLO_ZeroDummyCells( region,var )
    USE ModDataTypes
    USE ModDataStruct, ONLY : t_region
    TYPE(t_region)         :: region
    REAL(RFREAL), POINTER  :: var(:,:)
  END SUBROUTINE RFLO_ZeroDummyCells

  END INTERFACE

END MODULE RFLO_ModInterfacesLibrary

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModInterfacesLibrary.F90,v $
! Revision 1.23  2008/12/06 08:44:16  mtcampbe
! Updated license.
!
! Revision 1.22  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.21  2006/03/18 11:07:55  wasistho
! removed RFLO_GridFlatPatch
!
! Revision 1.20  2006/03/12 22:00:54  wasistho
! added RFLO_GridFlatPatch
!
! Revision 1.19  2006/03/12 20:39:40  wasistho
! added RFLO_EdgeDeformationStraight
!
! Revision 1.18  2005/11/11 07:16:27  wasistho
! removed RFLO_WriteDegeneratEC
!
! Revision 1.17  2005/10/27 05:55:56  wasistho
! removed RFLO_LaplaceGrid... routines
!
! Revision 1.16  2005/10/20 06:54:13  wasistho
! added RFLO_CalcCellCentroids and RFLO_CalcFaceCentroids
!
! Revision 1.15  2005/05/27 01:52:55  wasistho
! added rflo_gridremesh
!
! Revision 1.14  2004/11/17 16:28:43  haselbac
! Adapted interface for RFLO_SetMStageCoeffs
!
! Revision 1.13  2004/11/09 10:53:23  wasistho
! added RFLO_readStatRegion
!
! Revision 1.12  2004/09/30 16:57:20  wasistho
! added RFLO_extrapIntCellVec/Scal routines
!
! Revision 1.11  2004/08/21 00:34:43  wasistho
! added RFLO_writeDegeneratEC
!
! Revision 1.10  2004/07/26 19:29:48  wasistho
! changed POINTER to TARGET in RFLO_CheckValidity
!
! Revision 1.9  2004/07/26 19:04:28  wasistho
! add RFLO_CheckValidity
!
! Revision 1.8  2003/11/21 22:36:37  fnajjar
! Update Random Number Generator
!
! Revision 1.7  2003/08/25 21:51:24  jblazek
! Full version of global grid motion scheme.
!
! Revision 1.6  2003/08/11 21:51:17  jblazek
! Added basic global grid smoothing scheme.
!
! Revision 1.5  2003/07/03 21:48:45  jblazek
! Implemented dual-time stepping.
!
! Revision 1.4  2003/03/14 22:05:11  jblazek
! Improved mesh motion algorithm - node movement exchaged between blocks.
!
! Revision 1.3  2003/02/17 19:31:12  jferry
! Implemented portable random number generator ModRandom
!
! Revision 1.2  2003/02/03 19:20:46  jblazek
! Added treatment of edge and corner cells for one processor.
!
! Revision 1.1  2002/12/27 22:07:14  jblazek
! Splitted up RFLO_ModInterfaces and ModInterfaces.
!
!******************************************************************************






