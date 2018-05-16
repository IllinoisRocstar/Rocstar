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
! ******************************************************************************
!
! Purpose: define derived data type related to boundary patches.
!
! Description: none
!
! Notes: none
!
! ******************************************************************************
!
! $Id: ModBndPatch.F90,v 1.70 2010/02/18 21:47:40 juzhang Exp $
!
! Copyright: (c) 2001-2006 by the University of Illinois
!
! ******************************************************************************

MODULE ModBndPatch

  USE ModDataTypes
  USE ModParameters
  USE ModPartLag, ONLY : t_tile_plag, t_buffer_plag, t_surfstats_plag
#ifdef PEUL
  USE ModPartEul, ONLY : t_buffer_peul
#endif

  USE ModStencil, ONLY: t_stencil,t_stencilInfo

  IMPLICIT NONE

! ******************************************************************************
! Boundary condition values data structure 
! ******************************************************************************

  TYPE t_bcvalues
#ifdef RFLO
    LOGICAL :: bcSet, setMotion
    INTEGER :: nSendBuff, nRecvBuff, iRequest
#endif
    INTEGER :: distrib, nData, nSwitches

    INTEGER, POINTER      :: switches(:)
    REAL(RFREAL), POINTER :: vals(:,:)
#ifdef RFLO
    REAL(RFREAL)          :: maxChange, bndVel(3), amplitude
    REAL(RFREAL), POINTER :: sendBuff(:), recvBuff(:)
#endif

    TYPE(t_tbcvalues), POINTER :: tbcs(:)

    INTEGER :: cvState
    REAL(RFREAL), POINTER, DIMENSION(:,:) :: cv,cvOld,dv,gv,tv,rhs,rhsSum
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: gradFace
  END TYPE t_bcvalues

! ******************************************************************************
! Time-dependent boundary condition data structure
! ******************************************************************************

  TYPE t_tbcvalues
    INTEGER :: tbcType
    INTEGER, POINTER      :: switches(:)
    REAL(RFREAL), POINTER :: params(:), mean(:), svals(:), bvals(:,:)
  END TYPE t_tbcvalues

! ******************************************************************************
! Patch data structure 
! ******************************************************************************

  TYPE t_patch
    INTEGER      :: bcKind,reflect
    REAL(RFREAL) :: nscbcK

    INTEGER      :: bcType, bcCoupled, bcMotion
    REAL(RFREAL) :: periodAngle     ! angle between rotat. periodic boundaries
#ifdef RFLO
    INTEGER :: srcRegion, srcPatch
    INTEGER :: lbound, l1beg, l1end, l2beg, l2end
    INTEGER :: srcLbound, srcL1beg, srcL1end, srcL2beg, srcL2end
    INTEGER :: corns(4), position(4)
    INTEGER :: dirFlat
    LOGICAL :: align, thrustCalc, bndFlat, globalAeroCoeffs
    REAL(RFREAL) :: l1VecSrc(3), l2VecSrc(3)
    REAL(RFREAL), POINTER :: surfCoord(:,:,:)
    REAL(RFREAL), POINTER :: arclen1(:), arclen2(:)
    REAL(RFREAL), POINTER :: cp(:), ch(:)    
    REAL(RFREAL), POINTER :: cf(:,:), forceCoeffs(:,:), momentCoeffs(:,:)    

    REAL(RFREAL), POINTER :: st(:,:,:), stOld(:,:,:)             ! prmgrd (EPDE)
    REAL(RFREAL), POINTER :: sti(:,:,:), stj(:,:,:)              ! 1st deriv. st
    REAL(RFREAL), POINTER :: stii(:,:,:),stjj(:,:,:),stij(:,:,:) ! 2nd deriv. st
    REAL(RFREAL), POINTER :: pfun(:,:,:,:)                       ! ctrl function

    REAL(RFREAL), POINTER :: aij(:,:), aimj(:,:), aipj(:,:)      ! matrix coeffs
    REAL(RFREAL), POINTER :: aimjm(:,:), aijm(:,:), aipjm(:,:)   ! of EPDE
    REAL(RFREAL), POINTER :: aimjp(:,:), aijp(:,:), aipjp(:,:)

#ifdef GENX
    INTEGER, POINTER      :: bcFlag(:), bFlag(:,:)
    REAL(RFREAL), POINTER :: duAlp(:,:,:), mdotAlp(:,:), &
                             rhofvfAlp(:,:,:), tflmAlp(:,:), pf(:,:)     , &
                             qc(:,:), qr(:,:), rhofAlp(:,:), nfAlp(:,:,:), &
                             tracf(:,:,:), tempf(:,:), tempv(:,:), dnml(:,:)
#endif
#endif
#ifdef RFLU
    CHARACTER(CHRLEN) :: bcName 
    LOGICAL :: flatFlag,movePatch,plotFlag,plotStatsFlag,renumFlag,smoothGrid, &
               thrustFlag,transformFlag
    LOGICAL, DIMENSION(:), POINTER :: nbMap
    INTEGER :: axisRelated,cReconst,iBorder,iPatchGlobal,iPatchLocal, &
               iPatchRelated,moveBcType,movePatchDir,spaceOrder
    INTEGER :: nBCellsVirt,nBFaces,nBFacesMax,nBFacesTot,nBTris,nBTrisMax, &
               nBTrisTot,nBQuads,nBQuadsMax,nBQuadsTot,nBVert,nBVertEst, &
               nBVertMax,nBVertTot
    INTEGER, DIMENSION(:), POINTER :: bf2c,bf2cSorted,bf2cSortedKeys,bv,bvc, &
                                      bvcSorted,bvTemp,bv2bv
    INTEGER, DIMENSION(:,:), POINTER :: bf2ct,bf2v,bf2vSorted,bTri2v,bTri2vLoc, & 
                                        bQuad2v,bQuad2vLoc    
    TYPE(t_stencil), DIMENSION(:), POINTER :: bf2cs,bf2cs1D  
    TYPE(t_stencilInfo) :: bf2csInfo,bf2cs1DInfo
    REAL(RFREAL) :: angleRelated          
    REAL(RFREAL) :: pc(XCOORD:ZCOORD),pn(XCOORD:ZCOORD)
    REAL(RFREAL) :: tm(XCOORD:XYZMAG,XCOORD:XYZMAG)
    REAL(RFREAL), DIMENSION(:), POINTER :: gs,mfMixt,vfMixt
    REAL(RFREAL), DIMENSION(:), POINTER :: cofgDist    
    REAL(RFREAL), DIMENSION(:,:), POINTER :: bvn,dXyz,fc,fn
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: bfgwt
    REAL(RFREAL), DIMENSION(:), POINTER :: cp,ch,cmass    
    REAL(RFREAL), DIMENSION(:,:), POINTER :: cf,cmom,forceCoeffs,forceVacCoeffs,momentCoeffs    
    REAL(RFREAL), DIMENSION(:), POINTER :: massCoeffs
#ifdef GENX
    INTEGER :: bFlagInit
    INTEGER, DIMENSION(:), POINTER :: bFlag,bcFlag,cnstrType,patchNo    
    REAL(RFREAL), DIMENSION(:), POINTER :: mdotAlp,pf,qc,qr,rhofAlp,tbAlp, & 
                                           tempf,tflmAlp  
    REAL(RFREAL), DIMENSION(:,:), POINTER :: duAlp,nfAlp,rhofvfAlp,tracf, & 
                                             xyz
#endif
#ifndef NO_TECPLOT
    REAL(RFREAL), DIMENSION(:,:), POINTER :: varFaceTEC,varVertTEC
#endif   
#endif 

#ifdef RADI
    INTEGER :: radBcType
#endif
#ifdef PEUL
    TYPE(t_buffer_peul) :: bufferPeul
#endif

    TYPE(t_bcvalues) :: mixt,turb,spec,peul,valRadi,valBola
    TYPE(t_bcvalues), POINTER :: plag(:)
    TYPE(t_tile_plag)         :: tilePlag
    TYPE(t_buffer_plag)       :: bufferPlag
    TYPE(t_surfstats_plag), POINTER, DIMENSION(:) :: statsPlag
  END TYPE t_patch

END MODULE ModBndPatch

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModBndPatch.F90,v $
! Revision 1.70  2010/02/18 21:47:40  juzhang
! Heat transfer bc for non-propellant surface documented in Rocburn_PY_HT.pdf in Rocburn_PY directory is implemented within Rocburn_PY. Major changes were made to Rocburn, Rocman3, RocfluidMP/genx, RocfluidMP/modflo directories.
!
! Revision 1.69  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.68  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.67  2006/10/20 21:28:13  mparmar
! Added thrustFlag, cmass, cmom, forceVacCoeffs, momentCoeffs, massCoeffs
!
! Revision 1.66  2006/08/19 15:52:21  mparmar
! Renamed patch variables, rm bg2bf, bf2bgTot and modified t_bcvalues, t_patch
!
! Revision 1.65  2006/08/18 13:59:38  haselbac
! Added transformFlag and bvcSorted
!
! Revision 1.64  2006/04/15 16:57:18  haselbac
! Added cReconst and spaceOrder to patch data type
!
! Revision 1.63  2006/04/07 14:43:50  haselbac
! Added iPatchLocal and bf2cs1D
!
! Revision 1.62  2006/03/25 21:44:29  haselbac
! Added vars bcos of sype patches
!
! Revision 1.61  2006/03/18 08:21:47  wasistho
! added arclen1,2
!
! Revision 1.60  2006/03/17 06:36:28  wasistho
! added dirFlat
!
! Revision 1.59  2006/03/10 08:54:28  wasistho
! changed dimensions of aero coeffs
!
! Revision 1.58  2006/03/10 00:57:04  wasistho
! added globalAeroCoeffs
!
! Revision 1.57  2005/12/07 08:46:34  wasistho
! added stuff for surface mesh motion EPDE
!
! Revision 1.56  2005/12/05 05:53:19  wasistho
! added bndFlat
!
! Revision 1.55  2005/11/18 07:17:45  wasistho
! added amplitude in t_bcvalues
!
! Revision 1.54  2005/10/14 14:03:09  haselbac
! Added tbAlp
!
! Revision 1.53  2005/09/23 18:56:58  haselbac
! Added plotStatsFlag
!
! Revision 1.52  2005/09/09 03:17:44  wasistho
! added valMixt%setMotion and valMixt%bndVel(3)
!
! Revision 1.51  2005/08/09 00:57:35  haselbac
! Added plotFlag
!
! Revision 1.50  2005/06/16 01:26:35  wasistho
! added patch%position
!
! Revision 1.49  2005/06/13 21:44:02  wasistho
! added new patch variable patch%bcMotion
!
! Revision 1.48  2005/06/10 22:32:04  wasistho
! define patch%corns(4)
!
! Revision 1.47  2005/06/09 20:18:27  haselbac
! Added movePatchDir and cnstrType
!
! Revision 1.46  2005/05/04 03:32:29  haselbac
! Removed writeGrid, was not required for some time
!
! Revision 1.45  2005/05/02 18:05:10  wasistho
! move surfCoord outside ifdef Genx
!
! Revision 1.44  2004/12/21 15:01:37  fnajjar
! Added surface statistics infrastructure for PLAG
!
! Revision 1.43  2004/12/19 15:43:20  haselbac
! Added vfMixt and cofgDist
!
! Revision 1.42  2004/10/19 19:28:35  haselbac
! Added GENX variables, cosmetics
!
! Revision 1.41  2004/09/30 18:04:00  wasistho
! added radBcType within ifdef RADI
!
! Revision 1.40  2004/06/16 20:00:46  haselbac
! Added patch coefficients and Tecplot variables
!
! Revision 1.39  2004/05/25 01:33:33  haselbac
! Added bf2bgTot
!
! Revision 1.38  2004/01/29 22:57:19  haselbac
! Removed vfMixt
!
! Revision 1.37  2003/12/04 03:28:22  haselbac
! Added boundary stencil arrays
!
! Revision 1.36  2003/11/25 21:03:07  haselbac
! Added mfMixt and vfMixt for rocspecies support
!
! Revision 1.35  2003/06/02 17:11:32  jblazek
! Added computation of thrust.
!
! Revision 1.34  2003/05/31 01:41:54  wasistho
! installed turb. wall layer model
!
! Revision 1.33  2003/05/19 21:18:21  jblazek
! Automated switch to 0th-order extrapolation at slip walls and injection 
! boundaries.
!
! Revision 1.32  2003/05/07 00:23:33  haselbac
! Changed surfGridFlag to writeGrid
!
! Revision 1.31  2003/04/29 22:42:05  haselbac
! Added surfGridFlag
!
! Revision 1.30  2003/04/24 15:38:31  haselbac
! Added bFlagInit
!
! Revision 1.29  2003/04/10 00:56:45  jblazek
! Added ifdef PEUL around TYPE(t_buffer_peul) :: bufferPeul
!
! Revision 1.28  2003/04/09 22:51:44  jferry
! removed peul_save and peul_verify structures
!
! Revision 1.27  2003/04/09 14:06:43  fnajjar
! Included Pointer Type of Buffers for PEUL
!
! Revision 1.26  2003/03/15 17:40:29  haselbac
! Added bvFlag, renumFlag, nB*Tot, removed bf2ct
!
! Revision 1.25  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
! Revision 1.24  2003/01/28 16:37:52  haselbac
! Removed locally-numbered lists
!
! Revision 1.23  2002/12/06 22:29:26  jblazek
! Corrected bug for geometry exchange between minimal patches.
!
! Revision 1.22  2002/11/26 15:24:38  haselbac
! Added moveBcType
!
! Revision 1.21  2002/10/27 18:59:30  haselbac
! Added bvn, smoothGrid, and movePatch, moved dXyz & bTri2vl
!
! Revision 1.20  2002/10/25 14:02:21  f-najjar
! Define Tile and Buffer types for PLAG datastructure
!
! Revision 1.19  2002/10/12 19:11:03  haselbac
! Added bTri2vl for RFLU
!
! Revision 1.18  2002/10/05 18:57:03  haselbac
! Integration of RFLU into GENX
!
! Revision 1.17  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.16  2002/09/24 23:16:46  jblazek
! Changed bcflag to a pointer.
!
! Revision 1.15  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.14  2002/09/17 22:51:23  jferry
! Removed Fast Eulerian particle type
!
! Revision 1.13  2002/09/17 13:43:00  jferry
! Added Time-dependent boundary conditions
!
! Revision 1.12  2002/09/09 14:50:47  haselbac
! Added access arrays and face gradient weights
!
! Revision 1.11  2002/06/22 01:13:37  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.10  2002/06/14 20:14:03  haselbac
! Added iPatchGlobal and bf2ct array (for Charm)
!
! Revision 1.9  2002/03/29 23:15:22  jblazek
! Corrected bug in MPI send.
!
! Revision 1.8  2002/03/26 19:13:25  haselbac
! Enclosed bcSet in conditional compilation for ROCFLO
!
! Revision 1.7  2002/03/18 23:07:19  jblazek
! Finished multiblock and MPI.
!
! Revision 1.6  2002/03/14 19:05:16  haselbac
! Added entries for face centroid and face normal
!
! Revision 1.5  2002/03/01 16:39:10  haselbac
! Added nBFaces,bf2c,bf2v arrays
!
! Revision 1.4  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.3  2001/12/19 23:09:21  jblazek
! Added routines to read grid and solution.
!
! Revision 1.2  2001/12/04 17:15:38  haselbac
! Added ROCFLU variables
!
! Revision 1.1  2001/12/04 00:07:00  jblazek
! Modules BndPatch, Global and Grid moved to modfloflu directory.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
! ******************************************************************************






