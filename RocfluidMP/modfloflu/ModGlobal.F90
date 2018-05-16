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
! Purpose: Define global variables and global input data.
!
! Description: none
!
! Notes: none
!
! ******************************************************************************
!
! $Id: ModGlobal.F90,v 1.133 2009/03/02 00:19:33 mtcampbe Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! ******************************************************************************

MODULE ModGlobal

  USE ModDataTypes
  USE ModMaterials, ONLY : t_material

  IMPLICIT NONE

! ******************************************************************************
! Type definition
! ******************************************************************************

  TYPE t_global

! ******************************************************************************
!   Input
! ******************************************************************************

! ==============================================================================
!   Case name, grid & solution format, directories
! ==============================================================================

    LOGICAL :: calcCellCtr, calcFaceCtr
#ifdef RFLO
    LOGICAL :: degenrtEc
#endif
    CHARACTER(CHRLEN) :: winName,caseName
    CHARACTER(CHRLEN) :: inDir,outDir
    INTEGER :: gridFormat,solutFormat
    INTEGER :: gridSource,initFlowFlag

#ifdef GENX
    CHARACTER(CHRLEN) :: outDirHDF
#ifdef PLAG
    CHARACTER(CHRLEN) :: winp
#endif
#endif

! ==============================================================================
!   Reference values
! ==============================================================================

    REAL(RFREAL) :: refVelocity,refPressure,refDensity,refCp,refGamma, &
                    refLength,refREnum,refVisc,prLam,prTurb,scnLam,scnTurb


! ==============================================================================
!   Reference values for Gas-Liquid Mixture
! ==============================================================================

    REAL(RFREAL) :: refBetaPLiq,refBetaTLiq,refCvLiq,refDensityLiq, &
                    refPressLiq,refTempLiq

! ==============================================================================
!   Acceleration
! ==============================================================================

    LOGICAL      :: accelOn
    REAL(RFREAL) :: accelX,accelY,accelZ
    REAL(RFREAL) :: gravity

! ==============================================================================
!   Time stepping
! ==============================================================================

    LOGICAL :: dualTstSource,predCorrIter,predictSol,dtFixed
    INTEGER :: flowType,maxSubIter,solverType,tstepOrder
    INTEGER :: nrkSteps,currentIter,maxIter,printIter,rkScheme,writeIter
    REAL(RFREAL) :: currentTime,dTimeSystem,dtImposed,maxTime,printTime, &
                    resTol,stopRun,timeStamp,timeStampPrep,tolSubIter, &
                    writeTime
                    
#ifdef RFLU

    LOGICAL :: restartFromScratch
    INTEGER :: iterSinceRestart,restartIter
    REAL(RFREAL) :: currentTimeRK,restartTime,timeSincePrint,timeSinceProbe, &
                    timeSinceRestart,timeSinceWrite,ZoomFactor

! ==============================================================================
!   Time Zooming
! ==============================================================================

    REAL(RFREAL) :: tzMinPlane, tzMaxPlane, tzRhos, tzA, tzN, tzThroatRad, tzLenChamb, &
                    tzEpsNozz, tzVolBulk, tzDvolBulkDt, tzRadChamber, &
                    tzNozInlet, tzNozRad
    REAL(RFREAL), POINTER :: tzCvBulk(:),tzResBarQBulk(:),tzCvdVdtBulk(:)
    LOGICAL :: tzSubNoz
    INTEGER :: tzCoordLong, tzCoordTrans1, tzCoordTrans2

!     REAL(RFREAL),Dimension(CV_MIXT_DENS:CV_MIXT_ENER) :: CvBulk,resbar_Qbulk,CvdVdtBulk

! ==============================================================================
!   Constraints
! ==============================================================================

    REAL(RFREAL) :: cnstrCaseRad, cnstrTol1, cnstrTol2, cnstrEllipsL, cnstrEllipsT
    REAL(RFREAL) :: cnstrHeadEnd, cnstrAftEnd, cnstrNozY
    REAL(RFREAL) :: cnstrLMinPlane, cnstrLMaxPlane, cnstrT1MinPlane
    REAL(RFREAL) :: cnstrT1MaxPlane, cnstrT2MinPlane, cnstrT2MaxPlane
    INTEGER :: cnstrCoordL, cnstrCoordT1, cnstrCoordT2


#endif

! ==============================================================================
!   Multigrid
! ==============================================================================

    INTEGER :: startLevel,cycleType,refineIter

! ==============================================================================
!   Grid motion
! ==============================================================================

    INTEGER :: moveGridScheme, moveGridNiter, moveGridViter, moveGridSiter
    INTEGER :: moveGridNbour, moveGridRegNc, moveGridNsmatch, moveGridNsharedMax
    INTEGER :: moveGridOrthDir
    REAL(RFREAL) :: moveGridWeight, moveGridPower
    REAL(RFREAL) :: moveGridAmplifX, moveGridAmplifY, moveGridAmplifZ
    REAL(RFREAL) :: moveGridOrthWghtX, moveGridOrthWghtY, moveGridOrthWghtZ
    REAL(RFREAL) :: moveGridOrthCell

! ==============================================================================
!   Probe
! ==============================================================================

    LOGICAL :: probeOpenClose
    INTEGER :: nProbes, probeSaveIter
    INTEGER, POINTER :: probePos(:,:)
    REAL(RFREAL) :: probeSaveTime
    REAL(RFREAL), POINTER :: probeXYZ(:,:)

! ==============================================================================
!   Thrust
! ==============================================================================

    LOGICAL :: thrustOpenClose
    INTEGER :: thrustType,thrustPlane,thrustSaveIter
    REAL(RFREAL) :: thrustCoord,thrustMom,thrustPamb,thrustPress, & 
                    thrustSaveTime,thrustTotal

! ==============================================================================
!   Boundary conditions
! ==============================================================================

#ifdef RFLO
    INTEGER :: infloNijk, internDeform
    REAL(RFREAL), POINTER :: infloPlanEdges(:,:,:), xyzMinmax(:,:)
#endif

#ifdef STATS
! ==============================================================================
!   Statistics
! ==============================================================================

#ifdef GENX
    CHARACTER(CHRLEN), POINTER :: mixtStatNm(:,:,:), turbStatNm(:,:,:), &
                                  plagStatNm(:,:,:), peulStatNm(:,:,:)
#endif
    INTEGER :: doStat,reStat,mixtNStat,turbNStat,plagNStat,peulNStat,statBc
    INTEGER, POINTER :: mixtStatId(:,:),mixtStatCode(:,:,:)
    INTEGER, POINTER :: turbStatId(:,:),turbStatCode(:,:,:)
    INTEGER, POINTER :: plagStatId(:,:),plagStatCode(:,:,:)
    INTEGER, POINTER :: peulStatId(:,:),peulStatCode(:,:,:)
    REAL(RFREAL) :: integrTime
#endif

! ==============================================================================
!   Checking and verbosity level, warning counter
! ==============================================================================

    INTEGER :: checkLevel,verbLevel,warnCounter
#ifdef GENX
    INTEGER :: verbLevelCOM
#endif

! ==============================================================================
!   Transformation of grid and solution, enforcement of patch coordinates
! ==============================================================================

    LOGICAL  :: distortFlag,enforceFlag,transformFlag
    REAL(RFREAL) :: angleX,angleY,angleZ,distortX,distortY,distortZ,scaleX, &
                    scaleY,scaleZ


! ******************************************************************************
!   Variables
! ******************************************************************************

! ==============================================================================
!   MPI related data
! ==============================================================================

    INTEGER :: mpiComm,mpiTagMax,nProcs
    INTEGER :: nProcAlloc,myProcid,nRegionsProc,nRequests,iRequest
    INTEGER, POINTER :: requests(:)    
!    INTEGER :: nProcs
    INTEGER, DIMENSION(:), POINTER :: proc2RegMap,regMap
    INTEGER, DIMENSION(:,:), POINTER :: proc2RegMapInfo    

! ==============================================================================
!   GENX - NOTE always need to define
! ==============================================================================

    CHARACTER(CHRLEN) :: surfWinName,surfWinNameInput,volWinName, &
                         volWinNameInput,winNameIn,winNameOut
    INTEGER :: communicator,handleObtain
    INTEGER, DIMENSION(:), POINTER :: reg2PaneMap    

! ==============================================================================
!   Overall problem dimensions
! ==============================================================================

    INTEGER :: nRegions,nRegionsLocal
    INTEGER :: nLevels,nPatches

! ==============================================================================
!   Time stepping
! ==============================================================================

    REAL(RFREAL) :: dtMin,dtMinLimit,residual,resInit

! ==============================================================================
!   Limiter
! ==============================================================================

    REAL(RFREAL) :: limRef(3),limVolRef

! ==============================================================================
!   Optimal LES
! ==============================================================================

    REAL(RFREAL) :: dissOLES,enerOLES,uVarOLES,vVarOLES,wVarOLES

! ==============================================================================
!   Error codes, function tree
! ==============================================================================

    INTEGER :: error,mpierr,nFunTree
    CHARACTER(CHRLEN) :: functionTree(2,30)

! ==============================================================================
!   Grid quality
! ==============================================================================

    REAL(RFREAL) :: skewness, minVol

! ==============================================================================
!   Forces
! ==============================================================================

    REAL(RFREAL) :: forceX,forceY,forceZ
    REAL(RFREAL) :: forceRefArea,forceRefLength,forceRefXCoord,forceRefYCoord, &
                    forceRefZCoord
#ifdef RFLO
    INTEGER :: forcesOn, aeroCoeffs
    REAL(RFREAL) :: forceCoeffs(3,2), momentCoeffs(3,2)
    REAL(RFREAL) :: acBndBoxXmin, acBndBoxYmin, acBndBoxZmin
    REAL(RFREAL) :: acBndBoxXmax, acBndBoxYmax, acBndBoxZmax
#endif
#ifdef RFLU
    LOGICAL :: forceFlag,patchCoeffFlag
    INTEGER :: forceWriteCntr,thrustWriteCntr
#endif

! ==============================================================================
!   Mass flow
! ==============================================================================

    REAL(RFREAL) :: massIn,massOut

! ==============================================================================
!   Total volume and mass
! ==============================================================================

    REAL(RFREAL) :: totalMass,totalVol

! ==============================================================================
!   Material properties
! ==============================================================================

    INTEGER :: nMaterials
    TYPE(t_material), POINTER :: materials(:)

! ==============================================================================
!   Portable random number generator
! ==============================================================================

    INTEGER :: randSeedOffset,randSeedType

! ==============================================================================
!   Helpful constants
! ==============================================================================

    REAL(RFREAL) :: pi,rad
#ifdef RFLU
    REAL(RFREAL) :: deg2rad,rad2deg
#endif

! ==============================================================================
!   Work-variables reserved for common use by physical modules
! ==============================================================================

    REAL(RFREAL) :: moduleVar(9)

! ==============================================================================
!   Miscellaneous
! ==============================================================================

    INTEGER :: progressCounter
    INTEGER :: genxHandleBc,genxHandleGm
    INTEGER :: moduleType
                    
! ==============================================================================
!   Multiphysics modules
! ==============================================================================

    LOGICAL :: inrtUsed,peulUsed,plagUsed,plagUsedSave,specUsed

#ifdef PLAG
! ------------------------------------------------------------------------------
!   Particle module
! ------------------------------------------------------------------------------

    INTEGER :: initPlagFlag,nPclsCommTot,nPclsMax
#endif

#ifdef PEUL
! ------------------------------------------------------------------------------
!   Smoke module
! ------------------------------------------------------------------------------

    REAL(RFREAL) :: peulResInit,peulResidual
#endif

#ifdef RADI
! ------------------------------------------------------------------------------
!   Radiation module
! ------------------------------------------------------------------------------

    LOGICAL :: radiActive
#endif

#ifdef SPEC
! ------------------------------------------------------------------------------
!   Species
! ------------------------------------------------------------------------------

    INTEGER :: nSpecies ! Needed for postprocessing
#endif

#ifdef TURB
! ------------------------------------------------------------------------------
!   Turbulence module
! ------------------------------------------------------------------------------

    LOGICAL :: turbActive,turbCalcWDist,turbWorkUnused
    INTEGER :: turbWorkDim,turbWallDim,turbCalcWDistFreq
    REAL(RFREAL) :: esg1Sum,esg4Sum,esg1Psum,esg4Psum
    REAL(RFREAL), POINTER :: turbWork1D(:)
#endif

! ==============================================================================
!   Preprocessing
! ==============================================================================

    LOGICAL, POINTER :: prepBcDefined(:)    
#ifdef RFLU
    LOGICAL :: syPePatchesFlag
    INTEGER :: prepPartMode
#endif

! ==============================================================================
!   Postprocessing
! ==============================================================================

    INTEGER :: postPlotType

#ifdef RFLO
    INTEGER :: postIter,postOutFmt
    LOGICAL :: postStatsFlag,postTurbFlag,postPlagFlag,postRadiFlag, & 
               postSpecFlag
    REAL(RFREAL) :: postTime
#endif
#ifdef RFLU
    LOGICAL :: postCompErrFlag,postDiscFlag,postExtractFlag,postGradFlag, &
               postLag2EulFlag,postMergeFlag,postPlotPatchFlag, &
               postPlotVolFlag,postSpecFlag,postVortFlag,postVortCoreFlag, &
               postWriteMergeFlag
    INTEGER :: postInterpOrder,postInterpType,postNFringes,postNServers, &
               postOutputFormat,postPartNumber,postPartNumberSave,postSchType
    REAL(RFREAL) :: postSchExp
#endif

#ifdef RFLU
! ==============================================================================
!   Picking
! ==============================================================================

    LOGICAL :: pickCoordFlag
    REAL(RFREAL) :: pickXCoordLow,pickXCoordUpp,pickYCoordLow, & 
                    pickYCoordUpp,pickZCoordLow,pickZCoordUpp
#endif

#ifdef RFLO
! ==============================================================================
!   Flo2Flu conversion
! ==============================================================================

    INTEGER :: tofluNPatches,tofluNHexs,tofluNVerts,tofluNbfMax,tofluNbnMax
    INTEGER :: tofluNFaces,tofluNEdges, tofluMaxBind
    INTEGER, POINTER :: tofluNbVerts(:),tofluNbFaces(:)
    INTEGER, POINTER :: tofluHex2v(:,:),tofluQuad2v(:,:,:),tofluBLoc2g(:,:)
    INTEGER, POINTER :: tofluIq(:),tofluBType(:,:)
    REAL(RFREAL), POINTER :: tofluXyz(:,:)
#endif

  END TYPE t_global

END MODULE ModGlobal

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModGlobal.F90,v $
! Revision 1.133  2009/03/02 00:19:33  mtcampbe
! Added some ifdefs around Rocflo to disable particle injection on INFLOW
! boundaries and added some checks around MPI tags utilizing a new global
! data item, global%mpiTagMax.
!
! Revision 1.132  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.131  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.130  2008/01/14 22:08:57  mtcampbe
! added planar constr
!
! Revision 1.129  2007/04/14 22:36:57  mtcampbe
! Updated for TZ w/generalized geometry and support for submerged nozzle
!
! Revision 1.128  2007/04/14 14:29:10  mtcampbe
! Updated for TZ and rocket case constraints
!
! Revision 1.127  2007/03/06 23:11:38  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.126  2006/10/20 21:30:14  mparmar
! Added gravity and thrustWriteCntr
!
! Revision 1.125  2006/08/04 03:03:57  haselbac
! Added grid distortion parameters
!
! Revision 1.124  2006/05/09 23:36:56  wasistho
! added dtFixed
!
! Revision 1.123  2006/05/05 17:18:52  haselbac
! Cosmetics only
!
! Revision 1.122  2006/03/26 20:21:50  haselbac
! Added reference values for GL model
!
! Revision 1.121  2006/03/25 21:45:58  haselbac
! Added syPePatchesFlag
!
! Revision 1.120  2006/03/24 04:52:55  wasistho
! added forceCoeffs and momentCoeffs in RFLO
!
! Revision 1.119  2006/03/18 13:27:38  wasistho
! added orthDir and orthWghtX,Y,Z
!
! Revision 1.118  2006/03/18 11:04:46  wasistho
! added minvol
!
! Revision 1.117  2006/03/15 06:41:18  wasistho
! added global skewness
!
! Revision 1.116  2006/03/10 01:44:21  wasistho
! changed coeffsBox to acBndBox
!
! Revision 1.115  2006/03/10 00:56:40  wasistho
! added aeroCoeffs and coeffBox
!
! Revision 1.114  2006/03/09 19:38:15  wasistho
! made forceVariables available for RFLO
!
! Revision 1.113  2006/03/08 06:37:35  wasistho
! added moveGridSiter and Viter
!
! Revision 1.112  2006/01/06 22:10:26  haselbac
! Added entry for postGradFlag
!
! Revision 1.111  2005/12/10 23:28:51  haselbac
! Renamed geom post variables to pick variables
!
! Revision 1.110  2005/12/10 16:55:09  haselbac
! Added postLag2EulFlag
!
! Revision 1.109  2005/12/01 17:10:51  fnajjar
! Added randSeedType in global
!
! Revision 1.108  2005/10/28 22:49:48  wasistho
! added moveGridOrthCell
!
! Revision 1.107  2005/10/28 19:17:49  haselbac
! Added postPlotPatchFlag
!
! Revision 1.106  2005/10/20 06:49:08  wasistho
! added calcFaceCtr
!
! Revision 1.105  2005/10/05 20:04:39  haselbac
! Added vars for ENSIGHT filter
!
! Revision 1.104  2005/09/30 01:01:19  wasistho
! added internDeform for rflo bc
!
! Revision 1.103  2005/09/20 23:17:42  wasistho
! added infloNijk for Rocflo
!
! Revision 1.102  2005/08/28 23:47:38  wasistho
! added orthoWght for block orthogonality of RFLO global-gridmotion
!
! Revision 1.101  2005/08/25 23:05:15  wasistho
! added moveGridNsharedMax
!
! Revision 1.100  2005/08/18 19:47:04  wasistho
! added moveGridNsmatch
!
! Revision 1.99  2005/08/10 00:33:48  haselbac
! Added postVortCoreFlag
!
! Revision 1.98  2005/08/09 00:58:05  haselbac
! Added patchCoeffFlag, postCompErrFlag, postWriteMergeFlag
!
! Revision 1.97  2005/07/25 12:21:39  haselbac
! Added postVortFlag, changed postSchExp to float
!
! Revision 1.96  2005/07/05 19:27:21  haselbac
! Added postSchType and postSchExp
!
! Revision 1.95  2005/07/01 15:12:39  haselbac
! Added verbLevelCOM
!
! Revision 1.94  2005/06/25 03:16:28  wasistho
! enabled nRegions /= nProcs in type 2 gridmotion
!
! Revision 1.93  2005/06/23 01:34:48  wasistho
! added moveGridNbour
!
! Revision 1.92  2005/06/04 01:01:45  wasistho
! distinguished to AMPLIFX,Y,Z
!
! Revision 1.91  2005/06/02 22:57:29  wasistho
! added moveGridAmplif and moveGridPower
!
! Revision 1.90  2005/05/21 01:44:20  wasistho
! added statBc
!
! Revision 1.89  2005/05/18 22:06:17  fnajjar
! Added vars for initial parallelization of PLAG
!
! Revision 1.88  2005/05/03 08:14:55  wasistho
! added xyzMinmax for rflo
!
! Revision 1.87  2005/05/02 18:03:57  wasistho
! added infloPlanEdges
!
! Revision 1.86  2005/05/01 14:19:30  haselbac
! Added postDiscFlag and postNFringes
!
! Revision 1.85  2005/04/15 15:06:29  haselbac
! Removed Charm/FEM variables, added winNameIn
!
! Revision 1.84  2005/03/31 16:52:32  haselbac
! Added currentTimeRK
!
! Revision 1.83  2004/12/29 23:26:39  wasistho
! prepared statistics for PLAG and PEUL
!
! Revision 1.82  2004/11/17 16:29:01  haselbac
! Added rkScheme
!
! Revision 1.81  2004/10/26 15:17:12  haselbac
! Added postExtractFlag
!
! Revision 1.80  2004/10/21 15:32:32  haselbac
! Added reg2PaneMap
!
! Revision 1.79  2004/10/19 19:28:42  haselbac
! Added GENX variables, clean-up
!
! Revision 1.78  2004/10/09 16:34:28  fnajjar
! Added initialization flag for PLAG
!
! Revision 1.77  2004/08/21 00:30:01  wasistho
! added logical degenrtEc
!
! Revision 1.76  2004/08/18 02:08:56  wasistho
! added new toflu variables
!
! Revision 1.75  2004/08/17 00:55:08  wasistho
! prepared for utilities/rocflo/toflu
!
! Revision 1.74  2004/07/28 15:29:20  jferry
! created global variable for spec use
!
! Revision 1.73  2004/07/27 20:39:59  wasistho
! added pointer prepBcDefined(:)
!
! Revision 1.72  2004/07/25 05:03:15  wasistho
! fixed bug postSpecFlag already made general
!
! Revision 1.71  2004/07/24 03:46:51  wasistho
! define global variables for Rocflo-post
!
! Revision 1.70  2004/07/21 14:55:12  haselbac
! Added postInterpType
!
! Revision 1.69  2004/07/02 22:03:05  fnajjar
! Added winp for PLAG running in Gen3
!
! Revision 1.68  2004/06/16 20:00:51  haselbac
! Added variables, cosmetics
!
! Revision 1.67  2004/06/07 23:06:53  wasistho
! added mixtStatNm and turbStatNm
!
! Revision 1.66  2004/04/20 20:44:51  wasistho
! added calcTurbWDistFreq within ifdef TURB
!
! Revision 1.65  2004/03/05 22:09:01  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.64  2004/02/26 21:16:17  wasistho
! added esg1Sum, esg4Sum, esg1Psum, esg4Psum
!
! Revision 1.63  2004/01/29 22:57:21  haselbac
! Added deg2rad and rad2deg
!
! Revision 1.62  2003/11/21 22:38:21  fnajjar
! Added Random Seed Offset and Active flags for PLAG and PEUL
!
! Revision 1.61  2003/10/29 21:38:52  haselbac
! Added new global variables for time stepping
!
! Revision 1.60  2003/10/15 02:41:26  haselbac
! Added dtMinLimit and field flag for min dt location
!
! Revision 1.59  2003/10/07 20:30:08  wasistho
! rocturb work space made from 2D to 1D
!
! Revision 1.57  2003/08/28 20:05:39  jblazek
! Added acceleration terms.
!
! Revision 1.56  2003/08/11 21:50:00  jblazek
! Splitted ModInterfaces into 4 sections.
!
! Revision 1.55  2003/08/07 15:31:55  haselbac
! Added and changed var names
!
! Revision 1.54  2003/08/01 22:11:48  wasistho
! radiWrite/turbWrite to radiActive/turbActive
!
! Revision 1.53  2003/07/22 03:22:24  wasistho
! Added logical write-parameter for RADI and TURB
!
! Revision 1.52  2003/07/22 02:00:30  haselbac
! Added postInterOrder and warnCounter
!
! Revision 1.51  2003/07/03 21:48:45  jblazek
! Implemented dual-time stepping.
!
! Revision 1.50  2003/06/02 17:11:32  jblazek
! Added computation of thrust.
!
! Revision 1.49  2003/05/29 17:28:43  jblazek
! Implemented Roe scheme.
!
! Revision 1.48  2003/05/05 18:40:39  haselbac
! Added plotType
!
! Revision 1.47  2003/05/01 14:08:30  haselbac
! Added fieldFlagCntr
!
! Revision 1.46  2003/04/29 21:49:36  haselbac
! Added surfDiverFlag
!
! Revision 1.45  2003/04/28 22:41:32  haselbac
! Added vars for post and prep modules
!
! Revision 1.44  2003/04/12 21:37:08  haselbac
! FEMRocfluGrid now under global
!
! Revision 1.43  2003/04/07 14:22:57  haselbac
! Added field flag for probe positions
!
! Revision 1.42  2003/04/04 21:05:00  jblazek
! Corrected bug in dumping out the solution.
!
! Revision 1.41  2003/03/29 03:27:26  wasistho
! install ROCPERI
!
! Revision 1.40  2003/03/15 17:45:54  haselbac
! Added more field flags, removed splitFace
!
! Revision 1.39  2003/03/11 15:57:46  jferry
! Created data type for material properties
!
! Revision 1.38  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
! Revision 1.37  2003/02/05 21:07:30  jblazek
! Coordinated stop of a run works now for MPI.
!
! Revision 1.36  2003/02/01 00:29:36  haselbac
! Added enforceFlag, transformFlag now LOGICAL
!
! Revision 1.35  2003/01/30 19:07:05  haselbac
! Added timeStampPrep variable
!
! Revision 1.34  2003/01/23 17:48:53  jblazek
! Changed algorithm to dump convergence, solution and probe data.
!
! Revision 1.33  2002/12/02 20:10:41  jblazek
! Moved RFLU_ModGrid inside ifdef RFLU in ScaleGridSpeeds.
!
! Revision 1.32  2002/11/15 21:24:41  haselbac
! Changed field flag for integrals and added totalMass, totalVol
!
! Revision 1.31  2002/11/08 21:23:50  haselbac
! Added fieldFlagMass
!
! Revision 1.30  2002/11/02 01:53:57  wasistho
! Added TURB statistics
!
! Revision 1.29  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.28  2002/09/11 17:04:39  jblazek
! Added number of current predictor-corrector cycle.
!
! Revision 1.27  2002/09/11 16:29:00  jblazek
! Integrated into GENX.
!
! Revision 1.26  2002/09/09 14:53:16  haselbac
! Added several variables for OLES
!
! Revision 1.25  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.24  2002/08/30 01:47:58  jblazek
! Added support for moving grids.
!
! Revision 1.23  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
! Revision 1.22  2002/07/25 15:14:36  haselbac
! Added progress counter
!
! Revision 1.21  2002/07/12 21:50:07  jblazek
! Added tool to split single grid into multiple regions.
!
! Revision 1.20  2002/06/14 21:34:32  wasistho
! Added time avg statistics
!
! Revision 1.19  2002/06/14 20:16:01  haselbac
! Added nRegionsLocal (as ModLocal deleted) and nPatches
!
! Revision 1.18  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.17  2002/05/07 18:52:53  haselbac
! Added transformation variables
!
! Revision 1.16  2002/04/11 18:51:16  haselbac
! Added flag for initialization of solution
!
! Revision 1.15  2002/03/29 23:15:22  jblazek
! Corrected bug in MPI send.
!
! Revision 1.14  2002/03/26 19:17:05  haselbac
! Added gridSource and changed dimensions of probeXYZ to 2
!
! Revision 1.13  2002/03/18 23:07:19  jblazek
! Finished multiblock and MPI.
!
! Revision 1.12  2002/03/01 16:48:06  haselbac
! Added nLevels
!
! Revision 1.11  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.10  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.9  2002/01/31 20:23:59  jblazek
! Added treatment of edge & corner cells.
!
! Revision 1.8  2002/01/16 22:03:34  jblazek
! Added time-stepping routines.
!
! Revision 1.7  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.6  2002/01/10 18:21:29  jblazek
! Added iteration number and initial residual to solution file.
!
! Revision 1.5  2002/01/02 16:20:19  jblazek
! Added flow initialization and dummy cell geometry.
!
! Revision 1.4  2001/12/21 23:06:15  haselbac
! Added nRegionsProc for ROCFLU
!
! Revision 1.3  2001/12/19 23:09:21  jblazek
! Added routines to read grid and solution.
!
! Revision 1.2  2001/12/04 16:43:27  jblazek
! Makefiles modified because the modules Global, BndPatch and Grid moved
! to the modfloflu directory.
!
! Revision 1.1  2001/12/04 00:07:00  jblazek
! Modules BndPatch, Global and Grid moved to modfloflu directory.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
! ******************************************************************************






