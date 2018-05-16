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
! Purpose: define data types related to Lagrangian particles.
!
! Description:
!
! Input:
!
! Output:
!
! Notes:
!
! ******************************************************************************
!
! $Id: ModPartLag.F90,v 1.42 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! ******************************************************************************

MODULE ModPartLag

  USE ModDataTypes
  
  IMPLICIT NONE

! ******************************************************************************
! Input 
! ******************************************************************************

  TYPE t_plag_PDF
    INTEGER              :: nbins,locmax   ! number of data, maxloc
    REAL(RFREAL)         :: valmax         ! max val
    REAL(RFREAL),POINTER :: pdfvalues(:,:) 
  END TYPE t_plag_PDF

  TYPE t_plag_input
    INTEGER :: nCont                ! Total Number of Constituents
    INTEGER :: nPclsMax             ! Maximum Number of Particles in Region
    INTEGER :: ejecModel            ! Ejection Model Type
    INTEGER :: injcDiamDist         ! Injection Model Type
    INTEGER :: intrplMixtModel      ! Interpolation Model Type for Mixture
    INTEGER :: nPclsBuffTot         ! Total Buffer Size for Patches
    INTEGER :: nPclsBuffCECellsMax  ! Maximum Buffer Size for Corner and Edge Cells
    INTEGER :: breakupModel         ! Breakup Model Type
    INTEGER :: breakupWebSwi        ! Weber Switch for Breakup Model
    INTEGER :: readStatus           ! Status of reading for input sections
    INTEGER :: nPclsIni             ! Number of Initial Particles
    INTEGER :: findPclMethod	    ! Method for finding particles

    INTEGER, POINTER, DIMENSION(:) :: materialIndex  ! index pointing to appropriate Material

    REAL(RFREAL) :: injcVelRatio      ! Injection Velocity Ratio based on Mixture
    REAL(RFREAL) :: spLoad            ! Superparticle Loading
    REAL(RFREAL) :: injcDiamMean, &   ! Mean Diameter
                    injcStdDev, &     ! Standard Deviation
                    injcDiamMax, &    ! Maximum Diameter
                    injcDiamMin       ! Minimum Diameter
    REAL(RFREAL) :: injcBeta          ! Beta Coefficient for Ejection
    REAL(RFREAL) :: breakupFac        ! Breakup Factor

    REAL(RFREAL) :: iniRandDiamMax,   &  ! Initial Maximum Diameter for Random State
                    iniRandDiamMin,   &  ! Initial Minimum Diameter for Random State
                    iniRandTempMax,   &  ! Initial Maximum Temperature for Random State
                    iniRandTempMin,   &  ! Initial Minimum Temperature for Random State
                    iniRandSpLoadMax, &  ! Initial Maximum Superparticle Loading for Random State
                    iniRandSpLoadMin, &  ! Initial Maximum Superparticle Loading for Random State
                    iniRandUMax,      &
                    iniRandUMin,      & 
                    iniRandVMax,      &
                    iniRandVMin,      & 
                    iniRandWMax,      &
                    iniRandWMin,      & 
                    iniRandXMax,      &
                    iniRandXMin,      & 
                    iniRandYMax,      &
                    iniRandYMin,      &                                        
                    iniRandZMax,      &
                    iniRandZMin                     
                    
    REAL(RFREAL), POINTER, DIMENSION(:) :: dens,spht,molw, & ! Input structural parameters
                                           injcMassFluxRatio,injcTemp,surftens
    REAL(RFREAL), POINTER, DIMENSION(:) :: iniPosX,   &  ! Initial Xposition
                                           iniPosY,   &  ! Initial Yposition
                                           iniPosZ,   &  ! Initial Zposition
                                           iniDiam,   &  ! Initial Diameter
                                           iniTemp,   &  ! Initial Temperature
                                           iniSpload, &  ! Initial Superparticle Loading
                                           iniVelX,   &  ! Initial Xvelocity
                                           iniVelY,   &  ! Initial Yvelocity
                                           iniVelZ       ! Initial Zvelocity

    CHARACTER(CHRLEN), POINTER, DIMENSION(:) :: materialName ! String represent Material

    TYPE(t_plag_PDF) :: PDF
  END TYPE t_plag_input

! ******************************************************************************
! Data 
! ******************************************************************************

  TYPE t_plag
    INTEGER :: nAiv, nArv, nCv, nDv, nEv, nTv
    INTEGER :: nPcls,nPclsPrev,nPclsMax,nPclsGlobal
    INTEGER :: nRequests,nRequestsMetrics,nRequestsCECells,nRequestsStat
    INTEGER :: nInrtSources
    INTEGER :: nextIdNumber
    INTEGER :: nCellsNzPcl,nCellsNzPclMax

    INTEGER, POINTER, DIMENSION(:) :: requests, requestsMetrics
    INTEGER, POINTER, DIMENSION(:) :: requestsI, requestsR
    INTEGER, POINTER, DIMENSION(:) :: requestsCECells
    INTEGER, POINTER, DIMENSION(:) :: requestsCECellsI,requestsCECellsR 
    INTEGER, POINTER, DIMENSION(:) :: requestsStat
    INTEGER, POINTER, DIMENSION(:) :: cvPlagMass, dvPlagVolu
    INTEGER, POINTER, DIMENSION(:) :: icgNzPcl,iPclPerCellCSR,iPclPerCellCSRInfo
    INTEGER, POINTER, DIMENSION(:,:) :: aiv,aivOld

    REAL(RFREAL), POINTER, DIMENSION(:,:) :: arv, arvOld
    REAL(RFREAL), POINTER, DIMENSION(:,:) :: cv,  cvOld, dv, tv
    REAL(RFREAL), POINTER, DIMENSION(:,:) :: rhs, rhsSum
    REAL(RFREAL), POINTER, DIMENSION(:,:) :: inrtSources

    REAL(RFREAL), POINTER, DIMENSION(:,:,:) :: fc
    REAL(RFREAL), POINTER, DIMENSION(:,:)   :: si, sj, sk
    REAL(RFREAL), POINTER, DIMENSION(:,:)   :: ev, tav
  END TYPE t_plag

! ******************************************************************************
! Tile data structure 
! ******************************************************************************

  TYPE t_tile_plag
    INTEGER :: nCv, nDv
    INTEGER,      POINTER, DIMENSION(:)   :: nPclsInjc
    INTEGER,      POINTER, DIMENSION(:)   :: cvTileMass

    REAL(RFREAL), POINTER, DIMENSION(:,:) :: cv, cvOld
    REAL(RFREAL), POINTER, DIMENSION(:,:) :: dv
    REAL(RFREAL), POINTER, DIMENSION(:,:) :: rhs, rhsSum
  END TYPE t_tile_plag

! ******************************************************************************
! Communication data structure 
! ******************************************************************************

  TYPE t_buffer_plag
    INTEGER :: iRequest
    INTEGER :: nBuffSize
    INTEGER :: nBuffSizeDes
    INTEGER :: nSendBuffI, nRecvBuffI
    INTEGER :: nSendBuffR, nRecvBuffR
    INTEGER :: nBuffSizeTot
    INTEGER :: nSendBuffTotI, nRecvBuffTotI
    INTEGER :: nSendBuffTotR, nRecvBuffTotR
    INTEGER :: iRequestStat
    INTEGER :: nSendBuffStat,nRecvBuffStat
    INTEGER,      POINTER, DIMENSION(:)   :: recvBuffI, sendBuffI
    INTEGER,      POINTER, DIMENSION(:,:) :: aiv, aivOld

    REAL(RFREAL), POINTER, DIMENSION(:)   :: recvBuffR, sendBuffR
    REAL(RFREAL), POINTER, DIMENSION(:,:) :: arv, arvOld, cv, cvOld, dv, tv
    REAL(RFREAL), POINTER, DIMENSION(:,:) :: rhs, rhsSum
    REAL(RFREAL), POINTER, DIMENSION(:)   :: recvBuffStat, sendBuffStat
  END TYPE t_buffer_plag

! ******************************************************************************
! Surface statistics 
! ******************************************************************************

  TYPE t_surfstats_plag
    INTEGER, DIMENSION(:), POINTER :: nHits
    REAL(RFREAL), DIMENSION(:,:), POINTER :: vars
  END TYPE t_surfstats_plag

END MODULE ModPartLag

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModPartLag.F90,v $
! Revision 1.42  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.41  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.40  2007/03/31 23:52:28  haselbac
! Renamed nPclsTotGlobal to nPclsGlobal
!
! Revision 1.39  2007/03/27 00:18:55  haselbac
! Added arrays for improved PLAG initialization
!
! Revision 1.38  2007/03/20 22:01:49  fnajjar
! Added nPclsTotGlobal
!
! Revision 1.37  2007/03/06 23:12:07  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.36  2006/10/26 15:01:55  fnajjar
! Added min-max random velocities
!
! Revision 1.35  2006/09/18 20:27:01  fnajjar
! Added injcTemp in datastructure
!
! Revision 1.34  2006/05/05 17:23:22  haselbac
! Cosmetics
!
! Revision 1.33  2006/04/07 15:19:18  haselbac
! Removed tabs
!
! Revision 1.32  2005/05/18 22:10:12  fnajjar
! Added nPclsPrev
!
! Revision 1.31  2005/04/25 18:39:08  luca1
! Imposed PDF from file option for random particle ejection
!
! Revision 1.30  2005/03/31 20:25:10  fnajjar
! Included definitions for initial particle velocities
!
! Revision 1.29  2005/02/16 14:41:08  fnajjar
! Included MPI-based datastructure for statistics
!
! Revision 1.28  2005/01/08 20:35:37  fnajjar
! Included statistics datastructure
!
! Revision 1.27  2004/12/21 15:02:32  fnajjar
! Added surface statistics infrastructure for PLAG
!
! Revision 1.26  2004/10/11 19:38:15  fnajjar
! Renamed ininPlag to nPclsIni to follow naming convention
!
! Revision 1.25  2004/10/09 16:35:23  fnajjar
! Added infrastructure for particle preptool in random state and removed 
! initFlag as obsolete
!
! Revision 1.24  2004/10/08 22:08:53  haselbac
! Added findPclMethod
!
! Revision 1.23  2004/08/20 23:26:11  fnajjar
! Included variables for Plag prep tool
!
! Revision 1.22  2004/07/28 18:54:27  fnajjar
! Added nPclsTot for dynamic memory reallocation
!
! Revision 1.21  2004/06/17 15:17:08  fnajjar
! Included ejecModel variable
!
! Revision 1.20  2004/06/16 22:55:33  fnajjar
! Rename injcModel and injcTimeCoeff to injcModel and injcDiamDist to 
! prepare for CRE kernel
!
! Revision 1.19  2004/03/10 23:08:34  fnajjar
! Added infrastructure for MPI corner-edge cells
!
! Revision 1.18  2004/03/05 22:09:01  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.17  2004/01/16 23:22:28  fnajjar
! Included requestsMetrics and nRequestsMetrics for MPI of geometry metrics
!
! Revision 1.16  2003/11/21 22:39:48  fnajjar
! Removed nPclsTot and added nextIdNumber for proper restart
!
! Revision 1.15  2003/11/03 21:18:43  fnajjar
! Added face vectors to PLAG infrastructure
!
! Revision 1.14  2003/09/17 21:06:53  fnajjar
! Added infrastructure for skewed Log distribution in injection model
!
! Revision 1.13  2003/09/13 20:16:26  fnajjar
! Added infrastructure for Breakup model
!
! Revision 1.12  2003/09/10 23:36:24  fnajjar
! Removed flags that are subsumed with Rocinteract
!
! Revision 1.11  2003/05/27 19:06:16  fnajjar
! Removed distPartBurning and all pertinent LOGICAL datastructure
!
! Revision 1.10  2003/03/12 21:20:11  fnajjar
! Include Material Datastructure
!
! Revision 1.9  2003/03/04 22:12:34  jferry
! Initial import of Rocinteract
!
! Revision 1.8  2003/02/06 16:16:12  f-najjar
! Added definition for requests of integer and real data buffers
!
! Revision 1.7  2003/01/22 23:59:36  f-najjar
! Remove iRequest for t_plag datastructure
!
! Revision 1.6  2003/01/22 23:43:06  f-najjar
! Included MPI-related request data
!
! Revision 1.5  2003/01/22 16:49:22  f-najjar
! Added iRequest for MPI communication
!
! Revision 1.4  2003/01/16 17:45:25  f-najjar
! Added nBuffSizeDes for buffer infrastructure
!
! Revision 1.3  2002/10/25 14:04:53  f-najjar
! Finalize PLAG datastructure
!
! Revision 1.2  2001/12/11 21:59:29  jblazek
! memory allocation added.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
! ******************************************************************************






