!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
MODULE ROCSTAR_RocFrac

!!****f* Rocfrac/Rocfrac/Source/ROCSTAR_RocFrac.f90
!!
!!  NAME
!!    ROCSTAR_RocFrac
!!
!!  FUNCTION
!!    Global variable module
!!
!!***

  USE ROCSTAR_RocFracComm 
  
  IMPLICIT NONE
  
  SAVE
  
  CHARACTER(*), PARAMETER :: surWin = "sfrac"
  CHARACTER(*), PARAMETER :: volWin = "vfrac"

! - Derive Type for the boundary conditions
  TYPE bcvalues
     INTEGER :: BCtypeX, BCtypeY, BCtypeZ
     REAL*8  :: BCvalueX, BCvalueY, BCvalueZ
  END TYPE bcvalues

  TYPE, PUBLIC :: ROCFRAC_GLOBAL
     LOGICAL   :: iDummyRocfrac


!----  ** PARALLEL PARAMETER **
     INTEGER :: MPI_COMM_ROCFRAC  ! Communicator of Rocfrac
! -- Processor id number character
     CHARACTER*4 :: MyIdChr
!--   Non-block receive, Non-block send request arrays
     INTEGER, POINTER, DIMENSION(:) :: ReqRcv,ReqSnd
!--   Non-block receive, Non-block send status arrays
     INTEGER, POINTER, DIMENSION(:,:) :: StatRcv, StatSnd
!--   Number of Volumetric elements on Partition Boundary
     INTEGER :: NumElPartBndry 
! -- Total Number of Nodes being communicated
     INTEGER :: TotNumNdComm
! -- Total Number of Neighoring Processor communicating with.
     INTEGER :: TotNumNeighProcs 
! -- Number of Nodes Communicating per-processor
     INTEGER, POINTER, DIMENSION(:) ::  NumNdComm
! -- Received Data buffer
     TYPE(rcv_buf), POINTER, DIMENSION(:) :: RecvDataFrm
! -- List of Processors Neighbors
     INTEGER, POINTER, DIMENSION(:) :: NeighProcList
! -- List of Nodes Communicated
     TYPE(send_buf), POINTER, DIMENSION(:) :: NdCommList


!---- ** FILE PRAMETERS **
     CHARACTER(LEN=200) :: prefx   ! I/O file prefix
     INTEGER :: prefx_lngth        ! length of prefix
     INTEGER :: io_input      ! unit id of control input deck
     INTEGER :: io_sum        ! unit id of summary control inpute deck


!---- ANALYSIS PARAMETERS
     LOGICAL :: ReStart      ! Is the analysis a continuation ?
                                    !   .FALSE. = NO *
                                    !   .TRUE. = YES
     REAL*8 :: CourantRatio    ! ratio of delta and Courant conditon

     REAL*8 :: DT           ! time step (seconds)
     REAL*8 :: DTInv 
 
! -- Current sytem time

     REAL*8 :: CurrTime            

     LOGICAL :: ALEenabled     ! Is the ALE enabled ?
                                    !   .FALSE. = NO *
                                    !   .TRUE. = YES

     INTEGER, POINTER, DIMENSION(:) :: iSolnType       ! Analysis Type
                                    !  -1 = non-linear geom/ neo-hookean
                                    !   0 = non-linear geom/ arruda Boyce
                                    !   1 = non-linear geom/ linear mat
                                    !   2 = linear geom / linear mat
                                    !  10 = 4node ANDG /arruda Boyce

     INTEGER, POINTER, DIMENSION(:) :: iSolnTypeHT       ! Analysis Type
                                    !   0 = no heat transfer solution
                                    !   1 = heat transfer solution

!--   heat conductivity lke variable for mesh motion 
     REAL*8 :: kappa

!-- Used to subtract the initial pressure
     REAL*8, POINTER, DIMENSION(:,:) :: pstatic
!-- Used to subtract the initial pressure
     REAL*8, POINTER, DIMENSION(:,:) :: pstaticnb
!-- Flag to subract the initial pressure
     LOGICAL :: ipstatic

!---- ** VOLUMETRIC MATERIAL PROPERTIES **
!
     INTEGER :: NumMatVol    ! number of volumetric materials
!--   Young's moduli | Poisson's ratios | density | thermal expansion coeffs
     REAL*8, POINTER, DIMENSION(:) :: E, xnu, rho, alpha
     REAL*8, POINTER, DIMENSION(:) :: E1, E2, E3, nu12, nu13, nu23, G12, G13, G23
! -- shear Moduls, Bulk Modulus, Lame's Constant
     REAL*8, POINTER, DIMENSION(:) :: xmu, xkappa, xlambda
!--   Fastest Dilation Wave Speed
     REAL*8 :: cd_fastest
!--   elastic stiffness const
     REAL*8, POINTER, DIMENSION(:,:) :: ci,cj

!---- ** Mesh Variables **
!

! --- Order of Tetradehral
        ! iElType = 4 ; Four Node Tetrahedral
        ! iElType = 10 ; Ten Node Tetrahedral
     INTEGER :: iElType
!  -- Corresponding order for face of Tetrahedral   
     INTEGER :: iElType2D
!  -- Number of Gaussian Integration Points
     INTEGER :: iStrGss

! Number of Nodes, Number of Elements     
     INTEGER :: NumNP, NumElVol
     INTEGER, POINTER, DIMENSION(:) ::  NumElVolMat, NumElPartBndryMat
!--   mat number for CST element
     INTEGER, POINTER, DIMENSION(:) :: MatIdVol
!--   connectivity table for CST elem
     INTEGER, POINTER, DIMENSION(:,:) :: ElConnVol
!--   Undeformed Coordinate array (i.e. Original Mesh Coordinates)
!     REAL*8, POINTER, DIMENSION(:,:) :: coor
!--   Mesh Coordinate array ( meshcoor = coor if ALEenabled = .False.)
     REAL*8, POINTER, DIMENSION(:,:) :: meshcoor
!  - Stores the boundary conditions data values
     TYPE(bcvalues), POINTER, DIMENSION(:) :: bcCond
!  - Stores the boundary conditions data values
     TYPE(bcvalues), POINTER, DIMENSION(:) :: bcCondHT

!- Physical Boundary Flags
     INTEGER :: NumNdsBC  ! number of boundary nodes w/loads
     INTEGER :: NumNdsBCmm  ! number of boundary nodes w/loads mesh motion
     INTEGER :: NumNdsBCHT  ! number of boundary nodes w/loads mesh motion
!--   displacement and force dof
     INTEGER, POINTER, DIMENSION(:,:) :: BCFlag
!--   displacement and force dof
     INTEGER, POINTER, DIMENSION(:,:) :: BCFlagHT
!--   applied displacements and loads
     REAL*8, POINTER, DIMENSION(:,:) :: BCvalue
!--   applied displacements and loads
     REAL*8, POINTER, DIMENSION(:,:) :: BCvalueHT

!  - Stores the boundary meshmotion conditions data values
     TYPE(bcvalues), POINTER, DIMENSION(:) :: bcCondmm
     INTEGER :: NumBcFlagsmm
!--   displacement and force dof mesh motion
     INTEGER, POINTER, DIMENSION(:,:) :: BCFlagmm
!--   applied displacements and loads
     REAL*8, POINTER, DIMENSION(:,:) :: BCvaluemm

! - Mass Matrix
     REAL*8, POINTER, DIMENSION(:) :: xmass

     REAL*8, POINTER, DIMENSION(:) :: CapctInv

!--- Additional Arrays for average nodal deformation gradient elements

!  - Number of Elements Associated with node
     INTEGER, POINTER, DIMENSION(:) :: NumELNeigh
!  - Element ids of Element Associated with node
     INTEGER, POINTER, DIMENSION(:,:) :: ElConnNd
!  - Volume Ratio Associated with node
     REAL*8, POINTER, DIMENSION(:,:) :: AlphaR
!  - Origonal undeformed Volume of node
     REAL*8, POINTER, DIMENSION(:) :: VolUndfmd

! ** Solid/Fluid Interface ** 

! -- Solid communication map, Maps Surface Node to Volume Node
     INTEGER, POINTER, DIMENSION(:) :: MapNodeSF
     INTEGER, POINTER, DIMENSION(:) :: MapNodeSFnb
     INTEGER, POINTER, DIMENSION(:) :: MapNodeS

! -- Number of Nodes and Elements  on the Solid/Fluid Interface
     INTEGER :: InterfaceSFNumNodes, InterfaceSFNumElems
! -- Number of Nodes and Elements  on the Solid/Fluid Non-Burning Interface
     INTEGER :: InterfaceSFnbNumNodes, InterfaceSFnbNumElems
! -- Number of Nodes and Elements  on the Non-Solid/Fluid Interface
     INTEGER :: InterfaceSNumNodes, InterfaceSNumElems

! -- Element connectivity array           
     INTEGER, POINTER, DIMENSION(:,:) :: InterfaceSFElemConn           
     INTEGER, POINTER, DIMENSION(:,:) :: InterfaceSFnbElemConn
     INTEGER, POINTER, DIMENSION(:,:) :: InterfaceSElemConn
! -- Nodal coordinates
     REAL*8 , POINTER, DIMENSION(:,:) :: InterfaceSFNodalCoors
     REAL*8 , POINTER, DIMENSION(:,:) :: InterfaceSFnbNodalCoors
     REAL*8 , POINTER, DIMENSION(:,:) :: InterfaceSNodalCoors
! -- Displacements (incremental)
     REAL*8, POINTER, DIMENSION(:,:) :: InterfaceSFNodalDisps
     REAL*8, POINTER, DIMENSION(:,:) :: InterfaceSFnbNodalDisps
! -- Displacements (total)
     REAL*8, POINTER, DIMENSION(:,:) :: InterfaceSFTotalNodalDisps
     REAL*8, POINTER, DIMENSION(:,:) :: InterfaceSFnbTotalNodalDisps
! -- Velocities
     REAL*8, POINTER, DIMENSION(:,:) :: InterfaceSFNodalVels
     REAL*8, POINTER, DIMENSION(:,:) :: InterfaceSFnbNodalVels
! -- Tractions
!!$ #OLD     REAL*8, POINTER, DIMENSION(:,:) :: InterfaceSFElemTract
     REAL*8, POINTER, DIMENSION(:) :: InterfaceSFElemTract
     REAL*8, POINTER, DIMENSION(:) :: InterfaceSFnbElemTract

! heat flux
     REAL*8, POINTER, DIMENSION(:) :: InterfaceSFHeatFlux
     REAL*8, POINTER, DIMENSION(:) :: InterfaceSFNodalTemp

! -- Vbar 
     REAL*8, POINTER, DIMENSION(:,:) :: InterfaceSFVbar
     REAL*8, POINTER, DIMENSION(:,:) :: InterfaceSFnbVbar
     REAL*8, POINTER, DIMENSION(:,:) :: InterfaceSVbar
!
! -- Element Integration (10 node)
! 0 - full integration
! 1 - reduced integration

     INTEGER :: iElIntgratn


!---- ** COHESIVE MATERIAL PROPERTIES **
! - Number of cohesive materials
     INTEGER :: NumMatCoh   
!--  normal characteristic lengths | tangent characteristic lengths
     REAL*8, POINTER, DIMENSION(:) :: deltan, deltat
!--  max. normal stresses |  max. shearing stresses | friction coeffs. 
     REAL*8, POINTER, DIMENSION(:) :: SigmaMax, TauMax
!--  initial Sthresholds      
     REAL*8, POINTER, DIMENSION(:) :: Sinit

!--  initial Sthresholds      
     REAL*8, POINTER, DIMENSION(:,:) :: Sthresh1,Sthresh2

! -- Solution Dynamic Variables

!--   displacement vector
     REAL*8, POINTER, DIMENSION(:) :: Disp


!--   Temperature vector
     REAL*8, POINTER, DIMENSION(:) :: Temperature

!--   previous time step displacement
     REAL*8, POINTER, DIMENSION(:) :: DispOld
!--   mesh displacement vector
     REAL*8, POINTER, DIMENSION(:) :: DispBar, DispTotal !d_bar, d_total
!--   velocity vector at t=t+dt/2 
     REAL*8, POINTER, DIMENSION(:) :: VeloHalf !vhalf
!--   mesh velocity vector at previous step
     REAL*8, POINTER, DIMENSION(:) :: VeloBarOld !v_bar_old
!--   mesh velocity vector
     REAL*8, POINTER, DIMENSION(:) :: VeloBar !v_bar
!--   velocity of boundary nodes
     REAL*8, POINTER, DIMENSION(:) :: VeloBndry !vb
!--   acceleration of boundary nodes
     REAL*8, POINTER, DIMENSION(:) :: AccelBndry !ab 
!--   mesh acceleration vector
     REAL*8, POINTER, DIMENSION(:) :: AccelBar !a_bar
!--   acceleration vector
     REAL*8, POINTER, DIMENSION(:) :: Accel !a 
!--   CST stress
     REAL*8, POINTER, DIMENSION(:,:) :: S11, S22, S33, S12, S23, S13 
!--   Von Mises stress
     REAL*8, POINTER, DIMENSION(:) :: SVonMises

! -- Dummy Variables
     
! 1-3 for 4 node tet (i.e. 3 node triangles)
! 4-6 for 10 node tet (i.e. 6 node triangles, mid-side nodes get traction)
     INTEGER :: LwrBnd,UppBnd

! Mass conservation
! single processors total mass for the volume
     REAL*8 :: TotalMassSolidp
! single processors total mass flux
     REAL*8 :: xmdot_totalp ! not used
! single processors total burning area
     REAL*8 :: areap ! not used

     REAL*8 :: TotalGeomVolp
     REAL*8 :: TotalGeomUndefVolp
! --- Node tracking
     INTEGER :: NumNodeIO, NumNodeIOpid
     INTEGER, POINTER, DIMENSION(:) :: NodeIO

! -- For StandAlone Rocfrac
     REAL*8 :: DummyTractVal
     REAL*8 :: DummyBurnRate
     REAL*8 :: DummyFlux

     LOGICAL :: IONEWER ! change the input format in 2.5

! -- Damping [C] matrix
     REAL*8 :: KappaDamp
     LOGICAL :: DampEnabled

! -- Nodal Elements Mass Lumping

     INTEGER :: NdMassLump
! -- To enforce traction where no fluid is present
!    good for stand-alone mode, or cases like the artery
     LOGICAL :: EnforceTractionS
     LOGICAL :: EnforceTractionSF
! -- Associates the surface mesh element with the volumitric element
     INTEGER, POINTER, DIMENSION(:) :: MapSFElVolEl, MapSFnbElVolEl, MapSElVolEl

     REAL*8, POINTER, DIMENSION(:,:) :: AmpTable
     INTEGER :: NumEntries 


!--   proportionality constant
     REAL*8 :: prop
!--   m%d(prop)/dt
     REAL*8 :: slope    

     INTEGER :: iAmpCnt

     LOGICAL :: NdBasedEl

     LOGICAL :: UnDefConfig

     INTEGER :: NumMatVolHT

     logical :: HeatTransSoln

     real*8 :: Temperature0 ! intial temperature of the model
     real*8 :: ThermalDiffusivity

     real*8, POINTER, DIMENSION(:) :: KappaHT, Cp

     REAL*8, POINTER, DIMENSION(:,:,:) :: mixed_map
     REAL*8, POINTER, DIMENSION(:,:,:)  :: enhanced_map
     REAL*8, POINTER, DIMENSION(:,:)  :: Aenh

     real*8, pointer, dimension(:,:,:) :: dmat
    
     logical :: ArtificialDamping
     REAL*8, POINTER, DIMENSION(:,:) :: DetF_old

     
     REAL*8, POINTER, DIMENSION(:) :: BCValueGlb

     INTEGER :: NumNdsBCcrypt
     INTEGER, POINTER, DIMENSION(:,:) :: BCFlagCrypt

     INTEGER :: NumProbesEl, NumProbesNd
     
     REAL*8, POINTER, DIMENSION(:,:) :: ProbeCoorNd, ProbeCoorEl
     INTEGER, POINTER, DIMENSION(:) :: ProbeNd

     LOGICAL, POINTER, DIMENSION(:) :: PointOnProc

     INTEGER :: NSTATEV
     REAL*8, POINTER, DIMENSION(:) :: STATEV_Part1
     
     REAL*8, POINTER, DIMENSION(:) :: STATEV_Part2
     
     INTEGER :: NMATRIX
     REAL*8, POINTER, DIMENSION(:) :: MATRIX

     integer :: NPARTICLE, NPARTICLETYPE 
     REAL*8, POINTER, DIMENSION(:,:) :: PARTICLE

     integer :: NINTERFAC
     REAL*8, POINTER, DIMENSION(:) :: INTERFAC
     
     
     LOGICAL :: DebondPart,DebondPart_Matous
     LOGICAL :: ThermalExpansion

     REAL*8, POINTER, DIMENSION(:) :: StrainTrace

     LOGICAL :: AmplitudeTable

     REAL*8 :: alpha1, alpha2, c3, p1, p2,Yin, a_eta, a_zeta, cm, cb, c2

     REAL*8, POINTER, DIMENSION(:,:,:) :: L_tensor, M_tensor
     REAL*8, DIMENSION(1:6,1:6) :: L_bar, M_bar, Lo
     REAL*8, POINTER, DIMENSION(:,:) :: StrainOld, SoftParam, cd
     INTEGER :: NumMatVol_Part

     REAL*8, DIMENSION(1:2) :: ShrMod, BulkMod, PoisRat

     LOGICAL :: debug_state

     INTEGER :: NumNpOverlay,NumElOverlay
     INTEGER, POINTER, DIMENSION(:,:) :: ElConnOverLay
     REAL*8, POINTER, DIMENSION(:,:) :: CoorOverlay
!     INTEGER, POINTER, DIMENSION(:) :: Map2VolNdOverlay

     REAL*8, POINTER, DIMENSION(:,:) :: etaOverlay, nuOverlay

     INTEGER, POINTER, DIMENSION(:) :: MapFaceEl2Vol1, FaceOfVolEL1

     INTEGER, POINTER, DIMENSION(:) :: MapFaceEl2Vol2, FaceOfVolEL2

     INTEGER :: nf1, nf2
     
     INTEGER:: nsubn1, nsubf1, nsubn2, nsubf2, nn, nf
     INTEGER,  POINTER, DIMENSION(:,:) :: sd_subfaces1,sd_subfaces2
     INTEGER,  POINTER, DIMENSION(:) :: sd_subface_parentsB
     REAL*8, POINTER, DIMENSION(:,:) :: sd_coor1,sd_coor2
     INTEGER,  POINTER, DIMENSION(:) :: sd_subface_parents1,sd_subface_parents2
     REAL*4, POINTER, DIMENSION(:,:) :: sd_subface_nat_coors1,sd_subface_nat_coors2
     INTEGER,  POINTER, DIMENSION(:) :: sd_subface_counterparts1,sd_subface_counterparts2
      
     INTEGER :: Verb

     LOGICAL :: OverlayExist

  END TYPE ROCFRAC_GLOBAL


CONTAINS

  FUNCTION getTimeString(time)
    REAL*8        :: time
    CHARACTER*9   :: getTimeString
    CHARACTER*15  :: ichrstring
    
    WRITE(ichrstring,'(e13.6)') time*1.e9

    getTimeString = ichrstring(11:12)//'.'//ichrstring(3:8)
  END FUNCTION getTimeString
  
  SUBROUTINE associate_pointer( attr, ptr)
    TYPE(ROCFRAC_GLOBAL), POINTER :: attr, ptr
    ptr => attr
  END SUBROUTINE ASSOCIATE_POINTER

END MODULE ROCSTAR_RocFrac

