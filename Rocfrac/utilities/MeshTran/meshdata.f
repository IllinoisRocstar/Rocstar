      MODULE meshdata

      IMPLICIT NONE

      INTEGER :: io_input       ! unit id of control deck
      PARAMETER(io_input = 10)

!  - DECK FILE PRAMETERS
        CHARACTER*20 :: prefx   ! I/O file prefx
        INTEGER :: prefx_lngth  ! length of prefx
!  - Fastest Dilatitial wave speed
        REAL*8 :: cd_fastest
!  - Derive Type for the boundary conditions
      TYPE bcvalues
          INTEGER :: b1, b2, b3
          REAL*8 :: bc1, bc2, bc3
      END TYPE
!  - Stores the boundary conditions data values
      TYPE(bcvalues), DIMENSION(32) :: bc_conditions

! PRIMARY MESH DATA (i.e. mesh created with ansys, patran, or mesh3d):

!   - Number of nodes 
      INTEGER :: numnp_prmry
!   - Number of volumetric elements  
      INTEGER :: numelv_prmry
!   - meshing software flag
      INTEGER :: iansys         ! 0- no, 1-yes
      INTEGER :: ipatran        ! 0- no, 1-yes
      INTEGER :: itetmesh       ! 0- no, 1-yes
      INTEGER :: ipatcohin         ! 0- no, 1-yes
      INTEGER :: itetcohin      ! 0- no, 1-yes

      INTEGER :: numvertx

      REAL*8 :: ConvertUnit
      

!   - Coordinate array
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: coor
!   - Coordinate array
      REAL*8, ALLOCATABLE, DIMENSION(:) :: press_nodal
!   - Connectivity array for volumetric elements
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lmelv_prmry
!   - Element neighbor array
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: neigh
!   - Nodal boundary flag
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ibcflg
!   - Nodal mesh motion boundary flag
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ibcflg_mm
!   - Nodal corner boundary flag
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ibcaxi
!   - Partioned element's processor (from METIS)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: epart

!   - Number of surface triangles with applied tractions
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: numel_2d
!   - Number of surface nodes with applied mesh motion velocity
      INTEGER, ALLOCATABLE, DIMENSION(:) :: numnp_2d
! sub
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: lmtri
      INTEGER, DIMENSION(:), ALLOCATABLE :: elm_2D
      INTEGER, DIMENSION(:), ALLOCATABLE :: elm_2D_flag
      INTEGER :: numel_tri
      INTEGER :: numnp_tri
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: neigh_2d
      INTEGER, DIMENSION(:), ALLOCATABLE :: epart_2d
! sub

!   - Node tracking
        INTEGER :: NumNodeIO
        INTEGER, ALLOCATABLE, DIMENSION(:) :: NodeIO

!   - Node 

        INTEGER :: IOformat     ! 0 = unformatted (default), 1 = formatted 
        ! 0 = no output, 1 = output
        INTEGER :: iopmvis  ! PMVIS software to view the partition 0=no, 1 = yes
        INTEGER :: ipress  ! flag for if we have pressure loading 0=no, 1 = yes

        integer :: IntFaceFlag
        integer :: itype6,itype7,itype8

!
! -- derived surface type
!
c$$$        TYPE request
c$$$            INTEGER, DIMENSION(1:6) :: Conn ! list of nodes
c$$$            INTEGER :: NumEl
c$$$            INTEGER :: NumNp
c$$$            TYPE(request), POINTER :: next
c$$$        END TYPE request
c$$$      
c$$$        TYPE(request), POINTER :: head, tail
c$$$        TYPE(request), POINTER :: item, first

        integer, allocatable, dimension(:) :: iNdsBurnFlg

        integer :: NumMat

      INTEGER, DIMENSION(:), ALLOCATABLE :: MatId

      ! limit of 10 nodes to monitor history

      INTEGER, dimension(1:10) :: NdHistory
      integer :: NumNdHistory
        
      END MODULE meshdata 
