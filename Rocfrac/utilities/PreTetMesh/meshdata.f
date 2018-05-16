      MODULE meshdata

      IMPLICIT NONE

      INTEGER :: io_input       ! unit id of control deck
      PARAMETER(io_input = 10)

!  - DECK FILE PRAMETERS
      CHARACTER*20 :: prefx     ! I/O file prefx
      INTEGER :: prefx_lngth    ! length of prefx

! PRIMARY MESH DATA (i.e. 2D mesh datat created with patran)

!   - Number of nodes 
      INTEGER :: numnp
!   - Number of triangle elements  
      INTEGER :: numelv
!   - Coordinate array
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: coor
!   - Connectivity array for volumetric elements
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lmelv
!   - Face Flag
      INTEGER, ALLOCATABLE, DIMENSION(:) :: npress
!   - Nodal Flag
      INTEGER, ALLOCATABLE, DIMENSION(:) ::ibcflg


      END MODULE meshdata
