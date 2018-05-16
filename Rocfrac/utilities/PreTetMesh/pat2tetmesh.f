      PROGRAM pat2tetmesh

      USE meshdata

      IMPLICIT NONE

      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: refnum

      CALL read_deck()

      CALL read_patran()

! edge flags, not used

      ALLOCATE(refnum(1:3,1:numelv))
      refnum(:,:) = 0

      CALL TetMeshOut(numelv,lmelv,refnum,npress,numnp,coor,ibcflg,
     $     prefx,prefx_lngth)

      DEALLOCATE(lmelv,refnum,coor,ibcflg,npress)

      END
