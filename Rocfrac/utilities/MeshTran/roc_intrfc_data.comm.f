      SUBROUTINE rocface_data(lmtri,nprocs,numnp_old,NumElSurf,epart,
     $     coor,prefix,jpref,MapSurfEl2VolEl,ik13d,numelv_prmry,numvertx,numel_2d,elm_2D_flag,
     $     jk1, maxdim,NumNp3D,jk13d,maxdim3d,iNdsBurnFlg,NumBorderVol,NumMat,MapGlbEl2LocEl,ElOnPartBndry,MatId)

      USE linked_list_2d

      IMPLICIT NONE

      CHARACTER*20 :: prefix
      INTEGER jpref
      integer :: maxdim,maxdim3d
      integer, dimension(1:numnp_old) :: iNdsBurnFlg
      INTEGER :: numelv_prmry
      INTEGER, DIMENSION(1:numnp_old*3,1:nprocs) :: ik13d
      INTEGER, DIMENSION(1:nprocs) :: NumNp3D

C-- number of original node points & elements
      INTEGER :: numnp_old, NumElSurf
C-- number of nodes for triangle
      INTEGER :: numvertx
C-- dummy variables
      INTEGER :: i,ii,j,jj,k,kk
      INTEGER :: iaux,iaux1,iaux2,iaux3,iaux4,iaux5,iaux6,n,mm

      INTEGER :: nk
      INTEGER :: flag,ios

      INTEGER :: nn,nprocs,iunit

C-- Coordinate array from triangle
      REAL*8 :: coor(1:3,1:numnp_old)
C-- Connectivity array from triangle for linear strain triangles
      INTEGER, DIMENSION(1:numvertx,1:NumElSurf) :: lmtri
C-- Connectivity array (OUTPUT) for cohesive elements
C-- Partioned element's processor (from METIS)
      INTEGER, DIMENSION(1:NumElSurf) :: epart
C-- Number of nodal points & elements on each processor
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  numel
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NumNp2D
C-- For the lst: relates the old node numbering to the new nodes
C----- dimension:  <1>  old node number   <2> processor
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: ik1
      INTEGER, DIMENSION(1:maxdim,1:nprocs,1:2) :: jk1
      INTEGER, DIMENSION(1:maxdim3d,1:nprocs) :: jk13d
      CHARACTER*4 :: ai3

      INTEGER :: ip
      TYPE(vol_info_type), POINTER :: vol_item
C --  List of volumetric elements
      TYPE(vol_list_type), TARGET, ALLOCATABLE, DIMENSION(:,:) :: vol_list
! global Tetrahedral element element
      INTEGER, DIMENSION(1:NumElSurf) :: MapSurfEl2VolEl
      INTEGER, DIMENSION(1:nprocs,1:2) :: numel_2d
      INTEGER, DIMENSION(1:NumElSurf) :: elm_2D_flag
      INTEGER :: iflag,Ndglb

      INTEGER :: LocVolElId, VolElId

      INTEGER :: NumMat
      INTEGER, DIMENSION(1:numelv_prmry) :: MapGlbEl2LocEl
      LOGICAL, DIMENSION(1:numelv_prmry) :: ElOnPartBndry
      INTEGER, DIMENSION(1:numelv_prmry) :: MatId
      integer, DIMENSION(1:NumMat,1:nprocs) :: NumBorderVol

      integer :: imat

      ALLOCATE(NumNp2D(1:nprocs,1:2))
      NumNp2D(1:nprocs,1:2) = 0
      ALLOCATE(vol_item)
      ALLOCATE(vol_list(1:nprocs,1:2))
      ALLOCATE(numel(1:nprocs,1:2))
      numel(1:nprocs,1:2) = 0
      ALLOCATE(ik1(1:numnp_old,1:nprocs,1:2))
      ik1 = 0

      DO i = 1, NumElSurf
         ip = epart(i)
         iflag = elm_2D_flag(i)

! Global Volume element id associated with surface
         
         VolElId = MapSurfEl2VolEl(i)

! Material Type associated with Volume element

         imat = MatId(VolElId)
!
! flag =1 S/F interface
! flag = 2 S interface
         numel(ip,iflag) = numel(ip,iflag) + 1
         DO k = 1, numvertx
            nk = lmtri(k,i)
            IF(ik1(nk,ip,iflag).EQ.0)THEN
               
               NumNp2D(ip,iflag) = NumNp2D(ip,iflag) + 1
               jk1(NumNp2D(ip,iflag), ip, iflag) = nk
               ik1(nk,ip,iflag) = NumNp2D(ip,iflag)
            ENDIF
            vol_item%lmvol(k) = ik1(nk,ip,iflag)
         ENDDO
         
         vol_item%mat_vol = 1 ! fix matvol(i)

         IF( ElOnPartBndry( VolElId )  )THEN
            LocVolElId = NumBorderVol(imat,ip) - MapGlbEl2LocEl(VolElId) + 1
         ELSE
            LocVolElId = MapGlbEl2LocEl(VolElId) + NumBorderVol(imat,ip)
         ENDIF

         vol_item%mapelem =  LocVolElId
         vol_item%flag = elm_2D_flag(i) ! not correct

C     Add item to volumetric element list
         CALL vol_insert_tail(vol_list(ip,iflag),vol_item)
      ENDDO
C     WRITE THE COORDINATES
      iunit = 4000
      
      OPEN(iunit,FILE=prefix(1:jpref)//'/fracSF.im',STATUS='replace',FORM='formatted')
      OPEN(iunit+1,FILE=prefix(1:jpref)//'/fracS.im',STATUS='replace',FORM='formatted')

      WRITE(iunit,*) nprocs,numvertx
      WRITE(iunit+1,*) nprocs,numvertx

      DO i = 1, nprocs
         WRITE(iunit,*) i
         WRITE(iunit+1,*) i
         DO k = 1, 2

            IF(k.EQ.1)THEN
               DO kk = 1, NumNp3D(i)
                  Ndglb = jk13d(kk,i) ! unpartioned node number (tet)
                  IF(iNdsBurnFlg(Ndglb).EQ.1)THEN
                     IF(ik1(Ndglb,i,1).EQ.0.AND.ik1(Ndglb,i,2).EQ.0)THEN
                        NumNp2D(i,k) = NumNp2D(i,k) + 1
                        print*,'Processor',i,'Node',NumNp2D(i,k),Ndglb
                        jk1(NumNp2D(i,k),i,k) = Ndglb ! original unpartioned node number
                        ik1(Ndglb,i,k) = NumNp2D(i,k)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

!
! First surface list is S/F interface
! Second surface list is S interface

            WRITE(iunit+k-1,'(2i10)') NumNp2D(i,k), numel(i,k)

            IF(k.EQ.1)THEN
               IF(NumNp2D(i,k).GT.0.AND.NumNp2D(i,k).LE.2.AND.numel(i,k).EQ.0)THEN
                  PRINT*,'*********'
                  print*,'POSSIBLE ERROR: SURFACE MESH HAS ONLY NODES, NO ELEMENTS'
                  PRINT*,'**********'
                  PRINT*,'PROCESSOR file = ', i-1
               ENDIF
            ENDIF
                  
            DO j = 1, NumNp2D(i,k) ! jk1 maps the triangle surface nodes to tetrahedral nodes
!               print*,i,j,jk1(j,i,k),ik13d( jk1(j,i,k),i )
               WRITE(iunit+k-1,'(3(x,e16.9),1i10)') coor(1:3,jk1(j,i,k)),ik13d( jk1(j,i,k),i )
            ENDDO
            CALL print_vol_list(vol_list(i,k),numvertx,iunit+k-1)
         ENDDO

      ENDDO

      CLOSE(iunit)
      CLOSE(iunit+1)

      DEALLOCATE(ik1)
      RETURN
  
      END
