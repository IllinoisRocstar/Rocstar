      PROGRAM meshtran3D

      USE linked_list

      USE meshdata

      IMPLICIT NONE

!     Type definition 

      TYPE nodebc
          INTEGER :: node,bc
      END TYPE

      TYPE nodebc_lst
          INTEGER :: node_lst,bc_lst
      END TYPE

      TYPE sndrcv_buf
          INTEGER :: sndn,rcvn,nodes
      END TYPE

      TYPE(sndrcv_buf), ALLOCATABLE, DIMENSION(:) :: sndrcvnod
      TYPE(sndrcv_buf), ALLOCATABLE, DIMENSION(:) :: sndrcvnod_lst

!-- number of nodes for triangle
      INTEGER :: numvertx2d
!-- dummy variables
      INTEGER :: i,ii,j,jj,k,kk
      INTEGER :: iaux,iaux1,iaux2,iaux3,n,mm
      REAL*8 :: aux
      INTEGER :: ntime

      INTEGER :: nk
C-- Stores the material number associated with lst
      INTEGER, ALLOCATABLE, DIMENSION(:) :: matvol,iaux89
C-- Stores the cohesive material number associated with cohesive
      INTEGER :: flag

      INTEGER :: numbcss
      INTEGER :: indx,jndx
      INTEGER :: edgecut,nn,nprocs,iunit
      CHARACTER*3 :: ai4
      CHARACTER*4 :: ichr4,ai1
      INTEGER :: i1,i2,i3


       INTEGER, ALLOCATABLE, DIMENSION(:) :: ncoor

C-- Element Neighbor array from triangle
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lmtri_2d
C-- Connectivity array (OUTPUT) for cohesive elements
C----- dimension:  <1>  6 nodes   <2> local element id  <3> processor
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: lmcoh
C-- Tempory holding array of cohesive and no-cohesive (input to METIS)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: elmnts
C-- Tempory holding array of partitioned 'elmnts' array
      INTEGER, ALLOCATABLE, DIMENSION(:) :: epart_p
C-- Number of nodal points & elements on each processor
      INTEGER, ALLOCATABLE, DIMENSION(:) :: numnp
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: numel
C-- Number of cohesive elements on each processor
      INTEGER, ALLOCATABLE, DIMENSION(:) :: numclst
C-- For the lst: relates the old node numbering to the new nodes
C----- dimension:  <1>  old node number   <2> processor
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ik1
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jk1
C-- Number of nodal b.c. for each processor 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: numbc
C-- Number of nodal mesh motionb.c. for each processor 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: numbc_mm
C-- Number of cohesive elements on processor boundary for each processor
      INTEGER, ALLOCATABLE, DIMENSION(:) :: num_border_coh
C-- Number of neighboring processors for R_co
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nproc_neigh
C-- Number of neighboring processors for R_in
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nproc_neigh_lst
C-- Number of neighboring processors for R_co on No/cohesive boundary
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nproc_neigh_rco
C-- Stores the neighboring processors that I'm communcating with.
      INTEGER, ALLOCATABLE, DIMENSION(:) :: my_neigh
C-- Stores if node is on the processor boundary
      INTEGER, ALLOCATABLE, DIMENSION(:) :: bord_node_org
C-- Stores the node for R_co communication
      INTEGER, ALLOCATABLE, DIMENSION(:) :: bord_node_Rco
C--
      INTEGER :: err
      INTEGER, ALLOCATABLE, DIMENSION(:) :: matclst,lmtemp
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nnntemp_lst
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nnntemp_rco
    
      TYPE(nodebc), ALLOCATABLE, DIMENSION(:,:) :: ibc
      TYPE(nodebc), ALLOCATABLE, DIMENSION(:,:) :: ibc_mm
      TYPE(nodebc_lst), ALLOCATABLE, DIMENSION(:,:) :: ibc_lst
      TYPE proclist
         INTEGER, DIMENSION(:), POINTER :: proc_list
      END TYPE
      TYPE proclistrco
         INTEGER, DIMENSION(:), POINTER :: proc_list_rco
      END TYPE
      TYPE(proclistrco), ALLOCATABLE, DIMENSION(:) :: node_rco
      TYPE(proclist), ALLOCATABLE, DIMENSION(:) :: node
C-- Number of outer boundary nodes
      INTEGER :: ibcglobal
      INTEGER, DIMENSION(4) :: iboundx,iboundy
      REAL*8, DIMENSION(4) :: boundx,boundy
C-- wave properties and stability varibles
      REAL*8 :: xx,yy,zz,size1,size2,size3,size4,size5,size6
      REAL*8 :: dhmin,dt_courant
      REAL*8 :: dh_courant,cd_courant
      integer :: maxdim
C-

      INTEGER :: numnp_total,nface

      INTEGER, allocatable,DIMENSION(:,:) :: NumNp2D
      INTEGER, allocatable,DIMENSION(:,:,:) :: jk1_2D
      INTEGER, allocatable,DIMENSION(:,:,:) :: ik1_2D

      CHARACTER*3 :: ai3
      character*1 ichr1
      character*12 ichr12
      INTEGER :: num_zones,numcoh_zone,numelv_prmry_zone
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: izones
      INTEGER, ALLOCATABLE, DIMENSION(:) :: matzone,iizones,nbcele
      INTEGER, ALLOCATABLE, DIMENSION(:) :: numelv_prmry_zones, ik1_z
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MapGlbEl2LocEl

      INTEGER :: numel_z,numnp_z,p2max_lst,num_rco_border,ip,icount

      INTEGER iaux23
      INTEGER :: num_vol,num_coh
      INTEGER :: iflag
      INTEGER, DIMENSION(1:4) :: iflag89
      INTEGER :: ntri
      INTEGER :: ios
      INTEGER :: jk1_size

      TYPE(coh_info_type), POINTER :: coh_item
      TYPE(coh_list_type), TARGET, ALLOCATABLE, DIMENSION(:) :: coh_list
      TYPE(coh_comm_info_type), POINTER :: coh_comm_item
      TYPE(coh_comm_list_type), TARGET, ALLOCATABLE, DIMENSION(:,:) :: coh_comm_list
      TYPE(vol_info_type), POINTER :: vol_item
      TYPE(vol_list_type), TARGET, ALLOCATABLE, DIMENSION(:,:) :: vol_list
      TYPE(vol_box_point), ALLOCATABLE, DIMENSION(:) :: ik1_c4
      TYPE(comm_list_type), TARGET, ALLOCATABLE, DIMENSION(:,:) :: ID_sendto

      INTEGER ni,node1,node2,node3,node4,gnode1,gnode2,gnode3,gnode4

      TYPE(bcvalues), DIMENSION(16) :: bc_mshmtn

      INTEGER, ALLOCATABLE, DIMENSION(:) :: iaux_p

      INTEGER :: numploadelem
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: idpressload
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ipressflag
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: pressload

      LOGICAL, allocatable, dimension(:) :: ElOnPartBndry

      integer, allocatable, dimension(:,:,:) :: icnt
      integer, allocatable, dimension(:,:) :: NumBorderVol

      INTEGER :: myid
      REAL*8 :: shift
      INTEGER :: iargc
      CHARACTER*3 chr_arg
      CHARACTER*4 chr_procs
      CHARACTER*20 chr_units

      integer :: imat, itmp1, itmp2

      integer, dimension(1:10) :: NdHistoryFlag
      integer, allocatable, dimension(:) :: NumNdHistoryP
      integer, allocatable, dimension(:,:) :: NdHistoryP
      ConvertUnit = 1.d0

      n = iargc()
      IF(n.LE.1) GOTO 554

      j = 1
      DO i = 1, n/2
         CALL getarg( j, chr_arg)
         IF(chr_arg.EQ.'-np')THEN
            CALL getarg( j+1, chr_procs)
            READ(chr_procs,*) nprocs
            j = j +2
         ELSE IF(chr_arg.EQ.'-un')THEN ! units
            CALL getarg( j+1, chr_units)
            READ(chr_units,*) ConvertUnit
            j = j + 2
            PRINT*,'UNITS CONVERSION'
            PRINT*,' multiply by ', ConvertUnit
            
         ELSE
            GOTO 554
         ENDIF
      ENDDO

      GOTO 555

 554  CONTINUE
      PRINT*,'Usage:'
      PRINT*,'   MeshTran -np #'
      PRINT*,'         - where # is the number of processors'
      STOP

 555  CONTINUE

      myid = 0

      CALL readinp(ntime)

      PRINT*,'numvertx=',numvertx

c      CALL system('\rm -f -r '//prefx(1:prefx_lngth))

      num_zones = 1
      ALLOCATE(izones(1:2,1:num_zones))

c      DO i=1,num_zones
c         READ(3,*) izones(1,i),izones(2,i)
c      ENDDO
c      READ(3,'()')
C -- mesh motion boundary conditions
cc      DO i = 1, 16
cc         READ(3,*) iaux, bc_mshmtn(i)%b1, bc_mshmtn(i)%b2,bc_mshmtn(i)%b3,
cc     $        bc_mshmtn(i)%bc1, bc_mshmtn(i)%bc2,bc_mshmtn(i)%bc3
cc      ENDDO
cc      CLOSE(3)
            
      dhmin = 1000000000.d0

      IF(iansys.EQ.1)THEN
         PRINT*,'ANSYS NO LONGER SUPPORTED'
         PRINT*,'STOPPING'
         STOP
!        CALL read_ansys(numbcss, dhmin)
      ELSE IF(ipatran.EQ.1)THEN
         CALL read_patran(numbcss,numvertx2d,dhmin,nprocs)
      ELSE IF(ipatcohin.EQ.1)THEN
         CALL read_patran_cohin(numbcss,numvertx2d,dhmin,nprocs)
      ELSE IF(itetcohin.EQ.1)THEN
         CALL read_tetmesh_cohin(numbcss,numvertx2d,dhmin,nprocs)
      ELSE IF(itetmesh.EQ.1)THEN
!NOTE:
!IF both ascii and binary temesh input files exist,then
! the binary file will be the one that gets read.
         CALL read_tetmesh(numbcss,numvertx2d,dhmin,nprocs)
      ENDIF
c
c -- Checking Courant condition for time step
c
      dt_courant = dhmin/cd_fastest
      cd_courant = cd_fastest
      dh_courant = dhmin

      PRINT*,' --------------------------------------'
      PRINT*,' ---- COURANT STABILITY CONDITION -----'
      PRINT*,' --------------------------------------'
      PRINT*,'   dt_courant =',dt_courant
      PRINT*,'   cd_courant =',cd_courant
      PRINT*,'    h_courant =',dh_courant
      PRINT*,' --------------------------------------'

c      READ(19,*) numploadelem
c
      ALLOCATE(ipressflag(1:numelv_prmry))
      ipressflag(:) = 0
c     ALLOCATE(idpressload(1:2,1:numploadelem), pressload(1,1:numelv_prmry))
c     DO j=1,numploadelem
c     READ(19,*) iaux, ipressflag(iaux),pressload(1,iaux)
c     ENDDO

      ibcglobal = 0
      DO i = 1,numnp_prmry
ccccc    READ(11,*) ibcflg(i)
         IF(ibcflg(i).NE.0.AND.ibcflg(i).LT.1000) 
     $        ibcglobal = ibcglobal+1
      ENDDO
      ALLOCATE(ibc(1:30*ibcglobal,1:nprocs))
      ALLOCATE(ibc_mm(1:30*ibcglobal,1:nprocs))

      ALLOCATE(numelv_prmry_zones(1:2))
      ALLOCATE(iizones(1:numelv_prmry))
      ALLOCATE(matzone(1:numelv_prmry))
      ALLOCATE(matvol(1:numelv_prmry))

      numcoh_zone = 0
      numelv_prmry_zones(:) = 0
      ! izones is not used anymore
      izones(2,:) = 1
      izones(1,:) = 0 ! 1 ! 0 !0 !1 !!! 0  !!! 1 !!!!! 1 !!!!!!!!! 0
      iizones(:) =  0 ! 1 ! 0 ! 0 ! 1 ! 0 ! 1 ! 0 !1 !0 ! fix

      matvol(:) = 1 ! fix if want different materials. could remove
      matzone(:) = 1 ! currently not working, could remove
c
c-- LST
      ALLOCATE(bord_node_org(1:numnp_prmry))
      ALLOCATE(node(1:numnp_prmry))


      ALLOCATE(ibc_lst(1:2*ibcglobal,1:nprocs))
      ALLOCATE(bord_node_Rco(1:numnp_prmry))
      bord_node_org = 0
      bord_node_Rco = 0
      iaux = 1
      ALLOCATE(numel(1:NumMat,1:nprocs))
      numel(:,:) = 0
      ALLOCATE(nbcele(1:nprocs))
      nbcele(:) = 0

      DO n = 1, numelv_prmry          ! Loop over the tetrahedra elements
         imat = MatId(n)
         numel(imat,epart(n)) = numel(imat,epart(n)) + 1 
         DO i = 1, 4            ! Loop over the faces of the tetrahedral
C     Check to see that the side is not on the outside boundary
            IF(neigh(i,n).GT.0)THEN
C     test to see if the side is located between partitions and make sure 
C     it is a partion between similar elements
               IF(epart(neigh(i,n)).NE.epart(n).AND.
     $              iizones(n).EQ.iizones(neigh(i,n)))THEN
ccc               IF(epart(neigh(i,n)).NE.epart(n).AND.
ccc     $              iizones(n).EQ.0)THEN
                  IF(i.EQ.1)THEN
c     assign a value of one to nodes that are on the border
                     bord_node_org(lmelv_prmry(2,n)) = 1
                     bord_node_org(lmelv_prmry(3,n)) = 1
                     bord_node_org(lmelv_prmry(4,n)) = 1
C     create arrays that store what processor that node belongs to
                     ALLOCATE(node(lmelv_prmry(2,n))%proc_list(1:nprocs))
                     node(lmelv_prmry(2,n))%proc_list(1:nprocs) = -1
                     ALLOCATE(node(lmelv_prmry(3,n))%proc_list(1:nprocs))
                     node(lmelv_prmry(3,n))%proc_list(1:nprocs) = -1
                     ALLOCATE(node(lmelv_prmry(4,n))%proc_list(1:nprocs))
                     node(lmelv_prmry(4,n))%proc_list(1:nprocs) = -1
                     IF(numvertx.GT.4)THEN
                        bord_node_org(lmelv_prmry(6,n)) = 1
                        bord_node_org(lmelv_prmry(9,n)) = 1
                        bord_node_org(lmelv_prmry(10,n)) = 1
                        ALLOCATE(node(lmelv_prmry(6,n))%proc_list(1:nprocs))
                        node(lmelv_prmry(6,n))%proc_list(1:nprocs) = -1
                        ALLOCATE(node(lmelv_prmry(9,n))%proc_list(1:nprocs))
                        node(lmelv_prmry(9,n))%proc_list(1:nprocs) = -1
                        ALLOCATE(node(lmelv_prmry(10,n))%proc_list(1:nprocs))
                        node(lmelv_prmry(10,n))%proc_list(1:nprocs) = -1
                     ENDIF
                  ELSE IF(i.EQ.2)THEN
                     bord_node_org(lmelv_prmry(3,n)) = 1
                     bord_node_org(lmelv_prmry(1,n)) = 1
                     bord_node_org(lmelv_prmry(4,n)) = 1
                     ALLOCATE(node(lmelv_prmry(1,n))%proc_list(1:nprocs))
                     node(lmelv_prmry(1,n))%proc_list(1:nprocs) = -1
                     ALLOCATE(node(lmelv_prmry(3,n))%proc_list(1:nprocs))
                     node(lmelv_prmry(3,n))%proc_list(1:nprocs) = -1
                     ALLOCATE(node(lmelv_prmry(4,n))%proc_list(1:nprocs))
                     node(lmelv_prmry(4,n))%proc_list(1:nprocs) = -1
                     IF(numvertx.GT.4)THEN
                        bord_node_org(lmelv_prmry(7,n)) = 1
                        bord_node_org(lmelv_prmry(8,n)) = 1
                        bord_node_org(lmelv_prmry(10,n)) = 1
                        ALLOCATE(node(lmelv_prmry(7,n))%proc_list(1:nprocs))
                        node(lmelv_prmry(7,n))%proc_list(1:nprocs) = -1
                        ALLOCATE(node(lmelv_prmry(8,n))%proc_list(1:nprocs))
                        node(lmelv_prmry(8,n))%proc_list(1:nprocs) = -1
                        ALLOCATE(node(lmelv_prmry(10,n))%proc_list(1:nprocs))
                        node(lmelv_prmry(10,n))%proc_list(1:nprocs) = -1
                     ENDIF
                  ELSE IF(i.EQ.3)THEN
                     bord_node_org(lmelv_prmry(1,n)) = 1
                     bord_node_org(lmelv_prmry(2,n)) = 1
                     bord_node_org(lmelv_prmry(4,n)) = 1
                     ALLOCATE(node(lmelv_prmry(1,n))%proc_list(1:nprocs))
                     node(lmelv_prmry(1,n))%proc_list(1:nprocs) = -1
                     ALLOCATE(node(lmelv_prmry(2,n))%proc_list(1:nprocs))
                     node(lmelv_prmry(2,n))%proc_list(1:nprocs) = -1
                     ALLOCATE(node(lmelv_prmry(4,n))%proc_list(1:nprocs))
                     node(lmelv_prmry(4,n))%proc_list(1:nprocs) = -1
                    IF(numvertx.GT.4)THEN
                        bord_node_org(lmelv_prmry(5,n)) = 1
                        bord_node_org(lmelv_prmry(8,n)) = 1
                        bord_node_org(lmelv_prmry(9,n)) = 1
                        
                        ALLOCATE(node(lmelv_prmry(5,n))%proc_list(1:nprocs))
                        node(lmelv_prmry(5,n))%proc_list(1:nprocs) = -1
                        ALLOCATE(node(lmelv_prmry(8,n))%proc_list(1:nprocs))
                        node(lmelv_prmry(8,n))%proc_list(1:nprocs) = -1
                        ALLOCATE(node(lmelv_prmry(9,n))%proc_list(1:nprocs))
                        node(lmelv_prmry(9,n))%proc_list(1:nprocs) = -1
                     ENDIF
                  ELSE IF(i.EQ.4)THEN
                     bord_node_org(lmelv_prmry(1,n)) = 1
                     bord_node_org(lmelv_prmry(2,n)) = 1
                     bord_node_org(lmelv_prmry(3,n)) = 1
                     ALLOCATE(node(lmelv_prmry(1,n))%proc_list(1:nprocs))
                     node(lmelv_prmry(1,n))%proc_list(1:nprocs) = -1
                     ALLOCATE(node(lmelv_prmry(2,n))%proc_list(1:nprocs))
                     node(lmelv_prmry(2,n))%proc_list(1:nprocs) = -1
                     ALLOCATE(node(lmelv_prmry(3,n))%proc_list(1:nprocs))
                     node(lmelv_prmry(3,n))%proc_list(1:nprocs) = -1
                     IF(numvertx.GT.4)THEN
                        bord_node_org(lmelv_prmry(5,n)) = 1
                        bord_node_org(lmelv_prmry(6,n)) = 1
                        bord_node_org(lmelv_prmry(7,n)) = 1
                        ALLOCATE(node(lmelv_prmry(5,n))%proc_list(1:nprocs))
                        node(lmelv_prmry(5,n))%proc_list(1:nprocs) = -1
                        ALLOCATE(node(lmelv_prmry(6,n))%proc_list(1:nprocs))
                        node(lmelv_prmry(6,n))%proc_list(1:nprocs) = -1
                        ALLOCATE(node(lmelv_prmry(7,n))%proc_list(1:nprocs))
                        node(lmelv_prmry(7,n))%proc_list(1:nprocs) = -1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      num_rco_border = iaux - 1

      ALLOCATE(node_rco(1:num_rco_border))

      DO i=1,num_rco_border
         ALLOCATE(node_rco(i)%proc_list_rco(1:nprocs))
      ENDDO

      ALLOCATE(numnp(1:nprocs),numbc(1:nprocs),numbc_mm(1:nprocs))

      ALLOCATE(ik1_c4(1:numelv_prmry))
      DO ii = 1, NumMat
         jk1_size = MAXVAL(numel(ii,1:nprocs))*numvertx
      ENDDO
      ALLOCATE(jk1(1:jk1_size,1:nprocs))
      ALLOCATE(ik1(1:numnp_prmry*3,1:nprocs))
      ik1(1:numnp_prmry*3,1:nprocs) = 0
      jk1 = 0
      ii  = 0

      numnp(1:nprocs) = 0
      numel(1:NumMat,1:nprocs) = 0
      numbc(1:nprocs) = 0
      numbc_mm(1:nprocs) = 0
      ALLOCATE(vol_item)
      ALLOCATE(vol_list(1:NumMat,1:nprocs))
C
C -- Initialize link list variables

      vol_list(1:NumMat,1:nprocs)%num_border_vol = 0

      ALLOCATE(iaux89(1:nprocs))
      iaux89 = 0 ! fix no longer correct if iopmvis = 0
C
C     Open output file for Analysis Code for boundary information
C
      IF(iopmvis.EQ.1)THEN
      IF(IOformat.EQ.0)THEN
         DO i = 1, nprocs
            WRITE(ai1,'(i4.4)') i-1
            iunit = 2000 + i-1
            OPEN(iunit,FILE=prefx(1:prefx_lngth)//'/'//prefx(1:prefx_lngth)//'.'//ai1//'.bndry',
     $           STATUS='replace',FORM='unformatted', POSITION = 'append')
         ENDDO
      ELSE
         DO i = 1, nprocs
            WRITE(ai1,'(i4.4)') i-1
            iunit = 2000 + i-1
            OPEN(iunit,FILE=prefx(1:prefx_lngth)//'/'//prefx(1:prefx_lngth)//'.'//ai1//'.bndry',
     $           STATUS='replace',FORM='formatted', POSITION = 'append')
         ENDDO
      ENDIF
      ENDIF
C
C     Open output file for POV-RAY to display the boundary
C
      IF(iopmvis.EQ.1)THEN
         OPEN(88,FILE=prefx(1:prefx_lngth)//'/'//prefx(1:prefx_lngth)//'.inc',
     $        STATUS='replace',FORM='formatted')
         WRITE(88,'(a6)') 'mesh {'
      ENDIF
C
C -- RENUMBER THE VOLUMETRIC ELEMENTS USING LOCAL NUMBERING
      PRINT*,'RENUMBERING VOLUMETRIC ELEMENTS'

      allocate(MapGlbEl2LocEl(1:numelv_prmry))
      allocate(icnt(1:NumMat,1:2,1:nprocs))
      icnt(:,:,:) = 0
      allocate(ElOnPartBndry(1:numelv_prmry))
      DO i = 1, numelv_prmry
         ip = epart(i)
         imat = MatId(i)
         IF(iizones(i).EQ.1)THEN ! ELEMENT IN COHESIVE ADD ZONE
            ii = numnp(ip) + 1
            iflag = 0
            iflag89(1:numvertx) = 0
            DO k = 1,numvertx
C
C--   Number the lst according; assuming cohesive elements are present
               nk = lmelv_prmry(k,i)
               vol_item%lmvol(k) = ii + k - 1
               vol_item%mat_vol = matvol(i)
               jk1(ii + k - 1,ip) = nk
               ik1(nk,ip) = ii + k - 1 ! numnp(ip)
c     
c --  Transfer boundary conditions to local nodes 
c
               IF(ibcflg(lmelv_prmry(k,i)).GT.0
     $              .AND.ibcflg(lmelv_prmry(k,i)).LT.1000)THEN
                  numbc(ip) = numbc(ip) + 1
                  ibc(numbc(ip),ip)%node = ii + k - 1
                  ibc(numbc(ip),ip)%bc = ibcflg(lmelv_prmry(k,i))
               ENDIF
               IF(neigh(k,i).LE.0)THEN
                  iflag89(k) = k
               ENDIF
            ENDDO

            iunit = 2000 + ip - 1
            WRITE(ai4,'(i3.3)') ip
            DO k = 1, 4
               shift = .05 
               IF(iflag89(k).EQ.1)THEN
                  IF(iopmvis.EQ.1)THEN
                  WRITE(88,203)   'triangle { <',
     $                 coor(1,jk1(vol_item%lmvol(2),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(2),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(2),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(3),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(3),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(3),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(4),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(4),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(4),ip))+shift,
     $                 '> texture {T'//ai4//'} }'
                  IF(ioformat.EQ.0)THEN
                     WRITE(iunit) vol_item%lmvol(2),vol_item%lmvol(3),
     $                    vol_item%lmvol(4)
                  ELSE
                     WRITE(iunit,'(3i10)') vol_item%lmvol(2),vol_item%lmvol(3),
     $                    vol_item%lmvol(4)
                  ENDIF
                  ENDIF
                  iaux89(ip) = iaux89(ip) + 1
               ELSE IF(iflag89(k).EQ.2)THEN
                  IF(iopmvis.EQ.1)THEN
                  WRITE(88,203)   'triangle { <',
     $                 coor(1,jk1(vol_item%lmvol(1),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(1),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(1),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(3),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(3),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(3),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(4),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(4),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(4),ip))+shift,
     $                 '> texture {T'//ai4//'} }'
                  IF(ioformat.EQ.0)THEN
                     WRITE(iunit) vol_item%lmvol(1),vol_item%lmvol(3),
     $                    vol_item%lmvol(4)
                  ELSE
                     WRITE(iunit,'(3i10)') vol_item%lmvol(1),vol_item%lmvol(3),
     $                    vol_item%lmvol(4)
                  ENDIF
                  ENDIF
                  iaux89(ip) = iaux89(ip) + 1
               ELSEIF(iflag89(k).EQ.3)THEN
                  IF(iopmvis.EQ.1)THEN
                  WRITE(88,203)   'triangle { <',
     $                 coor(1,jk1(vol_item%lmvol(1),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(1),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(1),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(2),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(2),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(2),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(4),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(4),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(4),ip))+shift,
     $                 '> texture {T'//ai4//'} }'
                  IF(ioformat.EQ.0)THEN
                     WRITE(iunit) vol_item%lmvol(1),vol_item%lmvol(2),
     $                    vol_item%lmvol(4)
                  ELSE
                     WRITE(iunit,'(3i10)') vol_item%lmvol(1),vol_item%lmvol(2),
     $                    vol_item%lmvol(4)
                  ENDIF
                  ENDIF
                  iaux89(ip) = iaux89(ip) + 1
               ELSEIF(iflag89(k).EQ.4)THEN
                  IF(iopmvis.EQ.1)THEN
                  WRITE(88,203)   'triangle { <',
     $                 coor(1,jk1(vol_item%lmvol(1),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(1),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(1),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(2),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(2),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(2),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(3),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(3),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(3),ip))+shift,
     $                 '> texture {T'//ai4//'} }'
                  IF(ioformat.EQ.0)THEN
                     WRITE(iunit) vol_item%lmvol(1),vol_item%lmvol(2),
     $                    vol_item%lmvol(3)
                  ELSE
                     
                     WRITE(iunit,'(3i10)') vol_item%lmvol(1),vol_item%lmvol(2),
     $                    vol_item%lmvol(3)
                  ENDIF
                  ENDIF
                  iaux89(ip) = iaux89(ip) + 1
               ENDIF
            ENDDO
            IF(ipressflag(i).NE.0)THEN
               vol_item%iface = ipressflag(i)
               vol_item%press = pressload(1,i)
               nbcele(ip) = nbcele(ip) + 1
            ELSE
               vol_item%iface = 0
               vol_item%press = 0.d0
            ENDIF

C     Add item to volumetric element list
                        
            CALL vol_insert_tail(vol_list(imat,ip),vol_item)

            CALL vol_insert_ptr(vol_list(imat,ip),ik1_c4(i))
 
            numel(imat,ip) = numel(imat,ip) + 1
            numnp(ip) = ii + 3
         ELSE IF(iizones(i).EQ.0)THEN ! nocohesive elements
            numel(imat,ip) = numel(imat,ip) + 1
            ElOnPartBndry(i) = .FALSE.
            iflag89(1:4) = 0
            IF( SUM( bord_node_org( lmelv_prmry(:,i) ) ).GT.0) ElOnPartBndry(i) = .TRUE.
            DO k = 1, numvertx
               nk = lmelv_prmry(k,i)
               IF(k.LE.4)THEN
                  IF(neigh(k,i).LE.0)THEN
                     iflag89(k) = k
                  ENDIF
               ENDIF
               IF(ik1(nk,ip).EQ.0)THEN
                  numnp(ip) = numnp(ip) + 1
                  IF(bord_node_org(nk).EQ.1)THEN
                     node(nk)%proc_list(epart(i)) = numnp(ip)
                  ENDIF
                  DO jj=1,num_rco_border
                     IF(bord_node_Rco(jj).EQ.nk)THEN
                        node_rco(jj)%proc_list_rco(ip) = numnp(ip)
                     ENDIF
                  ENDDO
                  jk1(numnp(ip),ip) = nk
                  ik1(nk,ip) = numnp(ip)
C     
C -- Transfer the bc to the new node

                  IF(ibcflg(lmelv_prmry(k,i)).GT.0.AND.
     $                 ibcflg(lmelv_prmry(k,i)).LT.1000)THEN
                     numbc(ip)=numbc(ip)+1
                     ibc(numbc(ip),ip)%node = numnp(ip)
                     ibc(numbc(ip),ip)%bc = ibcflg(lmelv_prmry(k,i))
                  ENDIF
                  IF(ibcflg_mm(lmelv_prmry(k,i)).GT.0.AND.
     $                 ibcflg_mm(lmelv_prmry(k,i)).LT.1000)THEN
                     numbc_mm(ip)=numbc_mm(ip)+1 
                     ibc_mm(numbc_mm(ip),ip)%node = numnp(ip)
                     ibc_mm(numbc_mm(ip),ip)%bc = ibcflg_mm(lmelv_prmry(k,i))
                  ENDIF
               ENDIF
               vol_item%mat_vol = matvol(i)
               vol_item%lmvol(k) = ik1(nk,ip)
            ENDDO
            IF(ipressflag(i).NE.0)THEN
               vol_item%iface = ipressflag(i)
               vol_item%press = pressload(1,i)
               nbcele(ip) = nbcele(ip) + 1
            ELSE
               vol_item%iface = 0
               vol_item%press = 0.d0
            ENDIF
            WRITE(ai4,'(i3.3)') ip
            IF(iopmvis.EQ.1)THEN
            iunit = 2000 + ip - 1
            DO k = 1,4
               IF(iflag89(k).EQ.1)THEN
                  shift = .5
                  WRITE(88,203)   'triangle { <',
     $                 coor(1,jk1(vol_item%lmvol(2),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(2),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(2),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(3),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(3),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(3),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(4),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(4),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(4),ip))+shift,
     $                 '> texture {T'//ai4//'} }'
                  IF(ioformat.EQ.0)THEN
                     WRITE(iunit) vol_item%lmvol(2),vol_item%lmvol(3),
     $                    vol_item%lmvol(4)
                  ELSE
                     WRITE(iunit,'(3i10)') vol_item%lmvol(2),vol_item%lmvol(3),
     $                    vol_item%lmvol(4)  
                  ENDIF
                  iaux89(ip) = iaux89(ip) + 1
               ELSE IF(iflag89(k).EQ.2)THEN
                  WRITE(88,203)   'triangle { <',
     $                 coor(1,jk1(vol_item%lmvol(1),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(1),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(1),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(3),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(3),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(3),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(4),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(4),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(4),ip))+shift,
     $                 '> texture {T'//ai4//'} }'
                  IF(ioformat.EQ.0)THEN
                     WRITE(iunit) vol_item%lmvol(1),vol_item%lmvol(3),
     $                    vol_item%lmvol(4)
                  ELSE
                     WRITE(iunit,'(3i10)') vol_item%lmvol(1),vol_item%lmvol(3),
     $                    vol_item%lmvol(4)
                  ENDIF
                  iaux89(ip) = iaux89(ip) + 1
               ELSEIF(iflag89(k).EQ.3)THEN
                  WRITE(88,203)   'triangle { <',
     $                 coor(1,jk1(vol_item%lmvol(1),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(1),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(1),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(2),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(2),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(2),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(4),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(4),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(4),ip))+shift,
     $                 '> texture {T'//ai4//'} }'
                  IF(IOformat.EQ.0)THEN
                     WRITE(iunit) vol_item%lmvol(1),vol_item%lmvol(2),
     $                    vol_item%lmvol(4)
                  ELSE
                     WRITE(iunit,'(3i10)') vol_item%lmvol(1),vol_item%lmvol(2),
     $                    vol_item%lmvol(4)
                  ENDIF
                  iaux89(ip) = iaux89(ip) + 1
               ELSEIF(iflag89(k).EQ.4)THEN
                  WRITE(88,203)   'triangle { <',
     $                 coor(1,jk1(vol_item%lmvol(1),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(1),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(1),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(2),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(2),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(2),ip))+shift,'>, <',
     $                 coor(1,jk1(vol_item%lmvol(3),ip))+shift,
     $                 coor(2,jk1(vol_item%lmvol(3),ip))+shift,
     $                 coor(3,jk1(vol_item%lmvol(3),ip))+shift,
     $                 '> texture {T'//ai4//'} }'
                  IF(IOformat.EQ.0)THEN
                     WRITE(iunit) vol_item%lmvol(1),vol_item%lmvol(2),
     $                    vol_item%lmvol(3)
                  ELSE
                     WRITE(iunit,'(3i10)') vol_item%lmvol(1),vol_item%lmvol(2),
     $                    vol_item%lmvol(3)
                  ENDIF
                  iaux89(ip) = iaux89(ip) + 1
               ENDIF
            ENDDO
            ENDIF
!
!     Add item to volumetric element list
!
            
            IF(ElOnPartBndry(i))THEN ! Element contains a node that is on the partition boundary
               CALL vol_insert_head(vol_list(imat,ip),vol_item)
               icnt(imat,1,ip) = icnt(imat,1,ip) + 1
               MapGlbEl2LocEl(i) = icnt(imat,1,ip)
            ELSE
               CALL vol_insert_tail(vol_list(imat,ip),vol_item)
               icnt(imat,2,ip) = icnt(imat,2,ip) + 1
               MapGlbEl2LocEl(i) = icnt(imat,2,ip)
            ENDIF
         ENDIF
      ENDDO

      allocate(NumBorderVol(1:NumMat,1:nprocs)) ! made array so that I don't have to pass to 2d surface routine

      DO i = 1 , nprocs
         DO ii = 1, NumMat
            NumBorderVol(ii,i) = vol_list(ii,i)%num_border_vol
         ENDDO
      ENDDO

      IF(iopmvis.EQ.1)THEN
      DO i = 1, nprocs
         iunit = 2000 + i - 1
         CLOSE(iunit)
      ENDDO
      ENDIF

c$$$      OPEN(302,FILE='fort.301')
c$$$      REWIND(302)
c$$$      DO i = 1, numel_tri
c$$$         READ(302,'(6i10)') iaux,iaux,lmtri(1:3,i),iaux1
c$$$c         READ(302,'(6i10)') iaux,iaux,jk1(lmtri(1:3,i)),iaux1
c$$$         elm_2D(i) = iaux1
c$$$      ENDDO
c$$$      CLOSE(302)
c$$$      ALLOCATE(neigh_2d(1:3,1:numel_tri))
c$$$      neigh_2d = 0


      
! find the spike situation
      maxdim = (MAXVAL(numel_2d(1:nprocs,1))+MAXVAL(numel_2d(1:nprocs,2)))*6
      ALLOCATE(jk1_2D(1:maxdim,1:nprocs,1:2))
      jk1_2D = 0
      
c$$$ 
c$$$      DO i = 1, nprocs
c$$$         CALL FindSpike(vol_list(i),ip,NumNp2D,nprocs,jk1_2D,ik1_2D,ik1,maxdim)
c$$$      ENDDO

      IF(numel_tri.NE.0)THEN
         print*,'Rocface data'
         CALL rocface_data(lmtri,nprocs,numnp_prmry,numel_tri,epart_2d,
     $        coor,prefx,prefx_lngth,elm_2D,ik1,numelv_prmry,numvertx2d,
     $        numel_2d,elm_2D_flag,
     $        jk1_2D,maxdim,NumNp,jk1,jk1_size,iNdsBurnFlg,NumBorderVol,NumMat,MapGlbEl2LocEl,ElOnPartBndry,MatId )
         

         print*,'Rocface data finished'
      ENDIF
      deallocate(NumBorderVol)
      DEALLOCATE(iNdsBurnFlg)
      DEALLOCATE(jk1_2D)
      DEALLOCATE(ik1)
      DEALLOCATE(ibcaxi)
      DEALLOCATE(press_nodal)

C     WRITE THE COORDINATES
c
      ALLOCATE(iaux_p(1:nprocs))

      IF(IOFORMAT.EQ.0)THEN
         PRINT*,'Creating Binary output files'
      ELSE
         PRINT*,'Creating Ascii output files'
      ENDIF
      print*,'NumNdHistory',NumNdHistory
      IF(NumNdHistory.NE.0)THEN
         allocate(NumNdHistoryP(1:nprocs))
         NumNdHistoryP(1:nprocs) = 0
         allocate(NdHistoryP(1:10,1:nprocs))
         NdHistoryP(:,:) = 0
         NdHistoryFlag(:)=0
      ENDIF 
      DO i = 1, nprocs
         WRITE(ichr4,'(i4.4)') i - 1
         IF(IOFORMAT.EQ.0)THEN
            OPEN(4000,FILE=
     $           prefx(1:prefx_lngth)//'/'//prefx(1:prefx_lngth)//'.'//ichr4//'.inp',
     $           STATUS='replace',FORM='unformatted')         
!!            WRITE(4000) SUM(numnp(1:nprocs)),SUM(numel(1:nprocs))
            WRITE(4000) dt_courant ! '.005' ! dt_courant
            WRITE(4000) 4
            WRITE(4000) 1, 1.d0 ! 0.d0
            WRITE(4000) 30, 1.d0
            WRITE(4000) 200000, 1.d0 ! ntime, 1.d0
            WRITE(4000) 2000000, 1.d0 ! 10*ntime, 1.d0         
c
c -- Output the nodal cordinates
            iaux_p(:) = 0
            WRITE(4000) numnp(i)
            DO j=1,numnp(i)
               WRITE(4000) j,coor(1:3,jk1(j,i))
               DO k = 1, NumNodeIO
                  iaux_p(i) = iaux_p(i) + 1
                  IF(jk1(j,i).EQ.NodeIO(k)) 
     $                 PRINT*,'Old Node Number',NodeIO(k),' =', j,'proc ', i-1
               ENDDO 
            ENDDO
            
         ELSE

! Output To Each Processors Files

            OPEN(4000,FILE=
     $           prefx(1:prefx_lngth)//'/'//prefx(1:prefx_lngth)//'.'//ichr4//'.inp',
     $           STATUS='replace',FORM='formatted')         
!            WRITE(4000,*) SUM(numnp(1:nprocs)),SUM(numel(1:nprocs))
!            WRITE(4000,*) dt_courant ! '.005' ! dt_courant
!            WRITE(4000,*) 4
!            WRITE(4000,*) 1, 1.d0 ! 0.d0
!            WRITE(4000,*) 30, 1.d0
!            WRITE(4000,*) 200000, 1.d0 ! ntime, 1.d0
!            WRITE(4000,*) 2000000, 1.d0 ! 10*ntime, 1.d0

         
!
! -- Output the nodal cordinates
!
            ! Version of ROCSTAR_DATA
            WRITE(4000,*) 1
            WRITE(4000,*) 2.5 !,dt_courant, dhmin

            ! Number
            WRITE(4000,*) 2
            iaux_p(:) = 0
            WRITE(4000,*) numnp(i),0,0,0,0
            DO j=1,numnp(i)
               WRITE(4000,'(i9,3(1x,e16.9),2i9)') j,coor(1:3,jk1(j,i)),0
               DO k = 1, NumNdHistory
                  iaux_p(i) = iaux_p(i) + 1
                  IF(jk1(j,i).EQ.NdHistory(k).AND.NdHistoryFlag(k).EQ.0)THEN
                     NumNdHistoryP(i) = NumNdHistoryP(i) + 1
                     NdHistoryP(NumNdHistoryP(i),i) = j
                     NdHistoryFlag(k) = 1
                  ENDIF
               ENDDO  
  
            ENDDO
         ENDIF
         IF(NumNdHistory.NE.0)THEN
            IF(NumNdHistoryP(i).NE.0)THEN
               WRITE(4000,*) 7
               WRITE(4000,*) NumNdHistoryP(i)
               DO j = 1, NumNdHistoryP(i)
                  WRITE(4000,*) NdHistoryP(j,i)
               ENDDO
            ENDIF
         ENDIF
         
         
         CLOSE(4000)
      ENDDO

c      DEALLOCATE(coor)

      DEALLOCATE(matvol)

      IF(iopmvis.EQ.1) WRITE(88,'(a1)') '}'

      CLOSE(88)

         
 203  FORMAT(a12,2(3(1x,f12.8),1x,a4),3(1x,f12.8),1x,a19)

C
C --  Read element nodes that make up the structure boundary

c-avs      ii = 0
c-avs      ALLOCATE(ncoor(1:numnp_prmry))
c-avs      ncoor(:) = -1
! fix for ansys option
c$$$      ALLOCATE(lmtri_2d(1:3,1:numel_2d))
c$$$      IF(iansys.EQ.1)THEN
c$$$         DO i = 1, numel_2d
c$$$            READ(19,*) nn,nface,lmtri(1)
c$$$            READ(19,*) lmtri(2)
c$$$            READ(19,*) lmtri(3) 
c$$$cc            WRITE(111,'(5i10)') nn,nface,ik1(lmtri(1:3),1)
c$$$            lmtri_2d(1:3,i) = ik1(lmtri(1:3),1)
c$$$            ncoor(ik1(lmtri(1:3),1)) = 1
c$$$            ii = ii + 1
c$$$            IF(ii.EQ.7 .AND.i.NE.numel_2d)THEN
c$$$               READ(19,'()')
c$$$               READ(19,'()')
c$$$               READ(19,'()')
c$$$               ii = 0
c$$$            ENDIF
c$$$         ENDDO
c$$$      ENDIF
c-avs      WRITE(2001,'(5i8)')  numnp_2d,numel_2d, 1, 0, 0  
c-avs
c-avs      jj = 0
c-avs      DO i = 1, numnp(1)
c-avs         IF(ncoor(i).GT.0)THEN
c-avs            jj = jj + 1
ccc            WRITE(111,'(1i10,3(1x,f16.12))') i,coor(1:3,jk1(i,1))
c-avsc-avs            WRITE(2001,'(i8,3e14.6))') i,coor(1:3,jk1(i,1))
c-avs         ENDIF
c-avs      ENDDO
c      print*,'number of boundary nodes=',jj,ik1(6,1)
c      DO j=1,numel_2d
c         WRITE(2001,'(2i8,x,1a3,x,10i8)') j,1,'tri',ik1(lmtri_2d(1:3,j),1)
c      ENDDO
c$$$      OPEN(302,FILE='fort.301')
c$$$      REWIND(302)
c$$$      DO j = 1, numel_2d
c$$$         READ(302,*) jj,jj,lmtri_2d(1:3,j)
c$$$         WRITE(2001,'(2i8,x,1a3,x,10i8)') j,1,'tri',lmtri_2d(1:3,j)
c$$$      ENDDO
c$$$      CLOSE(302)
c-avs 1400 FORMAT(e14.5,e14.5,e14.5,e14.5)
c-avs      WRITE(2001,'(6i8)') 1,1
c-avs      WRITE(2001,'(a18)') 'Displacement,units' 
c-avs      DO i = 1, numnp(1)
c-avs         IF(ncoor(i).GT.0)THEN
c-avs            WRITE(2001,'(i8,3e14.6))') i,1.d0
c-avs         ENDIF
c-avs      ENDDO
c
c  2) CREATE COHESIVE ELEMENTS
c     Create list such that elements that are on a parallel partition boundary 
C     are listed at the begining of the list, the interior elements then follow
C
      PRINT*,' .....creating cohesive element connectivity table....'

c***      ALLOCATE(lmcoh(1:6,1:numelv_prmry*4,1:nprocs))
      ALLOCATE(numclst(1:nprocs))
      ALLOCATE(num_border_coh(1:nprocs))
      num_border_coh = 0
      numclst(1:nprocs) = 0

      ALLOCATE(coh_item)
      ALLOCATE(coh_list(1:nprocs))
      ALLOCATE(coh_comm_list(1:nprocs,1:nprocs))
      ALLOCATE(coh_comm_item)
      coh_comm_list(:,:)%num_comm_nodes = 0
      PRINT*,' .....creating cohesive element connectivity table....'
      DO i = 1, numelv_prmry
         ip = epart(i)
         IF(iizones(i).EQ.1)THEN
            DO j = 1, 4
               IF(neigh(j,i).GT.0)THEN
                  flag = 0
                  IF(neigh(j,i).GT.i.OR.ip.NE.epart(neigh(j,i)).OR.
     $                 (iizones(neigh(j,i)).EQ.0 ))THEN
                     IF(ip.NE.epart(neigh(j,i))) flag=epart(neigh(j,i))
                     numclst(ip) = numclst(ip) + 1
                     IF(j.EQ.1)THEN
                        IF(matzone(i).NE.matzone(neigh(j,i)))THEN ! NOTWORK
                           coh_item%mat_coh = ibcflg(lmelv_prmry(5,i))/1000
                        ELSE
                           coh_item%mat_coh = matzone(i)
                        ENDIF
                        coh_item%lmcoh(1) = ik1_c4(i)%ptr%vol_info%lmvol(4)
                        coh_item%lmcoh(2) = ik1_c4(i)%ptr%vol_info%lmvol(3)
                        coh_item%lmcoh(3) = ik1_c4(i)%ptr%vol_info%lmvol(2)
                        ni = jk1(ik1_c4(i)%ptr%vol_info%lmvol(2),ip)
                        node1 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(1)
                        node2 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(2)
                        node3 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(3)
                        node4 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(4)
                        gnode1 = jk1(node1,epart(neigh(j,i)))
                        gnode2 = jk1(node2,epart(neigh(j,i)))
                        gnode3 = jk1(node3,epart(neigh(j,i)))
                        gnode4 = jk1(node4,epart(neigh(j,i)))

                        CALL coh_opposite_side(numelv_prmry,ni,node1,node2,node3,node4,
     $                       gnode1,gnode2,gnode3,gnode4,neigh,j,i,coh_item)
                     ELSE IF(j.EQ.2)THEN
                        IF(matzone(i).NE.matzone(neigh(j,i)))THEN
                           coh_item%mat_coh = ibcflg(lmelv_prmry(6,i))/1000
                        ELSE
                           coh_item%mat_coh = matzone(i)
                        ENDIF
                        coh_item%lmcoh(1) =  ik1_c4(i)%ptr%vol_info%lmvol(4)
                        coh_item%lmcoh(2) =  ik1_c4(i)%ptr%vol_info%lmvol(1)
                        coh_item%lmcoh(3) =  ik1_c4(i)%ptr%vol_info%lmvol(3)
                        ni = jk1(ik1_c4(i)%ptr%vol_info%lmvol(3),ip)
                        node1 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(1)
                        node2 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(2)
                        node3 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(3)
                        node4 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(4)
                        gnode1 = jk1(node1,epart(neigh(j,i)))
                        gnode2 = jk1(node2,epart(neigh(j,i)))
                        gnode3 = jk1(node3,epart(neigh(j,i)))
                        gnode4 = jk1(node4,epart(neigh(j,i)))

                        CALL coh_opposite_side(numelv_prmry,ni,node1,node2,node3,node4,
     $                       gnode1,gnode2,gnode3,gnode4,neigh,j,i,coh_item)

                     ELSE IF(j.EQ.3)THEN
                        IF(matzone(i).NE.matzone(neigh(j,i)))THEN
                           coh_item%mat_coh = ibcflg(lmelv_prmry(4,i))/1000
                        ELSE
                           coh_item%mat_coh = matzone(i)
                        ENDIF
                        coh_item%lmcoh(1) = ik1_c4(i)%ptr%vol_info%lmvol(4)
                        coh_item%lmcoh(2) = ik1_c4(i)%ptr%vol_info%lmvol(2)
                        coh_item%lmcoh(3) = ik1_c4(i)%ptr%vol_info%lmvol(1)
                        ni = jk1(ik1_c4(i)%ptr%vol_info%lmvol(1),ip)
                        node1 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(1)
                        node2 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(2)
                        node3 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(3)
                        node4 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(4)
                        gnode1 = jk1(node1,epart(neigh(j,i)))
                        gnode2 = jk1(node2,epart(neigh(j,i)))
                        gnode3 = jk1(node3,epart(neigh(j,i)))
                        gnode4 = jk1(node4,epart(neigh(j,i)))

                        CALL coh_opposite_side(numelv_prmry,ni,node1,node2,node3,node4,
     $                       gnode1,gnode2,gnode3,gnode4,neigh,j,i,coh_item)

                     ELSE IF(j.EQ.4)THEN
                        IF(matzone(i).NE.matzone(neigh(j,i)))THEN
                           coh_item%mat_coh = ibcflg(lmelv_prmry(4,i))/1000
                        ELSE
                           coh_item%mat_coh = matzone(i)
                        ENDIF
                        coh_item%lmcoh(1) =  ik1_c4(i)%ptr%vol_info%lmvol(1)
                        coh_item%lmcoh(2) =  ik1_c4(i)%ptr%vol_info%lmvol(2)
                        coh_item%lmcoh(3) =  ik1_c4(i)%ptr%vol_info%lmvol(3)
                        ni = jk1(ik1_c4(i)%ptr%vol_info%lmvol(3),ip)
                        node1 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(1)
                        node2 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(2)
                        node3 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(3)
                        node4 = ik1_c4(neigh(j,i))%ptr%vol_info%lmvol(4)
                        gnode1 = jk1(node1,epart(neigh(j,i)))
                        gnode2 = jk1(node2,epart(neigh(j,i)))
                        gnode3 = jk1(node3,epart(neigh(j,i)))
                        gnode4 = jk1(node4,epart(neigh(j,i)))

                        CALL coh_opposite_side(numelv_prmry,ni,node1,node2,node3,node4,
     $                       gnode1,gnode2,gnode3,gnode4,neigh,j,i,coh_item)
                     ENDIF
                     coh_item%clst_type = flag - 1
C
C add cohesive element to the list
C
                     coh_comm_item%lmcoh_comm(1) = coh_item%lmcoh(4)
                     coh_comm_item%lmcoh_comm(2) = coh_item%lmcoh(5)
                     coh_comm_item%lmcoh_comm(3) = coh_item%lmcoh(6)
                     IF(coh_item%clst_type.GE.0)THEN
                        CALL coh_insert_head(coh_list(ip),coh_item)
                        CALL coh_comm(coh_comm_list(epart(neigh(j,i)),ip),coh_comm_item)
                     ELSE
                        CALL coh_insert_tail(coh_list(ip),coh_item)
                     ENDIF

                     IF(flag-1.GE.0)
     $                    num_border_coh(ip)=num_border_coh(ip)+1
c
c     for the special case of a zone boundary
c
                     IF(iizones(neigh(j,i)).NE.iizones(i).AND.
     $                    epart(neigh(j,i)).NE.ip)THEN
                        numclst(epart(neigh(j,i)))=
     $                       numclst(epart(neigh(j,i))) + 1
                        IF(j.EQ.1)THEN
                           coh_item%mat_coh = ibcflg(lmelv_prmry(5,i))/1000
c                           matcoh(numclst(epart(neigh(j,i))),
c     $                          epart(neigh(j,i))) = 
c     $                          ibcflg(lmelv_prmry(5,i))/100
                        ELSE IF(j.EQ.2)THEN
                           coh_item%mat_coh = ibcflg(lmelv_prmry(6,i))/1000
c                           matcoh(numclst(epart(neigh(j,i))),
c     $                          epart(neigh(j,i))) = 
c     $                          ibcflg(lmelv_prmry(6,i))/100  
                        ELSE
                           coh_item%mat_coh = ibcflg(lmelv_prmry(4,i))/1000
c                           matcoh(numclst(epart(neigh(j,i))),
c     $                          epart(neigh(j,i))) =
c     $                          ibcflg(lmelv_prmry(4,i))/100
                        ENDIF
ccc***                        coh_item%lmcoh(1) = lmcoh(3,numclst(ip),ip)
ccc***                        coh_item%lmcoh(2) = lmcoh(4,numclst(ip),ip)
ccc***                        coh_item%lmcoh(3) = lmcoh(1,numclst(ip),ip)
ccc***                        coh_item%lmcoh(4) = lmcoh(2,numclst(ip),ip)
ccc***                        coh_item%lmcoh(5) = lmcoh(6,numclst(ip),ip)
ccc***                        coh_item%lmcoh(6) = lmcoh(5,numclst(ip),ip)

                        coh_item%clst_type = ip-1

                        coh_comm_item%lmcoh_comm(1) = coh_item%lmcoh(4)
                        coh_comm_item%lmcoh_comm(2) = coh_item%lmcoh(5)
                        coh_comm_item%lmcoh_comm(3) = coh_item%lmcoh(6)

                        IF(coh_item%clst_type.GE.0)THEN
                           CALL coh_insert_head(coh_list(epart(neigh(j,i))),coh_item)
                           CALL coh_comm(coh_comm_list(epart(neigh(j,i)),ip),coh_comm_item)
                        ELSE
                           CALL coh_insert_tail(coh_list(epart(neigh(j,i))),coh_item)
                        ENDIF
                        IF(coh_item%clst_type.GE.0)
     $                       num_border_coh(epart(neigh(j,i))) =
     $                       num_border_coh(epart(neigh(j,i))) + 1
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      PRINT*,' Finished  cohesive element connectivity table....'
      DEALLOCATE(ik1_c4)
      DEALLOCATE(jk1)
      DEALLOCATE(matzone)
      DEALLOCATE(iizones)
ccc***      DEALLOCATE(lmcoh)
c      DEALLOCATE(node) !!! rm if have no cohesive elements

      PRINT*,'Setting up communication'
C
C Cohesive elements communciation

      ALLOcate(ID_sendto(1:nprocs,1:nprocs))
      ID_sendto(1:nprocs,1:nprocs)%num_border_comm = 0
      DO i = 1,numnp_prmry
         IF(bord_node_org(i).EQ.1)THEN         
            DO j = 1,nprocs
               DO k = 1,nprocs
                  IF(k.NE.j.AND.node(i)%proc_list(k).NE.-1)THEN
                     IF(node(i)%proc_list(j).GE.0)THEN
                        CALL addcommNd(ID_sendto(j,k),node(i)%proc_list(j))
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ENDDO

c$$$      DO i = 1,numnp_prmry
c$$$         IF(bord_node_org(i).EQ.1)THEN         
c$$$            DO j = 1,nprocs
c$$$               DO k = 1,nprocs
c$$$                  IF(k.NE.j.AND.node(i)%proc_list(k).NE.-1)THEN
c$$$                     p2max_lst = p2max_lst + 1
c$$$                     sndrcvnod_lst(p2max_lst)%sndn=j
c$$$                     sndrcvnod_lst(p2max_lst)%rcvn=k
c$$$                     sndrcvnod_lst(p2max_lst)%nodes=node(i)%proc_list(j)
c$$$                     IF(sndrcvnod_lst(p2max_lst)%nodes.GE.0)
c$$$     $                    nnntemp_lst(j,k) = nnntemp_lst(j,k) + 1
c$$$                     DO jj=1,num_rco_border
c$$$                        IF(sndrcvnod_lst(p2max_lst)%nodes.EQ.
c$$$     $                       node_rco(jj)%proc_list_rco(j))THEN
c$$$                           nnntemp_rco(j,k) = nnntemp_rco(j,k) + 1
c$$$                        ENDIF
c$$$                     ENDDO
c$$$                  ENDIF
c$$$               ENDDO
c$$$            ENDDO
c$$$         ENDIF
c$$$      ENDDO

c ----------------------------------------------------------------------
c -- OUTPUT  ----------------------------------------------------
c ---------------------------------------------------------

      WRITE(*,*)' ------------------------------------------------'
      WRITE(*,*)' ------------- MESH DESCRIPTION -----------------'
      WRITE(*,*)' ------------------------------------------------'
      WRITE(*,*)' Number of volumetric lst elements= ',SUM(numel)
      WRITE(*,*)' Number of cohesive 6-node elements= ',SUM(numclst)
      WRITE(*,*)' Number on Nodes= ',SUM(numnp)
      WRITE(*,*)' Average nodes per processor= ',
     $     DBLE(SUM(numnp))/DBLE(nprocs)
      WRITE(*,*)' Average elements per processor= ',
     $     DBLE(SUM(numel))/DBLE(nprocs)

c
c    Write each processor's input file
c
      ALLOCATE(nproc_neigh_lst(1:nprocs))
      ALLOCATE(nproc_neigh_rco(1:nprocs))
      ALLOCATE(nproc_neigh(1:nprocs))
      nproc_neigh = 0
      nproc_neigh_lst = 0
      nproc_neigh_rco = 0
C
C-- Make an input file for each processor to read
C     Loop over all the processors

!!      OPEN(113,FILE=prefx(1:prefx_lngth)//'/max_arrybnds.inp')
      DO i=1,nprocs
C
C     Determine the neighbor of processors 'i' is communicating with.
C     
         DO j = 1,nprocs
            IF(coh_comm_list(i,j)%num_comm_nodes.NE.0) 
     $           nproc_neigh(i) = nproc_neigh(i) + 1
            IF(ID_sendto(i,j)%num_border_comm.NE.0) 
     $           nproc_neigh_lst(i) = nproc_neigh_lst(i) + 1
!            IF(nnntemp_rco(i,j).NE.0) 
!     $           nproc_neigh_rco(i) = nproc_neigh_rco(i) + 1
         ENDDO
         WRITE(ichr4,'(i4.4)') i - 1
         IF(IOformat.EQ.0)THEN

            OPEN(4000,FILE=
     $           prefx(1:prefx_lngth)//'/'//prefx(1:prefx_lngth)//'.'//ichr4//'.inp',
     $           STATUS='OLD',FORM='unformatted', POSITION = 'append')
C
C --  Write Nodal Boundary Conditions
C
c
c     3) Boundary conditions
c
c      DO i=1,4
c         READ(3,20)
c         READ(3,*)iboundx(i),iboundy(i),boundx(i),boundy(i)
c      ENDDO
c
c      WRITE(*,*)'     top surface : '
c      IF(iboundx(1).EQ.1)THEN
c         WRITE(*,*)'           x-traction set at ',boundx(1)
c      ELSE
c         WRITE(*,*)'       x-displacement set at ',boundx(1)
c      ENDIF
c      IF(iboundy(1).EQ.1)THEN
c         WRITE(*,*)'           y-traction set at ',boundy(1)
c      ELSE
c         WRITE(*,*)'       y-displacement set at ',boundy(1)
c      ENDIF
c         print*,numbc(i)
            print*,' Write Nodal Boundary Conditions'	
            WRITE(4000) numbc(i)
            DO j = 1,numbc(i)
               k = ibc(j,i)%bc
               WRITE(4000) ibc(j,i)%node,k
            ENDDO
            print*,' Write MeshMotionNodal Boundary Conditions'
            WRITE(4000) numbc_mm(i)
            DO j = 1,numbc_mm(i)
               k = ibc_mm(j,i)%bc
               WRITE(4000) ibc_mm(j,i)%node,k 
            ENDDO
c$$$         WRITE(4000) numbc(i)
c$$$         DO j = 1,numbc(i)
c$$$            k = ibc(j,i)%bc
c$$$            WRITE(4000) ibc(j,i)%node,bc_mshmtn(k)
c$$$         ENDDO
C
C --  WRITE COHESIVE ELEMENT CONNECTIVITY ARRAY
C
C No. of 6-node cohesive  elements
C No. of 12-node cohesive elements
C No. of cohesvie elements being communicated
C No. of neighboring processors
C No. of cohesive elements on the partition border

            WRITE(4000) numclst(i),0,SUM(coh_comm_list(i,1:nprocs)%num_comm_nodes),
     $           nproc_neigh(i),num_border_coh(i)

            CALL print_coh_list(coh_list(i))

            WRITE(4000) 0       ! fix fix
C
C --  WRITE VOLUMETRIC ELEMENT CONNECTIVITY ARRAY 
C
C No. of 4-node tetrahedral
C No. of 10-node tetrahedral
C No. of lst on the partioned mesh boundary
            DO ii = 1, NumMat

               WRITE(4000) numel(ii,i),0,vol_list(ii,i)%num_border_vol,numvertx ! ,numel_2d(i) ! fix nbcele(i)

               CALL print_vol_list(vol_list(ii,i),ii)
            ENDDO
            PRINT*,'Communication information'
C
C --  WRITE MPI COMMUNICATION INFORMATION
C
C     1) Cohesive Force calculation communciation

            ii = 1
            ALLOCATE(my_neigh(1:nproc_neigh(i)))
            DO j=1,nprocs
               IF(coh_comm_list(i,j)%num_comm_nodes.GT.0)THEN
                  my_neigh(ii) = j - 1
                  ii = ii + 1
               ENDIF
               CALL print_coh_comm(coh_comm_list(i,j))  
            ENDDO

C     Number of neighboring processors
       WRITE(4000) nproc_neigh(i)
C     List of neighboring processors
       IF(nproc_neigh(i).GT.0)
     $      WRITE(4000) my_neigh(1:nproc_neigh(i))

       DEALLOCATE(my_neigh)

c     2) Internal Force calculaton communciation

       ii=1
       ALLOCATE(my_neigh(1:nproc_neigh_lst(i)))
       DO j=1,nprocs            ! receiving processor
          IF(nnntemp_lst(i,j).GT.0)THEN
             my_neigh(ii)=j-1
             ii=ii+1
          ENDIF
C     Number of nodes that need to be communicated for R_in calculation
            WRITE(4000) nnntemp_lst(i,j) ! common
C     List of nodes that need to be communicated for R_in calculation
            DO k=1,p2max_lst
               IF(sndrcvnod_lst(k)%sndn.EQ.i.AND.
     $              sndrcvnod_lst(k)%rcvn.EQ.j
     $              .AND.sndrcvnod_lst(k)%nodes.GE.0) 
     $              WRITE(4000) sndrcvnod_lst(k)%nodes
            ENDDO
         ENDDO
C     Number of neighboring proc. involved in R_in calculation
         WRITE(4000) nproc_neigh_lst(i)
C     List these neighboring processors
         IF(nproc_neigh_lst(i).GT.0)
     $        WRITE(4000) my_neigh(1:nproc_neigh_lst(i))

         DEALLOCATE(my_neigh)
C
C -Cohesive/Non-Cohesive Boundary communication
         ii = 1
         ALLOCATE(my_neigh(1:nproc_neigh_rco(i)))
         DO j=1,nprocs          ! receiving processor
            IF(nnntemp_rco(i,j).GT.0)THEN
               my_neigh(ii) = j - 1
               ii = ii + 1
            ENDIF
C     Number of nodes needed on a cohesive/no-cohesive boundary
            WRITE(4000) nnntemp_rco(i,j) ! common
C     List of nodes that are on a cohesive/no-cohesive boundary 
C     (needed for R_co calulations)
            DO k=1,p2max_lst
               DO jj=1,num_rco_border
                  IF(sndrcvnod_lst(k)%sndn.EQ.i.AND.
     $                 sndrcvnod_lst(k)%rcvn.EQ.j.AND.
     $                 sndrcvnod_lst(k)%nodes.EQ.
     $                 node_rco(jj)%proc_list_rco(i))THEN
                     WRITE(4000) node_rco(jj)%proc_list_rco(i)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
C     Number of neighboring processors
         WRITE(4000) nproc_neigh_rco(i)
C     List of neighboring processors
         IF(nproc_neigh_rco(i).GT.0)
     $        WRITE(4000) my_neigh(1:nproc_neigh_rco(i))

         DEALLOCATE(my_neigh)

         WRITE(4000) iaux89(i)
         CLOSE(4000)
C
C --- WRITE OUT PRESSURE SURFACE

c         IF(iansys.EQ.1)THEN
c            OPEN(19,FILE=prefx(1:prefx_lngth)//'.press', FORM='unformatted')
c         ENDIF
c         IF(iansys.EQ.1)THEN
c            REWIND(19)
! fix for ansys and pressure
c$$$            DO j=1,5
c$$$               READ(19,'()')
c$$$            ENDDO
c$$$            ntri = 0
c$$$            ii = 0
c$$$            DO j = 1, numel_2d
c$$$               READ(19,*) nn,nface,lmtri(1)
c$$$               READ(19,*) lmtri(2)
c$$$               READ(19,*) lmtri(3)
c$$$               IF(ik1_el(nn,i).NE.0)THEN
c$$$                  WRITE(112,'(5i10)') ik1_el(nn,i),nface,ik1(lmtri(1:3),i)
c$$$                  ntri = ntri + 1
c$$$               ENDIF
c$$$               ii = ii + 1
c$$$               IF(ii.EQ.7 .AND.j.NE.numel_2d)THEN
c$$$                  READ(19,'()')
c$$$                  READ(19,'()')
c$$$                  READ(19,'()')
c$$$                  ii = 0
c$$$               ENDIF
c$$$            ENDDO
c$$$            WRITE(113,'(2i10)') i-1,ntri
c            CLOSE(19)
c         ELSE IF(ipatran.EQ.1.OR.itetmesh.EQ.1)THEN
            IF(ipress.EQ.1)THEN !fix this for pressure
ccnotneeded               OPEN(112,FILE=
ccnotneeded    $              prefx(1:prefx_lngth)//'/'//prefx(1:prefx_lngth)//'.'//ichr4//'.press', 
ccnotneeded     $              FORM='unformatted',STATUS='unknown')
ccnotneeded               REWIND(112)
               OPEN(302,FILE='fort.301')
               DO j = 1, numelv_prmry
                  READ(302,'(5i10)',IOSTAT=ios) ii,nface,i1,i2,i3
                                ! Check for errors and end of file
                  IF(ios.LT.0) EXIT ! end of file
                                ! nodes that make up surface triangle
ccnotneeded                  IF(ii.EQ.i) WRITE(112) nface,ik1(i1,i),ik1(i2,i),ik1(i3,i)
                  READ(302,'()') ! regression flag 0 -no regression, 1- regression 
               ENDDO
               CLOSE(302)
               CLOSE(301)
            ENDIF
!!            WRITE(113,'(4i10)') i-1,numel_2d(i),numnp_2d(i),iaux89(i) ! not correct for more then 1 pr.
               

c$$$        OPEN(111,FILE=
c$$$     $           prefx(1:prefx_lngth)//'/'//prefx(1:prefx_lngth)//'.'//ichr4//'.intrfc', 
c$$$     $           FORM='formatted',STATUS='unknown')
c$$$
c$$$            WRITE(111,*) numnp_2d(i)
c$$$
c$$$            OPEN(505,FILE='fort.505')
c$$$            kk = 1
c$$$            DO
c$$$               READ(505,*) jj, ii
c$$$               IF(jj.EQ.i)THEN
c$$$                  WRITE(111,'(2i10)') ik1(ii,i)
c$$$                  kk = kk + 1
c$$$               ENDIF
c$$$               IF(kk.GT.numnp_2d(i)) EXIT
c$$$            ENDDO
c$$$            CLOSE(111)
c$$$            CLOSE(505)
               
c         ENDIF


         ELSE

            OPEN(4000,FILE=
     $           prefx(1:prefx_lngth)//'/'//prefx(1:prefx_lngth)//'.'//ichr4//'.inp',
     $           STATUS='OLD',FORM='formatted', POSITION = 'append')
C
C --  Write Nodal Boundary Conditions
C
c
c     3) Boundary conditions
c
c      DO i=1,4
c         READ(3,20)
c         READ(3,*)iboundx(i),iboundy(i),boundx(i),boundy(i)
c      ENDDO
c
c      WRITE(*,*)'     top surface : '
c      IF(iboundx(1).EQ.1)THEN
c         WRITE(*,*)'           x-traction set at ',boundx(1)
c      ELSE
c         WRITE(*,*)'       x-displacement set at ',boundx(1)
c      ENDIF
c      IF(iboundy(1).EQ.1)THEN
c         WRITE(*,*)'           y-traction set at ',boundy(1)
c      ELSE
c         WRITE(*,*)'       y-displacement set at ',boundy(1)
c      ENDIF
c         print*,numbc(i)

         WRITE(4000,*) 3

         WRITE(4000,*) numbc(i),0
         DO j = 1,numbc(i)
            k = ibc(j,i)%bc
            WRITE(4000,'(4i9)') ibc(j,i)%node, k, 0
         ENDDO
          WRITE(4000,*) 4
         WRITE(4000,*) numbc_mm(i),0
         DO j = 1,numbc_mm(i)
            k = ibc_mm(j,i)%bc
            WRITE(4000,'(4i9)') ibc_mm(j,i)%node, k, 0
            IF(k.EQ.0)THEN
               PRINT*,'Zero boundary Condition FLAG'
               PRINT*,'processor =',i,j
               PRINT*,'Local node number',ibc_mm(j,i)%node
               STOP
            ENDIF
         ENDDO
         

c$$$         WRITE(4000,*) numbc(i)
c$$$         DO j = 1,numbc(i)
c$$$            k = ibc(j,i)%bc
c$$$            WRITE(4000,'(i9,3i9,3e14.5)') ibc(j,i)%node,bc_mshmtn(k)
c$$$         ENDDO
C
C --  WRITE COHESIVE ELEMENT CONNECTIVITY ARRAY
C
C No. of 6-node cohesive  elements
C No. of 12-node cohesive elements
C No. of cohesvie elements being communicated
C No. of neighboring processors
C No. of cohesive elements on the partition border

c$$$         WRITE(4000,'(5i9)') numclst(i),0,SUM(coh_comm_list(i,1:nprocs)%num_comm_nodes),
c$$$     $        nproc_neigh(i),num_border_coh(i)
c$$$
c$$$         CALL print_coh_list(coh_list(i))
C
C --  WRITE VOLUMETRIC ELEMENT CONNECTIVITY ARRAY 
C
C No. of 4-node tetrahedral
C No. of 10-node tetrahedral
C No. of lst on the partioned mesh boundary

         itmp1 = SUM(numel(1:NumMat,i))
         itmp2 = SUM(vol_list(1:NumMat,i)%num_border_vol)

         WRITE(4000,*) 5

         DO ii = 1, NumMat

            WRITE(4000,'(6i9)') itmp1,itmp2,numel(ii,i),vol_list(ii,i)%num_border_vol,numvertx,0

            CALL print_vol_list(vol_list(ii,i),ii)
         ENDDO
C
C --  WRITE MPI COMMUNICATION INFORMATION
C
C     1) Cohesive Force calculation communciation

c$$$       ii = 1
c$$$       ALLOCATE(my_neigh(1:nproc_neigh(i)))
c$$$       DO j=1,nprocs
c$$$          IF(coh_comm_list(i,j)%num_comm_nodes.GT.0)THEN
c$$$             my_neigh(ii) = j - 1
c$$$             ii = ii + 1
c$$$          ENDIF
c$$$          CALL print_coh_comm(coh_comm_list(i,j))  
c$$$       ENDDO
c$$$
c$$$C     Number of neighboring processors
c$$$       WRITE(4000,*) nproc_neigh(i)
c$$$C     List of neighboring processors
c$$$       IF(nproc_neigh(i).GT.0)
c$$$     $      WRITE(4000,'(100000i9)') my_neigh(1:nproc_neigh(i))
c$$$
c$$$       DEALLOCATE(my_neigh)

c     2) Internal Force calculaton communciation

         write(4000,*) 6

C     Number of neighboring proc. involved in R_in calculation
         WRITE(4000,*) nproc_neigh_lst(i)

!     List these neighboring processors

       ALLOCATE(my_neigh(1:nproc_neigh_lst(i)))
       DO j=1,nprocs            ! receiving processor
          IF(ID_sendto(i,j)%num_border_comm.NE.0)THEN
C     Number of nodes that need to be communicated for R_in calculation
             WRITE(4000,*) j-1,ID_sendto(i,j)%num_border_comm ! common
C     List of nodes that need to be communicated for R_in calculation
             CALL print_comm_list(ID_sendto(i,j))
          ENDIF
       ENDDO

!       IF(nproc_neigh_lst(i).GT.0)
!     $      WRITE(4000,'(10000i9)') my_neigh(1:nproc_neigh_lst(i))

       DEALLOCATE(my_neigh)
c$$$C
c$$$C -Cohesive/Non-Cohesive Boundary communication
c$$$         ii = 1
c$$$         ALLOCATE(my_neigh(1:nproc_neigh_rco(i)))
c$$$         DO j=1,nprocs          ! receiving processor
c$$$            IF(nnntemp_rco(i,j).GT.0)THEN
c$$$               my_neigh(ii) = j - 1
c$$$               ii = ii + 1
c$$$            ENDIF
c$$$C     Number of nodes needed on a cohesive/no-cohesive boundary
c$$$            WRITE(4000,*) nnntemp_rco(i,j) ! common
c$$$C     List of nodes that are on a cohesive/no-cohesive boundary 
c$$$C     (needed for R_co calulations)
c$$$            DO k=1,p2max_lst
c$$$               DO jj=1,num_rco_border
c$$$                  IF(sndrcvnod_lst(k)%sndn.EQ.i.AND.
c$$$     $                 sndrcvnod_lst(k)%rcvn.EQ.j.AND.
c$$$     $                 sndrcvnod_lst(k)%nodes.EQ.
c$$$     $                 node_rco(jj)%proc_list_rco(i))THEN
c$$$                     WRITE(4000,*) node_rco(jj)%proc_list_rco(i)
c$$$                  ENDIF
c$$$               ENDDO
c$$$            ENDDO
c$$$         ENDDO
c$$$C     Number of neighboring processors
c$$$         WRITE(4000,*) nproc_neigh_rco(i)
c$$$C     List of neighboring processors
c$$$         IF(nproc_neigh_rco(i).GT.0)
c$$$     $        WRITE(4000,*) my_neigh(1:nproc_neigh_rco(i))
c$$$
c$$$         DEALLOCATE(my_neigh)

!         WRITE(4000,*) iaux89(i)

       write(4000,*) 99
         CLOSE(4000)
C
C --- WRITE OUT PRESSURE SURFACE

c         IF(iansys.EQ.1)THEN
c            OPEN(19,FILE=prefx(1:prefx_lngth)//'.press', FORM='formatted')
c         ENDIF
c         IF(iansys.EQ.1)THEN
c            REWIND(19)
! fix for ansys and pressure
c$$$            DO j=1,5
c$$$               READ(19,'()')
c$$$            ENDDO
c$$$            ntri = 0
c$$$            ii = 0
c$$$            DO j = 1, numel_2d
c$$$               READ(19,*) nn,nface,lmtri(1)
c$$$               READ(19,*) lmtri(2)
c$$$               READ(19,*) lmtri(3)
c$$$               IF(ik1_el(nn,i).NE.0)THEN
c$$$                  WRITE(112,'(5i10)') ik1_el(nn,i),nface,ik1(lmtri(1:3),i)
c$$$                  ntri = ntri + 1
c$$$               ENDIF
c$$$               ii = ii + 1
c$$$               IF(ii.EQ.7 .AND.j.NE.numel_2d)THEN
c$$$                  READ(19,'()')
c$$$                  READ(19,'()')
c$$$                  READ(19,'()')
c$$$                  ii = 0
c$$$               ENDIF
c$$$            ENDDO
c$$$            WRITE(113,'(2i10)') i-1,ntri
c            CLOSE(19)
c         ELSE IF(ipatran.EQ.1.OR.itetmesh.EQ.1)THEN
!            IF(ipress.EQ.1)THEN
!notneeded               OPEN(112,FILE=
!notneeded     $              prefx(1:prefx_lngth)//'/'//prefx(1:prefx_lngth)//'.'//ichr4//'.press', 
!notneeded     $              FORM='formatted',STATUS='unknown')
!notneeded               REWIND(112)
!               OPEN(302,FILE='fort.301')
!               DO j = 1, numelv_prmry
!                  READ(302,'(5i10)',IOSTAT=ios) ii,nface,i1,i2,i3
                                ! Check for errors and end of file
!                  IF(ios.LT.0) EXIT ! end of file
                                ! nodes that make up surface triangle
!notneeded                  IF(ii.EQ.i) WRITE(112,'(5i10)') nface,ik1(i1,i),ik1(i2,i),ik1(i3,i)
!                  READ(302,'()') ! regression flag 0 -no regression, 1- regression 
!               ENDDO
!               CLOSE(302)
!!               WRITE(113,'(4i10)') i-1,numel_2d(i),numnp_2d(i),iaux89(i) ! not correct for more then 1 pr.
!               CLOSE(301)
!            ENDIF
               

c$$$        OPEN(111,FILE=
c$$$     $           prefx(1:prefx_lngth)//'/'//prefx(1:prefx_lngth)//'.'//ichr4//'.intrfc', 
c$$$     $           FORM='formatted',STATUS='unknown')
c$$$
c$$$            WRITE(111,*) numnp_2d(i)
c$$$
c$$$            OPEN(505,FILE='fort.505')
c$$$            kk = 1
c$$$            DO
c$$$               READ(505,*) jj, ii
c$$$               IF(jj.EQ.i)THEN
c$$$                  WRITE(111,'(2i10)') ik1(ii,i)
c$$$                  kk = kk + 1
c$$$               ENDIF
c$$$               IF(kk.GT.numnp_2d(i)) EXIT
c$$$            ENDDO
c$$$            CLOSE(111)
c$$$            CLOSE(505)
               
c         ENDIF
      ENDIF

      ENDDO
      CALL system('rm -f fort.302')

      STOP



      END
