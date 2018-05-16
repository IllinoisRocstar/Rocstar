      SUBROUTINE read_tetmesh(numbc_prmry,numvertx2d,dhmin,nprocs)

      USE meshdata

      IMPLICIT NONE

      INTEGER :: numbc_prmry
      INTEGER :: numvertx2d
      REAL*8 :: dhmin
      INTEGER :: numbc_prmry_mm

! Local
      REAL*8 :: xx,yy,zz,size1,size2,size3,size4,size5,size6

      INTEGER :: id
      INTEGER :: i, j
      INTEGER :: itype
      INTEGER :: iface
      REAL*8 :: value
      REAL*8 :: press
C-- Tempory holding array of partitioned metis  arrays
      INTEGER, ALLOCATABLE, DIMENSION(:) :: npart
      INTEGER :: edgecut
      INTEGER :: nprocs
      INTEGER :: ip
c
      INTEGER :: iaux,iaux1,iaux2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: imin,nmin,imax,nmax,ninc
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lmelv_prmry_1d
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lmelv_prmry_2d
      REAL*4,  ALLOCATABLE, DIMENSION(:) :: coor_1d
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ipressflg
      INTEGER :: numnp2D, numelv2D
      INTEGER :: nnd1, nnd2, nnd3, nnd4
      INTEGER :: n1, n2, n3, n4
      INTEGER :: iloc

      INTEGER :: ios

      INTEGER :: NBELE,NBELEF,LOELE,LOELEF,IDEB,IFIN
      INTEGER :: NBPOI, LOPOI, NBPOIF, LOPOIF

      INTEGER :: TetForm ! format of Tetmesh input files 0=binary, 1=ascii

      PRINT*,'MESH OPTION:'
      PRINT*,'  READING TETMESH MESH'
c$$$
c$$$      TetForm = 0
c$$$      ipress = 0 ! flag for if we have pressure loading 0 = no 1 = yes
c$$$
c$$$! check to see if it a binary file
c$$$      OPEN(io_input,FILE=prefx(1:prefx_lngth)//'.noboiteb',
c$$$     $     FORM='unformatted',STATUS='old',IOSTAT=ios )
c$$$! If binary file does not exist, check to see that ascii file exist
c$$$      IF(ios.NE.0)THEN
c$$$         OPEN(io_input,FILE=prefx(1:prefx_lngth)//'.noboite',
c$$$     $        FORM='formatted',STATUS='old',IOSTAT=ios )
c$$$         TetForm = 1
c$$$      ENDIF
c$$$! Neither file exist, so there is a problem
c$$$      IF(ios.NE.0)THEN
c$$$         PRINT*,'ERROR: Problem opening ',prefx(1:prefx_lngth)//'.noboite'
c$$$         PRINT*,' .. STOPPING'
c$$$         STOP
c$$$      ENDIF
c$$$
c$$$! Record 1:
c$$$
c$$$      IF(TetForm.EQ.0)THEN
c$$$         PRINT*,' ...... Binary input/output'
c$$$         READ(10) numelv_prmry,numnp_prmry,iaux,iaux,iaux,
c$$$     $        NBELE,LOELE,NBELEF,LOELEF, ! for element records
c$$$     $        NBPOI,LOPOI,NBPOIF,LOPOIF, ! for coordinate records
c$$$     $        iaux,iaux,iaux,iaux ! for sub
c$$$
c$$$      ELSE
c$$$         PRINT*,' ...... ASCII input/output'
c$$$         OPEN(10,file=prefx(1:prefx_lngth)//'.noboite', FORM='formatted')
c$$$! Record 1:
c$$$         READ(10,*) numelv_prmry,numnp_prmry,iaux,iaux,iaux,
c$$$     $        NBELE,LOELE,NBELEF,LOELEF, ! for element records
c$$$     $        NBPOI,LOPOI,NBPOIF,LOPOIF, ! for coordinate records
c$$$     $        iaux,iaux,iaux,iaux ! for sub
c$$$
c$$$      ENDIF
c$$$
c$$$      PRINT*,' '
c$$$
c$$$! Record 2:
c$$$
c$$$      PRINT*,' Number of elements = ', numelv_prmry
c$$$      PRINT*,' Number of nodes = ', numnp_prmry
c$$$
c$$$      ALLOCATE(lmelv_prmry(1:4,1:numelv_prmry))
c$$$C
c$$$C --  Read element connectivity array
c$$$
c$$$      ALLOCATE(lmelv_prmry_1d(1:4*numelv_prmry))
c$$$      IF(TetForm.EQ.0)THEN
c$$$
c$$$         IDEB = 1
c$$$         IFIN = LOELE
c$$$         DO j = 1, NBELE
c$$$            READ(10) (lmelv_prmry_1d(i),i = IDEB, IFIN)
c$$$            IDEB = IDEB + LOELE
c$$$            IFIN = IFIN + LOELE
c$$$         ENDDO
c$$$         IF(NBELEF.NE.0)THEN
c$$$            IFIN = IDEB + LOELEF - 1
c$$$            READ(10) (lmelv_prmry_1d(i),i = IDEB, IFIN)
c$$$         ENDIF
c$$$      ELSE
c$$$         READ(10,*) lmelv_prmry_1d(1:numelv_prmry*4)
c$$$      ENDIF
c$$$      lmelv_prmry = RESHAPE(lmelv_prmry_1d,(/4,numelv_prmry/))
c$$$
c$$$      DEALLOCATE(lmelv_prmry_1d)
c$$$
c$$$      PRINT*,'   Finished element connectivity'
c$$$! Record 3
c$$$
c$$$      ALLOCATE(coor(1:3,1:numnp_prmry))
c$$$      
c$$$      ALLOCATE(coor_1d(3*numnp_prmry))
c$$$      IF(TetForm.EQ.0)THEN
c$$$
c$$$c         READ(10) ( coor_1d(i),i= 1, 3*numnp_prmry)
c$$$
c$$$         IDEB = 1
c$$$         IFIN = LOPOI
c$$$         DO j = 1, NBPOI
c$$$            READ(10) (coor_1d(i),i = IDEB, IFIN)
c$$$            IDEB = IDEB + LOPOI
c$$$            IFIN = IFIN + LOPOI
c$$$         ENDDO
c$$$         IF(NBPOIF.NE.0)THEN
c$$$            IFIN = IDEB + LOPOIF - 1
c$$$            READ(10) (coor_1d(i),i = IDEB, IFIN)
c$$$         ENDIF         
c$$$      ELSE
c$$$         READ(10,*) coor_1d(1:3*numnp_prmry)
c$$$      ENDIF
c$$$      coor = RESHAPE(coor_1d,(/3,numnp_prmry/))
c$$$
c$$$      DEALLOCATE(coor_1d)
c$$$
c$$$      PRINT*,'   Finished Nodal coordinates'
c$$$
c$$$      CLOSE(10)
c$$$
c$$$      ALLOCATE(press_nodal(1:numnp_prmry))
c$$$      press_nodal(:) = 0.d0
c$$$C
c$$$C --  Read element connectivity array
c$$$ 
c$$$      DO i = 1, numelv_prmry
c$$$!
c$$$! -- Find the size of the smallest element
c$$$!
c$$$         xx = coor(1,lmelv_prmry(1,i)) - coor(1,lmelv_prmry(2,i))
c$$$         yy = coor(2,lmelv_prmry(1,i)) - coor(2,lmelv_prmry(2,i))
c$$$         zz = coor(3,lmelv_prmry(1,i)) - coor(3,lmelv_prmry(2,i))
c$$$         size1 = SQRT(xx*xx+yy*yy+zz*zz)
c$$$         xx = coor(1,lmelv_prmry(2,i)) - coor(1,lmelv_prmry(3,i))
c$$$         yy = coor(2,lmelv_prmry(2,i)) - coor(2,lmelv_prmry(3,i))
c$$$         zz = coor(3,lmelv_prmry(2,i)) - coor(3,lmelv_prmry(3,i))
c$$$         size2 = SQRT(xx*xx+yy*yy+zz*zz)
c$$$         xx = coor(1,lmelv_prmry(3,i)) - coor(1,lmelv_prmry(1,i))
c$$$         yy = coor(2,lmelv_prmry(3,i)) - coor(2,lmelv_prmry(1,i))
c$$$         zz = coor(3,lmelv_prmry(3,i)) - coor(3,lmelv_prmry(1,i))
c$$$         size3 = SQRT(xx*xx+yy*yy+zz*zz)
c$$$         xx = coor(1,lmelv_prmry(4,i)) - coor(1,lmelv_prmry(1,i))
c$$$         yy = coor(2,lmelv_prmry(4,i)) - coor(2,lmelv_prmry(1,i))
c$$$         zz = coor(3,lmelv_prmry(4,i)) - coor(3,lmelv_prmry(1,i))
c$$$         size4 = SQRT(xx*xx+yy*yy+zz*zz)
c$$$         xx = coor(1,lmelv_prmry(4,i)) - coor(1,lmelv_prmry(2,i))
c$$$         yy = coor(2,lmelv_prmry(4,i)) - coor(2,lmelv_prmry(2,i))
c$$$         zz = coor(3,lmelv_prmry(4,i)) - coor(3,lmelv_prmry(2,i))
c$$$         size5 = SQRT(xx*xx+yy*yy+zz*zz)
c$$$         xx = coor(1,lmelv_prmry(4,i)) - coor(1,lmelv_prmry(3,i))
c$$$         yy = coor(2,lmelv_prmry(4,i)) - coor(2,lmelv_prmry(3,i))
c$$$         zz = coor(3,lmelv_prmry(4,i)) - coor(3,lmelv_prmry(3,i))
c$$$         size6 = SQRT(xx*xx+yy*yy+zz*zz)
c$$$         dhmin = MIN(size1,size2,size3,size4,size5,size6,dhmin)
c$$$      ENDDO
c$$$
c$$$C -- Partition the finite element mesh
c$$$
c$$$      PRINT*,' REQUESTED NUMBER OF PARTITIONS =',nprocs
c$$$      
c$$$      ALLOCATE(epart(1:numelv_prmry))
c$$$      ALLOCATE(npart(1:numnp_prmry))
c$$$C
c$$$C     Call METIS using the partition the mesh
c$$$C
c$$$      IF(nprocs.GT.1)THEN
c$$$C-ALT       CALL METIS_PartMeshNodal(numel_old_z,nn,
c$$$C-ALT$          elmnts,1,1,nprocs,edgecut,epart_p,npart_p)
c$$$            
c$$$         CALL METIS_PartMeshDual(numelv_prmry,numnp_prmry,
c$$$     $        lmelv_prmry,2,1,nprocs,edgecut,epart,npart)
c$$$      ELSE
c$$$         epart(:) = 1
c$$$      ENDIF
c$$$      DEALLOCATE(npart)
c$$$
c$$$
c$$$C -- create a view of the partitioned mesh
c$$$C -- PMVIS software
c$$$C -- Command: pmvis.sgi.bin -n pmvis.nod -c pmvis.ele -p pmvis.part -o 1
c$$$
c$$$      iopmvis = 0
c$$$
c$$$      IF(iopmvis.EQ.1)THEN
c$$$         
c$$$         OPEN(13,FILE='pmvis.nod')
c$$$         DO i = 1, numnp_prmry
c$$$            WRITE(13,'(3e16.9)') coor(1:3,i)
c$$$         ENDDO
c$$$         CLOSE(13)
c$$$         OPEN(13,FILE='pmvis.ele')
c$$$         DO i = 1, numelv_prmry
c$$$            WRITE(13,'(4i10)') lmelv_prmry(1:4,i)
c$$$         ENDDO
c$$$         CLOSE(13)
c$$$
c$$$         OPEN(13,FILE='pmvis.part')
c$$$         DO i = 1, numelv_prmry
c$$$            WRITE(13,'(1i10)') epart(i)
c$$$         ENDDO
c$$$         CLOSE(13)
c$$$      ENDIF
c$$$
c$$$      OPEN(10,FILE=prefx(1:prefx_lngth)//'.pat', FORM='formatted')
c$$$!
c$$$! - Packet Type 25: Title Card, Format(i2,8i8)
c$$$      READ(10,'()')
c$$$      READ(10,'()')
c$$$
c$$$! -  Packet Type 26: Summary Data, Format(i2,8i8)
c$$$! -      26 ID IV KC N1 N2 N3 N4 N5
c$$$! -  N1 = number of nodes
c$$$! -  N2 = number of elements
c$$$! -  N3 = number of materials
c$$$! -  N4 = number of Element Properties
c$$$! -  N5 = number of Coordinate Frames
c$$$
c$$$      READ(10,*) i,i,i,i,numnp2D,numelv2D
c$$$
c$$$! - part of title card
c$$$      READ(10,'()')
c$$$
c$$$C --  Read Nodal coordinates
c$$$
c$$$      DO i = 1, numnp2D
c$$$         READ(10,'()')
c$$$
c$$$C - format(3e16.9)
c$$$C - caresian coordinate of Nodes x, y, z
c$$$         READ(10,'()') 
c$$$C - ICF GTYPE NDF CONFIG CID PSPC
c$$$         READ(10,'()')
c$$$      ENDDO
c$$$C
c$$$
c$$$C --  Read element connectivity array
c$$$
c$$$
c$$$      ALLOCATE(lmelv_prmry_2d(1:3,1:numelv2D))
c$$$      DO i = 1, numelv2D
c$$$
c$$$C - Type 02: Element Data, format(i2,8i8)
c$$$C        2 ID IV KC N1
c$$$C     ID = element ID
c$$$C     IV = Shape (5 = tet)
c$$$         READ(10,*) j, id
c$$$C     NODES CONFIG CID CEID 
c$$$         READ(10,'()')
c$$$C     LNODES
c$$$C        LNODES = Element corner nodes followed by additional nodes
c$$$         READ(10,*) lmelv_prmry_2d(1:3,id)
c$$$      ENDDO
c$$$
c$$$C
c$$$C
c$$$C -- Continue to read in the rest of the input parameters
c$$$
c$$$      numbc_prmry = 0
c$$$
c$$$      ALLOCATE(ibcflg(1:numnp_prmry))
c$$$      ALLOCATE(ibcflg_mm(1:numnp_prmry))
c$$$      ALLOCATE(ibcaxi(1:numnp_prmry))
c$$$      ibcflg(:) = 0
c$$$      ibcflg_mm(:) = 0
c$$$      ibcaxi(:) = 0
c$$$
c$$$      ALLOCATE(ipressflg(1:numnp_prmry))
c$$$      ipressflg(:) = 0 
c$$$      numel_tri = 0
c$$$      DO
c$$$      READ(10,*) itype,id
c$$$
c$$$      IF(itype.EQ.99) EXIT
c$$$      IF(itype.EQ.4)THEN
c$$$         READ(10,'()')     
c$$$      ELSE IF(itype .EQ. 6)THEN      ! pressure loading
c$$$         ipress = 1
c$$$         READ(10,'()')
c$$$         nnd1 = lmelv_prmry_2d(1,id)
c$$$         nnd2 = lmelv_prmry_2d(2,id)
c$$$         nnd3 = lmelv_prmry_2d(3,id)
c$$$         ipressflg(nnd1) = 1
c$$$         ipressflg(nnd2) = 1
c$$$         ipressflg(nnd3) = 1
c$$$         numel_tri = numel_tri + 1
c$$$C Non-zero load value 
c$$$         READ(10,*) press
c$$$         press = ABS(press)
c$$$
c$$$         IF(press_nodal(nnd1).EQ.0.d0)THEN
c$$$            press_nodal(nnd1) = press
c$$$         ELSE
c$$$            press_nodal(nnd1) = MIN(press,press_nodal(nnd1))
c$$$         ENDIF
c$$$         IF(press_nodal(nnd2).EQ.0.d0)THEN
c$$$            press_nodal(nnd2) = press
c$$$         ELSE
c$$$            press_nodal(nnd2) = MIN(press,press_nodal(nnd2))
c$$$         ENDIF
c$$$         IF(press_nodal(nnd3).EQ.0.d0)THEN
c$$$            press_nodal(nnd3) = press
c$$$         ELSE
c$$$            press_nodal(nnd3) = MIN(press,press_nodal(nnd3))
c$$$         ENDIF
c$$$         press_nodal(nnd1) = press
c$$$         press_nodal(nnd2) = press
c$$$         press_nodal(nnd3) = press
c$$$
c$$$
c$$$      ELSE IF(itype .EQ. 8)THEN ! node displacements
c$$$C 
c$$$C     CID = Coordinate frame ID
c$$$C     ICOMP = 6 displacement component flags (0 or 1)
c$$$C     note the flag is passed to node number one (i.e. fixed displacement)
c$$$
c$$$         READ(10,'()')
c$$$C
c$$$C non-blank displacement components
c$$$C note: important!!! when entering the displacement values leave all other items blank
c$$$         READ(10,*) value
c$$$         IF(value.LT.1000)THEN
c$$$            numbc_prmry = numbc_prmry + 1
c$$$            ibcflg(id) = INT(value)
c$$$         ELSE IF(value.GE.1000)THEN
c$$$            numbc_prmry_mm = numbc_prmry_mm + 1
c$$$            ibcflg_mm(id) = INT(value/1000.d0)
c$$$         ENDIF
c$$$      ELSE IF(itype .EQ. 7)THEN ! mesh motion boundary conditions
c$$$
c$$$         READ(10,'()')
c$$$C 
c$$$C     CID = Coordinate frame ID
c$$$C     ICOMP = 6 displacement component flags (0 or 1)
c$$$C     note the flag is passed to node number one (i.e. fixed displacement)
c$$$         READ(10,'(1e16.9)') value
c$$$         numbc_prmry_mm = numbc_prmry_mm + 1
c$$$C
c$$$C non-blank displacement components
c$$$C note: important!!! when entering the displacement values leave all other items blank
c$$$         ibcflg_mm(id) = INT(value)
c$$$      ELSE IF(itype .EQ. 5)THEN ! coordinate frames
c$$$         READ(10,'()')
c$$$         READ(10,'()')
c$$$         READ(10,'()')
c$$$         READ(10,'()')
c$$$      ELSE IF(itype .EQ. 10)THEN ! nodal temperature (used to flag axisymetric)
c$$$         READ(10,*) value
c$$$         ibcaxi(id) = INT(value)
c$$$      ELSE
c$$$         PRINT*,'packet type',itype,' not accounted for in program'
c$$$         PRINT*,'Source of possible error, STOPPING...'
c$$$         STOP
c$$$      ENDIF
c$$$      ENDDO
c$$$      IF(numvertx.EQ.4)THEN
c$$$         numvertx2d = 3
c$$$      ELSE
c$$$         numvertx2d = 6
c$$$      ENDIF
c$$$      IF(numel_tri.GT.0)THEN
c$$$         ALLOCATE(lmtri(1:numvertx2d,1:numel_tri))
c$$$         ALLOCATE(elm_2D(1:numel_tri))
c$$$         ALLOCATE(elm_2D_flag(1:numel_tri))
c$$$         ALLOCATE(epart_2d(1:numel_tri))
c$$$         elm_2D(:) = 0
c$$$         elm_2D_flag(:) = 0
c$$$         ALLOCATE(numel_2D(1:nprocs))
c$$$         numel_2D(:) = 0
c$$$      ! Associate the surface elements with the volume elements
c$$$      ! because tetmesh does not keep the same element numbering
c$$$         iloc = 0
c$$$         DO i = 1, numelv_prmry
c$$$            nnd1 = lmelv_prmry(1,i)
c$$$            nnd2 = lmelv_prmry(2,i)
c$$$            nnd3 = lmelv_prmry(3,i)
c$$$            nnd4 = lmelv_prmry(4,i)
c$$$            n1 = ipressflg(lmelv_prmry(1,i))
c$$$            n2 = ipressflg(lmelv_prmry(2,i))
c$$$            n3 = ipressflg(lmelv_prmry(3,i))
c$$$            n4 = ipressflg(lmelv_prmry(4,i))
c$$$            id = i
c$$$            IF(n1+n2+n3+n4.GE.3)THEN
c$$$               ip = epart(id)
c$$$               IF(n1.EQ.1.AND.n2.EQ.1.AND.n3.EQ.1)THEN
c$$$                  iloc = iloc + 1
c$$$                  elm_2D(iloc) = id
c$$$                  epart_2d(iloc) = ip
c$$$                  lmtri(1,iloc) = nnd1
c$$$                  lmtri(2,iloc) = nnd3
c$$$                  lmtri(3,iloc) = nnd2
c$$$                  numel_2d(ip) = numel_2d(ip) + 1
c$$$               ENDIF
c$$$               IF(n1.EQ.1.AND.n2.EQ.1.AND.n4.EQ.1)THEN
c$$$                  iloc = iloc + 1
c$$$                  elm_2D(iloc) = id
c$$$                  epart_2d(iloc) = ip
c$$$                  lmtri(1,iloc) = nnd1
c$$$                  lmtri(2,iloc) = nnd2
c$$$                  lmtri(3,iloc) = nnd4
c$$$                  numel_2d(ip) = numel_2d(ip) + 1
c$$$               ENDIF
c$$$               IF(n2.EQ.1.AND.n3.EQ.1.AND.n4.EQ.1)THEN
c$$$                  iloc = iloc + 1
c$$$                  elm_2D(iloc) = id
c$$$                  epart_2d(iloc) = ip
c$$$                  lmtri(1,iloc) = nnd2
c$$$                  lmtri(2,iloc) = nnd3
c$$$                  lmtri(3,iloc) = nnd4
c$$$                  numel_2d(ip) = numel_2d(ip) + 1
c$$$               ENDIF
c$$$               IF(n1.EQ.1.AND.n3.EQ.1.AND.n4.EQ.1)THEN
c$$$                  iloc = iloc + 1
c$$$                  elm_2D(iloc) = id
c$$$                  epart_2d(iloc) = ip
c$$$                  lmtri(1,iloc) = nnd4
c$$$                  lmtri(2,iloc) = nnd3
c$$$                  lmtri(3,iloc) = nnd1
c$$$                  numel_2d(ip) = numel_2d(ip) + 1
c$$$               ENDIF
c$$$            ENDIF
c$$$         ENDDO
c$$$         PRINT*,'iloc count, numel_tri should match', iloc,numel_tri
c$$$      ENDIF
c$$$
c$$$! subroutine
c$$$      IF(numel_tri.GT.0)THEN
c$$$C
c$$$C - find 2D neighbors
c$$$
c$$$         ALLOCATE(neigh_2d(1:3,1:numel_tri))
c$$$         neigh_2d = 0
c$$$         DO i = 1, numel_tri
c$$$            iaux1 = lmtri(2,i)
c$$$            iaux2 = lmtri(3,i)
c$$$            DO j=1,numel_tri
c$$$               IF(j.NE.i)THEN
c$$$                  IF( ANY(lmtri(1:3,j).EQ.iaux1) .AND.
c$$$     $                 ANY(lmtri(1:3,j).EQ.iaux2) )THEN
c$$$                     neigh_2d(1,i)=j
c$$$                     GOTO 100
c$$$                  ENDIF
c$$$               ENDIF
c$$$            ENDDO
c$$$ 100        CONTINUE
c$$$            iaux1 = lmtri(1,i)
c$$$            iaux2 = lmtri(3,i)
c$$$            DO j = 1, numel_tri
c$$$               IF(j.NE.i)THEN
c$$$                  IF( ANY(lmtri(1:3,j).EQ.iaux1) .AND.
c$$$     $                 ANY(lmtri(1:3,j).EQ.iaux2) )THEN
c$$$                     neigh_2d(2,i)=j
c$$$                     GOTO 101
c$$$                  ENDIF
c$$$               ENDIF
c$$$            ENDDO
c$$$ 101        CONTINUE
c$$$            iaux1 = lmtri(1,i)
c$$$            iaux2 = lmtri(2,i)
c$$$            DO j = 1, numel_tri
c$$$               IF(j.NE.i)THEN
c$$$                  IF( ANY(lmtri(1:3,j).EQ.iaux1) .AND.
c$$$     $                 ANY(lmtri(1:3,j).EQ.iaux2) )THEN
c$$$                     neigh_2d(3,i) = j
c$$$                     GOTO 102
c$$$                  ENDIF
c$$$               ENDIF
c$$$            ENDDO
c$$$ 102        CONTINUE
c$$$         ENDDO
c$$$
c$$$!         DO i=1,numel_tri
c$$$!            WRITE(899,*) neigh_2d(1:3,1)
c$$$!         ENDDO
c$$$
c$$$         ALLOCATE(imin(1:numnp_prmry),nmin(1:numel_tri),imax(1:numnp_prmry),nmax(1:numel_tri))
c$$$         ALLOCATE(ninc(1:numel_tri))
c$$$
c$$$         CALL STRUCT_tri(numnp_prmry,numel_tri,lmtri,neigh_2d,IMIN,NMIN,IMAX,NMAX,NINC)
c$$$
c$$$         DEALLOCATE(imin,nmin,imax,nmax,ninc)
c$$$
c$$$!         DO i=1,numel_tri
c$$$!            WRITE(900,*) neigh_2d(1:3,i)
c$$$!         ENDDO
c$$$      ENDIF         
c$$$      DEALLOCATE(ipressflg)
c$$$      CLOSE(10)
c$$$
c$$$      PRINT*,' Number of faces load by a distributed load =', numel_2d
c$$$      PRINT*,' Number of nodes with a boundary flag =',  numbc_prmry
c$$$C
c$$$C --  READ THE ELEMENT NEIGHBOR ARRAY
c$$$      IF(IOformat.EQ.0)THEN
c$$$         OPEN(10,FILE=prefx(1:prefx_lngth)//'.neigh', FORM='unformatted')
c$$$         
c$$$         ALLOCATE(neigh(1:4,1:numelv_prmry))
c$$$         
c$$$         neigh(1:4,1:numelv_prmry) = -1
c$$$         DO i = 1, numelv_prmry
c$$$            READ(10) neigh(1:4,i)
c$$$         ENDDO
c$$$         CLOSE(10)
c$$$      ELSE
c$$$         OPEN(10,FILE=prefx(1:prefx_lngth)//'.neigh', FORM='formatted')
c$$$         
c$$$         ALLOCATE(neigh(1:4,1:numelv_prmry))
c$$$         
c$$$         neigh(1:4,1:numelv_prmry) = -1
c$$$         DO i = 1, numelv_prmry
c$$$            READ(10,*) neigh(1:4,i)
c$$$         ENDDO
c$$$         CLOSE(10)
c$$$      ENDIF

      END
    
