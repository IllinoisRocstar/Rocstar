      SUBROUTINE read_patran(numbc_prmry,numvertx2d,dhmin,nprocs)

      USE meshdata

      IMPLICIT NONE

      INTEGER :: numbc_prmry
      INTEGER :: numvertx2d
      REAL*8 :: dhmin
      INTEGER :: numbc_prmry_mm

! Local
      REAL*8 :: xx,yy,zz,size1,size2,size3,size4,size5,size6
      REAL*8 :: dt_courant

      INTEGER :: id
      INTEGER :: i, j
      INTEGER :: itype
      INTEGER :: n1,n2,n3,n4
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
      INTEGER, ALLOCATABLE, DIMENSION(:) :: imap
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lm4node
      INTEGER :: nnd1, nnd2, nnd3, nnd4, nnd5, nnd6, nnd7, nnd8, nnd9, nnd10
      INTEGER :: numnp4node

      integer :: iFaceId, itmp

      integer :: NumSFtri
      integer :: NumStri

      PRINT*,'MESH OPTION:'
      PRINT*,'  READING PATRAN MESH'
      PRINT*,' '

      itype6 = 0
      itype7 = 0
      itype8 = 0
      itmp = 0

      NumSFtri = 0
      NumStri = 0
      NumNdHistory = 0

      ipress = 0 ! flag for if we have pressure loading 0 = no 1 = yes

      OPEN(10,file=prefx(1:prefx_lngth)//'.pat', FORM='formatted')
C
C - Packet Type 25: Title Card, Format(i2,8i8)
      READ(10,'()')
      READ(10,'()')

C -  Packet Type 26: Summary Data, Format(i2,8i8)
C -      26 ID IV KC N1 N2 N3 N4 N5
C -  N1 = number of nodes
C -  N2 = number of elements
C -  N3 = number of materials
C -  N4 = number of Element Properties
C -  N5 = number of Coordinate Frames

      READ(10,*) i,i,i,i,numnp_prmry,numelv_prmry,NumMat

      IF(NumMat.LE.0) NumMat = 1

C - part of title card
      READ(10,'()')

      PRINT*,' Number of elements = ', numelv_prmry
      PRINT*,' Number of nodes = ', numnp_prmry
      PRINT*,' Number of Materials = ', NumMat
      PRINT*,' Number of vertex=',numvertx
C
C --  Read Nodal coordinates

      ALLOCATE(coor(1:3,1:numnp_prmry),lmelv_prmry(1:numvertx,1:numelv_prmry))
      ALLOCATE(press_nodal(1:numnp_prmry))
      
      IF(numvertx.NE.4)THEN
         ALLOCATE(lm4node( 1:4,1:numelv_prmry))
         ALLOCATE(imap(1:numnp_prmry))
         imap(:) = 0
         numnp4node = 0
      ELSE
         numnp4node = numnp_prmry
      ENDIF

      press_nodal(:) = 0.d0

!      ConvertUnit = .0254d0 ! for titan, units in engles
      
      DO i = 1, numnp_prmry

C - Type 01: Node Data,  format(i2,8i8)
C        1 ID IV KC
C    ID = node id
C    IV = 0 n/a
C    KC = 2
         READ(10,*) j,id
C - format(3e16.9)
C - caresian coordinate of Nodes x, y, z
         READ(10,'(3e16.9)') coor(1:3,id)
         coor(1:3,id) = coor(1:3,id)*ConvertUnit
C - ICF GTYPE NDF CONFIG CID PSPC
         READ(10,'()')
      ENDDO
C
C --  Read element connectivity array

      ALLOCATE(MatId(1:numelv_prmry))
 
      DO i = 1, numelv_prmry

C - Type 02: Element Data, format(i2,8i8)
C        2 ID IV KC N1
C     ID = element ID
C     IV = Shape (5 = tet)
         READ(10,*) j, id
C     NODES CONFIG CID CEID 
         READ(10,'(3i8)') numvertx,j, MatId(i)
         IF(MatId(i).LE.0) MatId(i) = 1
C     LNODES
C        LNODES = Element corner nodes followed by additional nodes
         READ(10,*) lmelv_prmry(1:numvertx,id)

         IF(numvertx.EQ.10)THEN
            DO j = 1, 4
               IF(imap(lmelv_prmry(j,id)).EQ.0)THEN
                  numnp4node = numnp4node + 1
                  imap(lmelv_prmry(j,id)) = numnp4node
               ENDIF
               lm4node(j,i) = imap(lmelv_prmry(j,id))
            ENDDO
         ENDIF

!
! -- Find the size of the smallest element
!
         xx = coor(1,lmelv_prmry(1,id)) - coor(1,lmelv_prmry(2,id))
         yy = coor(2,lmelv_prmry(1,id)) - coor(2,lmelv_prmry(2,id))
         zz = coor(3,lmelv_prmry(1,id)) - coor(3,lmelv_prmry(2,id))
         size1 = SQRT(xx*xx+yy*yy+zz*zz)
         xx = coor(1,lmelv_prmry(2,id)) - coor(1,lmelv_prmry(3,id))
         yy = coor(2,lmelv_prmry(2,id)) - coor(2,lmelv_prmry(3,id))
         zz = coor(3,lmelv_prmry(2,id)) - coor(3,lmelv_prmry(3,id))
         size2 = SQRT(xx*xx+yy*yy+zz*zz)
         xx = coor(1,lmelv_prmry(3,id)) - coor(1,lmelv_prmry(1,id))
         yy = coor(2,lmelv_prmry(3,id)) - coor(2,lmelv_prmry(1,id))
         zz = coor(3,lmelv_prmry(3,id)) - coor(3,lmelv_prmry(1,id))
         size3 = SQRT(xx*xx+yy*yy+zz*zz)
         xx = coor(1,lmelv_prmry(4,id)) - coor(1,lmelv_prmry(1,id))
         yy = coor(2,lmelv_prmry(4,id)) - coor(2,lmelv_prmry(1,id))
         zz = coor(3,lmelv_prmry(4,id)) - coor(3,lmelv_prmry(1,id))
         size4 = SQRT(xx*xx+yy*yy+zz*zz)
         xx = coor(1,lmelv_prmry(4,id)) - coor(1,lmelv_prmry(2,id))
         yy = coor(2,lmelv_prmry(4,id)) - coor(2,lmelv_prmry(2,id))
         zz = coor(3,lmelv_prmry(4,id)) - coor(3,lmelv_prmry(2,id))
         size5 = SQRT(xx*xx+yy*yy+zz*zz)
         xx = coor(1,lmelv_prmry(4,id)) - coor(1,lmelv_prmry(3,id))
         yy = coor(2,lmelv_prmry(4,id)) - coor(2,lmelv_prmry(3,id))
         zz = coor(3,lmelv_prmry(4,id)) - coor(3,lmelv_prmry(3,id))
         size6 = SQRT(xx*xx+yy*yy+zz*zz)
         dhmin = MIN(size1,size2,size3,size4,size5,size6,dhmin)
      ENDDO

      IF(numvertx.NE.4) DEALLOCATE(imap)

C -- Partition the finite element mesh

      PRINT*,' REQUESTED NUMBER OF PARTITIONS =',nprocs

      ALLOCATE(epart(1:numelv_prmry))
      ALLOCATE(npart(1:numnp4node))
C
C     Call METIS using the partition the mesh
C
      PRINT*,'CALLING METIS'
      IF(nprocs.GT.1)THEN
C-ALT       CALL METIS_PartMeshNodal(numel_old_z,nn,
C-ALT$          elmnts,1,1,nprocs,edgecut,epart_p,npart_p)
         IF(numvertx.EQ.10)THEN
            CALL METIS_PartMeshDual(numelv_prmry,numnp4node,
     $           lm4node,2,1,nprocs,edgecut,epart,npart)
            DEALLOCATE(lm4node)
         ELSE
            CALL METIS_PartMeshDual(numelv_prmry,numnp_prmry,
     $           lmelv_prmry,2,1,nprocs,edgecut,epart,npart)
         ENDIF
      ELSE
         epart(:) = 1
      ENDIF
      DEALLOCATE(npart)

      PRINT*,' Called METIS'

C -- create a view of the partitioned mesh
C -- PMVIS software
C -- Command: pmvis.sgi.bin -n pmvis.nod -c pmvis.ele -p pmvis.part -o 1

      iopmvis = 0
      
      IF(iopmvis.EQ.1)THEN

         OPEN(13,FILE='pmvis.nod')
         DO i = 1, numnp_prmry
            WRITE(13,'(3e16.9)') coor(1:3,i)
         ENDDO
         CLOSE(13)
         OPEN(13,FILE='pmvis.ele')
         DO i = 1, numelv_prmry
            WRITE(13,'(4i10)') lmelv_prmry(1:4,i)
         ENDDO
         CLOSE(13)

         OPEN(13,FILE='pmvis.part')
         DO i = 1, numelv_prmry
            WRITE(13,'(1i10)') epart(i)
         ENDDO
         CLOSE(13)
      ENDIF
C
C
C -- Continue to read in the rest of the input parameters

      numbc_prmry = 0

      ALLOCATE(numel_2d(1:nprocs,1:2))
      numel_2d(:,:) = 0

      ALLOCATE(ibcflg(1:numnp_prmry))
      ALLOCATE(ibcflg_mm(1:numnp_prmry))
      ALLOCATE(ibcaxi(1:numnp_prmry))
      ALLOCATE(iNdsBurnFlg(1:numnp_prmry))
      iNdsBurnFlg(:) = 0
      ibcflg(:) = 0
      ibcflg_mm(:) = 0
      ibcaxi(:) = 0
      
      PRINT*,'FINISH READING PATRAN FILE'
      numel_tri = 0
      DO
      READ(10,*) itype,id
      IF(itype.EQ.99) EXIT
      IF(itype.EQ.4)THEN
         READ(10,'()')

! Marks if IntFaceFlag :
!
!   SolidFluid Interface, Non-burning = -1
!   SolidFluid Interface, burning = 1
!   NotASolidFluid Interface = 0

      ELSE IF(itype .EQ. 6)THEN      ! pressure loading
c$$$         IF(itype6.EQ.0)THEN
c$$$            Allocate(MeshLT0(1:nprocs))
c$$$            Allocate(MeshGT0(1:nprocs))
c$$$            Allocate(MeshEQ0(1:nprocs))
c$$$            ALLOCATE(Mesh2dItem)
c$$$            itype6 = 1
c$$$         ENDIF
         
         ipress = 1
         ip = epart(id)
         READ(10,'(i1,i1,i1,i6,4i1)')i,i,i,i,n1,n2,n3,n4
         READ(10,*) press
         IntFaceFlag = INT(press)

         IF(IntFaceFlag.NE.0)THEN
            IntFaceFlag = 1     ! S/F interface
            NumSFtri = NumSFtri + 1
         ELSE
            IntFaceFlag = 2     ! S interface
            NumStri = NumStri + 1
         ENDIF


         numel_tri = numel_tri + 1
         numel_2d(ip,IntFaceFlag) = numel_2d(ip,IntFaceFlag) + 1
         
         IF(n1.EQ.1.AND.n2.EQ.1.AND.n3.EQ.1)THEN
            IF(numvertx.EQ.4)THEN
               nnd1 = lmelv_prmry(1,id)
               nnd2 = lmelv_prmry(3,id)
               nnd3 = lmelv_prmry(2,id)
               itmp = itmp + 1
               iFaceID = 1
            ELSE
               nnd1 = lmelv_prmry(1,id)
               nnd2 = lmelv_prmry(3,id)
               nnd3 = lmelv_prmry(2,id)
               nnd4 = lmelv_prmry(7,id)
               nnd5 = lmelv_prmry(6,id)
               nnd6 = lmelv_prmry(5,id)
               iFaceID = 1
            ENDIF
         ELSE IF(n1.EQ.1.AND.n2.EQ.1.AND.n4.EQ.1)THEN
            IF(numvertx.EQ.4)THEN
               nnd1 = lmelv_prmry(1,id)
               nnd2 = lmelv_prmry(2,id)
               nnd3 = lmelv_prmry(4,id)
               itmp = itmp + 1
               iFaceID = 2
            ELSE
               nnd1 = lmelv_prmry(1,id)
               nnd2 = lmelv_prmry(2,id)
               nnd3 = lmelv_prmry(4,id)
               nnd4 = lmelv_prmry(5,id)
               nnd5 = lmelv_prmry(9,id)
               nnd6 = lmelv_prmry(8,id)
               iFaceID = 2
            ENDIF
         ELSE IF(n2.EQ.1.AND.n3.EQ.1.AND.n4.EQ.1)THEN
            IF(numvertx.EQ.4)THEN
               nnd1 = lmelv_prmry(2,id)
               nnd2 = lmelv_prmry(3,id)
               nnd3 = lmelv_prmry(4,id)
               itmp = itmp + 1
               iFaceID = 3
            ELSE
               nnd1 = lmelv_prmry(2,id)
               nnd2 = lmelv_prmry(3,id)
               nnd3 = lmelv_prmry(4,id)
               nnd4 = lmelv_prmry(6,id)
               nnd5 = lmelv_prmry(10,id)
               nnd6  = lmelv_prmry(9,id)
               iFaceID = 3
            ENDIF
         ELSE IF(n1.EQ.1.AND.n3.EQ.1.AND.n4.EQ.1)THEN
            IF(numvertx.EQ.4)THEN
               nnd1 = lmelv_prmry(4,id)
               nnd2 = lmelv_prmry(3,id)
               nnd3 = lmelv_prmry(1,id)
               itmp = itmp + 1
               iFaceID = 4
            ELSE
               nnd1 = lmelv_prmry(4,id)
               nnd2 = lmelv_prmry(3,id)
               nnd3 = lmelv_prmry(1,id)
               nnd4 = lmelv_prmry(10,id)
               nnd5 = lmelv_prmry(7,id)
               nnd6 = lmelv_prmry(8,id)
               iFaceID = 4
            ENDIF
         ELSE
            PRINT*,'Error in pressure face numbering'
            STOP
         ENDIF
         IF(numvertx.EQ.4)THEN
            WRITE(301,'(6i10)') ip,nnd1,nnd2,nnd3,id,IntFaceFlag
            IF(IntFaceFlag.EQ.1)THEN
               IF(iNdsBurnFlg(nnd1).EQ.0) iNdsBurnFlg(nnd1) = 1
               IF(iNdsBurnFlg(nnd2).EQ.0) iNdsBurnFlg(nnd2) = 1
               IF(iNdsBurnFlg(nnd3).EQ.0) iNdsBurnFlg(nnd3) = 1
            ENDIF
         ELSE
            WRITE(301,'(9i10)') ip,nnd1,nnd2,nnd3,nnd4,nnd5,nnd6,id,IntFaceFlag
            IF(IntFaceFlag.EQ.1)THEN
               IF(iNdsBurnFlg(nnd1).EQ.0) iNdsBurnFlg(nnd1) = 1
               IF(iNdsBurnFlg(nnd2).EQ.0) iNdsBurnFlg(nnd2) = 1
               IF(iNdsBurnFlg(nnd3).EQ.0) iNdsBurnFlg(nnd3) = 1
               IF(iNdsBurnFlg(nnd4).EQ.0) iNdsBurnFlg(nnd4) = 1
               IF(iNdsBurnFlg(nnd5).EQ.0) iNdsBurnFlg(nnd5) = 1
               IF(iNdsBurnFlg(nnd6).EQ.0) iNdsBurnFlg(nnd6) = 1
            ENDIF
         ENDIF

      ELSE IF(itype .EQ. 8)THEN ! node displacements
C 
C     CID = Coordinate frame ID
C     ICOMP = 6 displacement component flags (0 or 1)
C     note the flag is passed to node number one (i.e. fixed displacement)

         READ(10,'()')
         READ(10,*) value
         numbc_prmry = numbc_prmry + 1
C
C non-blank displacement components
C note: important!!! when entering the displacement values leave all other items blank

         ibcflg(id) = INT(value)
         
      ELSE IF(itype .EQ. 7)THEN ! mesh motion boundary conditions
C 
C     CID = Coordinate frame ID
C     ICOMP = 6 displacement component flags (0 or 1)
C     note the flag is passed to node number one (i.e. fixed displacement)

         READ(10,'()')
         READ(10,*) value
         numbc_prmry_mm = numbc_prmry_mm + 1
C
C non-blank displacement components
C note: important!!! when entering the displacement values leave all other items blank
         ibcflg_mm(id) = INT(value)

      ELSE IF(itype .EQ. 5)THEN ! coordinate frames
         READ(10,'()')
         READ(10,'()')
         READ(10,'()')
         READ(10,'()')
      ELSE IF(itype .EQ. 10)THEN ! nodal temperature (corners)
         READ(10,*) value
         NumNdHistory = NumNdHistory + 1
         NdHistory(NumNdHistory) = id
         print*,'Monitoring Node =', id
      ELSE
         PRINT*,'packet type',itype,' not accounted for in program'
         PRINT*,'Source of possable error, STOPPING...'
         STOP
      ENDIF
      ENDDO
      CLOSE(301)

      PRINT*,'finished reading patran input'
      PRINT*,'  Number of SF element =', NumSFtri,itmp
      PRINT*,'  Number of S element =', NumStri


      IF(numvertx.EQ.4)THEN
         numvertx2d = 3
      ELSE
         numvertx2d = 6
      ENDIF
      IF(numel_tri.GT.0)THEN
         ALLOCATE(lmtri(1:numvertx2d,1:numel_tri))
         ALLOCATE(elm_2D(1:numel_tri))
         ALLOCATE(elm_2D_flag(1:numel_tri))
         ALLOCATE(epart_2d(1:numel_tri))
         elm_2D(:) = 0
         elm_2D_flag(:) = 0
         OPEN(302,FILE='fort.301')
         REWIND(302)
         DO i = 1, numel_tri
            READ(302,*) ip,lmtri(:,i),iaux1,iaux2 ! iaux1 is the volumetic elment id
               
            elm_2D(i) = iaux1
            epart_2d(i) = ip
            elm_2D_flag(i) = iaux2
         ENDDO
         CLOSE(302)
c$$$C
c$$$C - find 2D neighbors
c$$$
c$$$         PRINT*,'Finding 2D neighbors'
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
c$$$         ALLOCATE(imin(1:numnp_prmry),nmin(1:numel_tri),imax(1:numnp_prmry),nmax(1:numel_tri))
c$$$         ALLOCATE(ninc(1:numel_tri))
c$$$
c$$$         PRINT*,'Calling Struct_tri'
c$$$
c$$$         CALL STRUCT_tri(numnp_prmry,numel_tri,lmtri,neigh_2d,IMIN,NMIN,IMAX,NMAX,NINC)
c$$$
c$$$         DEALLOCATE(imin,nmin,imax,nmax,ninc)

!        DO i=1,numel_tri
!           WRITE(900,*) neigh_2d(1:3,i)
!        ENDDO
      ENDIF


!      PRINT*,' Number of faces load by a distributed load =', numel_2d
      PRINT*,' Number of nodes with a boundary flag =',  numbc_prmry

      CLOSE(10)
C
C --  READ THE ELEMENT NEIGHBOR ARRAY

      OPEN(10,FILE=prefx(1:prefx_lngth)//'.neigh', FORM='formatted')

      ALLOCATE(neigh(1:4,1:numelv_prmry))

      neigh(1:4,1:numelv_prmry) = -1
      DO i = 1, numelv_prmry
         READ(10,*) neigh(1:4,i)
      ENDDO
      CLOSE(10)

      RETURN

      END
    
