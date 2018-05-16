      SUBROUTINE read_patran_cohin(numbc_prmry,numvertx2d,dhmin,nprocs)

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
      INTEGER, DIMENSION(:,:),  ALLOCATABLE :: v_numnp
C-- Tempory holding array of partitioned metis  arrays
      INTEGER, ALLOCATABLE, DIMENSION(:) :: npart
      INTEGER :: edgecut
      INTEGER :: nprocs
      INTEGER :: ip
c
      INTEGER :: iaux,iaux1,iaux2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: imin,nmin,imax,nmax,ninc

      INTEGER :: nnd1, nnd2, nnd3, nnd4, nnd5, nnd6, nnd7, nnd8, nnd9, nnd10
      REAL*8 :: nface

      PRINT*,'MESH OPTION:'
      PRINT*,'  READING PATRAN MESH'
      PRINT*,' '
c$$$
c$$$      OPEN(10,file=prefx(1:prefx_lngth)//'.pat', FORM='formatted')
c$$$
c$$$      OPEN(12,file=prefx(1:prefx_lngth)//'.msh', FORM='formatted')
c$$$
c$$$C
c$$$C - Packet Type 25: Title Card, Format(i2,8i8)
c$$$      READ(10,'()')
c$$$      READ(10,'()')
c$$$
c$$$C -  Packet Type 26: Summary Data, Format(i2,8i8)
c$$$C -      26 ID IV KC N1 N2 N3 N4 N5
c$$$C -  N1 = number of nodes
c$$$C -  N2 = number of elements
c$$$C -  N3 = number of materials
c$$$C -  N4 = number of Element Properties
c$$$C -  N5 = number of Coordinate Frames
c$$$
c$$$      READ(10,*) i,i,i,i,numnp_prmry,numelv_prmry, i, nummat
c$$$      WRITE(12,'(2i10)') numelv_prmry, numnp_prmry
c$$$C - part of title card
c$$$      READ(10,'()')
c$$$
c$$$      PRINT*,' Number of elements = ', numelv_prmry
c$$$      PRINT*,' Number of nodes = ', numnp_prmry
c$$$      PRINT*,' Number of materials ', nummat
c$$$C
c$$$C --  Read Nodal coordinates
c$$$
c$$$      ALLOCATE(coor(1:3,1:numnp_prmry),lmelv_prmry(1:numvertx,1:numelv_prmry))
c$$$      ALLOCATE(press_nodal(1:numnp_prmry))
c$$$      press_nodal(:) = 0.d0
c$$$
c$$$!      ConvertUnit = .0254d0 ! for titan, units in engles
c$$$      
c$$$      DO i = 1, numnp_prmry
c$$$
c$$$C - Type 01: Node Data,  format(i2,8i8)
c$$$C        1 ID IV KC
c$$$C    ID = node id
c$$$C    IV = 0 n/a
c$$$C    KC = 2
c$$$         READ(10,*) j,id
c$$$C - format(3e16.9)
c$$$C - caresian coordinate of Nodes x, y, z
c$$$         READ(10,'(3e16.9)') coor(1:3,id)
c$$$         coor(1:3,id) = coor(1:3,id)*ConvertUnit
c$$$         WRITE(12,'(i10,3f16.9)') id,coor(1:3,id)
c$$$C - ICF GTYPE NDF CONFIG CID PSPC
c$$$         READ(10,'()')
c$$$      ENDDO
c$$$C
c$$$C --  Read element connectivity array
c$$$ 
c$$$      DO i = 1, numelv_prmry
c$$$
c$$$C - Type 02: Element Data, format(i2,8i8)
c$$$C        2 ID IV KC N1
c$$$C     ID = element ID
c$$$C     IV = Shape (5 = tet)
c$$$         READ(10,*) j, id
c$$$C     NODES CONFIG CID CEID 
c$$$         READ(10,'(3i8)') numvertx,j,matid
c$$$C     LNODES
c$$$C        LNODES = Element corner nodes followed by additional nodes
c$$$         READ(10,*) lmelv_prmry(1:numvertx,id)
c$$$!
c$$$! -- Find the size of the smallest element
c$$$!
c$$$         xx = coor(1,lmelv_prmry(1,id)) - coor(1,lmelv_prmry(2,id))
c$$$         yy = coor(2,lmelv_prmry(1,id)) - coor(2,lmelv_prmry(2,id))
c$$$         zz = coor(3,lmelv_prmry(1,id)) - coor(3,lmelv_prmry(2,id))
c$$$         size1 = SQRT(xx*xx+yy*yy+zz*zz)
c$$$         xx = coor(1,lmelv_prmry(2,id)) - coor(1,lmelv_prmry(3,id))
c$$$         yy = coor(2,lmelv_prmry(2,id)) - coor(2,lmelv_prmry(3,id))
c$$$         zz = coor(3,lmelv_prmry(2,id)) - coor(3,lmelv_prmry(3,id))
c$$$         size2 = SQRT(xx*xx+yy*yy+zz*zz)
c$$$         xx = coor(1,lmelv_prmry(3,id)) - coor(1,lmelv_prmry(1,id))
c$$$         yy = coor(2,lmelv_prmry(3,id)) - coor(2,lmelv_prmry(1,id))
c$$$         zz = coor(3,lmelv_prmry(3,id)) - coor(3,lmelv_prmry(1,id))
c$$$         size3 = SQRT(xx*xx+yy*yy+zz*zz)
c$$$         xx = coor(1,lmelv_prmry(4,id)) - coor(1,lmelv_prmry(1,id))
c$$$         yy = coor(2,lmelv_prmry(4,id)) - coor(2,lmelv_prmry(1,id))
c$$$         zz = coor(3,lmelv_prmry(4,id)) - coor(3,lmelv_prmry(1,id))
c$$$         size4 = SQRT(xx*xx+yy*yy+zz*zz)
c$$$         xx = coor(1,lmelv_prmry(4,id)) - coor(1,lmelv_prmry(2,id))
c$$$         yy = coor(2,lmelv_prmry(4,id)) - coor(2,lmelv_prmry(2,id))
c$$$         zz = coor(3,lmelv_prmry(4,id)) - coor(3,lmelv_prmry(2,id))
c$$$         size5 = SQRT(xx*xx+yy*yy+zz*zz)
c$$$         xx = coor(1,lmelv_prmry(4,id)) - coor(1,lmelv_prmry(3,id))
c$$$         yy = coor(2,lmelv_prmry(4,id)) - coor(2,lmelv_prmry(3,id))
c$$$         zz = coor(3,lmelv_prmry(4,id)) - coor(3,lmelv_prmry(3,id))
c$$$         size6 = SQRT(xx*xx+yy*yy+zz*zz)
c$$$         dhmin = MIN(size1,size2,size3,size4,size5,size6,dhmin)
c$$$
c$$$         WRITE(12,'(7i10)') id, 0, matid, lmelv_prmry(1:numvertx,id)
c$$$
c$$$      ENDDO
c$$$
c$$$      PRINT*,'Minimum Sized element =', dhmin
c$$$      dt_courant = dhmin/cd_fastest
c$$$      PRINT*,'Courant condition element =', dt_courant
c$$$
c$$$C -- Partition the finite element mesh
c$$$
c$$$!      PRINT*,'Enter the number of processors'
c$$$        nprocs = 1
c$$$!       READ*, nprocs 
c$$$
c$$$      ALLOCATE(epart(1:numelv_prmry))
c$$$cc      ALLOCATE(npart(1:numnp_prmry))
c$$$C
c$$$C     Call METIS using the partition the mesh
c$$$C
c$$$cc      IF(nprocs.GT.1)THEN
c$$$C-ALT       CALL METIS_PartMeshNodal(numel_old_z,nn,
c$$$C-ALT$          elmnts,1,1,nprocs,edgecut,epart_p,npart_p)
c$$$cc            
c$$$cc         CALL METIS_PartMeshDual(numelv_prmry,numnp_prmry,
c$$$cc     $        lmelv_prmry,2,1,nprocs,edgecut,epart,npart)
c$$$cc      ELSE
c$$$      epart(:) = 1
c$$$cc      ENDIF
c$$$cc      DEALLOCATE(npart)
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
c$$$C
c$$$C
c$$$C -- Continue to read in the rest of the input parameters
c$$$
c$$$      numbc_prmry = 0
c$$$      nface = 0
c$$$
c$$$      ALLOCATE(numel_2d(1:nprocs))
c$$$      numel_2d(:) = 0
c$$$      ALLOCATE(numnp_2d(1:nprocs))
c$$$      numnp_2d(:) = 0
c$$$      ALLOCATE(v_numnp(1:nprocs,1:numnp_prmry))
c$$$      v_numnp(:,:) = 0
c$$$
c$$$      ALLOCATE(ibcflg(1:numnp_prmry))
c$$$      ALLOCATE(ibcflg_mm(1:numnp_prmry))
c$$$      ALLOCATE(ibcaxi(1:numnp_prmry))
c$$$      ibcflg(:) = 0
c$$$      ibcflg_mm(:) = 0
c$$$      ibcaxi(:) = 0
c$$$      DO
c$$$      READ(10,*) itype,id
c$$$
c$$$      IF(itype.EQ.99) EXIT
c$$$      IF(itype.EQ.4)THEN
c$$$         READ(10,'()')         
c$$$      ELSE IF(itype .EQ. 6)THEN  ! pressure loading
c$$$         ipress = 1
c$$$         ip = epart(id)
c$$$         READ(10,'(i1,i1,i1,i6,4i1)')i,i,i,i,n1,n2,n3,n4
c$$$
c$$$         numel_2d(ip) = numel_2d(ip) + 1
c$$$
c$$$         nface = nface + 1
c$$$         
c$$$         IF(n1.EQ.1.AND.n2.EQ.1.AND.n3.EQ.1)THEN
c$$$            iface = 1
c$$$            IF(numvertx.EQ.4)THEN
c$$$               nnd1 = lmelv_prmry(1,id)
c$$$               nnd2 = lmelv_prmry(2,id)
c$$$               nnd3 = lmelv_prmry(3,id)
c$$$            
c$$$c-plate        if(coor(3,nnd1).GT.0.d0)THEN
c$$$               v_numnp(ip,nnd1) = 1
c$$$               v_numnp(ip,nnd2) = 1
c$$$               v_numnp(ip,nnd3) = 1
c$$$               WRITE(301,'(6i10)') ip, iface,nnd1,nnd3,nnd2,id
c$$$c-plate        endif
c$$$               WRITE(501,'(2i8,x,1a3,x,10i8)') numel_2d(ip),1,'tri', lmelv_prmry(1:3,id)
c$$$            ELSE
c$$$               nnd1 = lmelv_prmry(1,id)
c$$$               nnd2 = lmelv_prmry(2,id)
c$$$               nnd3 = lmelv_prmry(3,id)
c$$$               nnd5 = lmelv_prmry(5,id)
c$$$               nnd6 = lmelv_prmry(6,id)
c$$$               nnd7 = lmelv_prmry(7,id)
c$$$
c$$$               v_numnp(ip,nnd1) = 1
c$$$               v_numnp(ip,nnd2) = 1
c$$$               v_numnp(ip,nnd3) = 1
c$$$               v_numnp(ip,nnd5) = 1
c$$$               v_numnp(ip,nnd6) = 1
c$$$               v_numnp(ip,nnd7) = 1
c$$$
c$$$               WRITE(301,'(9i10)') ip, iface,nnd1,nnd3,nnd2,nnd7,nnd6,nnd5,id
c$$$
c$$$            ENDIF
c$$$         ELSE IF(n1.EQ.1.AND.n2.EQ.1.AND.n4.EQ.1)THEN
c$$$            iface = 2
c$$$            IF(numvertx.EQ.4)THEN
c$$$               nnd1 = lmelv_prmry(1,id)
c$$$               nnd2 = lmelv_prmry(2,id)
c$$$               nnd4 = lmelv_prmry(4,id)
c$$$c-plate        if(coor(3,nnd1).GT.0.d0)THEN
c$$$               v_numnp(ip,nnd1) = 1
c$$$               v_numnp(ip,nnd2) = 1
c$$$               v_numnp(ip,nnd4) = 1
c$$$               WRITE(301,'(6i10)') ip,iface,nnd1,nnd2,nnd4,id
c$$$c-plate        ENDIF
c$$$               WRITE(501,'(2i8,x,1a3,x,10i8)') numel_2d(ip),1,'tri',lmelv_prmry(1:2,id),lmelv_prmry(4,id)
c$$$
c$$$            ELSE
c$$$               nnd1 = lmelv_prmry(1,id)
c$$$               nnd2 = lmelv_prmry(2,id)
c$$$               nnd4 = lmelv_prmry(4,id)
c$$$               nnd5 = lmelv_prmry(5,id)
c$$$               nnd9 = lmelv_prmry(9,id)
c$$$               nnd8 = lmelv_prmry(8,id)
c$$$
c$$$               v_numnp(ip,nnd1) = 1
c$$$               v_numnp(ip,nnd2) = 1
c$$$               v_numnp(ip,nnd4) = 1
c$$$               v_numnp(ip,nnd5) = 1
c$$$               v_numnp(ip,nnd9) = 1
c$$$               v_numnp(ip,nnd8) = 1
c$$$
c$$$               WRITE(301,'(9i10)') ip,iface,nnd1,nnd2,nnd4,nnd5,nnd9,nnd8,id
c$$$
c$$$            ENDIF
c$$$
c$$$         ELSE IF(n2.EQ.1.AND.n3.EQ.1.AND.n4.EQ.1)THEN
c$$$            iface = 3
c$$$            IF(numvertx.EQ.4)THEN
c$$$ 
c$$$               nnd2 = lmelv_prmry(2,id)
c$$$               nnd3 = lmelv_prmry(3,id)
c$$$               nnd4 = lmelv_prmry(4,id)
c$$$c-plate        if(coor(3,nnd2).GT.0.d0)THEN
c$$$               v_numnp(ip,nnd2) = 1
c$$$               v_numnp(ip,nnd3) = 1
c$$$               v_numnp(ip,nnd4) = 1
c$$$               WRITE(301,'(6i10)') ip, iface,nnd2,nnd3,nnd4,id
c$$$c-plate        endif
c$$$               WRITE(501,'(2i8,x,1a3,x,10i8)') numel_2d(ip),1,'tri',lmelv_prmry(2:3,id),lmelv_prmry(4,id)
c$$$            ELSE
c$$$               nnd2 = lmelv_prmry(2,id)
c$$$               nnd3 = lmelv_prmry(3,id)
c$$$               nnd4 = lmelv_prmry(4,id)
c$$$               nnd6 = lmelv_prmry(6,id)
c$$$               nnd10 = lmelv_prmry(10,id)
c$$$               nnd9  = lmelv_prmry(9,id)
c$$$
c$$$               v_numnp(ip,nnd2) = 1
c$$$               v_numnp(ip,nnd3) = 1
c$$$               v_numnp(ip,nnd4) = 1
c$$$               v_numnp(ip,nnd6) = 1
c$$$               v_numnp(ip,nnd10) = 1
c$$$               v_numnp(ip,nnd9)  = 1
c$$$
c$$$               WRITE(301,'(9i10)') ip, iface,nnd2,nnd3,nnd4,nnd6,nnd10,nnd9,id
c$$$            ENDIF
c$$$
c$$$
c$$$         ELSE IF(n1.EQ.1.AND.n3.EQ.1.AND.n4.EQ.1)THEN
c$$$            iface = 4
c$$$            IF(numvertx.EQ.4)THEN
c$$$               nnd1 = lmelv_prmry(1,id)
c$$$               nnd3 = lmelv_prmry(3,id)
c$$$               nnd4 = lmelv_prmry(4,id)
c$$$c-plate        if(coor(3,nnd1).GT.0.d0)THEN
c$$$               v_numnp(ip,nnd1) = 1
c$$$               v_numnp(ip,nnd3) = 1
c$$$               v_numnp(ip,nnd4) = 1
c$$$               WRITE(301,'(6i10)') ip, iface,nnd4,nnd3,nnd1,id
c$$$c-plate        ENDIF
c$$$               WRITE(501,'(2i8,x,1a3,x,10i8)') numel_2d(ip),1,'tri',lmelv_prmry(1,id),lmelv_prmry(3:4,id)
c$$$            ELSE
c$$$               nnd1 = lmelv_prmry(1,id)
c$$$               nnd3 = lmelv_prmry(3,id)
c$$$               nnd4 = lmelv_prmry(4,id)
c$$$               nnd7 = lmelv_prmry(7,id)
c$$$               nnd10 = lmelv_prmry(10,id)
c$$$               nnd8  = lmelv_prmry(8,id)
c$$$
c$$$               v_numnp(ip,nnd1) = 1
c$$$               v_numnp(ip,nnd3) = 1
c$$$               v_numnp(ip,nnd4) = 1
c$$$               v_numnp(ip,nnd7) = 1
c$$$               v_numnp(ip,nnd10) = 1
c$$$               v_numnp(ip,nnd8) = 1
c$$$
c$$$               WRITE(301,'(9i10)') ip, iface,nnd4,nnd3,nnd1,nnd10,nnd7,nnd8,id
c$$$            ENDIF
c$$$         ELSE
c$$$            PRINT*,'Error in pressure face numbering'
c$$$            STOP
c$$$         ENDIF
c$$$C Non-zero load value 
c$$$         READ(10,*) press
c$$$         IF(press.GT.10.d0)THEN
c$$$            WRITE(301,*) 0 ! no regression
c$$$         ELSE
c$$$            WRITE(301,*) 1 ! regression
c$$$         ENDIF
c$$$! NOTE:: IMPORTANT the value that has the highest priority will get the value
c$$$         IF(iface.EQ.1)THEN
c$$$            IF(press_nodal(nnd1).EQ.0.d0)THEN
c$$$               press_nodal(nnd1) = press
c$$$            ELSE
c$$$               press_nodal(nnd1) = MIN(press,press_nodal(nnd1))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd2).EQ.0.d0)THEN
c$$$               press_nodal(nnd2) = press
c$$$            ELSE
c$$$               press_nodal(nnd2) = MIN(press,press_nodal(nnd2))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd3).EQ.0.d0)THEN
c$$$               press_nodal(nnd3) = press
c$$$            ELSE
c$$$               press_nodal(nnd3) = MIN(press,press_nodal(nnd3))
c$$$            ENDIF
c$$$            IF(numvertx.GT.4)THEN
c$$$            IF(press_nodal(nnd5).EQ.0.d0)THEN
c$$$               press_nodal(nnd5) = press
c$$$            ELSE
c$$$               press_nodal(nnd5) = MIN(press,press_nodal(nnd5))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd6).EQ.0.d0)THEN
c$$$               press_nodal(nnd6) = press
c$$$            ELSE
c$$$               press_nodal(nnd6) = MIN(press,press_nodal(nnd6))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd7).EQ.0.d0)THEN
c$$$               press_nodal(nnd7) = press
c$$$            ELSE
c$$$               press_nodal(nnd7) = MIN(press,press_nodal(nnd7))
c$$$            ENDIF
c$$$            ENDIF
c$$$         ELSE IF(iface.EQ.2)THEN
c$$$            IF(press_nodal(nnd1).EQ.0.d0)THEN
c$$$               press_nodal(nnd1) = press
c$$$            ELSE
c$$$               press_nodal(nnd1) = MIN(press,press_nodal(nnd1))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd2).EQ.0.d0)THEN
c$$$               press_nodal(nnd2) = press
c$$$            ELSE
c$$$               press_nodal(nnd2) = MIN(press,press_nodal(nnd2))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd4).EQ.0.d0)THEN
c$$$               press_nodal(nnd4) = press
c$$$            ELSE
c$$$               press_nodal(nnd4) = MIN(press,press_nodal(nnd4))
c$$$            ENDIF
c$$$            IF(numvertx.GT.4)THEN
c$$$            IF(press_nodal(nnd5).EQ.0.d0)THEN
c$$$               press_nodal(nnd5) = press
c$$$            ELSE
c$$$               press_nodal(nnd5) = MIN(press,press_nodal(nnd5))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd9).EQ.0.d0)THEN
c$$$               press_nodal(nnd9) = press
c$$$            ELSE
c$$$               press_nodal(nnd9) = MIN(press,press_nodal(nnd9))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd8).EQ.0.d0)THEN
c$$$               press_nodal(nnd8) = press
c$$$            ELSE
c$$$               press_nodal(nnd8) = MIN(press,press_nodal(nnd8))
c$$$            ENDIF
c$$$            ENDIF
c$$$         ELSE IF(iface.EQ.3)THEN
c$$$            IF(press_nodal(nnd3).EQ.0.d0)THEN
c$$$               press_nodal(nnd3) = press
c$$$            ELSE
c$$$               press_nodal(nnd3) = MIN(press,press_nodal(nnd3))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd2).EQ.0.d0)THEN
c$$$               press_nodal(nnd2) = press
c$$$            ELSE
c$$$               press_nodal(nnd2) = MIN(press,press_nodal(nnd2))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd4).EQ.0.d0)THEN
c$$$               press_nodal(nnd4) = press
c$$$            ELSE
c$$$               press_nodal(nnd4) = MIN(press,press_nodal(nnd4))
c$$$            ENDIF
c$$$            IF(numvertx.GT.4)THEN
c$$$            IF(press_nodal(nnd6).EQ.0.d0)THEN
c$$$               press_nodal(nnd6) = press
c$$$            ELSE
c$$$               press_nodal(nnd6) = MIN(press,press_nodal(nnd6))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd10).EQ.0.d0)THEN
c$$$               press_nodal(nnd10) = press
c$$$            ELSE
c$$$               press_nodal(nnd10) = MIN(press,press_nodal(nnd10))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd9).EQ.0.d0)THEN
c$$$               press_nodal(nnd9) = press
c$$$            ELSE
c$$$               press_nodal(nnd9) = MIN(press,press_nodal(nnd9))
c$$$            ENDIF
c$$$            ENDIF
c$$$         ELSE IF(iface.EQ.4)THEN
c$$$            IF(press_nodal(nnd3).EQ.0.d0)THEN
c$$$               press_nodal(nnd3) = press
c$$$            ELSE
c$$$               press_nodal(nnd3) = MIN(press,press_nodal(nnd3))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd1).EQ.0.d0)THEN
c$$$               press_nodal(nnd1) = press
c$$$            ELSE
c$$$               press_nodal(nnd1) = MIN(press,press_nodal(nnd1))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd4).EQ.0.d0)THEN
c$$$               press_nodal(nnd4) = press
c$$$            ELSE
c$$$               press_nodal(nnd4) = MIN(press,press_nodal(nnd4))
c$$$            ENDIF
c$$$            IF(numvertx.GT.4)THEN
c$$$            IF(press_nodal(nnd7).EQ.0.d0)THEN
c$$$               press_nodal(nnd7) = press
c$$$            ELSE
c$$$               press_nodal(nnd7) = MIN(press,press_nodal(nnd7))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd10).EQ.0.d0)THEN
c$$$               press_nodal(nnd10) = press
c$$$            ELSE
c$$$               press_nodal(nnd10) = MIN(press,press_nodal(nnd10))
c$$$            ENDIF
c$$$            IF(press_nodal(nnd8).EQ.0.d0)THEN
c$$$               press_nodal(nnd8) = press
c$$$            ELSE
c$$$               press_nodal(nnd8) = MIN(press,press_nodal(nnd8))
c$$$            ENDIF
c$$$            ENDIF
c$$$         ENDIF
c$$$
c$$$      ELSE IF(itype .EQ. 8)THEN ! node displacements
c$$$C 
c$$$C     CID = Coordinate frame ID
c$$$C     ICOMP = 6 displacement component flags (0 or 1)
c$$$C     note the flag is passed to node number one (i.e. fixed displacement)
c$$$
c$$$         READ(10,'()')
c$$$         READ(10,*) value
c$$$         numbc_prmry = numbc_prmry + 1
c$$$C
c$$$C non-blank displacement components
c$$$C note: important!!! when entering the displacement values leave all other items blank
c$$$
c$$$         ibcflg(id) = INT(value)
c$$$         
c$$$      ELSE IF(itype .EQ. 7)THEN ! mesh motion boundary conditions
c$$$C 
c$$$C     CID = Coordinate frame ID
c$$$C     ICOMP = 6 displacement component flags (0 or 1)
c$$$C     note the flag is passed to node number one (i.e. fixed displacement)
c$$$
c$$$         READ(10,'()')
c$$$         READ(10,*) value
c$$$         numbc_prmry_mm = numbc_prmry_mm + 1
c$$$C
c$$$C non-blank displacement components
c$$$C note: important!!! when entering the displacement values leave all other items blank
c$$$         ibcflg_mm(id) = INT(value)
c$$$
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
c$$$         PRINT*,'Source of possable error, STOPPING...'
c$$$         STOP
c$$$      ENDIF
c$$$      ENDDO
c$$$
c$$$      numnp_2d(:) = SUM( v_numnp(:,:), 2)
c$$$      CLOSE(301)
c$$$
c$$$      WRITE(12,*) numbc_prmry
c$$$      IF(numbc_prmry.NE.0)THEN
c$$$         DO i = 1, numnp_prmry
c$$$            IF(ibcflg(i).NE.0)THEN
c$$$               WRITE(12,'(2i10)') i,ibcflg(i)
c$$$            ENDIF
c$$$         ENDDO
c$$$      ENDIF
c$$$
c$$$! subroutine
c$$$      itype = 20
c$$$      IF(itype.EQ.20)THEN
c$$$         numel_tri = SUM(numel_2d(:))
c$$$         WRITE(12,*) numel_tri
c$$$         numnp_tri = SUM(numnp_2d)
c$$$         IF(numvertx.EQ.4)THEN
c$$$            numvertx2d = 3
c$$$         ELSE
c$$$            numvertx2d = 6
c$$$         ENDIF
c$$$         ALLOCATE(lmtri(1:numvertx2d,1:numel_tri))
c$$$         ALLOCATE(elm_2D(1:numel_tri))
c$$$         ALLOCATE(elm_2D_flag(1:numel_tri))
c$$$         ALLOCATE(epart_2d(1:numel_tri))
c$$$         elm_2D(:) = 0
c$$$         elm_2D_flag(:) = 0
c$$$         OPEN(302,FILE='fort.301')
c$$$         REWIND(302)
c$$$         DO i = 1, numel_tri
c$$$            READ(302,*) ip,iaux,lmtri(:,i),iaux1 ! iaux1 is the volumetic elment id
c$$$            elm_2D(i) = iaux1
c$$$            epart_2d(i) = ip
c$$$            READ(302,*) iaux1
c$$$            elm_2D_flag(i) = iaux1
c$$$            WRITE(12,'(5i10)') 0,lmtri(:,i), iaux1
c$$$         ENDDO
c$$$         CLOSE(302)
c$$$
c$$$         CLOSE(10)
c$$$         CLOSE(12)
c$$$
c$$$         STOP
c$$$
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
c$$$c         DO i=1,numel_tri
c$$$c            WRITE(899,*) neigh_2d(1:3,1)
c$$$c         ENDDO
c$$$
c$$$         ALLOCATE(imin(1:numnp_prmry),nmin(1:numel_tri),imax(1:numnp_prmry),nmax(1:numel_tri))
c$$$         ALLOCATE(ninc(1:numel_tri))
c$$$
c$$$         CALL STRUCT_tri(numnp_prmry,numel_tri,lmtri,neigh_2d,IMIN,NMIN,IMAX,NMAX,NINC)
c$$$
c$$$         DEALLOCATE(imin,nmin,imax,nmax,ninc)
c$$$
c$$$         DO i=1,numel_tri
c$$$            WRITE(900,*) neigh_2d(1:3,i)
c$$$         ENDDO
c$$$
c$$$c$$$      CALL rocface_data(lmtri,nprocs,numnp_prmry,numel_tri,epart_2d,
c$$$c$$$     $     neigh_2d,coor,press_nodal,prefx,prefx_lngth,elm_2D,ik1_el,numelv_prmry)
c$$$c$$$
c$$$c$$$      ALLOCATE(neigh_2d(1:3,1:numel_tri))
c$$$c$$$      
c$$$c$$$      ALLOCATE(imin(1:numnp_tri),nmin(1:numel_tri))
c$$$c$$$      ALLOCATE(imax(1:numnp_tri),nmax(1:numel_tri))
c$$$c$$$      ALLOCATE(ninc(1:numel_tri))
c$$$c$$$
c$$$c$$$      CALL STRUCT_tri(numnp_tri,numel_tri,lmtri,neigh_2d,IMIN,NMIN,IMAX,NMAX,NINC)
c$$$
c$$$
c$$$      ENDIF
c$$$
c$$$
c$$$!subroutine
c$$$c$$$C
c$$$c$$$  C --- write surface to avs of original mesh
c$$$c$$$
c$$$c$$$      WRITE(2001,'(5i8)')  numnp_2d,numel_2d, 1, 0, 0
c$$$c$$$      DO i = 1, numnp
c$$$c$$$         IF((v_numnp(i).GT.0)THEN
c$$$c$$$            WRITE(2001,'(i8,3e14.6))') i,coor(1:3,i)
c$$$c$$$         ENDIF
c$$$c$$$      ENDDO
c$$$c$$$      DO j=1,numel_2d
c$$$c$$$         WRITE(2001,'(2i8,x,1a3,x,10i8)') j,1,'tri',ik1(lmtri_2d(1:3,j),1)
c$$$c$$$      ENDDO
c$$$
c$$$      WRITE(500,'(5i8)')  SUM(v_numnp(1,:)), numel_2d, 1, 0, 0
c$$$      
c$$$      DO i = 1, numnp_prmry
c$$$         IF(v_numnp(1,i).EQ.1)THEN
c$$$            WRITE(500,'(i8,3e14.6))') i,coor(1:3,i)
c$$$         ENDIF
c$$$      ENDDO
c$$$      WRITE(501,'(6i8)') 1,1
c$$$      WRITE(501,'(a18)') 'Displacement,units' 
c$$$      DO i = 1, numnp_prmry
c$$$         IF(v_numnp(1,i).EQ.1)THEN
c$$$            WRITE(501,'(i8,3e14.6))') i,1.d0
c$$$         ENDIF
c$$$      ENDDO
c$$$      CLOSE(501)
c$$$      CLOSE(500)
c$$$
c$$$      CALL system('rm -f press.surf.inp')
c$$$      CALL system('cat fort.500 fort.501 > press.surf.inp')
c$$$C
c$$$C -mesh motion nodes
c$$$c$$$      DO j = 1 , nprocs
c$$$c$$$         DO i=1, numnp_prmry
c$$$c$$$            IF(v_numnp(j,i).EQ.1)THEN
c$$$c$$$               WRITE(505,'(2i10,)') j,i
c$$$c$$$            ENDIF
c$$$c$$$         ENDDO
c$$$c$$$      ENDDO
c$$$c$$$      CLOSE(505)
c$$$      DEALLOCATE(v_numnp)
c$$$      
c$$$      PRINT*,' Number of faces load by a distributed load =', numel_2d
c$$$      PRINT*,' Number of nodes with a boundary flag =',  numbc_prmry
c$$$
c$$$      CLOSE(10)
c$$$C
c$$$C --  READ THE ELEMENT NEIGHBOR ARRAY
c$$$
c$$$      OPEN(10,FILE=prefx(1:prefx_lngth)//'.neigh', FORM='formatted')
c$$$
c$$$      ALLOCATE(neigh(1:4,1:numelv_prmry))
c$$$
c$$$      neigh(1:4,1:numelv_prmry) = -1
c$$$      DO i = 1, numelv_prmry
c$$$         READ(10,*) neigh(1:4,i)
c$$$      ENDDO
c$$$      CLOSE(10)

      RETURN

      END
    
