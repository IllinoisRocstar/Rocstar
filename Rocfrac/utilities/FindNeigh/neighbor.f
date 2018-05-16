      PROGRAM neighbor
C
C     guidef90 -O3 -o neigh neighbor.f
C     setenv OMP_NUM_THREADS 1
C     
      IMPLICIT REAL*8 (a-h,o-z)

      CHARACTER*8 ichr
      CHARACTER*1 ichr1
      CHARACTER*12 ichr12

      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: neigh
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lm
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lm_1d
      INTEGER, ALLOCATABLE, DIMENSION(:) :: imin,nmin,imax,nmax,ninc
      INTEGER :: IOERROR

      CHARACTER*200 :: keywd    ! keyword parameter

      INTEGER :: ios            ! io error
      INTEGER :: io_input
      CHARACTER*20 :: prefx
      INTEGER :: prefx_lngth

      INTEGER :: TetForm

      INTEGER :: norder
      INTEGER :: icheck

      INTEGER :: NBELE,NBELEF,LOELE,LOELEF,IDEB,IFIN

      INTEGER :: IOFORMAT

      TetForm = 0
      IOformat = 1 ! acsii by default

C
C -- Open Analysis Deck File

      io_input = 10

      OPEN(io_input,FILE='fractography3d.inp',STATUS='old',IOSTAT=ios)
      IF(ios.NE.0)THEN
         PRINT*, 'Unable to find fractography3d.inp'
         PRINT*, ' ...Trying RocfracControl.txt'

         OPEN(io_input,FILE='RocfracControl.txt',STATUS='old',IOSTAT=ios)
         IF(ios.NE.0)THEN
            PRINT*, 'Unable to find RocfracControl.txt'
            PRINT*, ' ...STOPPING'
            STOP
         ENDIF
      ENDIF
      
C
C -- Read Analysis Deck File
      icheck = 0

      REWIND io_input
 10   READ(io_input,'(A)',IOSTAT=ios) keywd
!      print*,keywd
      IF(ios.LT.0) GOTO 20
      IF(keywd(1:7).EQ.'*PREFIX') THEN
         icheck = icheck + 1
         CALL PREFIX_SUB(io_input,prefx,prefx_lngth)
         PRINT*,'  PREFIX = ', prefx(1:prefx_lngth)
         IF(icheck.EQ.2)THEN
            GOTO 30
         ELSE
            GOTO 10
         ENDIF 
      ELSE IF(keywd(1:9).EQ.'*IOFORMAT') THEN
         icheck = icheck + 1
         CALL IOFORMAT_SUB(IOformat,io_input)
         IF(icheck.EQ.2)THEN
            GOTO 30
         ELSE
            GOTO 10
         ENDIF
      ELSE IF(keywd(1:9).EQ.'*MESHSOFT') THEN
         icheck = icheck + 1
         CALL MESHSOFT_SUB(io_input,iansys,ipatran,itetmesh) 
         IF(icheck.EQ.2)THEN
            GOTO 30
         ELSE
            GOTO 10
         ENDIF
      ELSE
         GOTO 10
      ENDIF

 20   PRINT*,'ERROR: *PREFIX AND/OR MESHSOFT not found in fractography3d.inp'
      PRINT*,'... STOPPING'
      STOP

 30   CONTINUE
      CLOSE(io_input)

      io_input = 11
C
C -- Open the element file
      IF(iansys.EQ.1)THEN
         OPEN(io_input,FILE=prefx(1:prefx_lngth)//'.ele',
     $        FORM='formatted',STATUS='old')
         PRINT*,'ENTER NUMBER OF TETRAHEDRAL ELEMENTS'
         READ*,numel
         ALLOCATE(lm(1:4,1:numel),neigh(1:4,1:numel))
         icount=0
         DO i=1,5
            READ(io_input,'(a1)')
         ENDDO
         DO i = 1,numel
            READ(io_input,*) iaux1,iaux,iaux,iaux,iaux,iaux,
     $           lm(1,i),lm(2,i),lm(3,i),lm(4,i)
            icount = icount+1
            IF(icount.EQ.20.AND.i.NE.numel)THEN
               READ(io_input,'(a1)') ichr1
               READ(io_input,'(a12)') ichr12
               READ(io_input,'(a1)') ichr1
               icount = 0
            ENDIF
         ENDDO

      ELSE IF(ipatran.EQ.1)THEN
         OPEN(io_input,FILE=prefx(1:prefx_lngth)//'.pat',
     $        FORM='formatted',STATUS='old')

C
C - Packet Type 25: Title Card, Format(i2,8i8)
         READ(io_input,'()')
         READ(io_input,'()')

C -  Packet Type 26: Summary Data, Format(i2,8i8)
C -      26 ID IV KC N1 N2 N3 N4 N5
C -  N1 = number of nodes
C -  N2 = number of elements
C -  N3 = number of materials
C -  N4 = number of Element Properties
C -  N5 = number of Coordinate Frames

      READ(io_input,*) i,i,i,i,numnp, numel
C - part of title card
      READ(io_input,'()')

      PRINT*,' Number of elements = ', numel
      PRINT*,' Number of nodes = ', numnp
C
C --  Read Nodal coordinates

      ALLOCATE(neigh(1:4,1:numel))
      
      DO i = 1, numnp

C - Type 01: Node Data,  format(i2,8i8)
C        1 ID IV KC
C    ID = node id
C    IV = 0 n/a
C    KC = 2
         READ(io_input,'()')
C - caresian coordinate of Nodes x, y, z
         READ(io_input,'()')
C - ICF GTYPE NDF CONFIG CID PSPC
         READ(io_input,'()')
      ENDDO
C
C --  Read element connectivity array

C - Type 02: Element Data, format(i2,8i8)
C        2 ID IV KC N1
C     ID = element ID
C     IV = Shape (5 = tet)
      READ(io_input,*) j, id
C     NODES CONFIG CID CEID 
      READ(io_input,'(i8)') norder
      ALLOCATE(lm(1:norder,1:numel))
C     LNODES
C        LNODES = Element corner nodes followed by additional nodes
      READ(io_input,*) lm(1:norder,id)

      DO i = 2, numel

C - Type 02: Element Data, format(i2,8i8)
C        2 ID IV KC N1
C     ID = element ID
C     IV = Shape (5 = tet)
         READ(io_input,*) j, id
C     NODES CONFIG CID CEID 
         READ(io_input,'(i8)') norder
C     LNODES
C        LNODES = Element corner nodes followed by additional nodes
         READ(io_input,*) lm(1:norder,id)
      ENDDO
C------------------------------------------------------------------------
C-- TETMESH INPUT
C------------------------------------------------------------------------

      ELSE IF(itetmesh.EQ.1) THEN
         PRINT*,'Tetmesh file'
! check to see if it a binary file
         OPEN(io_input,FILE=prefx(1:prefx_lngth)//'.noboiteb',
     $           FORM='unformatted',STATUS='old',IOSTAT=ios )
         TetForm = 0
! Binary file does not exist, check to see that ascii file exist
         IF(ios.NE.0)THEN
            OPEN(io_input,FILE=prefx(1:prefx_lngth)//'.noboite',
     $           FORM='formatted',STATUS='old',IOSTAT=ios )
            TetForm = 1
         ENDIF
! Neither file exist, so there is a problem
         IF(ios.NE.0)THEN
            PRINT*,'ERROR: Problem opening ',prefx(1:prefx_lngth)//'.noboite'
            PRINT*,' .. STOPPING'
            STOP
         ENDIF

         norder = 4
         IF(TetForm.EQ.0)THEN
            PRINT*, '.. Binary format'
! Record 1:
            READ(io_input) numel,nnode,iaux,iaux,iaux,
     $           NBELE,LOELE,NBELEF,LOELEF, ! for element records
     $           iaux,iaux,iaux,iaux, ! for coordinate records
     $           iaux,iaux,iaux,iaux ! for sub
         ELSE
            PRINT*,'.. ASCII format'
! Record 1:
            READ(io_input,*) numel,nnode,iaux,iaux,iaux,
     $           NBELE,LOELE,NBELEF,LOELEF, ! for element records
     $           iaux,iaux,iaux,iaux, ! for coordinate records
     $           iaux,iaux,iaux,iaux  ! for sub
         ENDIF
            
! Record 2:

         PRINT*,' Number of elements = ', numel
         PRINT*,' Number of nodes = ', nnode

         ALLOCATE(lm(1:4,1:numel))
C
C --  Read element connectivity array

         ALLOCATE(lm_1d(1:4*numel))
         IF(TetForm.EQ.0)THEN
            IDEB = 1
            IFIN = LOELE
            DO j = 1, NBELE
               READ(io_input) (lm_1d(i),i = IDEB, IFIN)
               IDEB = IDEB + LOELE
               IFIN = IFIN + LOELE
            ENDDO
            IF(NBELEF.NE.0)THEN
               IFIN = IDEB + LOELEF - 1
               READ(io_input) (lm_1d(i),i = IDEB, IFIN)
            ENDIF
         ELSE
            READ(io_input,*) lm_1d(1:numel*4)
         ENDIF

         lm = RESHAPE(lm_1d,(/4,numel/))

         DEALLOCATE(lm_1d)

         PRINT*,'   Finished element connectivity'
         ALLOCATE(neigh(1:4,1:numel))
      ENDIF

      CLOSE(io_input)

      PRINT*,'  Finished reading mesh data ...'

C
C     FIND MAXIMUM VERTEX ADDRESS
C
      nnode       = lm(1,1)
      IF(norder.EQ.4)THEN
         DO l = 1, numel 
            nnode    = MAX(nnode,lm(1,l),lm(2,l),lm(3,l),lm(4,l))
         ENDDO
      ELSE
         DO l = 1, numel 
            nnode    = MAX(nnode,lm(1,l),lm(2,l),lm(3,l),lm(4,l),
     $           lm(5,l),lm(6,l),lm(7,l),lm(8,l),lm(9,l),lm(10,l))
         ENDDO
      ENDIF

      PRINT*,'     Maximum vertex address ...', nnode

      ALLOCATE(imin(1:nnode),nmin(1:numel),imax(1:nnode),nmax(1:numel))
      ALLOCATE(ninc(1:numel))

      CALL STRUCT (nnode,numel,lm(1:4,1:numel),neigh,imin,nmin,imax,nmax,ninc)
      DEALLOCATE(imin,nmin,imax,nmax,ninc)

      io_input = 12

      IF(IOformat.EQ.0)THEN
         OPEN(io_input,FILE=prefx(1:prefx_lngth)//'.neigh',
     $        FORM='unformatted',STATUS='unknown')
         DO i=1,numel
            WRITE(io_input) neigh(2,i),neigh(3,i),neigh(4,i),neigh(1,i)
         ENDDO
      ELSE
         OPEN(io_input,FILE=prefx(1:prefx_lngth)//'.neigh',
     $        FORM='formatted',STATUS='unknown')
         
         DO i=1,numel
            WRITE(io_input,'(4I10)') neigh(2,i),neigh(3,i),neigh(4,i),neigh(1,i)
         ENDDO
      ENDIF

      CLOSE(io_input)

      END

      SUBROUTINE PREFIX_SUB(io_input,prefx,prefx_lngth)

      IMPLICIT NONE

      INTEGER :: io_input, prefx_lngth
      CHARACTER*20 :: prefx
         
      READ(io_input,'(a20)') prefx
      prefx_lngth = LEN_TRIM(prefx)
      
      RETURN
      END

      SUBROUTINE MESHSOFT_SUB(io_input,iansys,ipatran,itetmesh)

      IMPLICIT NONE

!   - meshing software flag
      INTEGER :: io_input
      INTEGER :: iansys         ! 0- no, 1-yes
      INTEGER :: ipatran        ! 0- no, 1-yes
      INTEGER :: itetmesh       ! 0- no, 1-yes

      CHARACTER*1 :: chr

      READ(io_input,*) chr

      iansys   = 0 ! set default -no-
      ipatran  = 0 
      itetmesh = 0

      IF(chr.EQ.'T'.OR.chr.EQ.'t')THEN
         itetmesh = 1
      ELSE IF(chr.EQ.'A'.OR.chr.EQ.'a')THEN
         iansys = 1
      ELSE IF(chr.EQ.'P'.OR.chr.EQ.'p')THEN
         ipatran = 1
      ELSE
         PRINT*,' ERROR: MESHING PACKAGE NOT SUPPORTED'
         PRINT*,'STOPPING'
         STOP
      ENDIF

      RETURN
      END

      SUBROUTINE IOFORMAT_SUB(IOformat,io_input)

      IMPLICIT NONE

      INTEGER :: IOformat,io_input

! 0 = binary
! 1 = ascii

      READ(io_input,*) IOformat

      RETURN
      END


     
