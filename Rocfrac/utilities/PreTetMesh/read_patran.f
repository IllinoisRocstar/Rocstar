      SUBROUTINE read_patran()

      USE meshdata

      IMPLICIT NONE

      INTEGER :: i,j,id,numvertx,itype

      REAL*8 :: value
      REAL*8 :: press

      PRINT*,'MESH OPTION:'
      PRINT*,'  READING PATRAN MESH'
      PRINT*,' '

      OPEN(io_input,file=prefx(1:prefx_lngth)//'.pat', FORM='formatted')
!
! - Packet Type 25: Title Card, Format(i2,8i8)
      READ(io_input,'()')
      READ(io_input,'()')

! -  Packet Type 26: Summary Data, Format(i2,8i8)
! -      26 ID IV KC N1 N2 N3 N4 N5
! -  N1 = number of nodes
! -  N2 = number of elements
! -  N3 = number of materials
! -  N4 = number of Element Properties
! -  N5 = number of Coordinate Frames

      READ(io_input,*) i,i,i,i,numnp,numelv
! - part of title card
      READ(io_input,'()')

      PRINT*,' Number of 2D skin elements = ', numelv
      PRINT*,' Number of skin nodes = ', numnp
C
C --  Read Nodal coordinates

      ALLOCATE(coor(1:3,1:numnp),lmelv(1:3,1:numelv))
      ALLOCATE(npress(1:numelv))
      npress(:) = 0

      DO i = 1, numnp

C - Type 01: Node Data,  format(i2,8i8)
C        1 ID IV KC
C    ID = node id
C    IV = 0 n/a
C    KC = 2
         READ(io_input,*) j,id

C - format(3e16.9)
C - caresian coordinate of Nodes x, y, z
         READ(io_input,'(3e16.9)') coor(1:3,id)
C - ICF GTYPE NDF CONFIG CID PSPC
         READ(io_input,'()')
      ENDDO
C
C --  Read element connectivity array
 
      DO i = 1, numelv

C - Type 02: Element Data, format(i2,8i8)
C        2 ID IV KC N1
C     ID = element ID
C     IV = Shape (5 = tet)
         READ(10,*) j, id
C     NODES CONFIG CID CEID 
         READ(10,'(i8)') numvertx
C     LNODES
C        LNODES = Element corner nodes followed by additional nodes
         READ(10,*) lmelv(1:numvertx,id)
      ENDDO
C
C
C -- Continue to read in the rest of the input parameters

      ALLOCATE(ibcflg(1:numnp))
      ibcflg(:) = 0

      DO
         READ(10,*) itype,id

         IF(itype.EQ.99) EXIT
         IF(itype .EQ. 4)THEN   ! element properties
            READ(10,'()')
         ELSE IF(itype .EQ. 6)THEN ! pressure loading
            READ(10,'()')
C Non-zero load value 
            READ(10,*) press
            npress(id) = INT(press)
         ELSE IF(itype .EQ. 8)THEN ! node displacements

! USED FOR THE PHYSICAL BOUNDARY CONDITION FLAGS AT NODES

C 
C     CID = Coordinate frame ID
C     ICOMP = 6 displacement component flags (0 or 1)
C     note the flag is passed to node number one (i.e. fixed displacement)

            READ(10,'()')
C
C non-blank displacement components
C note: important!!! when entering the displacement values leave all other items blank
            READ(10,*) value
            ibcflg(id) = ibcflg(id) + INT(value)*100
         ELSE IF(itype .EQ. 7)THEN ! mesh motion boundary conditions

! USED FOR THE MESHMOTION BOUNDARY CONDITION FLAGS AT NODES

C 
C     CID = Coordinate frame ID
C     ICOMP = 6 displacement component flags (0 or 1)
C     note the flag is passed to node number one (i.e. fixed displacement)

            READ(10,'()')
            READ(10,'(1e16.9)') value
C
C non-blank displacement components
C note: important!!! when entering the displacement values leave all other items blank
            ibcflg(id) = ibcflg(id) + INT(value)
         ELSE IF(itype .EQ. 5)THEN ! coordinate frames
            READ(10,'()')
            READ(10,'()')
            READ(10,'()')
            READ(10,'()')
         ELSE IF(itype .EQ. 10)THEN ! nodal temperature (used to flag axisymetric)
            READ(10,*) value
         ELSE
            PRINT*,'packet type',itype,' not accounted for in program'
            PRINT*,'Source of possable error, STOPPING...'
            STOP
         ENDIF
      ENDDO

      CLOSE(10)

      END
    
