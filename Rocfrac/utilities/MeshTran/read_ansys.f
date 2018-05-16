      SUBROUTINE read_ansys(numbc_prmry, dhmin)

      USE meshdata

      IMPLICIT NONE

      INTEGER :: numbc_prmry
      REAL*8 :: dhmin
! Local
      REAL*8 :: xx,yy,zz,size1,size2,size3,size4,size5,size6
      INTEGER :: ii, i, iaux,iaux1, ios
      INTEGER :: id
      CHARACTER*9 :: ichr
      REAL*8 :: value

      PRINT*,'MESH OPTION:'
      PRINT*,'  READING ANSYS MESH'
      PRINT*,' '
!
! Summary input file of the mesh
!   NOTE: The file is not automatically created with ansys, user must
!         edit the file and provide the information.

      OPEN(10,FILE=prefx(1:prefx_lngth)//'.ansys', FORM='formatted')
      
      READ(10,*) numelv_prmry
      READ(10,*) numnp_prmry
      READ(10,*) numbc_prmry
      READ(10,*) numel_2d
      READ(10,*) numnp_2d

      PRINT*,' Number of elements = ', numelv_prmry
      PRINT*,' Number of nodes = ', numnp_prmry

      CLOSE(10)
C
C --  Read Nodal coordinates

      OPEN(10,file=prefx(1:prefx_lngth)//'.nod', FORM='formatted')

      ALLOCATE(coor(1:3,1:numnp_prmry))

      DO i=1,5
         READ(10,'()') 
      ENDDO

      ii = 0
      DO i = 1, numnp_prmry
         READ(10,*) iaux,coor(1,i), coor(2,i), coor(3,i)
         ii = ii + 1
         IF(ii.EQ.20 .AND. i.NE.numnp_prmry)THEN
            READ(10,'(a1)')
            READ(10,'(a1)')
            ii = 0
         ENDIF
      ENDDO
      CLOSE(10)
C
C --  Read element connectivity array

      OPEN(10,file=prefx(1:prefx_lngth)//'.ele', FORM='formatted')

      ALLOCATE(lmelv_prmry(1:numvertx,1:numelv_prmry))

      DO i=1,5
         READ(10,'()')
      ENDDO
      ii = 0
      IF(numvertx.EQ.4)THEN
         DO i = 1, numelv_prmry
            READ(10,*) iaux1,iaux,iaux,iaux,iaux,iaux,lmelv_prmry(1:4,i)
            ii = ii + 1
            IF(ii.EQ.20.AND.i.NE.numelv_prmry)THEN
               READ(10,'()')
               READ(10,'()')
               READ(10,'()')
               ii = 0
            ENDIF
c
c -- Find the size of the smallest element
c
            xx = coor(1,lmelv_prmry(1,i)) - coor(1,lmelv_prmry(2,i))
            yy = coor(2,lmelv_prmry(1,i)) - coor(2,lmelv_prmry(2,i))
            zz = coor(3,lmelv_prmry(1,i)) - coor(3,lmelv_prmry(2,i))
            size1 = SQRT(xx*xx+yy*yy+zz*zz)
            xx = coor(1,lmelv_prmry(2,i)) - coor(1,lmelv_prmry(3,i))
            yy = coor(2,lmelv_prmry(2,i)) - coor(2,lmelv_prmry(3,i))
            zz = coor(3,lmelv_prmry(2,i)) - coor(3,lmelv_prmry(3,i))
            size2 = SQRT(xx*xx+yy*yy+zz*zz)
            xx = coor(1,lmelv_prmry(3,i)) - coor(1,lmelv_prmry(1,i))
            yy = coor(2,lmelv_prmry(3,i)) - coor(2,lmelv_prmry(1,i))
            zz = coor(3,lmelv_prmry(3,i)) - coor(3,lmelv_prmry(1,i))
            size3 = SQRT(xx*xx+yy*yy+zz*zz)
            xx = coor(1,lmelv_prmry(4,i)) - coor(1,lmelv_prmry(1,i))
            yy = coor(2,lmelv_prmry(4,i)) - coor(2,lmelv_prmry(1,i))
            zz = coor(3,lmelv_prmry(4,i)) - coor(3,lmelv_prmry(1,i))
            size4 = SQRT(xx*xx+yy*yy+zz*zz)
            xx = coor(1,lmelv_prmry(4,i)) - coor(1,lmelv_prmry(2,i))
            yy = coor(2,lmelv_prmry(4,i)) - coor(2,lmelv_prmry(2,i))
            zz = coor(3,lmelv_prmry(4,i)) - coor(3,lmelv_prmry(2,i))
            size5 = SQRT(xx*xx+yy*yy+zz*zz)
            xx = coor(1,lmelv_prmry(4,i)) - coor(1,lmelv_prmry(3,i))
            yy = coor(2,lmelv_prmry(4,i)) - coor(2,lmelv_prmry(3,i))
            zz = coor(3,lmelv_prmry(4,i)) - coor(3,lmelv_prmry(3,i))
            size6 = SQRT(xx*xx+yy*yy+zz*zz)
            dhmin = MIN(size1,size2,size3,size4,size5,size6,dhmin)
         ENDDO
      ELSE IF(numvertx.EQ.10)THEN
         DO i = 1, numelv_prmry
            READ(10,*) iaux1,iaux,iaux,iaux,iaux,iaux,lmelv_prmry(1:10,i)
            ii = ii + 1
            IF(ii.EQ.20.AND.i.NE.numelv_prmry)THEN
               READ(10,'()')
               READ(10,'()')
               READ(10,'()')
               ii = 0
            ENDIF
c
c -- Find the size of the smallest element
c
            xx = coor(1,lmelv_prmry(1,i)) - coor(1,lmelv_prmry(2,i))
            yy = coor(2,lmelv_prmry(1,i)) - coor(2,lmelv_prmry(2,i))
            zz = coor(3,lmelv_prmry(1,i)) - coor(3,lmelv_prmry(2,i))
            size1 = SQRT(xx*xx+yy*yy+zz*zz)
            xx = coor(1,lmelv_prmry(2,i)) - coor(1,lmelv_prmry(3,i))
            yy = coor(2,lmelv_prmry(2,i)) - coor(2,lmelv_prmry(3,i))
            zz = coor(3,lmelv_prmry(2,i)) - coor(3,lmelv_prmry(3,i))
            size2 = SQRT(xx*xx+yy*yy+zz*zz)
            xx = coor(1,lmelv_prmry(3,i)) - coor(1,lmelv_prmry(1,i))
            yy = coor(2,lmelv_prmry(3,i)) - coor(2,lmelv_prmry(1,i))
            zz = coor(3,lmelv_prmry(3,i)) - coor(3,lmelv_prmry(1,i))
            size3 = SQRT(xx*xx+yy*yy+zz*zz)
            xx = coor(1,lmelv_prmry(4,i)) - coor(1,lmelv_prmry(1,i))
            yy = coor(2,lmelv_prmry(4,i)) - coor(2,lmelv_prmry(1,i))
            zz = coor(3,lmelv_prmry(4,i)) - coor(3,lmelv_prmry(1,i))
            size4 = SQRT(xx*xx+yy*yy+zz*zz)
            xx = coor(1,lmelv_prmry(4,i)) - coor(1,lmelv_prmry(2,i))
            yy = coor(2,lmelv_prmry(4,i)) - coor(2,lmelv_prmry(2,i))
            zz = coor(3,lmelv_prmry(4,i)) - coor(3,lmelv_prmry(2,i))
            size5 = SQRT(xx*xx+yy*yy+zz*zz)
            xx = coor(1,lmelv_prmry(4,i)) - coor(1,lmelv_prmry(3,i))
            yy = coor(2,lmelv_prmry(4,i)) - coor(2,lmelv_prmry(3,i))
            zz = coor(3,lmelv_prmry(4,i)) - coor(3,lmelv_prmry(3,i))
            size6 = SQRT(xx*xx+yy*yy+zz*zz)
            dhmin = MIN(size1,size2,size3,size4,size5,size6,dhmin)
         ENDDO
      ENDIF
      CLOSE(10)
C
C -- Read the boundary condition file

      OPEN(10,FILE=prefx(1:prefx_lngth)//'.bc',  FORM='formatted',IOSTAT=ios)
      IF(ios .NE. 0) THEN
         PRINT*,' .... No boundary condition file'
      ELSE
         ALLOCATE(ibcflg(1:numnp_prmry))
         ibcflg(:) = 0
         DO i=1,5 
            READ(10,'()') 
         ENDDO
         ii = 0
         DO i = 1, numbc_prmry
            READ(10,'(i9,a8,e14.5)') id,ichr,value ! node, ,value
            ii = ii + 1
            IF(ii.EQ.20.AND.numbc_prmry.NE.i)THEN
               READ(10,'()')
               READ(10,'()')
               ii = 0
            ENDIF
            ibcflg(id) = INT(value)
         ENDDO
      ENDIF
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
