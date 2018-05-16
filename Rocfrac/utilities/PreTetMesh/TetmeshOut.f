      SUBROUTINE TetMeshOut(NBELEM,NUMSOM,REFNUM,NSD,NBSINI,COOR,NRS,
     $     prefx,prefx_lngth)

      IMPLICIT NONE

      INTEGER :: i

      CHARACTER*20 :: prefx     ! I/O file prefx
      INTEGER :: prefx_lngth    ! length of prefx

! The number of elements: faces of the skin mesh, plus the enforced internal edges and faces, if any.      
      INTEGER :: NBELEM
! A reserved parameter (not sure what this is)
      INTEGER :: NRESER = 0 
! The number of vertices of the element (3 = triangle)
      INTEGER :: NTYPE = 3
! The vetex numbers of the face
      INTEGER, DIMENSION(1:3,1:NBELEM) :: NUMSOM
! The face attribute, characterizing the physical property of the face
      INTEGER,DIMENSION(1:NBELEM) :: NSD
! The attributes of the face and edges
! edge 1 is the edge joining vertex 1 to vertex 2, 
! edge 2 is the edge joining vertex 2 to vertex 3,
! edge 3 is the edge joining vertex 1 to vertex 3,
      INTEGER, DIMENSION(1:3,1:NBELEM) :: REFNUM
! The number of intial vertices of the skin mesh
      INTEGER :: NBSINI
! The Vertex Coordinates
      REAL*8, DIMENSION(1:3,1:NBSINI) :: COOR
! The local absolute stepsize of a specified point (local size of the tetrahedra holding this point)
!     REAL*8, DIMENSION(1:NBSINI) :: SIZE
! The vertex attribure, characteizing the physical property associated to the vertex
      INTEGER, DIMENSION(1:NBSINI) :: NRS


      OPEN(10,FILE=prefx(1:prefx_lngth)//'.faces',STATUS='unknown')

      WRITE(10,'(2i10)') NBELEM, NRESER

      DO i = 1, NBELEM
         WRITE(10,'(10i10)') NTYPE,NUMSOM(1:NTYPE,i),NSD(i),REFNUM(1:NTYPE,i)
      ENDDO

      CLOSE(10)

      OPEN(12,FILE=prefx(1:prefx_lngth)//'.points',STATUS='unknown')
      
      WRITE(12,'(1i10)') NBSINI

      DO i = 1, NBSINI
         WRITE(12,'(3e16.9,1i10)') COOR(1:3,i),NRS(i)
      ENDDO

      CLOSE(12)

      RETURN 
      END

      

      
