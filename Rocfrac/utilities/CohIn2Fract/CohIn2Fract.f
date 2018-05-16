      PROGRAM cohin2fract

      USE meshdata

! f90 -O3 -o CohIn2Fract Cohin2fract.f
      INTEGER :: numelc,numel

      REAL*8 :: dt_courant
      REAL*8, DIMENSION(1:3) :: coor3
      INTEGER :: iaux, iaux1, iaux2
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lmvol

      INTEGER, DIMENSION(1:10) :: lm

      REAL*8 :: xminedge

      CALL read_deck()

      OPEN(10,file=prefx(1:prefx_lngth)//'.msh.coh')
      OPEN(12,FILE=
     $     prefx(1:prefx_lngth)//'/'//prefx(1:prefx_lngth)//'.0000.inp',
     $     STATUS='replace',FORM='formatted')     
c$$$  WRITE(12,*) 1, 1 ! fix? still us? SUM(numnp(1:nprocs)),SUM(numel(1:nprocs))
c
c -- Output the nodal cordinates

      READ(10,*) numelv, numnp, numelc, xminedge

      WRITE(12,*) xminedge ! minimum edge length

      numel = numelv + numelc

      WRITE(12,*) numnp
      DO j = 1, numnp
         READ(10,*) iaux,coor3(:)
         WRITE(12,'(i9,3(1x,e14.5))') iaux,coor3(:)
      ENDDO
      print*,' Write Nodal Boundary Conditions'
      READ(10,*) numbc
      WRITE(12,*) numbc
      DO j = 1, numbc
         READ(10,*) iaux1, iaux2
         WRITE(12,*) iaux1, iaux2
      ENDDO
      WRITE(12,*) 0 ! mesh motion bcs.
      
      ALLOCATE(lmvol(1:5,1:numelv))
      jj = 0
       WRITE(12,'(12i10)') numelc,0,0,0,0
      DO j = 1, numel
!              id  | id flag | prop id | connectivity
         READ(10,*) iaux1, iaux2, iaux3, lm(1:iaux2+4)
         IF(iaux2.EQ.2)THEN
            WRITE(12,'(12i10)') iaux3, lm(1:iaux2+4),-1
         ELSE
            jj = jj + 1
            lmvol(1,jj) = iaux3
            lmvol(2:5,jj) = lm(1:iaux2+4)
         ENDIF

      ENDDO

      WRITE(12,*) 0
      WRITE(12,'(12i10)') jj,0,0,4,2 ! not sure what this last 2 are
      DO j = 1, jj
         WRITE(12,'(12i10)') lmvol(1:5,j),0
      ENDDO 

      DO i = 1, 6
         WRITE(12,*) 0
      ENDDO
      WRITE(12,*) 0
! last one was used for pov-ray i think, will not work
! number of boundary triangles
      
      
      CLOSE(4000)

      END
