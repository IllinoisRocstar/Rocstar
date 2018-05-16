
!!****
!!
!!  NAME
!!     LocThermCap_v3d8
!!
!!  FUNCTION
!!     This subroutine constructs the process local capacitance matrix
!!     for volumetric stess 8 node, 3d elements as a consistent matrix.
!!     The resulting matrix is put into the variables in comp_row_global.
!!
!!  INPUTS
!!     NumEl -- The number of elements that this processor knows about
!!     NumNp -- The number of nodes that this processor knows about
!!     NumMat -- The global number of materials in the structure
!!     coor -- The coordinates of each of the nodes
!!     nodes -- The connectivity table
!!     MatType -- Mapping of material types to elements
!!     ri -- Gauss point positions within an element
!!     rho -- Density of the materials
!!     cap -- Specific heat capacity of the materials
!!     ElConnVol -- The connectivity table (redundant)
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE LocThermCap_v3d8( NumEl, NumNP, NumMat, &
     coor, nodes, MatType, ri, rho, cap, ElConnVol )

  USE Precision
  USE implicit_global
  USE comp_row_global
  USE ROCSTAR_RocFracComm

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  ! ... Input variables
  INTEGER :: NumEl, NumNP, NumMat
  INTEGER, DIMENSION(1:NumEl) :: MatType
  INTEGER, DIMENSION(1:8,1:NumEl) :: nodes
  REAL(KIND=wp), DIMENSION(1:3,1:NumNP) :: coor
  REAL(KIND=wp), DIMENSION(1:3,1:8) :: ri
  REAL(KIND=wp), DIMENSION(1:NumMat) :: rho, cap
  INTEGER, DIMENSION(1:8,1:NumEl) :: ElConnVol

  ! ... Internal variables
  REAL(kind=wp), DIMENSION(1:8) :: N
  REAL(kind=wp), DIMENSION(1:8,1:3) :: dN
  REAL(kind=wp), DIMENSION(1:8,1:8) :: NtN
  INTEGER :: igpt
  INTEGER :: ngpts = 8
  INTEGER :: ielem
  INTEGER :: imat
  INTEGER :: i, j, k, l, tempi, tempj
  INTEGER :: dof1, dof2, gdof1, gdof2, ldof1, ldof2
  REAL(kind=wp) :: tempval, tempval2, tempval3
  INTEGER :: nnzold
  REAL(kind=wp) :: element_volume
  REAL(kind=wp), DIMENSION(1:3,1:3) :: jac, jacinv
  REAL(kind=wp) :: detj
  ! ... diag_width is the max number of non-zero elements expected in a row 
  ! ... in the global capacitance matrix. There might be some algorithm to determine that a
  ! ... priori. But right now its set to 30
  INTEGER :: diag_width = 30
  REAL(kind=wp) :: t1,t2,t3,t4
  LOGICAL :: error
  INTEGER :: n1,n2,n3,n4,n5,n6,n7,n8
  REAL(kind=wp),DIMENSION(1:3,1:8) :: coordtmp
  INTEGER,DIMENSION(1:8) :: L2G_temp

  ! ... Output variables
  INTEGER, ALLOCATABLE, DIMENSION(:) :: rp, cval
  REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: aval
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:) :: GlMassMat

  ! ... 'values' and 'col_index' are the vectors 'aval' and 'cval' 
  ! ... respectively. Their size is also allocated using diag_width
  REAL(kind=wp), ALLOCATABLE, DIMENSION(:) :: values
  INTEGER, ALLOCATABLE, DIMENSION(:) :: col_index
  INTEGER :: nnz


  ! ... Initialize the output variables

  !  ALLOCATE(rp(1:3*GNumNp+1))
  !  ALLOCATE(cval(1:1))
  !  ALLOCATE(aval(1:1))
  !  ALLOCATE(GlMassMat(1:3*GNumNp,1:3*GNumNp))

  ! ... allocate enough memory of the values, col_index (=aval,cval)
  ! ... vectors to hold diag_width values for each row of the
  ! ... explicit mass matrix.
  ALLOCATE(values(1:GNumNp*diag_width))
  ALLOCATE(col_index(1:GNumNp*diag_width))



  nnz = 1
  ! rp(1) = 0
  !rp(2:3*GNumNp+1) = 1
  ! cval(1) = 0
  ! aval(1) = 0.0
  !GlMassMat(:,:)=0.0
  col_index(:)=-1
  values(:)=0.0

  !  PRINT*,'Entering Do loop', myid
  CALL MPI_BARRIER(ROCSTAR_COMMUNICATOR,ierr)
  CALL CPU_time(t1)
  ! Loop through the elements
  DO ielem = 1, NumEl

     ! ... Node positions
     n1 = ElConnVol(1,ielem)
     n2 = ElConnVol(2,ielem)
     n3 = ElConnVol(3,ielem)
     n4 = ElConnVol(4,ielem)
     n5 = ElConnVol(5,ielem)
     n6 = ElConnVol(6,ielem)
     n7 = ElConnVol(7,ielem)
     n8 = ElConnVol(8,ielem)
     ! ... x,y,z coordtmpinates of node n#.
     coordtmp(1,1) = coor(1,n1)
     coordtmp(2,1) = coor(2,n1)
     coordtmp(3,1) = coor(3,n1)
     coordtmp(1,2) = coor(1,n2)
     coordtmp(2,2) = coor(2,n2)
     coordtmp(3,2) = coor(3,n2)
     coordtmp(1,3) = coor(1,n3)
     coordtmp(2,3) = coor(2,n3)
     coordtmp(3,3) = coor(3,n3)
     coordtmp(1,4) = coor(1,n4)
     coordtmp(2,4) = coor(2,n4)
     coordtmp(3,4) = coor(3,n4)
     coordtmp(1,5) = coor(1,n5)
     coordtmp(2,5) = coor(2,n5)
     coordtmp(3,5) = coor(3,n5)
     coordtmp(1,6) = coor(1,n6)
     coordtmp(2,6) = coor(2,n6)
     coordtmp(3,6) = coor(3,n6)
     coordtmp(1,7) = coor(1,n7)
     coordtmp(2,7) = coor(2,n7)
     coordtmp(3,7) = coor(3,n7)
     coordtmp(1,8) = coor(1,n8)
     coordtmp(2,8) = coor(2,n8)
     coordtmp(3,8) = coor(3,n8)

     ! Find which material this element is made of
     imat = MatType(ielem)

     ! ... Initialize stuff
     NtN(:,:) = 0.0
     element_volume = 0.0

     ! ... Loop throught the gauss points in this element
     DO igpt = 1, ngpts

        ! ... Find the shape functions N and the derivative dN
        ! ... at gauss point = igpt
        N(:) = 0.0
        dN(:,:) = 0.0
        CALL get_shape(ri,N,dN,igpt)

        ! ... Find the volume of the element
        jac(1:3,1:3) = 0.0
        jacinv(1:3,1:3) = 0.0
        CALL get_jacobien(coordtmp,3,8,dN,jac,jacinv,detj,error)
        element_volume = element_volume + detj ! * w

        ! ... Integrate and assemble into NtN
        DO i = 1, 8
           DO j = 1, 8
              NtN(i,j) = NtN(i,j) + cap(imat)*rho(imat)*N(i)*N(j)*detj  ! * w
           ENDDO
        ENDDO

     ENDDO


     ! ... Local element capacitance matrix computed
     ! ... now assemble

     DO i = 1, 8  ! node1 loop
        gdof1 = Local2Global(ElConnVol(i,ielem))
        ldof1 = i
        DO j = 1, 8  ! node2 loop
           ldof2 = j
           gdof2 = Local2Global(ElConnVol(j,ielem))

           ! ... Make sure the matrix is exactly symmetric
           IF ( ldof1 <= ldof2 ) THEN
              tempval = NtN(ldof1,ldof2)
           ELSE
              tempval = NtN(ldof2,ldof1)
           ENDIF

           IF ( tempval == 0.0 ) THEN
              print*,'ZERO DETECTED IN LOCAL CAPACITANCE MATRIX!',gdof1,gdof2,myid
           END IF
           tempi = (gdof1-1)*diag_width
           DO k=1,diag_width
              ! ... if index k in col_index is empty, then fill index k 
              ! ... in col_index and values with the corresponding column
              ! ... index and capacitance matix entry
              IF(col_index(tempi+k).EQ.-1)THEN
                 !col_index(tempi+k)=gdof2-1
                 col_index(tempi+k)=gdof2
                 values(tempi+k)=tempval

                 EXIT
                 ! ... Ex., if, for a given row, column 5 needs be inserted, but
                 ! ... col_val(...k,k+1,..)=<...4,6...>, shift k+1,to k+2,... then
                 ! ... insert 5 at k+1.
              ELSEIF((col_index(tempi+k).LT.gdof2).AND.(col_index(tempi+k+1).GT.gdof2))THEN
                 DO l= diag_width-1, k+1, -1
                    tempj = col_index(tempi+l)
                    col_index(tempi+l+1)=tempj
                    tempval2 = values(tempi+l)
                    values(tempi+l+1) = tempval2
                 ENDDO
                 !col_index(tempi+k+1) = gdof2-1
                 col_index(tempi+k+1) = gdof2
                 values(tempi+k+1) = tempval
                 EXIT
                 ! ... if entry in k is less than given column number, but no entry in k+1
                 ! ... then enter column value into k+1
              ELSEIF((col_index(tempi+k).LT.gdof2).AND.(col_index(tempi+k+1).EQ.-1))THEN
                 values(tempi+k+1)=tempval
                 col_index(tempi+k+1)=gdof2
                 EXIT
                 ! ... if the column value entered into k is the same as the given
                 ! ... column value, then add the new entry into k in the values (=aval)
                 ! ... vector.
              ELSEIF(col_index(tempi+k).EQ.gdof2)THEN
                 values(tempi+k)= values(tempi+k) + tempval
                 EXIT
                 ! ... if column value in k is greater than the given column value, 
                 ! ... shift values k,k+1,... to k+1,k+2... and enter given column
                 ! ... value and capacitance value into col_index and values vectors.
              ELSEIF(col_index(tempi+k).GT.gdof2)THEN
                 DO l= diag_width-1, k, -1
                    tempj = col_index(tempi+l)
                    col_index(tempi+l+1)=tempj
                    tempval2 = values(tempi+l)
                    values(tempi+l+1) = tempval2
                 ENDDO
                 col_index(tempi+k) = gdof2
                 values(tempi+k) = tempval
                 EXIT
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  if(myid.eq.0) print*,'FINISHED FORMING GLB CAPACITANCE MATRIX FOR ELEMENT',ielem

  ! ... Determining the exact number of non-zero elements in the global mass matrix
  tempi=0
  DO i= 1,GNumNp*diag_width
     IF (col_index(i).NE.-1)THEN
        tempi = tempi + 1
     ENDIF
  ENDDO

  ! ... Allocating cval_temp, aval_temp and rp_temp
  ALLOCATE(rp_temp(1:GNumNp+1))
  ALLOCATE(cval_temp(1:tempi))
  ALLOCATE(aval_temp(1:tempi))
  tempj=tempi

  ! ... Populating aval-temp, cval_temp
  tempi=1
  DO i= 1, GNumNp*diag_width
     IF (col_index(i).NE.-1)THEN
        cval_temp(tempi) = col_index(i)
        aval_temp(tempi) = values(i)
        tempi = tempi + 1
     ENDIF
  ENDDO

  ! ... Populating rp_temp
  rp_temp(1) = 0
  DO i = 1, GNumNp
     tempi = 0
     ! ... loop through row i, look for non-zeros
     DO j = 1, diag_width
        IF (col_index((i-1)*diag_width+j).NE.-1)THEN
           tempi = tempi + 1
        ENDIF
     ENDDO
     rp_temp(i+1) = rp_temp(i) + tempi
  ENDDO

  ! ... cval_temp(i) = cval_temp(i)-1 per C standards for BLOCKSOLVE95
  DO i=1,tempj
     cval_temp(i)=cval_temp(i)-1
  ENDDO
  nnz = tempj
  nnz_temp = nnz

  DEALLOCATE(col_index)
  DEALLOCATE(values)

END SUBROUTINE LocThermCap_v3d8
