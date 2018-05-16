
!!****
!!
!!  NAME
!!     LocThermStiff_v3d8
!!
!!  FUNCTION
!!     This subroutine constructs the process local thermal stiffness
!!     matrix for volumetric stess 8 node, 3d elements.  The resulting
!!     matrix is put into the variables aval_kt, cval_kt and rp_kt
!!     in comp_row_global.
!!
!!  INPUTS
!!     coor -- Coordinates of the nodes
!!     ElConnVol -- The connectivity table
!!     dmat -- Material stiffness matrix
!!     numnp -- Number of nodes
!!     NumEl -- Number of elements
!!     MatType -- Mapping of materials to elements
!!     NumMatType -- Total number of materials
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE LocThermStiff_v3d8(coor,ElConnVol,dmat, & 
     numnp,NumEL,MatType,NumMatType)


  USE comp_row_global 
  USE implicit_global 
  USE Precision 

  IMPLICIT NONE 

  ! ... number of node points that this process knows of
  INTEGER :: numnp       
  ! ... number of elements that this process knows of
  INTEGER :: NumEl       
  ! ... number of materials in the model  
  INTEGER :: NumMatType 
  ! ... coordinate array 
  REAL*8, DIMENSION(1:3,1:numnp) :: coor 
  ! ... connectivity table for Brick  
  INTEGER, DIMENSION(1:8,1:NumEl) :: ElConnVol 
  INTEGER, DIMENSION(1:NumEl) :: MatType 

  ! ... Local variables 
  ! ... node numbers 
  INTEGER :: n1,n2,n3,n4,n5,n6,n7,n8
  ! ... coordinate array
  REAL*8, DIMENSION(1:3,1:8) :: coord 
  
  REAL*8 :: element_volume
  ! ... element index, material index, Gauss point index
  INTEGER :: ielem, imat, igpt
  ! ... shape function, derivatives, jacobian and products
  REAL*8 :: N(8), dN(8,3), jac(3,3), jacinv(3,3), NkN(8,8), dmat(3,3),Bmat(8,3)

  ! ... Gauss point coordinates
  REAL*8, DIMENSION(1:3,1:8) :: ri = RESHAPE( & 
       (/-0.577350269189626,-0.577350269189626,-0.577350269189626, & 
       0.577350269189626,-0.577350269189626,-0.577350269189626, & 
       0.577350269189626, 0.577350269189626,-0.577350269189626, & 
       -0.577350269189626, 0.577350269189626,-0.577350269189626, & 
       -0.577350269189626,-0.577350269189626, 0.577350269189626, & 
       0.577350269189626,-0.577350269189626, 0.577350269189626, & 
       0.577350269189626, 0.577350269189626, 0.577350269189626, & 
       -0.577350269189626, 0.577350269189626, 0.577350269189626/),(/3,8/) ) 

  ! ... determinant of the jacobian
  REAL*8 :: detj
  LOGICAL :: error, debug 
  INTEGER,parameter :: ngpts = 8 
  INTEGER :: i,k,j,l

  INTEGER :: m 
  INTEGER :: dof1, dof2 
  INTEGER :: counter 
  REAL(kind=wp) :: tempval

  ! ... Initialize the K matrix 
  ! ... find number of non-zeros in the process local thermal stiffness matrix
  CALL KtNumNnz(NumNp,NumEl,ElConnVol,nnz_kt) 
  ALLOCATE(rp_kt(1:GNumNp+1)) 
  ALLOCATE(cval_kt(1:nnz_kt)) 
  ALLOCATE(aval_kt(1:nnz_kt))
  rp_kt(:) = 0
  cval_kt(:) = 0
  aval_kt(:) = 0.0

  ! ... fill the cval_kt and rp_kt vectors
  CALL KtInitCRS(NumNp,NumEl,ElConnVol,nnz_kt,rp_kt,cval_kt,aval_kt) 

  
  ! ... Loop through the elements
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

     ! ... x,y,z coordinates of node n#.
     ! ... local node 1
     coord(1,1) = coor(1,n1)
     coord(2,1) = coor(2,n1)
     coord(3,1) = coor(3,n1)
     ! ... local node 2
     coord(1,2) = coor(1,n2)
     coord(2,2) = coor(2,n2)
     coord(3,2) = coor(3,n2)
     ! ... local node 3
     coord(1,3) = coor(1,n3)
     coord(2,3) = coor(2,n3)
     coord(3,3) = coor(3,n3)
     ! ... local node 4
     coord(1,4) = coor(1,n4)
     coord(2,4) = coor(2,n4)
     coord(3,4) = coor(3,n4)
     ! ... local node 5
     coord(1,5) = coor(1,n5)
     coord(2,5) = coor(2,n5)
     coord(3,5) = coor(3,n5)
     ! ... local node 6
     coord(1,6) = coor(1,n6)
     coord(2,6) = coor(2,n6)
     coord(3,6) = coor(3,n6)
     ! ... local node 7
     coord(1,7) = coor(1,n7)
     coord(2,7) = coor(2,n7)
     coord(3,7) = coor(3,n7)
     ! ... local node 8
     coord(1,8) = coor(1,n8)
     coord(2,8) = coor(2,n8)
     coord(3,8) = coor(3,n8)

     ! ... Find which material this element is made of
     imat = MatType(ielem)

     ! ... Initialize stuff
     NkN(:,:) = 0.0
     element_volume = 0.0

     ! ... Loop throught the gauss points in this element
     DO igpt = 1, 8

        ! ... Find the shape functions N and the derivative dN
        ! ... at gauss point = igpt
        N(:) = 0.0
        dN(:,:) = 0.0
        CALL get_shape(ri,N,dN,igpt)

        ! ... Find the volume of the element
        jac(1:3,1:3) = 0.0
        jacinv(1:3,1:3) = 0.0
        detj = 0

        CALL get_jacobien(coord,3,8,dN,jac,jacinv,detj,error)
        element_volume = element_volume + detj ! * w


        Bmat(:,:) = 0.0
        do i =1,8
           Bmat(i,:) = matmul(dN(i,:),jacinv)
        enddo
           

        ! ... Integrate and assemble into NkN
        DO i = 1, 8
           DO j = 1, 8
              DO k = 1, 3
                 ! ... conductivity tensor 3x3
                 DO l = 1,3 
                    NkN(i,j) = NkN(i,j) + dmat(k,l) * Bmat(i,k) * Bmat(j,k) *detj  ! * w


                 ENDDO
              ENDDO
           ENDDO
       ENDDO

     ENDDO




     !************************************************************** 
     ! Assemble local K (i.e. matrix) into Global stiffness matrix 
     !************************************************************** 
     ! ... Loop over the first node 
     DO i = 1, 8
        ! ... Find the corresponding global node
        dof1 = Local2Global(ElConnVol(i,ielem))

        ! ... Loop over second node
        DO j = 1, 8

           ! ... Find the corresponding global node
           dof2 = Local2Global(ElConnVol(j,ielem))

           ! Make sure the matrix is exactly symmetric by using only values in the upper triangle of NkN
           IF ( i <= j ) THEN
              tempval = NkN(i,j)
           ELSE
              tempval = NkN(j,i)
           ENDIF

           ! ... Place the nonzeros into the K matrix
           IF (tempval /= 0.0) THEN  
              counter = 0
              DO m = rp_kt(dof1)+1, rp_kt(dof1+1)
                 IF (cval_kt(m) == dof2-1) THEN
                       aval_kt(m) = aval_kt(m) + tempval

                    counter = 1
                 ENDIF

              ENDDO
              if(counter==0) print*,'WARNING:  Unable to add value to K matrix at (',dof1,',',dof2,') on processor ',myid,rp_kt(dof1)+1,rp_kt(dof1+1)
           ELSE
              print*,myid,'ZERO DETECTED IN LOCAL STIFFNESS MATRIX!',dof1,dof2
           END IF
           
        ENDDO
     ENDDO
  ENDDO

  ! ... the arrays are stored in global temporary arrays because pointers of unknown
  ! ... size can't be passed as arguments
  ALLOCATE(rp_temp(1:GNumNp+1),cval_temp(1:nnz_kt),aval_temp(1:nnz_kt)) 
  nnz_temp = nnz_kt 
  rp_temp = rp_kt 
  cval_temp = cval_kt 
  aval_temp = aval_kt 

  DEALLOCATE(rp_kt,cval_kt,aval_kt) 

END SUBROUTINE LocThermStiff_v3d8
      
    
