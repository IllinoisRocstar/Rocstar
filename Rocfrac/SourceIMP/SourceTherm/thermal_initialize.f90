

!!****
!!
!!  NAME
!!     thermal_initialize
!!
!!  FUNCTION
!!     Creates all of the data structures and matrices needed
!!     by the implicit thermal solver
!!
!!  INPUTS
!!     none
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     MPI
!!
!!****

SUBROUTINE thermal_initialize(global)

  USE implicit_global
  USE ROCSTAR_RocFracComm

  IMPLICIT NONE

  TYPE(ROCFRAC_GLOBAL) :: global
  INTEGER :: i, j, m

  ! ... Set up the global to local node mapping
  ! ... GNumNP is the number of nodes in the entire model
  ! ... global%NumNP is the number of node points that
  ! ... this process knows about.
  ! ... Node points that are represented by GNumNP but 
  ! ... are not known ot this process (NumNP) are represented
  ! ... in the Global2Local array by the value -1.

  ALLOCATE(Global2Local(1:GNumNP))
  Global2Local(:) = -1
  DO i = 1, global%NumNP
     Global2Local(Local2Global(i)) = i
  ENDDO



  ! ... Initialize MPI communications

  ! ... InitComm1 sets up variables that aid in communication
  ! ... of boundary nodes on one process to the process that 
  ! ... owns that node.  Usefull in constructing capacitance and 
  ! ... stiffness matrices.

  ! ... InitComm2 sets up variables that aid in communication
  ! ... of nodes owned by one process to the processes that contain 
  ! ... those nodes in their partition.  Usefull in constructing  
  ! ... internal load vector.

  IF(myid==0) PRINT*,'INITIALIZING COMMUNICATIONS'
  IF (nprocs > 1) THEN
     CALL InitComm1(global)
     CALL InitComm2(global)
  ENDIF



  ! ... Set up some constants
  ! ... nstart_ktc is the 1st node number that this process is
  ! ... responsible for. nrows_ktc is the number of rows of the 
  ! ... capacitance and stiffness matrices owned by this process.

  nstart_ktc = GNumNp
  DO m = 1, GNumNp
     DO i = 1, global%NumNp
        IF (Local2Global(i) == m) THEN
           IF (NodeProc(i) == myid) THEN
              nstart_ktc = MIN(nstart_ktc,Local2Global(i))
           ENDIF
        ENDIF
     ENDDO
  ENDDO

  nrows_ktc = LNumNp
  
  ! ... Set up the nodal DOF vectors
  ! ... Need to do more research on this

  ALLOCATE(node_flag(1:global%NumNp,1:3))
  ALLOCATE(boundary_value(1:global%NumNp,1:3))
  node_flag(1:global%NumNp,1) = 0
  DO i = 1, global%NumNdsBCHT
     IF ( global%BCFlagHT(2,i) == 0 ) THEN
        node_flag(global%BCFlagHT(1,i),1:3) = 8
        boundary_value(global%BCFlagHT(1,i),1:3) = global%BCvalueHT(1,i)
     ENDIF
     IF ( global%BCFlagHT(2,i) == 1 ) THEN
        node_flag(global%BCFlagHT(1,i),1:3) = 7
        boundary_value(global%BCFlagHT(1,i),1:3) = global%BCvalueHT(1,i)
     ENDIF
  ENDDO


  ! ... Create the global thermal capacitance matrix
  CALL GlbThermCap(global)
  
  ! ... Create the global thermal stiffness matrix  
  CALL GlbThermStiff(global) 



END SUBROUTINE thermal_initialize
