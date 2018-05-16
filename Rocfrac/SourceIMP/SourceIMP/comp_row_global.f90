
!!****
!!
!!  NAME
!!     comp_row_global
!!
!!  FUNCTION
!!     The variables within this module can store an entire matrix in
!!     compressed row storage.  These variables are obviously
!!     inherantly allocatable and constantly change size depending on
!!     the number of nonzeros in the matrix.  Since FORTRAN90 does not
!!     allow allocatable arrays to be passed to functions or subroutines,
!!     they must be places as such in a global module.
!!
!!  INPUTS
!!     none
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     none
!!
!!****

MODULE comp_row_global
  
  USE Precision
  
  IMPLICIT NONE


  REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: aval_temp
  INTEGER,ALLOCATABLE,DIMENSION(:) :: cval_temp
  INTEGER,ALLOCATABLE,DIMENSION(:) :: rp_temp
  INTEGER :: nnz_temp

END MODULE comp_row_global
