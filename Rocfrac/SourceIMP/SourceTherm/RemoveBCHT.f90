
!!****
!!
!!  NAME
!!     RemoveBCHT_reff
!!
!!  FUNCTION
!!     This subroutine removes values from the load vector
!!     that are associated with perscribed boundary
!!     conditions.
!!
!!  INPUTS
!!     ndim -- The size of one dimension of the global Meff matrix
!!     nstart -- The global index of the first row assigned to this processor
!!     pbar -- The part of the global load vector that is assigned to this processor
!!     newndim -- The size of the global Meff matrix after BCs have been removed
!!
!!  OUTPUTS
!!     newpbar -- The part of the global load vector that is assigned to this processor after BCs have been removed
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE RemoveBCHT_reff(nstart,ndim,pbar,newndim,newpbar)

  USE removeBCs_global
  USE comp_row_global
  USE implicit_global
  USE Precision

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: i, m, counter1, counter2

  INTEGER :: ndim, newndim, nstart
  REAL(kind=wp), DIMENSION(ndim) :: pbar
  REAL(kind=wp), DIMENSION(newndim) :: newpbar
  

  ! Construct new vector
  counter2 = 0
  DO i = 1, ndim
     counter1 = 0
     DO m = 1, GNumTemp
        IF(GTempBC(m) == i + nstart - 1) THEN
           counter1 = 1
        ENDIF
     ENDDO
     IF(counter1 == 0) THEN
        counter2 = counter2 + 1
        newpbar(counter2) = pbar(i)
     ENDIF
  ENDDO
  

END SUBROUTINE RemoveBCHT_reff









!!****
!!
!!  NAME
!!     removeBCs_newa
!!
!!  FUNCTION
!!     This subroutine takes the acceleration vector that
!!     does not include boundary conditions and expands it
!!     back into the acceleration vector with rows associated
!!     with boundary conditions in it.
!!
!!  INPUTS
!!     ndim -- The size of one dimension of the global Meff matrix
!!     nstart -- The global index of the first row assigned to this processor
!!     newa -- The part of the global acceleration vector assigned to this proc without rows associated with boundary conditions
!!     newndim -- The size of the global Meff matrix after BCs have been removed
!!
!!  OUTPUTS
!!     a -- The part of the global acceleration vector assigned to this proc with rows associated with boundary conditions
!!
!!  USES
!!     none
!!
!!****

SUBROUTINE RemoveBCHT_newTemp(nstart,ndim,a,newndim,newa)

  USE removeBCs_global
  USE comp_row_global
  USE implicit_global
  USE Precision

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: i, m, counter1, counter2, counter3

  INTEGER :: ndim, newndim, nstart
  REAL(kind=wp), DIMENSION(ndim) :: a
  REAL(kind=wp), DIMENSION(newndim) :: newa
  

  ! Construct new vector
  counter2 = 0
  DO i = 1, newndim
     counter3 = 0
     DO m = 1, GNumTemp
        IF((GTempBC(m) <= i + counter3 + nstart - 1).AND.(GTempBC(m) >= nstart)) THEN
           counter3 = counter3 + 1
        ENDIF
     ENDDO
     a(i+counter3) = newa(i)
     !if(myid==2) print*,myid,i,counter3,i+counter3
  ENDDO
  

END SUBROUTINE RemoveBCHT_newTemp
