!-----------------------------------------------------GET_MAT_CONDUCTION
SUBROUTINE ConductivityTensor(NumMatVol,k,dmat)
! *--------------------------------------------------------------------*
! |                                                                    |
! |    Returns the material thermal conduction matrx, given            |
! |    the material properties.  Currently written for a material      |
! |    with isotropic thermal properties.                              | 
! |                                                                    |
! |    <NumMatVol> : number of materials in the model                  |
! |    <k>         : isotropic thermal conduction coefficient          |
! |    <dmat>      : 3 x 3 matrix of the conduction tensor             |
! |                                                                    |
! *--------------------------------------------------------------------*

  IMPLICIT NONE
  
  REAL*8 , DIMENSION(NumMatVol,1:3,1:3) :: dmat
  REAL*8  :: k
  REAL*8  :: d1111,d1122,d1212
  INTEGER :: NumMatVol, i
  
  
  dmat(:,:,:) = 0.d0

  DO i = 1,NumMatVol 
     dmat(i,1,1) = k*1.d0
     dmat(i,2,2) = k*1.d0
     dmat(i,3,3) = k*1.d0
  ENDDO

END SUBROUTINE ConductivityTensor
