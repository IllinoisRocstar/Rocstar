!!****
!!
!!  NAME
!!     RockOut_8node
!!
!!  FUNCTION
!!     This subroutine writes an ASCII file of 8 node elements and associated temperature
!!     data at every timestep.
!!
!!  INPUTS
!!     glb             -- The global pointer variable
!!     CurrentTimeStep -- The current time step (for file naming)
!!
!!  OUTPUTS
!!     newnrows -- The number of rows assigned to this proc after BCs have been removed
!!     newnstart -- The global index of the first row assigned to this proc after BCs have been removed
!!     newndim -- The size of the global Meff matrix after BCs have been removed
!!
!!  USES
!!     none
!!
!!****
  SUBROUTINE RockOut_8node(glb,CurrentTimeStep)
     
    USE implicit_global
    USE comp_row_global
    USE Precision
    
    IMPLICIT NONE

    TYPE(ROCFRAC_GLOBAL) :: glb
    CHARACTER(LEN=4) :: chr1
    INTEGER :: CurrentTimeStep
    INTEGER :: i,j,k
    
    PRINT*,'ROCK OUT!!!! MEEDLEY-MEEDLEY!!!'

       WRITE(chr1,'(i4.4)') 0!CurrentTimeStep
!       i = CurrentTimeStep
 i=0      
       print*,i,glb%NumNP
       ! WRITE(ichr2,'(I4.4)') myid
    PRINT*,'ROCK OUT!!!! before IO'
       
       !          OPEN(9940+i,FILE='Rocfrac/Rocout/RocTherm.'//ichr1//'.'//ichr2,POSITION='APPEND')
       OPEN(900+i,FILE='Rocfrac/Rocout/RocTherm.'//chr1//'.dat',POSITION='APPEND')
!        print*,'1'
!        WRITE(900+i,*) 'TITLE= "Temperature"'
!        print*,'2'
!        WRITE(900+i,*) 'VARIABLES= "X", "Y", "Z", "Temperature"'
!        print*,'3'
!        WRITE(900+i,*) 'ZONE N=',glb%NumNP,', E=',glb%NumElVol,', DATAPACKING=POINT, ZONETYPE=FEBRICK'
!     PRINT*,'ROCK OUT!!!! before meshcoor'
!        DO j = 1, glb%NumNP
!           write(900+i,*) glb%meshcoor(j,1),glb%meshcoor(j,2),glb%meshcoor(j,3),glb%Temperature(j)
!        ENDDO
!        WRITE(900+i,*) " "
!        DO j = 1, glb%NumElVol
!           write(900+i,*) (glb%ElConnVol(j,k),k=1,8)
!       ENDDO
        WRITE(900+i,*) ' '
        WRITE(900+i,*) CurrentTimeStep
        do j=1,glb%NumNP
           WRITE(900+i,*) j,glb%Temperature(j)
        end do
       close(900+i)
    
    
  END SUBROUTINE RockOut_8node
