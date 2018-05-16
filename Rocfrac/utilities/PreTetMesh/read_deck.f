      SUBROUTINE read_deck()

      USE meshdata

! -- local
      CHARACTER*200 :: keywd    ! keyword parameter

      INTEGER :: ios            ! io error
      INTEGER :: ierr           ! mpi error

      INTEGER :: i              ! loop counter
!
! -- Open Analysis Deck File

      ios = 0

!      OPEN(io_input,FILE='fractography3d.inp',STATUS='old',IOSTAT=ios)
      OPEN(io_input,FILE='RocfracControl.txt',STATUS='old',IOSTAT=ios,FORM='formatted')
      IF(ios.NE.0)THEN
         PRINT*, 'Warning: Unable to find RocfracControl.txt'
         PRINT*, ' '
         PRINT*,' Entering manual mode...'
         PRINT*,'  (1) Enter the prefix of the patran neutral file'
         PRINT*,'    NOTE: the patran neutral file should be called: <prefix>.pat'
         READ(*,'(a20)') prefx
         prefx_lngth = LEN_TRIM(prefx)
         GOTO 1
      ENDIF
!
! -- Read Analysis Deck File

      REWIND io_input

 1    READ(io_input,'(A)',IOSTAT=ios) keywd
      IF(ios.LT.0) THEN ! Negative ios means end-of-file
         PRINT*,' *END parameter not found - STOPPING'
         STOP
      ENDIF
!
! Comment field

      IF(keywd(1:2).EQ.'**'.OR.keywd(1:1).NE.'*')THEN
         GOTO 1
      ENDIF      

      IF(keywd(1:4).EQ.'*END') THEN
         GOTO 3
      ELSE IF(keywd(1:7).EQ.'*PREFIX') THEN
         CALL PREFIX_SUB()
         GOTO 1
      ELSE
         GOTO 1
      ENDIF

 3    CONTINUE

      CLOSE(io_input)

      RETURN

      END

      SUBROUTINE PREFIX_SUB()

      USE meshdata

      IMPLICIT NONE
         
      READ(io_input,'(a20)') prefx
      print*,prefx
      prefx_lngth = LEN_TRIM(prefx)
      
      RETURN
      END
