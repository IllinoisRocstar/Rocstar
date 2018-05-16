!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
SUBROUTINE readinp(ntime)
  
  USE meshdata
  
  IMPLICIT NONE
  
  INTEGER :: i
  INTEGER :: ntime

! -- local
  CHARACTER*200 :: keywd    ! keyword parameter
  INTEGER :: ios            ! io error

! -- Default

  NumNodeIO = 0
  IOformat = 1
  numvertx = 4
!
! -- Open Analysis Deck File

  OPEN(io_input,FILE='fractography3d.inp',STATUS='old',IOSTAT=ios)


  IF(ios.NE.0)THEN
     PRINT*, 'Unable to find fractography3d.inp'
     PRINT*, ' ...Trying RocfracControl.txt'
     
     OPEN(io_input,FILE='RocfracControl.txt',STATUS='old',IOSTAT=ios)
     IF(ios.NE.0)THEN
        PRINT*, 'Unable to find RocfracControl.txt'
        PRINT*, ' ...STOPPING'
        STOP
     ENDIF
  ENDIF

!  CALL COM_call_system('mkdir -p Rocin')
!
! -- Read Analysis Deck File

  REWIND io_input
10 READ(io_input,'(A)',IOSTAT=ios) keywd
  IF(ios.LT.0) THEN ! Negative ios means end-of-file
     PRINT*,' *END parameter not found - STOPPING'
     STOP
  ENDIF
  IF(keywd(1:4).EQ.'*END') THEN
     GOTO 40
  ELSE IF(keywd(1:7).EQ.'*PREFIX') THEN
     CALL PREFIX_SUB()
!!$     CALL system ('mkdir '//prefx(1:prefx_lngth))
     GOTO 10
  ELSE IF(keywd(1:7).EQ.'*MATVOL'.OR.keywd(1:13).EQ.'*HYPERELASTIC'.OR. &
       keywd(1:8).EQ.'*ELASTIC') THEN ! Read volumetric material props.
     CALL MATVOL_SUB()
     GOTO 10  
  ELSE IF(keywd(1:8).EQ.'*ELEMENT'.AND.keywd(1:15).NE.'*ELEMENT OUTPUT')THEN
     CALL ELEMENT_SUB(keywd)
     GOTO 10
!      ELSE IF(keywd(1:5).EQ.'*NODE') THEN
!         CALL NODEIO_SUB()
!         GOTO 10
  ELSE IF(keywd(1:9).EQ.'*IOFORMAT') THEN
     CALL IOFORMAT_SUB()
     GOTO 10
  ELSE
     GOTO 10   
  ENDIF

40 READ(io_input,'(A)',IOSTAT=ios) keywd
  IF(ios.LT.0) THEN ! Negative ios means end-of-file
     PRINT*,' *END parameter not found - STOPPING'
     STOP
  ENDIF
  IF(keywd(1:3).EQ.'*END') THEN
     GOTO 50
!      ELSE IF(keywd(1:3).EQ.'*BC') THEN
!         CALL BC_SUB()
!         GOTO 40
  ELSE IF(keywd(1:9).EQ.'*MESHSOFT') THEN
     CALL MESHSOFT_SUB()
     GOTO 40
  ENDIF
50 CONTINUE

  CLOSE(io_input)

  RETURN
END SUBROUTINE readinp

SUBROUTINE PREFIX_SUB()
  
  USE meshdata
  
  IMPLICIT NONE
  
  READ(io_input,'(a20)') prefx
  prefx_lngth = LEN_TRIM(prefx)
  
  RETURN
END SUBROUTINE PREFIX_SUB

SUBROUTINE NRUN_SUB(ntime)
  
  USE meshdata
  
  IMPLICIT NONE
  
  INTEGER :: ntime,ii
  REAL*8 :: iaux
  READ(io_input,*) iaux, numvertx
  
  RETURN
END SUBROUTINE NRUN_SUB

SUBROUTINE MATVOL_SUB()
  
  USE meshdata
  
  IMPLICIT NONE
  
  INTEGER :: i              ! loop counter
  
  INTEGER :: numat_vol      ! number of volumetric materials
  
  REAL*8 :: E, xnu, rho, alpha
  
  READ(io_input,*) numat_vol
  
  cd_fastest = 0.d0
  DO i = 1, numat_vol
     READ(io_input,*) E, xnu, rho, alpha
     cd_fastest = MAX( cd_fastest, &
          SQRT(E*(1.d0-xnu)/rho/(1.d0+xnu)/(1.d0-2.d0*xnu)) )
  ENDDO
  
  RETURN
END SUBROUTINE MATVOL_SUB

SUBROUTINE BC_SUB()
  
  USE meshdata
  
  IMPLICIT NONE
  
  INTEGER :: iaux
  INTEGER :: i
  
  DO i = 1, 32
     READ(io_input,*) iaux, bc_conditions(i)%b1, bc_conditions(i)%b2, bc_conditions(i)%b3, &
          bc_conditions(i)%bc1, bc_conditions(i)%bc2, bc_conditions(i)%bc3
  ENDDO
  
  RETURN
END SUBROUTINE BC_SUB

SUBROUTINE NODEIO_SUB()
  
  USE meshdata
  
  IMPLICIT NONE
  
  INTEGER :: i              ! loop counter
  
  READ(io_input,*) NumNodeIO
  
  ALLOCATE(NodeIO(1:NumNodeIO))
  
  DO i = 1, NumNodeIO
     READ(io_input,*) NodeIO(i)
  ENDDO
  
  RETURN
END SUBROUTINE NODEIO_SUB

SUBROUTINE MESHSOFT_SUB()
  
  USE meshdata
  
  IMPLICIT NONE
  
  CHARACTER*1 :: chr
  
  READ(io_input,*) chr
  
  iansys   = 0 ! set default -no-
  ipatran  = 0 
  itetmesh = 0
  ipatcohin   = 0
  itetcohin   = 0
  
  IF(chr.EQ.'T'.OR.chr.EQ.'t')THEN
     itetmesh = 1
  ELSE IF(chr.EQ.'A'.OR.chr.EQ.'a')THEN
     iansys = 1
  ELSE IF(chr.EQ.'P'.OR.chr.EQ.'p')THEN
     ipatran = 1
  ELSE IF(chr.EQ.'C'.OR.chr.EQ.'c')THEN
     READ(io_input,*) chr
     IF(chr.EQ.'P'.OR.chr.EQ.'p')THEN
        ipatcohin = 1
     ELSE IF(chr.EQ.'T'.OR.chr.EQ.'t')THEN
        itetcohin = 1
     ENDIF
  ELSE
     PRINT*,' ERROR: MESHING PACKAGE NOT SUPPORTED'
     PRINT*,'STOPPING'
     STOP
  ENDIF
  
  RETURN
END SUBROUTINE MESHSOFT_SUB

SUBROUTINE IOFORMAT_SUB()
  
  USE meshdata
  
  IMPLICIT NONE

! 0 = binary
! 1 = ascii

  READ(io_input,*) IOformat
  
  RETURN
END SUBROUTINE IOFORMAT_SUB
SUBROUTINE ELEMENT_SUB(keywd)
  
  USE meshdata
  
  IMPLICIT NONE
  
  CHARACTER(len=200) :: keywd
  INTEGER :: k1, k2
  
  CHARACTER(len=16) :: ElType
  
  CALL locchr(keywd,'TYPE                      ',4,8,k1,k2)
  
  ElType = keywd(k1:k2)
  
  SELECT CASE (TRIM(ElType))
  CASE ('V3D4')
     numvertx = 4
  CASE ('V3D4NCC')
     numvertx = 4
  CASE ('V3D4N')
     numvertx = 4
  CASE ('V3D10R')
     numvertx = 10
  CASE ('V3D10BBAR')
     numvertx = 10
  CASE ('V3D10')
     numvertx = 10
  CASE ('V3D8')
     numvertx = 8
  CASE ('V3D8ME')
     numvertx = 8
  CASE default
     PRINT*,' ERROR:'
     PRINT*,'*ELEMENT TYPE NOT FOUND'
     STOP
  END SELECT
  RETURN
END SUBROUTINE ELEMENT_SUB

