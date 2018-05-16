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
   SUBROUTINE locchr(text,varna,lvari,kpos0,kpos1,kpos2)

!!****f* Rocfrac/Rocfrac/Source/locchr.f90/locchr
!!
!!  NAME
!!     locchr
!!
!!  FUNCTION
!!
!!      Locates the keyword value after a keyword
!!
!!  INPUTS
!!    text -- character string
!!    varna -- variable name to search for
!!    lvari -- length of the variable name
!!    kpos0 -- initial position in 'text' so start looking for varna
!!  
!!  OUTPUT
!!    kpos1 -- Start of keyword value in string
!!    kpos2 -- End of keywork value in string
!!
!!****

    IMPLICIT NONE

    CHARACTER(len=200) :: text
    CHARACTER(len=26) :: varna

    INTEGER :: kpos0,kpos1,kpos2
! length of the input variable
    integer :: lvari

    INTEGER :: kposi
    INTEGER :: key,k

    kposi = kpos0-1
    IF(kposi.LT.0) kposi = 0
    DO
       kposi=kposi+1
       IF(kposi+lvari-1.GT.200) THEN
          kpos2 = -1
          RETURN
       ENDIF
       IF(text(kposi:kposi+lvari-1).EQ.varna(1:lvari)) EXIT
    ENDDO

    kpos1 = kposi
    DO
       kpos1=kpos1+1
       IF(text(kpos1:kpos1).EQ.'=') EXIT
    ENDDO

    DO
       kpos1 = kpos1 + 1
       IF(text(kpos1:kpos1).NE.' ') EXIT
    ENDDO
    kpos2 = kpos1

    key = 0

    DO
       kpos2 = kpos2 + 1
       IF(kpos2.GE.120) EXIT
!!! if(text(kpos2:kpos2).ne.' ') key=1
       IF(text(kpos2:kpos2).EQ.',') EXIT
       IF(text(kpos2:kpos2).EQ.' '.AND.&
            text(kpos2+1:kpos2+1).EQ.' ') EXIT
    ENDDO
!!! kpos1=kpos1+1
    kpos2=kpos2-1
    DO k=1,10
       IF(text(kpos1-1+k:kpos1-1+k).NE.' ') EXIT
    END DO
    kpos1 = kpos1 + k - 1
    DO k=1,10
       IF(text(kpos2+1-k:kpos2+1-k).NE.' '.AND. &
            text(kpos2+1-k:kpos2+1-k).NE.',') EXIT
    END DO

    kpos2=kpos2-k+1
    RETURN

  END SUBROUTINE locchr

  SUBROUTINE conchr(text,varna,lvari,kpos0,key)

!!****f* Rocfrac/Rocfrac/Source/locchr.f90/conchr
!!
!!  NAME
!!     conchr
!!
!!  FUNCTION
!!
!!      To determine if a control deck keyword is specified
!!
!!  INPUTS
!!    text -- character string
!!    varna -- variable name to search for
!!    lvari -- length of the variable name
!!    kpos0 -- initial position in 'text' so start looking for varna
!!  
!!  OUTPUT
!!     key --  0 = no, 1 = yes
!!
!!****


    IMPLICIT NONE
! In
    INTEGER :: kpos0, lvari
    CHARACTER(len=200) :: text
    CHARACTER(len=26 ) :: varna
! Out
    INTEGER :: key ! 1=found, 0=not found

    INTEGER :: lll, kposi

    key = 0

    CALL dtext(text,lll)

    kposi = kpos0 - 1

    IF(kposi.LT.0) kposi = 0

    DO
       kposi = kposi + 1

! Found keyword
!
       IF(text(kposi:kposi+lvari-1).EQ.varna(1:lvari)) EXIT
!
! Keyword not found

       IF(kposi+lvari-1.GT.lll) RETURN

    ENDDO

    key=1
    RETURN
  END SUBROUTINE conchr

  SUBROUTINE dtext(text,lll)

!!****f* Rocfrac/Rocfrac/Source/locchr.f90/dtext
!!
!!  NAME
!!     dtext
!!
!!  FUNCTION
!!
!!     To determine the string length
!!
!!  INPUTS
!!     text -- character string
!!  
!!  OUTPUT
!!     lll -- length of string
!!
!!****

    IMPLICIT NONE
! In
    CHARACTER(len=200) :: text
! Out
    INTEGER :: lll

    INTEGER :: mlen

    mlen = lll
    IF(lll.LT.1) mlen = 200
    mlen = 200
    DO lll = mlen, 1, -1
       IF(text(lll:lll).NE.' ') RETURN
    END DO
    lll=1
    RETURN
  END SUBROUTINE dtext

  SUBROUTINE dchar(char,key)

!!****f* Rocfrac/Rocfrac/Source/locchr.f90/dtext
!!
!!  NAME
!!     dtext
!!
!!  FUNCTION
!!
!!     Converts a character string to an integer
!!
!!  INPUTS
!!     char -- character string
!!  
!!  OUTPUT
!!     key -- integer 
!!
!!****

    
    IMPLICIT NONE
    
    INTEGER :: key
    CHARACTER(len=16) :: char
    
    INTEGER :: k
    
    key = 0
    DO k = 1, 16
       IF(char(k:k).NE.'0'.AND.char(k:k).NE.'1'.AND. &
            char(k:k).NE.'2'.AND.char(k:k).NE.'3'.AND. &
            char(k:k).NE.'4'.AND.char(k:k).NE.'5'.AND. &
            char(k:k).NE.'6'.AND.char(k:k).NE.'7'.AND. &
            char(k:k).NE.'8'.AND.char(k:k).NE.'9'.AND. &
            char(k:k).NE.' ') RETURN
    END DO
    READ(char,'(I16)') key
    RETURN
    
  END SUBROUTINE dchar

  SUBROUTINE rchar(char,key)

!!****f* Rocfrac/Rocfrac/Source/locchr.f90/rchar
!!
!!  NAME
!!     rchar
!!
!!  FUNCTION
!!
!!     Converts a character string to a real
!!
!!  INPUTS
!!     char -- character string
!!  
!!  OUTPUT
!!     key -- real
!!
!!****

    
    IMPLICIT NONE
    
    real*8 :: key
    CHARACTER(len=16) :: char
    
    INTEGER :: k
    
    key = 0
!
! Don't enforce that key includes a '.' real number,
! allow for integer input.

!!$    DO k = 1, 16
!!$       IF(char(k:k).NE.'0'.AND.char(k:k).NE.'1'.AND. &
!!$            char(k:k).NE.'2'.AND.char(k:k).NE.'3'.AND. &
!!$            char(k:k).NE.'4'.AND.char(k:k).NE.'5'.AND. &
!!$            char(k:k).NE.'6'.AND.char(k:k).NE.'7'.AND. &
!!$            char(k:k).NE.'8'.AND.char(k:k).NE.'9'.AND. &
!!$            char(k:k).NE.' '.AND.char(k:k).NE.'.') RETURN
!!$    END DO

    READ(char,*) key
    RETURN
    
  END SUBROUTINE rchar

