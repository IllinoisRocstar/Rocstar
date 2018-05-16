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
!******************************************************************************
!
! Purpose: store x-,y-,z-coordinates of grid nodes in binary format.
!
! Notes: only the finest grid is stored.
!
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE WRITE_GRDB(iblock,ni,nj,nk,x,y,z,F_name)

  IMPLICIT NONE
  INTEGER :: iblock,ni,nj,nk
  DOUBLE PRECISION :: time,x(ni,nj,nk),y(ni,nj,nk),z(ni,nj,nk)
  INTEGER :: i,j,k
  CHARACTER (LEN=500) :: F_name,filnam
  INTEGER :: IF_GRID = 17

! open file; write time stamp

  IF (ABS(iblock) == 1) THEN
    WRITE(filnam,11) F_name(1:INDEX(F_name,' ')-1)
    OPEN (IF_GRID,file=filnam,form='unformatted',status='unknown')
    time = 0.0
    WRITE(IF_GRID) time
  ENDIF

! write block number and dimensions

  WRITE(IF_GRID) ABS(iblock),ni-1,nj-1,nk-1

! write coordinates

  WRITE(IF_GRID) &
    (((x(i,j,k), i=1,ni), j=1,nj), k=1,nk),&
    (((y(i,j,k), i=1,ni), j=1,nj), k=1,nk),&
    (((z(i,j,k), i=1,ni), j=1,nj), k=1,nk)

! close file

  IF (iblock < 0) THEN
    CLOSE(IF_GRID)
  ENDIF

11 FORMAT(a)
 
END SUBROUTINE


