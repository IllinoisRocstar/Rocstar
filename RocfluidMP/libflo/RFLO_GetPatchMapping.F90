! *********************************************************************
! * Rocstar Simulation Suite                                          *
! * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
! *                                                                   *
! * Illinois Rocstar LLC                                              *
! * Champaign, IL                                                     *
! * www.illinoisrocstar.com                                           *
! * sales@illinoisrocstar.com                                         *
! *                                                                   *
! * License: See LICENSE file in top level of distribution package or *
! * http://opensource.org/licenses/NCSA                               *
! *********************************************************************
! *********************************************************************
! * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
! * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
! * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
! * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
! * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
! * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
! * Arising FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
! * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
! *********************************************************************
!******************************************************************************
!
! Purpose: build mapping matrix from current to source patch.
!
! Description: the mapping is from i,j,k of the current patch (and its
!              dummy cells/nodes to is,js,ks of the source region and patch.
!              Hence:
!
!              is = i*mapMat(1,1) + j*mapMat(1,2) + k*mapMat(1,3) + mapMat(1,4)
!              js = i*mapMat(2,1) + j*mapMat(2,2) + k*mapMat(2,3) + mapMat(2,4)
!              ks = i*mapMat(3,1) + j*mapMat(3,2) + k*mapMat(3,3) + mapMat(3,4)
!
! Input: lb          = block face (current patch)
!        lbs         = block face (source patch)
!        l1SrcDir    = orientation of l1-direction on source patch wrp.
!                      to current patch (= 1 or -1)
!        l2SrcDir    = orientation of l2-direction on source patch wrp.
!                      to current patch (= 1 or -1)
!        align       = alignment between l1- and l2-directions on source
!                      and on current patch (true or false)
!        i/j/kdir    = direction: +-1 for coordinate normal to current patch,
!                      0 otherwise
!        i/j/kdirSrc = direction: +-1 for coordinate normal to source patch,
!                      0 otherwise
!        i/j/kbeg    = start indices of current patch in i-,j-,k-direction
!        i/j/kend    = end indices of current patch in i-,j-,k-direction
!        i/j/kbegSrc = start indices of source patch in i-,j-,k-direction
!        i/j/kendSrc = end indices of source patch in i-,j-,k-direction.
!
! Output: mapMat(3,4) = mapping matrix.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_GetPatchMapping.F90,v 1.3 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_GetPatchMapping( lb,lbs,l1SrcDir,l2SrcDir,align, &
                                 idir,jdir,kdir,idirSrc,jdirSrc,kdirSrc, &
                                 ibeg,iend,jbeg,jend,kbeg,kend, &
                                 ibegSrc,iendSrc,jbegSrc,jendSrc, &
                                 kbegSrc,kendSrc,mapMat )

  USE ModDataTypes
  USE ModError
  IMPLICIT NONE

! ... parameters
  INTEGER :: lb, lbs, l1SrcDir, l2SrcDir
  INTEGER :: idir, jdir, kdir, idirSrc, jdirSrc, kdirSrc
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: ibegSrc, iendSrc, jbegSrc, jendSrc, kbegSrc, kendSrc
  INTEGER :: mapMat(3,4)

  LOGICAL :: align

! ... local variables
  INTEGER :: l1, l2, lf, l1Src, l2Src, lfSrc, dir, dirSrc
  INTEGER :: ijkBeg(3), ijkBegSrc(3)

!******************************************************************************
! set matrix to zero

  mapMat(:,:) = 0

! l1/l2-direction on source patch

  IF      (lbs==1 .OR. lbs==2) THEN    ! l1Src = j, l2Src = k
    l1Src  = 2
    l2Src  = 3
    lfSrc  = 1
    dirSrc = idirSrc

                      ijkBegSrc(1) = ibegSrc
                      ijkBegSrc(2) = jbegSrc
    IF (l1SrcDir < 0) ijkBegSrc(2) = jendSrc
                      ijkBegSrc(3) = kbegSrc
    IF (l2SrcDir < 0) ijkBegSrc(3) = kendSrc
  ELSE IF (lbs==3 .OR. lbs==4) THEN    ! l1Src = k, l2Src = i
    l1Src  = 3
    l2Src  = 1
    lfSrc  = 2
    dirSrc = jdirSrc

                      ijkBegSrc(2) = jbegSrc
                      ijkBegSrc(3) = kbegSrc
    IF (l1SrcDir < 0) ijkBegSrc(3) = kendSrc
                      ijkBegSrc(1) = ibegSrc
    IF (l2SrcDir < 0) ijkBegSrc(1) = iendSrc
  ELSE IF (lbs==5 .OR. lbs==6) THEN    ! l1Src = i, l2Src = j
    l1Src  = 1
    l2Src  = 2
    lfSrc  = 3
    dirSrc = kdirSrc

                      ijkBegSrc(3) = kbegSrc
                      ijkBegSrc(1) = ibegSrc
    IF (l1SrcDir < 0) ijkBegSrc(1) = iendSrc
                      ijkBegSrc(2) = jbegSrc
    IF (l2SrcDir < 0) ijkBegSrc(2) = jendSrc
  ENDIF

! l1/l2-direction on current patch

  IF      (lb==1 .OR. lb==2) THEN      ! l1 = j, l2 = k
    IF (align) THEN
      l1 = 2
      l2 = 3
    ELSE
      l1 = 3
      l2 = 2
    ENDIF
    lf        = 1
    dir       = idir
    ijkBeg(1) = ibeg - idir
    ijkbeg(2) = jbeg
    ijkbeg(3) = kbeg
  ELSE IF (lb==3 .OR. lb==4) THEN      ! l1 = k, l2 = i
    IF (align) THEN
      l1 = 3
      l2 = 1
    ELSE
      l1 = 1
      l2 = 3
    ENDIF
    lf        = 2
    dir       = jdir
    ijkBeg(1) = ibeg
    ijkbeg(2) = jbeg - jdir
    ijkbeg(3) = kbeg
  ELSE IF (lb==5 .OR. lb==6) THEN      ! l1 = i, l2 = j
    IF (align) THEN
      l1 = 1
      l2 = 2
    ELSE
      l1 = 2
      l2 = 1
    ENDIF
    lf        = 3
    dir       = kdir
    ijkBeg(1) = ibeg
    ijkbeg(2) = jbeg
    ijkbeg(3) = kbeg - kdir
  ENDIF

! assemble matrix

  mapMat(l1Src,l1) = l1SrcDir
  mapMat(l2Src,l2) = l2SrcDir
  mapMat(lfSrc,lf) = -dir*dirSrc

! add offsets

  mapMat(1,4) = ijkBegSrc(1) - (mapMat(1,1)*ijkBeg(1)+mapMat(1,2)*ijkBeg(2)+ &
                                mapMat(1,3)*ijkBeg(3))
  mapMat(2,4) = ijkBegSrc(2) - (mapMat(2,1)*ijkBeg(1)+mapMat(2,2)*ijkBeg(2)+ &
                                mapMat(2,3)*ijkBeg(3))
  mapMat(3,4) = ijkBegSrc(3) - (mapMat(3,1)*ijkBeg(1)+mapMat(3,2)*ijkBeg(2)+ &
                                mapMat(3,3)*ijkBeg(3))

END SUBROUTINE RFLO_GetPatchMapping

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_GetPatchMapping.F90,v $
! Revision 1.3  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.3  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.2  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/03/18 21:56:39  jblazek
! Finished multiblock and MPI.
!
!******************************************************************************






