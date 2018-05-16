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
! Purpose: check if control volumes > 0 and if sum of face vectors = 0.
!
! Description: none.
!
! Input: iReg               = region number
!        region%levels%grid = dimensions, coordinates, face vectors
!                             (current region)
!
! Output: none.
!
! Notes: metrics is checked for physical cells only.
!
!******************************************************************************
!
! $Id: RFLO_CheckMetrics.F90,v 1.15 2009/10/17 01:17:16 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CheckMetrics( iReg,region )

  USE ModDataTypes
  USE ModGlobal,     ONLY : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset, &
                            RFLO_GetNodeOffset, RFLO_GetDimensDummy
  USE ModError
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: iReg

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iLev, i, j, k, l,myerr

! ... local variables
  CHARACTER(CHRLEN) :: msg

  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ijkCell, ijkNode(8), face(6)
  INTEGER :: im, jm, km

  REAL(RFREAL) :: sum(3), edgeLen, edgeMin2, edgeMin3, edge(3,4), dprod(8)
  REAL(RFREAL) :: ds, dprodm, dpmax, emag, smag1, smag2, skewness, minVol
  REAL(RFREAL), POINTER :: xyz(:,:), si(:,:), sj(:,:), sk(:,:), vol(:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global
  skewness = 0.0

  CALL RegisterFunction( global,'RFLO_CheckMetrics',&
  'RFLO_CheckMetrics.F90' )

! loop over all grid levels ---------------------------------------------------

  DO iLev=1,region%nGridLevels

    CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                             jpcbeg,jpcend,kpcbeg,kpcend )
    CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
    CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

    xyz => region%levels(iLev)%grid%xyz
    si  => region%levels(iLev)%grid%si
    sj  => region%levels(iLev)%grid%sj
    sk  => region%levels(iLev)%grid%sk
    vol => region%levels(iLev)%grid%vol

! - calculate the shortest distance

    edgeLen = 1.E+30_RFREAL

!    call MPI_BARRIER(MPI_COMM_WORLD,myerr)
!      if(global%myprocid .eq. 0) then
!       write(*,*) 'calling shortest distance'
!    endif

    DO k=kpcbeg,kpcend
      DO j=jpcbeg,jpcend
        DO i=ipcbeg,ipcend
          ijkNode(1) = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
          ijkNode(2) = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
          ijkNode(3) = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
          ijkNode(4) = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
          ds      = SQRT((xyz(XCOORD,ijkNode(2))-xyz(XCOORD,ijkNode(1)))**2+ &
                         (xyz(YCOORD,ijkNode(2))-xyz(YCOORD,ijkNode(1)))**2+ &
                         (xyz(ZCOORD,ijkNode(2))-xyz(ZCOORD,ijkNode(1)))**2)
          edgeLen = MIN(edgeLen,ds)
          ds      = SQRT((xyz(XCOORD,ijkNode(3))-xyz(XCOORD,ijkNode(1)))**2+ &
                         (xyz(YCOORD,ijkNode(3))-xyz(YCOORD,ijkNode(1)))**2+ &
                         (xyz(ZCOORD,ijkNode(3))-xyz(ZCOORD,ijkNode(1)))**2)
          edgeLen = MIN(edgeLen,ds)
          ds      = SQRT((xyz(XCOORD,ijkNode(4))-xyz(XCOORD,ijkNode(1)))**2+ &
                         (xyz(YCOORD,ijkNode(4))-xyz(YCOORD,ijkNode(1)))**2+ &
                         (xyz(ZCOORD,ijkNode(4))-xyz(ZCOORD,ijkNode(1)))**2)
          edgeLen = MIN(edgeLen,ds)
        ENDDO
      ENDDO
    ENDDO

!    call MPI_BARRIER(MPI_COMM_WORLD,myerr)
!      if(global%myprocid .eq. 0) then
!       write(*,*) 'doing closedness'
!    endif

    edgeMin2 = 0.1_RFREAL*edgeLen*edgeLen
    edgeMin3 = edgeLen*edgeMin2

! - check face vectors: closedness

    DO k=kpcbeg,kpcend
      DO j=jpcbeg,jpcend
        DO i=ipcbeg,ipcend
          face(1) = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
          face(2) = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
          face(4) = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
          face(6) = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)

          sum(1:3) = si(1:3,face(1)) - si(1:3,face(2)) + &
                     sj(1:3,face(1)) - sj(1:3,face(4)) + &
                     sk(1:3,face(1)) - sk(1:3,face(6))

          IF (MAXVAL(ABS(sum)) > edgeMin2) THEN
            WRITE(msg,1000) MAXVAL(ABS(sum)),iReg,iLev,i,j,k
            CALL ErrorStop( global,ERR_FACEVEC_SUM,__LINE__,msg )
          ENDIF
        ENDDO
      ENDDO
    ENDDO

!    call MPI_BARRIER(MPI_COMM_WORLD,myerr)
!!      if(global%myprocid .eq. 0) then
!       write(*,*) 'doing face vectors'
!    endif

! - check face vectors: inverted cells

    region%levels(iLev)%grid%remesh = 0
    dprodm =  1._RFREAL
    dpmax  = -1._RFREAL

    DO k=kpcbeg,kpcend
      DO j=jpcbeg,jpcend
        DO i=ipcbeg,ipcend
          ijkNode(1) = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
          ijkNode(2) = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
          ijkNode(3) = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
          ijkNode(4) = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
          ijkNode(5) = IndIJK(i+1,j+1,k+1,iNOff,ijNOff)
          ijkNode(6) = IndIJK(i+1,j+1,k  ,iNOff,ijNOff)
          ijkNode(7) = IndIJK(i+1,j  ,k+1,iNOff,ijNOff)
          ijkNode(8) = IndIJK(i  ,j+1,k+1,iNOff,ijNOff)

          face(1) = ijkNode(1)
          face(2) = ijkNode(2)
          face(4) = ijkNode(3)
          face(6) = ijkNode(4)

! ------- check for inverted i-faces

          edge(1,1) = xyz(XCOORD,ijkNode(1))-xyz(XCOORD,ijkNode(2))
          edge(2,1) = xyz(YCOORD,ijkNode(1))-xyz(YCOORD,ijkNode(2))
          edge(3,1) = xyz(ZCOORD,ijkNode(1))-xyz(ZCOORD,ijkNode(2))

          edge(1,2) = xyz(XCOORD,ijkNode(3))-xyz(XCOORD,ijkNode(6))
          edge(2,2) = xyz(YCOORD,ijkNode(3))-xyz(YCOORD,ijkNode(6))
          edge(3,2) = xyz(ZCOORD,ijkNode(3))-xyz(ZCOORD,ijkNode(6))

          edge(1,3) = xyz(XCOORD,ijkNode(8))-xyz(XCOORD,ijkNode(5))
          edge(2,3) = xyz(YCOORD,ijkNode(8))-xyz(YCOORD,ijkNode(5))
          edge(3,3) = xyz(ZCOORD,ijkNode(8))-xyz(ZCOORD,ijkNode(5))

          edge(1,4) = xyz(XCOORD,ijkNode(4))-xyz(XCOORD,ijkNode(7))
          edge(2,4) = xyz(YCOORD,ijkNode(4))-xyz(YCOORD,ijkNode(7))
          edge(3,4) = xyz(ZCOORD,ijkNode(4))-xyz(ZCOORD,ijkNode(7))

          DO l=1,4
            emag  = SQRT( edge(1,l)*edge(1,l) + &
                          edge(2,l)*edge(2,l) + &
                          edge(3,l)*edge(3,l) )
            smag1 = SQRT( si(1,face(1))*si(1,face(1)) + &
                          si(2,face(1))*si(2,face(1)) + &
                          si(3,face(1))*si(3,face(1)) )
            smag2 = SQRT( si(1,face(2))*si(1,face(2)) + &
                          si(2,face(2))*si(2,face(2)) + &
                          si(3,face(2))*si(3,face(2)) )

            dprod(l)  =  (si(1,face(1))*edge(1,l) + &
                          si(2,face(1))*edge(2,l) + &
                          si(3,face(1))*edge(3,l))/(smag1*emag)

            dprod(l+4) = (si(1,face(2))*edge(1,l) + &
                          si(2,face(2))*edge(2,l) + &
                          si(3,face(2))*edge(3,l))/(smag2*emag)
          ENDDO

          IF (MINVAL(dprod) < dprodm) THEN
            dprodm = MINVAL( dprod )
            im = i
            jm = j
            km = k
          ENDIF

          IF (MAXVAL(dprod) > dpmax) THEN
            dpmax = MAXVAL( dprod )
          ENDIF

! ------- check for inverted j-faces

          edge(1,1) = xyz(XCOORD,ijkNode(1))-xyz(XCOORD,ijkNode(3))
          edge(2,1) = xyz(YCOORD,ijkNode(1))-xyz(YCOORD,ijkNode(3))
          edge(3,1) = xyz(ZCOORD,ijkNode(1))-xyz(ZCOORD,ijkNode(3))

          edge(1,2) = xyz(XCOORD,ijkNode(4))-xyz(XCOORD,ijkNode(8))
          edge(2,2) = xyz(YCOORD,ijkNode(4))-xyz(YCOORD,ijkNode(8))
          edge(3,2) = xyz(ZCOORD,ijkNode(4))-xyz(ZCOORD,ijkNode(8))

          edge(1,3) = xyz(XCOORD,ijkNode(7))-xyz(XCOORD,ijkNode(5))
          edge(2,3) = xyz(YCOORD,ijkNode(7))-xyz(YCOORD,ijkNode(5))
          edge(3,3) = xyz(ZCOORD,ijkNode(7))-xyz(ZCOORD,ijkNode(5))

          edge(1,4) = xyz(XCOORD,ijkNode(2))-xyz(XCOORD,ijkNode(6))
          edge(2,4) = xyz(YCOORD,ijkNode(2))-xyz(YCOORD,ijkNode(6))
          edge(3,4) = xyz(ZCOORD,ijkNode(2))-xyz(ZCOORD,ijkNode(6))

          DO l=1,4
            emag  = SQRT( edge(1,l)*edge(1,l) + &
                          edge(2,l)*edge(2,l) + &
                          edge(3,l)*edge(3,l) )
            smag1 = SQRT( sj(1,face(1))*sj(1,face(1)) + &
                          sj(2,face(1))*sj(2,face(1)) + &
                          sj(3,face(1))*sj(3,face(1)) )
            smag2 = SQRT( sj(1,face(4))*sj(1,face(4)) + &
                          sj(2,face(4))*sj(2,face(4)) + &
                          sj(3,face(4))*sj(3,face(4)) )

            dprod(l)  =  (sj(1,face(1))*edge(1,l) + &
                          sj(2,face(1))*edge(2,l) + &
                          sj(3,face(1))*edge(3,l))/(smag1*emag)

            dprod(l+4) = (sj(1,face(4))*edge(1,l) + &
                          sj(2,face(4))*edge(2,l) + &
                          sj(3,face(4))*edge(3,l))/(smag2*emag)
          ENDDO

          IF (MINVAL(dprod) < dprodm) THEN
            dprodm = MINVAL( dprod )
            im = i
            jm = j
            km = k
          ENDIF

          IF (MAXVAL(dprod) > dpmax) THEN
            dpmax = MAXVAL( dprod )
          ENDIF

! ------- check for inverted k-faces

          edge(1,1) = xyz(XCOORD,ijkNode(1))-xyz(XCOORD,ijkNode(4))
          edge(2,1) = xyz(YCOORD,ijkNode(1))-xyz(YCOORD,ijkNode(4))
          edge(3,1) = xyz(ZCOORD,ijkNode(1))-xyz(ZCOORD,ijkNode(4))

          edge(1,2) = xyz(XCOORD,ijkNode(2))-xyz(XCOORD,ijkNode(7))
          edge(2,2) = xyz(YCOORD,ijkNode(2))-xyz(YCOORD,ijkNode(7))
          edge(3,2) = xyz(ZCOORD,ijkNode(2))-xyz(ZCOORD,ijkNode(7))

          edge(1,3) = xyz(XCOORD,ijkNode(6))-xyz(XCOORD,ijkNode(5))
          edge(2,3) = xyz(YCOORD,ijkNode(6))-xyz(YCOORD,ijkNode(5))
          edge(3,3) = xyz(ZCOORD,ijkNode(6))-xyz(ZCOORD,ijkNode(5))

          edge(1,4) = xyz(XCOORD,ijkNode(3))-xyz(XCOORD,ijkNode(8))
          edge(2,4) = xyz(YCOORD,ijkNode(3))-xyz(YCOORD,ijkNode(8))
          edge(3,4) = xyz(ZCOORD,ijkNode(3))-xyz(ZCOORD,ijkNode(8))

          DO l=1,4
            emag  = SQRT( edge(1,l)*edge(1,l) + &
                          edge(2,l)*edge(2,l) + &
                          edge(3,l)*edge(3,l) )
            smag1 = SQRT( sk(1,face(1))*sk(1,face(1)) + &
                          sk(2,face(1))*sk(2,face(1)) + &
                          sk(3,face(1))*sk(3,face(1)) )
            smag2 = SQRT( sk(1,face(6))*sk(1,face(6)) + &
                          sk(2,face(6))*sk(2,face(6)) + &
                          sk(3,face(6))*sk(3,face(6)) )

            dprod(l)  =  (sk(1,face(1))*edge(1,l) + &
                          sk(2,face(1))*edge(2,l) + &
                          sk(3,face(1))*edge(3,l))/(smag1*emag)

            dprod(l+4) = (sk(1,face(6))*edge(1,l) + &
                          sk(2,face(6))*edge(2,l) + &
                          sk(3,face(6))*edge(3,l))/(smag2*emag)
          ENDDO

          IF (MINVAL(dprod) < dprodm) THEN
            dprodm = MINVAL( dprod )
            im = i
            jm = j
            km = k
          ENDIF

          IF (MAXVAL(dprod) > dpmax) THEN
            dpmax = MAXVAL( dprod )
          ENDIF

        ENDDO
      ENDDO
    ENDDO

!    call MPI_BARRIER(MPI_COMM_WORLD,myerr)
!      if(global%myprocid .eq. 0) then
!       write(*,*) 'done doing face vectors'
!    endif
!    if(dpmax .eq. 0) then
!       write(*,*) 'dpmax == 0'
!       stop
!    endif
!    skewness = dprodm/dpmax
!    region%levels(iLev)%grid%skewness = skewness  ! block skewness
!    global%skewness                   = skewness  ! for global skewness

!    IF (skewness < 1.E-3_RFREAL) THEN     ! positive criterion
    IF (skewness < 0._RFREAL) THEN         ! negative criterion
!      WRITE(msg,2000) skewness, &
!                      vol(IndIJK(im,jm,km,iCOff,ijCOff)), &
!                      iReg,iLev,im,jm,km
      WRITE(STDOUT,'(A)') SOLVER_NAME//'Warning: inverted cell, '//msg
!      region%levels(iLev)%grid%remesh = 1
    ENDIF

! - check volumes
!    call MPI_BARRIER(MPI_COMM_WORLD,myerr)
!    write(*,*) 'checking volume'
    minVol = HUGE( 1._RFREAL )
!    DO k=kdcbeg,kdcend
!      DO j=jdcbeg,jdcend
!        DO i=idcbeg,idcend
    DO k=kpcbeg,kpcend
      DO j=jpcbeg,jpcend
        DO i=ipcbeg,ipcend
          ijkCell = IndIJK(i,j,k,iCOff,ijCOff)
          minVol  = MIN( minVol,vol(ijkCell) )
          
!          IF (vol(ijkCell) < edgeMin3) THEN   ! positive criterion
          IF (vol(ijkCell) < 0._RFREAL) THEN   ! negative criterion
            WRITE(msg,1000) vol(ijkCell),iReg,iLev,i,j,k
            CALL ErrorStop( global,ERR_VOLUME_SIZE,__LINE__,msg )
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    region%levels(iLev)%grid%minVol = minVol  ! block minimum cell-vol
    global%minVol                   = minVol  ! for global minimum cell-vol
!    call MPI_BARRIER(MPI_COMM_WORLD,myerr)
!    write(*,*) 'checking volume done'

  ENDDO   ! iLev

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1000 FORMAT('value = ',1PE12.5,', region ',I5,', lev ',I1,', ijk=',3I4)
2000 FORMAT('skewness, vol.= ',2(1PE10.3),', reg.',I5,', lev.',I1,', ijk=',3I4)

END SUBROUTINE RFLO_CheckMetrics

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CheckMetrics.F90,v $
! Revision 1.15  2009/10/17 01:17:16  mtcampbe
!
! Port to Mercury, George's skewness fix, MP NATIVE IO fix
!
! Revision 1.14  2009/08/12 04:15:58  mtcampbe
! Major update, bugfix from Abe development, more propagation compatibility,
! some Rocstar IO changes, Ju's temporary clipping fix for turbulence. A bug
! fix for initialization IO.
!
! Revision 1.13  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.12  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.11  2006/03/18 11:05:20  wasistho
! computed minvol
!
! Revision 1.10  2006/03/15 06:42:10  wasistho
! added region and global skewness
!
! Revision 1.9  2006/03/08 09:14:22  wasistho
! changed skewness report
!
! Revision 1.8  2006/03/04 04:27:54  wasistho
! normalized cell skewness parameter
!
! Revision 1.7  2005/05/28 05:41:59  wasistho
! switch positive to negative criteria for inverted cells and cell-volume
!
! Revision 1.6  2005/05/27 08:40:16  wasistho
! warning inverted cell after solver_name
!
! Revision 1.5  2005/05/27 01:51:02  wasistho
! added remeshing
!
! Revision 1.4  2005/05/06 07:36:07  wasistho
! modified warning for inverted cells
!
! Revision 1.3  2005/04/25 23:40:35  wasistho
! modified inverted cell criterion
!
! Revision 1.2  2005/04/25 04:59:26  wasistho
! added inverted cell detector
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.13  2004/08/04 21:33:27  wasistho
! insert missing ABS in face-vectors check
!
! Revision 1.12  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.8  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.7  2002/12/06 22:29:26  jblazek
! Corrected bug for geometry exchange between minimal patches.
!
! Revision 1.6  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.5  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/02/01 22:17:38  jblazek
! Change addressing of face vectors at block boundaries.
!
! Revision 1.3  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.2  2002/01/08 22:09:17  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.1  2001/12/19 23:09:21  jblazek
! Added routines to read grid and solution.
!
!******************************************************************************







