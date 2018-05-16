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
! Purpose: calculation of pointers to grid cells/nodes or of i,j,k indices
!          from a pointer.
!
! Description: - GetIJK    = i,j,k indices from a pointer
!              - IndIJKMap = pointer to grid cell/node based on i,j,k indices
!                            of a patch and mapping matrix to source patch
!
! Notes: subroutines RFLO_GetCellOffset or RFLO_GetNodeOffset must be called
!        BEFORE the functions are used in order to set the correct offsets.
!        Further, the subroutine RFLO_GetPatchMapping has to be called before
!        the function IndIJKMap can be invoked, in order to set the mapping
!        matrix.
!
!******************************************************************************
!
! $Id: ModIndexing.F90,v 1.4 2008/12/06 08:44:14 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModIndexing

  IMPLICIT NONE

  CONTAINS

! from a pointer to i, j, k ---------------------------------------------------

  SUBROUTINE GetIJK( ijk,iOffset,ijOffset,nDumCells,i,j,k )
    INTEGER :: ijk, iOffset, ijOffset, nDumCells, i, j, k
    INTEGER :: ijkmin, q

    ijkmin = 1 - nDumCells*(1+iOffset+ijOffset)
    q      = (ijk-ijkmin)/iOffset
    i      = MOD(ijk-ijkmin,iOffset) - nDumCells + 1
    j      = MOD(q,ijOffset/iOffset) - nDumCells + 1
    k      = (ijk-ijkmin)/ijOffset   - nDumCells + 1
  END SUBROUTINE GetIJK

! pointer to cell/node using patch mapping ------------------------------------

  INTEGER FUNCTION IndIJKMap( i,j,k,mapMat,iOffset,ijOffset )
    INTEGER :: i, j, k, mapMat(3,4), iOffset, ijOffset
    INTEGER :: is, js, ks

    is = i*mapMat(1,1) + j*mapMat(1,2) + k*mapMat(1,3) + mapMat(1,4)
    js = i*mapMat(2,1) + j*mapMat(2,2) + k*mapMat(2,3) + mapMat(2,4)
    ks = i*mapMat(3,1) + j*mapMat(3,2) + k*mapMat(3,3) + mapMat(3,4)

    IndIJKMap = is + (js-1)*iOffset + (ks-1)*ijOffset
  END FUNCTION IndIJKMap

END MODULE ModIndexing

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModIndexing.F90,v $
! Revision 1.4  2008/12/06 08:44:14  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.1  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.5  2002/03/18 22:46:07  jblazek
! Finished multiblock and MPI.
!
! Revision 1.4  2002/01/12 00:02:48  jblazek
! Added postprocessor.
!
! Revision 1.3  2001/12/11 21:59:28  jblazek
! memory allocation added.
!
! Revision 1.2  2001/12/08 00:18:41  jblazek
! Added routines to read BC input file.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
!******************************************************************************






