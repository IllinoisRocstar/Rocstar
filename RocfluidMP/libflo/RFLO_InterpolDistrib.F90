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
! Purpose: interpolate distribution of boundary values from finer to coarser
!          grid.
!
! Description: interpolation is based on simple arithmetic averaging. Area
!              weighted averaging would be more accurate, but ...
!
! Input: n1f, n2f = no. of values on fine grid in 1st and 2nd direction
!        n1c, n2c = no. of values on coarse grid in 1st and 2nd direction
!        nData    = number of equations per node/cell
!        valf     = values on fine grid
!
! Output: valc = values on coarse grid
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_InterpolDistrib.F90,v 1.3 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_InterpolDistrib( n1f,n2f,n1c,n2c,nData,valf,valc )

  USE ModDataTypes
  USE ModError
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: n1f, n2f, n1c, n2c, nData

  REAL(RFREAL), POINTER :: valf(:,:), valc(:,:)

! ... loop variables
  INTEGER :: i, j, n

! ... local variables
  INTEGER :: icOff, ifOff, ii, jj, ijcp

!******************************************************************************
! loop over no. of data and patch` dimensions

  icOff = n1c + 1
  ifOff = n1f + 1

  DO n=1,nData
    DO j=0,n2c
      jj = 2*j
      DO i=0,n1c
        ii   = 2*i
        ijcp = IndIJ(i,j,icOff)
        valc(n,ijcp) = 0.25_RFREAL*(valf(n,IndIJ(ii  ,jj  ,ifOff))+ &
                                    valf(n,IndIJ(ii+1,jj  ,ifOff))+ &
                                    valf(n,IndIJ(ii  ,jj+1,ifOff))+ &
                                    valf(n,IndIJ(ii+1,jj+1,ifOff)))
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE RFLO_InterpolDistrib

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_InterpolDistrib.F90,v $
! Revision 1.3  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.6  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.5  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/03/30 00:50:48  jblazek
! Cleaned up with flint.
!
! Revision 1.3  2002/01/11 17:13:30  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.2  2001/12/22 00:09:36  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1  2001/12/08 00:18:40  jblazek
! Added routines to read BC input file.
!
!******************************************************************************






