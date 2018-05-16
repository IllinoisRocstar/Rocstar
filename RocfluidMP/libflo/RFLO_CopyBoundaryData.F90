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
! Purpose: copy data associated with a patch from finer to coarser grid.
!
! Description: none.
!
! Input: patchPrev = patch on previous grid level.
!
! Output: patch = patch on current grid level.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_CopyBoundaryData.F90,v 1.4 2008/12/06 08:44:06 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CopyBoundaryData( global,patchPrev,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_InterpolDistrib
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_patch), POINTER  :: patchPrev, patch
  TYPE(t_global), POINTER :: global

! ... local variables
  INTEGER :: n1, n2, n1p, n2p, iOff, ijBeg, ijEnd, errorFlag

!******************************************************************************

  CALL RegisterFunction( global,'RFLO_CopyBoundaryData',&
  'RFLO_CopyBoundaryData.F90' )

! copy switches

  IF (patch%mixt%nSwitches > 0) THEN
    ALLOCATE( patch%mixt%switches(patch%mixt%nSwitches),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
    __LINE__ )
    patch%mixt%switches(:) = patchPrev%mixt%switches(:)
  ENDIF

! copy distributions

  IF (patch%mixt%distrib == BCDAT_DISTRIB) THEN
    n1p   = ABS(patchPrev%l1end-patchPrev%l1beg)
    n2p   = ABS(patchPrev%l2end-patchPrev%l2beg)
    n1    = ABS(patch%l1end-patch%l1beg)
    n2    = ABS(patch%l2end-patch%l2beg)
    iOff  = n1 + 1
    ijBeg = IndIJ( 0, 0,iOff)
    ijEnd = IndIJ(n1,n2,iOff)
    IF (patch%mixt%nData > 0) THEN
      ALLOCATE( patch%mixt%vals(patch%mixt%nData,ijBeg:ijEnd), &
                stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
      __LINE__ )
      CALL RFLO_InterpolDistrib( n1p,n2p,n1,n2,patch%mixt%nData, &
                                 patchPrev%mixt%vals,patch%mixt%vals )
    ENDIF

  ELSE
    IF (patch%mixt%nData > 0) THEN
      ALLOCATE( patch%mixt%vals(patch%mixt%nData,0:1),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,&
      __LINE__ )
      patch%mixt%vals(:,:) = patchPrev%mixt%vals(:,:)
    ENDIF
  ENDIF

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_CopyBoundaryData

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CopyBoundaryData.F90,v $
! Revision 1.4  2008/12/06 08:44:06  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:37:59  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.6  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.3  2003/05/15 02:57:01  jblazek
! Inlined index function.
!
! Revision 1.2  2002/10/25 18:36:19  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.1  2002/10/19 00:40:30  jblazek
! Added utility (rflosurf) to write out surface grids for GenX.
!
! Revision 1.2  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.1  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.5  2002/02/27 18:38:19  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
! Revision 1.4  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.3  2002/01/11 17:18:30  jblazek
! Updated description of I/O variables.
!
! Revision 1.2  2001/12/22 00:09:38  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1  2001/12/07 16:54:31  jblazek
! Added files to read user input.
!
!******************************************************************************







