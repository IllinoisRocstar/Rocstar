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
! Purpose: calculate start and end indices of a patch, as well as
!          the direction normal to the boundary of the patch.
!
! Description: file contains the following subroutines:
!
!  - GetPatchIndices      = cell indices in terms of physical cells
!  - GetPatchIndicesNodes = nodes indices in terms of physical nodes
!  - GetPatchDirection    = coordinate direction into the dummy/physical
!                           cell layers.
!                        
!
! Input: region = current region
!        patch  = current patch
!        iLev   = current grid level
!
! Output: ibeg, iend = indices in i-direction
!         jbeg, jend = indices in j-direction
!         kbeg, kend = indices in k-direction
!         i/j/kdir   = direction (+-1 for coordinate normal to patch,
!                      =0 otherwise)
!
! Notes: ordering of l1, l2 patch coordinates is i, j, k (cyclic).
!        The sign of i/j/kdir is such that the vector (idir,jdir,kdir)
!        always points INTO the physical domain. Thus, the n-th dummy
!        cell is: {i-n*idir, j-n*jdir, k-n*kdir}. The m-th physical cell
!        is: {i+(m-1)*idir, j+(m-1)*jdir, k+(m-1)*kdir}, with i,j,k being
!        the indices within a patch (and 1st physical cell).
!
!******************************************************************************
!
! $Id: RFLO_GetPatchDimens.F90,v 1.3 2008/12/06 08:44:07 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend,jbeg,jend, &
                                 kbeg,kend )

  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys
  USE ModError
  IMPLICIT NONE

! ... parameters
  INTEGER        :: iLev, ibeg, iend, jbeg, jend, kbeg, kend
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... local variables
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_GetPatchIndices',&
  'RFLO_GetPatchDimens.F90' )

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )

  SELECT CASE (patch%lbound)
    CASE (1)
      ibeg = ipcbeg
      iend = ipcbeg
      jbeg = patch%l1beg
      jend = patch%l1end
      kbeg = patch%l2beg
      kend = patch%l2end
    CASE (2)
      ibeg = ipcend
      iend = ipcend
      jbeg = patch%l1beg
      jend = patch%l1end
      kbeg = patch%l2beg
      kend = patch%l2end
    CASE (3)
      ibeg = patch%l2beg
      iend = patch%l2end
      jbeg = jpcbeg
      jend = jpcbeg
      kbeg = patch%l1beg
      kend = patch%l1end
    CASE (4)
      ibeg = patch%l2beg
      iend = patch%l2end
      jbeg = jpcend
      jend = jpcend
      kbeg = patch%l1beg
      kend = patch%l1end
    CASE (5)
      ibeg = patch%l1beg
      iend = patch%l1end
      jbeg = patch%l2beg
      jend = patch%l2end
      kbeg = kpcbeg
      kend = kpcbeg
    CASE (6)
      ibeg = patch%l1beg
      iend = patch%l1end
      jbeg = patch%l2beg
      jend = patch%l2end
      kbeg = kpcend
      kend = kpcend
  END SELECT

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_GetPatchIndices

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_GetPatchIndicesNodes( region,patch,iLev,ibeg,iend, &
                                      jbeg,jend,kbeg,kend )

  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes
  USE ModError
  IMPLICIT NONE

! ... parameters
  INTEGER        :: iLev, ibeg, iend, jbeg, jend, kbeg, kend
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... local variables
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_GetPatchIndicesNodes',&
  'RFLO_GetPatchDimens.F90' )

  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )

  SELECT CASE (patch%lbound)
    CASE (1)
      ibeg = ipnbeg
      iend = ipnbeg
      jbeg = patch%l1beg
      jend = patch%l1end + 1
      kbeg = patch%l2beg
      kend = patch%l2end + 1
    CASE (2)
      ibeg = ipnend
      iend = ipnend
      jbeg = patch%l1beg
      jend = patch%l1end + 1
      kbeg = patch%l2beg
      kend = patch%l2end + 1
    CASE (3)
      ibeg = patch%l2beg
      iend = patch%l2end + 1
      jbeg = jpnbeg
      jend = jpnbeg
      kbeg = patch%l1beg
      kend = patch%l1end + 1
    CASE (4)
      ibeg = patch%l2beg
      iend = patch%l2end + 1
      jbeg = jpnend
      jend = jpnend
      kbeg = patch%l1beg
      kend = patch%l1end + 1
    CASE (5)
      ibeg = patch%l1beg
      iend = patch%l1end + 1
      jbeg = patch%l2beg
      jend = patch%l2end + 1
      kbeg = kpnbeg
      kend = kpnbeg
    CASE (6)
      ibeg = patch%l1beg
      iend = patch%l1end + 1
      jbeg = patch%l2beg
      jend = patch%l2end + 1
      kbeg = kpnend
      kend = kpnend
  END SELECT

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_GetPatchIndicesNodes

! #############################################################################
! #############################################################################

SUBROUTINE RFLO_GetPatchDirection( patch,idir,jdir,kdir )

  USE ModBndPatch, ONLY : t_patch
  IMPLICIT NONE

! ... parameters
  INTEGER       :: idir, jdir, kdir
  TYPE(t_patch) :: patch

!******************************************************************************

  SELECT CASE (patch%lbound)
    CASE (1)
      idir = +1
      jdir =  0
      kdir =  0
    CASE (2)
      idir = -1
      jdir =  0
      kdir =  0
    CASE (3)
      idir =  0
      jdir = +1
      kdir =  0
    CASE (4)
      idir =  0
      jdir = -1
      kdir =  0
    CASE (5)
      idir =  0
      jdir =  0
      kdir = +1
    CASE (6)
      idir =  0
      jdir =  0
      kdir = -1
  END SELECT

END SUBROUTINE RFLO_GetPatchDirection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_GetPatchDimens.F90,v $
! Revision 1.3  2008/12/06 08:44:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:20  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 21:25:16  wasistho
! lower to upper case
!
! Revision 1.8  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.4  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.3  2002/03/18 21:56:39  jblazek
! Finished multiblock and MPI.
!
! Revision 1.2  2002/02/21 23:25:04  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/08 22:09:16  jblazek
! Added calculation of face vectors and volumes.
!
!******************************************************************************








