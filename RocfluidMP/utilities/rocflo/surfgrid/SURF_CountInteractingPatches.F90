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
! Purpose: count number of boundary patches interacting with GenX.
!
! Description: none.
!
! Input: regions = region dimensions and topology.
!
! Output: nInteract = no. of interacting boundary patches.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: SURF_CountInteractingPatches.F90,v 1.3 2008/12/06 08:44:52 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE CountInteractingPatches( regions,nInteract )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  INTEGER :: nInteract

  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'CountInteractingPatches', &
                         'SURF_CountInteractingPatches.F90' )

! loop over all regions

  nInteract = 0

  DO iReg=1,global%nRegions

! - loop over all patches of a region

    DO iPatch=1,regions(iReg)%nPatches
      patch => regions(iReg)%levels(1)%patches(iPatch)
      IF (patch%bcCoupled == BC_EXTERNAL) nInteract = nInteract + 1
    ENDDO   ! iPatch

  ENDDO     ! iReg

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE CountInteractingPatches

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SURF_CountInteractingPatches.F90,v $
! Revision 1.3  2008/12/06 08:44:52  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:02  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/03 02:47:00  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:49:09  wasistho
! lower to upper case
!
! Revision 1.3  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.2  2003/03/20 19:48:09  haselbac
! Corrected mistake in phased check-in
!
! Revision 1.1  2002/10/19 00:40:31  jblazek
! Added utility (rflosurf) to write out surface grids for GenX.
!
!******************************************************************************







