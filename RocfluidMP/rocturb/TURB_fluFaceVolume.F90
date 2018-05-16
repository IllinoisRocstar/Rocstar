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
! Purpose: Define face volume.
!
! Description: The face volume is obtained by averaging the two adjacent cell 
!              volumes
!
! Input: region  = data of current region
!
! Notes: Face volume is used in several routines, a.o LesMij & LesCalcEddyVis.
!
!******************************************************************************
!
! $Id: TURB_fluFaceVolume.F90,v 1.6 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FluFaceVolume( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iN, iPatch, ifl

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch

  INTEGER :: ijkC0, ijkC1, ifg, ifgBeg
  INTEGER, POINTER      :: f2c(:,:)
  REAL(RFREAL), POINTER :: fvol(:), bfVol(:), vol(:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_fluFaceVolume.F90,v $ $Revision: 1.6 $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FluFaceVolume',&
  'TURB_fluFaceVolume.F90' )

! get indices and pointers ---------------------------------------------------

  f2c   => region%grid%f2c
  vol   => region%grid%vol
  fvol  => region%turb%fvolI
  bfVol => region%turb%bfVolI

! perform 2-point averaging of volume from centers to faces

  DO iN = 1,region%grid%nFaces
    ijkC0 = f2c(1,iN)
    ijkC1 = f2c(2,iN)

    fvol(iN) = 0.5_RFREAL*(vol(ijkC1)+vol(ijkC0))
  ENDDO  ! iC

  DO iPatch = 1,region%grid%nPatches
    patch  => region%patches(iPatch)
! TEMPORARY : removing usage of bf2bg from everywhere
!    ifgBeg =  patch%bf2bg(BF2BG_BEG)

    DO ifl = 1,patch%nBFaces
      ijkC0 = patch%bf2c(ifl)
      ifg   = ifl + ifgBeg-1
      bfVol(ifg) = vol(ijkC0)
    ENDDO ! ifl
  ENDDO   ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FluFaceVolume

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_fluFaceVolume.F90,v $
! Revision 1.6  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/08/19 15:41:02  mparmar
! Removed bf2bg
!
! Revision 1.3  2005/12/29 19:48:58  wasistho
! added boundary/patch treatment
!
! Revision 1.2  2004/03/27 02:16:42  wasistho
! compiled with Rocflu
!
! Revision 1.1  2004/03/19 02:54:58  wasistho
! prepared for RFLU
!
!
!******************************************************************************







