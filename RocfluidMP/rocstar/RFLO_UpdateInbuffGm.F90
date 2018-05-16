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
! Purpose: dummy routine to fill the geometry buffer.
!
! Description: none.
!
! Input: globalGenx = pointer to global data
!        dAlpha     = relative time step.
!
! Output: nothing useful.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_UpdateInbuffGm.F90,v 1.3 2008/12/06 08:44:01 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_UpdateInbuffGm( globalGenx,dAlpha )

  USE ModDataTypes
  USE ModRocstar, ONLY       : t_globalGenx
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  DOUBLE PRECISION, INTENT(in) :: dAlpha

  TYPE(t_globalGenx), POINTER :: globalGenx

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  INTEGER :: iLev

  TYPE(t_global), POINTER  :: global
  TYPE(t_patch) , POINTER  :: patch
  TYPE (t_region), POINTER :: regions(:)

!******************************************************************************

  global  => globalGenx%global
  regions => globalGenx%regions

  CALL RegisterFunction( global,'RFLO_UpdateInbuffGm',&
  'RFLO_UpdateInbuffGm.F90' )

! fill geometry buffer

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
      iLev = regions(iReg)%currLevel
      DO iPatch=1,regions(iReg)%nPatches
        patch => regions(iReg)%levels(iLev)%patches(iPatch)
        IF (patch%bcCoupled == BC_EXTERNAL) THEN        ! data from outside
          patch%duAlp(:,:,:) = 0._RFREAL
        ENDIF
      ENDDO    ! iPatch
    ENDIF      ! region on this processor and active
  ENDDO        ! iReg

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_UpdateInbuffGm

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_UpdateInbuffGm.F90,v $
! Revision 1.3  2008/12/06 08:44:01  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:23:55  haselbac
! Initial revision after changing case
!
! Revision 1.4  2003/05/15 02:57:00  jblazek
! Inlined index function.
!
! Revision 1.3  2003/01/28 21:35:33  jblazek
! Corrected the loop over patches.
!
! Revision 1.2  2002/10/30 22:07:51  jiao
! Changed to initialize duAlp to 0 instead of to stop the code.
!
! Revision 1.1  2002/09/20 22:22:34  jblazek
! Finalized integration into GenX.
!
!******************************************************************************







