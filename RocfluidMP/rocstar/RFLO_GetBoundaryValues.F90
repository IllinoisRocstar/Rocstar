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
! Purpose: obtain boundary values from GenX.
!
! Description: none.
!
! Input: region = dimensions and topology.
!
! Output: region%levels%patch = input values for boundary conditions.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_GetBoundaryValues.F90,v 1.7 2008/12/06 08:44:00 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_GetBoundaryValues( region )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iPatch, i, j, k, ijkC, n1, n2, ng1, ng2

! ... local variables
  INTEGER :: iLev, bcType, lbound, iCOff, ijCOff, i2d, nOff
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend

  REAL(RFREAL), POINTER :: vals(:,:)

  TYPE(t_global), POINTER :: global
  TYPE(t_patch) , POINTER :: patch

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_GetBoundaryValues',&
  'RFLO_GetBoundaryValues.F90' )

! get offsets

  iLev = region%currLevel

  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

! loop over all cells of the patch (if an interface)

  DO iPatch=1,region%nPatches

    patch  => region%levels(iLev)%patches(iPatch)
    vals   => patch%mixt%vals
    bcType =  patch%bcType
    lbound =  patch%lbound
    nOff   =  ABS(patch%l1end-patch%l1beg) + 1

    IF (patch%bcCoupled == BC_EXTERNAL) THEN        ! data from outside

      CALL RFLO_GetPatchIndices( region,patch,iLev, &
                                 ibeg,iend,jbeg,jend,kbeg,kend )

      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkC = IndIJK(i,j,k,iCOff,ijCOff)
            IF      (lbound==1 .OR. lbound==2) THEN
              n1 = j - jbeg
              n2 = k - kbeg
              IF (lbound == 2) THEN
                ng1 = j - jbeg + 1
              ELSE
                ng1 = jend - j + 1
              ENDIF
              ng2 = k - kbeg + 1
            ELSE IF (lbound==3 .OR. lbound==4) THEN
              n1 = k - kbeg
              n2 = i - ibeg
              IF (lbound == 4) THEN
                ng2 = i - ibeg + 1
              ELSE
                ng2 = iend - i + 1
              ENDIF
              ng1 = k - kbeg + 1
            ELSE IF (lbound==5 .OR. lbound==6) THEN
              n1 = i - ibeg
              n2 = j - jbeg
              IF (lbound == 6) THEN
                ng1 = i - ibeg + 1
              ELSE
                ng1 = iend - i + 1
              ENDIF
              ng2 = j - jbeg + 1
            ENDIF
            i2d = IndIJ(n1,n2,nOff)
            IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) THEN
              IF (patch%bFlag(ng1,ng2) == 1) THEN ! burning
                vals(BCDAT_INJECT_MFRATE,i2d) = patch%mdotAlp(ng1,ng2)
              ELSE                                  ! not ignited
                vals(BCDAT_INJECT_MFRATE,i2d) = -1._RFREAL
              ENDIF
              vals(BCDAT_INJECT_TEMP ,i2d) = patch%tflmAlp  (  ng1,ng2)
              vals(BCDAT_INJECT_RFVFU,i2d) = patch%rhofvfAlp(1,ng1,ng2)
              vals(BCDAT_INJECT_RFVFV,i2d) = patch%rhofvfAlp(2,ng1,ng2)
              vals(BCDAT_INJECT_RFVFW,i2d) = patch%rhofvfAlp(3,ng1,ng2)
            ELSE IF (bcType>=BC_NOSLIPWALL .AND. &
                     bcType<=BC_NOSLIPWALL+BC_RANGE) THEN
              ! store Tflm as wall temperature ???
              ! wall movement deduced from grid deformation
              ! This block reserved for heat conduction later,
              ! for now set alpha-variables to zero
              patch%mdotAlp(:,:)     = 0._RFREAL
              patch%tflmAlp(:,:)     = 0._RFREAL
              patch%rhofvfAlp(:,:,:) = 0._RFREAL
            ELSE
              patch%mdotAlp(:,:)     = 0._RFREAL
              patch%tflmAlp(:,:)     = 0._RFREAL
              patch%rhofvfAlp(:,:,:) = 0._RFREAL
            ENDIF
          ENDDO
        ENDDO
      ENDDO

    ENDIF  ! external BC
  ENDDO    ! iPatch

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_GetBoundaryValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_GetBoundaryValues.F90,v $
! Revision 1.7  2008/12/06 08:44:00  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/08/19 15:37:48  mparmar
! Renamed patch variables
!
! Revision 1.4  2006/03/03 07:11:09  wasistho
! zero out incoming varsAlp for non-injection
!
! Revision 1.3  2005/02/01 21:18:09  wasistho
! reactivated bflag
!
! Revision 1.2  2005/01/21 00:21:26  wasistho
! temporary commented bflag condition
!
! Revision 1.1  2004/12/01 21:23:49  haselbac
! Initial revision after changing case
!
! Revision 1.9  2004/11/13 22:37:35  wasistho
! invert orientation of genx-surface-variables
!
! Revision 1.8  2003/11/20 16:40:33  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/05/15 02:57:00  jblazek
! Inlined index function.
!
! Revision 1.4  2002/10/12 20:14:35  jblazek
! Corrected bug in indices of bFlag.
!
! Revision 1.3  2002/10/10 23:49:48  jblazek
! Changed orientation of surface grid.
!
! Revision 1.2  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.1  2002/09/20 22:22:34  jblazek
! Finalized integration into GenX.
!
!******************************************************************************







