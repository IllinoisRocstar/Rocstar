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
! Purpose: store outgoing data into GenX buffers before updating boundary
!          conditions (for now, the routine only sends the density at
!          the interface).
!
! Description: none.
!
! Input: region = dimensions of patches, types of BC`s, flow variables
!
! Output: regions%levels%patches = data in buffers for GenX.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_SendBoundaryValuesAlpha.F90,v 1.7 2008/12/06 08:44:01 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_SendBoundaryValuesAlpha( region )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetCellOffset, &
                            MixtPerf_R_M, MixtPerf_D_PRT
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iPatch, i, j, k, ijkC, i2d, n1, n2, ng1, ng2

! ... local variables
  INTEGER :: iLev, bcType, lbound, iCOff, ijCOff, nOff
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, distrib, indMol

  REAL(RFREAL)          :: tBurn, rgas
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), gv(:,:), vals(:,:)

  TYPE(t_global), POINTER :: global
  TYPE(t_patch) , POINTER :: patch

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_SendBoundaryValuesAlpha',&
  'RFLO_SendBoundaryValuesAlpha.F90' )

! store pointers to variables -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  indMol = region%levels(iLev)%mixt%indMol

  cv  => region%levels(iLev)%mixt%cv
  dv  => region%levels(iLev)%mixt%dv
  gv  => region%levels(iLev)%mixt%gv

! loop over all cells of the patch (if an interface)

  DO iPatch=1,region%nPatches

    patch   => region%levels(iLev)%patches(iPatch)
    bcType  =  patch%bcType
    lbound  =  patch%lbound
    distrib =  patch%mixt%distrib
    nOff    =  ABS(patch%l1end-patch%l1beg) + 1
    vals    => patch%mixt%vals

    IF ( patch%bcCoupled == BC_EXTERNAL .OR. &   ! data from outside
        (patch%bcCoupled == BC_INTERNAL .AND. &  ! data from internal APN
         bcType == BC_INJECTION_APN)) THEN        

      CALL RFLO_GetPatchIndices( region,patch,iLev, &
                                 ibeg,iend,jbeg,jend,kbeg,kend )

! --- outgoing data

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
            IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) THEN
              i2d   = distrib * IndIJ(n1,n2,nOff)
              tBurn = vals(BCDAT_INJECT_TEMP,i2d)

              IF (bcType==BC_INJECTION_APN) THEN
                patch%mdotAlp(    ng1,ng2) = vals(BCDAT_INJECT_MFRATE,i2d)
                patch%tflmAlp(    ng1,ng2) = vals(BCDAT_INJECT_TEMP  ,i2d)
                patch%rhofvfAlp(1,ng1,ng2) = vals(BCDAT_INJECT_RFVFU ,i2d) 
                patch%rhofvfAlp(2,ng1,ng2) = vals(BCDAT_INJECT_RFVFV ,i2d) 
                patch%rhofvfAlp(3,ng1,ng2) = vals(BCDAT_INJECT_RFVFW ,i2d) 
              ENDIF
              IF (patch%bFlag(ng1,ng2) == 1) THEN    ! burning
                rgas                   = MixtPerf_R_M( gv(GV_MIXT_MOL,ijkC*indMol) )
                patch%rhofAlp(ng1,ng2) = MixtPerf_D_PRT( dv(DV_MIXT_PRES,ijkC), &
                                                         rgas,tBurn )
              ELSE                                   ! not burning
                patch%rhofAlp(ng1,ng2) = cv(CV_MIXT_DENS,ijkC)
              ENDIF
            ELSE                           ! not an injection boundary
              patch%rhofAlp(ng1,ng2) = cv(CV_MIXT_DENS,ijkC)
            ENDIF
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k

    ENDIF  ! external BC
  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_SendBoundaryValuesAlpha

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_SendBoundaryValuesAlpha.F90,v $
! Revision 1.7  2008/12/06 08:44:01  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/08/19 15:37:52  mparmar
! Renamed patch variables
!
! Revision 1.4  2006/01/25 04:31:33  wasistho
! bug fixed allow injectionAPN to send
!
! Revision 1.3  2006/01/25 01:37:04  wasistho
! output tflame_alp and rhofVf_alp to Genx
!
! Revision 1.2  2006/01/24 07:05:07  wasistho
! send mdotAlp for InjectionAPN case
!
! Revision 1.1  2004/12/01 21:23:53  haselbac
! Initial revision after changing case
!
! Revision 1.7  2004/11/13 22:37:50  wasistho
! invert orientation of genx-surface-variables
!
! Revision 1.6  2003/11/20 16:40:34  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.3  2003/05/15 02:57:00  jblazek
! Inlined index function.
!
! Revision 1.2  2002/10/15 00:49:29  jblazek
! Got rid of second nOff ...
!
! Revision 1.1  2002/10/15 00:38:18  jblazek
! Added routine to send newest fluids density to GenX.
!
!******************************************************************************







