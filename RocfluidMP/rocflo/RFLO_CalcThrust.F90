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
! Purpose: calculate amount of thrust for a region.
!
! Description: none.
!
! Input: region%levels%mixt        = flow variables
!        region%levels%grid%si/j/k = face vectors (at boundaries)
!
! Output: global%thrustTotal = thrust
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_CalcThrust.F90,v 1.3 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_CalcThrust( region )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetCellOffset, RFLO_GetNodeOffset, &
                            RFLO_GetPatchIndices, RFLO_GetPatchDirection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iPatch, i, j, k

! ... local variables
  INTEGER :: iLev, lbound, ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ijkC, ijkN
  INTEGER :: inode, jnode, knode

  REAL(RFREAL)          :: sgn, thrustMom, thrustPress, rhoua, rhova, rhowa, &
                           ua, va, wa, dS
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), sFace(:,:)

  TYPE(t_patch), POINTER :: patch

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_CalcThrust',&
  'RFLO_CalcThrust.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  cv => region%levels(iLev)%mixt%cv
  dv => region%levels(iLev)%mixt%dv

  thrustMom   = region%global%thrustMom
  thrustPress = region%global%thrustPress

! loop over patches -----------------------------------------------------------

  DO iPatch=1,region%nPatches

    patch => region%levels(iLev)%patches(iPatch)

    IF (patch%thrustCalc) THEN

! --- get dimensions of the patch

      lbound = patch%lbound

      CALL RFLO_GetPatchIndices( region,patch,iLev, &
                                 ibeg,iend,jbeg,jend,kbeg,kend )
      CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )

! --- to take the right face vector and make it point outwards

      sgn   = +1._RFREAL
      inode = 0
      jnode = 0
      knode = 0
      IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
        sgn   = -1._RFREAL
        inode = -idir
        jnode = -jdir
        knode = -kdir
      ENDIF

! --- get the appropriate face vector

      IF (lbound==1 .OR. lbound==2) sFace => region%levels(iLev)%grid%si
      IF (lbound==3 .OR. lbound==4) sFace => region%levels(iLev)%grid%sj
      IF (lbound==5 .OR. lbound==6) sFace => region%levels(iLev)%grid%sk

! --- loop over all cells of the patch

      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkC   = IndIJK(i,j,k,iCOff,ijCOff)
            ijkN   = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
            rhoua  = cv(CV_MIXT_XMOM,ijkC)
            rhova  = cv(CV_MIXT_YMOM,ijkC)
            rhowa  = cv(CV_MIXT_ZMOM,ijkC)
            ua     = dv(DV_MIXT_UVEL,ijkC)
            va     = dv(DV_MIXT_VVEL,ijkC)
            wa     = dv(DV_MIXT_WVEL,ijkC)
            thrustMom = sgn*sFace(XCOORD,ijkN)*rhoua*ua + &
                        sgn*sFace(YCOORD,ijkN)*rhova*va + &
                        sgn*sFace(ZCOORD,ijkN)*rhowa*wa + thrustMom
            IF (region%global%thrustType == THRUST_MOMP) THEN
              dS = SQRT(sFace(XCOORD,ijkN)*sFace(XCOORD,ijkN)+ &
                        sFace(YCOORD,ijkN)*sFace(YCOORD,ijkN)+ &
                        sFace(ZCOORD,ijkN)*sFace(ZCOORD,ijkN))
              thrustPress = (dv(DV_MIXT_PRES,ijkC)- &
                             region%global%thrustPamb)*dS + thrustPress
            ENDIF
          ENDDO
        ENDDO
      ENDDO

    ENDIF   ! thrustCalc

  ENDDO     ! iPatch

! finalize --------------------------------------------------------------------

  region%global%thrustMom   = thrustMom
  region%global%thrustPress = thrustPress

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_CalcThrust

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_CalcThrust.F90,v $
! Revision 1.3  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.5  2003/11/20 16:40:37  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.1  2003/06/02 17:12:00  jblazek
! Added computation of thrust.
!
!******************************************************************************







