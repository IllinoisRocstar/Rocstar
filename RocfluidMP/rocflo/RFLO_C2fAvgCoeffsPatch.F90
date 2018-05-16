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
! Purpose: Construct averaging coeffients at patches of physical and symmetry
!          boundaries
!
! Description: For symmetry and physical bc with linear dummy extrapolation
!              the averaging coefficients are half-half, while for physical
!              bc with constant dummy extrapolation, they are one-zero.
!
! Input: region        = info of current region data
!        patch         = current patch data
!  
! Output: Averaging coefficients c2fCoI, c2fCoJ, c2fCoK at boundary patches. 
!
! Notes: Mother routine = RFLO_C2fAvgCoeffs.
!
!******************************************************************************
!
! $Id: RFLO_C2fAvgCoeffsPatch.F90,v 1.4 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_C2fAvgCoeffsPatch( region,patch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices,RFLO_GetPatchDirection, &
                            RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  TYPE(t_patch)  :: patch

! ... loop variables
  INTEGER :: i, j, k 

! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: iLev, lbound, bcType
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: inode, jnode, knode, idir, jdir, kdir, coIndx(2) 
  INTEGER :: iNOff, ijNOff, ijkN
  INTEGER :: iExtrap_SlipW,iExtrap_Inject

  REAL(RFREAL), POINTER :: avgCo(:,:)

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RFLO_C2fAvgCoeffsPatch',&
  'RFLO_C2fAvgCoeffsPatch.F90' )

! get dimensions and pointers-------------------------------------------------

  iLev   = region%currLevel
  lbound = patch%lbound
  bcType = patch%bcType

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend,jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  inode = 0
  jnode = 0
  knode = 0
  coIndx(1) = 1
  coIndx(2) = 2
  IF( lbound==2 .OR. lbound==4 .OR. lbound==6 ) THEN
    inode = -idir
    jnode = -jdir
    knode = -kdir
    coIndx(1) = 2
    coIndx(2) = 1
  ENDIF

  IF (lbound==1 .OR. lbound==2) THEN
    avgCo => region%levels(iLev)%grid%c2fCoI
  ELSE IF (lbound==3 .OR. lbound==4) THEN
    avgCo => region%levels(iLev)%grid%c2fCoJ
  ELSE IF (lbound==5 .OR. lbound==6) THEN
    avgCo => region%levels(iLev)%grid%c2fCoK
  ENDIF

  IF (bcType>=BC_SLIPWALL .AND. bcType<=BC_SLIPWALL+BC_RANGE) &
     iExtrap_SlipW  = patch%mixt%switches(BCSWI_SLIPW_EXTRAP)
  IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) &
     iExtrap_Inject = patch%mixt%switches(BCSWI_INJECT_EXTRAP)

  IF ((bcType>=BC_SLIPWALL .AND. bcType<=BC_SLIPWALL+BC_RANGE &
                           .AND. iExtrap_Slipw==EXTRAPOL_LINEAR) .OR. &
      (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) .OR. &
      (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE &
                            .AND. iExtrap_Inject==EXTRAPOL_LINEAR) .OR. &
      (bcType>=BC_INFLOW .AND. bcType<=BC_INFLOW+BC_RANGE) .OR. &
      (bcType>=BC_OUTFLOW .AND. bcType<=BC_OUTFLOW+BC_RANGE) .OR. &
      (bcType>=BC_FARFIELD .AND. bcType<=BC_FARFIELD+BC_RANGE) .OR. &
      (bcType>=BC_SYMMETRY .AND. bcType<=BC_SYMMETRY+BC_RANGE)) THEN

    CALL AvgCoPatchLinear

  ELSE IF ((bcType>=BC_SLIPWALL .AND. bcType<=BC_SLIPWALL+BC_RANGE &
                                .AND. iExtrap_Slipw==EXTRAPOL_CONST ) .OR. &
           (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE &
                                 .AND. iExtrap_Inject==EXTRAPOL_CONST )) THEN

    CALL AvgCoPatchConstant

  ELSE IF (bcType>=BC_REGIONINT .AND. bcType<=BC_REGIONINT+BC_RANGE) THEN

  ELSE IF (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) THEN

! everything else 

  ELSE

  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

! ==============================================================================
! completion subroutines for linear and constant extrap. case
! ==============================================================================

CONTAINS

  SUBROUTINE AvgCoPatchLinear

! define averaging coefficient at patch with linearly extrapolated dummies

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend

          ijkN  = IndIJK(i+inode  ,j+jnode  ,k+knode  ,iNOff,ijNOff)

          avgCo(coIndx(2),ijkN) = 0.5_RFREAL
          avgCo(coIndx(1),ijkN) = 0.5_RFREAL

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

  END SUBROUTINE AvgCoPatchLinear

! ##############################################################################

  SUBROUTINE AvgCoPatchConstant

! define averaging coefficient at patch with linearly extrapolated dummies

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend

          ijkN  = IndIJK(i+inode  ,j+jnode  ,k+knode  ,iNOff,ijNOff)

          avgCo(coIndx(2),ijkN) = 0._RFREAL
          avgCo(coIndx(1),ijkN) = 1._RFREAL

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

  END SUBROUTINE AvgCoPatchConstant

END SUBROUTINE RFLO_C2fAvgCoeffsPatch

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_C2fAvgCoeffsPatch.F90,v $
! Revision 1.4  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/08/19 15:39:26  mparmar
! Renamed patch variables
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.3  2004/08/03 00:51:51  wasistho
! changed avgCo to c2fCo in the description
!
! Revision 1.2  2004/08/02 19:33:02  wasistho
! changed grid%avgCo to grid%c2fCo
!
! Revision 1.1  2004/07/30 17:30:31  wasistho
! initial import routines starting with RFLO_c2fAvg...
!
!
!
!******************************************************************************







