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
! Purpose: allocate memory for all variables associated with the Eulerian
!          particles for all active regions on current processor.
!
! Description: none.
!
! Input: region = current region
!
! Output: region%levels(iLev)%peul = Eulerian particle variables
!         region%levels(iLev)%grid = grid variables
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_AllocateMemory.F90,v 1.3 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_AllocateMemory( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal,     ONLY : t_global
  USE ModGrid,       ONLY : t_grid
  USE ModPartEul,    ONLY : t_peul
  USE ModError
  USE ModParameters

  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset, &
                            RFLO_GetDimensDummy, RFLO_GetCellOffset
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: iLev

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ibc, iec, ibn, ien, iCellOffset, ijCellOffset
  INTEGER :: iNodeOffset, ijNodeOffset, errorFlag

  TYPE(t_level),  POINTER :: level
  TYPE(t_grid) ,  POINTER :: grid, gridOld
  TYPE(t_peul) ,  POINTER :: peul
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_AllocateMemory.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_AllocateMemory',&
  'PEUL_AllocateMemory.F90' )

! begin -----------------------------------------------------------------------

! loop over all grid levels

  DO iLev=1,region%nGridLevels

    level   => region%levels(iLev)
    grid    => level%grid
    gridOld => level%gridOld
    peul    => level%peul

! - get cell and node dimensions

    CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( region,iLev,iCellOffset,ijCellOffset )
    ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCellOffset,ijCellOffset)
    iec = IndIJK(idcend,jdcend,kdcend,iCellOffset,ijCellOffset)

    CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetNodeOffset( region,iLev,iNodeOffset,ijNodeOffset )
    ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNodeOffset,ijNodeOffset)
    ien = IndIJK(idnend,jdnend,kdnend,iNodeOffset,ijNodeOffset)

! - spectral radii

    ALLOCATE( peul%srad(3,ibc:iec),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! - coefficients of implicit residual smoothing

    IF (global%flowType==FLOW_STEADY .AND. &
        region%peulInput%smoocf>0._RFREAL) THEN
      ALLOCATE( peul%epsIrs(3,ibc:iec), stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    ELSE
      NULLIFY( peul%epsIrs )
    ENDIF

! - Eulerian particle variables

    ALLOCATE( peul%cv   (peul%nCv,ibc:iec),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    ALLOCATE( peul%cvOld(peul%nCv,ibc:iec),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    ALLOCATE( peul%rhs  (peul%nCv,ibc:iec),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    ALLOCATE( peul%diss (peul%nCv,ibc:iec),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

    IF (peul%nDv > 0) THEN
      ALLOCATE( peul%dv(peul%nDv,ibc:iec),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    ELSE
      NULLIFY( peul%dv )
    ENDIF

    IF (peul%nTv > 0) THEN
      ALLOCATE( peul%tv(peul%nTv,ibc:iec),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    ELSE
      NULLIFY( peul%tv )
    ENDIF

    IF (iLev>1 .AND. global%cycleType/=MGCYCLE_NO) THEN
      ALLOCATE( peul%fterm(peul%nCv,ibc:iec),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    ELSE
      NULLIFY( peul%fterm )
    ENDIF

    IF (global%flowType == FLOW_UNSTEADY) THEN
      ALLOCATE( peul%rhsSum(peul%nCv,ibc:iec),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    ELSE
      NULLIFY( peul%rhsSum )
    ENDIF

  ENDDO   ! iLev

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_AllocateMemory

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_AllocateMemory.F90,v $
! Revision 1.3  2008/12/06 08:44:38  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:13  haselbac
! Initial revision after changing case
!
! Revision 1.5  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.4  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.3  2003/04/07 18:29:01  jferry
! added inflow boundary condition and initialization to a constant
!
! Revision 1.2  2003/02/12 23:34:48  jferry
! Replaced [io]stat=global%error with local errorFlag
!
! Revision 1.1  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







