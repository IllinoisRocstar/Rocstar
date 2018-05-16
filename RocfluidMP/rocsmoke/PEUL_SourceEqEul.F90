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
! Purpose: compute source terms associated with Equilibrium Eulerian method
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: region%levels%peul%rhs = source terms.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_SourceEqEul.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_SourceEqEul( region,ipt )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartEul,    ONLY : t_peul_ptype
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region
  INTEGER,        INTENT(IN)    :: ipt

! ... loop variables
  INTEGER :: ijkN,i,j,k

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iLev,ibn,ien,errorFlag,ic,ijkNpi,ijkNpj,ijkNpk
  INTEGER :: ipcbeg,ipcend,jpcbeg,jpcend,kpcbeg,kpcend,iCOff,ijCOff
  INTEGER :: idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend,iNOff,ijNOff

  REAL(RFREAL) :: ux,uy,uz,vx,vy,vz,wx,wy,wz,taufac,mass,tau,trg2av

  REAL(RFREAL), POINTER,     DIMENSION(:,:) :: gTv,gradi,gradj,gradk,sCv,sRhs
  REAL(RFREAL), POINTER,     DIMENSION(:)   :: vol
  REAL(RFREAL), ALLOCATABLE, DIMENSION(:)   :: trg2i,trg2j,trg2k

  TYPE(t_peul_ptype), POINTER :: ptype
  TYPE(t_global),     POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_SourceEqEul.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_SourceEqEul',&
  'PEUL_SourceEqEul.F90' )

! begin -----------------------------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

  vol   => region%levels(iLev)%grid%vol

  gTv   => region%levels(iLev)%mixt%tv

  gradi => region%levels(iLev)%mixt%gradi
  gradj => region%levels(iLev)%mixt%gradj
  gradk => region%levels(iLev)%mixt%gradk

  sCv   => region%levels(iLev)%peul%cv
  sRhs  => region%levels(iLev)%peul%rhs

  ptype => region%peulInput%ptypes(ipt)

! allocate temporary arrays ---------------------------------------------------

  ALLOCATE( trg2i(ibn:ien),stat=errorFlag )
  errorFlag = global%error
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( trg2j(ibn:ien),stat=errorFlag )
  errorFlag = global%error
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( trg2k(ibn:ien),stat=errorFlag )
  errorFlag = global%error
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  DO ijkN = ibn,ien

    ux  = gradi(GR_MIXT_UX,ijkN)
    uy  = gradi(GR_MIXT_UY,ijkN)
    uz  = gradi(GR_MIXT_UZ,ijkN)
    vx  = gradi(GR_MIXT_VX,ijkN)
    vy  = gradi(GR_MIXT_VY,ijkN)
    vz  = gradi(GR_MIXT_VZ,ijkN)
    wx  = gradi(GR_MIXT_WX,ijkN)
    wy  = gradi(GR_MIXT_WY,ijkN)
    wz  = gradi(GR_MIXT_WZ,ijkN)

    trg2i(ijkN) = ux**2 + vy**2 + wz**2 + 2._RFREAL*(uy*vx + vz*wy + wx*uz)

    ux  = gradj(GR_MIXT_UX,ijkN)
    uy  = gradj(GR_MIXT_UY,ijkN)
    uz  = gradj(GR_MIXT_UZ,ijkN)
    vx  = gradj(GR_MIXT_VX,ijkN)
    vy  = gradj(GR_MIXT_VY,ijkN)
    vz  = gradj(GR_MIXT_VZ,ijkN)
    wx  = gradj(GR_MIXT_WX,ijkN)
    wy  = gradj(GR_MIXT_WY,ijkN)
    wz  = gradj(GR_MIXT_WZ,ijkN)

    trg2j(ijkN) = ux**2 + vy**2 + wz**2 + 2._RFREAL*(uy*vx + vz*wy + wx*uz)

    ux  = gradk(GR_MIXT_UX,ijkN)
    uy  = gradk(GR_MIXT_UY,ijkN)
    uz  = gradk(GR_MIXT_UZ,ijkN)
    vx  = gradk(GR_MIXT_VX,ijkN)
    vy  = gradk(GR_MIXT_VY,ijkN)
    vz  = gradk(GR_MIXT_VZ,ijkN)
    wx  = gradk(GR_MIXT_WX,ijkN)
    wy  = gradk(GR_MIXT_WY,ijkN)
    wz  = gradk(GR_MIXT_WZ,ijkN)

    trg2k(ijkN) = ux**2 + vy**2 + wz**2 + 2._RFREAL*(uy*vx + vz*wy + wx*uz)

  END DO ! ijkN

  taufac = ptype%tauVcoef

  DO k=kpcbeg,kpcend
    DO j=jpcbeg,jpcend
      DO i=ipcbeg,ipcend

        ic = IndIJK(i,j,k,iCOff,ijCOff)
        mass = sCv(ipt,ic)*vol(ic)

        IF (mass > 0._RFREAL) THEN ! only effects positive masses

          tau = taufac / gTv(TV_MIXT_MUEL,ic)

          ijkN   = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
          ijkNpi = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
          ijkNpj = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
          ijkNpk = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)

          trg2av = (trg2i(ijkN) + trg2i(ijkNpi) + &
                    trg2j(ijkN) + trg2j(ijkNpj) + &
                    trg2k(ijkN) + trg2k(ijkNpk)) / 6._RFREAL

          sRhs(ipt,ic) = sRhs(ipt,ic) - mass * tau * trg2av

        END IF ! mass

      END DO ! i
    END DO   ! j
  END DO     ! k

! Deallocate temporary array --------------------------------------------------

  DEALLOCATE( trg2i,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  DEALLOCATE( trg2j,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  DEALLOCATE( trg2k,stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_SourceEqEul

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_SourceEqEul.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:59  haselbac
! Initial revision after changing case
!
! Revision 1.1  2004/05/04 20:27:37  jferry
! Implemented equilibrium Eulerian method
!
!******************************************************************************







