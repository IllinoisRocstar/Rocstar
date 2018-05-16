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
! Purpose: compute artificial numerical dissipation for smoke.
!
! Description: the dissipation is 4th-order only, with no 2nd-order part
!
! Input: region = data of current region.
!
! Output: region%levels%peul%diss = dissipative fluxes.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_CentralDissipation.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_CentralDissipation( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, RFLO_GetCellOffset
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: i, j, k, ii, jj, kk, iCv

! ... local variables
  INTEGER, PARAMETER :: NPEUL_MAX = 10

  CHARACTER(CHRLEN)  :: RCSIdentString

  INTEGER :: ipcbeg,ipcend,jpcbeg,jpcend,kpcbeg,kpcend
  INTEGER :: nCv,iLev,iCOff,ijCOff,ijkC0,ijkCm1,ijkCp1,ijkCp2

  REAL(RFREAL) :: beta,eval,eps2,eps4,pmax,fd,vis4(NPEUL_MAX)
  REAL(RFREAL), POINTER :: sCv(:,:),sDiss(:,:),srad(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_CentralDissipation.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_CentralDissipation',&
  'PEUL_CentralDissipation.F90' )

! begin -----------------------------------------------------------------------

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  sCv  => region%levels(iLev)%peul%cv
  nCv  =  region%levels(iLev)%peul%nCv

  sDiss => region%levels(iLev)%peul%diss
  srad => region%levels(iLev)%peul%srad

  beta =  region%mixtInput%betrk(region%irkStep)

! the implementation assumes that each sCv is a smoke density, so check that
! nCv is indeed the number of smoke particle types

  IF (nCv /= region%peulInput%nPtypes) &
    CALL ErrorStop( global,ERR_PEUL_NPMISMATCH,__LINE__ )

  IF (nCv > NPEUL_MAX) &
    CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

  vis4(:)     = 0._RFREAL
  vis4(1:nCv) = beta*region%peulInput%ptypes(1:nCv)%vis4

! dissipation in i-direction ------------------------------------------------

  DO k=kpcbeg,kpcend
    DO j=jpcbeg,jpcend

! --- dissipative fluxes at I+1/2

      ii = 0
      DO i=ipcbeg-1,ipcend
        ii     = ii + 1
        ijkC0  = IndIJK(i  ,j,k,iCOff,ijCOff)
        ijkCm1 = IndIJK(i-1,j,k,iCOff,ijCOff)
        ijkCp1 = IndIJK(i+1,j,k,iCOff,ijCOff)
        ijkCp2 = IndIJK(i+2,j,k,iCOff,ijCOff)
        eval   = 0.5_RFREAL*(srad(ICOORD,ijkC0)+srad(ICOORD,ijkCp1) + &
                         MAX(srad(JCOORD,ijkC0)+srad(JCOORD,ijkCp1), &
                             srad(KCOORD,ijkC0)+srad(KCOORD,ijkCp1)))

        DO iCv = 1,nCv
          eps4 = eval*vis4(iCv)
          fd   = eps4*( (sCv(iCv,ijkCm1) - sCv(iCv,ijkCp2)) + &
              3._RFREAL*(sCv(iCv,ijkCp1) - sCv(iCv,ijkC0 )) )

          sDiss(iCv,ijkC0 ) = sDiss(iCv,ijkC0 ) + fd
          sDiss(iCv,ijkCp1) = sDiss(iCv,ijkCp1) - fd
        ENDDO ! iCv

      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! dissipation in j-direction ------------------------------------------------

  DO k=kpcbeg,kpcend
    DO i=ipcbeg,ipcend

! --- dissipative fluxes at J+1/2

      jj = 0
      DO j=jpcbeg-1,jpcend
        jj     = jj + 1
        ijkC0  = IndIJK(i,j  ,k,iCOff,ijCOff)
        ijkCm1 = IndIJK(i,j-1,k,iCOff,ijCOff)
        ijkCp1 = IndIJK(i,j+1,k,iCOff,ijCOff)
        ijkCp2 = IndIJK(i,j+2,k,iCOff,ijCOff)
        eval   = 0.5_RFREAL*(srad(JCOORD,ijkC0)+srad(JCOORD,ijkCp1) + &
                         MAX(srad(ICOORD,ijkC0)+srad(ICOORD,ijkCp1), &
                             srad(KCOORD,ijkC0)+srad(KCOORD,ijkCp1)))

        DO iCv = 1,nCv
          eps4 = eval*vis4(iCv)
          fd   = eps4*( (sCv(iCv,ijkCm1) - sCv(iCv,ijkCp2)) + &
              3._RFREAL*(sCv(iCv,ijkCp1) - sCv(iCv,ijkC0 )) )

          sDiss(iCv,ijkC0 ) = sDiss(iCv,ijkC0 ) + fd
          sDiss(iCv,ijkCp1) = sDiss(iCv,ijkCp1) - fd
        ENDDO ! iCv

      ENDDO   ! j
    ENDDO     ! i
  ENDDO       ! k

! dissipation in k-direction --------------------------------------------------

  DO j=jpcbeg,jpcend
    DO i=ipcbeg,ipcend

! --- dissipative fluxes at K+1/2

      kk = 0
      DO k=kpcbeg-1,kpcend
        kk     = kk + 1
        ijkC0  = IndIJK(i,j,k  ,iCOff,ijCOff)
        ijkCm1 = IndIJK(i,j,k-1,iCOff,ijCOff)
        ijkCp1 = IndIJK(i,j,k+1,iCOff,ijCOff)
        ijkCp2 = IndIJK(i,j,k+2,iCOff,ijCOff)
        eval   = 0.5_RFREAL*(srad(KCOORD,ijkC0)+srad(KCOORD,ijkCp1) + &
                         MAX(srad(ICOORD,ijkC0)+srad(ICOORD,ijkCp1), &
                             srad(JCOORD,ijkC0)+srad(JCOORD,ijkCp1)))

        DO iCv = 1,nCv
          eps4 = eval*vis4(iCv)
          fd   = eps4*( (sCv(iCv,ijkCm1) - sCv(iCv,ijkCp2)) + &
              3._RFREAL*(sCv(iCv,ijkCp1) - sCv(iCv,ijkC0 )) )

          sDiss(iCv,ijkC0 ) = sDiss(iCv,ijkC0 ) + fd
          sDiss(iCv,ijkCp1) = sDiss(iCv,ijkCp1) - fd
        ENDDO ! iCv

      ENDDO   ! k
    ENDDO     ! i
  ENDDO       ! j

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_CentralDissipation

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_CentralDissipation.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:25  haselbac
! Initial revision after changing case
!
! Revision 1.8  2004/07/28 15:42:13  jferry
! deleted defunct constructs: useDetangle, useSmokeDrag, useSmokeHeatTransfer
!
! Revision 1.7  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2004/03/02 21:48:09  jferry
! First phase of replacing Detangle interaction
!
! Revision 1.5  2003/09/25 15:42:57  jferry
! Added mixture source terms due to active smoke (for Detangle interaction)
!
! Revision 1.4  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.3  2003/05/01 22:57:24  jferry
! substituted macro for IndIJK
!
! Revision 1.2  2003/04/09 15:12:04  jferry
! miscellaneous stylistic changes
!
! Revision 1.1  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







