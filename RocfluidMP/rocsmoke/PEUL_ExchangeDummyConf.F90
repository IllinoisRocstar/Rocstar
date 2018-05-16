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
! Purpose: exchange values between interior and dummy cells of the
!          corresponding patches of the adjacent regions (both regions
!          are on the same processor) specific to PEUL Module.
!
! Description: none.
!
! Input: region    = current region
!        regionSrc = source region
!        patch     = current patch of region
!        patchSrc  = source patch of regionSrc.
!
! Output: regions%levels%peul%cv = eulerian particle variables in dummy cells.
!
! Notes: intended for conforming grid boundaries only.
!
!******************************************************************************
!
! $Id: PEUL_ExchangeDummyConf.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_ExchangeDummyConf( region,regionSrc,patch,patchSrc )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModIndexing, ONLY   : IndIJKMap
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
        RFLO_GetCellOffset, RFLO_GetPatchMapping
  USE ModPartEul, ONLY    : t_peul
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), INTENT(INOUT) :: region
  TYPE(t_region), INTENT(IN)    :: regionSrc
  TYPE(t_patch),  INTENT(INOUT) :: patch
  TYPE(t_patch),  INTENT(IN )   :: patchSrc

! ... loop variables
  INTEGER :: iCv, idum, i, j, k, ii, jj, kk

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir, &
             iCOff, ijCOff, ijkD, iLev
  INTEGER :: ibegSrc, iendSrc, jbegSrc, jendSrc, kbegSrc, kendSrc, &
             idirSrc, jdirSrc, kdirSrc, iCOffSrc, ijCOffSrc, ijkCSrc
  INTEGER :: lb, lbs, l1SrcDir, l2SrcDir, mapMat(3,4)
  INTEGER :: nCv

  LOGICAL :: align

  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv, pCvSrc

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PEUL_ExchangeDummyConf.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PEUL_ExchangeDummyConf',&
  'PEUL_ExchangeDummyConf.F90' )

! check if the source region is active ----------------------------------------

  IF (regionSrc%active == OFF) THEN
    CALL ErrorStop( region%global,ERR_SRCREGION_OFF,__LINE__ )
  ENDIF

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend, &
                             jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchIndices( regionSrc,patchSrc,iLev,ibegSrc,iendSrc, &
                             jbegSrc,jendSrc,kbegSrc,kendSrc )
  CALL RFLO_GetPatchDirection( patch   ,idir   ,jdir,   kdir    )
  CALL RFLO_GetPatchDirection( patchSrc,idirSrc,jdirSrc,kdirSrc )
  CALL RFLO_GetCellOffset( region   ,iLev,iCOff   ,ijCOff    )
  CALL RFLO_GetCellOffset( regionSrc,iLev,iCOffSrc,ijCOffSrc )

  nCv  =  region%levels(iLev)%peul%nCv

  pCv    => region%levels(iLev)%peul%cv
  pCvSrc => regionSrc%levels(iLev)%peul%cv

! the implementation assumes that each cv is a smoke density, so check that
! nCv is indeed the number of smoke particle types

  IF (nCv /= region%peulInput%nPtypes) &
    CALL ErrorStop( region%global,ERR_PEUL_NPMISMATCH,__LINE__ )

! mapping between patches -----------------------------------------------------

                                       l1SrcDir =  1
  IF (patch%srcL1beg > patch%srcL1end) l1SrcDir = -1
                                       l2SrcDir =  1
  IF (patch%srcL2beg > patch%srcL2end) l2SrcDir = -1

  lb    = patch%lbound
  lbs   = patch%srcLbound
  align = patch%align

  CALL RFLO_GetPatchMapping( lb,lbs,l1SrcDir,l2SrcDir,align, &
                             idir,jdir,kdir,idirSrc,jdirSrc,kdirSrc, &
                             ibeg,iend,jbeg,jend,kbeg,kend, &
                             ibegSrc,iendSrc,jbegSrc,jendSrc,kbegSrc,kendSrc, &
                             mapMat )

! loop over dummy nodes of current patch --------------------------------------

  DO iCv=1,nCv
    DO idum=1,region%nDumCells
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ii      = i - idum*idir
            jj      = j - idum*jdir
            kk      = k - idum*kdir
            ijkD    = IndIJK(ii,jj,kk,iCOff,ijCOff)
            ijkCSrc = IndIJKMap(ii,jj,kk,mapMat,iCOffSrc,ijCOffSrc)

            pCv(iCv,ijkD) = pCvSrc(iCv,ijkCSrc)
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k
    ENDDO        ! idum
  ENDDO          ! iCv

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_ExchangeDummyConf

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_ExchangeDummyConf.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:35  haselbac
! Initial revision after changing case
!
! Revision 1.3  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.2  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.1  2003/04/09 14:33:58  fnajjar
! Initial Import of Multi-region rocsmoke
!
!******************************************************************************







