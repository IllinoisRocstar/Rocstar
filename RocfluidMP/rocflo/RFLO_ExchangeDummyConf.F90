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
!          are on the same processor).
!
! Description: none.
!
! Input: region    = current region
!        regionSrc = source region
!        patch     = current patch of region
!        patchSrc  = source patch of regionSrc.
!
! Output: regions%levels%mixt%cv = flow variables in dummy cells.
!
! Notes: intended for conforming grid boundaries only.
!
!******************************************************************************
!
! $Id: RFLO_ExchangeDummyConf.F90,v 1.4 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ExchangeDummyConf( region,regionSrc,patch,patchSrc )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModIndexing, ONLY   : IndIJKMap
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
        RFLO_GetCellOffset, RFLO_GetPatchMapping, MixtureProperties
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region, regionSrc
  TYPE(t_patch)  :: patch, patchSrc

! ... loop variables
  INTEGER :: idum, i, j, k, ii, jj, kk

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir, &
             iCOff, ijCOff, ijkD, iLev, gasModel
  INTEGER :: ibegSrc, iendSrc, jbegSrc, jendSrc, kbegSrc, kendSrc, &
             idirSrc, jdirSrc, kdirSrc, iCOffSrc, ijCOffSrc, ijkCSrc
  INTEGER :: lb, lbs, l1SrcDir, l2SrcDir, mapMat(3,4)

  LOGICAL :: align

  REAL(RFREAL), POINTER :: cv(:,:), cvSrc(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_ExchangeDummyConf',&
  'RFLO_ExchangeDummyConf.F90' )

! check if the source region is active

  IF (regionSrc%active == OFF) THEN
    CALL ErrorStop( region%global,ERR_SRCREGION_OFF,__LINE__ )
  ENDIF

! get dimensions and pointers

  iLev = region%currLevel

  CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend, &
                             jbeg,jend,kbeg,kend )
  CALL RFLO_GetPatchIndices( regionSrc,patchSrc,iLev,ibegSrc,iendSrc, &
                             jbegSrc,jendSrc,kbegSrc,kendSrc )
  CALL RFLO_GetPatchDirection( patch   ,idir   ,jdir,   kdir    )
  CALL RFLO_GetPatchDirection( patchSrc,idirSrc,jdirSrc,kdirSrc )
  CALL RFLO_GetCellOffset( region   ,iLev,iCOff   ,ijCOff    )
  CALL RFLO_GetCellOffset( regionSrc,iLev,iCOffSrc,ijCOffSrc )

  gasModel = region%mixtInput%gasModel

  cv    => region%levels(iLev)%mixt%cv
  cvSrc => regionSrc%levels(iLev)%mixt%cv

! mapping between patches

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

! loop over dummy nodes of current patch

  DO idum=1,region%nDumCells
    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ii      = i - idum*idir
          jj      = j - idum*jdir
          kk      = k - idum*kdir
          ijkD    = IndIJK(ii,jj,kk,iCOff,ijCOff)
          ijkCSrc = IndIJKMap(ii,jj,kk,mapMat,iCOffSrc,ijCOffSrc)

          cv(CV_MIXT_DENS,ijkD) = cvSrc(CV_MIXT_DENS,ijkCSrc)
          cv(CV_MIXT_XMOM,ijkD) = cvSrc(CV_MIXT_XMOM,ijkCSrc)
          cv(CV_MIXT_YMOM,ijkD) = cvSrc(CV_MIXT_YMOM,ijkCSrc)
          cv(CV_MIXT_ZMOM,ijkD) = cvSrc(CV_MIXT_ZMOM,ijkCSrc)
          cv(CV_MIXT_ENER,ijkD) = cvSrc(CV_MIXT_ENER,ijkCSrc)

          IF (gasModel == GAS_MODEL_TCPERF) THEN
            CALL MixtureProperties( region,ijkD,ijkD,.false. )
          ELSE
            CALL MixtureProperties( region,ijkD,ijkD,.true.  )
          ENDIF
        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! idum

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ExchangeDummyConf

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ExchangeDummyConf.F90,v $
! Revision 1.4  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/10/31 21:09:36  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.12  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.8  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.7  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.6  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.5  2002/08/16 21:33:48  jblazek
! Changed interface to MixtureProperties.
!
! Revision 1.4  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.3  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.2  2002/01/31 20:56:30  jblazek
! Added basic boundary conditions.
!
! Revision 1.1  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
!******************************************************************************







