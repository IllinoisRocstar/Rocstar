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
! Purpose: exchange Radiation transport variables between interior and dummy 
!          cells of the corresponding patches of the adjacent regions (both 
!          regions are on the same processor).
!
! Description: none.
!
! Input: region    = current region
!        regionSrc = source region
!        patch     = current patch of region
!        patchSrc  = source patch of regionSrc.
!
! Output: regions%levels%radi%cv = Radiation variables in dummy cells.
!
! Notes: intended for conforming grid boundaries only.
!
!******************************************************************************
!
! $Id: RADI_floFlimExchangeDummyConf.F90,v 1.3 2008/12/06 08:44:38 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FloFlimExchangeDummyConf( region,regionSrc,patch,patchSrc )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModIndexing, ONLY   : IndIJKMap
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset, RFLO_GetPatchMapping
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region, regionSrc
  TYPE(t_patch)  :: patch, patchSrc

! ... loop variables
  INTEGER :: idum, i, j, k, l, ii, jj, kk

! ... local variables
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idir, jdir, kdir, &
             iCOff, ijCOff, ijkD, iLev, nCv
  INTEGER :: ibegSrc, iendSrc, jbegSrc, jendSrc, kbegSrc, kendSrc, &
             idirSrc, jdirSrc, kdirSrc, iCOffSrc, ijCOffSrc, ijkCSrc
  INTEGER :: lb, lbs, l1SrcDir, l2SrcDir, mapMat(3,4)

  LOGICAL :: align

  REAL(RFREAL), POINTER :: rcv(:,:), rcvSrc(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RADI_FloFlimExchangeDummyConf', &
                         'RADI_floFlimExchangeDummyConf.F90' )

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

  nCv    =  region%radiInput%nCv

  rcv    => region%levels(iLev)%radi%cv
  rcvSrc => regionSrc%levels(iLev)%radi%cv

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

          DO l=1,nCv
            rcv(l,ijkD) = rcvSrc(l,ijkCSrc)
          ENDDO

        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k
  ENDDO        ! idum

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RADI_FloFlimExchangeDummyConf

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_floFlimExchangeDummyConf.F90,v $
! Revision 1.3  2008/12/06 08:44:38  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
!
!
!******************************************************************************







