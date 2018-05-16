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
! Purpose: coarsen dummy cells at patches associated with interblock
!          or periodic boundaries (where real geometry was copied to).
!
! Description: none.
!
! Input: region = current region
!        iPatch = current patch
!
! Output: region%levels%grid%xyz = geometry in dummy cells at coarse levels.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_ExchangeGeometryLevels.F90,v 1.3 2008/12/06 08:44:26 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_ExchangeGeometryLevels( region,iPatch )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetPatchDirection, &
                            RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

  INTEGER :: iPatch

! ... loop variables
  INTEGER :: idum, iLev, i, j, k, ii, jj, kk

! ... local variables
  INTEGER :: ibegF, iendF, jbegF, jendF, kbegF, kendF, idir, jdir, kdir, &
             iNOffFine, ijNOffFine, ijkF
  INTEGER :: ibegC, iendC, jbegC, jendC, kbegC, kendC, &
             iNOffCoarse, ijNOffCoarse, ijkC

  REAL(RFREAL), POINTER :: xyzFine(:,:), xyzCoarse(:,:)

  TYPE(t_patch), POINTER :: patchFine, patchCoarse

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_ExchangeGeometryLevels',&
  'RFLO_ExchangeGeometryLevels.F90' )

! loop over dummy nodes of current patch

  DO iLev=2,region%nGridLevels

    xyzFine     => region%levels(iLev-1)%grid%xyz
    xyzCoarse   => region%levels(iLev  )%grid%xyz
    patchFine   => region%levels(iLev-1)%patches(iPatch)
    patchCoarse => region%levels(iLev  )%patches(iPatch)

    CALL RFLO_GetPatchIndicesNodes( region,patchFine  ,iLev-1,ibegF,iendF, &
                                    jbegF,jendF,kbegF,kendF )
    CALL RFLO_GetPatchIndicesNodes( region,patchCoarse,iLev  ,ibegC,iendC, &
                                    jbegC,jendC,kbegC,kendC )
    CALL RFLO_GetPatchDirection( patchFine,idir,jdir,kdir )
    CALL RFLO_GetNodeOffset( region,iLev-1,iNOffFine  ,ijNOffFine   )
    CALL RFLO_GetNodeOffset( region,iLev  ,iNOffCoarse,ijNOffCoarse )

    DO idum=1,region%nDumCells
      DO k=kbegC,kendC
        kk = (k-kbegC)*(2-ABS(kdir)) + kbegF
        DO j=jbegC,jendC
          jj = (j-jbegC)*(2-ABS(jdir)) + jbegF
          DO i=ibegC,iendC
            ii   = (i-ibegC)*(2-ABS(idir)) + ibegF
            ijkF = IndIJK(ii-idum*idir,jj-idum*jdir,kk-idum*kdir,iNOffFine ,ijNOffFine )
            ijkC = IndIJK( i-idum*idir, j-idum*jdir, k-idum*kdir,iNOffCoarse,ijNOffCoarse)

            xyzCoarse(XCOORD,ijkC) = xyzFine(XCOORD,ijkF)
            xyzCoarse(YCOORD,ijkC) = xyzFine(YCOORD,ijkF)
            xyzCoarse(ZCOORD,ijkC) = xyzFine(ZCOORD,ijkF)
          ENDDO  ! i
        ENDDO    ! j
      ENDDO      ! k
    ENDDO        ! idum

  ENDDO          ! iLev

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_ExchangeGeometryLevels

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ExchangeGeometryLevels.F90,v $
! Revision 1.3  2008/12/06 08:44:26  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:37  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.9  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.4  2002/09/27 00:57:10  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.3  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/06/22 01:13:38  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.1  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
!******************************************************************************







