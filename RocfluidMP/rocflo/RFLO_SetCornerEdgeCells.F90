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
! Purpose: set values in edge and corner cells of a region.
!
! Description: this is for the case if the cells are not located
!              within another region.
!
! Input: region = region data.
!
! Output: new values of conservative variables.
!
! Notes: values in ALL edge and corner cells are averaged from neighboring
!        cells. Hence, this routine should be called first, before values
!        from adjacent regions are copied.
!
!******************************************************************************
!
! $Id: RFLO_SetCornerEdgeCells.F90,v 1.4 2008/12/06 08:44:28 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_SetCornerEdgeCells( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetDimensPhys, &
                            RFLO_GetCellOffset
  USE ModError
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iLev, gasModel, iCOff, ijCOff, ijkD, ijkC1, ijkC2

  REAL(RFREAL), POINTER :: cv(:,:)

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_SetCornerEdgeCells',&
  'RFLO_SetCornerEdgeCells.F90' )

! get dimensions and pointers -------------------------------------------------

  iLev = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetDimensPhys( region,iLev,ipcbeg,ipcend, &
                           jpcbeg,jpcend,kpcbeg,kpcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )

  gasModel = region%mixtInput%gasModel

  cv => region%levels(iLev)%mixt%cv

! edges 9, 10, 11, 12 ---------------------------------------------------------

  DO i=ipcbeg,ipcend
    DO j=jpcbeg-1,jdcbeg,-1         ! edges 9, 10
      DO k=kpcbeg-1,kdcbeg,-1
        ijkD  = IndIJK(i,j  ,k  ,iCOff,ijCOff)
        ijkC1 = IndIJK(i,j+1,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i,j  ,k+1,iCOff,ijCOff)
        CALL RFLO_AverageEdgeCells( ijkD,ijkC1,ijkC2,gasModel,cv,region )
      ENDDO
      DO k=kpcend+1,kdcend
        ijkD  = IndIJK(i,j,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i,j+1,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i,j  ,k-1,iCOff,ijCOff)
        CALL RFLO_AverageEdgeCells( ijkD,ijkC1,ijkC2,gasModel,cv,region )
      ENDDO
    ENDDO

    DO j=jpcend+1,jdcend           ! edges 11, 12
      DO k=kpcbeg-1,kdcbeg,-1
        ijkD  = IndIJK(i,j,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i,j-1,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i,j  ,k+1,iCOff,ijCOff)
        CALL RFLO_AverageEdgeCells( ijkD,ijkC1,ijkC2,gasModel,cv,region )
      ENDDO
      DO k=kpcend+1,kdcend
        ijkD  = IndIJK(i,j,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i,j-1,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i,j  ,k-1,iCOff,ijCOff)
        CALL RFLO_AverageEdgeCells( ijkD,ijkC1,ijkC2,gasModel,cv,region )
      ENDDO
    ENDDO
  ENDDO

! edges 2, 4, 6, 8 ------------------------------------------------------------

  DO j=jpcbeg,jpcend
    DO i=ipcbeg-1,idcbeg,-1        ! edges 2, 4
      DO k=kpcend+1,kdcend
        ijkD  = IndIJK(i  ,j,k  ,iCOff,ijCOff)
        ijkC1 = IndIJK(i+1,j,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j,k-1,iCOff,ijCOff)
        CALL RFLO_AverageEdgeCells( ijkD,ijkC1,ijkC2,gasModel,cv,region )
      ENDDO
      DO k=kpcbeg-1,kdcbeg,-1
        ijkD  = IndIJK(i  ,j,k  ,iCOff,ijCOff)
        ijkC1 = IndIJK(i+1,j,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j,k+1,iCOff,ijCOff)
        CALL RFLO_AverageEdgeCells( ijkD,ijkC1,ijkC2,gasModel,cv,region )
      ENDDO
    ENDDO

    DO i=ipcend+1,idcend           ! edges 6, 8
      DO k=kpcend+1,kdcend
        ijkD  = IndIJK(i  ,j,k  ,iCOff,ijCOff)
        ijkC1 = IndIJK(i-1,j,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j,k-1,iCOff,ijCOff)
        CALL RFLO_AverageEdgeCells( ijkD,ijkC1,ijkC2,gasModel,cv,region )
      ENDDO
      DO k=kpcbeg-1,kdcbeg,-1
        ijkD  = IndIJK(i  ,j,k  ,iCOff,ijCOff)
        ijkC1 = IndIJK(i-1,j,k  ,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j,k+1,iCOff,ijCOff)
        CALL RFLO_AverageEdgeCells( ijkD,ijkC1,ijkC2,gasModel,cv,region )
      ENDDO
    ENDDO
  ENDDO

! edges 1, 3, 5, 7 ------------------------------------------------------------

  DO k=kdcbeg,kdcend
    DO i=ipcbeg-1,idcbeg,-1        ! edges 1, 3
      DO j=jpcbeg-1,jdcbeg,-1
        ijkD  = IndIJK(i  ,j  ,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i+1,j  ,k,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j+1,k,iCOff,ijCOff)
        CALL RFLO_AverageEdgeCells( ijkD,ijkC1,ijkC2,gasModel,cv,region )
      ENDDO
      DO j=jpcend+1,jdcend
        ijkD  = IndIJK(i  ,j  ,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i+1,j  ,k,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j-1,k,iCOff,ijCOff)
        CALL RFLO_AverageEdgeCells( ijkD,ijkC1,ijkC2,gasModel,cv,region )
      ENDDO
    ENDDO

    DO i=ipcend+1,idcend           ! edges 5, 7
      DO j=jpcbeg-1,jdcbeg,-1
        ijkD  = IndIJK(i  ,j  ,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i-1,j  ,k,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j+1,k,iCOff,ijCOff)
        CALL RFLO_AverageEdgeCells( ijkD,ijkC1,ijkC2,gasModel,cv,region )
      ENDDO
      DO j=jpcend+1,jdcend
        ijkD  = IndIJK(i  ,j  ,k,iCOff,ijCOff)
        ijkC1 = IndIJK(i-1,j  ,k,iCOff,ijCOff)
        ijkC2 = IndIJK(i  ,j-1,k,iCOff,ijCOff)
        CALL RFLO_AverageEdgeCells( ijkD,ijkC1,ijkC2,gasModel,cv,region )
      ENDDO
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( region%global )

! =============================================================================
!   Averaging subroutine
! =============================================================================

  CONTAINS

    SUBROUTINE RFLO_AverageEdgeCells( ijkD,ijkC1,ijkC2,gasModel,cv,region )
      USE ModDataTypes
      USE ModDataStruct, ONLY : t_region
      USE ModInterfaces, ONLY : MixtureProperties
      USE ModParameters
      IMPLICIT NONE

      INTEGER               :: ijkD, ijkC1, ijkC2, gasModel
      REAL(RFREAL), POINTER :: cv(:,:)
      TYPE(t_region)        :: region

      cv(CV_MIXT_DENS,ijkD) = 0.5_RFREAL*(cv(CV_MIXT_DENS,ijkC1)+ &
                                          cv(CV_MIXT_DENS,ijkC2))
      cv(CV_MIXT_XMOM,ijkD) = 0.5_RFREAL*(cv(CV_MIXT_XMOM,ijkC1)+ &
                                          cv(CV_MIXT_XMOM,ijkC2))
      cv(CV_MIXT_YMOM,ijkD) = 0.5_RFREAL*(cv(CV_MIXT_YMOM,ijkC1)+ &
                                          cv(CV_MIXT_YMOM,ijkC2))
      cv(CV_MIXT_ZMOM,ijkD) = 0.5_RFREAL*(cv(CV_MIXT_ZMOM,ijkC1)+ &
                                          cv(CV_MIXT_ZMOM,ijkC2))
      cv(CV_MIXT_ENER,ijkD) = 0.5_RFREAL*(cv(CV_MIXT_ENER,ijkC1)+ &
                                          cv(CV_MIXT_ENER,ijkC2))

      IF (gasModel == GAS_MODEL_TCPERF) THEN
        CALL MixtureProperties( region,ijkD,ijkD,.false. )
      ELSE
        CALL MixtureProperties( region,ijkD,ijkD,.true.  )
      ENDIF

    END SUBROUTINE RFLO_AverageEdgeCells

END SUBROUTINE RFLO_SetCornerEdgeCells

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_SetCornerEdgeCells.F90,v $
! Revision 1.4  2008/12/06 08:44:28  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:39  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/10/31 21:09:36  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.1  2004/11/29 20:51:40  wasistho
! lower to upper case
!
! Revision 1.12  2003/11/20 16:40:40  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.8  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.7  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.6  2002/08/16 21:33:48  jblazek
! Changed interface to MixtureProperties.
!
! Revision 1.5  2002/03/18 23:11:33  jblazek
! Finished multiblock and MPI.
!
! Revision 1.4  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.3  2002/02/16 07:16:00  jblazek
! Added implicit residual smoothing.
!
! Revision 1.2  2002/02/06 00:15:39  jblazek
! Improved injection BC. Added pointers to gradients.
!
! Revision 1.1  2002/01/31 20:20:13  jblazek
! Added treatment of edge & corner cells.
!
!******************************************************************************







