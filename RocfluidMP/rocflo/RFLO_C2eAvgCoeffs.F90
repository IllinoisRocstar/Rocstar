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
! Purpose: Construct cell to edge (four points) averaging coefficients.
!
! Description: the coefficients are obtained following an i,j,k cycle:
!              coefficients surrounding edge-i are computed by first performing
!              cell2face averaging in the j direction and then face2edge in the
!              k direction. The same way for edge-j and edge-k.
!
! Input: region = info of current region data
!
! Output: c2eCoI, c2eCoJ, c2eCoK = averaging coefficients at i, j and k edges
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_C2eAvgCoeffs.F90,v 1.3 2008/12/06 08:44:25 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************
SUBROUTINE RFLO_C2eAvgCoeffs( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY      : t_region
  USE ModGlobal, ONLY          : t_global
  USE ModInterfaces, ONLY      : RFLO_GetDimensDummy, &
                                 RFLO_GetCellOffset, RFLO_GetNodeOffset
  USE RFLO_ModInterfacesSolver, ONLY: RFLO_C2eAvgCoeffsDegec
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: idcbeg,idcend,jdcbeg,jdcend,kdcbeg,kdcend
  INTEGER :: iLev,iCOff,ijCOff,iNOff,ijNOff
  INTEGER, PARAMETER :: IDIR=1, JDIR=2, KDIR=3

  REAL(RFREAL), POINTER :: avgCo(:,:), c2fCo2(:,:), c2fCo3(:,:)

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RFLO_C2eAvgCoeffs',&
  'RFLO_C2eAvgCoeffs.F90' )

! get indices and pointers ---------------------------------------------------

  iLev   =  region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                           jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! compute averaging coefficient at i, j, k edges

  CALL ComputeAvgCo( IDIR )
  CALL ComputeAvgCo( JDIR )
  CALL ComputeAvgCo( KDIR )

! correct averaging coefficients at degenerated edges/corners -----------------

  IF (global%degenrtEc) THEN
    CALL RFLO_C2eAvgCoeffsDegec( region )
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

! =============================================================================
!   Averaging coefficients computation subroutine
! =============================================================================

CONTAINS

  SUBROUTINE ComputeAvgCo( ijk )

! ... parameters
    INTEGER   :: ijk 

! ... local variables
    INTEGER   :: ibeg,iend,jbeg,jend,kbeg,kend
    INTEGER   :: i2,j2,k2, i3,j3,k3
    INTEGER   :: corner(3)
    REAL(RFREAL) :: ec0, ec1

! - Set limits and pointers ---------------------------------------------------
    IF (ijk==IDIR) THEN
      ibeg = idcbeg
      iend = idcend
      jbeg = jdcbeg+1
      jend = jdcend
      kbeg = kdcbeg+1
      kend = kdcend
      i2   = 0
      j2   = 0
      k2   = -1
      i3   = 0
      j3   = -1
      k3   = 0
      avgCo  => region%levels(iLev)%grid%c2eCoI
      c2fCo2 => region%levels(iLev)%grid%c2fCoJ
      c2fCo3 => region%levels(iLev)%grid%c2fCoK
    ELSEIF (ijk==JDIR) THEN
      ibeg = idcbeg+1
      iend = idcend
      jbeg = jdcbeg
      jend = jdcend
      kbeg = kdcbeg+1
      kend = kdcend
      i2   = -1
      j2   = 0
      k2   = 0
      i3   = 0
      j3   = 0
      k3   = -1
      avgCo  => region%levels(iLev)%grid%c2eCoJ
      c2fCo2 => region%levels(iLev)%grid%c2fCoK
      c2fCo3 => region%levels(iLev)%grid%c2fCoI
    ELSEIF (ijk==KDIR) THEN
      ibeg = idcbeg+1
      iend = idcend
      jbeg = jdcbeg+1
      jend = jdcend
      kbeg = kdcbeg
      kend = kdcend
      i2   = 0
      j2   = -1
      k2   = 0
      i3   = -1
      j3   = 0
      k3   = 0
      avgCo  => region%levels(iLev)%grid%c2eCoK
      c2fCo2 => region%levels(iLev)%grid%c2fCoI
      c2fCo3 => region%levels(iLev)%grid%c2fCoJ
    ENDIF

! - define 4-points averaging coefficients for non-uniform grid

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          corner(1) = IndIJK(i       ,j       ,k       ,iNOff,ijNOff)
          corner(2) = IndIJK(i+i2    ,j+j2    ,k+k2    ,iNOff,ijNOff)
          corner(3) = IndIJK(i+i3    ,j+j3    ,k+k3    ,iNOff,ijNOff)

          ec0 = 0.5_RFREAL*(c2fCo3(1,corner(1))+c2fCo3(1,corner(3)))
          ec1 = 1._RFREAL-ec0

          avgCo(1,corner(1)) = ec1*c2fCo2(2,corner(2))
          avgCo(2,corner(1)) = ec1*c2fCo2(1,corner(2))
          avgCo(3,corner(1)) = ec0*c2fCo2(1,corner(1))
          avgCo(4,corner(1)) = ec0*c2fCo2(2,corner(1))
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

  END SUBROUTINE ComputeAvgCo

END SUBROUTINE RFLO_C2eAvgCoeffs

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_C2eAvgCoeffs.F90,v $
! Revision 1.3  2008/12/06 08:44:25  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:36  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.3  2004/08/25 07:45:33  wasistho
! added RFLO_C2eAvgCoeffsDegec
!
! Revision 1.2  2004/08/09 21:44:00  wasistho
! removed TURB_ModParameters
!
! Revision 1.1  2004/08/03 22:46:50  wasistho
! added RFLO_c2eAvgCoeffs
!
!
!******************************************************************************







