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
! Purpose: Construct cell to face (two points) averaging coefficients.
!
! Description: Averaging coefficients are first computed for the interior 
!              i, j and k faces, then the patches and finally extrapolated
!              to edges and corners, depending on the patch bc.
!
! Input: region = info of current region data
!
! Output: c2fCoI, c2fCoJ, c2fCoK = averaging coefficients at i, j and k faces
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_C2fAvgCoeffs.F90,v 1.3 2008/12/06 08:44:25 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************
SUBROUTINE RFLO_C2fAvgCoeffs( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY      : t_region
  USE ModGlobal, ONLY          : t_global
  USE ModInterfaces, ONLY      : RFLO_GetDimensDummy, &
                                 RFLO_GetCellOffset, RFLO_GetNodeOffset, &
                                 RFLO_CopyEdgeFaceNorm, RFLO_CopyEdgeFaceParal
  USE RFLO_ModInterfacesSolver, ONLY: &
                                 RFLO_C2fAvgCoeffsPatch, &
                                 RFLO_C2fAvgCoeffsDummy, RFLO_C2fAvgCoeffsDegec
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, l, ipatch

! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: idcbeg,idcend,jdcbeg,jdcend,kdcbeg,kdcend
  INTEGER :: iLev,iCOff,ijCOff,iNOff,ijNOff
  INTEGER, PARAMETER :: IDIR=1, JDIR=2, KDIR=3

  REAL(RFREAL), POINTER :: xyz(:,:), avgCo(:,:)
  REAL(RFREAL), POINTER :: avgCoI(:,:), avgCoJ(:,:), avgCoK(:,:)
  REAL(RFREAL), PARAMETER :: REAL_SMALL = 1.E-16_RFREAL

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RFLO_C2fAvgCoeffs',&
  'RFLO_C2fAvgCoeffs.F90' )

! get indices and pointers ---------------------------------------------------

  iLev   =  region%currLevel
  xyz    => region%levels(iLev)%grid%xyz
  avgCoI => region%levels(iLev)%grid%c2fCoI
  avgCoJ => region%levels(iLev)%grid%c2fCoJ
  avgCoK => region%levels(iLev)%grid%c2fCoK

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                           jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! compute averaging coefficient at face i, j, k interior domain

  CALL ComputeAvgCo( IDIR )
  CALL ComputeAvgCo( JDIR )
  CALL ComputeAvgCo( KDIR )

! compute averaging coefficient at dummy points

  DO ipatch=1,region%nPatches
    CALL RFLO_C2fAvgCoeffsPatch( region,region%levels(ilev)%patches(ipatch) )
    CALL RFLO_C2fAvgCoeffsDummy( region,region%levels(ilev)%patches(ipatch) )
  ENDDO

! and copy to edge/corner dummies -------------------------------------------- 
! the edge length is extended that corners are covered to minimize work

  CALL RFLO_CopyEdgeFaceNorm( region,1,2,avgCoI,avgCoJ,avgCoK )
  CALL RFLO_CopyEdgeFaceParal( region,1,2,avgCoI,avgCoJ,avgCoK )

! correct averaging coefficients at degenerated edges/corners -----------------

  IF (global%degenrtEc) THEN
    CALL RFLO_C2fAvgCoeffsDegec( region )
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
    INTEGER   :: ibeg,iend,jbeg,jend,kbeg,kend, iadd,jadd,kadd
    INTEGER   :: i2,i3,i4,j2,j3,j4,k2,k3,k4
    INTEGER   :: corner(12)
    REAL(RFREAL) :: denom,segment(2),faceCenter(XCOORD:ZCOORD)
    REAL(RFREAL) :: cel0Center(XCOORD:ZCOORD),cel1Center(XCOORD:ZCOORD)

! - Set limits and pointers ---------------------------------------------------
    IF (ijk==IDIR) THEN
      ibeg = idcbeg+1
      iend = idcend
      jbeg = jdcbeg
      jend = jdcend
      kbeg = kdcbeg
      kend = kdcend
      iadd = -1
      jadd = 0
      kadd = 0
      i2   = 0
      i3   = 0
      i4   = 0
      j2   = 0
      j3   = 1
      j4   = 1
      k2   = 1
      k3   = 1
      k4   = 0
      avgCo => region%levels(iLev)%grid%c2fCoI
    ELSEIF (ijk==JDIR) THEN
      ibeg = idcbeg
      iend = idcend
      jbeg = jdcbeg+1
      jend = jdcend
      kbeg = kdcbeg
      kend = kdcend
      iadd = 0
      jadd = -1
      kadd = 0
      i2   = 0
      i3   = 1
      i4   = 1
      j2   = 0
      j3   = 0
      j4   = 0
      k2   = 1
      k3   = 1
      k4   = 0
      avgCo => region%levels(iLev)%grid%c2fCoJ
    ELSEIF (ijk==KDIR) THEN
      ibeg = idcbeg
      iend = idcend
      jbeg = jdcbeg
      jend = jdcend
      kbeg = kdcbeg+1
      kend = kdcend
      iadd = 0
      jadd = 0
      kadd = -1
      i2   = 0
      i3   = 1
      i4   = 1
      j2   = 1
      j3   = 1
      j4   = 0
      k2   = 0
      k3   = 0
      k4   = 0
      avgCo => region%levels(iLev)%grid%c2fCoK
    ENDIF

! define 2-point averaging coefficients for non-uniform grid

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          corner(1) = IndIJK(i    ,j    ,k    ,iNOff,ijNOff)
          corner(2) = IndIJK(i+i2 ,j+j2 ,k+k2 ,iNOff,ijNOff)
          corner(3) = IndIJK(i+i3 ,j+j3 ,k+k3 ,iNOff,ijNOff)
          corner(4) = IndIJK(i+i4 ,j+j4 ,k+k4 ,iNOff,ijNOff)
          corner(5) = IndIJK(i-iadd    ,j-jadd   ,k-kadd     ,iNOff,ijNOff)
          corner(6) = IndIJK(i+i2-iadd ,j+j2-jadd ,k+k2-kadd ,iNOff,ijNOff)
          corner(7) = IndIJK(i+i3-iadd ,j+j3-jadd ,k+k3-kadd ,iNOff,ijNOff)
          corner(8) = IndIJK(i+i4-iadd ,j+j4-jadd ,k+k4-kadd ,iNOff,ijNOff)
          corner(9) = IndIJK(i+iadd    ,j+jadd   ,k+kadd     ,iNOff,ijNOff)
          corner(10)= IndIJK(i+i2+iadd ,j+j2+jadd ,k+k2+kadd ,iNOff,ijNOff)
          corner(11)= IndIJK(i+i3+iadd ,j+j3+jadd ,k+k3+kadd ,iNOff,ijNOff)
          corner(12)= IndIJK(i+i4+iadd ,j+j4+jadd ,k+k4+kadd ,iNOff,ijNOff)
         
          DO l=1,3
            faceCenter(l) = 0.25_RFREAL*(xyz(l,corner(1))+ &
                                         xyz(l,corner(2))+ &
                                         xyz(l,corner(3))+ &
                                         xyz(l,corner(4)))
          ENDDO         
          DO l=1,3
            cel0Center(l) = 0.125_RFREAL*(xyz(l,corner(1))+ &
                                          xyz(l,corner(2))+ &
                                          xyz(l,corner(3))+ &
                                          xyz(l,corner(4))+ &
                                          xyz(l,corner(5))+ &
                                          xyz(l,corner(6))+ &
                                          xyz(l,corner(7))+ &
                                          xyz(l,corner(8)))
          ENDDO         
          DO l=1,3
            cel1Center(l) = 0.125_RFREAL*(xyz(l,corner(1)) + &
                                          xyz(l,corner(2)) + &
                                          xyz(l,corner(3)) + &
                                          xyz(l,corner(4)) + &
                                          xyz(l,corner(9)) + &
                                          xyz(l,corner(10))+ &
                                          xyz(l,corner(11))+ &
                                          xyz(l,corner(12)))
          ENDDO

          segment(1) = SQRT((faceCenter(XCOORD)-cel1Center(XCOORD))**2 + &
                            (faceCenter(YCOORD)-cel1Center(YCOORD))**2 + &
                            (faceCenter(ZCOORD)-cel1Center(ZCOORD))**2)

          segment(2) = SQRT((faceCenter(XCOORD)-cel0Center(XCOORD))**2 + &
                            (faceCenter(YCOORD)-cel0Center(YCOORD))**2 + &
                            (faceCenter(ZCOORD)-cel0Center(ZCOORD))**2)

!          IF (segment(1) < REAL_SMALL) THEN
!            segment(1) = 0._RFREAL
!          ENDIF
!          IF (segment(2) < REAL_SMALL) THEN
!            segment(2) = 0._RFREAL
!          ENDIF
!          IF ((segment(1)+segment(2)) < REAL_SMALL) THEN
!            denom = REAL_SMALL
!          ELSE
!            denom = (segment(1)+segment(2))        
!          ENDIF

          denom = MAX( segment(1)+segment(2),REAL_SMALL ) 
          avgCo(2,corner(1)) = segment(1)/denom
          avgCo(1,corner(1)) = 1._RFREAL-avgCo(2,corner(1))
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

  END SUBROUTINE ComputeAvgCo

END SUBROUTINE RFLO_C2fAvgCoeffs

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_C2fAvgCoeffs.F90,v $
! Revision 1.3  2008/12/06 08:44:25  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:36  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.6  2004/08/25 07:45:18  wasistho
! added RFLO_C2fAvgCoeffsDegec
!
! Revision 1.5  2004/08/09 21:43:51  wasistho
! removed TURB_ModParameters
!
! Revision 1.4  2004/08/03 03:17:49  wasistho
! replaced IF statement inside do loop by MAX statement
!
! Revision 1.3  2004/08/03 00:51:45  wasistho
! changed avgCo to c2fCo in the description
!
! Revision 1.2  2004/08/02 19:32:56  wasistho
! changed grid%avgCo to grid%c2fCo
!
! Revision 1.1  2004/07/30 17:30:31  wasistho
! initial import routines starting with RFLO_c2fAvg...
!
!
!******************************************************************************







