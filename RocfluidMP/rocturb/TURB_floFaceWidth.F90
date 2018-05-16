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
! Purpose: Obtained the width of all faces.
!
! Description: Face width is defined as distance between opposing edge-mids.
!              Each face has two face widths, e.g. face i has fw in j & k dir. 
!
! Input: region  = data of current region
!
! Output: Face width stored in workI,J,K(2,:)
!
! Notes: This routine is employed to construct non-uniform filter coefficients.
!        Hence is only relevant if non-uniform filter is selected.
!
!******************************************************************************
!
! $Id: TURB_floFaceWidth.F90,v 1.5 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************
SUBROUTINE TURB_FloFaceWidth( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY      : t_region
  USE ModGlobal, ONLY          : t_global
  USE ModInterfaces, ONLY      : RFLO_GetDimensDummy, RFLO_GetNodeOffset, &
                                 RFLO_CopyEdgeFaceNorm, RFLO_CopyEdgeFaceParal
  USE TURB_ModInterfaces, ONLY : TURB_FloFaceWidthDummy
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, ijkN, ipatch

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER           :: ibeg,iend,jbeg,jend,kbeg,kend
  INTEGER           :: idcbeg,idcend,jdcbeg,jdcend,kdcbeg,kdcend
  INTEGER           :: iLev,iNOff,ijNOff
  REAL(RFREAL), POINTER :: xyz(:,:), fw(:,:)
  REAL(RFREAL), POINTER :: fwI(:,:), fwJ(:,:), fwK(:,:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floFaceWidth.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloFaceWidth',&
  'TURB_floFaceWidth.F90' )

! get indices and pointers ---------------------------------------------------

  iLev =  region%currLevel
  xyz  => region%levels(iLev)%grid%xyz
  fwI  => region%levels(iLev)%turb%workI
  fwJ  => region%levels(iLev)%turb%workJ
  fwK  => region%levels(iLev)%turb%workK

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                           jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! compute averaging coefficient at face i, j, k interior domain

  CALL ComputeFaceWidth( DIRI )
  CALL ComputeFaceWidth( DIRJ )
  CALL ComputeFaceWidth( DIRK )

! compute averaging coefficient at dummy points

  DO ipatch=1,region%nPatches
    CALL TURB_FloFaceWidthDummy( region,region%levels(ilev)%patches(ipatch) )
  ENDDO

! and copy to edge/corner dummies -------------------------------------------- 
! the edge length is extended that corners are covered to minimize work

  CALL RFLO_CopyEdgeFaceNorm( region,1,2,fwI,fwJ,fwK )
  CALL RFLO_CopyEdgeFaceParal( region,1,2,fwI,fwJ,fwK )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

! =============================================================================
!   Face width computation subroutine
! =============================================================================

CONTAINS

  SUBROUTINE ComputeFaceWidth( ijk )

! ... parameters
    INTEGER   :: ijk

! ... local variables
    INTEGER   :: i2,i3,i4,j2,j3,j4,k2,k3,k4
    INTEGER   :: corner(4)
    REAL(RFREAL) :: fc1(XCOORD:ZCOORD), fc2(XCOORD:ZCOORD)

! - Set limits and pointers ---------------------------------------------------
    IF (ijk==DIRI) THEN
      ibeg = idcbeg
      iend = idcend+1
      jbeg = jdcbeg
      jend = jdcend
      kbeg = kdcbeg
      kend = kdcend
      i2   = 0
      i3   = 0
      i4   = 0
      j2   = 0
      j3   = 1
      j4   = 1
      k2   = 1
      k3   = 1
      k4   = 0
      fw   => region%levels(iLev)%turb%workI
    ELSEIF (ijk==DIRJ) THEN
      ibeg = idcbeg
      iend = idcend
      jbeg = jdcbeg
      jend = jdcend+1
      kbeg = kdcbeg
      kend = kdcend
      i2   = 0
      i3   = 1
      i4   = 1
      j2   = 0
      j3   = 0
      j4   = 0
      k2   = 1
      k3   = 1
      k4   = 0
      fw   => region%levels(iLev)%turb%workJ
    ELSEIF (ijk==DIRK) THEN
      ibeg = idcbeg
      iend = idcend
      jbeg = jdcbeg
      jend = jdcend
      kbeg = kdcbeg
      kend = kdcend+1
      i2   = 0
      i3   = 1
      i4   = 1
      j2   = 1
      j3   = 1
      j4   = 0
      k2   = 0
      k3   = 0
      k4   = 0
      fw   => region%levels(iLev)%turb%workK
    ENDIF

! define 2 face-widths for each face of non-uniform grid;
! face I has face-widths in J and K direction, and so on

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          corner(1) = IndIJK(i    ,j    ,k    ,iNOff,ijNOff)
          corner(2) = IndIJK(i+i2 ,j+j2 ,k+k2 ,iNOff,ijNOff)
          corner(3) = IndIJK(i+i3 ,j+j3 ,k+k3 ,iNOff,ijNOff)
          corner(4) = IndIJK(i+i4 ,j+j4 ,k+k4 ,iNOff,ijNOff)
         
          fc1(XCOORD:ZCOORD) = 0.5_RFREAL* &
                               (xyz(XCOORD:ZCOORD,corner(1))+ &
                                xyz(XCOORD:ZCOORD,corner(2)))
          fc2(XCOORD:ZCOORD) = 0.5_RFREAL* &
                               (xyz(XCOORD:ZCOORD,corner(3))+ &
                                xyz(XCOORD:ZCOORD,corner(4)))
          fw(1,ijkN) = SQRT((fc2(XCOORD)-fc1(XCOORD))**2 + &
                            (fc2(YCOORD)-fc1(YCOORD))**2 + &
                            (fc2(ZCOORD)-fc1(ZCOORD))**2)
         
          fc1(XCOORD:ZCOORD) = 0.5_RFREAL* &
                               (xyz(XCOORD:ZCOORD,corner(1))+ &
                                xyz(XCOORD:ZCOORD,corner(4)))
          fc2(XCOORD:ZCOORD) = 0.5_RFREAL* &
                               (xyz(XCOORD:ZCOORD,corner(2))+ &
                                xyz(XCOORD:ZCOORD,corner(3)))
          fw(2,ijkN) = SQRT((fc2(XCOORD)-fc1(XCOORD))**2 + &
                            (fc2(YCOORD)-fc1(YCOORD))**2 + &
                            (fc2(ZCOORD)-fc1(ZCOORD))**2)

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

  END SUBROUTINE ComputeFaceWidth

END SUBROUTINE TURB_FlofaceWidth

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floFaceWidth.F90,v $
! Revision 1.5  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/08/04 02:49:20  wasistho
! removed turb%avgCoI,J,K as it is defined as grid%c2fCoI,J,K
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.2  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************







