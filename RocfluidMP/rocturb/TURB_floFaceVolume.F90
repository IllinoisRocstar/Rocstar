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
! Purpose: Define face volume.
!
! Description: The face volume is obtained by averaging the two adjacent cell 
!              volumes then extrapolated to outerst faces.
!
! Input: region  = data of current region
!        ijk     = ijk-face is being treated
!
! Output: Face volume defined at all nodes.
!
! Notes: Face volume is used in several routines, a.o LesMij & LesCalcEddyVis.
!
!******************************************************************************
!
! $Id: TURB_floFaceVolume.F90,v 1.4 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloFaceVolume( region,ijk )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, &
                            RFLO_GetNodeOffset,RFLO_GetCellOffset
  USE ModError
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: ijk

! ... loop variables
  INTEGER :: i, j, k, ijkC0, ijkC1, ijkN, ijkN1

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend, iadd,jadd,kadd
  INTEGER :: idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend
  INTEGER :: iLev,iNOff,ijNOff,iCOff,ijCOff

  REAL(RFREAL), POINTER :: vol(:), fvol(:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floFaceVolume.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloFaceVolume',&
  'TURB_floFaceVolume.F90' )

! get indices and pointers --------------------------------------------------

  iLev = region%currLevel
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  vol => region%levels(ilev)%grid%vol

  IF (ijk==DIRI) THEN
    fvol => region%levels(iLev)%turb%fvolI
    ibeg = idnbeg
    iend = idnend
    jbeg = jdnbeg
    jend = jdnend-1
    kbeg = kdnbeg
    kend = kdnend-1
    iadd = -1
    jadd = 0
    kadd = 0
  ELSEIF (ijk==DIRJ) THEN
    fvol => region%levels(iLev)%turb%fvolJ
    ibeg = idnbeg
    iend = idnend-1
    jbeg = jdnbeg
    jend = jdnend
    kbeg = kdnbeg
    kend = kdnend-1
    iadd = 0
    jadd = -1
    kadd = 0
  ELSEIF (ijk==DIRK) THEN
    fvol => region%levels(iLev)%turb%fvolK
    ibeg = idnbeg
    iend = idnend-1
    jbeg = jdnbeg
    jend = jdnend-1
    kbeg = kdnbeg
    kend = kdnend
    iadd = 0
    jadd = 0
    kadd = -1
  ENDIF

! Define face volumes and store in fvol(:)

  DO k=kbeg-kadd,kend+kadd
    DO j=jbeg-jadd,jend+jadd
      DO i=ibeg-iadd,iend+iadd
        ijkC0 = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
        ijkC1 = ijkC0 + iadd + jadd*iCOff + kadd*ijCOff
        ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
        fvol(ijkN) = 0.5_RFREAL*(vol(ijkC1)+vol(ijkC0))
      ENDDO
    ENDDO
  ENDDO

  IF (ijk==DIRI) THEN

! - extrapolation to outerst layers

    DO k=kbeg,kend
      DO j=jbeg,jend
        i = ibeg
        ijkC0 = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
        ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
        fvol(ijkN) = vol(ijkC0)
        i = iend
        ijkC1 = IndIJK(i+iadd,j+jadd,k+kadd,iCOff,ijCOff)
        ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
        fvol(ijkN) = vol(ijkC1)
      ENDDO
    ENDDO

    DO k=kdnbeg,kdnend
      DO i=idnbeg,idnend
        ijkN  = IndIJK(i    ,jdnend ,k     ,iNOff,ijNOff)
        ijkN1 = ijkN - iNOff
        fvol(ijkN) = fvol(ijkN1)
        fvol(ijkN) = fvol(ijkN1)
        fvol(ijkN) = fvol(ijkN1) 
        fvol(ijkN) = fvol(ijkN1)
      ENDDO     ! i
    ENDDO       ! k
    DO j=jdnbeg,jdnend
      DO i=idnbeg,idnend
        ijkN  = IndIJK(i    ,j    ,kdnend  ,iNOff,ijNOff)
        ijkN1 = ijkN - ijNOff
        fvol(ijkN) = fvol(ijkN1)
        fvol(ijkN) = fvol(ijkN1)
        fvol(ijkN) = fvol(ijkN1) 
        fvol(ijkN) = fvol(ijkN1)
      ENDDO     ! i
    ENDDO       ! j

  ELSEIF (ijk==DIRJ) THEN

! - extrapolation to outerst layers

    DO k=kbeg,kend
      DO i=ibeg,iend
        j = jbeg
        ijkC0 = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
        ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
        fvol(ijkN) = vol(ijkC0)
        j = jend
        ijkC1 = IndIJK(i+iadd,j+jadd,k+kadd,iCOff,ijCOff)
        ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
        fvol(ijkN) = vol(ijkC1)
      ENDDO
    ENDDO

    DO k=kdnbeg,kdnend
      DO j=jdnbeg,jdnend
        ijkN  = IndIJK(idnend  ,j     ,k     ,iNOff,ijNOff)
        ijkN1 = ijkN - 1
        fvol(ijkN) = fvol(ijkN1)
        fvol(ijkN) = fvol(ijkN1)
        fvol(ijkN) = fvol(ijkN1) 
        fvol(ijkN) = fvol(ijkN1)
      ENDDO     ! j
    ENDDO       ! k
    DO j=jdnbeg,jdnend
      DO i=idnbeg,idnend
        ijkN  = IndIJK(i    ,j    ,kdnend  ,iNOff,ijNOff)
        ijkN1 = ijkN - ijNOff
        fvol(ijkN) = fvol(ijkN1)
        fvol(ijkN) = fvol(ijkN1)
        fvol(ijkN) = fvol(ijkN1) 
        fvol(ijkN) = fvol(ijkN1)
      ENDDO     ! i
    ENDDO       ! j

  ELSEIF (ijk==DIRK) THEN

! - extrapolation to outerst layers

    DO j=jbeg,jend
      DO i=ibeg,iend
        k = kbeg
        ijkC0 = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
        ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
        fvol(ijkN) = vol(ijkC0)
        k = kend
        ijkC1 = IndIJK(i+iadd,j+jadd,k+kadd,iCOff,ijCOff)
        ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
        fvol(ijkN) = vol(ijkC1)
      ENDDO
    ENDDO

    DO k=kdnbeg,kdnend
      DO j=jdnbeg,jdnend
        ijkN  = IndIJK(idnend  ,j     ,k     ,iNOff,ijNOff)
        ijkN1 = ijkN - 1
        fvol(ijkN) = fvol(ijkN1)
        fvol(ijkN) = fvol(ijkN1)
        fvol(ijkN) = fvol(ijkN1) 
        fvol(ijkN) = fvol(ijkN1)
      ENDDO     ! j
    ENDDO       ! k
    DO k=kdnbeg,kdnend
      DO i=idnbeg,idnend
        ijkN  = IndIJK(i    ,jdnend  ,k     ,iNOff,ijNOff)
        ijkN1 = ijkN - iNOff
        fvol(ijkN) = fvol(ijkN1)
        fvol(ijkN) = fvol(ijkN1)
        fvol(ijkN) = fvol(ijkN1) 
        fvol(ijkN) = fvol(ijkN1)
      ENDDO     ! i
    ENDDO       ! k

  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FlofaceVolume

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floFaceVolume.F90,v $
! Revision 1.4  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.3  2003/10/09 23:07:08  wasistho
! renamed CalcEddyVis to LesCalcEddyVis
!
! Revision 1.2  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************







