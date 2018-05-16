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
! Purpose: Get face values of cv for non-uniform grid.
!
! Description: It is obtained by two point averaging the adjacent cell values
!              using the averaging coefficients avgCoI,J,K.
!
! Input: region = data of current region
!        ijk    = ijk-face is being treated
!
! Output: fVar = face values of cv.
!
! Notes: This routine is currently only relevant if non-uniform filter is 
!        selected.
!
!******************************************************************************
!
! $Id: TURB_floLesGenC2F.F90,v 1.7 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloLesGenC2F( region,ijk )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensPhys, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: ijk

! ... loop variables
  INTEGER :: i, j, k, m, ijkC0, ijkC1, ijkN, ijkN1

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER           :: ibeg,iend,jbeg,jend,kbeg,kend, iadd,jadd,kadd
  INTEGER           :: idcbeg,idcend,jdcbeg,jdcend,kdcbeg,kdcend
  INTEGER           :: iLev,iCOff,ijCOff,iNOff,ijNOff
  INTEGER           :: indxC0(2),indxN(2),indxN1(2)
  REAL(RFREAL), POINTER :: avgCo(:,:),cv(:,:),fVar(:,:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floLesGenC2F.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloLesGenC2F',&
  'TURB_floLesGenC2F.F90' )

! get indices and pointers ---------------------------------------------------

  iLev   =  region%currLevel
  cv     => region%levels(iLev)%mixt%cv
  fVar   => region%levels(iLev)%turb%fVar

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  IF (ijk==DIRI) THEN
    ibeg = idcbeg+1
    iend = idcend
    jbeg = jdcbeg
    jend = jdcend
    kbeg = kdcbeg
    kend = kdcend
    iadd = -1
    jadd = 0
    kadd = 0
    avgCo => region%levels(iLev)%grid%c2fCoI
  ELSEIF (ijk==DIRJ) THEN
    ibeg = idcbeg
    iend = idcend
    jbeg = jdcbeg+1
    jend = jdcend
    kbeg = kdcbeg
    kend = kdcend
    iadd = 0
    jadd = -1
    kadd = 0
    avgCo => region%levels(iLev)%grid%c2fCoJ
  ELSEIF (ijk==DIRK) THEN
    ibeg = idcbeg
    iend = idcend
    jbeg = jdcbeg
    jend = jdcend
    kbeg = kdcbeg+1
    kend = kdcend
    iadd = 0
    jadd = 0
    kadd = -1
    avgCo => region%levels(iLev)%grid%c2fCoK
  ENDIF

! perform 2-point averaging of rho,rhou1,rhou2 and rhou3 from centers to faces

  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,iend

          ijkC0 = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          ijkC1 = ijkC0 + iadd + jadd*iCOff + kadd*ijCOff
          ijkN  = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)

          fVar(CV_TURB_DENS,ijkN)=avgCo(2,ijkN)*cv(CV_MIXT_DENS,ijkC0)+ &
                                  avgCo(1,ijkN)*cv(CV_MIXT_DENS,ijkC1)
          fVar(CV_TURB_XMOM,ijkN)=avgCo(2,ijkN)*cv(CV_MIXT_XMOM,ijkC0)+ &
                                  avgCo(1,ijkN)*cv(CV_MIXT_XMOM,ijkC1)
          fVar(CV_TURB_YMOM,ijkN)=avgCo(2,ijkN)*cv(CV_MIXT_YMOM,ijkC0)+ &
                                  avgCo(1,ijkN)*cv(CV_MIXT_YMOM,ijkC1)
          fVar(CV_TURB_ZMOM,ijkN)=avgCo(2,ijkN)*cv(CV_MIXT_ZMOM,ijkC0)+ &
                                  avgCo(1,ijkN)*cv(CV_MIXT_ZMOM,ijkC1)

        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k

!$startcopy TURB_FloLesUniC2F
! extrapolation to outerst layers

  IF (ijk==DIRI) THEN
    indxN1(1)=ibeg-1
    indxC0(1)=ibeg-1
    indxN(1) =ibeg
    indxN1(2)=iend+1
    indxC0(2)=iend
    indxN(2) =iend
    DO m=1,2
      DO k=kbeg,kend
        DO j=jbeg,jend
          ijkN  = IndIJK(indxN(m)   ,j     ,k     ,iNOff,ijNOff)
          ijkN1 = ijkN + indxN1(m)-indxN(m)
          ijkC0 = IndIJK(indxC0(m)  ,j     ,k     ,iCOff,ijCOff)

          fVar(CV_TURB_DENS,ijkN1)=2._RFREAL*cv(CV_MIXT_DENS,ijkC0)- &
                                           fVar(CV_TURB_DENS,ijkN)
          fVar(CV_TURB_XMOM,ijkN1)=2._RFREAL*cv(CV_MIXT_XMOM,ijkC0)- &
                                           fVar(CV_TURB_XMOM,ijkN)
          fVar(CV_TURB_YMOM,ijkN1)=2._RFREAL*cv(CV_MIXT_YMOM,ijkC0)- &
                                           fVar(CV_TURB_YMOM,ijkN)
          fVar(CV_TURB_ZMOM,ijkN1)=2._RFREAL*cv(CV_MIXT_ZMOM,ijkC0)- &
                                           fVar(CV_TURB_ZMOM,ijkN)
        ENDDO     ! j
      ENDDO       ! k
    ENDDO         ! m

    DO k=kdcbeg,kdcend+1
      DO i=idcbeg,idcend+1
        ijkN  = IndIJK(i    ,jdcend+1 ,k     ,iNOff,ijNOff)
        ijkN1 = ijkN - iNOff
        fVar(CV_TURB_DENS,ijkN) = fVar(CV_TURB_DENS,ijkN1)
        fVar(CV_TURB_XMOM,ijkN) = fVar(CV_TURB_XMOM,ijkN1)
        fVar(CV_TURB_YMOM,ijkN) = fVar(CV_TURB_YMOM,ijkN1) 
        fVar(CV_TURB_ZMOM,ijkN) = fVar(CV_TURB_ZMOM,ijkN1)
      ENDDO     ! i
    ENDDO       ! k
    DO j=jdcbeg,jdcend+1
      DO i=idcbeg,idcend+1
        ijkN  = IndIJK(i    ,j    ,kdcend+1  ,iNOff,ijNOff)
        ijkN1 = ijkN - ijNOff
        fVar(CV_TURB_DENS,ijkN) = fVar(CV_TURB_DENS,ijkN1)
        fVar(CV_TURB_XMOM,ijkN) = fVar(CV_TURB_XMOM,ijkN1)
        fVar(CV_TURB_YMOM,ijkN) = fVar(CV_TURB_YMOM,ijkN1) 
        fVar(CV_TURB_ZMOM,ijkN) = fVar(CV_TURB_ZMOM,ijkN1)
      ENDDO     ! i
    ENDDO       ! j

  ELSEIF (ijk==DIRJ) THEN
    indxN1(1)=jbeg-1
    indxC0(1)=jbeg-1
    indxN(1) =jbeg
    indxN1(2)=jend+1
    indxC0(2)=jend
    indxN(2) =jend
    DO m=1,2
      DO k=kbeg,kend
        DO i=ibeg,iend
          ijkN  = IndIJK(i     ,indxN(m)     ,k     ,iNOff,ijNOff)
          ijkN1 = ijkN + (indxN1(m)-indxN(m))*iNOff
          ijkC0 = IndIJK(i     ,indxC0(m)    ,k     ,iCOff,ijCOff)

          fVar(CV_TURB_DENS,ijkN1)=2._RFREAL*cv(CV_MIXT_DENS,ijkC0)- &
                                           fVar(CV_TURB_DENS,ijkN)
          fVar(CV_TURB_XMOM,ijkN1)=2._RFREAL*cv(CV_MIXT_XMOM,ijkC0)- &
                                           fVar(CV_TURB_XMOM,ijkN)
          fVar(CV_TURB_YMOM,ijkN1)=2._RFREAL*cv(CV_MIXT_YMOM,ijkC0)- &
                                           fVar(CV_TURB_YMOM,ijkN)
          fVar(CV_TURB_ZMOM,ijkN1)=2._RFREAL*cv(CV_MIXT_ZMOM,ijkC0)- &
                                           fVar(CV_TURB_ZMOM,ijkN)
        ENDDO     ! i
      ENDDO       ! k
    ENDDO         ! m

    DO k=kdcbeg,kdcend+1
      DO j=jdcbeg,jdcend+1
        ijkN  = IndIJK(idcend+1 ,j     ,k     ,iNOff,ijNOff)
        ijkN1 = ijkN - 1
        fVar(CV_TURB_DENS,ijkN) = fVar(CV_TURB_DENS,ijkN1)
        fVar(CV_TURB_XMOM,ijkN) = fVar(CV_TURB_XMOM,ijkN1)
        fVar(CV_TURB_YMOM,ijkN) = fVar(CV_TURB_YMOM,ijkN1) 
        fVar(CV_TURB_ZMOM,ijkN) = fVar(CV_TURB_ZMOM,ijkN1)
      ENDDO     ! j
    ENDDO       ! k
    DO j=jdcbeg,jdcend+1
      DO i=idcbeg,idcend+1
        ijkN  = IndIJK(i    ,j    ,kdcend+1  ,iNOff,ijNOff)
        ijkN1 = ijkN - ijNOff
        fVar(CV_TURB_DENS,ijkN) = fVar(CV_TURB_DENS,ijkN1)
        fVar(CV_TURB_XMOM,ijkN) = fVar(CV_TURB_XMOM,ijkN1)
        fVar(CV_TURB_YMOM,ijkN) = fVar(CV_TURB_YMOM,ijkN1) 
        fVar(CV_TURB_ZMOM,ijkN) = fVar(CV_TURB_ZMOM,ijkN1)
      ENDDO     ! i
    ENDDO       ! j

  ELSEIF (ijk==DIRK) THEN
    indxN1(1)=kbeg-1
    indxC0(1)=kbeg-1
    indxN(1) =kbeg
    indxN1(2)=kend+1
    indxC0(2)=kend
    indxN(2) =kend
    DO m=1,2
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkN  = IndIJK(i     ,j     ,indxN(m)     ,iNOff,ijNOff)
          ijkN1 = ijkN + (indxN1(m)-indxN(m))*ijNOff
          ijkC0 = IndIJK(i     ,j     ,indxC0(m)    ,iCOff,ijCOff)

          fVar(CV_TURB_DENS,ijkN1)=2._RFREAL*cv(CV_MIXT_DENS,ijkC0)- &
                                           fVar(CV_TURB_DENS,ijkN)
          fVar(CV_TURB_XMOM,ijkN1)=2._RFREAL*cv(CV_MIXT_XMOM,ijkC0)- &
                                           fVar(CV_TURB_XMOM,ijkN)
          fVar(CV_TURB_YMOM,ijkN1)=2._RFREAL*cv(CV_MIXT_YMOM,ijkC0)- &
                                           fVar(CV_TURB_YMOM,ijkN)
          fVar(CV_TURB_ZMOM,ijkN1)=2._RFREAL*cv(CV_MIXT_ZMOM,ijkC0)- &
                                           fVar(CV_TURB_ZMOM,ijkN)
        ENDDO     ! i
      ENDDO       ! j
    ENDDO         ! m

    DO k=kdcbeg,kdcend+1
      DO j=jdcbeg,jdcend+1
        ijkN  = IndIJK(idcend+1 ,j     ,k     ,iNOff,ijNOff)
        ijkN1 = ijkN - 1
        fVar(CV_TURB_DENS,ijkN) = fVar(CV_TURB_DENS,ijkN1)
        fVar(CV_TURB_XMOM,ijkN) = fVar(CV_TURB_XMOM,ijkN1)
        fVar(CV_TURB_YMOM,ijkN) = fVar(CV_TURB_YMOM,ijkN1) 
        fVar(CV_TURB_ZMOM,ijkN) = fVar(CV_TURB_ZMOM,ijkN1)
      ENDDO     ! j
    ENDDO       ! k
    DO k=kdcbeg,kdcend+1
      DO i=idcbeg,idcend+1
        ijkN  = IndIJK(i    ,jdcend+1  ,k     ,iNOff,ijNOff)
        ijkN1 = ijkN - iNOff
        fVar(CV_TURB_DENS,ijkN) = fVar(CV_TURB_DENS,ijkN1)
        fVar(CV_TURB_XMOM,ijkN) = fVar(CV_TURB_XMOM,ijkN1)
        fVar(CV_TURB_YMOM,ijkN) = fVar(CV_TURB_YMOM,ijkN1) 
        fVar(CV_TURB_ZMOM,ijkN) = fVar(CV_TURB_ZMOM,ijkN1)
      ENDDO     ! i
    ENDDO       ! k

  ENDIF         ! ijk
!$endcopy

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloLesGenC2F

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floLesGenC2F.F90,v $
! Revision 1.7  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2004/08/02 19:33:49  wasistho
! changed grid%avgCo to grid%c2fCo
!
! Revision 1.4  2004/07/30 17:46:48  wasistho
! replaced turb face averaging by rflo face averaging
!
! Revision 1.3  2004/03/19 02:52:51  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
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







