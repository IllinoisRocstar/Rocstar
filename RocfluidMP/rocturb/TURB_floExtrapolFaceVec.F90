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
! Purpose: Extrapolate face vector variable to outerst dummy faces.
!
! Description: The values are already available up to second outerst dummy faces.
!              Here the values are extrapolated to the outers dummy faces.
!
! Input: region  = data of current region
!        intDIR  = identifier for integration direction
!        idBeg,idEnd = begin and end variable index of the vector variable
!        fVec        = vector variable to be extrapolated
!
! Output: fVec = values of fVec at the remaining dummy cells are defined.
!
! Notes: This routine is used by LesUniFiltFF and LesGenFiltFF.
!
!******************************************************************************
!
! $Id: TURB_floExtrapolFaceVec.F90,v 1.4 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloExtrapolFaceVec( region,intDIR,idBeg,idEnd,fVec )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset
  USE ModTurbulence
  USE ModError
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region)        :: region
  INTEGER               :: intDIR,idBeg,idEnd
  REAL(RFREAL), POINTER :: fVec(:,:)

! ... loop variables
  INTEGER :: i, j, k, l, m, ijkN, ijkN0, ijkN1, ijkN2
  INTEGER :: indxN, indxN0, indxN1, indxN2

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER           :: ibeg,iend,jbeg,jend,kbeg,kend
  INTEGER           :: idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend
  INTEGER           :: iLev,iNOff,ijNOff,nDum,nadd,sign(2),indx(2)
  INTEGER, PARAMETER:: TWO=2

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floExtrapolFaceVec.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloExtrapolFaceVec',&
  'TURB_floExtrapolFaceVec.F90' )

! get indices ---------------------------------------------------------------

  sign(1)= 1
  sign(2)=-1
  nDum   = region%nDumCells
  nadd   = -nDum
  iLev   = region%currLevel

  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                            jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  ibeg = idnbeg
  iend = idnend
  jbeg = jdnbeg
  jend = jdnend
  kbeg = kdnbeg
  kend = kdnend

! extrapolation to outerst layers

  IF (intDir==DIRI) THEN
      indx(1)=ibeg
      indx(2)=iend
      DO m=1,2
        indxN2 = indx(m)-TWO*sign(m)*nadd
        indxN1 = indx(m)-sign(m)*nadd
        indxN  = indx(m)
        DO k=kbeg,kend
          DO j=jbeg,jend
            ijkN  = IndIJK(indxN   ,j   ,k   ,iNOff,ijNOff)
            ijkN1 = ijkN + indxN1-indxN
            ijkN2 = ijkN + indxN2-indxN

            DO l=idBeg,idEnd
              fVec(l,ijkN)=2._RFREAL*fVec(l,ijkN1)-fVec(l,ijkN2) 
            ENDDO   ! l
          ENDDO     ! j
        ENDDO       ! k
      ENDDO         ! m
  ELSEIF (intDir==DIRJ) THEN
      indx(1)=jbeg
      indx(2)=jend
      DO m=1,2
        indxN2 = indx(m)-TWO*sign(m)*nadd
        indxN1 = indx(m)-sign(m)*nadd
        indxN  = indx(m)
        DO k=kbeg,kend
          DO i=ibeg,iend
            ijkN  = IndIJK(i   ,indxN    ,k   ,iNOff,ijNOff)
            ijkN1 = ijkN + (indxN1-indxN)*iNOff
            ijkN2 = ijkN + (indxN2-indxN)*iNOff

            DO l=idBeg,idEnd
              fVec(l,ijkN)=2._RFREAL*fVec(l,ijkN1)-fVec(l,ijkN2) 
            ENDDO   ! l
          ENDDO     ! i
        ENDDO       ! k
      ENDDO         ! m
  ELSEIF (intDir==DIRK) THEN
      indx(1)=kbeg
      indx(2)=kend
      DO m=1,2
        indxN2 = indx(m)-TWO*sign(m)*nadd
        indxN1 = indx(m)-sign(m)*nadd
        indxN  = indx(m)
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkN  = IndIJK(i   ,j   ,indxN     ,iNOff,ijNOff)
            ijkN1 = ijkN + (indxN1-indxN)*ijNOff
            ijkN2 = ijkN + (indxN2-indxN)*ijNOff

            DO l=idBeg,idEnd
              fVec(l,ijkN)=2._RFREAL*fVec(l,ijkN1)-fVec(l,ijkN2) 
            ENDDO   ! l
          ENDDO     ! i
        ENDDO       ! j
      ENDDO         ! m
  ENDIF             ! intDir

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloExtrapolFaceVec

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floExtrapolFaceVec.F90,v $
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
! Revision 1.2  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************







