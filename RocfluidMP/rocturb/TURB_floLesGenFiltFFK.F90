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
! Purpose: Perform face to face filtering in K-direction.
!
! Description: The filtering is performed using the filter coefficients stored
!              in turb%ffCofk1,2,4. If filterwidth is zero in this direction
!              the cell-variable is simply copied to the filtered cell-variable.
!
! Input: region  = data of current region 
!        ijk     = ijk-face is being treated
!        i,j,kbeg, i,j,kend = node index boundaries
!        iNOff, ijNOff      = node Offset    
!        nDel    = three components filter width parameter
!        idBeg   = begin variable index to be filtered
!        idEnd   = end variable index to be filtered
!        fVar    = face variable to be filtered
!
! Output: filtVar  = filtered face variable
!
! Notes: Mother routine = LesGenFiltFF.
!
!******************************************************************************
!
! $Id: TURB_floLesGenFiltFFK.F90,v 1.4 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloLesGenFiltFFK( region,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
                                  iNOff,ijNOff,nDel,idBeg,idEnd,fVar,filtVar )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModTurbulence, ONLY : t_turb
  USE ModError
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), TARGET :: region
  INTEGER                :: ijk,ibeg,iend,jbeg,jend,kbeg,kend,iNOff,ijNOff
  INTEGER                :: nDel(DIRI:DIRK),idBeg,idEnd
  REAL(RFREAL), POINTER  :: fVar(:,:),filtVar(:,:)

! ... loop variables
  INTEGER :: i, j, k, l, ijkN, ijkNb, ijkNe

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_turb) , POINTER  :: turb

  INTEGER                 :: iLev
  REAL(RFREAL)            :: tmpb(idBeg:idEnd),tmpe(idBeg:idEnd)
  REAL(RFREAL),  POINTER  :: ffCof(:,:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floLesGenFiltFFK.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloLesGenFiltFFK',&
  'TURB_floLesGenFiltFFK.F90' )

! get parameters and pointers ------------------------------------------------

  iLev =  region%currLevel
  turb => region%levels(iLev)%turb

! integration over K-direction

  IF (nDel(DIRK)==FILWIDTH_ZERO) THEN

! - no filtering in K-direction is needed and we copy fVar into filtVar

    filtVar = fVar

  ELSEIF ((nDel(DIRK)==FILWIDTH_ONE) .OR. &
          (nDel(DIRK)==FILWIDTH_TWO)) THEN

    IF (nDel(DIRK)==FILWIDTH_ONE) THEN
      IF (ijk==DIRI) THEN
        ffCof => turb%ffCofk1I
      ELSEIF (ijk==DIRJ) THEN
        ffCof => turb%ffCofk1J
      ELSEIF (ijk==DIRK) THEN
        ffCof => turb%ffCofk1K
      ENDIF
    ELSEIF (nDel(DIRK)==FILWIDTH_TWO) THEN
      IF (ijk==DIRI) THEN
        ffCof => turb%ffCofk2I
      ELSEIF (ijk==DIRJ) THEN
        ffCof => turb%ffCofk2J
      ELSEIF (ijk==DIRK) THEN
        ffCof => turb%ffCofk2K
      ENDIF
    ENDIF

! - perform integration

    DO k=kbeg,kend
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkN = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          DO l=idBeg,idEnd
            filtVar(l,ijkN)=ffCof(1,ijkN)*fVar(l,ijkN-ijNOff) + &
                            ffCof(2,ijkN)*fVar(l,ijkN)        + &
                            ffCof(3,ijkN)*fVar(l,ijkN+ijNOff)
          END DO ! l
        ENDDO    ! i
      ENDDO      ! j
    ENDDO        ! k

  ELSEIF (nDel(DIRK)==FILWIDTH_FOUR) THEN

    IF (ijk==DIRI) THEN
      ffCof => turb%ffCofk4I
    ELSEIF (ijk==DIRJ) THEN
      ffCof => turb%ffCofk4J
    ELSEIF (ijk==DIRK) THEN
      ffCof => turb%ffCofk4K
    ENDIF
  
! - perform integration

    DO k=kbeg+1,kend-1           
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijkN = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          DO l=idBeg,idEnd
            filtVar(l,ijkN)=ffCof(1,ijkN)*fVar(l,ijkN-2*ijNOff) + &
                            ffCof(2,ijkN)*fVar(l,ijkN-ijNOff)   + &
                            ffCof(3,ijkN)*fVar(l,ijkN)          + &
                            ffCof(4,ijkN)*fVar(l,ijkN+ijNOff)   + &
                            ffCof(5,ijkN)*fVar(l,ijkN+2*ijNOff)
          ENDDO ! l
        ENDDO   ! i
      ENDDO     ! j
    ENDDO       ! k
  
! - perform integration at the edges

    DO j=jbeg,jend
      DO i=ibeg,iend
        ijkNb = IndIJK(i     ,j     ,kbeg  ,iNOff,ijNOff)
        ijkNe = IndIJK(i     ,j     ,kend  ,iNOff,ijNOff)
        DO l=idBeg,idEnd
          tmpb(l) = 2._RFREAL*fVar(l,ijkNb-ijNOff)-fVar(l,ijkNb)
          tmpe(l) = 2._RFREAL*fVar(l,ijkNe+ijNOff)-fVar(l,ijkNe)
          filtVar(l,ijkNb)=ffCof(1,ijkNb)*tmpb(l)               + &
                           ffCof(2,ijkNb)*fVar(l,ijkNb-ijNOff)  + &
                           ffCof(3,ijkNb)*fVar(l,ijkNb)         + &
                           ffCof(4,ijkNb)*fVar(l,ijkNb+ijNOff)  + &
                           ffCof(5,ijkNb)*fVar(l,ijkNb+2*ijNOff)
          filtVar(l,ijkNe)=ffCof(1,ijkNe)*fVar(l,ijkNe-2*ijNOff)+ &
                           ffCof(2,ijkNe)*fVar(l,ijkNe-ijNOff)  + &
                           ffCof(3,ijkNe)*fVar(l,ijkNe)         + &
                           ffCof(4,ijkNe)*fVar(l,ijkNe+ijNOff)  + &
                           ffCof(5,ijkNe)*tmpe(l) 
        ENDDO ! l
      ENDDO   ! i
    ENDDO     ! j
  ENDIF       ! nDel

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloLesGenFiltFFK

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floLesGenFiltFFK.F90,v $
! Revision 1.4  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.3  2003/08/29 01:42:14  wasistho
! Added TARGET attribute to region variable, since pointers are cached into it
!
! Revision 1.2  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/14 23:55:30  wasistho
! Install Rocturb
!
!
!******************************************************************************







