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
! Purpose: Perform face to face filtering in J-direction.
!
! Description: The filtering is performed using the filter coefficients stored
!              in turb%ffCofj1,2,4. If filterwidth is zero in this direction
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
! $Id: TURB_floLesGenFiltFFJ.F90,v 1.4 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloLesGenFiltFFJ( region,ijk,ibeg,iend,jbeg,jend,kbeg,kend, &
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
  TYPE(t_turb) , POINTER :: turb

  INTEGER                :: iLev
  REAL(RFREAL)           :: tmpb(idBeg:idEnd),tmpe(idBeg:idEnd)
  REAL(RFREAL),  POINTER :: ffCof(:,:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floLesGenFiltFFJ.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloLesGenFiltFFJ',&
  'TURB_floLesGenFiltFFJ.F90' )

! get parameters and pointers ------------------------------------------------

  iLev =  region%currLevel
  turb => region%levels(iLev)%turb

! integration over J-direction

  IF (nDel(DIRJ)==FILWIDTH_ZERO) THEN

! - no filtering in J-direction is needed and we copy fVar into filtVar

    filtVar = fVar

  ELSEIF ((nDel(DIRJ)==FILWIDTH_ONE) .OR. &
          (nDel(DIRJ)==FILWIDTH_TWO)) THEN

    IF (nDel(DIRJ)==FILWIDTH_ONE) THEN
      IF (ijk==DIRI) THEN
        ffCof => turb%ffCofj1I
      ELSEIF (ijk==DIRJ) THEN
        ffCof => turb%ffCofj1J
      ELSEIF (ijk==DIRK) THEN
        ffCof => turb%ffCofj1K
      ENDIF
    ELSEIF (nDel(DIRJ)==FILWIDTH_TWO) THEN
      IF (ijk==DIRI) THEN
        ffCof => turb%ffCofj2I
      ELSEIF (ijk==DIRJ) THEN
        ffCof => turb%ffCofj2J
      ELSEIF (ijk==DIRK) THEN
        ffCof => turb%ffCofj2K
      ENDIF
    ENDIF

    DO k=kbeg,kend
      DO i=ibeg,iend

! ----- perform integration

        DO j=jbeg,jend
          ijkN = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          DO l=idBeg,idEnd
            filtVar(l,ijkN)=ffCof(1,ijkN)*fVar(l,ijkN-iNOff) + &
                            ffCof(2,ijkN)*fVar(l,ijkN)       + &
                            ffCof(3,ijkN)*fVar(l,ijkN+iNOff)
          END DO ! l
        ENDDO    ! j
         
      ENDDO      ! i
    ENDDO        ! k

  ELSEIF (nDel(DIRJ)==FILWIDTH_FOUR) THEN

    IF (ijk==DIRI) THEN
      ffCof => turb%ffCofj4I
    ELSEIF (ijk==DIRJ) THEN
      ffCof => turb%ffCofj4J
    ELSEIF (ijk==DIRK) THEN
      ffCof => turb%ffCofj4K
    ENDIF

    DO k=kbeg,kend
      DO i=ibeg,iend
  
! ----- perform integration

        DO j=jbeg+1,jend-1           
          ijkN = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          DO l=idBeg,idEnd
            filtVar(l,ijkN)=ffCof(1,ijkN)*fVar(l,ijkN-2*iNOff) + &
                            ffCof(2,ijkN)*fVar(l,ijkN-iNOff)   + &
                            ffCof(3,ijkN)*fVar(l,ijkN)         + &
                            ffCof(4,ijkN)*fVar(l,ijkN+iNOff)   + &
                            ffCof(5,ijkN)*fVar(l,ijkN+2*iNOff)
          ENDDO ! l
        ENDDO   ! j
  
! ----- perform integration at the edges

        ijkNb = IndIJK(i     ,jbeg  ,k     ,iNOff,ijNOff)
        ijkNe = IndIJK(i     ,jend  ,k     ,iNOff,ijNOff)
        DO l=idBeg,idEnd
          tmpb(l) = 2._RFREAL*fVar(l,ijkNb-iNOff)-fVar(l,ijkNb)
          tmpe(l) = 2._RFREAL*fVar(l,ijkNe+iNOff)-fVar(l,ijkNe)
          filtVar(l,ijkNb)=ffCof(1,ijkNb)*tmpb(l)              + &
                           ffCof(2,ijkNb)*fVar(l,ijkNb-iNOff)  + &
                           ffCof(3,ijkNb)*fVar(l,ijkNb)        + &
                           ffCof(4,ijkNb)*fVar(l,ijkNb+iNOff)  + &
                           ffCof(5,ijkNb)*fVar(l,ijkNb+2*iNOff)
          filtVar(l,ijkNe)=ffCof(1,ijkNe)*fVar(l,ijkNe-2*iNOff)+ &
                           ffCof(2,ijkNe)*fVar(l,ijkNe-iNOff)  + &
                           ffCof(3,ijkNe)*fVar(l,ijkNe)        + &
                           ffCof(4,ijkNe)*fVar(l,ijkNe+iNOff)  + &
                           ffCof(5,ijkNe)*tmpe(l)
        ENDDO   ! l
      ENDDO     ! i
    ENDDO       ! k
  ENDIF         ! nDel

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloLesGenFiltFFJ

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floLesGenFiltFFJ.F90,v $
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
! Revision 1.3  2003/08/29 01:42:10  wasistho
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







