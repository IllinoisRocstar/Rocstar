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
! Purpose: Perform face to face uniform filtering in J-direction.
!
! Description: The filtering is performed using fixed filter coefficients.
!              If filterwidth is zero in this direction the cell-variable is 
!              simply copied to the filtered cell-variable.
!
! Input: global  = global data, only used in register function 
!        i,j,kbeg, i,j,kend = node index boundaries
!        iNOff, ijNOff      = node Offset    
!        nDel    = three components filter width parameter
!        idBeg   = begin variable index to be filtered
!        idEnd   = end variable index to be filtered
!        fact1,2 = fixed filter coefficients for uniform grid
!        fVar    = face variable to be filtered
!
! Output: filtVar  = filtered face variable
!
! Notes: Mother routine = LesUniFiltFF.
!
!******************************************************************************
!
! $Id: TURB_floLesUniFiltFFJ.F90,v 1.4 2008/12/06 08:44:44 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloLesUniFiltFFJ( global,ibeg,iend,jbeg,jend,kbeg,kend,iNOff, &
                                  ijNOff,nDel,idBeg,idEnd,fact1,fact2,fVar, &
                                  filtVar )
  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModTurbulence
  USE ModError
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_global), POINTER :: global
  INTEGER               :: ibeg,iend,jbeg,jend,kbeg,kend,iNOff,ijNOff
  INTEGER               :: nDel(DIRI:DIRK),idBeg,idEnd
  REAL(RFREAL)          :: fact1(FILWIDTH_FOUR),fact2(FILWIDTH_FOUR)
  REAL(RFREAL), POINTER :: fVar(:,:),filtVar(:,:)

! ... loop variables
  INTEGER :: i, j, k, l, ijkN, ijkNb, ijkNe

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  REAL(RFREAL)      :: tmpb(idBeg:idEnd),tmpe(idBeg:idEnd)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floLesUniFiltFFJ.F90,v $'

  CALL RegisterFunction( global,'TURB_FloLesUniFiltFFJ',&
  'TURB_floLesUniFiltFFJ.F90' )

! integration over J-direction -----------------------------------------------

  IF (nDel(DIRJ)==FILWIDTH_ZERO) THEN

! - no filtering in J-direction is needed and we copy fVar into filtVar

    filtVar = fVar

  ELSEIF ((nDel(DIRJ)==FILWIDTH_ONE) .OR. &
          (nDel(DIRJ)==FILWIDTH_TWO)) THEN

    DO k=kbeg,kend
      DO i=ibeg,iend

! ----- perform integration

        DO j=jbeg,jend
          ijkN = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          DO l=idBeg,idEnd
            filtVar(l,ijkN)=fact1(nDel(DIRJ))*(fVar(l,ijkN-iNOff)+ &
                                               fVar(l,ijkN+iNOff))+ &
                            fact2(nDel(DIRJ))*fVar(l,ijkN)
          END DO ! l
        ENDDO    ! j
         
      ENDDO      ! i
    ENDDO        ! k

  ELSEIF (nDel(DIRJ)==FILWIDTH_FOUR) THEN

    DO k=kbeg,kend
      DO i=ibeg,iend
  
! ----- perform integration

        DO j=jbeg+1,jend-1           
          ijkN = IndIJK(i     ,j     ,k     ,iNOff,ijNOff)
          DO l=idBeg,idEnd
            filtVar(l,ijkN)=fact1(nDel(DIRJ))*(fVar(l,ijkN-2*iNOff)+ &
                                               fVar(l,ijkN+2*iNOff))+ &
                            fact2(nDel(DIRJ))*(fVar(l,ijkN+iNOff)+ &
                                               fVar(l,ijkN)+fVar(l,ijkN-iNOff))
          ENDDO ! l
        ENDDO   ! j
  
! ----- perform integration at the edges

        ijkNb = IndIJK(i     ,jbeg  ,k     ,iNOff,ijNOff)
        ijkNe = IndIJK(i     ,jend  ,k     ,iNOff,ijNOff)
        DO l=idBeg,idEnd
          tmpb(l) = 2._RFREAL*fVar(l,ijkNb-iNOff)-fVar(l,ijkNb)
          tmpe(l) = 2._RFREAL*fVar(l,ijkNe+iNOff)-fVar(l,ijkNe)
          filtVar(l,ijkNb)=fact1(nDel(DIRJ))*(tmpb(l)+fVar(l,ijkNb+2*iNOff))+ &
                           fact2(nDel(DIRJ))*(fVar(l,ijkNb+iNOff)+ &
                                 fVar(l,ijkNb)+fVar(l,ijkNb-iNOff))
          filtVar(l,ijkNe)=fact1(nDel(DIRJ))*(tmpe(l)+fVar(l,ijkNe-2*iNOff))+ &
                           fact2(nDel(DIRJ))*(fVar(l,ijkNe+iNOff)+ &
                                fVar(l,ijkNe)+fVar(l,ijkNe-iNOff))
        ENDDO   ! l
      ENDDO     ! i
    ENDDO       ! k
  ENDIF         ! nDel

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloLesUniFiltFFJ

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floLesUniFiltFFJ.F90,v $
! Revision 1.4  2008/12/06 08:44:44  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
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
! Revision 1.1  2002/10/14 23:55:30  wasistho
! Install Rocturb
!
!******************************************************************************







