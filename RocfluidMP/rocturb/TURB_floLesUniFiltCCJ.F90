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
! Purpose: Perform cell to cell uniform filtering in J-direction.
!
! Description: The filtering is performed using fixed filter coefficients.
!              If filterwidth is zero in this direction the cell-variable is 
!              simply copied to the filtered cell-variable.
!
! Input: global  = global data, only used in register function 
!        nDum    = number of dummy cells
!        i,j,kbeg, i,j,kend = cell index boundaries
!        iCOff, ijCOff      = cell Offset    
!        nDel    = three components filter width parameter
!        idBeg   = begin variable index to be filtered
!        idEnd   = end variable index to be filtered
!        fVar    = cell variable to be filtered
!
! Output: filtVar  = filtered cell variable
!
! Notes: Mother routine = LesUniFiltCC.
!
!******************************************************************************
!
! $Id: TURB_floLesUniFiltCCJ.F90,v 1.5 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloLesUniFiltCCJ( global,nDum,ibeg,iend,jbeg,jend,kbeg,kend, &
                              iCOff,ijCOff,nDel,idBeg,idEnd,fact1,fact2,fVar, &
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
  INTEGER               :: nDum,ibeg,iend,jbeg,jend,kbeg,kend,iCOff,ijCOff
  INTEGER               :: nDel(DIRI:DIRK),idBeg,idEnd
  REAL(RFREAL)          :: fact1(FILWIDTH_FOUR),fact2(FILWIDTH_FOUR)
  REAL(RFREAL), POINTER :: fVar(:,:),filtVar(:,:)

! ... loop variables
  INTEGER :: i, j, k, l, ijkC, ijkCb, ijkCe, ijkCb1, ijkCe1

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  REAL(RFREAL)      :: tmpb(idBeg:idEnd,2),tmpe(idBeg:idEnd,2)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floLesUniFiltCCJ.F90,v $'

  CALL RegisterFunction( global,'TURB_FloLesUniFiltCCJ',&
  'TURB_floLesUniFiltCCJ.F90' )

! integration over J-direction -----------------------------------------------

  IF (nDel(DIRJ)==FILWIDTH_ZERO) THEN

! - no filtering in J-direction is needed and we copy fVar into filtVar

    filtVar(idBeg:idEnd,:)=fVar(idBeg:idEnd,:) ! no simple copy filtvar=fvar
                                               ! array sizes can be different

  ELSEIF ((nDel(DIRJ)==FILWIDTH_ONE) .OR. &
          (nDel(DIRJ)==FILWIDTH_TWO)) THEN

    DO k=kbeg,kend
      DO i=ibeg,iend

! ----- perform integration

        DO j=jbeg+1,jend-1
          ijkC = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          DO l=idBeg,idEnd
            filtVar(l,ijkC)=fact1(nDel(DIRJ))*(fVar(l,ijkC-iCOff) + &
                                               fVar(l,ijkC+iCOff))+ &
                            fact2(nDel(DIRJ))*fVar(l,ijkC)
          END DO ! l
        ENDDO    ! j
  
! ----- perform integration at the edges

        ijkCb = IndIJK(i    ,jbeg  ,k     ,iCOff,ijCOff)
        ijkCe = IndIJK(i    ,jend  ,k     ,iCOff,ijCOff)

        DO l=idBeg,idEnd
          tmpb(l,1)=2._RFREAL*fVar(l,ijkCb)-fVar(l,ijkCb+iCOff)
          tmpe(l,1)=2._RFREAL*fVar(l,ijkCe)-fVar(l,ijkCe-iCOff)
          filtVar(l,ijkCb)=fact1(nDel(DIRJ))*(tmpb(l,1)+fVar(l,ijkCb+iCOff))+ &
                           fact2(nDel(DIRJ))*fVar(l,ijkCb)
          filtVar(l,ijkCe)=fact1(nDel(DIRJ))*(tmpe(l,1)+fVar(l,ijkCe-iCOff))+ &
                           fact2(nDel(DIRJ))*fVar(l,ijkCe)
        ENDDO   ! l
         
      ENDDO     ! i
    ENDDO       ! k

  ELSEIF (nDel(DIRJ)==FILWIDTH_FOUR) THEN

    DO k=kbeg,kend
      DO i=ibeg,iend
  
! ----- perform integration

        DO j=jbeg+2,jend-2           
          ijkC = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          DO l=idBeg,idEnd
            filtVar(l,ijkC)=fact1(nDel(DIRJ))*(fVar(l,ijkC-2*iCOff) + &
                                               fVar(l,ijkC+2*iCOff))+ &
                            fact2(nDel(DIRJ))*(fVar(l,ijkC+iCOff) + &
                                               fVar(l,ijkC)+fVar(l,ijkC-iCOff))
          ENDDO ! l
        ENDDO   ! j
  
! ----- perform integration at the edges

        ijkCb = IndIJK(i    ,jbeg   ,k     ,iCOff,ijCOff)
        ijkCe = IndIJK(i    ,jend   ,k     ,iCOff,ijCOff)
        ijkCb1= IndIJK(i    ,jbeg+1 ,k     ,iCOff,ijCOff)
        ijkCe1= IndIJK(i    ,jend-1 ,k     ,iCOff,ijCOff)

        DO l=idBeg,idEnd
          tmpb(l,1)=2._RFREAL*fVar(l,ijkCb)-fVar(l,ijkCb+iCOff)
          tmpe(l,1)=2._RFREAL*fVar(l,ijkCe)-fVar(l,ijkCe-iCOff)
          filtVar(l,ijkCb1)=fact1(nDel(DIRJ))*(tmpb(l,1)+fVar(l,ijkCb1+2*iCOff))+ &
                           fact2(nDel(DIRJ))*(fVar(l,ijkCb1+iCOff)+fVar(l,ijkCb1)+ &
                                              fVar(l,ijkCb1-iCOff))
          filtVar(l,ijkCe1)=fact1(nDel(DIRJ))*(tmpe(l,1)+fVar(l,ijkCe1-2*iCOff))+ &
                           fact2(nDel(DIRJ))*(fVar(l,ijkCe1+iCOff)+fVar(l,ijkCe1)+ &
                                              fVar(l,ijkCe1-iCOff))
        ENDDO   ! l

        DO l=idBeg,idEnd
          tmpb(l,2)=2._RFREAL*fVar(l,ijkCb)-fVar(l,ijkCb+2*iCOff)
          tmpe(l,2)=2._RFREAL*fVar(l,ijkCe)-fVar(l,ijkCe-2*iCOff)
          filtVar(l,ijkCb)=fact1(nDel(DIRJ))*(tmpb(l,2)+fVar(l,ijkCb+2*iCOff))+ &
                           fact2(nDel(DIRJ))*(fVar(l,ijkCb+iCOff)+fVar(l,ijkCb)+ &
                                              tmpb(l,1))
          filtVar(l,ijkCe)=fact1(nDel(DIRJ))*(tmpe(l,2)+fVar(l,ijkCe-2*iCOff))+ &
                           fact2(nDel(DIRJ))*(fVar(l,ijkCe-iCOff)+fVar(l,ijkCe)+ &
                                              tmpe(l,1))
        ENDDO   ! l
      ENDDO     ! i
    ENDDO       ! k
  ENDIF         ! nDel

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloLesUniFiltCCJ

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floLesUniFiltCCJ.F90,v $
! Revision 1.5  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/06/03 02:12:49  wasistho
! expand CC-filtering to all dummy layers
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.3  2003/05/16 05:44:12  wasistho
! modified array range of CC-filtered
!
! Revision 1.2  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/14 23:55:30  wasistho
! Install Rocturb
!
!******************************************************************************







