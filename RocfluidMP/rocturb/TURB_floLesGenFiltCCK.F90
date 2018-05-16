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
! Purpose: Perform cell to cell non-unifrom filtering in K-direction.
!
! Description: The filtering is performed using the filter coefficients stored
!              in turb%ccCofk1,2,4. If filterwidth is zero in this direction
!              the cell-variable is simply copied to the filtered cell-variable.
!
! Input: region  = data of current region 
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
! Notes: Mother routine = LesGenFiltCC.
!
!******************************************************************************
!
! $Id: TURB_floLesGenFiltCCK.F90,v 1.5 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloLesGenFiltCCK( region,nDum,ibeg,iend,jbeg,jend,kbeg,kend, &
                                  iCOff,ijCOff,nDel,idBeg,idEnd,fVar,filtVar )

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
  INTEGER                :: nDum,ibeg,iend,jbeg,jend,kbeg,kend,iCOff,ijCOff
  INTEGER                :: nDel(DIRI:DIRK),idBeg,idEnd
  REAL(RFREAL), POINTER  :: fVar(:,:),filtVar(:,:)

! ... loop variables
  INTEGER :: i, j, k, l, ijkC, ijkCb, ijkCe, ijkCb1, ijkCe1

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_turb) , POINTER :: turb

  INTEGER                :: iLev
  REAL(RFREAL)           :: tmpb(idBeg:idEnd,2),tmpe(idBeg:idEnd,2)
  REAL(RFREAL),  POINTER :: ccCof(:,:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floLesGenFiltCCK.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloLesGenFiltCCK',&
  'TURB_floLesGenFiltCCK.F90' )

! get parameters and pointers ------------------------------------------------

  iLev =  region%currLevel
  turb => region%levels(iLev)%turb

! integration over K-direction

  IF (nDel(DIRK)==FILWIDTH_ZERO) THEN

! - no filtering in K-direction is needed and we copy fVar into filtVar

    filtVar(idBeg:idEnd,:)=fVar(idBeg:idEnd,:) ! no simple copy filtvar=fvar
                                               ! array sizes can be different

  ELSEIF ((nDel(DIRK)==FILWIDTH_ONE) .OR. &
          (nDel(DIRK)==FILWIDTH_TWO)) THEN

    IF (nDel(DIRI)==FILWIDTH_ONE) THEN
      ccCof => turb%ccCofk1
    ELSEIF (nDel(DIRI)==FILWIDTH_TWO) THEN
      ccCof => turb%ccCofk2
    ENDIF

    DO j=jbeg,jend
      DO i=ibeg,iend

! ----- perform integration

        DO k=kbeg+1,kend-1
          ijkC = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          DO l=idBeg,idEnd
            filtVar(l,ijkC)=ccCof(1,ijkC)*fVar(l,ijkC-ijCOff) + &
                            ccCof(2,ijkC)*fVar(l,ijkC)        + &
                            ccCof(3,ijkC)*fVar(l,ijkC+ijCOff)
          END DO ! l
        ENDDO    ! k
  
! ----- perform integration at the edges

        ijkCb = IndIJK(i     ,j     ,kbeg  ,iCOff,ijCOff)
        ijkCe = IndIJK(i     ,j     ,kend  ,iCOff,ijCOff)

        DO l=idBeg,idEnd
          tmpb(l,1)=2._RFREAL*fVar(l,ijkCb)-fVar(l,ijkCb+ijCOff)
          tmpe(l,1)=2._RFREAL*fVar(l,ijkCe)-fVar(l,ijkCe-ijCOff)
          filtVar(l,ijkCb)=ccCof(1,ijkCb)*tmpb(l,1)     + &
                           ccCof(2,ijkCb)*fVar(l,ijkCb) + &
                           ccCof(3,ijkCb)*fVar(l,ijkCb+ijCOff)
          filtVar(l,ijkCe)=ccCof(1,ijkCe)*fVar(l,ijkCe-ijCOff) + &
                           ccCof(2,ijkCe)*fVar(l,ijkCe)        + &
                           ccCof(3,ijkCe)*tmpe(l,1)
        ENDDO    ! l
      ENDDO      ! i
    ENDDO        ! j

  ELSEIF (nDel(DIRK)==FILWIDTH_FOUR) THEN

    ccCof => turb%ccCofk4

    DO j=jbeg,jend
      DO i=ibeg,iend
  
! ----- perform integration

        DO k=kbeg+2,kend-2
          ijkC = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)
          DO l=idBeg,idEnd
            filtVar(l,ijkC)=ccCof(1,ijkC)*fVar(l,ijkC-2*ijCOff) + &
                            ccCof(2,ijkC)*fVar(l,ijkC-ijCOff)   + &
                            ccCof(3,ijkC)*fVar(l,ijkC)          + &
                            ccCof(4,ijkC)*fVar(l,ijkC+ijCOff)   + &
                            ccCof(5,ijkC)*fVar(l,ijkC+2*ijCOff)
          ENDDO ! l
        ENDDO   ! k
  
! ----- perform integration at the edges

        ijkCb = IndIJK(i     ,j     ,kbeg   ,iCOff,ijCOff)
        ijkCe = IndIJK(i     ,j     ,kend   ,iCOff,ijCOff)
        ijkCb1= IndIJK(i     ,j     ,kbeg+1 ,iCOff,ijCOff)
        ijkCe1= IndIJK(i     ,j     ,kend-1 ,iCOff,ijCOff)

        DO l=idBeg,idEnd
          tmpb(l,1)=2._RFREAL*fVar(l,ijkCb)-fVar(l,ijkCb+ijCOff)
          tmpe(l,1)=2._RFREAL*fVar(l,ijkCe)-fVar(l,ijkCe-ijCOff)
          filtVar(l,ijkCb1)=ccCof(1,ijkCb1)*tmpb(l,1)               + &
                            ccCof(2,ijkCb1)*fVar(l,ijkCb1-ijCOff)   + &
                            ccCof(3,ijkCb1)*fVar(l,ijkCb1)          + &
                            ccCof(4,ijkCb1)*fVar(l,ijkCb1+ijCOFf)   + &
                            ccCof(5,ijkCb1)*fVar(l,ijkCb1+2*ijCOff) 
          filtVar(l,ijkCe1)=ccCof(1,ijkCe1)*fVar(l,ijkCe1-2*ijCOff) + &
                            ccCof(2,ijkCe1)*fVar(l,ijkCe1-ijCOff)   + &
                            ccCof(3,ijkCe1)*fVar(l,ijkCe1)          + &
                            ccCof(4,ijkCe1)*fVar(l,ijkCe1+ijCOff)   + &
                            ccCof(5,ijkCe1)*tmpe(l,1)
        ENDDO   ! l
  
        DO l=idBeg,idEnd
          tmpb(l,2)=2._RFREAL*fVar(l,ijkCb)-fVar(l,ijkCb+2*ijCOff)
          tmpe(l,2)=2._RFREAL*fVar(l,ijkCe)-fVar(l,ijkCe-2*ijCOff)
          filtVar(l,ijkCb)=ccCof(1,ijkCb)*tmpb(l,2)              + &
                           ccCof(2,ijkCb)*tmpb(l,1)              + &
                           ccCof(3,ijkCb)*fVar(l,ijkCb)          + &
                           ccCof(4,ijkCb)*fVar(l,ijkCb+ijCOff)   + &
                           ccCof(5,ijkCb)*fVar(l,ijkCb+2*ijCOff) 
          filtVar(l,ijkCe)=ccCof(1,ijkCe)*fVar(l,ijkCe-2*ijCOff) + &
                           ccCof(2,ijkCe)*fVar(l,ijkCe-ijCOff)   + &
                           ccCof(3,ijkCe)*fVar(l,ijkCe)          + &
                           ccCof(4,ijkCe)*tmpe(l,1)              + &
                           ccCof(5,ijkCe)*tmpe(l,2)
        ENDDO   ! l
      ENDDO     ! i
    ENDDO       ! j
  ENDIF         ! nDel

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_FloLesGenFiltCCK

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floLesGenFiltCCK.F90,v $
! Revision 1.5  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/06/03 02:13:28  wasistho
! expand CC-filtering to all dummy layers
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.4  2003/08/29 01:41:53  wasistho
! Added TARGET attribute to region variable, since pointers are cached into it
!
! Revision 1.3  2003/05/16 05:44:59  wasistho
! modified array range of CC-filtered
!
! Revision 1.2  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/14 23:55:30  wasistho
! Install Rocturb
!
!******************************************************************************







