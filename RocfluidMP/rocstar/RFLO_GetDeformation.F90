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
! Purpose: obtain deformation to boundary nodes from GenX.
!
! Description: none.
!
! Input: region     = grid dimensions and topology
!        boundMoved = flag for boundaries of region which have moved.
!
! Output: dNode = deformations at all six boundaries.
!
! Notes: variable dNode contains the whole 3-D field.
!
!******************************************************************************
!
! $Id: RFLO_GetDeformation.F90,v 1.8 2009/08/12 04:15:57 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_GetDeformation( region,boundMoved,dNode )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset, &
                            RFLO_GetDimensPhysNodes
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  LOGICAL :: boundMoved(6)

  REAL(RFREAL), POINTER :: dNode(:,:)

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iPatch, i, j, k, ng1, ng2, ijkN

! ... local variables
  INTEGER :: iLev, lbound, ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff
  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend, nfixed
  LOGICAL, POINTER :: allExternal(:)

  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch

  REAL(RFREAL) :: dextrem(1:2),dist

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_GetDeformation',&
  'RFLO_GetDeformation.F90' )

! initialize variables --------------------------------------------------------

  iLev = 1     ! finest grid
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )

! zero out displacements

  dNode(:,:) = 0._RFREAL

! reset boundary movement flags

  allExternal => region%levels(iLev)%grid%allExternal
  boundMoved(:)  = .false.
  allExternal(:) = .false. 

! obtain displacements --------------------------------------------------------

  IF ( global%myProcid == 0 .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A,1X,I3)') SOLVER_NAME, &
                                   'Getting displacements from Rocstar....'
  END IF ! global%myProcid


  DO iPatch=1,region%nPatches
    patch  => region%levels(iLev)%patches(iPatch)
    lbound =  patch%lbound
    dextrem(2) = 10.0_RFREAL
    dextrem(1) = 0._RFREAL
    nfixed = 0
    IF (patch%bcMotion == BC_EXTERNAL) THEN
       CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
            ibeg,iend,jbeg,jend,kbeg,kend )
       boundMoved(lbound) = .true.
       
       DO k=kbeg,kend
          DO j=jbeg,jend
             DO i=ibeg,iend
                ijkN = IndIJK(i,j,k,iNOff,ijNOff)
                IF      (lbound==1 .OR. lbound==2) THEN
                   IF (lbound == 2) THEN
                      ng1 = j - jbeg + 1
                   ELSE
                      ng1 = jend - j + 1
                   ENDIF
                   ng2 = k - kbeg + 1
                ELSE IF (lbound==3 .OR. lbound==4) THEN
                   ng1 = k - kbeg + 1
                   IF (lbound == 4) THEN
                      ng2 = i - ibeg + 1
                   ELSE
                      ng2 = iend - i + 1
                   ENDIF
                ELSE IF (lbound==5 .OR. lbound==6) THEN
                   IF (lbound == 6) THEN
                      ng1 = i - ibeg + 1
                   ELSE
                      ng1 = iend - i + 1
                   ENDIF
                   ng2 = j - jbeg + 1
                ENDIF
                dNode(XCOORD,ijkN) = patch%duAlp(1,ng1,ng2)
                dNode(YCOORD,ijkN) = patch%duAlp(2,ng1,ng2)
                dNode(ZCOORD,ijkN) = patch%duAlp(3,ng1,ng2)
                dist = SQRT(dNode(XCOORD,ijkN)*dNode(XCOORD,ijkN) +&
                     dNode(YCOORD,ijkN)*dNode(YCOORD,ijkN) + &
                     dNode(ZCOORD,ijkN)*dNode(ZCOORD,ijkN))
                if(dist > dextrem(2)) dextrem(2) = dist
                if(dist < dextrem(1)) dextrem(1) = dist
                if(dist < 1E-12_RFREAL) nfixed = nfixed + 1
             ENDDO
          ENDDO
       ENDDO
!       WRITE(*,*) 'MOTION EXTREMA FOR PATCH: ',dextrem(1),dextrem(2),nfixed
      IF (lbound==1 .OR. lbound==2) THEN
        IF (jbeg==jpnbeg .AND. jend==jpnend .AND. &
            kbeg==kpnbeg .AND. kend==kpnend) THEN
          allExternal(lbound) = .true.
        ENDIF
      ELSEIF (lbound==3 .OR. lbound==4) THEN
        IF (ibeg==ipnbeg .AND. iend==ipnend .AND. &
            kbeg==kpnbeg .AND. kend==kpnend) THEN
          allExternal(lbound) = .true.
        ENDIF
      ELSEIF (lbound==5 .OR. lbound==6) THEN
        IF (ibeg==ipnbeg .AND. iend==ipnend .AND. &
            jbeg==jpnbeg .AND. jend==jpnend) THEN
          allExternal(lbound) = .true.
        ENDIF
      ENDIF

    ELSE  ! external BC
 !     IF(patch%bcMotion .ne. BC_EXTERNAL) THEN
      IF (patch%mixt%setMotion .AND. &
         global%verbLevel >= VERBOSE_HIGH) THEN
         WRITE(*,*) 'SETTING MOTION ON PATCH '
        CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                        ibeg,iend,jbeg,jend,kbeg,kend )
        boundMoved(lbound) = .true.
      
        DO k=kbeg,kend
          DO j=jbeg,jend
            DO i=ibeg,iend
              ijkN = IndIJK(i,j,k,iNOff,ijNOff)
              if(dNode(XCOORD,ijkN) .eq. 0) then
                 dNode(XCOORD,ijkN) = patch%mixt%bndVel(XCOORD)*global%dtMin
              endif
              if(dNode(YCOORD,ijkN) .eq. 0) then
                 dNode(YCOORD,ijkN) = patch%mixt%bndVel(YCOORD)*global%dtMin
              endif
              if(dNode(XCOORD,ijkN) .eq. 0) then
                 dNode(ZCOORD,ijkN) = patch%mixt%bndVel(ZCOORD)*global%dtMin
              endif
            ENDDO
          ENDDO
        ENDDO

      ENDIF  ! motion prescribed internally
    ENDIF    ! external BC
  ENDDO      ! iPatch

! finalize --------------------------------------------------------------------
  IF ( global%myProcid == 0 .AND. & 
       global%verbLevel >= VERBOSE_HIGH ) THEN
    WRITE(STDOUT,'(A,1X,A,1X,I3)') SOLVER_NAME, &
                                   'Getting displacements from Rocstar done.'
  END IF ! global%myProcid

999  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_GetDeformation

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_GetDeformation.F90,v $
! Revision 1.8  2009/08/12 04:15:57  mtcampbe
! Major update, bugfix from Abe development, more propagation compatibility,
! some Rocstar IO changes, Ju's temporary clipping fix for turbulence. A bug
! fix for initialization IO.
!
! Revision 1.7  2008/12/06 08:44:00  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/08/19 15:37:49  mparmar
! Renamed patch variables
!
! Revision 1.4  2005/09/28 20:49:26  wasistho
! modified for non-external moving boundary
!
! Revision 1.3  2005/06/13 21:45:33  wasistho
! changed patch%bcCoupled to patch%bcMotion
!
! Revision 1.2  2005/06/05 23:03:02  wasistho
! distinguish external boundary to be fully and partly external
!
! Revision 1.1  2004/12/01 21:23:50  haselbac
! Initial revision after changing case
!
! Revision 1.8  2004/11/13 22:37:41  wasistho
! invert orientation of genx-surface-variables
!
! Revision 1.7  2003/11/20 16:40:33  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.4  2003/05/15 02:57:00  jblazek
! Inlined index function.
!
! Revision 1.3  2002/10/10 23:49:48  jblazek
! Changed orientation of surface grid.
!
! Revision 1.2  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.1  2002/09/20 22:22:34  jblazek
! Finalized integration into GenX.
!
!******************************************************************************







