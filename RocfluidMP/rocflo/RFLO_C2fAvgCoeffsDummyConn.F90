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
! Purpose: Extrapolate coefficient from the interior domain or at the patch 
!          to the dummy points at connecting boundaries.
!
! Description: Dummy values of connecting bnd has been computed in prev. step 
!              up to second outerst layer. They are copied to the
!              outerst layer in this routine.
!
! Input: region   = current region data
!        lbound   = patch boundary ID
!        i,j,kdir = direction identifier
!        i,j,kndBeg = i,j,k begin index
!        i,j,kndEnd = i,j,k end index
!  
! Output: Averaging coefficients c2fCoI, c2fCoJ, c2fCoK at outerst layer of
!         conn. bc.
!
! Notes: Mother routine = RFLO_C2fAvgCoeffsDummy.
!
!******************************************************************************
!
! $Id: RFLO_C2fAvgCoeffsDummyConn.F90,v 1.3 2008/12/06 08:44:25 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_C2fAvgCoeffsDummyConn( region,lbound,idir,jdir,kdir, &
                                 indBeg,indEnd,jndBeg,jndEnd,kndBeg,kndEnd )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: lbound,idir,jdir,kdir
  INTEGER        :: indBeg, indEnd, jndBeg, jndEnd, kndBeg, kndEnd

! ... loop variables
  INTEGER :: i, j, k, idum, jdum, kdum

! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: iLev, iNOff, ijNOff
  INTEGER :: ijkNI, ijkNJ, ijkNK, ijkND
  INTEGER :: iG(3), jG(3), kG(3)
  INTEGER :: iDumB(3), jDumB(3), kDumB(3), iDumE(3), jDumE(3), kDumE(3)

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RFLO_C2fAvgCoeffsDummyConn',&
  'RFLO_C2fAvgCoeffsDummyConn.F90' )

! get dimensions -------------------------------------------------------------

  iLev = region%currLevel
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  iDumB(:) = -idir*region%nDumCells
  jDumB(:) = -jdir*region%nDumCells
  kDumB(:) = -kdir*region%nDumCells
  iDumE(:) = iDumB(:)
  jDumE(:) = jDumB(:)
  kDumE(:) = kDumB(:)

  iG(:) = -idir*(region%nDumCells-1)
  jG(:) = -jdir*(region%nDumCells-1)
  kG(:) = -kdir*(region%nDumCells-1)

! 2D loop over patch nodes ----------------------------------------------------

  DO k=kndBeg,kndEnd
    DO j=jndBeg,jndEnd
      DO i=indBeg,indEnd

        ijkNI = IndIJK(i+iG(1),j+jG(1),k+kG(1),iNOff,ijNOff)  ! i reference
        ijkNJ = IndIJK(i+iG(2),j+jG(2),k+kG(2),iNOff,ijNOff)  ! j reference
        ijkNK = IndIJK(i+iG(3),j+jG(3),k+kG(3),iNOff,ijNOff)  ! k reference

! ----- 1D loop in direction normal to patch to define grads at dummy faces

        DO idum=iDumB(1),iDumE(1)
          DO jdum=jDumB(1),jDumE(1)
            DO kdum=kDumB(1),kDumE(1) 
              ijkND = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff)
              region%levels(iLev)%grid%c2fCoI(2,ijkND) = &
              region%levels(iLev)%grid%c2fCoI(2,ijkNI)
              region%levels(iLev)%grid%c2fCoI(1,ijkND) = &
              1._RFREAL-region%levels(iLev)%grid%c2fCoI(2,ijkND)
            ENDDO
          ENDDO
        ENDDO

        DO idum=iDumB(2),iDumE(2)
          DO jdum=jDumB(2),jDumE(2)
            DO kdum=kDumB(2),kDumE(2) 
              ijkND  = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff)
              region%levels(iLev)%grid%c2fCoJ(2,ijkND) = &
              region%levels(iLev)%grid%c2fCoJ(2,ijkNJ)
              region%levels(iLev)%grid%c2fCoJ(1,ijkND) = &
              1._RFREAL-region%levels(iLev)%grid%c2fCoJ(2,ijkND)
            ENDDO
          ENDDO
        ENDDO

        DO idum=iDumB(3),iDumE(3)
          DO jdum=jDumB(3),jDumE(3)
            DO kdum=kDumB(3),kDumE(3) 
              ijkND  = IndIJK(i+idum,j+jdum,k+kdum,iNOff,ijNOff)
              region%levels(iLev)%grid%c2fCoK(2,ijkND) = &
              region%levels(iLev)%grid%c2fCoK(2,ijkNK)
              region%levels(iLev)%grid%c2fCoK(1,ijkND) = &
              1._RFREAL-region%levels(iLev)%grid%c2fCoK(2,ijkND)
            ENDDO
          ENDDO
        ENDDO

      ENDDO   ! i
    ENDDO     ! j
  ENDDO       ! k

! copy coefficients from region boundary to patch edge

  CALL AvgCoPatchEdge

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

! ==============================================================================
! subroutine for copying to patch edges 
! ==============================================================================

CONTAINS
 
  SUBROUTINE AvgCoPatchEdge

! ... local variables

  INTEGER :: inbeg, inend, jnbeg, jnend, knbeg, knend
  INTEGER :: ijkN, ijkN1, iref, jref, kref
  REAL(RFREAL), POINTER :: avgCo(:,:)
      
  IF (lbound==1 .OR. lbound==3 .OR. lbound==5) THEN
    inbeg = indBeg-idir*region%nDumCells
    jnbeg = jndBeg-jdir*region%nDumCells
    knbeg = kndBeg-kdir*region%nDumCells
    inend = indEnd-idir
    jnend = jndEnd-jdir
    knend = kndEnd-kdir
  ELSE
    inbeg = indBeg
    jnbeg = jndBeg
    knbeg = kndBeg
    inend = indEnd-idir*region%nDumCells
    jnend = jndEnd-jdir*region%nDumCells
    knend = kndEnd-kdir*region%nDumCells
  ENDIF

  IF ((lbound==1).OR.(lbound==2)) THEN
    avgCo => region%levels(iLev)%grid%c2fCoI
    IF (lbound==1) THEN
      iref = indEnd
    ELSE
      iref = indBeg-1
    ENDIF 
    DO i=inbeg,inend
      DO k=knbeg,knend-1
        ijkN    = IndIJK(i    ,jnbeg  ,k,iNOff,ijNOff)
        ijkN1   = IndIJK(iref ,jnbeg  ,k,iNOff,ijNOff)
        IF (avgCo(1,ijkN)==0._RFREAL .OR. avgCo(2,ijkN)==0._RFREAL) THEN
          avgCo(:,ijkN)=avgCo(:,ijkN1)
        ENDIF
        ijkN    = IndIJK(i    ,jnend  ,k,iNOff,ijNOff)
        ijkN1   = IndIJK(iref ,jnend  ,k,iNOff,ijNOff)
        IF (avgCo(1,ijkN)==0._RFREAL .OR. avgCo(2,ijkN)==0._RFREAL) THEN
          avgCo(:,ijkN)=avgCo(:,ijkN1)
        ENDIF
      END DO  
      DO j=jnbeg,jnend-1
        ijkN    = IndIJK(i    ,j,knbeg  ,iNOff,ijNOff)
        ijkN1   = IndIJK(iref ,j,knbeg  ,iNOff,ijNOff)
        IF (avgCo(1,ijkN)==0._RFREAL .OR. avgCo(2,ijkN)==0._RFREAL) THEN
          avgCo(:,ijkN)=avgCo(:,ijkN1)
        ENDIF
        ijkN    = IndIJK(i    ,j,knend  ,iNOff,ijNOff)
        ijkN1   = IndIJK(iref ,j,knend  ,iNOff,ijNOff)
        IF (avgCo(1,ijkN)==0._RFREAL .OR. avgCo(2,ijkN)==0._RFREAL) THEN
          avgCo(:,ijkN)=avgCo(:,ijkN1)
        ENDIF
      END DO  
    END DO
  ELSEIF ((lbound==3).OR.(lbound==4)) THEN
    avgCo => region%levels(iLev)%grid%c2fCoJ
    IF (lbound==3) THEN
      jref = jndEnd
    ELSE
      jref = jndBeg-1
    ENDIF 
    DO j=jnbeg,jnend
      DO k=knbeg,knend-1
        ijkN    = IndIJK(inbeg ,j    ,k,iNOff,ijNOff)
        ijkN1   = IndIJK(inbeg ,jref ,k,iNOff,ijNOff)
        IF (avgCo(1,ijkN)==0._RFREAL .OR. avgCo(2,ijkN)==0._RFREAL) THEN
          avgCo(:,ijkN)=avgCo(:,ijkN1)
        ENDIF
        ijkN    = IndIJK(inend ,j    ,k,iNOff,ijNOff)
        ijkN1   = IndIJK(inend ,jref ,k,iNOff,ijNOff)
        IF (avgCo(1,ijkN)==0._RFREAL .OR. avgCo(2,ijkN)==0._RFREAL) THEN
          avgCo(:,ijkN)=avgCo(:,ijkN1)
        ENDIF
      END DO  
      DO i=inbeg,inend-1
        ijkN    = IndIJK(i,j    ,knbeg ,iNOff,ijNOff)
        ijkN1   = IndIJK(i,jref ,knbeg ,iNOff,ijNOff)
        IF (avgCo(1,ijkN)==0._RFREAL .OR. avgCo(2,ijkN)==0._RFREAL) THEN
          avgCo(:,ijkN)=avgCo(:,ijkN1)
        ENDIF
        ijkN    = IndIJK(i,j    ,knend ,iNOff,ijNOff)
        ijkN1   = IndIJK(i,jref ,knend ,iNOff,ijNOff)
        IF (avgCo(1,ijkN)==0._RFREAL .OR. avgCo(2,ijkN)==0._RFREAL) THEN
          avgCo(:,ijkN)=avgCo(:,ijkN1)
        ENDIF
      END DO  
    END DO  
  ELSEIF ((lbound==5).OR.(lbound==6)) THEN
    avgCo => region%levels(iLev)%grid%c2fCoK
    IF (lbound==5) THEN
      kref = kndEnd
    ELSE
      kref = kndBeg-1
    ENDIF 
    DO k=knbeg,knend
      DO i=inbeg,inend-1
        ijkN    = IndIJK(i,jnbeg ,k    ,iNOff,ijNOff)
        ijkN1   = IndIJK(i,jnbeg ,kref ,iNOff,ijNOff)
        IF (avgCo(1,ijkN)==0._RFREAL .OR. avgCo(2,ijkN)==0._RFREAL) THEN
          avgCo(:,ijkN)=avgCo(:,ijkN1)
        ENDIF
        ijkN    = IndIJK(i,jnend ,k    ,iNOff,ijNOff)
        ijkN1   = IndIJK(i,jnend ,kref ,iNOff,ijNOff)
        IF (avgCo(1,ijkN)==0._RFREAL .OR. avgCo(2,ijkN)==0._RFREAL) THEN
          avgCo(:,ijkN)=avgCo(:,ijkN1)
        ENDIF
      END DO  
      DO j=jnbeg,jnend-1
        ijkN    = IndIJK(inbeg ,j,k    ,iNOff,ijNOff)
        ijkN1   = IndIJK(inbeg ,j,kref ,iNOff,ijNOff)
        IF (avgCo(1,ijkN)==0._RFREAL .OR. avgCo(2,ijkN)==0._RFREAL) THEN
          avgCo(:,ijkN)=avgCo(:,ijkN1)
        ENDIF
        ijkN    = IndIJK(inend ,j,k    ,iNOff,ijNOff)
        ijkN1   = IndIJK(inend ,j,kref ,iNOff,ijNOff)
        IF (avgCo(1,ijkN)==0._RFREAL .OR. avgCo(2,ijkN)==0._RFREAL) THEN
          avgCo(:,ijkN)=avgCo(:,ijkN1)
        ENDIF
      END DO  
    END DO  
  END IF

  END SUBROUTINE AvgCoPatchEdge

END SUBROUTINE RFLO_C2fAvgCoeffsDummyConn

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: RFLO_C2fAvgCoeffsDummyConn.F90,v $
! Revision 1.3  2008/12/06 08:44:25  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:36  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.3  2004/08/03 00:52:08  wasistho
! changed avgCo to c2fCo in the description
!
! Revision 1.2  2004/08/02 19:32:43  wasistho
! changed grid%avgCo to grid%c2fCo
!
! Revision 1.1  2004/07/30 17:30:31  wasistho
! initial import routines starting with RFLO_c2fAvg...
!
!
!
!******************************************************************************







