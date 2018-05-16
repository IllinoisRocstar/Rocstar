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
! ******************************************************************************
!
! Purpose: Suite for extrapolation routines.
!
! Description: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModExtrapolation.F90,v 1.3 2008/12/06 08:44:16 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModExtrapolation

  USE ModGlobal, ONLY    : t_global 
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY      : t_grid
  USE ModBndPatch, ONLY  : t_patch
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModMPI
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLO_ExtrapRegDummyNode, &
            RFLO_ExtrapRegLastNode

! private :
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModExtrapolation.F90,v $ $Revision: 1.3 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS

!******************************************************************************
!
! Purpose: extrapolate to region dummy nodes
!
! Description: based on uniform grid linear interpolation
!
! Input: ibeg, iend, jbeg, jend, kbeg, kend = region dimension incl. dummies
!        ndum     = number of dummy points/cells
!        iNOff    = i-stride
!        ijNOff   = ij-stride
!        idb, ide = begin and end indexing of first dimension of var(:,:)
!        var      = real variable to be extrapolated
!
! Output: var = dummy points of var assigned extrapolated values
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_ExtrapRegDummyNode( ibeg,iend,jbeg,jend,kbeg,kend,ndum, &
                                    iNOff,ijNOff,idb,ide,var )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, ndum, idb, ide, iNOff, ijNOff
  REAL(RFREAL), POINTER :: var(:,:)

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: nelm, ndim, ijkN, ijkN1, ijkN2, error

!******************************************************************************

  nelm = ide-idb+1
  ndim = (iend-ibeg+1)*(jend-jbeg+1)*(kend-kbeg+1)

  IF ((SIZE( var,1 ) < nelm) .OR. (SIZE( var,2 ) /= ndim)) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR in RFLO_ExtrapRegDummyNode: '
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'inconsistent 1st or 2nd dimension'
    WRITE(STDERR,'(A)') SOLVER_NAME
#ifdef MPI
    CALL MPI_Abort( error )
#endif
    STOP
  ENDIF

  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg+ndum-1,ibeg,-1
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i+1 ,j   ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i+2 ,j   ,k   ,iNOff,ijNOff)
        var(idb:ide,ijkN) = 2*var(idb:ide,ijkN1) - var(idb:ide,ijkN2)
      ENDDO
      DO i=iend-ndum+1,iend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i-1 ,j   ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i-2 ,j   ,k   ,iNOff,ijNOff)
        var(idb:ide,ijkN) = 2*var(idb:ide,ijkN1) - var(idb:ide,ijkN2)
      ENDDO
    ENDDO
  ENDDO

  DO i=ibeg,iend
    DO k=kbeg,kend
      DO j=jbeg+ndum-1,jbeg,-1
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j+1 ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j+2 ,k   ,iNOff,ijNOff)
        var(idb:ide,ijkN) = 2*var(idb:ide,ijkN1) - var(idb:ide,ijkN2)
      ENDDO
      DO j=jend-ndum+1,jend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j-1 ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j-2 ,k   ,iNOff,ijNOff)
        var(idb:ide,ijkN) = 2*var(idb:ide,ijkN1) - var(idb:ide,ijkN2)
      ENDDO
    ENDDO
  ENDDO

  DO j=jbeg,jend
    DO i=ibeg,iend
      DO k=kbeg+ndum-1,kbeg,-1
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j   ,k+1 ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j   ,k+2 ,iNOff,ijNOff)
        var(idb:ide,ijkN) = 2*var(idb:ide,ijkN1) - var(idb:ide,ijkN2)
      ENDDO
      DO k=kend-ndum+1,kend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j   ,k-1 ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j   ,k-2 ,iNOff,ijNOff)
        var(idb:ide,ijkN) = 2*var(idb:ide,ijkN1) - var(idb:ide,ijkN2)
      ENDDO
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

END SUBROUTINE RFLO_ExtrapRegDummyNode

!******************************************************************************
!
! Purpose: extrapolate to last layer of dummy nodes
!
! Description: based on uniform grid linear interpolation
!
! Input: ibeg, iend, jbeg, jend, kbeg, kend = region dimension incl. dummies
!        ndum     = number of dummy points/cells
!        iNOff    = i-stride
!        ijNOff   = ij-stride
!        idb, ide = begin and end indexing of first dimension of var(:,:)
!        var      = real variable to be extrapolated
!
! Output: var = dummy points of var assigned extrapolated values
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_ExtrapRegLastNode( ibeg,iend,jbeg,jend,kbeg,kend, &
                                   iNOff,ijNOff,idb,ide,var )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, idb, ide, iNOff, ijNOff
  REAL(RFREAL), POINTER :: var(:,:)

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: nelm, ndim, ijkN, ijkN1, ijkN2, error

!******************************************************************************

  nelm = ide-idb+1
  ndim = (iend-ibeg+1)*(jend-jbeg+1)*(kend-kbeg+1)

  IF ((SIZE( var,1 ) < nelm) .OR. (SIZE( var,2 ) /= ndim)) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR in RFLO_ExtrapRegDummyNode: '
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'inconsistent 1st or 2nd dimension'
    WRITE(STDERR,'(A)') SOLVER_NAME
#ifdef MPI
    CALL MPI_Abort( error )
#endif
    STOP
  ENDIF

  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg,ibeg
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i+1 ,j   ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i+2 ,j   ,k   ,iNOff,ijNOff)
        var(idb:ide,ijkN) = 2*var(idb:ide,ijkN1) - var(idb:ide,ijkN2)
      ENDDO
      DO i=iend,iend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i-1 ,j   ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i-2 ,j   ,k   ,iNOff,ijNOff)
        var(idb:ide,ijkN) = 2*var(idb:ide,ijkN1) - var(idb:ide,ijkN2)
      ENDDO
    ENDDO
  ENDDO

  DO i=ibeg,iend
    DO k=kbeg,kend
      DO j=jbeg,jbeg
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j+1 ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j+2 ,k   ,iNOff,ijNOff)
        var(idb:ide,ijkN) = 2*var(idb:ide,ijkN1) - var(idb:ide,ijkN2)
      ENDDO
      DO j=jend,jend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j-1 ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j-2 ,k   ,iNOff,ijNOff)
        var(idb:ide,ijkN) = 2*var(idb:ide,ijkN1) - var(idb:ide,ijkN2)
      ENDDO
    ENDDO
  ENDDO

  DO j=jbeg,jend
    DO i=ibeg,iend
      DO k=kbeg,kbeg
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j   ,k+1 ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j   ,k+2 ,iNOff,ijNOff)
        var(idb:ide,ijkN) = 2*var(idb:ide,ijkN1) - var(idb:ide,ijkN2)
      ENDDO
      DO k=kend,kend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j   ,k-1 ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j   ,k-2 ,iNOff,ijNOff)
        var(idb:ide,ijkN) = 2*var(idb:ide,ijkN1) - var(idb:ide,ijkN2)
      ENDDO
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

END SUBROUTINE RFLO_ExtrapRegLastNode

! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModExtrapolation

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModExtrapolation.F90,v $
! Revision 1.3  2008/12/06 08:44:16  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2005/12/03 09:39:47  wasistho
! initial import
!
!
!
! ******************************************************************************






