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
! Purpose: Suite for finite difference routines.
!
! Description: None.
!
! ******************************************************************************
!
! $Id: RFLO_ModFiniteDifference.F90,v 1.4 2008/12/06 08:44:16 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModFiniteDifference

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
  PUBLIC :: RFLO_FinDiffCompI, &
            RFLO_FinDiffCompJ, &
            RFLO_FinDiffCompK, &
            RFLO_FinDiffCompII, &
            RFLO_FinDiffCompJJ, &
            RFLO_FinDiffCompKK, &
            RFLO_FinDiffCompIs, &
            RFLO_FinDiffCompJs, &
            RFLO_FinDiffCompIIs, &
            RFLO_FinDiffCompJJs

! private :
    
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLO_ModFiniteDifference.F90,v $ $Revision: 1.4 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS

!******************************************************************************
!
! Purpose: perform 1st derivative finite difference in I-direction in
!          computational space
!
! Description: based on 2nd order central differencing
!
! Input: ibeg, iend, jbeg, jend, kbeg, kend = region range indices
!        iNOff    = i-stride
!        ijNOff   = ij-stride
!        idb, ide = begin and end indexing of first dimension of var(:,:)
!        var      = real variable to be extrapolated
!
! Output: dvar = resulting derivative
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_FinDiffCompI( ibeg,iend,jbeg,jend,kbeg,kend,ndum, &
                              iNOff,ijNOff,idb,ide,var,dvar )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, ndum, idb, ide, iNOff, ijNOff
  REAL(RFREAL), POINTER :: var(:,:), dvar(:,:)

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: nelm, ndim, ijkN, ijkNp, ijkNm, ijkN1, ijkN2, error
  REAL(RFREAL) :: dh, r2dh

!******************************************************************************

  nelm = ide-idb+1
  ndim = (iend-ibeg+1)*(jend-jbeg+1)*(kend-kbeg+1)

  IF ((SIZE( var,1 ) < nelm) .OR. (SIZE( var,2 ) /= ndim)) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR in RFLO_FinDiffCompI: '
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'inconsistent 1st or 2nd dimension'
    WRITE(STDERR,'(A)') SOLVER_NAME
#ifdef MPI
    CALL MPI_Abort( error )
#endif
    STOP
  ENDIF

  dh   = 1._RFREAL/REAL( iend-ibeg-2*ndum )
  r2dh = 0.5_RFREAL/dh

  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg+1,iend-1
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkNp = IndIJK(i+1 ,j   ,k   ,iNOff,ijNOff)
        ijkNm = IndIJK(i-1 ,j   ,k   ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = (var(idb:ide,ijkNp) - var(idb:ide,ijkNm))*r2dh
      ENDDO
      DO i=ibeg,ibeg
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i+1 ,j   ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i+2 ,j   ,k   ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = 2*dvar(idb:ide,ijkN1) - dvar(idb:ide,ijkN2)
      ENDDO
      DO i=iend,iend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i-1 ,j   ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i-2 ,j   ,k   ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = 2*dvar(idb:ide,ijkN1) - dvar(idb:ide,ijkN2)
      ENDDO
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

END SUBROUTINE RFLO_FinDiffCompI


!******************************************************************************
!
! Purpose: perform 1st derivative finite difference in J-direction in 
!          computational space
!
! Description: based on 2nd order central differencing
!
! Input: ibeg, iend, jbeg, jend, kbeg, kend = region range indices
!        iNOff    = i-stride
!        ijNOff   = ij-stride
!        idb, ide = begin and end indexing of first dimension of var(:,:)
!        var      = real variable to be extrapolated
!
! Output: dvar = resulting derivative
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_FinDiffCompJ( ibeg,iend,jbeg,jend,kbeg,kend,ndum, &
                              iNOff,ijNOff,idb,ide,var,dvar )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, ndum, idb, ide, iNOff, ijNOff
  REAL(RFREAL), POINTER :: var(:,:), dvar(:,:)

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: nelm, ndim, ijkN, ijkNp, ijkNm, ijkN1, ijkN2, error
  REAL(RFREAL) :: dh, r2dh

!******************************************************************************

  nelm = ide-idb+1
  ndim = (iend-ibeg+1)*(jend-jbeg+1)*(kend-kbeg+1)

  IF ((SIZE( var,1 ) < nelm) .OR. (SIZE( var,2 ) /= ndim)) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR in RFLO_FinDiffCompJ: '
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'inconsistent 1st or 2nd dimension'
    WRITE(STDERR,'(A)') SOLVER_NAME
#ifdef MPI
    CALL MPI_Abort( error )
#endif
    STOP
  ENDIF

  dh   = 1._RFREAL/REAL( jend-jbeg-2*ndum )
  r2dh = 0.5_RFREAL/dh

  DO k=kbeg,kend
    DO j=jbeg+1,jend-1
      DO i=ibeg,iend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkNp = IndIJK(i   ,j+1 ,k   ,iNOff,ijNOff)
        ijkNm = IndIJK(i   ,j-1 ,k   ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = (var(idb:ide,ijkNp) - var(idb:ide,ijkNm))*r2dh
      ENDDO
    ENDDO
  ENDDO

  DO i=ibeg,iend
    DO k=kbeg,kend
      DO j=jbeg,jbeg
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j+1 ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j+2 ,k   ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = 2*dvar(idb:ide,ijkN1) - dvar(idb:ide,ijkN2)
      ENDDO
      DO j=jend,jend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j-1 ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j-2 ,k   ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = 2*dvar(idb:ide,ijkN1) - dvar(idb:ide,ijkN2)
      ENDDO
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

END SUBROUTINE RFLO_FinDiffCompJ


!******************************************************************************
!
! Purpose: perform 1st derivative finite difference in K-direction in 
!          computational space
!
! Description: based on 2nd order central differencing
!
! Input: ibeg, iend, jbeg, jend, kbeg, kend = region range indices
!        iNOff    = i-stride
!        ijNOff   = ij-stride
!        idb, ide = begin and end indexing of first dimension of var(:,:)
!        var      = real variable to be extrapolated
!
! Output: dvar = resulting derivative
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_FinDiffCompK( ibeg,iend,jbeg,jend,kbeg,kend,ndum, &
                              iNOff,ijNOff,idb,ide,var,dvar )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, ndum, idb, ide, iNOff, ijNOff
  REAL(RFREAL), POINTER :: var(:,:), dvar(:,:)

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: nelm, ndim, ijkN, ijkNp, ijkNm, ijkN1, ijkN2, error
  REAL(RFREAL) :: dh, r2dh

!******************************************************************************

  nelm = ide-idb+1
  ndim = (iend-ibeg+1)*(jend-jbeg+1)*(kend-kbeg+1)

  IF ((SIZE( var,1 ) < nelm) .OR. (SIZE( var,2 ) /= ndim)) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR in RFLO_FinDiffCompK: '
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'inconsistent 1st or 2nd dimension'
    WRITE(STDERR,'(A)') SOLVER_NAME
#ifdef MPI
    CALL MPI_Abort( error )
#endif
    STOP
  ENDIF

  dh   = 1._RFREAL/REAL( kend-kbeg-2*ndum )
  r2dh = 0.5_RFREAL/dh

  DO k=kbeg+1,kend-1
    DO j=jbeg,jend
      DO i=ibeg,iend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkNp = IndIJK(i   ,j   ,k+1 ,iNOff,ijNOff)
        ijkNm = IndIJK(i   ,j   ,k-1 ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = (var(idb:ide,ijkNp) - var(idb:ide,ijkNm))*r2dh
      ENDDO
    ENDDO
  ENDDO

  DO i=ibeg,iend
    DO j=jbeg,jend
      DO k=kbeg,kbeg
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j   ,k+1 ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j   ,k+2 ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = 2*dvar(idb:ide,ijkN1) - dvar(idb:ide,ijkN2)
      ENDDO
      DO k=kend,kend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j   ,k-1 ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j   ,k-2 ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = 2*dvar(idb:ide,ijkN1) - dvar(idb:ide,ijkN2)
      ENDDO
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

END SUBROUTINE RFLO_FinDiffCompK


!******************************************************************************
!
! Purpose: perform 2nd derivative finite difference in I-direction in
!          computational space
!
! Description: based on 2nd order central differencing
!
! Input: ibeg, iend, jbeg, jend, kbeg, kend = region range indices
!        iNOff    = i-stride
!        ijNOff   = ij-stride
!        idb, ide = begin and end indexing of first dimension of var(:,:)
!        var      = real variable to be extrapolated
!
! Output: dvar = resulting derivative
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_FinDiffCompII( ibeg,iend,jbeg,jend,kbeg,kend,ndum, &
                               iNOff,ijNOff,idb,ide,var,dvar )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, ndum, idb, ide, iNOff, ijNOff
  REAL(RFREAL), POINTER :: var(:,:), dvar(:,:)

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: nelm, ndim, ijkN, ijkNp, ijkNm, ijkN1, ijkN2, error
  REAL(RFREAL) :: dh, rdh2

!******************************************************************************

  nelm = ide-idb+1
  ndim = (iend-ibeg+1)*(jend-jbeg+1)*(kend-kbeg+1)

  IF ((SIZE( var,1 ) < nelm) .OR. (SIZE( var,2 ) /= ndim)) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR in RFLO_FinDiffCompII: '
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'inconsistent 1st or 2nd dimension'
    WRITE(STDERR,'(A)') SOLVER_NAME
#ifdef MPI
    CALL MPI_Abort( error )
#endif
    STOP
  ENDIF

  dh   = 1._RFREAL/REAL( iend-ibeg-2*ndum )
  rdh2 = 1._RFREAL/(dh*dh)

  DO k=kbeg,kend
    DO j=jbeg,jend
      DO i=ibeg+1,iend-1
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkNp = IndIJK(i+1 ,j   ,k   ,iNOff,ijNOff)
        ijkNm = IndIJK(i-1 ,j   ,k   ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = (var(idb:ide,ijkNp) -2._RFREAL*var(idb:ide,ijkN)+ &
                              var(idb:ide,ijkNm))*rdh2
      ENDDO
      DO i=ibeg,ibeg
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i+1 ,j   ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i+2 ,j   ,k   ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = 2*dvar(idb:ide,ijkN1) - dvar(idb:ide,ijkN2)
      ENDDO
      DO i=iend,iend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i-1 ,j   ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i-2 ,j   ,k   ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = 2*dvar(idb:ide,ijkN1) - dvar(idb:ide,ijkN2)
      ENDDO
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

END SUBROUTINE RFLO_FinDiffCompII


!******************************************************************************
!
! Purpose: perform 2nd derivative finite difference in J-direction in 
!          computational space
!
! Description: based on 2nd order central differencing
!
! Input: ibeg, iend, jbeg, jend, kbeg, kend = region range indices
!        iNOff    = i-stride
!        ijNOff   = ij-stride
!        idb, ide = begin and end indexing of first dimension of var(:,:)
!        var      = real variable to be extrapolated
!
! Output: dvar = resulting derivative
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_FinDiffCompJJ( ibeg,iend,jbeg,jend,kbeg,kend,ndum, &
                               iNOff,ijNOff,idb,ide,var,dvar )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, ndum, idb, ide, iNOff, ijNOff
  REAL(RFREAL), POINTER :: var(:,:), dvar(:,:)

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: nelm, ndim, ijkN, ijkNp, ijkNm, ijkN1, ijkN2, error
  REAL(RFREAL) :: dh, rdh2

!******************************************************************************

  nelm = ide-idb+1
  ndim = (iend-ibeg+1)*(jend-jbeg+1)*(kend-kbeg+1)

  IF ((SIZE( var,1 ) < nelm) .OR. (SIZE( var,2 ) /= ndim)) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR in RFLO_FinDiffCompJJ: '
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'inconsistent 1st or 2nd dimension'
    WRITE(STDERR,'(A)') SOLVER_NAME
#ifdef MPI
    CALL MPI_Abort( error )
#endif
    STOP
  ENDIF

  dh   = 1._RFREAL/REAL( jend-jbeg-2*ndum )
  rdh2 = 1._RFREAL/(dh*dh)

  DO k=kbeg,kend
    DO j=jbeg+1,jend-1
      DO i=ibeg,iend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkNp = IndIJK(i   ,j+1 ,k   ,iNOff,ijNOff)
        ijkNm = IndIJK(i   ,j-1 ,k   ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = (var(idb:ide,ijkNp)-2._RFREAL*var(idb:ide,ijkN)+ & 
                              var(idb:ide,ijkNm))*rdh2
      ENDDO
    ENDDO
  ENDDO

  DO i=ibeg,iend
    DO k=kbeg,kend
      DO j=jbeg,jbeg
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j+1 ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j+2 ,k   ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = 2*dvar(idb:ide,ijkN1) - dvar(idb:ide,ijkN2)
      ENDDO
      DO j=jend,jend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j-1 ,k   ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j-2 ,k   ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = 2*dvar(idb:ide,ijkN1) - dvar(idb:ide,ijkN2)
      ENDDO
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

END SUBROUTINE RFLO_FinDiffCompJJ


!******************************************************************************
!
! Purpose: perform 2nd derivative finite difference in K-direction in 
!          computational space
!
! Description: based on 2nd order central differencing
!
! Input: ibeg, iend, jbeg, jend, kbeg, kend = region range indices
!        iNOff    = i-stride
!        ijNOff   = ij-stride
!        idb, ide = begin and end indexing of first dimension of var(:,:)
!        var      = real variable to be extrapolated
!
! Output: dvar = resulting derivative
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_FinDiffCompKK( ibeg,iend,jbeg,jend,kbeg,kend,ndum, &
                               iNOff,ijNOff,idb,ide,var,dvar )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, ndum, idb, ide, iNOff, ijNOff
  REAL(RFREAL), POINTER :: var(:,:), dvar(:,:)

! ... loop variables
  INTEGER :: i, j, k

! ... local variables
  INTEGER :: nelm, ndim, ijkN, ijkNp, ijkNm, ijkN1, ijkN2, error
  REAL(RFREAL) :: dh, rdh2

!******************************************************************************

  nelm = ide-idb+1
  ndim = (iend-ibeg+1)*(jend-jbeg+1)*(kend-kbeg+1)

  IF ((SIZE( var,1 ) < nelm) .OR. (SIZE( var,2 ) /= ndim)) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR in RFLO_FinDiffCompKK: '
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'inconsistent 1st or 2nd dimension'
    WRITE(STDERR,'(A)') SOLVER_NAME
#ifdef MPI
    CALL MPI_Abort( error )
#endif
    STOP
  ENDIF

  dh   = 1._RFREAL/REAL( kend-kbeg-2*ndum )
  rdh2 = 1._RFREAL/(dh*dh)

  DO k=kbeg+1,kend-1
    DO j=jbeg,jend
      DO i=ibeg,iend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkNp = IndIJK(i   ,j   ,k+1 ,iNOff,ijNOff)
        ijkNm = IndIJK(i   ,j   ,k-1 ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = (var(idb:ide,ijkNp)-2._RFREAL*var(idb:ide,ijkN)+ & 
                              var(idb:ide,ijkNm))*rdh2
      ENDDO
    ENDDO
  ENDDO

  DO i=ibeg,iend
    DO j=jbeg,jend
      DO k=kbeg,kbeg
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j   ,k+1 ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j   ,k+2 ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = 2*dvar(idb:ide,ijkN1) - dvar(idb:ide,ijkN2)
      ENDDO
      DO k=kend,kend
        ijkN  = IndIJK(i   ,j   ,k   ,iNOff,ijNOff)
        ijkN1 = IndIJK(i   ,j   ,k-1 ,iNOff,ijNOff)
        ijkN2 = IndIJK(i   ,j   ,k-2 ,iNOff,ijNOff)
        dvar(idb:ide,ijkN) = 2*dvar(idb:ide,ijkN1) - dvar(idb:ide,ijkN2)
      ENDDO
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

END SUBROUTINE RFLO_FinDiffCompKK


!******************************************************************************
!
! Purpose: perform 1st derivative finite difference in I-direction in
!          2D computational space (patch surface)
!
! Description: based on 2nd order central differencing
!
! Input: ni, nj   = patch dimensions
!        idb, ide = begin and end indexing of first dimension of var(:,:)
!        var      = real variable to be extrapolated
!
! Output: dvar = resulting derivative
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_FinDiffCompIs( ni,nj,idb,ide,var,dvar )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: ni, nj, idb, ide
  REAL(RFREAL), POINTER :: var(:,:,:), dvar(:,:,:)

! ... loop variables
  INTEGER :: i, j

! ... local variables
  INTEGER :: nelm, ndim1, ndim2, error
  REAL(RFREAL) :: dh, r2dh

!******************************************************************************

  nelm  = ide-idb+1
  ndim1 = ni
  ndim2 = nj

  IF ((SIZE( var,1 ) <  nelm ) .OR. &
      (SIZE( var,2 ) /= ndim1) .OR. &
      (SIZE( var,3 ) /= ndim2)) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR in RFLO_FinDiffCompIs: '
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'inconsistent 1st or 2nd dimension'
    WRITE(STDERR,'(A)') SOLVER_NAME
#ifdef MPI
    CALL MPI_Abort( error )
#endif
    STOP
  ENDIF

  dh   = 1._RFREAL/REAL( ni-1 )
  r2dh = 0.5_RFREAL/dh

  DO j=1,nj
    DO i=2,ni-1
      dvar(idb:ide,i,j) = (var(idb:ide,i+1,j) - var(idb:ide,i-1,j))*r2dh
    ENDDO
    DO i=1,1
      dvar(idb:ide,i,j) = 2*dvar(idb:ide,i+1,j) - dvar(idb:ide,i+2,j)
    ENDDO
    DO i=ni,ni
      dvar(idb:ide,i,j) = 2*dvar(idb:ide,i-1,j) - dvar(idb:ide,i-2,j)
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

END SUBROUTINE RFLO_FinDiffCompIs


!******************************************************************************
!
! Purpose: perform 1st derivative finite difference in J-direction in
!          2D computational space (patch surface)
!
! Description: based on 2nd order central differencing
!
! Input: ni, nj   = patch dimensions
!        idb, ide = begin and end indexing of first dimension of var(:,:)
!        var      = real variable to be extrapolated
!
! Output: dvar = resulting derivative
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_FinDiffCompJs( ni,nj,idb,ide,var,dvar )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: ni, nj, idb, ide
  REAL(RFREAL), POINTER :: var(:,:,:), dvar(:,:,:)

! ... loop variables
  INTEGER :: i, j

! ... local variables
  INTEGER :: nelm, ndim1, ndim2, error
  REAL(RFREAL) :: dh, r2dh

!******************************************************************************

  nelm  = ide-idb+1
  ndim1 = ni
  ndim2 = nj

  IF ((SIZE( var,1 ) <  nelm ) .OR. &
      (SIZE( var,2 ) /= ndim1) .OR. &
      (SIZE( var,3 ) /= ndim2)) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR in RFLO_FinDiffCompJs: '
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'inconsistent 1st or 2nd dimension'
    WRITE(STDERR,'(A)') SOLVER_NAME
#ifdef MPI
    CALL MPI_Abort( error )
#endif
    STOP
  ENDIF

  dh   = 1._RFREAL/REAL( nj-1 )
  r2dh = 0.5_RFREAL/dh

  DO j=2,nj-1
    DO i=1,ni
      dvar(idb:ide,i,j) = (var(idb:ide,i,j+1) - var(idb:ide,i,j-1))*r2dh
    ENDDO
  ENDDO
  DO j=1,1
    DO i=1,ni
      dvar(idb:ide,i,j) = 2*dvar(idb:ide,i,j+1) - dvar(idb:ide,i,j+2)
    ENDDO
  ENDDO
  DO j=nj,nj
    DO i=1,ni
      dvar(idb:ide,i,j) = 2*dvar(idb:ide,i,j-1) - dvar(idb:ide,i,j-2)
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

END SUBROUTINE RFLO_FinDiffCompJs


!******************************************************************************
!
! Purpose: perform 2nd derivative finite difference in I-direction in
!          2D computational space (patch surface)
!
! Description: based on 2nd order central differencing
!
! Input: ni, nj   = patch dimensions
!        idb, ide = begin and end indexing of first dimension of var(:,:)
!        var      = real variable to be extrapolated
!
! Output: dvar = resulting derivative
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_FinDiffCompIIs( ni,nj,idb,ide,var,dvar )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: ni, nj, idb, ide
  REAL(RFREAL), POINTER :: var(:,:,:), dvar(:,:,:)

! ... loop variables
  INTEGER :: i, j

! ... local variables
  INTEGER :: nelm, ndim1, ndim2, error
  REAL(RFREAL) :: dh, rdh2

!******************************************************************************

  nelm  = ide-idb+1
  ndim1 = ni
  ndim2 = nj

  IF ((SIZE( var,1 ) <  nelm ) .OR. &
      (SIZE( var,2 ) /= ndim1) .OR. &
      (SIZE( var,3 ) /= ndim2)) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR in RFLO_FinDiffCompIIs: '
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'inconsistent 1st or 2nd dimension'
    WRITE(STDERR,'(A)') SOLVER_NAME
#ifdef MPI
    CALL MPI_Abort( error )
#endif
    STOP
  ENDIF

  dh   = 1._RFREAL/REAL( ni-1 )
  rdh2 = 1._RFREAL/(dh*dh)

  DO j=1,nj
    DO i=2,ni-1
      dvar(idb:ide,i,j) = (var(idb:ide,i+1,j) - 2._RFREAL*var(idb:ide,i,j) + &
                           var(idb:ide,i-1,j))*rdh2
    ENDDO
    DO i=1,1
      dvar(idb:ide,i,j) = 2*dvar(idb:ide,i+1,j) - dvar(idb:ide,i+2,j)
    ENDDO
    DO i=ni,ni
      dvar(idb:ide,i,j) = 2*dvar(idb:ide,i-1,j) - dvar(idb:ide,i-2,j)
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

END SUBROUTINE RFLO_FinDiffCompIIs


!******************************************************************************
!
! Purpose: perform 2nd derivative finite difference in J-direction in
!          2D computational space (patch surface)
!
! Description: based on 2nd order central differencing
!
! Input: ni, nj   = patch dimensions
!        idb, ide = begin and end indexing of first dimension of var(:,:)
!        var      = real variable to be extrapolated
!
! Output: dvar = resulting derivative
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE RFLO_FinDiffCompJJs( ni,nj,idb,ide,var,dvar )

  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: ni, nj, idb, ide
  REAL(RFREAL), POINTER :: var(:,:,:), dvar(:,:,:)

! ... loop variables
  INTEGER :: i, j

! ... local variables
  INTEGER :: nelm, ndim1, ndim2, error
  REAL(RFREAL) :: dh, rdh2

!******************************************************************************

  nelm  = ide-idb+1
  ndim1 = ni
  ndim2 = nj

  IF ((SIZE( var,1 ) <  nelm ) .OR. &
      (SIZE( var,2 ) /= ndim1) .OR. &
      (SIZE( var,3 ) /= ndim2)) THEN 
    WRITE(STDERR,'(A)') SOLVER_NAME
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'ERROR in RFLO_FinDiffCompJJs: '
    WRITE(STDERR,'(A,1X,A)') SOLVER_NAME,'inconsistent 1st or 2nd dimension'
    WRITE(STDERR,'(A)') SOLVER_NAME
#ifdef MPI
    CALL MPI_Abort( error )
#endif
    STOP
  ENDIF

  dh   = 1._RFREAL/REAL( nj-1 )
  rdh2 = 1._RFREAL/(dh*dh)

  DO j=2,nj-1
    DO i=1,ni
      dvar(idb:ide,i,j) = (var(idb:ide,i,j+1) - 2._RFREAL*var(idb:ide,i,j) + &
                           var(idb:ide,i,j-1))*rdh2
    ENDDO
  ENDDO
  DO j=1,1
    DO i=1,ni
      dvar(idb:ide,i,j) = 2*dvar(idb:ide,i,j+1) - dvar(idb:ide,i,j+2)
    ENDDO
  ENDDO
  DO j=nj,nj
    DO i=1,ni
      dvar(idb:ide,i,j) = 2*dvar(idb:ide,i,j-1) - dvar(idb:ide,i,j-2)
    ENDDO
  ENDDO

! finalize --------------------------------------------------------------------

END SUBROUTINE RFLO_FinDiffCompJJs


! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModFiniteDifference

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModFiniteDifference.F90,v $
! Revision 1.4  2008/12/06 08:44:16  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/12/07 08:45:51  wasistho
! added stuff for surface mesh motion EPDE
!
! Revision 1.1  2005/12/03 09:39:47  wasistho
! initial import
!
!
!
! ******************************************************************************






