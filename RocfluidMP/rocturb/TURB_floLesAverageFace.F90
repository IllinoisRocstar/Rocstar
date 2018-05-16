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
! Purpose: Average face variable to another face variable whose size depends
!          on the number of homogeneous directions. The field to be filtered 
!          is contained in faceVar and on exit the averaged data is stored in 
!          avgFaceVar
!
! Description: Only the inner field (including boundaries) is averaged.
!              If the averaging direction (denoted by homDir) is the same as 
!              the working face direction (denoted by ijk) the averaging 
!              interval ranges from the first face up to the last face, 
!              otherwise from the first center to the last center, for instance 
!              if homDir=(1,0,0) and ijk=DIRI the averaging interval is 
!              [ipnbeg:ipnend], whereas if homDir=(0,1,0) and ijk=DIRI it is 
!              [jpcbeg:jpcend]. If there is no homogeneous direction, a local
!              averaging by applying a test filter procedure and the averaged
!              data is contained in a volume array. If there is only one
!              averaging (homogeneous) direction, the averaged data is
!              contained in a surface array. For two averaging directions, it 
!              is a line array, whereas for three averaging directions 
!              (isotropic flows) the averaged data has one single value. On 
!              exit the averaged data, eiter volume, surface, line array or 
!              point (scalar) is stored in 1D array avgFaceVar. 
!              For each ijk case we perform line averaging-operation in three 
!              consecutive sweeps corresponding to the three computational 
!              directions, depending on whether a particular direction 
!              homogeneous (homDir(i)=1) or not (=0). Averaging is only 
!              performed in homogeneous direction. We use trapezoidal rule if 
!              the averaging direction is normal to the faces and ensemble 
!              rule otherwise. homDir(i)=0 for all three directions implies no 
!              homogeneous direction is present and local face-face filtering 
!              is used in this case
!
! Input: region  = current region data
!        ijk     = ijk-face is being treated
!        faceVar = face variables to be averaged 
!
! Output: avgFaceVar = the resulting averaged variables.
!
! Notes: Several generic averaging routines contained at the end of this file.
!
!******************************************************************************
!
! $Id: TURB_floLesAverageFace.F90,v 1.6 2008/12/06 08:44:43 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_FloLesAverageFace( region,ijk,faceVar,avgFaceVar )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensPhysNodes,RFLO_GetDimensDummyNodes, &
                            RFLO_GetNodeOffset
  USE TURB_ModInterfaces, ONLY : TURB_FloLesUniFiltFF, TURB_FloLesGenFiltFF
  USE ModTurbulence
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region)        :: region
  INTEGER               :: ijk
  REAL(RFREAL), POINTER :: faceVar(:),avgFaceVar(:)

! ... loop variables
  INTEGER :: i, j, k, ijkN, ijN, ikN, jkN

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER           :: ibeg,iend,jbeg,jend,kbeg,kend
  INTEGER           :: ipnbeg,ipnend,jpnbeg,jpnend,kpnbeg,kpnend
  INTEGER           :: idnbeg,idnend,jdnbeg,jdnend,kdnbeg,kdnend
  INTEGER           :: iLev,iNOff,jNOff,kNOff,ijNOff,ikNOff,jkNOff
  INTEGER           :: ibn,ien,homSum,minL,maxL,actSize,refSize
  INTEGER           :: tNDel(DIRI:DIRK),homDir(DIRI:DIRK)

  REAL(RFREAL)      :: avgPoint
  REAL(RFREAL), POINTER :: avgLine(:),avgSurf(:)
  REAL(RFREAL), POINTER :: oneDVec(:,:),avgOneDVec(:,:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_floLesAverageFace.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'TURB_FloLesAverageFace',&
  'TURB_floLesAverageFace.F90' )

! get indices and pointers --------------------------------------------------

  iLev   = region%currLevel
  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibn    = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien    = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)
  jNOff  = jdnend-jdnbeg+1
  kNOff  = kdnend-kdnbeg+1
  ikNOff = iNOff*(kdnend-kdnbeg+1)
  jkNOff = jNOff*(kdnend-kdnbeg+1)
  minL   = MIN(idnbeg,jdnbeg,kdnbeg)
  maxL   = MAX(idnend,jdnend,kdnend)

  IF (ijk==DIRI) THEN
    ibeg = ipnbeg
    iend = ipnend
    jbeg = jpnbeg
    jend = jpnend-1
    kbeg = kpnbeg
    kend = kpnend-1
  ELSEIF (ijk==DIRJ) THEN
    ibeg = ipnbeg
    iend = ipnend-1
    jbeg = jpnbeg
    jend = jpnend
    kbeg = kpnbeg
    kend = kpnend-1
  ELSEIF (ijk==DIRK) THEN
    ibeg = ipnbeg
    iend = ipnend-1
    jbeg = jpnbeg
    jend = jpnend-1
    kbeg = kpnbeg
    kend = kpnend
  ENDIF

  homDir(:) = region%turbInput%homDir(:)
  homSum    = homDir(DIRI) + homDIR(DIRJ) + homDIR(DIRK)

! for local point-averaging we use test filter width twice grid-spacing

! - Volume filtering with possibly more terms filtered across regions
!  tNDel(DIRI)=2
!  tNDel(DIRJ)=2
!  tNDel(DIRK)=2

! - Surface filtering for least terms filtered across region boundaries
  IF (ijk==DIRI) THEN
    tNDel(DIRI)=0
    tNDel(DIRJ)=2
    tNDel(DIRK)=2
  ELSEIF (ijk==DIRJ) THEN
    tNDel(DIRI)=2
    tNDel(DIRJ)=0
    tNDel(DIRK)=2
  ELSEIF (ijk==DIRK) THEN
    tNDel(DIRI)=2
    tNDel(DIRJ)=2
    tNDel(DIRK)=0
  ENDIF

  IF (region%turbInput%filterType == FILTYPE_NONUNIF) THEN
    IF ((tNDel(DIRI)/=2 .AND. tNDel(DIRI)/=0) .OR. &
        (tNDel(DIRJ)/=2 .AND. tNDel(DIRJ)/=0) .OR. &
        (tNDel(DIRK)/=2 .AND. tNDel(DIRK)/=0)) THEN 
      CALL ErrorStop( global,ERR_TURB_FIXPARAM,__LINE__, &
      'filter coeff. possibly not allocated if tNDel =/ 2 for nonunif. filter.')
    ENDIF
  ENDIF

! allocate temporary arrays
    
  ALLOCATE( avgLine(minL:maxL),avgSurf(ibn:ien) )
  ALLOCATE( oneDVec(1,ibn:ien ),avgOneDVec(1,ibn:ien) )

! perform averaging process depending on number of homogeneous directions

  IF (homSum == 3) THEN

! - all three directions are homogeneous

    actSize=1

    IF (homDir(DIRI)==ACTIVE .AND. homDir(DIRJ)==ACTIVE) THEN

      IF (ijk==DIRI) THEN
! ----- averaging in all direction; trapezoidal in i-direction and ensemble 
!       in other directions
        CALL Trapezoid3Lijk( ibeg,iend,jbeg,jend,kbeg,kend,faceVar,avgSurf )
        CALL Ensemble2Loops( jbeg,jend,kbeg,kend,jNOff,avgSurf,avgLine )
        CALL Ensemble1Loop( kbeg,kend,avgLine,avgPoint )

      ELSEIF (ijk==DIRJ) THEN
! ----- averaging in all direction; trapezoidal in j-direction and ensemble 
!       in other directions
        CALL Ensemble3Lijk( ibeg,iend,jbeg,jend,kbeg,kend,faceVar,avgSurf )
        CALL Trapezoid2Loops( jbeg,jend,kbeg,kend,jNOff,avgSurf,avgLine )
        CALL Ensemble1Loop( kbeg,kend,avgLine,avgPoint )

      ELSEIF (ijk==DIRK) THEN
! ----- averaging in all direction; trapezoidal in k-direction and ensemble 
!       in other directions
        CALL Ensemble3Lijk( ibeg,iend,jbeg,jend,kbeg,kend,faceVar,avgSurf )
        CALL Ensemble2Loops( jbeg,jend,kbeg,kend,jNOff,avgSurf,avgLine )
        CALL Trapezoid1Loop( kbeg,kend,avgLine,avgPoint )

      ENDIF 
    ELSE
      CALL ErrorStop( global,ERR_TURB_LESINPUT,__LINE__,'homDir is not consistent.')
    ENDIF       

    avgFaceVar = avgPoint

  ELSEIF (homSum == 2) THEN

    IF (homDir(DIRJ)==ACTIVE .AND. homDir(DIRK)==ACTIVE) THEN
      actSize = iNOff

      IF (ijk==DIRI) THEN
! ----- averaging in j and k direction; ensemble in both direcions
        CALL Ensemble3Ljki( jbeg,jend,kbeg,kend,ibeg,iend,faceVar,avgSurf )
        CALL Ensemble2Loops( kbeg,kend,ibeg,iend,kNOff,avgSurf,avgLine )

      ELSEIF (ijk==DIRJ) THEN
! ----- averaging in j and k direction; trapezoidal in j direction and
!       ensemble in k direcions
        CALL Trapezoid3Ljki( jbeg,jend,kbeg,kend,ibeg,iend,faceVar,avgSurf )
        CALL Ensemble2Loops( kbeg,kend,ibeg,iend,kNOff,avgSurf,avgLine )

      ELSEIF (ijk==DIRK) THEN
! ----- averaging in j and k direction; trapezoidal in k direction and
!       ensemble in j direcions
        CALL Ensemble3Ljki( jbeg,jend,kbeg,kend,ibeg,iend,faceVar,avgSurf )
        CALL Trapezoid2Loops( kbeg,kend,ibeg,iend,kNOff,avgSurf,avgLine )

      ENDIF

! --- copy to avgFaceVar
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff)
            avgFaceVar(ijkN) = avgLine(i)
          ENDDO
        ENDDO
      ENDDO

    ELSEIF (homDir(DIRI)==ACTIVE .AND. homDir(DIRK)==ACTIVE) THEN
      actSize = jNOff

      IF (ijk==DIRI) THEN
! ----- averaging in i and k direction; trapezoidal in i direction and
!       ensemble in k direcions:
        CALL Trapezoid3Likj( ibeg,iend,kbeg,kend,jbeg,jend,faceVar,avgSurf )
        CALL Ensemble2Loops( kbeg,kend,jbeg,jend,kNOff,avgSurf,avgLine )

      ELSEIF (ijk==DIRJ) THEN
! ----- averaging in i and k direction; ensemble in both direcions
        CALL Ensemble3Likj( ibeg,iend,kbeg,kend,jbeg,jend,faceVar,avgSurf )
        CALL Ensemble2Loops( kbeg,kend,jbeg,jend,kNOff,avgSurf,avgLine )

      ELSEIF (ijk==DIRK) THEN
! ----- averaging in i and k direction; trapezoidal in k direction and
!       ensemble in i direcions
        CALL Ensemble3Likj( ibeg,iend,kbeg,kend,jbeg,jend,faceVar,avgSurf )
        CALL Trapezoid2Loops( kbeg,kend,jbeg,jend,kNOff,avgSurf,avgLine )

      ENDIF

! --- copy to avgFaceVar
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff)
            avgFaceVar(ijkN) = avgLine(j)
          ENDDO
        ENDDO
      ENDDO

    ELSEIF (homDir(DIRI)==ACTIVE .AND. homDir(DIRJ)==ACTIVE) THEN
      actSize = kNOff

      IF (ijk==DIRI) THEN
! ----- averaging in i and j direction; trapezoidal in i direction and
!       ensemble in j direcions
        CALL Trapezoid3Lijk( ibeg,iend,jbeg,jend,kbeg,kend,faceVar,avgSurf )
        CALL Ensemble2Loops( jbeg,jend,kbeg,kend,jNOff,avgSurf,avgLine )

      ELSEIF (ijk==DIRJ) THEN
! ----- averaging in i and j direction; trapezoidal in j direction and
!       ensemble in i direcions
        CALL Ensemble3Lijk( ibeg,iend,jbeg,jend,kbeg,kend,faceVar,avgSurf )
        CALL Trapezoid2Loops( jbeg,jend,kbeg,kend,jNOff,avgSurf,avgLine )

      ELSEIF (ijk==DIRK) THEN
! ----- averaging in i and j direction; ensemble in both directions
        CALL Ensemble3Lijk( ibeg,iend,jbeg,jend,kbeg,kend,faceVar,avgSurf )
        CALL Ensemble2Loops( jbeg,jend,kbeg,kend,jNOff,avgSurf,avgLine )

      ENDIF 

! --- copy to avgFaceVar
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff)
            avgFaceVar(ijkN) = avgLine(k)
          ENDDO
        ENDDO
      ENDDO

    ELSE
      CALL ErrorStop( global,ERR_TURB_LESINPUT,__LINE__,'homDir is not consistent.')
    ENDIF       

  ELSEIF (homSum == 1) THEN

    IF (homDir(DIRJ)/=ACTIVE .AND. homDir(DIRK)/=ACTIVE) THEN
      actSize = jkNOff

      IF (ijk==DIRI) THEN
! ----- trapezoidal averaging in i direction
        CALL Trapezoid3Lijk( ibeg,iend,jbeg,jend,kbeg,kend,faceVar,avgSurf )

      ELSEIF (ijk==DIRJ) THEN
! ----- ensemble averaging in i direction
        CALL Ensemble3Lijk( ibeg,iend,jbeg,jend,kbeg,kend,faceVar,avgSurf )

      ELSEIF (ijk==DIRK) THEN
! ----- ensemble averaging in i direction.
        CALL Ensemble3Lijk( ibeg,iend,jbeg,jend,kbeg,kend,faceVar,avgSurf )

      ENDIF 

! --- copy to avgFaceVar
      DO k=kbeg,kend
        DO j=jbeg,jend
          jkN = IndIJ(j ,k ,jNOff)
          DO i=ibeg,iend
            ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff)
            avgFaceVar(ijkN) = avgSurf(jkN)
          ENDDO
        ENDDO
      ENDDO

    ELSEIF (homDir(DIRI)/=ACTIVE .AND. homDir(DIRK)/=ACTIVE) THEN
      actSize = ikNOff

      IF (ijk==DIRI) THEN
! ----- ensemble averaging in j direction.
        CALL Ensemble3Ljik( jbeg,jend,ibeg,iend,kbeg,kend,faceVar,avgSurf )

      ELSEIF (ijk==DIRJ) THEN
! ----- trapezoidal averaging in j direction
        CALL Trapezoid3Ljik( jbeg,jend,ibeg,iend,kbeg,kend,faceVar,avgSurf )

      ELSEIF (ijk==DIRK) THEN
! ----- ensemble averaging in j direction
        CALL Ensemble3Ljik( jbeg,jend,ibeg,iend,kbeg,kend,faceVar,avgSurf )

      ENDIF 

! --- copy to avgFaceVar
      DO k=kbeg,kend
        DO i=ibeg,iend
          ikN = IndIJ(i ,k ,iNOff)
          DO j=jbeg,jend
            ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff)
            avgFaceVar(ijkN) = avgSurf(ikN)
          ENDDO
        ENDDO
      ENDDO

    ELSEIF (homDir(DIRI)/=ACTIVE .AND. homDir(DIRJ)/=ACTIVE) THEN
      actSize = ijNOff

      IF (ijk==DIRI) THEN
! ----- ensemble averaging in k direction
        CALL Ensemble3Lkij( kbeg,kend,ibeg,iend,jbeg,jend,faceVar,avgSurf )

      ELSEIF (ijk==DIRJ) THEN
! ----- ensemble averaging in k direction.
        CALL Ensemble3Lkij( kbeg,kend,ibeg,iend,jbeg,jend,faceVar,avgSurf )

      ELSEIF (ijk==DIRK) THEN
! ----- trapezoidal averaging in k direction.
        CALL Trapezoid3Lkij( kbeg,kend,ibeg,iend,jbeg,jend,faceVar,avgSurf )

      ENDIF 

! --- copy to avgFaceVar
      DO j=jbeg,jend
        DO i=ibeg,iend
          ijN = IndIJ(i ,j ,iNOff)
          DO k=kbeg,kend
            ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff)
            avgFaceVar(ijkN) = avgSurf(ijN)
          ENDDO
        ENDDO
      ENDDO

    ELSE
      CALL ErrorStop( global,ERR_TURB_LESINPUT,__LINE__,'homDir is not consistent.')
    ENDIF       

  ELSEIF (homSum == 0) THEN

    IF (homDir(DIRI)/=ACTIVE .AND. homDir(DIRJ)/=ACTIVE) THEN

      actSize = ijNOff*(kdnend-kdnbeg+1)

      DO i = ibn,ien
        oneDVec(1,i) = faceVar(i) 
      ENDDO
      IF (region%turbInput%filterType == FILTYPE_UNIFORM) THEN
        CALL TURB_FloLesUniFiltFF( region,ijk,tNDel,1,1,oneDVec,avgOneDVec )
      ELSE
        CALL TURB_FloLesGenFiltFF( region,ijk,tNDel,1,1,oneDVec,avgOneDVec )
      ENDIF
      DO i = ibn,ien
        avgFaceVar(i) = avgOneDVec(1,i)
      ENDDO
    ELSE
      CALL ErrorStop( global,ERR_TURB_LESINPUT,__LINE__,'homDir is not consistent.')
    ENDIF       
  ELSE
    CALL ErrorStop( global,ERR_TURB_LESINPUT,__LINE__,'homDir is not consistent.')
  ENDIF

! check size of averaged variable array

!  IF (refSize /= actSize) THEN
!    CALL ErrorStop( global,ERR_TURB_LESINPUT,__LINE__,'coefDim is not consistent.')
!  ENDIF 

! deallocate temporary arrays

  DEALLOCATE( avgLine,avgSurf,oneDVec,avgOneDVec )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

! ==============================================================================
! summation subroutines
! ==============================================================================

CONTAINS

  SUBROUTINE Trapezoid1Loop( ib,ie,lineVar,pointVar )

! ... parameters
    INTEGER               :: ib,ie
    REAL(RFREAL), POINTER :: lineVar(:)
    REAL(RFREAL)          :: pointVar  

! ... local variables
    INTEGER               :: i
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    pointVar = 0._RFREAL

    sum = 0._RFREAL
    DO i=ib+1,ie-1
      sum = sum + lineVar(i)
    ENDDO
    sum = sum + 0.5_RFREAL*(lineVar(ib)+lineVar(ie))
    pointVar = sum/REAL(ie-ib)

  END SUBROUTINE Trapezoid1Loop

!#####################################################################

  SUBROUTINE Ensemble1Loop( ib,ie,lineVar,pointVar )

! ... parameters
    INTEGER               :: ib,ie
    REAL(RFREAL), POINTER :: lineVar(:)
    REAL(RFREAL)          :: pointVar  

! ... local variables
    INTEGER               :: i
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    pointVar = 0._RFREAL

    sum = 0._RFREAL
    DO i=ib,ie
      sum = sum + lineVar(i)
    ENDDO
    pointVar = sum/REAL(ie-ib+1)

  END SUBROUTINE Ensemble1Loop

!#####################################################################

  SUBROUTINE Trapezoid2Loops( ib,ie,jb,je,nOff,surfVar,lineVar )

! ... parameters
    INTEGER               :: ib,ie,jb,je,nOff
    REAL(RFREAL), POINTER :: surfVar(:),lineVar(:)

! ... local variables
    INTEGER               :: i,j,ijN,ijNb,ijNe
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    lineVar = 0._RFREAL

    DO j=jb,je
      sum = 0._RFREAL
      DO i=ib+1,ie-1
        ijN = IndIJ(i ,j ,nOff)
        sum = sum + surfVar(ijN)
      ENDDO
      ijNb = IndIJ(ib ,j ,nOff) 
      ijNe = IndIJ(ie ,j ,nOff)
      sum  = sum + 0.5_RFREAL*(surfVar(ijNb)+surfVar(ijNe))
      lineVar(j) = sum/REAL(ie-ib)
    ENDDO

  END SUBROUTINE Trapezoid2Loops

!#####################################################################

  SUBROUTINE Ensemble2Loops( ib,ie,jb,je,nOff,surfVar,lineVar )

! ... parameters
    INTEGER               :: ib,ie,jb,je,nOff
    REAL(RFREAL), POINTER :: surfVar(:),lineVar(:)

! ... local variables
    INTEGER               :: i,j,ijN
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    lineVar = 0._RFREAL

    DO j=jb,je
      sum = 0._RFREAL
      DO i=ib,ie
        ijN = IndIJ(i ,j ,nOff) 
        sum = sum+surfVar(ijN)
      ENDDO
      lineVar(j) = sum/REAL(ie-ib+1)
    ENDDO

  END SUBROUTINE Ensemble2Loops

!#####################################################################

  SUBROUTINE Trapezoid3Lijk( ib,ie,jb,je,kb,ke,volVar,surfVar )

! ... parameters
    INTEGER               :: ib,ie,jb,je,kb,ke
    REAL(RFREAL), POINTER :: volVar(:),surfVar(:)

! ... local variables
    INTEGER               :: i,j,k,jkN,ijkN,ijkNb,ijkNe
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    surfVar = 0._RFREAL

    DO k=kb,ke
      DO j=jb,je
        jkN = IndIJ(j ,k ,jNOff) 
        sum = 0._RFREAL
        DO i=ib+1,ie-1
          ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff)
          sum  = sum + volVar(ijkN)
        ENDDO
        ijkNb = IndIJK(ib ,j ,k ,iNOff,ijNOff) 
        ijkNe = ijkNb + ie-ib
        sum   = sum + 0.5_RFREAL*(volVar(ijkNb)+volVar(ijkNe))
        surfVar(jkN) = sum/REAL(ie-ib)
      ENDDO
    ENDDO          

  END SUBROUTINE Trapezoid3Lijk
!===============================================================
  SUBROUTINE Trapezoid3Likj( ib,ie,kb,ke,jb,je,volVar,surfVar )

! ... parameters
    INTEGER               :: ib,ie,jb,je,kb,ke
    REAL(RFREAL), POINTER :: volVar(:),surfVar(:)

! ... local variables
    INTEGER               :: i,j,k,kjN,ijkN,ijkNb,ijkNe
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    surfVar = 0._RFREAL

    DO j=jb,je
      DO k=kb,ke
        kjN = IndIJ(k ,j ,kNOff) 
        sum = 0._RFREAL
        DO i=ib+1,ie-1
          ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff) 
          sum  = sum + volVar(ijkN)
        ENDDO
        ijkNb = IndIJK(ib ,j ,k ,iNOff,ijNOff) 
        ijkNe = ijkNb + ie-ib
        sum   = sum + 0.5_RFREAL*(volVar(ijkNb)+volVar(ijkNe))
        surfVar(kjN) = sum/REAL(ie-ib)
      ENDDO
    ENDDO          

  END SUBROUTINE Trapezoid3Likj
!===============================================================
  SUBROUTINE Trapezoid3Ljik( jb,je,ib,ie,kb,ke,volVar,surfVar )

! ... parameters
    INTEGER               :: ib,ie,jb,je,kb,ke
    REAL(RFREAL), POINTER :: volVar(:),surfVar(:)

! ... local variables
    INTEGER               :: i,j,k,ikN,ijkN,ijkNb,ijkNe
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    surfVar = 0._RFREAL

    DO k=kb,ke
      DO i=ib,ie
        ikN = IndIJ(i ,k ,iNOff) 
        sum = 0._RFREAL
        DO j=jb+1,je-1
          ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff) 
          sum  = sum + volVar(ijkN)
        ENDDO
        ijkNb = IndIJK(i ,jb ,k ,iNOff,ijNOff) 
        ijkNe = ijkNb + (je-jb)*iNOff
        sum   = sum + 0.5_RFREAL*(volVar(ijkNb)+volVar(ijkNe))
        surfVar(ikN) = sum/REAL(je-jb)
      ENDDO
    ENDDO          

  END SUBROUTINE Trapezoid3Ljik
!===============================================================
  SUBROUTINE Trapezoid3Ljki( jb,je,kb,ke,ib,ie,volVar,surfVar )

! ... parameters
    INTEGER               :: ib,ie,jb,je,kb,ke
    REAL(RFREAL), POINTER :: volVar(:),surfVar(:)

! ... local variables
    INTEGER               :: i,j,k,kiN,ijkN,ijkNb,ijkNe
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    surfVar = 0._RFREAL

    DO i=ib,ie
      DO k=kb,ke
        kiN = IndIJ(k ,i ,kNOff) 
        sum = 0._RFREAL
        DO j=jb+1,je-1
          ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff) 
          sum  = sum + volVar(ijkN)
        ENDDO
        ijkNb = IndIJK(i ,jb ,k ,iNOff,ijNOff) 
        ijkNe = ijkNb + (je-jb)*iNOff
        sum   = sum + 0.5_RFREAL*(volVar(ijkNb)+volVar(ijkNe))
        surfVar(kiN) = sum/REAL(je-jb)
      ENDDO
    ENDDO          

  END SUBROUTINE Trapezoid3Ljki
!===============================================================
  SUBROUTINE Trapezoid3Lkij( kb,ke,ib,ie,jb,je,volVar,surfVar )

! ... parameters
    INTEGER               :: ib,ie,jb,je,kb,ke
    REAL(RFREAL), POINTER :: volVar(:),surfVar(:)

! ... local variables
    INTEGER               :: i,j,k,ijN,ijkN,ijkNb,ijkNe
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    surfVar = 0._RFREAL

    DO j=jb,je
      DO i=ib,ie
        ijN = IndIJ(i ,j ,iNOff) 
        sum = 0._RFREAL
        DO k=kb+1,ke-1
          ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff) 
          sum  = sum + volVar(ijkN)
        ENDDO
        ijkNb = IndIJK(i ,j ,kb ,iNOff,ijNOff) 
        ijkNe = ijkNb + (ke-kb)*ijNOff
        sum   = sum + 0.5_RFREAL*(volVar(ijkNb)+volVar(ijkNe))
        surfVar(ijN) = sum/REAL(ke-kb)
      ENDDO
    ENDDO          

  END SUBROUTINE Trapezoid3Lkij
!===============================================================
  SUBROUTINE Trapezoid3Lkji( kb,ke,jb,je,ib,ie,volVar,surfVar )

! ... parameters
    INTEGER               :: ib,ie,jb,je,kb,ke
    REAL(RFREAL), POINTER :: volVar(:),surfVar(:)

! ... local variables
    INTEGER               :: i,j,k,jiN,ijkN,ijkNb,ijkNe
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    surfVar = 0._RFREAL

    DO i=ib,ie
      DO j=jb,je
        jiN = IndIJ(j ,i ,jNOff) 
        sum = 0._RFREAL
        DO k=kb+1,ke-1
          ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff) 
          sum  = sum + volVar(ijkN)
        ENDDO
        ijkNb = IndIJK(i ,j ,kb ,iNOff,ijNOff) 
        ijkNe = ijkNb + (ke-kb)*ijNOff
        sum   = sum + 0.5_RFREAL*(volVar(ijkNb)+volVar(ijkNe))
        surfVar(jiN) = sum/REAL(ke-kb)
      ENDDO
    ENDDO          

  END SUBROUTINE Trapezoid3Lkji

!#####################################################################

  SUBROUTINE Ensemble3Lijk( ib,ie,jb,je,kb,ke,volVar,surfVar )

! ... parameters
    INTEGER               :: ib,ie,jb,je,kb,ke
    REAL(RFREAL), POINTER :: volVar(:),surfVar(:)

! ... local variables
    INTEGER               :: i,j,k,jkN,ijkN
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    surfVar = 0._RFREAL

    DO k=kb,ke
      DO j=jb,je
        jkN = IndIJ(j ,k ,jNOff) 
        sum = 0._RFREAL
        DO i=ib,ie
          ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff) 
          sum  = sum + volVar(ijkN)
        ENDDO
        surfVar(jkN) = sum/REAL(ie-ib+1)
      ENDDO
    ENDDO          

  END SUBROUTINE Ensemble3Lijk
!===============================================================
  SUBROUTINE Ensemble3Likj( ib,ie,kb,ke,jb,je,volVar,surfVar )

! ... parameters
    INTEGER               :: ib,ie,jb,je,kb,ke
    REAL(RFREAL), POINTER :: volVar(:),surfVar(:)

! ... local variables
    INTEGER               :: i,j,k,kjN,ijkN
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    surfVar = 0._RFREAL

    DO j=jb,je
      DO k=kb,ke
        kjN = IndIJ(k ,j ,kNOff) 
        sum = 0._RFREAL
        DO i=ib,ie
          ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff) 
          sum  = sum + volVar(ijkN)
        ENDDO
        surfVar(kjN) = sum/REAL(ie-ib+1)
      ENDDO
    ENDDO          

  END SUBROUTINE Ensemble3Likj
!===============================================================
  SUBROUTINE Ensemble3Ljik( jb,je,ib,ie,kb,ke,volVar,surfVar )

! ... parameters
    INTEGER               :: ib,ie,jb,je,kb,ke
    REAL(RFREAL), POINTER :: volVar(:),surfVar(:)

! ... local variables
    INTEGER               :: i,j,k,ikN,ijkN
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    surfVar = 0._RFREAL

    DO k=kb,ke
      DO i=ib,ie
        ikN = IndIJ(i ,k ,iNOff) 
        sum = 0._RFREAL
        DO j=jb,je
          ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff) 
          sum  = sum + volVar(ijkN)
        ENDDO
        surfVar(ikN) = sum/REAL(je-jb+1)
      ENDDO
    ENDDO          

  END SUBROUTINE Ensemble3Ljik
!===============================================================
  SUBROUTINE Ensemble3Ljki( jb,je,kb,ke,ib,ie,volVar,surfVar )

! ... parameters
    INTEGER               :: ib,ie,jb,je,kb,ke
    REAL(RFREAL), POINTER :: volVar(:),surfVar(:)

! ... local variables
    INTEGER               :: i,j,k,kiN,ijkN
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    surfVar = 0._RFREAL

    DO i=ib,ie
      DO k=kb,ke
        kiN = IndIJ(k ,i ,kNOff) 
        sum = 0._RFREAL
        DO j=jb,je
          ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff) 
          sum  = sum + volVar(ijkN)
        ENDDO
        surfVar(kiN) = sum/REAL(je-jb+1)
      ENDDO
    ENDDO          

  END SUBROUTINE Ensemble3Ljki
!===============================================================
  SUBROUTINE Ensemble3Lkij( kb,ke,ib,ie,jb,je,volVar,surfVar )

! ... parameters
    INTEGER               :: ib,ie,jb,je,kb,ke
    REAL(RFREAL), POINTER :: volVar(:),surfVar(:)

! ... local variables
    INTEGER               :: i,j,k,ijN,ijkN
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    surfVar = 0._RFREAL

    DO j=jb,je
      DO i=ib,ie
        ijN = IndIJ(i ,j ,iNOff) 
        sum = 0._RFREAL
        DO k=kb,ke
          ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff) 
          sum  = sum + volVar(ijkN)
        ENDDO
        surfVar(ijN) = sum/REAL(ke-kb+1)
      ENDDO
    ENDDO          

  END SUBROUTINE Ensemble3Lkij
!===============================================================
  SUBROUTINE Ensemble3Lkji( kb,ke,jb,je,ib,ie,volVar,surfVar )

! ... parameters
    INTEGER               :: ib,ie,jb,je,kb,ke
    REAL(RFREAL), POINTER :: volVar(:),surfVar(:)

! ... local variables
    INTEGER               :: i,j,k,jiN,ijkN
    REAL(RFREAL)          :: sum

! - start procedure ---------------------------------------

    surfVar = 0._RFREAL

    DO i=ib,ie
      DO j=jb,je
        jiN = IndIJ(j ,i ,jNOff) 
        sum = 0._RFREAL
        DO k=kb,ke
          ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff) 
          sum  = sum + volVar(ijkN)
        ENDDO
        surfVar(jiN) = sum/REAL(ke-kb+1)
      ENDDO
    ENDDO          

  END SUBROUTINE Ensemble3Lkji

END SUBROUTINE TURB_FloLesAverageFace

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_floLesAverageFace.F90,v $
! Revision 1.6  2008/12/06 08:44:43  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:55  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2004/06/03 02:14:54  wasistho
! provided surface and volume (commented) local filtering
!
! Revision 1.3  2004/05/29 03:39:55  wasistho
! used local surface filtering instead of volume filtering in case of no hom. direction
!
! Revision 1.2  2004/03/12 02:55:36  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/08 23:35:45  wasistho
! changed turb nomenclature
!
! Revision 1.4  2004/02/28 01:19:10  wasistho
! reset filterwidth ratio for local averaging to 2
!
! Revision 1.3  2004/02/18 21:27:51  wasistho
! used smaller stencil for local averaging by filtering
!
! Revision 1.2  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!******************************************************************************







