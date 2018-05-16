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
! Purpose: Collection of routines to compute force and moment coefficients.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLO_ModForcesMoments.F90,v 1.6 2008/12/06 08:44:16 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLO_ModForcesMoments

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal,     ONLY: t_global
  USE ModGrid,       ONLY: t_grid
  USE ModBndPatch,   ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLO_ComputeIntegralForceMomCo, &  
            RFLO_ComputePatchForceMomCo, &
            RFLO_FindPatchCoeffsGlo, &
            RFLO_OpenForceMomCoFile, &
            RFLO_WriteIntegralForceMomCo, &
            RFLO_WritePatchCoeffsInfo

! private :

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLO_ModForcesMoments.F90,v $ $Revision: 1.6 $'

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  



! *******************************************************************************
!
! Purpose: Compute patch force and moment coefficients.
!
! Description: none.
!
! Input: region = current region data
!
! Output: patch%forceCoeffs, patch%momentCoeffs = 
!         patch force and moment coefficients computed
!
! Notes: none.
!
!*******************************************************************************

SUBROUTINE RFLO_ComputePatchForceMomCo( region )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, &
                            RFLO_GetPatchIndicesNodes, &
                            RFLO_GetPatchDirection, RFLO_GetNodeOffset
  IMPLICIT NONE

#include "Indexing.h"
        
! ... parameters
  TYPE (t_region) :: region

! ... loop variables
  INTEGER :: iPatch, i, j, k

! ... local variables
  INTEGER :: iLev, lbound, bcType, ibeg, iend, jbeg, jend, kbeg, kend
  INTEGER :: inbeg, inend, jnbeg, jnend, knbeg, knend, nOff, iNOff, ijNOff
  INTEGER :: idir, jdir, kdir, inode, jnode, knode, n1, n2, i2d, ijkN
  REAL(RFREAL) :: boxXmin, boxXmax, boxYmin, boxYmax, boxZmin, boxZmax
  REAL(RFREAL) :: xRef, yRef, zRef, xMin, xMax, yMin, yMax, zMin, zMax
  REAL(RFREAL) :: fpx, fpy, fpz, fvx, fvy, fvz, mpx, mpy, mpz, mvx, mvy, mvz
  REAL(RFREAL) :: sfx, sfy, sfz, sfm, xc, yc, zc, cp, cfx, cfy, cfz, ch
  REAL(RFREAL), POINTER :: xyz(:,:), sFace(:,:), cFace(:,:)
  TYPE (t_global), POINTER :: global
  TYPE (t_patch),  POINTER :: patch

! ******************************************************************************
    
  global => region%global
  CALL RegisterFunction(global,'RFLO_ComputePatchForceMomCo',&
       'RFLO_ModForcesMoments.F90')

! get dimensions, parameters and pointers --------------------------------------

  iLev = region%currLevel
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  xRef = global%forceRefXCoord
  yRef = global%forceRefYCoord
  zRef = global%forceRefZCoord

  boxXmin = global%acBndBoxXmin
  boxXmax = global%acBndBoxXmax
  boxYmin = global%acBndBoxYmin
  boxYmax = global%acBndBoxYmax
  boxZmin = global%acBndBoxZmin
  boxZmax = global%acBndBoxZmax

  xyz => region%levels(iLev)%grid%xyz

! initialize global values

  global%forceCoeffs  = 0._RFREAL
  global%momentCoeffs = 0._RFREAL
              
! loop over patches ------------------------------------------------------------
        
  DO iPatch = 1,region%nPatches
    patch  => region%levels(iLev)%patches(iPatch)
    lbound =  patch%lbound
    bcType =  patch%bcType
    nOff   =  ABS(patch%l1end-patch%l1beg) + 1

! - initialized local variables
            
    fpx = 0.0_RFREAL
    fpy = 0.0_RFREAL
    fpz = 0.0_RFREAL
                        
    fvx = 0.0_RFREAL
    fvy = 0.0_RFREAL
    fvz = 0.0_RFREAL                  
      
    mpx = 0.0_RFREAL
    mpy = 0.0_RFREAL
    mpz = 0.0_RFREAL
                     
    mvx = 0.0_RFREAL
    mvy = 0.0_RFREAL
    mvz = 0.0_RFREAL      

! - get patch dimensions and appropriate face vector

    CALL RFLO_GetPatchIndices( region,patch,iLev,ibeg,iend,jbeg,jend,kbeg,kend )
    CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                    inbeg,inend,jnbeg,jnend,knbeg,knend )
    CALL RFLO_GetPatchDirection( patch,idir,jdir,kdir )

    inode = 0
    jnode = 0
    knode = 0
    IF (lbound==2 .OR. lbound==4 .OR. lbound==6) THEN
      inode = -idir
      jnode = -jdir
      knode = -kdir
    ENDIF

    IF (lbound==1 .OR. lbound==2) THEN
      sFace => region%levels(iLev)%grid%si
      cFace => region%levels(iLev)%grid%cfcI
    ELSE IF (lbound==3 .OR. lbound==4) THEN
      sFace => region%levels(iLev)%grid%sj
      cFace => region%levels(iLev)%grid%cfcJ
    ELSE IF (lbound==5 .OR. lbound==6) THEN
      sFace => region%levels(iLev)%grid%sk
      cFace => region%levels(iLev)%grid%cfcK
    ENDIF

! - search for min/max xyz

    xMin =  HUGE( 1._RFREAL )
    xMax = -HUGE( 1._RFREAL )
    yMin =  HUGE( 1._RFREAL )
    yMax = -HUGE( 1._RFREAL )
    zMin =  HUGE( 1._RFREAL )
    zMax = -HUGE( 1._RFREAL )

    DO k=knbeg,knend
      DO j=jnbeg,jnend
        DO i=inbeg,inend
          ijkN = IndIJK(i,j,k,iNOff,ijNOff)
          xMin = MIN( xMin,xyz(XCOORD,ijkN) )
          xMax = MAX( xMax,xyz(XCOORD,ijkN) )
          yMin = MIN( yMin,xyz(YCOORD,ijkN) )
          yMax = MAX( yMax,xyz(YCOORD,ijkN) )
          zMin = MIN( zMin,xyz(ZCOORD,ijkN) )
          zMax = MAX( zMax,xyz(ZCOORD,ijkN) )
        ENDDO
      ENDDO
    ENDDO

! - loop over faces ----------------------------------------------------------

    IF (((bcType>=BC_SLIPWALL   .AND. bcType<=BC_SLIPWALL+BC_RANGE)   .OR. &
         (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) .OR. &
         (bcType>=BC_INJECTION  .AND. bcType<=BC_INJECTION+BC_RANGE)) .AND. &
        ((xMin > boxXmin .AND. xMax < boxXmax) .AND. &
         (yMin > boxYmin .AND. yMax < boxYmax) .AND. & 
         (zMin > boxZmin .AND. zMax < boxZmax))) THEN

      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend

            ijkN = IndIJK(i+inode,j+jnode,k+knode,iNOff,ijNOff)
      
            sfx = sFace(XCOORD,ijkN)
            sfy = sFace(YCOORD,ijkN)
            sfz = sFace(ZCOORD,ijkN)                
            sfm = sFace(XYZMAG,ijkN)
      
            xc  = cFace(XCOORD,ijkN)
            yc  = cFace(YCOORD,ijkN)
            zc  = cFace(ZCOORD,ijkN)                
      
! --------- get coefficients

            IF      (lbound==1 .OR. lbound==2) THEN
              n1 = j - jbeg
              n2 = k - kbeg
            ELSE IF (lbound==3 .OR. lbound==4) THEN
              n1 = k - kbeg
              n2 = i - ibeg
            ELSE IF (lbound==5 .OR. lbound==6) THEN
              n1 = i - ibeg
              n2 = j - jbeg
            ENDIF
            i2d = IndIJ(n1,n2,nOff)
      
            cp  = patch%cp(i2d)
            cfx = patch%cf(XCOORD,i2d)
            cfy = patch%cf(YCOORD,i2d)
            cfz = patch%cf(ZCOORD,i2d)
            ch  = patch%ch(i2d)                

! --------- compute contributions to force and moment coefficients
      
            fpx = fpx + cp*sfx*sfm
            fpy = fpy + cp*sfy*sfm
            fpz = fpz + cp*sfz*sfm

            fvx = fvx + cfx*sfm
            fvy = fvy + cfy*sfm
            fvz = fvz + cfz*sfm
                        
            mpx = mpx - cp*(sfy*(zc - zRef) + sfz*(yc - yRef))*sfm
            mpy = mpy + cp*(sfx*(zc - zRef) - sfz*(xc - xRef))*sfm
            mpz = mpz + cp*(sfy*(xc - xRef) - sfx*(yc - yRef))*sfm

            mvx = mvx - (cfy*(zc - zRef) + cfz*(yc - yRef))*sfm
            mvy = mvy + (cfx*(zc - zRef) - cfz*(xc - xRef))*sfm
            mvz = mvz + (cfy*(xc - xRef) - cfx*(yc - yRef))*sfm
          ENDDO ! i
        ENDDO   ! j
      ENDDO     ! k

! --- normalize and store coefficients ----------------------------------------
                    
      patch%forceCoeffs(XCOORD,FORCES_PRESS)  = fpx
      patch%forceCoeffs(YCOORD,FORCES_PRESS)  = fpy
      patch%forceCoeffs(ZCOORD,FORCES_PRESS)  = fpz
       
      patch%forceCoeffs(XCOORD,FORCES_VISC)   = fvx
      patch%forceCoeffs(YCOORD,FORCES_VISC)   = fvy
      patch%forceCoeffs(ZCOORD,FORCES_VISC)   = fvz   

      patch%momentCoeffs(XCOORD,FORCES_PRESS) = mpx
      patch%momentCoeffs(YCOORD,FORCES_PRESS) = mpy
      patch%momentCoeffs(ZCOORD,FORCES_PRESS) = mpz
       
      patch%momentCoeffs(XCOORD,FORCES_VISC)  = mvx
      patch%momentCoeffs(YCOORD,FORCES_VISC)  = mvy
      patch%momentCoeffs(ZCOORD,FORCES_VISC)  = mvz

      global%forceCoeffs(XCOORD,1)  = global%forceCoeffs(XCOORD,1)  + fpx
      global%forceCoeffs(YCOORD,1)  = global%forceCoeffs(YCOORD,1)  + fpy
      global%forceCoeffs(ZCOORD,1)  = global%forceCoeffs(ZCOORD,1)  + fpz
      global%forceCoeffs(XCOORD,2)  = global%forceCoeffs(XCOORD,2)  + fvx
      global%forceCoeffs(YCOORD,2)  = global%forceCoeffs(YCOORD,2)  + fvy
      global%forceCoeffs(ZCOORD,2)  = global%forceCoeffs(ZCOORD,2)  + fvz
      global%momentCoeffs(XCOORD,1) = global%momentCoeffs(XCOORD,1) + mpx
      global%momentCoeffs(YCOORD,1) = global%momentCoeffs(YCOORD,1) + mpy
      global%momentCoeffs(ZCOORD,1) = global%momentCoeffs(ZCOORD,1) + mpz
      global%momentCoeffs(XCOORD,2) = global%momentCoeffs(XCOORD,2) + mvx
      global%momentCoeffs(YCOORD,2) = global%momentCoeffs(YCOORD,2) + mvy
      global%momentCoeffs(ZCOORD,2) = global%momentCoeffs(ZCOORD,2) + mvz

    ENDIF ! bcType and boxBnd
  ENDDO   ! iPatch

! finalize ---------------------------------------------------------------------

  CALL DeregisterFunction(global)

END SUBROUTINE RFLO_ComputePatchForceMomCo





! *******************************************************************************
!
! Purpose: Compute global force and moment coefficients.
!
! Description: none.
!
! Input: global = global data
!
! Output: global%forceCoeffs, global%momentCoeffs = 
!         global force and moment coefficients computed
!
! Notes: none.
!
!*******************************************************************************

SUBROUTINE RFLO_ComputeIntegralForceMomCo( global )

  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... loop variables

! ... local variables
  REAL(RFREAL) :: fFact, mFact
  REAL(RFREAL), DIMENSION(12) :: globalVals, localVals

! *****************************************************************************
    
  CALL RegisterFunction(global,'RFLO_ComputeIntegralForceMomCo',&
       'RFLO_ModForcesMoments.F90')

! Get constants ---------------------------------------------------------------

  fFact = 1.0_RFREAL/global%forceRefArea
  mFact = 1.0_RFREAL/(global%forceRefArea*global%forceRefLength)
  
! Gather data -----------------------------------------------------------------

  localVals(1:3)   = global%forceCoeffs(XCOORD:ZCOORD,1)
  localVals(4:6)   = global%forceCoeffs(XCOORD:ZCOORD,2)
  localVals(7:9)   = global%momentCoeffs(XCOORD:ZCOORD,1)
  localVals(10:12) = global%momentCoeffs(XCOORD:ZCOORD,2)

! Perform reduction operation -------------------------------------------------

#ifdef MPI
  CALL MPI_AllReduce( localVals,globalVals,SIZE(localVals),MPI_RFREAL,MPI_SUM,&
                      global%mpiComm,global%mpierr )
  IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
       __LINE__ )
#else
  globalVals = localVals
#endif   

! Assign data -----------------------------------------------------------------

  global%forceCoeffs(:,1)  = fFact*globalVals(1:3)
  global%forceCoeffs(:,2)  = fFact*globalVals(4:6)
  global%momentCoeffs(:,1) = mFact*globalVals(7:9)  
  global%momentCoeffs(:,2) = mFact*globalVals(10:12)  

! Finalize --------------------------------------------------------------------

  CALL DeregisterFunction(global)

END SUBROUTINE RFLO_ComputeIntegralForceMomCo





! *******************************************************************************
!
! Purpose: Write time dependent global force and moment coefficients to file.
!
! Description: None.
!
! Input:  global = global force and moment coefficients
!
! Output: None.
!
! Notes: only master processor writes to file
!
! ******************************************************************************

  SUBROUTINE RFLO_WriteIntegralForceMomCo( global )

    IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... loop variables

! ... local variables
  CHARACTER(CHRLEN+9) :: fname

  INTEGER :: errorFlag
  REAL(RFREAL) :: fpCo(3), fvCo(3), mpCo(3), mvCo(3)

!******************************************************************************
    
  CALL RegisterFunction(global,'RFLO_WriteIntegralForceMomCo',&
       'RFLO_ModForcesMoments.F90')
                
! Write global force and moment coefficients ----------------------------------

  IF (global%myProcid==MASTERPROC) THEN

    fpCo(1:3) = global%forceCoeffs(XCOORD:ZCOORD,1)
    fvCo(1:3) = global%forceCoeffs(XCOORD:ZCOORD,2)
    mpCo(1:3) = global%momentCoeffs(XCOORD:ZCOORD,1)
    mvCo(1:3) = global%momentCoeffs(XCOORD:ZCOORD,2)

    IF (global%flowType == FLOW_STEADY) THEN
      WRITE(IF_FORMOM,1000,IOSTAT=errorFlag) global%currentIter, &
                                       fpCo(1:3),fvCo(1:3),mpCo(1:3),mvCo(1:3)
    ELSE
      WRITE(IF_FORMOM,1005,IOSTAT=errorFlag) global%currentTime, &
                                       fpCo(1:3),fvCo(1:3),mpCo(1:3),mvCo(1:3)
    ENDIF

    global%error = errorFlag
    IF (global%error /= 0) THEN
      CALL ErrorStop( global,ERR_FILE_WRITE,&
           __LINE__, &
                      'force-moment coeffs file' )
    ENDIF

! - Close and open .fom file (instead of fflush) ------------------------------

    IF (global%probeOpenClose) THEN
      WRITE(fname,'(A)') TRIM(global%outDir)//TRIM(global%casename)//'.fom'
      CLOSE(IF_FORMOM)
      OPEN( IF_FORMOM,FILE=fname,FORM='FORMATTED',STATUS='OLD', &
            POSITION='APPEND')
    ENDIF

  ENDIF   ! only masterproc

! Finalize --------------------------------------------------------------------

  CALL DeregisterFunction(global)

1000 FORMAT(I6,12(1PE13.5))
1005 FORMAT(1PE14.7,12(1PE13.5))

END SUBROUTINE RFLO_WriteIntegralForceMomCo



!******************************************************************************
!
! Purpose: open file for force and moment coefficients, write header.
!
! Description: none.
!
! Input: global  = global force/moment coeffs, flow type, restart.
!
! Output: to file.
!
! Notes:
!      This routine should only be called for the master process. It will not 
!      do anything if called by any other process (safeguard).
!
!******************************************************************************

SUBROUTINE RFLO_OpenForceMomCoFile( global )

  USE ModTools, ONLY: FloatLess
  IMPLICIT NONE

! ... parameters
  TYPE (t_global), POINTER :: global

! ... loop variables

! ... local variables
  CHARACTER(CHRLEN+9) :: fname

  LOGICAL :: fileExists, fileAppend
  INTEGER :: iFile, lastIter, errorFlag
  REAL(RFREAL) :: lastTime

!******************************************************************************

  CALL RegisterFunction( global,'RFLO_OpenForceMomCoFile',&
       'RFLO_ModForcesMoments.F90' )

! store id and generate file name ---------------------------------------------

  iFile = IF_FORMOM

  WRITE(fname,'(A)') TRIM(global%outDir)//TRIM(global%casename)//'.fom'

! append to existing file (restart) or create new file ------------------------

  IF ((global%flowType==FLOW_UNSTEADY .AND. &
       global%currentTime>0._RFREAL)  .OR.  &
      (global%flowType==FLOW_STEADY   .AND. &
       global%currentIter>1)) THEN

    INQUIRE( FILE=fname,EXIST=fileExists )
    IF (fileExists) THEN
      fileAppend = .TRUE.
      PRINT *,SOLVER_NAME,' Appending to ',TRIM(fname)
      OPEN(iFile,FILE=fname,FORM='formatted',STATUS='old',POSITION='append', &
           IOSTAT=errorFlag)
    ELSE
      fileAppend = .FALSE.
      PRINT *,SOLVER_NAME,' Overwriting ',TRIM(fname)
      OPEN(iFile,FILE=fname,FORM='formatted',STATUS='unknown',IOSTAT=errorFlag)
    ENDIF
  ELSE  ! new file

    fileAppend = .FALSE.
    PRINT *,SOLVER_NAME,' Creating new ',TRIM(fname)
    OPEN(iFile,FILE=fname,FORM='formatted',STATUS='unknown',IOSTAT=errorFlag)
  ENDIF ! fileStatus

  global%error = errorFlag
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,&
    __LINE__,'File: '//TRIM(fname) )

! write header ----------------------------------------------------------------

  IF (.NOT. fileAppend) THEN

    IF (global%flowType==FLOW_STEADY) THEN
      WRITE(iFile,1000,IOSTAT=errorFlag)
    ELSE
      WRITE(iFile,1005,IOSTAT=errorFlag)
    ENDIF

    global%error = errorFlag
    IF (global%error /= ERR_NONE) &
      CALL ErrorStop( global,ERR_FILE_WRITE,&
      __LINE__,'File: '//TRIM(fname) )
  ELSE

! - read the last line to get the last dump time.  If the current time is
!   earlier, back up to a time prior to the current one.

    IF ( global%flowType == FLOW_UNSTEADY ) THEN
      lastTime = HUGE( 1.0_RFREAL )

      LoopUnsteady: DO 
        BACKSPACE( iFile,IOSTAT=errorFlag )
        IF (errorFlag /= ERR_NONE) EXIT LoopUnsteady
        READ(iFile, FMT=*, IOSTAT=errorFlag) lastTime
        IF (errorFlag /= ERR_NONE) EXIT LoopUnsteady
        IF (FLOATLESS( lastTime,global%currentTime )) THEN
          EXIT LoopUnsteady
        ELSE
          BACKSPACE( iFile,IOSTAT=errorFlag )
          IF (errorFlag /= ERR_NONE) EXIT LoopUnsteady
        ENDIF
      ENDDO LoopUnsteady

      PRINT *,SOLVER_NAME,' positioned ',TRIM(fname),' at time ',lastTime

    ELSE
      lastIter = HUGE( 1 )

      LoopSteady: DO 
        BACKSPACE( iFile,IOSTAT=errorFlag )
        IF (errorFlag /= ERR_NONE) EXIT LoopSteady
        READ(iFile, FMT=*, IOSTAT=errorFlag) lastIter
        IF (errorFlag /= ERR_NONE) EXIT LoopSteady
        IF (lastIter < global%currentIter) THEN
          EXIT LoopSteady
        ELSE
          BACKSPACE( iFile,IOSTAT=errorFlag )
          IF (errorFlag /= ERR_NONE) EXIT LoopSteady
        ENDIF
      ENDDO LoopSteady

      PRINT *,SOLVER_NAME,' positioned ',TRIM(fname),' at iteration ',lastIter

    ENDIF  ! flowType
  ENDIF    ! append or new

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1000 FORMAT('iter, fpCoef(X,Y,Z), fvCoef(X,Y,Z), mpCoef(X,Y,Z), mvCoef(X,Y,Z)')
1005 FORMAT('time, fpCoef(X,Y,Z), fvCoef(X,Y,Z), mpCoef(X,Y,Z), mvCoef(X,Y,Z)')

END SUBROUTINE RFLO_OpenForceMomCoFile



! *******************************************************************************
!
! Purpose: Find region and patches contributing to global force and moment 
!          coefficients.
!
! Description: none.
!
! Input: regions = global input and patches%bcType in all regions
!
! Output: patch%globalAeroCoeffs==.true. = contributing regions and patches
!
! Notes: the search and write to file hereafter is done in parallel manner
!        to mimic the actual process in computing the coefficients, hence
!        minimizing generation of incorrect info.
!
!*******************************************************************************

SUBROUTINE RFLO_FindPatchCoeffsGlo( regions )

  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset
  IMPLICIT NONE

#include "Indexing.h"
        
! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch, i, j, k

! ... local variables
  INTEGER :: iLev, bcType, ijkN, iNOff, ijNOff
  INTEGER :: inbeg, inend, jnbeg, jnend, knbeg, knend
  REAL(RFREAL) :: boxXmin, boxXmax, boxYmin, boxYmax, boxZmin, boxZmax
  REAL(RFREAL) :: xMin, xMax, yMin, yMax, zMin, zMax
  REAL(RFREAL), POINTER :: xyz(:,:)
  TYPE (t_global), POINTER :: global
  TYPE (t_patch),  POINTER :: patch

! ******************************************************************************
    
  global => regions(1)%global
  CALL RegisterFunction(global,'RFLO_FindPatchCoeffsGlo',&
       'RFLO_ModForcesMoments.F90')

! get needed parameters --------------------------------------------------------

  iLev = 1

  boxXmin = global%acBndBoxXmin
  boxXmax = global%acBndBoxXmax
  boxYmin = global%acBndBoxYmin
  boxYmax = global%acBndBoxYmax
  boxZmin = global%acBndBoxZmin
  boxZmax = global%acBndBoxZmax

! loop over all regions and patches in parallel --------------------------------

  DO iReg = 1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor


      CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )
      xyz => regions(iReg)%levels(iLev)%grid%xyz
        
      DO iPatch = 1,regions(iReg)%nPatches
        patch  => regions(iReg)%levels(iLev)%patches(iPatch)
        bcType =  patch%bcType

! ----- get patch dimensions

        CALL RFLO_GetPatchIndicesNodes( regions(iReg),patch,iLev, &
                                        inbeg,inend,jnbeg,jnend,knbeg,knend )

! ----- search for min/max xyz

        xMin =  HUGE( 1._RFREAL )
        xMax = -HUGE( 1._RFREAL )
        yMin =  HUGE( 1._RFREAL )
        yMax = -HUGE( 1._RFREAL )
        zMin =  HUGE( 1._RFREAL )
        zMax = -HUGE( 1._RFREAL )

        DO k=knbeg,knend
          DO j=jnbeg,jnend
            DO i=inbeg,inend
              ijkN = IndIJK(i,j,k,iNOff,ijNOff)
              xMin = MIN( xMin,xyz(XCOORD,ijkN) )
              xMax = MAX( xMax,xyz(XCOORD,ijkN) )
              yMin = MIN( yMin,xyz(YCOORD,ijkN) )
              yMax = MAX( yMax,xyz(YCOORD,ijkN) )
              zMin = MIN( zMin,xyz(ZCOORD,ijkN) )
              zMax = MAX( zMax,xyz(ZCOORD,ijkN) )
            ENDDO
          ENDDO
        ENDDO

! ----- test for globalAeroCoeffs

        IF (((bcType>=BC_SLIPWALL   .AND. bcType<=BC_SLIPWALL+BC_RANGE)   .OR. &
             (bcType>=BC_NOSLIPWALL .AND. bcType<=BC_NOSLIPWALL+BC_RANGE) .OR. &
             (bcType>=BC_INJECTION  .AND. bcType<=BC_INJECTION+BC_RANGE)) .AND.&
            ((xMin > boxXmin .AND. xMax < boxXmax) .AND. &
             (yMin > boxYmin .AND. yMax < boxYmax) .AND. & 
             (zMin > boxZmin .AND. zMax < boxZmax))) THEN

          patch%globalAeroCoeffs = .TRUE.

        ENDIF ! bcType and boxBnd
      ENDDO   ! iPatch

    ENDIF     ! myProcid
  ENDDO       ! iReg

! finalize ---------------------------------------------------------------------

  CALL DeregisterFunction(global)

END SUBROUTINE RFLO_FindPatchCoeffsGlo



! ******************************************************************************
!
! Purpose: write location of region and patches where global force and moment
!          coefficients are defined.
!
! Description: only MASTERPROC write the information.
!
! Notes: none.
!
! ******************************************************************************

SUBROUTINE RFLO_WritePatchCoeffsInfo( regions )

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iPatch
   
! ... local variables
  CHARACTER(CHRLEN) :: fname
  TYPE(t_global), POINTER :: global
  TYPE(t_patch),  POINTER :: patch

#ifdef MPI
  INTEGER :: status(MPI_STATUS_SIZE)
#endif
  INTEGER :: iLev, errorFlag, iFile, tag
  
! ******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLO_WritePatchCoeffsInfo',&
       'RFLO_ModForcesMoments.F90')

! open file --------------------------------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN

    iFile = IF_PATCH_COEF
    WRITE(fname,'(A)') TRIM(global%outDir)//TRIM(global%casename)//'.pcoin'

    OPEN(iFile,FILE=fname,FORM="formatted",STATUS="unknown",IOSTAT=errorFlag)   
    global%error = errorFlag        

    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop(global,ERR_FILE_OPEN,&
      __LINE__,fname)

  ENDIF  ! MASTERPROC

! header and general information ---------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    WRITE(iFile,1005)' Regions and patches of force and moment coefficients'  
    WRITE(iFile,1000)
    WRITE(iFile,1005)' region-number    patch-number'  
    WRITE(iFile,1000)
  ENDIF

! print region and patches of global aero coeffs -----------------------------

  DO iReg=1,global%nRegions
    iLev  =  regions(iReg)%currLevel

    DO iPatch=1,regions(iReg)%nPatches
      patch => regions(iReg)%levels(iLev)%patches(iPatch)
      tag   =  regions(iReg)%localNumber + MPI_PATCHOFF*iPatch

! --- master receives and writes data, others send them

      IF (global%myProcid == MASTERPROC) THEN
#ifdef MPI
        IF (regions(iReg)%procid /= MASTERPROC) THEN
          CALL MPI_Recv( patch%globalAeroCoeffs,1,MPI_LOGICAL, &
                         regions(iReg)%procid,tag,global%mpiComm,status, &
                         global%mpierr )
          IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
               __LINE__ )
        ENDIF
#endif
        IF (patch%globalAeroCoeffs) WRITE(iFile,'(2I10)') iReg, iPatch

      ELSE  ! not the master
#ifdef MPI
        IF (regions(iReg)%procid == global%myProcid) THEN
          CALL MPI_Send( patch%globalAeroCoeffs,1,MPI_LOGICAL,MASTERPROC,tag, &
                         global%mpiComm,global%mpierr )
          IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,&
               __LINE__ )
        ENDIF
#endif
      ENDIF

    ENDDO ! iPatch
  ENDDO   ! iReg

! close file -------------------------------------------------------------------

  IF (global%myProcid == MASTERPROC) THEN
    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) &
      CALL ErrorStop(global,ERR_FILE_CLOSE,&
      __LINE__,fname)
  ENDIF

! finalize ---------------------------------------------------------------------

  CALL DeregisterFunction( global )

1000 FORMAT(/,1X,40('-'))
1005 FORMAT(/,A)
  
END SUBROUTINE RFLO_WritePatchCoeffsInfo

! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLO_ModForcesMoments

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_ModForcesMoments.F90,v $
! Revision 1.6  2008/12/06 08:44:16  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:27  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/03/25 04:34:49  wasistho
! only masterproc write .pcoin and .fom files
!
! Revision 1.3  2006/03/25 01:14:46  wasistho
! cosmetics
!
! Revision 1.2  2006/03/24 23:32:47  wasistho
! added RFLO_FindPatchCoeffsGlo
!
! Revision 1.1  2006/03/24 05:05:51  wasistho
! initial import
!
!
! ******************************************************************************












