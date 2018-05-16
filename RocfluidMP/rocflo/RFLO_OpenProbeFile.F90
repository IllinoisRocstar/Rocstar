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
! Purpose: open file(s) for probe data, write header(s).
!
! Description: none.
!
! Input: regions = region data (processor number, active flag)
!        global  = probe location, flow type, restart.
!
! Output: to file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_OpenProbeFile.F90,v 1.7 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_OpenProbeFile( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, CentroidHexa
  USE ModError
  USE ModParameters
  USE ModMPI
  USE ModTools, ONLY: FloatLess
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE (t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iprobe

! ... local variables
  CHARACTER(CHRLEN+9) :: fname

  INTEGER :: iReg, iLev, i, j, k, iNOff, ijNOff, errorFlag
  INTEGER :: corner(8)
  INTEGER :: ipcbeg, ipcend, jpcbeg, jpcend, kpcbeg, kpcend
  INTEGER :: iCoff, ijCOff
  INTEGER :: iprobemax
  INTEGER :: probeIter

  LOGICAL :: fileExists, fileAppend

  REAL(RFREAL)          :: xc, yc, zc, xyzHexa(3,8)
  REAL(RFREAL), POINTER :: xyz(:,:)
  REAL(RFREAL)          :: xmin, xmax, ymin, ymax, zmin, zmax
  REAL(RFREAL)          :: xmn, xmx, ymn, ymx, zmn, zmx
  REAL(RFREAL)          :: probeTime

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_OpenProbeFile',&
  'RFLO_OpenProbeFile.F90' )

! open file

  IF (global%nProbes > 0) THEN
    DO iprobe=1,global%nProbes

! --- check if region`s number within range

      IF (global%probePos(iprobe,1)>global%nRegions) &
        CALL ErrorStop( global,ERR_PROBE_LOCATION,__LINE__ )

! --- support entering 0 and coordinates

      IF (global%probePos(iprobe,1)<1) THEN

        IF (global%myProcid == 0 .AND. &
            global%verbLevel >= VERBOSE_MED) THEN
          WRITE(STDOUT,'(A,1X,A,I1,A)')SOLVER_NAME,&
                         'Trying to find cell containing probe '&
                         ,iProbe,' at'
          WRITE(STDOUT,'(A,1X,3(E14.5))')&
                         SOLVER_NAME,global%probeXYZ(iprobe,2),&
                         global%probeXYZ(iprobe,3), &
                         global%probeXYZ(iprobe,4)
        ENDIF

regs:   DO iReg = 1,global%nRegions

          IF (regions(iReg)%procid==global%myProcid .AND. &
              regions(iReg)%active==ACTIVE) THEN

            iLev = regions(iReg)%currLevel
            CALL RFLO_GetDimensPhys( regions(iReg),iLev,ipcbeg,ipcend, &
                                     jpcbeg,jpcend,kpcbeg,kpcend )
            CALL RFLO_GetCellOffset( regions(iReg),iLev,iCOff,ijCOff )
            CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )

            xyz  => regions(iReg)%levels(iLev)%grid%xyz

! --------- see if x, y, and z is within bounding box of any of my cells

            xmn = 1.0e+30
            xmx = -1.0e+30
            ymn = 1.0e+30
            ymx = -1.0e+30
            zmn = 1.0e+30
            zmx = -1.0e+30
            DO k=kpcbeg,kpcend
              DO j=jpcbeg,jpcend
                DO i=ipcbeg,ipcend
                  corner(1) = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
                  corner(2) = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
                  corner(3) = IndIJK(i  ,j+1,k+1,iNOff,ijNOff)
                  corner(4) = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
                  corner(5) = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
                  corner(6) = IndIJK(i+1,j  ,k+1,iNOff,ijNOff)
                  corner(7) = IndIJK(i+1,j+1,k+1,iNOff,ijNOff)
                  corner(8) = IndIJK(i+1,j+1,k  ,iNOff,ijNOff)

                  xmin = MIN(xyz(1,corner(1)),xyz(1,corner(2)),xyz(1,corner(3)), &
                             xyz(1,corner(4)),xyz(1,corner(5)),xyz(1,corner(6)), &
                             xyz(1,corner(7)),xyz(1,corner(8)))
                  xmax = MAX(xyz(1,corner(1)),xyz(1,corner(2)),xyz(1,corner(3)), &
                             xyz(1,corner(4)),xyz(1,corner(5)),xyz(1,corner(6)), &
                             xyz(1,corner(7)),xyz(1,corner(8)))
                  ymin = MIN(xyz(2,corner(1)),xyz(2,corner(2)),xyz(2,corner(3)), &
                             xyz(2,corner(4)),xyz(2,corner(5)),xyz(2,corner(6)), &
                             xyz(2,corner(7)),xyz(2,corner(8)))
                  ymax = MAX(xyz(2,corner(1)),xyz(2,corner(2)),xyz(2,corner(3)), &
                             xyz(2,corner(4)),xyz(2,corner(5)),xyz(2,corner(6)), &
                             xyz(2,corner(7)),xyz(2,corner(8)))
                  zmin = MIN(xyz(3,corner(1)),xyz(3,corner(2)),xyz(3,corner(3)), &
                             xyz(3,corner(4)),xyz(3,corner(5)),xyz(3,corner(6)), &
                             xyz(3,corner(7)),xyz(3,corner(8)))
                  zmax = MAX(xyz(3,corner(1)),xyz(3,corner(2)),xyz(3,corner(3)), &
                             xyz(3,corner(4)),xyz(3,corner(5)),xyz(3,corner(6)), &
                             xyz(3,corner(7)),xyz(3,corner(8)))

                  IF (xmin < xmn) xmn = xmin
                  IF (xmax > xmx) xmx = xmax
                  IF (ymin < ymn) ymn = ymin
                  IF (ymax > ymx) ymx = ymax
                  IF (zmin < zmn) zmn = zmin
                  IF (zmax > zmx) zmx = zmax

                  IF ((xmin <= global%probeXYZ(iprobe,2)) .AND.  &
                      (xmax >= global%probeXYZ(iprobe,2)) .AND.  &
                      (ymin <= global%probeXYZ(iprobe,3)) .AND.  &
                      (ymax >= global%probeXYZ(iprobe,3)) .AND.  &
                      (zmin <= global%probeXYZ(iprobe,4)) .AND.  &
                      (zmax >= global%probeXYZ(iprobe,4)) ) THEN

! ----------------- yes; assign block number and cell indices to probePos

                    IF (global%verbLevel >= VERBOSE_MED) THEN
                      WRITE(STDOUT,'(A,1X,A,2(I2,A),3(I2),A)')&
                           SOLVER_NAME,'Found probe ',iprobe, &
                           ' in block ',iReg,' cell ',i,j,k,&
                           ' with bounding box'
                      WRITE(STDOUT,'(A,1X,6(E14.5))')&
                            SOLVER_NAME, xmin, xmax, ymin, ymax,&
                            zmin, zmax
                    ENDIF
                    global%probePos(iprobe,1) = iReg
                    global%probePos(iprobe,2) = i
                    global%probePos(iprobe,3) = j
                    global%probePos(iprobe,4) = k

! ----------------- jump to named loop to stop searching over regions

                    EXIT regs
                  ENDIF
                ENDDO  ! i
              ENDDO    ! j
            ENDDO      ! k

          ENDIF     ! is this region on my proc
        ENDDO regs  ! iReg

#ifdef MPI
! ----- must broadcast block with probe to the others.  Use
!       MAX probePos(iprobe,1); eliminates redundant finds.
!       Other do not need the 2:4 elements.

        CALL MPI_Allreduce(global%probePos(iprobe,1),iprobemax,1,MPI_INTEGER,  &
                           MPI_MAX,global%mpiComm,global%mpierr)
        IF (global%mpierr /= 0) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
        global%probePos(iprobe,1) = iprobemax
#endif

      ENDIF  ! probePOS(iprobe,1) lower than range

! --- check if probe`s region on current processor

      IF (regions(global%probePos(iprobe,1))%procid==global%myProcid .AND. &
          regions(global%probePos(iprobe,1))%active==ACTIVE) THEN

! ----- generate file name

        WRITE(fname,'(A,I4.4)') TRIM(global%outDir)//TRIM(global%casename)//'.prb_',iprobe

! ----- append to existing file (restart) or create new file

        IF ((global%flowType==FLOW_UNSTEADY .AND. &
             global%currentTime>0._RFREAL)  .OR.  &
            (global%flowType==FLOW_STEADY   .AND. global%currentIter>1)) THEN
          INQUIRE(file=fname,exist=fileExists)
          IF (fileExists) THEN
            fileAppend = .TRUE.
            IF (global%verbLevel >= VERBOSE_MED) THEN
              PRINT *,SOLVER_NAME,' Appending to ',TRIM(fname)
            ENDIF
            OPEN(IF_PROBE+iprobe-1,file=fname,form='formatted',status='old', &
                                   position='append',iostat=errorFlag)
          ELSE
            fileAppend = .FALSE.
            IF (global%verbLevel >= VERBOSE_MED) THEN
              PRINT *,SOLVER_NAME,' Overwriting ',TRIM(fname)
            ENDIF
            OPEN(IF_PROBE+iprobe-1,file=fname,form='formatted', &
                                   status='unknown',iostat=errorFlag)
          ENDIF
        ELSE
          fileAppend = .FALSE.
          IF (global%verbLevel >= VERBOSE_MED) THEN
            WRITE(STDOUT,'(A,A,A)')SOLVER_NAME,' Creating new ',TRIM(fname)
          ENDIF
          OPEN(IF_PROBE+iprobe-1,file=fname,form='formatted',status='unknown', &
                                 iostat=errorFlag)
        ENDIF
        global%error = errorFlag
        IF (global%error /= 0) &
          CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname) )

! ----- write header ...

        iReg = global%probePos(iprobe,1)
        iLev = regions(iReg)%currLevel
        xyz => regions(iReg)%levels(iLev)%grid%xyz

        CALL RFLO_GetNodeOffset( regions(iReg),iLev,iNOff,ijNOff )

        i         = global%probePos(iprobe,2)
        j         = global%probePos(iprobe,3)
        k         = global%probePos(iprobe,4)
        corner(1) = IndIJK(i  ,j  ,k  ,iNOff,ijNOff)
        corner(2) = IndIJK(i  ,j  ,k+1,iNOff,ijNOff)
        corner(3) = IndIJK(i  ,j+1,k+1,iNOff,ijNOff)
        corner(4) = IndIJK(i  ,j+1,k  ,iNOff,ijNOff)
        corner(5) = IndIJK(i+1,j  ,k  ,iNOff,ijNOff)
        corner(6) = IndIJK(i+1,j  ,k+1,iNOff,ijNOff)
        corner(7) = IndIJK(i+1,j+1,k+1,iNOff,ijNOff)
        corner(8) = IndIJK(i+1,j+1,k  ,iNOff,ijNOff)

        xyzHexa(1:3,1) = xyz(1:3,corner(1))
        xyzHexa(1:3,2) = xyz(1:3,corner(2))
        xyzHexa(1:3,3) = xyz(1:3,corner(3))
        xyzHexa(1:3,4) = xyz(1:3,corner(4))
        xyzHexa(1:3,5) = xyz(1:3,corner(5))
        xyzHexa(1:3,6) = xyz(1:3,corner(6))
        xyzHexa(1:3,7) = xyz(1:3,corner(7))
        xyzHexa(1:3,8) = xyz(1:3,corner(8))

        CALL CentroidHexa( xyzHexa,xc,yc,zc )

! ----- only if we created a new probe file (not appending to an existing one)

        IF (.NOT. fileAppend) THEN
          WRITE(IF_PROBE+iprobe-1,1000,iostat=errorFlag) &
                global%probePos(iprobe,1),global%probePos(iprobe,2), &
                global%probePos(iprobe,3),global%probePos(iprobe,4), &
                xc,yc,zc
          global%error = errorFlag
          IF (global%error /= ERR_NONE) &
            CALL ErrorStop( global,ERR_FILE_WRITE,__LINE__,'File: '//TRIM(fname) )
        ELSE

! ----- read the last line to get the last probe dump time.  If the initial
! ----- time is earlier, back up to a time prior to the initial one.

          IF ( global%flowType == FLOW_UNSTEADY ) THEN
            probeTime = HUGE(1.0_RFREAL)
            LoopUnsteady: DO 
              BACKSPACE (IF_PROBE+iprobe-1,IOSTAT=errorFlag)
              IF (errorFlag /= ERR_NONE) EXIT LoopUnsteady
              READ (IF_PROBE+iprobe-1,FMT=*,iostat=errorFlag) probeTime
              IF (errorFlag /= ERR_NONE) EXIT LoopUnsteady
              IF (FloatLess(probeTime,global%currentTime)) THEN
                EXIT LoopUnsteady
              ELSE
                BACKSPACE (IF_PROBE+iprobe-1,IOSTAT=errorFlag)
                IF (errorFlag /= ERR_NONE) EXIT LoopUnsteady
              ENDIF
            ENDDO LoopUnsteady
            IF (global%verbLevel >= VERBOSE_MED) THEN
              PRINT *,SOLVER_NAME,' positioned ',TRIM(fname),' at time ',probeTime
            ENDIF
          ELSE
            probeIter = HUGE(probeIter)
            LoopSteady: DO 
              BACKSPACE (IF_PROBE+iprobe-1,IOSTAT=errorFlag)
              IF (errorFlag /= ERR_NONE) EXIT LoopSteady
              READ (IF_PROBE+iprobe-1,FMT=*,iostat=errorFlag) probeIter
              IF (errorFlag /= ERR_NONE) EXIT LoopSteady
              IF (probeIter < global%currentIter) THEN
                EXIT LoopSteady
              ELSE
                BACKSPACE (IF_PROBE+iprobe-1,IOSTAT=errorFlag)
                IF (errorFlag /= ERR_NONE) EXIT LoopSteady
              ENDIF
            ENDDO LoopSteady
          ENDIF

        ENDIF ! Append or new

      ENDIF   ! probe located within region

    ENDDO     ! iprobe
  ENDIF       ! nProbes > 0

! finalize

  CALL DeregisterFunction( global )

! format

1000 FORMAT('# probe data (iteration/time, density, u, v, w, p, T)',/, &
            '# region ',I5,', icell ',I5,', jcell ',I5,', kcell ',I5,/, &
            '# x=',E13.5,', y=',E13.5,', z=',E13.5)

END SUBROUTINE RFLO_OpenProbeFile

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_OpenProbeFile.F90,v $
! Revision 1.7  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2005/06/28 21:58:08  wasistho
! set it back to currentTime
!
! Revision 1.4  2005/06/28 08:51:22  rfiedler
! Remove local currentTime; use timeStamp in place of currentTime to open probes.
!
! Revision 1.2  2005/02/25 01:50:41  rfiedler
! Probe files now not affected by restarting at an older dump time.  Less verbose.
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.14  2004/07/22 21:00:02  wasistho
! fixed missing ifdef MPI around MPI calls
!
! Revision 1.13  2004/07/21 21:11:59  wasistho
! allow probes input by coordinates
!
! Revision 1.12.2.2  2004/07/02 21:28:37  rfiedler
! Bug fix: use MPI_Allreduce to tell all processes who has the probes.  RAF
!
! Revision 1.12.2.1  2004/07/02 04:11:25  rfiedler
! Allows Rocflo probes to be specified by coordinates.  This routine finds the
! first region containing a cell whose bounding box contains the probe.  RAF
!
! Revision 1.12  2003/11/20 16:40:40  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.8  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.7  2003/02/06 01:22:29  jblazek
! Added check for presence of an old probe file.
!
! Revision 1.6  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.5  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.4  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.3  2002/06/22 01:13:38  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.2  2002/04/03 02:28:52  jblazek
! Added x,y,z location to probe file header.
!
! Revision 1.1  2002/02/25 22:36:53  jblazek
! Simplified solver initialization routine.
!
!******************************************************************************







