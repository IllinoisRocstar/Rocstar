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
! Purpose: write grid (and solution) data to plot file in ASCII format.
!
! Description: none.
!
! Input: iReg     = region number
!        iLev     = grid level
!        plotType = 1 - grid + smoke only, = 2 - grid + all fields
!        region   = region data (dimensions, flow variables)
!
! Output: to plot file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: WriteTecplotAscii.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE WriteTecplotAscii( iReg,iLev,plotType,region )

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE RFLO_ModInterfacesPost, ONLY : RFLO_GetDimensPhysNodes, &
                                     RFLO_GetNodeOffset, &
                                     RFLO_GetCellOffset, Aver, AverDiv
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: iReg, iLev, plotType

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i, j, k, iPeul

! ... local variables
  CHARACTER(CHRLEN+4) :: fname
  CHARACTER(256)      :: varStr
  CHARACTER(16)       :: concStr

  INTEGER :: ipnbeg, ipnend, jpnbeg, jpnend, kpnbeg, kpnend
  INTEGER :: iCOff, ijCOff, iNOff, ijNOff, ijkN, cell(8), errorFlag, nPeul

  REAL(RFREAL)          :: rho, u, v, w, press, temp, mach, c, qq
  REAL(RFREAL), POINTER :: xyz(:,:), cv(:,:), dv(:,:), tv(:,:), gv(:,:)
  REAL(RFREAL), POINTER :: peulCv(:,:), conc(:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'WriteTecplotAscii',&
  'WriteTecplotAscii.F90' )

! set parameters --------------------------------------------------------------

  CALL RFLO_GetDimensPhysNodes( region,iLev,ipnbeg,ipnend, &
                                jpnbeg,jpnend,kpnbeg,kpnend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  nPeul = region%levels(iLev)%peul%nCv

  IF (nPeul > 0) ALLOCATE( conc(nPeul),stat=errorFlag )

  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! open file and write the header ----------------------------------------------

  IF (iReg == 1) THEN
    WRITE(fname,'(A,ES11.5,A)') &
       TRIM(global%casename)//'.peul_',global%currentTime,'.dat'
    OPEN(IF_PLOT,FILE=fname,status='unknown',form='formatted',iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,fname )

    IF (global%flowType == FLOW_STEADY) THEN
      WRITE(IF_PLOT,1000,err=10) TRIM(global%casename),global%currentIter
    ELSE
      IF (global%currentTime <= 0._RFREAL) THEN
        WRITE(IF_PLOT,1005,err=10) TRIM(global%casename),global%timeStamp
      ELSE
        WRITE(IF_PLOT,1005,err=10) TRIM(global%casename),global%currentTime
      ENDIF
    ENDIF

    varStr = 'x y z'
    IF (plotType == PLOT_GRID_FLOW) varStr = 'x y z rho u v w p T M'

    DO iPeul=1,nPeul
      SELECT CASE (iPeul)
      CASE ( 0: 9)
        WRITE(concStr,'(A,I1)') 'c_', iPeul
      CASE (10:99)
        WRITE(concStr,'(A,I2)') 'c_', iPeul
      CASE DEFAULT
        WRITE(concStr,'(A)') 'c_?'
      END SELECT ! iPeul
      varStr = TRIM(varStr)//' '//TRIM(concStr)
    ENDDO ! iPeul

    WRITE(IF_PLOT,1010,err=10) TRIM(varStr)

  ENDIF   ! iReg=1

! write zone header

  WRITE(IF_PLOT,1015) iReg,ipnend-ipnbeg+1,jpnend-jpnbeg+1,kpnend-kpnbeg+1

! write data ------------------------------------------------------------------
! pointer to variables

  xyz => region%levels(iLev)%grid%xyz
  cv  => region%levels(iLev)%mixt%cv
  dv  => region%levels(iLev)%mixt%dv
  tv  => region%levels(iLev)%mixt%tv
  gv  => region%levels(iLev)%mixt%gv
  IF (nPeul > 0) peulCv => region%levels(iLev)%peul%cv

  DO k=kpnbeg,kpnend
    DO j=jpnbeg,jpnend
      DO i=ipnbeg,ipnend
        ijkN = IndIJK(i,j,k,iNOff,ijNOff)
        cell(1) = IndIJK(i  ,j  ,k  ,iCOff,ijCOff)
        cell(2) = IndIJK(i-1,j  ,k  ,iCOff,ijCOff)
        cell(3) = IndIJK(i  ,j-1,k  ,iCOff,ijCOff)
        cell(4) = IndIJK(i-1,j-1,k  ,iCOff,ijCOff)
        cell(5) = IndIJK(i  ,j  ,k-1,iCOff,ijCOff)
        cell(6) = IndIJK(i-1,j  ,k-1,iCOff,ijCOff)
        cell(7) = IndIJK(i  ,j-1,k-1,iCOff,ijCOff)
        cell(8) = IndIJK(i-1,j-1,k-1,iCOff,ijCOff)

        DO iPeul=1,nPeul
          conc(iPeul) = Aver(cell,iPeul,peulCv)
        ENDDO

        IF (plotType == PLOT_GRID_ONLY) THEN

! ------- write grid coordinates and smoke fields only

          WRITE(IF_PLOT,1020,err=10) xyz(XCOORD,ijkN), &
                                     xyz(YCOORD,ijkN), &
                                     xyz(ZCOORD,ijkN), &
                                     (conc(iPeul),iPeul=1,nPeul)

        ELSE

! ------- write grid and all fields

          rho   = Aver(cell,CV_MIXT_DENS,cv)
          u     = AverDiv(cell,CV_MIXT_XMOM,cv,CV_MIXT_DENS,cv)
          v     = AverDiv(cell,CV_MIXT_YMOM,cv,CV_MIXT_DENS,cv)
          w     = AverDiv(cell,CV_MIXT_ZMOM,cv,CV_MIXT_DENS,cv)
          press = Aver(cell,DV_MIXT_PRES,dv)
          temp  = Aver(cell,DV_MIXT_TEMP,dv)
          c     = Aver(cell,DV_MIXT_SOUN,dv)
          qq    = u*u + v*v + w*w
          mach  = SQRT(qq)/c

          WRITE(IF_PLOT,1020,err=10) xyz(XCOORD,ijkN), &
                                     xyz(YCOORD,ijkN), &
                                     xyz(ZCOORD,ijkN), &
                                     rho,u,v,w,press,temp,mach, &
                                     (conc(iPeul),iPeul=1,nPeul)

        ENDIF ! plotType

      ENDDO
    ENDDO
  ENDDO

! close file, handle errors ---------------------------------------------------

  IF (iReg == global%nRegions) THEN
    CLOSE(IF_PLOT,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,fname )
  ENDIF

  IF (nPeul > 0) DEALLOCATE( conc,stat=errorFlag )

  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_WRITE,__LINE__,fname )

! formats ---------------------------------------------------------------------

1000 FORMAT('TITLE="',A,'. Iteration: ',I8,'."')
1005 FORMAT('TITLE="',A,'. Time: ',ES11.5,'."')
1010 FORMAT('VARIABLES= ',A)
1015 FORMAT('ZONE T="',I5.5,'", I=',I6,', J=',I6,', K=',I6,', F=POINT')
1020 FORMAT(999(1X,ES13.6))

999  CONTINUE
END SUBROUTINE WriteTecplotAscii

!******************************************************************************
!
! RCS Revision history:
!
! $Log: WriteTecplotAscii.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:29:09  fnajjar
! Initial revision after changing case
!
! Revision 1.2  2004/03/02 21:50:57  jferry
! Added timestamp to name of peulpost output file
!
! Revision 1.1  2003/09/25 15:40:22  jferry
! Implented Rocsmoke post-processing
!
!******************************************************************************







