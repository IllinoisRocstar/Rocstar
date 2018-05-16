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
! Purpose: sort particle entering nozzle by bins.
!
! Description: none.
!
! Input: iReg    = region number
!        iLev    = current level
!        region  = region data (dimensions, plag variables)
!        iRegBin = region number to bin
!
! Output: to plot file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_BinSortSpatialDist.F90,v 1.3 2008/12/06 08:45:07 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_BinSortSpatialDist ( iReg, iLev, region, iRegBin )

  USE ModDataTypes
  USE ModError
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global

  USE ModMPI
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  INTEGER :: iReg, iLev, iRegBin

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i,j

! ... local variables
  CHARACTER(CHRLEN)   :: RCSIdentString
  CHARACTER(CHRLEN+4) :: fname

  INTEGER :: errorFlag, nDimPlag
  INTEGER :: pidini,regini,icell,indexi,indexj,indexk
  INTEGER :: iBin,nBins
  INTEGER :: iContAl,iContAlOx
  INTEGER, POINTER :: cvMass(:)
  INTEGER, POINTER :: aiv(:,:)
  INTEGER, DIMENSION(0:500) :: sizeBin

  REAL(RFREAL)          :: diamMicron,deltaY,xPosPlag,yPosPlag,yPosPlagShift
  REAL(RFREAL)          :: massL,yHeight,yHeightChan,yHeightNozz,yHeightMin
  REAL(RFREAL)          :: xMinRange
  REAL(RFREAL), DIMENSION(0:500) :: compAlBin,compAlOxBin,diam32Bin,diam43Bin,&
                                    diam2Bin,diam3Bin,diam4Bin,massAlBin,     &
				    massAlOxBin,massTotBin,uVelBin,vVelBin,   &
				    volBin,wVelBin,xMomBin,yBin,yMomBin,zMomBin
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_BinSortSpatialDist.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PLAG_BinSortSpatialDist', 'PLAG_BinSortSpatialDist.F90' )

! set parameters --------------------------------------------------------------

  nDimPlag = region%levels(iLev)%plag%nPcls
 
  sizeBin(0:500) = 0
  diam2Bin(0:500) = 0.0_RFREAL
  diam3Bin(0:500) = 0.0_RFREAL
  diam4Bin(0:500) = 0.0_RFREAL
  diam32Bin(0:500) = 0.0_RFREAL
  diam43Bin(0:500) = 0.0_RFREAL

  nBins=20
  
  iContAl   = 1
  iContAlOx = 2  
  compAlBin(0:500)   = 0.0_RFREAL
  compAlOxBin(0:500) = 0.0_RFREAL
  
  massAlBin(0:500)   = 0.0_RFREAL
  massAlOxBin(0:500) = 0.0_RFREAL
  massTotBin(0:500)  = 0.0_RFREAL
  volBin(0:500)      = 0.0_RFREAL
  
  uVelBin(0:500) = 0.0_RFREAL
  vVelBin(0:500) = 0.0_RFREAL
  wVelBin(0:500) = 0.0_RFREAL

  xMomBin(0:500) = 0.0_RFREAL
  ymomBin(0:500) = 0.0_RFREAL
  zmomBin(0:500) = 0.0_RFREAL

  IF (iReg /= iRegBin) GOTO 1999
  
  PRINT*,'PLAG_BinSortSpatialDist: iRegBin, nDimPlag = ',iRegBin,nDimPlag

! set values ------------------------------------------------------------------

! TEMPORARY: Data Pertinent to c1

  yHeightChan = 0.09_RFREAL
  yHeightNozz = 0.0481267_RFREAL
  
  SELECT CASE (iRegBin)
    CASE(41)
      yHeight = yHeightChan
      yHeightMin = 0.0_RFREAL
      xMinRange  = 0.0_RFREAL

    CASE(52)
      yHeight = yHeightNozz
      yHeightMin = 0.0209367_RFREAL
      xMinRange = 0.465_RFREAL

    CASE DEFAULT
      yHeight = 1.0_RFREAL
  END SELECT ! iRegBin
 
  deltaY = yHeight/REAL(nBins,KIND=RFREAL)
  yBin(1) = yHeightMin
  
  DO i = 2, nBins+1
    yBin(i) = yBin(i-1) +deltaY
  END DO ! i

! set pointers ----------------------------------------------------------------

  aiv => region%levels(iLev)%plag%aiv
  cv  => region%levels(iLev)%plag%cv
  dv  => region%levels(iLev)%plag%dv
  cvMass  => region%levels(iLev)%plag%cvPlagMass

! determine size ---------------------------------------------------------------

  DO i = 1, nDimPlag
    xPosPlag = cv(CV_PLAG_XPOS,i)
    
    IF ( xPosPlag >= xMinRange ) THEN
      yPosPlag = cv(CV_PLAG_YPOS,i)
      yPosPlagShift = yPosPlag -yHeightMin
    
      iBin=NINT(yPosPlagShift/deltaY)
      sizeBin(iBin) = sizeBin(iBin)+1
    
      diamMicron = dv(DV_PLAG_DIAM,i)*1.0E+06_RFREAL
      massL = SUM( cv(cvMass(:),i) )    
    
      massAlBin(iBin)   = massAlBin(iBin) +cv(cvMass(iContAl),i)
      massAlOxBin(iBin) = massAlOxBin(iBin) +cv(cvMass(iContAlOx),i)
      massTotBin(iBin)  = massTotBin(iBin) +massL
      diam2Bin(iBin)    = diam2Bin(iBin) +diamMicron**2
      diam3Bin(iBin)    = diam3Bin(iBin) +diamMicron**3
      diam4Bin(iBin)    = diam4Bin(iBin) +diamMicron**4

      xMomBin(iBin) = xMomBin(iBin) +cv(CV_PLAG_XMOM,i)
      yMomBin(iBin) = yMomBin(iBin) +cv(CV_PLAG_YMOM,i)
      zMomBin(iBin) = zMomBin(iBin) +cv(CV_PLAG_ZMOM,i)
    ENDIF ! xPosPlag
  ENDDO ! i 

  DO iBin = 1, nBins
    IF ( sizeBin(iBin) > 0 ) THEN
      compAlBin(iBin)   = massAlBin(iBin)/massTotBin(iBin)
      compAlOxBin(iBin) = massAlOxBin(iBin)/massTotBin(iBin)

      diam43Bin(iBin) = diam4Bin(iBin)/diam3Bin(iBin)
      diam32Bin(iBin) = diam3Bin(iBin)/diam2Bin(iBin)

      uVelBin(iBin) = xMomBin(iBin)/massTotBin(iBin)
      vVelBin(iBin) = yMomBin(iBin)/massTotBin(iBin)
      wVelBin(iBin) = zMomBin(iBin)/massTotBin(iBin)
    END IF ! sizeBin
  ENDDO ! iBin

! open file and write the header ----------------------------------------------

  WRITE(fname,'(A,I2,A,1PE11.5,A)') &
   TRIM(global%casename)//'.plag_ydist_bin_reg_',iRegBin,'_',global%currentTime,'.dat'
  OPEN(IF_PLOT,FILE=fname,status='unknown',form='formatted',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,fname )

  IF (global%currentTime <= 0._RFREAL) THEN
    WRITE(IF_PLOT,1005,err=10) TRIM(global%casename),global%timeStamp
  ELSE
    WRITE(IF_PLOT,1005,err=10) TRIM(global%casename),global%currentTime
  ENDIF

  WRITE(IF_PLOT,1010,err=10) 'iBin yBin size compAl compAlOx diam32 diam43 u v w'

  WRITE(IF_PLOT,1015) iRegBin, nBins

  DO iBin = 1, nBins
    WRITE(IF_PLOT,1020,err=10) iBin,(yBin(iBin)-yHeightMin)/yHeight,sizeBin(iBin),&
                               NINT(compAlBin(iBin)*100.0_RFREAL),&
			       NINT(compAlOxBin(iBin)*100.0_RFREAL), &
			       NINT(diam32Bin(iBin)),NINT(diam43Bin(iBin)),&
			       uVelBin(iBin),vVelBin(iBin),wVelBin(iBin)
  ENDDO ! iBin

! close file, handle errors ---------------------------------------------------

  CLOSE(IF_PLOT,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,fname )

! open file and write the header for time-averaging ---------------------------

  WRITE(fname,'(A,I2,A,1PE11.5,A)') &
   TRIM(global%casename)//'.plag_cv_ydist_bin_reg_',iRegBin,'_',global%currentTime,'.dat'
  OPEN(IF_PLOT,FILE=fname,status='unknown',form='formatted',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,fname )

  IF (global%currentTime <= 0._RFREAL) THEN
    WRITE(IF_PLOT,1005,err=10) TRIM(global%casename),global%timeStamp
  ELSE
    WRITE(IF_PLOT,1005,err=10) TRIM(global%casename),global%currentTime
  ENDIF

  WRITE(IF_PLOT,1010,err=10) 'iBin yBin size massAl massAlOxBin massTot diam2 diam3 diam4 xMom yMom zMom '

  WRITE(IF_PLOT,1015) iRegBin, nBins

  DO iBin = 1, nBins
    WRITE(IF_PLOT,1030,err=10) iBin,(yBin(iBin)-yHeightMin)/yHeight,sizeBin(iBin),&
                               massAlBin(iBin),massAlOxBin(iBin),massTotBin(iBin),&
			       diam2Bin(iBin),diam3Bin(iBin),diam4Bin(iBin),&
			       xMomBin(iBin),yMomBin(iBin),zMomBin(iBin)
  ENDDO ! iBin

! close file, handle errors ---------------------------------------------------

  CLOSE(IF_PLOT,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,fname )

1999 CONTINUE
  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_WRITE,__LINE__,fname )

! formats ---------------------------------------------------------------------

1005 FORMAT('TITLE="',A,'. Time: ',1PE11.5,'."')
1010 FORMAT('VARIABLES= ',A)
1015 FORMAT('ZONE T="',I5.5,'", I=',I10,', F=POINT')
1020 FORMAT(1X,I5,1X,1PE12.5,5(1X,I5),3(1X,1PE12.5))
1030 FORMAT(1X,I5,1X,1PE12.5,1X,I5,9(1X,1E23.16))

999  CONTINUE
END SUBROUTINE PLAG_BinSortSpatialDist

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_BinSortSpatialDist.F90,v $
! Revision 1.3  2008/12/06 08:45:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/13 22:00:00  fnajjar
! Initial import
!
!******************************************************************************







