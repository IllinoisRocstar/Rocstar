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
! $Id: PLAG_BinSortNozzleInlet.F90,v 1.6 2008/12/06 08:45:07 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_BinSortNozzleInlet ( iReg, iLev, region, iRegBin )

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
  INTEGER :: iBinLog,nBinsLog
  INTEGER :: iContAl,iContAlOx
  INTEGER, POINTER :: cvMass(:)
  INTEGER, POINTER :: aiv(:,:)
  INTEGER, DIMENSION(0:500) :: sizeBin,sizeBinLog

  REAL(RFREAL)          :: diamBinMicron,x,y,z,diamMicron,diamMicronLog
  REAL(RFREAL)          :: diamMaxLog,diamMinLog,diamBinLogDelta,riBin
  REAL(RFREAL)          :: dLogjph,dLogjmh
  REAL(RFREAL)          :: massL,xMinRange,xPosPlag
  REAL(RFREAL)          :: diamMax,diamMin
  REAL(RFREAL), DIMENSION(0:500) :: compAlBin,compAlOxBin,compAlBinLog,&
                                    compAlOxBinLog
  REAL(RFREAL), DIMENSION(0:500) :: diamBin, diamBinLog,weightBin
  REAL(RFREAL), DIMENSION(0:500) :: diam32Bin,diam43Bin,&
                                    diam2Bin,diam3Bin,diam4Bin,massAlBin,&
				    massAlOxBin,massTotBin
  REAL(RFREAL), DIMENSION(0:500) :: xMomBin,yMomBin,zMomBin
  REAL(RFREAL), DIMENSION(0:500) :: uVelBin,vVelBin,wVelBin
  REAL(RFREAL), DIMENSION(0:500) :: diam32BinLog,diam43BinLog,&
                                    diam2BinLog,diam3BinLog,diam4BinLog,&
				    massAlBinLog,massAlOxBinLog,massTotBinLog
  REAL(RFREAL), DIMENSION(0:500) :: xMomBinLog,yMomBinLog,zMomBinLog
  REAL(RFREAL), DIMENSION(0:500) :: uVelBinLog,vVelBinLog,wVelBinLog
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:)

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_BinSortNozzleInlet.F90,v $ $Revision: 1.6 $'

  global => region%global

  CALL RegisterFunction( global,'PLAG_BinSortNozzleInlet', 'PLAG_BinSortNozzleInlet.F90' )

! set parameters --------------------------------------------------------------

  nDimPlag = region%levels(iLev)%plag%nPcls
 
  sizeBin(0:500) = 0
  sizeBinLog(0:500) = 0
  diamBin(0:500) = 0.0_RFREAL
  diamBinLog(0:500) = 0.0_RFREAL
  weightBin(0:500) = 1.0_RFREAL

  diamMinLog = LOG(  5.0_RFREAL)
  diamMaxLog = LOG(250.0_RFREAL)
  diamBinMicron = 5.0_RFREAL
  
  diamMin = 5.0_RFREAL
  diamMax = 250.0_RFREAL
  
  nBins=50
  nBinsLog=25
  
  iContAl   = 1
  iContAlOx = 2
  compAlBin(0:500)      = 0.0_RFREAL
  compAlOxBin(0:500)    = 0.0_RFREAL
  compAlBinLog(0:500)   = 0.0_RFREAL
  compAlOxBinLog(0:500) = 0.0_RFREAL

  massAlBin(0:500)   = 0.0_RFREAL
  massAlOxBin(0:500) = 0.0_RFREAL
  massTotBin(0:500)  = 0.0_RFREAL  

  xMomBin(0:500) = 0.0_RFREAL
  yMomBin(0:500) = 0.0_RFREAL
  zMomBin(0:500) = 0.0_RFREAL
  
  uVelBin(0:500) = 0.0_RFREAL
  vVelBin(0:500) = 0.0_RFREAL
  wVelBin(0:500) = 0.0_RFREAL

  massAlBinLog(0:500)   = 0.0_RFREAL
  massAlOxBinLog(0:500) = 0.0_RFREAL
  massTotBinLog(0:500)  = 0.0_RFREAL  

  xMomBinLog(0:500) = 0.0_RFREAL
  yMomBinLog(0:500) = 0.0_RFREAL
  zMomBinLog(0:500) = 0.0_RFREAL
  
  uVelBinLog(0:500) = 0.0_RFREAL
  vVelBinLog(0:500) = 0.0_RFREAL
  wVelBinLog(0:500) = 0.0_RFREAL

  IF (iReg /= iRegBin) GOTO 1999
  
  PRINT*,'PLAG_BinSortNozzleInlet: iRegBin, nDimPlag = ',iRegBin,nDimPlag

! set values ------------------------------------------------------------------

! TEMPORARY: Data Pertinent to c1

  xMinRange = 0.0_RFREAL

  SELECT CASE (iRegBin)
    CASE(41)
      xMinRange = 0.39_RFREAL

    CASE(52)
      xMinRange = 0.465_RFREAL

    CASE DEFAULT
      xMinRange = 0.0_RFREAL
  END SELECT ! iRegBin

  diamBinLogDelta = (diamMaxLog-diamMinLog)/REAL((nBinsLog-1),KIND=RFREAL)

  DO  iBin = 1, nBinsLog
    riBin = REAL(iBin,KIND=RFREAL)
    diamBinLog(iBin) = ((nBinsLog-riBin)*diamMinLog +(riBin-1.0_RFREAL)*diamMaxLog)/&
                        REAL((nBinsLog-1),KIND=RFREAL)

    dLogjph = ( (nBinsLog-riBin-0.5_RFREAL)*diamMinLog  &
               +(riBin+0.5_RFREAL-1)*diamMaxLog)/       &
                REAL((nBinsLog-1),KIND=RFREAL)

    dLogjmh = ( (nBinsLog-riBin+0.5_RFREAL)*diamMinLog &
               +(riBin-0.5_RFREAL-1)*diamMaxLog)/      &
                 REAL((nBinsLog-1),KIND=RFREAL)

    weightBin(iBin) = EXP(dLogjph) -EXP(dLogjmh)

!   PRINT*,' iBin, diamBin, weightBin = ', iBin,EXP(diamBinLog(iBin)),weightBin(iBin)
  END DO

! set pointers ----------------------------------------------------------------

  aiv => region%levels(iLev)%plag%aiv
  cv  => region%levels(iLev)%plag%cv
  dv  => region%levels(iLev)%plag%dv
  cvMass  => region%levels(iLev)%plag%cvPlagMass

! determine size ---------------------------------------------------------------

  DO i = 1, nDimPlag
    xPosPlag = cv(CV_PLAG_XPOS,i)
    
    IF ( xPosPlag >= xMinRange ) THEN
      diamMicron = dv(DV_PLAG_DIAM,i)*1.0E+06_RFREAL

      iBin=NINT(diamMicron/diamBinMicron)

      IF(diamMicron < diamMin ) iBin = 1
      IF(diamMicron > diamMax ) iBin = nBins

      sizeBin(iBin) = sizeBin(iBin)+1

      massL = SUM(cv(cvMass(:),i))   

      massAlBin(iBin)   = massAlBin(iBin)   +cv(cvMass(iContAl),i)
      massAlOxBin(iBin) = massAlOxBin(iBin) +cv(cvMass(iContAlOx),i)
      massTotBin(iBin)  = massTotBin(iBin)  +massL

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
   TRIM(global%casename)//'.plag_bin_reg_',iRegBin,'_',global%currentTime,'.dat'
  OPEN(IF_PLOT,FILE=fname,status='unknown',form='formatted',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,fname )

  IF (global%currentTime <= 0._RFREAL) THEN
    WRITE(IF_PLOT,1005,err=10) TRIM(global%casename),global%timeStamp
  ELSE
    WRITE(IF_PLOT,1005,err=10) TRIM(global%casename),global%currentTime
  ENDIF

  WRITE(IF_PLOT,1010,err=10) 'diamBin sizeBin compAl compAlOx diam32 diam43 u v w'

  WRITE(IF_PLOT,1015) iRegBin, nBins

  DO iBin = 1, nBins
    WRITE(IF_PLOT,1020,err=10) iBin*5,sizeBin(iBin),&
                               NINT(compAlBin(iBin)*100.0_RFREAL),&
			       NINT(compAlOxBin(iBin)*100.0_RFREAL), &
			       NINT(diam32Bin(iBin)),NINT(diam43Bin(iBin)),&
			       uVelBin(iBin),vVelBin(iBin),wVelBin(iBin)
  ENDDO ! iBin

  CLOSE(IF_PLOT,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,fname )

! write conserved variables for time-averaging ----------------------------------

  WRITE(fname,'(A,I2,A,1PE11.5,A)') &
   TRIM(global%casename)//'.plag_cv_bin_reg_',iRegBin,'_',global%currentTime,'.dat'
  OPEN(IF_PLOT,FILE=fname,status='unknown',form='formatted',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,fname )

  IF (global%currentTime <= 0._RFREAL) THEN
    WRITE(IF_PLOT,1005,err=10) TRIM(global%casename),global%timeStamp
  ELSE
    WRITE(IF_PLOT,1005,err=10) TRIM(global%casename),global%currentTime
  ENDIF

  WRITE(IF_PLOT,1010,err=10) 'diamBin sizeBin massAl massAlOxBin massTot diam2 diam3 diam4 xMom yMom zMom'

  WRITE(IF_PLOT,1015) iRegBin, nBins

  DO iBin = 1, nBins
    WRITE(IF_PLOT,1030,err=10) iBin*5,sizeBin(iBin),&
                               massAlBin(iBin),massAlOxBin(iBin),massTotBin(iBin),&
			       diam2Bin(iBin),diam3Bin(iBin),diam4Bin(iBin),&
			       xMomBin(iBin),yMomBin(iBin),zMomBin(iBin)
  ENDDO ! iBin

  CLOSE(IF_PLOT,iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,fname )

! compute log-binned distribution -----------------------------------------------

  DO i = 1, nDimPlag
    xPosPlag = cv(CV_PLAG_XPOS,i)
    
    IF ( xPosPlag >= xMinRange ) THEN
      diamMicron = dv(DV_PLAG_DIAM,i)*1.0E+06_RFREAL
      diamMicronLog = LOG(diamMicron) 

      iBinLog=NINT(1.0_RFREAL + (diamMicronLog-diamMinLog)/diamBinLogDelta)

      IF(diamMicronLog < diamMinLog ) iBinLog = 1
      IF(diamMicronLog > diamMaxLog ) iBinLog = nBinsLog

      sizeBinLog(iBinLog) = sizeBinLog(iBinLog)+1

      massL = SUM(cv(cvMass(:),i))    

      massAlBinLog(iBinLog)   = massAlBinLog(iBinLog)   +cv(cvMass(iContAl),i)
      massAlOxBinLog(iBinLog) = massAlOxBinLog(iBinLog) +cv(cvMass(iContAlOx),i)
      massTotBinLog(iBinLog)  = massTotBinLog(iBinLog)  +massL

      diam2BinLog(iBinLog) = diam2BinLog(iBinLog) +diamMicron**2
      diam3BinLog(iBinLog) = diam3BinLog(iBinLog) +diamMicron**3
      diam4BinLog(iBinLog) = diam4BinLog(iBinLog) +diamMicron**4

      xMomBinLog(iBinLog) = xMomBinLog(iBinLog) +cv(CV_PLAG_XMOM,i)
      yMomBinLog(iBinLog) = yMomBinLog(iBinLog) +cv(CV_PLAG_YMOM,i)
      zMomBinLog(iBinLog) = zMomBinLog(iBinLog) +cv(CV_PLAG_ZMOM,i)
    ENDIF ! xPosPlag
  ENDDO ! i

  DO iBinLog = 1, nBinsLog
    IF ( sizeBinLog(iBinLog) > 0 ) THEN
      compAlBinLog(iBinLog)   = massAlBinLog(iBinLog)/massTotBinLog(iBinLog)
      compAlOxBinLog(iBinLog) = massAlOxBinLog(iBinLog)/massTotBinLog(iBinLog)

      diam43BinLog(iBinLog) = diam4BinLog(iBinLog)/diam3BinLog(iBinLog)
      diam32BinLog(iBinLog) = diam3BinLog(iBinLog)/diam2BinLog(iBinLog)

      uVelBinLog(iBinLog) = xMomBinLog(iBinLog)/massTotBinLog(iBinLog)
      vVelBinLog(iBinLog) = yMomBinLog(iBinLog)/massTotBinLog(iBinLog)
      wVelBinLog(iBinLog) = zMomBinLog(iBinLog)/massTotBinLog(iBinLog)
    END IF ! sizeBinLog
  ENDDO ! iBinLog

! open file and write the header ----------------------------------------------

  WRITE(fname,'(A,I2,A,1PE11.5,A)') &
   TRIM(global%casename)//'.plag_logbin_reg_',iRegBin,'_',global%currentTime,'.dat'
  OPEN(IF_PLOT,FILE=fname,status='unknown',form='formatted',iostat=errorFlag)
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,fname )

  IF (global%currentTime <= 0._RFREAL) THEN
    WRITE(IF_PLOT,1005,err=10) TRIM(global%casename),global%timeStamp
  ELSE
    WRITE(IF_PLOT,1005,err=10) TRIM(global%casename),global%currentTime
  ENDIF

  WRITE(IF_PLOT,1010,err=10) 'diamBinLog sizeBinLog sizeBinWeight compAl compAlOx diam32 diam43 u v w'
  WRITE(IF_PLOT,1015) iRegBin, nBinsLog

  DO iBinLog = 1, nBinsLog
    WRITE(IF_PLOT,1025,err=10) EXP(diamBinLog(iBinLog)),sizeBinLog(iBinLog),&
                               sizeBinLog(iBinLog)/weightBin(iBinLog),&
                               NINT(compAlBinLog(iBinLog)*100.0_RFREAL),&
                               NINT(compAlOxBinLog(iBinLog)*100.0_RFREAL), &
                               NINT(diam32BinLog(iBinLog)),NINT(diam43BinLog(iBinLog)),&
                               uVelBinLog(iBinLog),vVelBinLog(iBinLog),wVelBinLog(iBinLog)
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
1020 FORMAT(1X,I5,1X,I6,4(1X,I5),3(1X,1PE12.5))
1025 FORMAT((1X,1PE12.5),1X,I5,(1X,1PE12.5),4(1X,I5),3(1X,1PE12.5))
1030 FORMAT(1X,I4,1X,I6,9(1X,E23.16))

999  CONTINUE
END SUBROUTINE PLAG_BinSortNozzleInlet

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_BinSortNozzleInlet.F90,v $
! Revision 1.6  2008/12/06 08:45:07  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2004/12/13 18:07:15  fnajjar
! Bug fix for out-of-bound bin indices
!
! Revision 1.3  2004/11/17 22:14:15  fnajjar
! Cosmetic changes and bug fixes
!
! Revision 1.2  2004/11/13 22:01:04  fnajjar
! Included log binning
!
! Revision 1.1  2004/05/24 14:25:19  fnajjar
! Initial import of binning routine
!
!******************************************************************************







