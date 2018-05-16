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
! Purpose: write solution data to plot file in ASCII format.
!
! Description: none.
!
! Input: iReg    = region number
!        iLev    = current level
!        region  = region data (dimensions, plag variables)
!
! Output: to plot file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_WriteTecplotAscii.F90,v 1.3 2008/12/06 08:45:07 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_WriteTecplotAscii( iReg, iLev, region )

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
  INTEGER :: iReg, iLev

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: i

! ... local variables
  CHARACTER(CHRLEN)   :: RCSIdentString
  CHARACTER(CHRLEN+4) :: fname

  INTEGER :: errorFlag, nDimPlag
  INTEGER :: pidini,regini,icell,indexi,indexj,indexk
  INTEGER, POINTER :: cvPlagMass(:)
  INTEGER, POINTER :: aiv(:,:)

  REAL(RFREAL)          :: x,y,z,diam,evapor,temp,up,vp,wp
  REAL(RFREAL), POINTER :: cv(:,:), dv(:,:)
  REAL(RFREAL)          :: massp,muem,taup

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_WriteTecplotAscii.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PLAG_WriteTecplotAscii', 'PLAG_WriteTecplotAscii.F90' )

! set parameters --------------------------------------------------------------

  nDimPlag = region%levels(iLev)%plag%nPcls
  
!  PRINT*,'PLAG_writeTecplotAscii: iReg, nDimPlag = ',iReg,nDimPlag

! open file and write the header ----------------------------------------------

  IF (iReg == 1) THEN
    WRITE(fname,'(A,1PE11.5,A)') &
       TRIM(global%casename)//'.plag_',global%currentTime,'.dat'
    OPEN(IF_PLOT,FILE=fname,status='unknown',form='formatted',iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,fname )

    IF (global%currentTime <= 0._RFREAL) THEN
      WRITE(IF_PLOT,1005,err=10) TRIM(global%casename),global%timeStamp
    ELSE
      WRITE(IF_PLOT,1005,err=10) TRIM(global%casename),global%currentTime
    ENDIF
    
    WRITE(IF_PLOT,1010,err=10) 'x y z diam Evapor temp up vp wp tau pidini regini icellp ip jp kp'

  ENDIF   ! iReg=1

! for zero nDimPlag exit ------------------------------------------------------
 
  IF ( nDimPlag <= 0 ) GOTO 1999

! write zone header -----------------------------------------------------------

  WRITE(IF_PLOT,1015) iReg, nDimPlag

! write data ------------------------------------------------------------------
! pointer to variables

  aiv => region%levels(iLev)%plag%aiv
  cv  => region%levels(iLev)%plag%cv
  dv  => region%levels(iLev)%plag%dv
  cvPlagMass  => region%levels(iLev)%plag%cvPlagMass

  DO i = 1, nDimPlag
    x      = cv(CV_PLAG_XPOS,     i)
    y      = cv(CV_PLAG_YPOS,     i)
    z      = cv(CV_PLAG_ZPOS,     i)
    diam   = dv(DV_PLAG_DIAM,     i)
    evapor = cv(CV_PLAG_ENERVAPOR,i)
    temp   = dv(DV_PLAG_TEMP,     i)
    up     = dv(DV_PLAG_UVEL,     i)
    vp     = dv(DV_PLAG_VVEL,     i)
    wp     = dv(DV_PLAG_WVEL,     i)

    pidini = aiv(AIV_PLAG_PIDINI, i)
    regini = aiv(AIV_PLAG_REGINI, i)
    icell  = aiv(AIV_PLAG_ICELLS, i)
    indexi = aiv(AIV_PLAG_INDEXI, i)
    indexj = aiv(AIV_PLAG_INDEXJ, i)
    indexk = aiv(AIV_PLAG_INDEXK, i)

    massp = SUM( cv(cvPlagMass(:),i) )
    muem  = 3.6E-04_RFREAL
    taup  = massp/(3.0_RFREAL*global%pi*muem*diam)

    WRITE(IF_PLOT,1020,err=10) x,y,z,diam,evapor,temp,up,vp,wp,taup, &
       pidini,regini,icell,indexi,indexj,indexk

  ENDDO ! i 

! close file, handle errors ---------------------------------------------------

  IF (iReg == global%nRegions) THEN
    CLOSE(IF_PLOT,iostat=errorFlag)
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_FILE_CLOSE,__LINE__,fname )
  ENDIF

1999 CONTINUE
  CALL DeregisterFunction( global )
  GOTO 999

10   CONTINUE
  CALL ErrorStop( global,ERR_FILE_WRITE,__LINE__,fname )

! formats ---------------------------------------------------------------------

1005 FORMAT('TITLE="',A,'. Time: ',1PE11.5,'."')
1010 FORMAT('VARIABLES= ',A)
1015 FORMAT('ZONE T="',I5.5,'", I=',I10,', F=POINT')
1020 FORMAT(1P,10(1X,E13.6),1X,I10,5(1X,I8))

999  CONTINUE
END SUBROUTINE PLAG_WriteTecplotAscii

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_WriteTecplotAscii.F90,v $
! Revision 1.3  2008/12/06 08:45:07  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:18  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:00:47  fnajjar
! Initial revision after changing case
!
! Revision 1.7  2004/11/18 22:06:30  fnajjar
! Modified formatting to properly write high particle count index
!
! Revision 1.6  2004/11/13 22:03:55  fnajjar
! Included particle timescale in Tecplot output file
!
! Revision 1.5  2004/06/29 14:07:03  fnajjar
! Removed mixture properties from Tecplot file since mixture data not allocated and read
!
! Revision 1.4  2004/05/24 14:23:24  fnajjar
! Added variables to Tecplot file
!
! Revision 1.3  2004/03/05 22:52:00  fnajjar
! Added particle temperature to Tecplot output file and changed extension to .dat
!
! Revision 1.2  2004/03/02 21:51:15  jferry
! Added output of vapor energy to rplagpost output file
!
! Revision 1.1.1.1  2003/05/06 16:14:38  fnajjar
! Import of postprocessing tool for Rocpart
!
!******************************************************************************







