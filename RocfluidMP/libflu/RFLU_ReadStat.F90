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
! Purpose: read in time averaged solution
!
! Description: file to be read in contains statistics solution in 
!              interior cells; only binary format is supported
!
! Input: region      = region data (dimensions)
!        integrTime  = integrated time from file
!        mixttav     = time averaged mixture data from file
!        turbtav     = time averaged TURB data from file
!
! Output: region%mixt%tav     = time averaged mixture variables
!         region%turb%tav     = time averaged TURB variables
!         global%integrTime   = integrated time
!
! Notes: each region has statistics solution file
!
!******************************************************************************
!
! $Id: RFLU_ReadStat.F90,v 1.11 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_ReadStat( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: ijk, ind, l

! ... local variables
  CHARACTER(CHRLEN+23) :: fname
  CHARACTER(CHRLEN)    :: RCSIdentString,msg,timeString1,timeString2

  INTEGER :: errorFlag,iReg, ijkbeg, ijkend, nCells
  INTEGER :: nTavgVar
  INTEGER :: mixtVarId(2,region%global%mixtNStat)
  INTEGER :: turbVarId(2,region%global%turbNStat)

  REAL(RFREAL)          :: currentTime
  REAL(RFREAL), POINTER :: mixttav(:,:), turbtav(:,:) 

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ReadStat.F90,v $'

  global => region%global

  CALL RegisterFunction(global,'RFLU_ReadStat',&
  'RFLU_ReadStat.F90')

! open file -----------------------------------------------------------------
  
  iReg = region%iRegionGlobal
  WRITE(fname,'(A,I5.5,A,1PE11.5)') TRIM(global%outDir)// & 
                                    TRIM(global%casename)//'.statb_',iReg, &
                                    '_',global%currentTime

  OPEN(IF_STAT,file=fname,form='unformatted',status='old',iostat=errorFlag)
  global%error = errorFlag   
  IF (global%error /= 0) &
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname))

! read physical time and integrated time ------------------------------------

  READ(IF_STAT,err=10,end=10) currentTime
  READ(IF_STAT,err=10,end=10) global%integrTime
  WRITE(timeString1,'(1PE11.5)') global%currentTime
  WRITE(timeString2,'(1PE11.5)') currentTime          

  IF ( TRIM(timeString1) /= TRIM(timeString2) ) THEN
    WRITE(msg,1010) iReg,currentTime,global%currentTime
    CALL ErrorStop(global,ERR_TIME_SOLUTION,__LINE__, & 
                   msg//' File: '//TRIM(fname) )
  ENDIF

! read dimensions & check them ----------------------------------------------

  READ(IF_STAT,err=10,end=10) nCells

  IF (nCells /= region%grid%nCells) THEN
    WRITE(msg,1000) iReg, nCells, region%grid%nCells
    CALL ErrorStop(global,ERR_GRID_DIMENSIONS,__LINE__,msg)
  ENDIF

! mixture NSTAT and ID

  IF (global%mixtNStat > 0) THEN
    READ(IF_STAT,err=10,end=10) nTavgVar,mixtVarId(1,:)
    mixtVarId(2,:) = mod(mixtVarId(1,:),10)
    mixtVarId(1,:) = (mixtVarId(1,:)-mixtVarId(2,:))/10

    IF (nTavgVar /= global%mixtNStat) THEN
      CALL ErrorStop(global,ERR_STATS_RESTART,__LINE__)
    ENDIF

    DO ind=1,2
      DO l=1,global%mixtNStat
        IF (mixtVarId(ind,l)  /= global%mixtStatId(ind,l)) THEN
          CALL ErrorStop(global,ERR_STATS_RESTART,__LINE__)
        ENDIF
      ENDDO
    ENDDO
  ENDIF

! turbulence NSTAT and ID

#ifdef TURB
  IF (global%turbNStat > 0) THEN
    READ(IF_STAT,err=10,end=10) nTavgVar,turbVarId(1,:)
    turbVarId(2,:) = mod(turbVarId(1,:),10)
    turbVarId(1,:) = (turbVarId(1,:)-turbVarId(2,:))/10

    IF (nTavgVar /= global%turbNStat) THEN
      CALL ErrorStop(global,ERR_STATS_RESTART,__LINE__)
    ENDIF

    DO ind=1,2
      DO l=1,global%turbNStat
        IF (turbVarId(ind,l)  /= global%turbStatId(ind,l)) THEN
          CALL ErrorStop(global,ERR_STATS_RESTART,__LINE__)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
#endif

! read time averaged variables ----------------------------------------------

  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME,'  - read statistics'

  ijkbeg = 1
  ijkend = region%grid%nCells

  IF (global%mixtNStat > 0) THEN
    mixttav => region%mixt%tav
    READ(IF_STAT,err=10,end=10) ((mixttav(l,ijk), ijk=ijkbeg,ijkend), &
                                  l=1,global%mixtNStat)
  ENDIF
#ifdef TURB
  IF (global%turbNStat > 0) THEN
    turbtav => region%turb%tav
    READ(IF_STAT,err=10,end=10) ((turbtav(l,ijk), ijk=ijkbeg,ijkend), &
                                  l=1,global%turbNStat)
  ENDIF
#endif


! close file ----------------------------------------------------------------

  CLOSE(IF_STAT,iostat=errorFlag)
  global%error = errorFlag   
  IF (global%error /= 0) &
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname))

! finalization & error handling ---------------------------------------------

  CALL DeregisterFunction(global)
  GOTO 999

10   CONTINUE
  CALL ErrorStop(global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname))

999  CONTINUE
1000 FORMAT('Region ',I5,', nCells= ',I10,', nCellExpected= ',I10)
1010 FORMAT('Region ',I5,', time in file is= ',1PE12.5, &
            ' but it should be= ',E12.5,'.')

END SUBROUTINE RFLU_ReadStat

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ReadStat.F90,v $
! Revision 1.11  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2006/01/02 11:22:09  wasistho
! changed timeStamp to currentTime
!
! Revision 1.8  2002/11/02 01:51:40  wasistho
! Added TURB statistics
!
! Revision 1.7  2002/10/08 15:48:57  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.6  2002/10/05 18:52:34  haselbac
! Cosmetic changes only
!
! Revision 1.5  2002/09/09 14:15:01  haselbac
! global now under regions
!
! Revision 1.4  2002/07/22 15:45:50  wasistho
! Cleaned-up conforming Coding Rule
!
! Revision 1.3  2002/06/18 00:37:07  wasistho
! Added prefix SOLVER NAME to satistics STDOutput
!
! Revision 1.2  2002/06/17 18:33:34  wasistho
! modified times matching check
!
! Revision 1.1  2002/06/14 21:23:56  wasistho
! Added time avg statistics
!
!
!******************************************************************************







