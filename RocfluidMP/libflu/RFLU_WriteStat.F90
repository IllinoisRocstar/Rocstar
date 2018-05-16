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
! Purpose: store time averaged solution
!
! Description: solution in interior cells is stored;
!              only binary format is supported
!
! Input: region = region data (dimensions)
!
! Output: global%currentTime           = physical time
!         global%integrTime            = integrated time
!         region%grid%nCells           = region number of cells
!         global%mixtNStat             = number of mixture statistics variables
!         global%mixtStatId            = mixturestatistics variable Id
!         region%mixt%tav              = time averaged mixture variables
!         global%turbNStat             = number of TURB statistics variables
!         global%turbStatId            = TURB statistics variable Id
!         region%turb%tav              = time averaged TURB variables
!         output to file
!
! Notes: each region has statistics solution file
!
!******************************************************************************
!
! $Id: RFLU_WriteStat.F90,v 1.9 2008/12/06 08:44:13 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_WriteStat( region )

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
  INTEGER :: ijk, l

! ... local variables
  CHARACTER(CHRLEN+23) :: fname
  CHARACTER(CHRLEN)    :: RCSIdentString

  INTEGER :: errorFlag,iReg, ijkbeg, ijkend
  INTEGER :: mixtStatId(region%global%mixtNStat)
  INTEGER :: turbStatId(region%global%turbNStat)

  REAL(RFREAL), POINTER :: mixttav(:,:), turbtav(:,:) 

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_WriteStat.F90,v $'

  global => region%global

  CALL RegisterFunction(global,'RFLU_WriteStat',&
  'RFLU_WriteStat.F90')

! open file ----------------------------------------------------------------

  iReg = region%iRegionGlobal
  WRITE(fname,'(A,I5.5,A,1PE11.5)') TRIM(global%outDir)// & 
                                    TRIM(global%casename)//'.statb_',iReg, &
                                    '_',global%currentTime

  OPEN(IF_STAT,file=fname,form='unformatted',status='unknown', &
                                             iostat=errorFlag)
  global%error = errorFlag                                              
  IF (global%error /= 0) &
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname))

! write physical time and integrated time ----------------------------------

  WRITE(IF_STAT,err=10) global%currentTime
  WRITE(IF_STAT,err=10) global%integrTime

! write dimensions ---------------------------------------------------------

  WRITE(IF_STAT,err=10) region%grid%nCells

! mixture NSTAT and ID

  IF (global%mixtNStat > 0) THEN
    mixtStatId(:)=global%mixtStatId(1,:)*10+global%mixtStatId(2,:)
    WRITE(IF_STAT,err=10) global%mixtNStat,mixtStatId
  ENDIF

! turbulence NSTAT and ID

#ifdef TURB
  IF (global%turbNStat > 0) THEN
    turbStatId(:)=global%turbStatId(1,:)*10+global%turbStatId(2,:)
    WRITE(IF_STAT,err=10) global%turbNStat,turbStatId
  ENDIF
#endif

! write time averaged variables --------------------------------------------

  IF (global%myProcid==MASTERPROC .AND. global%verbLevel>=VERBOSE_HIGH) &
    WRITE(STDOUT,'(A)') SOLVER_NAME,'  - write statistics'

  ijkbeg = 1
  ijkend = region%grid%nCells

  IF (global%mixtNStat > 0) THEN
    mixttav => region%mixt%tav
    WRITE(IF_STAT,err=10) ((mixttav(l,ijk), ijk=ijkbeg,ijkend), &
                              l=1,global%mixtNStat)
  ENDIF
#ifdef TURB
  IF (global%turbNStat > 0) THEN
    turbtav => region%turb%tav
    WRITE(IF_STAT,err=10) ((turbtav(l,ijk), ijk=ijkbeg,ijkend), &
                              l=1,global%turbNStat)
  ENDIF
#endif

! close file ---------------------------------------------------------------

  CLOSE(IF_STAT,iostat=errorFlag)
  global%error = errorFlag   
  IF (global%error /= 0) &
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname))

! finalization & error handling --------------------------------------------

  CALL DeregisterFunction(global)
  GOTO 999

10   CONTINUE
  CALL ErrorStop(global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname))

999  CONTINUE

END SUBROUTINE RFLU_WriteStat

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_WriteStat.F90,v $
! Revision 1.9  2008/12/06 08:44:13  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2002/11/02 01:51:49  wasistho
! Added TURB statistics
!
! Revision 1.6  2002/10/08 15:48:57  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.5  2002/10/05 18:52:50  haselbac
! Cosmetic changes only
!
! Revision 1.4  2002/09/09 14:15:02  haselbac
! global now under regions
!
! Revision 1.3  2002/07/22 15:46:09  wasistho
! Cleaned-up conforming Coding Rule
!
! Revision 1.2  2002/06/18 00:37:16  wasistho
! Added prefix SOLVER NAME to satistics STDOutput
!
! Revision 1.1  2002/06/14 21:23:56  wasistho
! Added time avg statistics
!
!
!******************************************************************************











