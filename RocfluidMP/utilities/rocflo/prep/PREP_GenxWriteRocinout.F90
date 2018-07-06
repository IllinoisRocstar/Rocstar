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
! Purpose: mkdir Rocin/out if not exist and write Rocin control files.
!
! Description: regions mapping info and initial solution file names are 
!              written into mapping files
!
! Input: global = global data, to support MP prep.
!
! Output: mapping files.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PREP_GenxWriteRocinout.F90,v 1.8 2008/12/06 08:44:50 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE GenxWriteRocinout( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  IMPLICIT NONE
  INCLUDE "comf90.h"

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  CHARACTER(2*CHRLEN+17) :: mapfname, solfname
  LOGICAL :: dirExist
  INTEGER :: errFlg, testFlg

!******************************************************************************

  CALL RegisterFunction( global,'GenxWriteRocinout',&
  'PREP_GenxWriteRocinout.F90' )

! open files and write mapping info ------------------------------------------

! mixture volume mapping file:

  dirExist = .FALSE.
  testFlg  = COM_call_system( "cd ../Rocout")  ! Intel-Fortran INQUIRE problem
  INQUIRE (FILE=TRIM('../Rocout'), EXIST=dirExist)
  IF ((.NOT. dirExist) .AND. (testFlg /= 0)) THEN
    errFlg = COM_call_system( "mkdir ../Rocout")
    global%error = errFlg
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_SYSTEM_COMMAND,__LINE__,'mkdir ../Rocout' )
  ENDIF

  dirExist = .FALSE.
  testFlg  = COM_call_system( "cd ../Rocin")  ! Intel-Fortran INQUIRE problem
  INQUIRE (FILE=TRIM('../Rocin'), EXIST=dirExist)
  IF ((.NOT. dirExist) .AND. (testFlg /= 0)) THEN
    errFlg = COM_call_system( "mkdir ../Rocin")
    global%error = errFlg
    IF (global%error /= 0) &
      CALL ErrorStop( global,ERR_SYSTEM_COMMAND,__LINE__,'mkdir ../Rocin' )
  ENDIF

  WRITE(mapfname,1000)'../Rocin/fluid_in_00.000000.txt'
#ifndef USE_CGNS
  WRITE(solfname,1000)'fluid_%5b.hdf'
#else
  WRITE(solfname,1000)'fluid_%5b.cgns'
#endif
  OPEN(IF_REGMAP,file=mapfname,form='formatted',status='unknown',iostat=errFlg)
  global%error = errFlg
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(mapfname) )

  WRITE(IF_REGMAP,1000) '@Proc: *'
  WRITE(IF_REGMAP,1010) '@Files: '//TRIM(solfname)
  WRITE(IF_REGMAP,1000) '@Panes: @Block ',global%nRegions,' 100 100'
  CLOSE(IF_REGMAP)

! mixture surface mapping file:

  WRITE(mapfname,1000)'../Rocin/ifluid_in_00.000000.txt'
#ifndef USE_CGNS
  WRITE(solfname,1000)'ifluid_%5b.hdf'
#else
  WRITE(solfname,1000)'ifluid_%5b.cgns'
#endif
  OPEN(IF_REGMAP,file=mapfname,form='formatted',status='unknown',iostat=errFlg)
  global%error = errFlg
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(mapfname) )

  WRITE(IF_REGMAP,1000) '@Proc: *'
  WRITE(IF_REGMAP,1010) '@Files: '//TRIM(solfname)
  WRITE(IF_REGMAP,1000) '@Panes: @Block ',global%nRegions,' 100 100'
  CLOSE(IF_REGMAP)

! MP mapping files:

  WRITE(mapfname,1000)'../Rocin/fluid_plag_in_00.000000.txt'
  WRITE(solfname,1000)'fluid_plag_%5b.hdf'
  OPEN(IF_REGMAP,file=mapfname,form='formatted',status='unknown',iostat=errFlg)
  global%error = errFlg
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(mapfname) )

  WRITE(IF_REGMAP,1000) '@Proc: *'
  WRITE(IF_REGMAP,1010) '@Files:'
  WRITE(IF_REGMAP,1000) '@Panes: @Block ',global%nRegions,' 100 100'
  CLOSE(IF_REGMAP)

! notes

  WRITE(mapfname,1000)'../RocfloNotes.txt'
  OPEN(IF_REGMAP,file=mapfname,form='formatted',status='unknown',iostat=errFlg)
  global%error = errFlg
  IF (global%error /= 0) &
    CALL ErrorStop( global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(mapfname) )

  WRITE(IF_REGMAP,1000) 'Dataset generated using Rfloprep pre version 2.3.0.0'
  WRITE(IF_REGMAP,1000) 'must be run with Rocstar compiled with' 
  WRITE(IF_REGMAP,1000) 'PRE_RFLOPREP_V2300=1 flag'
  CLOSE(IF_REGMAP)
    
1000 FORMAT(A,I10,A)
1010 FORMAT(A,A)

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE GenxWriteRocinout

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PREP_GenxWriteRocinout.F90,v $
! Revision 1.8  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2005/05/11 19:56:32  wasistho
! modified RocfloNotes.txt
!
! Revision 1.5  2005/04/20 02:58:16  wasistho
! modified text in RocfloNotes.txt
!
! Revision 1.4  2005/04/20 02:51:17  wasistho
! added text in RocfloNotes.txt
!
! Revision 1.3  2005/04/19 18:41:47  wasistho
! added RocfloNotes.txt
!
! Revision 1.2  2005/01/11 01:35:17  wasistho
! added testFlag due to Intel-Fortran problem with INQUIRE to directory
!
! Revision 1.1  2004/12/03 02:20:08  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:40:49  wasistho
! lower to upper case
!
! Revision 1.5  2004/10/13 16:20:34  jiao
! Updated to write the Rocin control files using the %5b placeholder.
!
! Revision 1.4  2004/10/12 04:32:36  wasistho
! split to one block per HDF file
!
! Revision 1.3  2004/10/09 19:41:26  jiao
! Changed the definition of offset for BlockCyclic and BlockBlockCyclic mapping.
!
! Revision 1.2  2004/10/07 04:25:31  jiao
! Fixed Rocin control files.
!
! Revision 1.1  2004/07/27 03:39:05  wasistho
! initial import genxWriteRocinout and printPrepInput
!
!
!******************************************************************************







