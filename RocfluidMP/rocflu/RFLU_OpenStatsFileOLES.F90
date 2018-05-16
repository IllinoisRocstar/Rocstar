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
! Purpose: Open file for optimal LES statistics.
!
! Description: None.
!
! Input: 
!   global      Pointer to global type
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_OpenStatsFileOLES.F90,v 1.8 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_OpenStatsFileOLES(global)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY : t_global
  USE ModMPI
  USE ModParameters
  
  USE ModBuildFileNames, ONLY: BuildFileNamePlain   
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: fname,RCSIdentString
  INTEGER :: errorFlag

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_OpenStatsFileOLES.F90,v $ $Revision: 1.8 $'

  CALL RegisterFunction(global,'RFLU_OpenStatsFileOLES',&
  'RFLU_OpenStatsFileOLES.F90' )

! ==============================================================================
! Open file
! ==============================================================================

  IF ( global%myProcid == MASTERPROC ) THEN   
    CALL BuildFileNamePlain(global,FILEDEST_OUTDIR,'.oles',fname) 

! ==============================================================================
!   Append to existing file (restart) or create new file
! ==============================================================================

    IF ( global%currentTime > 0.0_RFREAL ) THEN
      OPEN(IF_STATS_OLES,FILE=fname,FORM='FORMATTED',STATUS='OLD', &
                         POSITION='APPEND',IOSTAT=errorFlag)
    ELSE
      OPEN(IF_STATS_OLES,FILE=fname,FORM='FORMATTED',STATUS='UNKNOWN', &
                         IOSTAT=errorFlag)
    END IF ! global
    
    global%error = errorFlag
    IF ( global%error /= 0 ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname))
    END IF ! global%error
  END IF ! global%myProcid

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_OpenStatsFileOLES

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_OpenStatsFileOLES.F90,v $
! Revision 1.8  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2004/06/16 20:01:13  haselbac
! Added use of ModBuildFileNames, cosmetics
!
! Revision 1.5  2003/01/28 14:45:23  haselbac
! Use common building of file name
!
! Revision 1.4  2002/10/27 19:16:15  haselbac
! Removed tabs
!
! Revision 1.3  2002/10/08 15:49:29  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.2  2002/10/05 19:24:26  haselbac
! GENX integration: added outDir to file name
!
! Revision 1.1  2002/09/09 16:28:02  haselbac
! Initial revision
! 
! ******************************************************************************







