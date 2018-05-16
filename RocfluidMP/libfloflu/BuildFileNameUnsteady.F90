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
! Purpose: Build file name for unsteady flow, that is, file name consisting of 
!   directory, case name, extension, domain id, and time stamp.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   id          Region index
!   tm          Time stamp
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
!******************************************************************************
!
! $Id: BuildFileNameUnsteady.F90,v 1.4 2008/12/06 08:44:08 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE BuildFileNameUnsteady(global,dest,ext,id,tm,fileName)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: destString,RCSIdentString
  
! ==============================================================================
! Arguments
! ==============================================================================

  CHARACTER(*), INTENT(IN) :: ext
  INTEGER, INTENT(IN) :: dest,id
  REAL(RFREAL), INTENT(IN) :: tm
  TYPE(t_global), POINTER :: global 
  
  CHARACTER(CHRLEN), INTENT(OUT) :: fileName  
  
! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: BuildFileNameUnsteady.F90,v $ $Revision: 1.4 $'

  CALL RegisterFunction(global,'BuildFileNameUnsteady',&
  'BuildFileNameUnsteady.F90')
   
  IF ( dest == FILEDEST_INDIR ) THEN 
    WRITE(destString,'(A)') global%inDir
  ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
    WRITE(destString,'(A)') global%outDir
  ELSE 
    CALL ErrorStop(global,ERR_FILEDEST_INVALID,__LINE__)
  END IF ! dest
           
  WRITE(fileName,'(A,A,I5.5,A,1PE11.5)') & 
    TRIM(destString)//TRIM(global%caseName)//TRIM(ext),'_',id,'_',tm
   
  CALL DeregisterFunction(global)     
    
! ******************************************************************************
! End
! ******************************************************************************
  
END SUBROUTINE BuildFileNameUnsteady


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: BuildFileNameUnsteady.F90,v $
!   Revision 1.4  2008/12/06 08:44:08  mtcampbe
!   Updated license.
!
!   Revision 1.3  2008/11/19 22:17:22  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.2  2006/04/07 15:19:15  haselbac
!   Removed tabs
!
!   Revision 1.1  2004/12/01 16:48:14  haselbac
!   Initial revision after changing case
!
!   Revision 1.1  2003/01/28 16:12:15  haselbac
!   Initial revision
!
! ******************************************************************************







