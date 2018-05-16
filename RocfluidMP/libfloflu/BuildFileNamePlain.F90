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
! Purpose: Build plain file name, that is, file name consisting only of 
!   directory, case name, and extension.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
!******************************************************************************
!
! $Id: BuildFileNamePlain.F90,v 1.4 2008/12/06 08:44:08 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE BuildFileNamePlain(global,dest,ext,fileName)

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
  INTEGER, INTENT(IN) :: dest
  TYPE(t_global), POINTER :: global 
  
  CHARACTER(CHRLEN), INTENT(OUT) :: fileName  
  
! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: BuildFileNamePlain.F90,v $ $Revision: 1.4 $'

  CALL RegisterFunction(global,'BuildFileNamePlain',&
  'BuildFileNamePlain.F90')
   
  IF ( dest == FILEDEST_INDIR ) THEN 
    WRITE(destString,'(A)') global%inDir
  ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
    WRITE(destString,'(A)') global%outDir
  ELSE 
    CALL ErrorStop(global,ERR_FILEDEST_INVALID,__LINE__)
  END IF ! dest
   
  WRITE(fileName,'(A)') TRIM(destString)//TRIM(global%caseName)//TRIM(ext) 
   
  CALL DeregisterFunction(global)     
    
! ******************************************************************************
! End
! ******************************************************************************
  
END SUBROUTINE BuildFileNamePlain


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: BuildFileNamePlain.F90,v $
!   Revision 1.4  2008/12/06 08:44:08  mtcampbe
!   Updated license.
!
!   Revision 1.3  2008/11/19 22:17:22  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.2  2006/04/07 15:19:15  haselbac
!   Removed tabs
!
!   Revision 1.1  2004/12/01 16:48:09  haselbac
!   Initial revision after changing case
!
!   Revision 1.1  2003/01/28 16:12:15  haselbac
!   Initial revision
!
! ******************************************************************************







