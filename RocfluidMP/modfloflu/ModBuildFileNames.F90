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
! Purpose: Collection of utility routines for building file names.
!
! Description: None
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ModBuildFileNames.F90,v 1.5 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE ModBuildFileNames

  USE ModParameters
  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModError

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  

! ==============================================================================
! Public data
! ==============================================================================

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: ModBuildFileNames.F90,v $ $Revision: 1.5 $'        

! ==============================================================================
! Public functions
! ==============================================================================

  PUBLIC :: BuildFileNameBasic, &
            BuildFileNamePlain, & 
            BuildFileNamePlainSteady, &
            BuildFileNamePlainUnsteady, &
            BuildFileNameSteady, &
            BuildFileNameUnsteady, & 
            BuildRegionIdString

! ==============================================================================
! Private functions
! ==============================================================================
                
! ******************************************************************************
! Routines
! ******************************************************************************
                
  CONTAINS




! ******************************************************************************
!
! Purpose: Build basic file name, that is, file name consisting only of 
!   directory, case name, extension, and region id.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   regId       Region index
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE BuildFileNameBasic(global,dest,ext,regId,fileName)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString,regIdString
  
! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest,regId
    TYPE(t_global), POINTER :: global 

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName  
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNameBasic',&
         'ModBuildFileNames.F90')

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN 
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE 
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,&
           __LINE__)
    END IF ! dest

    CALL BuildRegionIdString(global,regId,regIdString)

    WRITE(fileName,'(A)') TRIM(destString)//TRIM(global%caseName)// & 
                          TRIM(ext)//'_'//TRIM(regIdString)
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)     
  
  END SUBROUTINE BuildFileNameBasic








! ******************************************************************************
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
! ******************************************************************************

  SUBROUTINE BuildFileNamePlain(global,dest,ext,fileName)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString
  
! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest
    TYPE(t_global), POINTER :: global 

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName  
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNamePlain',&
         'ModBuildFileNames.F90')

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN 
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE 
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,&
           __LINE__)
    END IF ! dest

    WRITE(fileName,'(A)') TRIM(destString)//TRIM(global%caseName)//TRIM(ext)   
    
! ******************************************************************************
!   End
! ******************************************************************************
  
    CALL DeregisterFunction(global)     
  
  END SUBROUTINE BuildFileNamePlain








! ******************************************************************************
!
! Purpose: Build plain file name with stamp, that is, file name consisting 
!   only of directory, case name, iteration stamp, and extension.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   it          Iteration counter
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE BuildFileNamePlainSteady(global,dest,ext,it,fileName)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString
  
! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest,it
    TYPE(t_global), POINTER :: global 

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName  
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNamePlainSteady',&
         'ModBuildFileNames.F90')

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN 
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE 
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,&
           __LINE__)
    END IF ! dest

    WRITE(fileName,'(A,A,I6.6,A)') & 
      TRIM(destString)//TRIM(global%caseName),'_',it,TRIM(ext) 
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)     
  
  END SUBROUTINE BuildFileNamePlainSteady








! ******************************************************************************
!
! Purpose: Build plain file name with stamp, that is, file name consisting 
!   only of directory, case name, time stamp, and extension.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   tm          Time stamp
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE BuildFileNamePlainUnsteady(global,dest,ext,tm,fileName)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString
  
! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest
    REAL(RFREAL), INTENT(IN) :: tm
    TYPE(t_global), POINTER :: global 

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName  

! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNamePlainUnsteady',&
         'ModBuildFileNames.F90')

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN 
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE 
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,&
           __LINE__)
    END IF ! dest

    WRITE(fileName,'(A,A,1PE11.5,A)') & 
      TRIM(destString)//TRIM(global%caseName),'_',tm,TRIM(ext) 

! ******************************************************************************
! End
! ******************************************************************************

    CALL DeregisterFunction(global)     
  
  END SUBROUTINE BuildFileNamePlainUnsteady








! ******************************************************************************
!
! Purpose: Build file name for steady flow, that is, file name consisting of 
!   directory, case name, extension, region id, and iteration counter.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   regId       Region index
!   it          Iteration counter
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE BuildFileNameSteady(global,dest,ext,regId,it,fileName)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString,regIdString
  
! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest,regId,it
    TYPE(t_global), POINTER :: global 

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName  
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNameSteady',&
         'ModBuildFileNames.F90')

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN 
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE 
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,&
           __LINE__)
    END IF ! dest

    CALL BuildRegionIdString(global,regId,regIdString)

    WRITE(fileName,'(A,I6.6)') TRIM(destString)//TRIM(global%caseName)// & 
                               TRIM(ext)//'_'//TRIM(regIdString)//'_',it
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)     
  
  END SUBROUTINE BuildFileNameSteady







! ******************************************************************************
!
! Purpose: Build file name for unsteady flow, that is, file name consisting of 
!   directory, case name, extension, region id, and time stamp.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   dest        Integer parameter indicating destination directory
!   ext         Extension (including ".")
!   regId       Region index
!   tm          Time stamp
!
! Output: 
!   fileName    File name
!
! Notes: None.
!
! ******************************************************************************
 
  SUBROUTINE BuildFileNameUnsteady(global,dest,ext,regId,tm,fileName)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: destString,regIdString
  
! ==============================================================================
!   Arguments
! ==============================================================================

    CHARACTER(*), INTENT(IN) :: ext
    INTEGER, INTENT(IN) :: dest,regId
    REAL(RFREAL), INTENT(IN) :: tm
    TYPE(t_global), POINTER :: global 

    CHARACTER(CHRLEN), INTENT(OUT) :: fileName  
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildFileNameUnsteady',&
         'ModBuildFileNames.F90')

! ******************************************************************************
!   Build file name
! ******************************************************************************

    IF ( dest == FILEDEST_INDIR ) THEN 
      WRITE(destString,'(A)') global%inDir
    ELSE IF ( dest == FILEDEST_OUTDIR ) THEN
      WRITE(destString,'(A)') global%outDir
    ELSE 
      CALL ErrorStop(global,ERR_FILEDEST_INVALID,&
           __LINE__)
    END IF ! dest

    CALL BuildRegionIdString(global,regId,regIdString)

    WRITE(fileName,'(A,1PE11.5)') TRIM(destString)//TRIM(global%caseName)// & 
                                  TRIM(ext)//'_'//TRIM(regIdString)//'_',tm
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)     
  
  END SUBROUTINE BuildFileNameUnsteady








! ******************************************************************************
!
! Purpose: Build region id string.
!
! Description: None.
!
! Input: 
!   global      Global pointer
!   regId       Region index
!
! Output: 
!   regIdString Region string
!
! Notes: None.
!
! ******************************************************************************
 
  SUBROUTINE BuildRegionIdString(global,regId,regIdString)

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================
  
! ==============================================================================
!   Arguments
! ==============================================================================

    INTEGER, INTENT(IN) :: regId
    TYPE(t_global), POINTER :: global 

    CHARACTER(CHRLEN), INTENT(OUT) :: regIdString 
  
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'BuildRegionIdString',&
         'ModBuildFileNames.F90')

! ******************************************************************************
!   Write region id into region string
! ******************************************************************************

    WRITE(regIdString,'(I5.5)') regId
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)     
  
  END SUBROUTINE BuildRegionIdString







END MODULE ModBuildFileNames

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModBuildFileNames.F90,v $
! Revision 1.5  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/07 15:19:18  haselbac
! Removed tabs
!
! Revision 1.2  2004/10/19 19:28:38  haselbac
! Added BuildRegionIdString bcos needed in GENX modules, cosmetics
!
! Revision 1.1  2004/06/16 20:00:43  haselbac
! Initial revision
!
! ******************************************************************************













