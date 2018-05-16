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
! Purpose: Write solution file for turbulence in binary ROCFLU format.
!
! Description: None.
!
! Input: region = Pointer to region
!
! Output: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: TURB_rFLU_WriteSolutionBinary.F90,v 1.5 2008/12/06 08:44:45 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE TURB_RFLU_WriteSolutionBinary( region ) ! PUBLIC

  USE ModDataTypes
  USE ModGlobal, ONLY    : t_global
  USE ModDataStruct, ONLY: t_region 
  USE ModGrid, ONLY      : t_grid
  USE ModError
  USE ModMPI
  USE ModParameters
  
  USE ModBuildFileNames, ONLY: BuildFileNameSteady, & 
                               BuildFileNameUnsteady  
  
  USE TURB_ModParameters

  IMPLICIT NONE
  
! ... arguments
  TYPE(t_region), POINTER :: region  

! ... loop variables
  INTEGER :: i, j
   
! ... local variables
  CHARACTER(CHRLEN) :: iFileName,sectionString,RCSIdentString
  TYPE(t_grid),   POINTER :: grid
  TYPE(t_global), POINTER :: global

  INTEGER :: errorFlag,iFile,nVars
  REAL(RFREAL) :: currentTime
  REAL(RFREAL), DIMENSION(:,:), POINTER :: tv, tcv, vort
  REAL(RFREAL), DIMENSION(:),   POINTER :: lens
  
! ******************************************************************************
! Start, open file
! ******************************************************************************

  RCSIdentString = '$RCSfile: TURB_rFLU_WriteSolutionBinary.F90,v $ $Revision: 1.5 $'

  global => region%global

  CALL RegisterFunction(global,'TURB_RFLU_WriteSolutionBinary',&
  'TURB_rFLU_WriteSolutionBinary.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary turbulence file...'
  END IF ! global%verbLevel

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.turbb', & 
                               region%iRegionGlobal,global%currentTime, & 
                               iFileName)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN                                             
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       region%iRegionGlobal
      WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                          global%currentTime  
    END IF ! global%verbLevel
  ELSE
    CALL BuildFileNameSteady(global,FILEDEST_OUTDIR,'.turbb', & 
                             region%iRegionGlobal,global%currentIter,iFileName)  

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN                                      
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       region%iRegionGlobal  
      WRITE(STDOUT,'(A,3X,A,1X,I6.6)') SOLVER_NAME,'Current iteration '// &
                                       'number:',global%currentIter             
    END IF ! global%verbLevel
  ENDIF ! global%flowType

  iFile = IF_SOLUT
  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", & 
       IOSTAT=errorFlag)
  global%error = errorFlag          
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
  END IF ! global%error  

! ==============================================================================
! Header and general information
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
  END IF ! global%verbLevel

  sectionString = '# ROCTURB solution file'
  WRITE(iFile) sectionString
  
  sectionString = '# Precision and range'
  WRITE(iFile) sectionString
  WRITE(iFile) PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)

  sectionString = '# Initial residual'
  WRITE(iFile) sectionString
  WRITE(iFile) global%resInit
  
  sectionString = '# Physical time'
  WRITE(iFile) sectionString
  WRITE(iFile) global%currentTime 
  
  sectionString = '# Esg1Sum'
  WRITE(iFile) sectionString
  WRITE(iFile) global%esg1Sum
  
  sectionString = '# Esg4Sum'
  WRITE(iFile) sectionString
  WRITE(iFile) global%esg4Sum

! ==============================================================================
! Dimensions
! ==============================================================================
  
  nVars = region%turbInput%nOutField
  
  grid => region%grid  
  
  sectionString = '# Dimensions'
  WRITE(iFile) sectionString
  WRITE(iFile) grid%nCellsTot,nVars 
  
! ==============================================================================
! Eddy viscosity
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Eddy viscosity...'
  END IF ! global%verbLevel

  sectionString = '# Eddy viscosity'
  WRITE(iFile) sectionString

  IF (region%turbInput%modelClass == MODEL_LES) THEN
    tv => region%mixt%tv
    WRITE(iFile) (tv(TV_MIXT_MUET,j),j=1,grid%nCellsTot)
  ELSEIF (region%turbInput%modelClass == MODEL_RANS) THEN
    tcv => region%turb%cv
    WRITE(iFile) (tcv(CV_SA_NUTIL,j),j=1,grid%nCellsTot)
  ENDIF

! ==============================================================================
! Total vorticity
! ==============================================================================

  IF ( nVars > 1 ) THEN
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Total vorticity...'
    END IF ! global%verbLevel

    sectionString = '# Total vorticity'
    WRITE(iFile) sectionString
      
    vort => region%turb%vort
    WRITE(iFile) (SQRT( vort(XCOORD,j)**2 + &
                        vort(YCOORD,j)**2 + &
                        vort(ZCOORD,j)**2 ),j=1,grid%nCellsTot)
  ENDIF 
 
! ==============================================================================
! RaNS or DES length scale
! ==============================================================================

  IF ( nVars > 2 ) THEN
    IF (region%turbInput%modelClass == MODEL_RANS) THEN

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_LOW ) THEN  
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Model length scale...'
      END IF ! global%verbLevel

      sectionString = '# RANS length scale'
      WRITE(iFile) sectionString
      
      lens => region%turb%lens
      WRITE(iFile) (lens(j),j=1,grid%nCellsTot)
    ENDIF                                      
  ENDIF 
  
! ==============================================================================
! End marker
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
  END IF ! global%verbLevel

  sectionString = '# End'
  WRITE(iFile) sectionString  

! ==============================================================================
! Close file
! ==============================================================================

  CLOSE(iFile,IOSTAT=errorFlag)
  global%error = errorFlag      
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
  END IF ! global%error
   
  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN   
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary flow file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)
    
! ******************************************************************************
! End
! ******************************************************************************
 
END SUBROUTINE TURB_RFLU_WriteSolutionBinary


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_rFLU_WriteSolutionBinary.F90,v $
! Revision 1.5  2008/12/06 08:44:45  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:57  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/01/10 07:18:10  wasistho
! fixed screen print vort to lengthscale
!
! Revision 1.2  2004/06/16 20:01:34  haselbac
! Added use of ModBuildFileNames, cosmetics
!  
! Revision 1.1  2004/03/27 02:19:15  wasistho  
! added routines specific for Rocflu           
!
! ******************************************************************************







