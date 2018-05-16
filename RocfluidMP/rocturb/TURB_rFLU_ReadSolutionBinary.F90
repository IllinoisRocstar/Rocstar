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
! Purpose: Read file for turbulence quantities in binary ROCFLU format.
!
! Description: None.
!
! Input: pRegion = Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: TURB_rFLU_ReadSolutionBinary.F90,v 1.6 2009/08/26 12:28:53 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE TURB_RFLU_ReadSolutionBinary( region ) ! PUBLIC

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
  INTEGER :: iVars, j
   
! ... local variables
  CHARACTER(CHRLEN) :: errorString,iFileName,sectionString,RCSIdentString, & 
                       timeString1,timeString2
  TYPE(t_grid),   POINTER :: grid
  TYPE(t_global), POINTER :: global

  INTEGER :: errorFlag,iFile,loopCounter,nCellsTot,nCellsExpected,nVars, &
             nVarsExpected,precActual,precExpected,rangeActual,rangeExpected
  REAL(RFREAL) :: currentTime
  REAL(RFREAL), DIMENSION(:,:), POINTER :: tv, tcv, vort
  REAL(RFREAL), DIMENSION(:),   POINTER :: lens
  
! ******************************************************************************
! Start, open file
! ******************************************************************************

  RCSIdentString = '$RCSfile: TURB_rFLU_ReadSolutionBinary.F90,v $ $Revision: 1.6 $'

  global => region%global
  CALL RegisterFunction(global,'TURB_RFLU_ReadSolutionBinary',&
  'TURB_rFLU_ReadSolutionBinary.F90')

#ifdef GENX
  GOTO 999    ! Genx doesn't need to read turbulence solution file
#endif

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary turbulence file...'
  END IF ! global%verbLevel

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    currentTime = global%currentTime

    CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.turba', & 
                               region%iRegionGlobal,currentTime,iFileName)

    IF ( global%myProcid == MASTERPROC .AND. & 
         global%verbLevel > VERBOSE_NONE ) THEN                                         
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       region%iRegionGlobal  
      WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                          currentTime 
    END IF ! global%verbLevel                                       
  ELSE
    CALL BuildFileNameSteady(global,FILEDEST_INDIR,'.turba', & 
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
  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD", & 
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

  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# ROCTURB solution file' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM  
  
! -----------------------------------------------------------------------------
! Precision and range
! -----------------------------------------------------------------------------
    
  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# Precision and range' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM  
  
  precExpected  = PRECISION(1.0_RFREAL)
  rangeExpected = RANGE(1.0_RFREAL)
  
  READ(iFile) precActual,rangeActual
  IF ( precActual < precExpected .OR. rangeActual < rangeExpected ) THEN 
    CALL ErrorStop(global,ERR_PREC_RANGE,__LINE__)
  END IF ! precActual
  
! -----------------------------------------------------------------------------
! Initial residual (for RaNS), physical time, and energy subgrid quantities
! -----------------------------------------------------------------------------
  
  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# Initial residual' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
  END IF ! TRIM  
  
  READ(iFile) global%resInit
  
  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# Physical time' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
  END IF ! TRIM    
   
  READ(iFile) currentTime 

  IF ( global%flowType == FLOW_UNSTEADY ) THEN
    IF ( global%currentTime < 0.0_RFREAL ) THEN
      global%currentTime = currentTime
    ELSE
      WRITE(timeString1,'(1PE11.5)') global%currentTime
      WRITE(timeString2,'(1PE11.5)') currentTime          
      IF ( TRIM(timeString1) /= TRIM(timeString2) ) THEN
        CALL ErrorStop(global,ERR_TIME_SOLUTION,__LINE__,TRIM(iFileName))
      END IF ! global%currentTime 
    END IF ! global%currentTime
  END IF ! global%flowType  
  
  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# Esg1Sum' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
  END IF ! TRIM    
   
  READ(iFile) global%esg1Sum
  
  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# Esg4Sum' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
  END IF ! TRIM    
   
  READ(iFile) global%esg4Sum
  
! ==============================================================================
! Dimensions
! ==============================================================================
  
  grid => region%grid  
  
  nVarsExpected  = region%turbInput%nOutField
  nCellsExpected = grid%nCellsTot    
  
  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM
  
  READ(iFile) nCellsTot,nVars
  IF ( nCellsTot /= nCellsExpected ) THEN 
    WRITE(errorString,'(A,1X,I6,1X,A,1X,I6))') 'Specified:',nCellsTot, & 
                                               'but expected:',nCellsExpected
    CALL ErrorStop(global,ERR_INVALID_NCELLS,__LINE__,errorString)
  END IF ! nCellsExpected     
  
  IF ( nVars /= nVarsExpected ) THEN 
    WRITE(errorString,'(A,1X,I6,1X,A,1X,I6))') 'Specified:',nVars, & 
                                               'but expected:',nVarsExpected  
    CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
  END IF ! nVarsExpected   
  
! ==============================================================================
! Rest of file
! ==============================================================================

  iVars       = 0
  loopCounter = 0

  DO ! set up infinite loop
    loopCounter = loopCounter + 1
  
    READ(iFile) sectionString

    SELECT CASE ( TRIM(sectionString) ) 

! ------------------------------------------------------------------------------
!   Eddy viscosity
! ------------------------------------------------------------------------------

      CASE ( '# Eddy viscosity' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Eddy viscosity...'
        END IF ! global%verbLevel    

        IF (region%turbInput%modelClass == MODEL_LES) THEN
          tv => region%mixt%tv
      
          iVars = iVars + 1
          READ(iFile) (tv(TV_MIXT_MUET,j),j=1,grid%nCellsTot)

          IF (ASSOCIATED( region%turb%postv ) .eqv. .true.) THEN
            region%turb%postv(1,:) = tv(TV_MIXT_MUET,:)
          ENDIF
        ELSEIF (region%turbInput%modelClass == MODEL_RANS) THEN
          tcv => region%turb%cv
      
          iVars = iVars + 1
          READ(iFile) (tcv(CV_SA_NUTIL,j),j=1,grid%nCellsTot)

          IF (ASSOCIATED( region%turb%postv ).eqv..true.) THEN
            region%turb%postv(1,:) = tcv(CV_SA_NUTIL,:)
          ENDIF
        ENDIF

! ------------------------------------------------------------------------------
!   Total vorticity
! ------------------------------------------------------------------------------

      CASE ( '# Total vorticity' ) 
        IF ( global%myProcid == MASTERPROC .AND. & 
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Total vorticity...'
        END IF ! global%verbLevel         
      
        vort => region%turb%vort
      
        iVars = iVars + 1      
        READ(iFile) (vort(XCOORD,j),j=1,grid%nCellsTot)

! ------------------------------------------------------------------------------
!   RaNS or DES length scale
! ------------------------------------------------------------------------------

      CASE ( '# RANS length scale' ) 
        IF ( global%myProcid == MASTERPROC .AND. & 
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'RaNS length scale...'
        END IF ! global%verbLevel          
      
        lens => region%turb%lens
      
        iVars = iVars + 1      
        READ(iFile) (lens(j),j=1,grid%nCellsTot)

        IF (ASSOCIATED( region%turb%postv ).eqv..true.) THEN
          region%turb%postv(2,:) = lens(:)
        ENDIF

! ------------------------------------------------------------------------------
!     End marker
! ------------------------------------------------------------------------------ 
      
      CASE ( '# End' ) 
        IF ( global%myProcid == MASTERPROC .AND. & 
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
        END IF ! global%verbLevel           
      
        EXIT
      
! ------------------------------------------------------------------------------
!     Invalid section string
! ------------------------------------------------------------------------------ 
      
      CASE DEFAULT
        IF ( global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
        END IF ! verbosityLevel           
      
        CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)        
             
    END SELECT ! TRIM
  
! ------------------------------------------------------------------------------
!   Guard against infinite loop - might be unnecessary because of read errors?
! ------------------------------------------------------------------------------  
  
    IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
      CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
    END IF ! loopCounter
  
  END DO ! <empty>

! ==============================================================================
! Check and information about number of variables read
! ==============================================================================

  IF ( iVars /= nVars ) THEN 
    CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
  END IF ! iVar

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
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary turbulence file done.'
  END IF ! global%verbLevel

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction(global)
    
! ******************************************************************************
! End
! ******************************************************************************
 
END SUBROUTINE TURB_RFLU_ReadSolutionBinary


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_rFLU_ReadSolutionBinary.F90,v $
! Revision 1.6  2009/08/26 12:28:53  mtcampbe
! Ported to Hera.   Fixed logical expression syntax errors.  Replaced all
! IF (logical_variable)  with IF (logical_variable .eqv. .true.) as
! consistent with the specification.  Also changed: IF( ASSOCIATED(expr) )
! to IF ( ASSOCIATED(expr) .eqv. .true. ).   Intel compilers produce code
! which silently fails for some mal-formed expressions, so these changes
! are a net which should ensure that they are evaluated as intended.
!
! Revision 1.5  2008/12/06 08:44:45  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:56  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/01/12 09:44:14  wasistho
! copy to postv
!
! Revision 1.2  2004/06/16 20:01:31  haselbac
! Added use of ModBuildFileNames, cosmetics
!   
! Revision 1.1  2004/03/27 02:19:15  wasistho  
! added routines specific for Rocflu           
!
! ******************************************************************************







