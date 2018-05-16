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
! Purpose: Read flow file for particles in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: PLAG_RFLU_ReadSolutionBinary.F90,v 1.6 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_ReadSolutionBinary(pRegion)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region 
  USE ModPartLag, ONLY: t_plag,t_tile_plag     
  USE ModMPI

  USE PLAG_ModParameters

  USE ModBuildFileNames, ONLY: BuildFileNameUnsteady

  IMPLICIT NONE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: errorString,iFileName,sectionString,RCSIdentString, & 
                       timeString1,timeString2
  INTEGER :: errorFlag,iCont,iFile,ifl,iMass,iPatch,iVars,j,loopCounter,nCont, &
             nPcls,nPclsExpected,nVars,nVarsExpected,precActual,precExpected, &
             rangeActual,rangeExpected
  INTEGER, DIMENSION(:,:), POINTER :: pAiv
  REAL(RFREAL) :: currentTime
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pArv,pCv
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_plag), POINTER :: pPlag
  TYPE(t_tile_plag), POINTER :: pTilePlag  
  
! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion  
  
! ******************************************************************************
! Start, open file
! ******************************************************************************

  RCSIdentString = &
    '$RCSfile: PLAG_RFLU_ReadSolutionBinary.F90,v $ $Revision: 1.6 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_ReadSolutionBinary',&
  'PLAG_RFLU_ReadSolutionBinary.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary particle file...'
  END IF ! global%verbLevel

  CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.plag_sol', & 
                             pRegion%iRegionGlobal,global%currentTime, & 
                             iFileName)

  iFile = IF_SOLUT
  OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD",IOSTAT=errorFlag)
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
  IF ( TRIM(sectionString) /= '# ROCFLU particle file' ) THEN 
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
! Initial residual and physical time
! -----------------------------------------------------------------------------
    
  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# Physical time' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
  END IF ! TRIM    
   
  READ(iFile) currentTime 

#ifndef GENX  
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
#endif  

! ==============================================================================
! Dimensions
! ==============================================================================
  
  pGrid => pRegion%grid  
  pPlag => pRegion%plag
   
  nCont = pRegion%plagInput%nCont 
   
  nVarsExpected = 13 ! Hard-coded for now
  nPclsExpected = pPlag%nPcls    
  
  READ(iFile) sectionString
  IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
    CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
  END IF ! TRIM

  READ(iFile) nPcls,nVars
  
  IF ( nPcls /= nPclsExpected ) THEN 
    WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nPcls, & 
                                              'but expected:',nPclsExpected
    CALL ErrorStop(global,ERR_PLAG_INVALID_NPCLS,__LINE__,errorString)
  END IF ! nCellsExpected     
  
  IF ( nVars /= nVarsExpected ) THEN 
    WRITE(errorString,'(A,1X,I6,1X,A,1X,I6)') 'Specified:',nVars, & 
                                              'but expected:',nVarsExpected  
    CALL ErrorStop(global,ERR_INVALID_NVARS,__LINE__)
  END IF ! nVarsExpected    

! ==============================================================================
! Rest of file
! ==============================================================================

  iCont       = 0
  iVars       = 0
  loopCounter = 0

  DO ! set up infinite loop
    loopCounter = loopCounter + 1
  
    READ(iFile) sectionString

    SELECT CASE ( TRIM(sectionString) ) 

! ------------------------------------------------------------------------------
!     Particle x-momentum
! ------------------------------------------------------------------------------

      CASE ( '# Particle x-momentum' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle x-momentum...'
        END IF ! global%verbLevel    
      
        pCv => pRegion%plag%cv
      
        iVars = iVars + 1
        READ(iFile) (pCv(CV_PLAG_XMOM,j),j=1,pPlag%nPcls)

! ------------------------------------------------------------------------------
!     Particle y-momentum
! ------------------------------------------------------------------------------

      CASE ( '# Particle y-momentum' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle y-momentum...'
        END IF ! global%verbLevel    
      
        pCv => pRegion%plag%cv
      
        iVars = iVars + 1
        READ(iFile) (pCv(CV_PLAG_YMOM,j),j=1,pPlag%nPcls)

! ------------------------------------------------------------------------------
!     Particle z-momentum
! ------------------------------------------------------------------------------

      CASE ( '# Particle z-momentum' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle z-momentum...'
        END IF ! global%verbLevel    
      
        pCv => pRegion%plag%cv
      
        iVars = iVars + 1
        READ(iFile) (pCv(CV_PLAG_ZMOM,j),j=1,pPlag%nPcls)

! ------------------------------------------------------------------------------
!     Particle energy
! ------------------------------------------------------------------------------

      CASE ( '# Particle energy' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle energy...'
        END IF ! global%verbLevel    
      
        pCv => pRegion%plag%cv
      
        iVars = iVars + 1
        READ(iFile) (pCv(CV_PLAG_ENER,j),j=1,pPlag%nPcls)

! ------------------------------------------------------------------------------
!     Particle x-location
! ------------------------------------------------------------------------------

      CASE ( '# Particle x-location' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle x-location...'
        END IF ! global%verbLevel    
      
        pCv => pRegion%plag%cv
      
        iVars = iVars + 1
        READ(iFile) (pCv(CV_PLAG_XPOS,j),j=1,pPlag%nPcls)

! ------------------------------------------------------------------------------
!     Particle y-location
! ------------------------------------------------------------------------------

      CASE ( '# Particle y-location' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle y-location...'
        END IF ! global%verbLevel    
      
        pCv => pRegion%plag%cv
      
        iVars = iVars + 1
        READ(iFile) (pCv(CV_PLAG_YPOS,j),j=1,pPlag%nPcls)

! ------------------------------------------------------------------------------
!     Particle z-location
! ------------------------------------------------------------------------------

      CASE ( '# Particle z-location' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle z-location...'
        END IF ! global%verbLevel    
      
        pCv => pRegion%plag%cv
      
        iVars = iVars + 1
        READ(iFile) (pCv(CV_PLAG_ZPOS,j),j=1,pPlag%nPcls)

! ------------------------------------------------------------------------------
!     Particle vapor energy
! ------------------------------------------------------------------------------

      CASE ( '# Particle vapor energy' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle vapor energy...'
        END IF ! global%verbLevel    
      
        pCv => pRegion%plag%cv
      
        iVars = iVars + 1
        READ(iFile) (pCv(CV_PLAG_ENERVAPOR,j),j=1,pPlag%nPcls)

! ------------------------------------------------------------------------------
!     Particle mass
! ------------------------------------------------------------------------------

      CASE ( '# Particle mass' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle mass...'
        END IF ! global%verbLevel    
      
        pCv => pRegion%plag%cv
      
        IF ( iCont == 0 ) THEN 
          iVars = iVars + 1
        END IF ! iCont
              
        iCont = iCont + 1
        iMass = pPlag%cvPlagMass(iCont)
        
        READ(iFile) (pCv(iMass,j),j=1,pPlag%nPcls)

! ------------------------------------------------------------------------------
!     Particle superloading
! ------------------------------------------------------------------------------

      CASE ( '# Particle superloading' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle superloading...'
        END IF ! global%verbLevel    
      
        pArv => pRegion%plag%arv
      
        iVars = iVars + 1
        READ(iFile) (pArv(ARV_PLAG_SPLOAD,j),j=1,pPlag%nPcls)

! ------------------------------------------------------------------------------
!     Particle initial identifier
! ------------------------------------------------------------------------------

      CASE ( '# Particle initial identifier' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle initial identifier...'
        END IF ! global%verbLevel    
      
        pAiv => pRegion%plag%aiv
      
        iVars = iVars + 1
        READ(iFile) (pAiv(AIV_PLAG_PIDINI,j),j=1,pPlag%nPcls)

! ------------------------------------------------------------------------------
!     Particle initial region
! ------------------------------------------------------------------------------

      CASE ( '# Particle initial region' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle initial region...'
        END IF ! global%verbLevel    
      
        pAiv => pRegion%plag%aiv
      
        iVars = iVars + 1
        READ(iFile) (pAiv(AIV_PLAG_REGINI,j),j=1,pPlag%nPcls)
        
! ------------------------------------------------------------------------------
!     Particle cell
! ------------------------------------------------------------------------------

      CASE ( '# Particle cell' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Particle cell...'
        END IF ! global%verbLevel    
      
        pAiv => pRegion%plag%aiv
      
        iVars = iVars + 1
        READ(iFile) (pAiv(AIV_PLAG_ICELLS,j),j=1,pPlag%nPcls)        

! ------------------------------------------------------------------------------
!     Patch data
! ------------------------------------------------------------------------------

      CASE ( '# Patch data' )
        IF ( global%myProcid == MASTERPROC .AND. &
             global%verbLevel > VERBOSE_LOW ) THEN  
          WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Patch data...'
        END IF ! global%verbLevel                
      
        DO iPatch = 1,pGrid%nPatches
          pPatch => pRegion%patches(iPatch)

         IF ( (pPatch%bcType >= BC_INJECTION .AND. pPatch%bcType <= BC_INJECTION + BC_RANGE) .OR. &
              (pPatch%bcType >= BC_INFLOW    .AND. pPatch%bcType <= BC_INFLOW    + BC_RANGE)      )  THEN
            pTilePlag   => pPatch%tilePlag 

            DO ifl = 1,pPatch%nBFaces
              READ(iFile) pTilePlag%cv(CV_TILE_MOMNRM,ifl), & 
                          pTilePlag%cv(CV_TILE_ENER  ,ifl)
            END DO ! ifl 
            
            DO iCont = 1,nCont
              iMass = pTilePlag%cvTileMass(iCont)
              READ(iFile) (pTilePlag%cv(iMass,ifl), &
                           ifl=1,pPatch%nBFaces)
            END DO ! iCont            
            
            DO ifl = 1,pPatch%nBFaces
              READ(iFile) pTilePlag%dv(DV_TILE_COUNTDOWN,ifl), & 
                          pTilePlag%dv(DV_TILE_DIAM    ,ifl), & 
                          pTilePlag%dv(DV_TILE_SPLOAD  ,ifl)
            END DO ! ifl                                    
          END IF ! pPatch%bcType
        END DO ! iPatch

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
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary particle file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)
    
! ******************************************************************************
! End
! ******************************************************************************
 
END SUBROUTINE PLAG_RFLU_ReadSolutionBinary


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_ReadSolutionBinary.F90,v $
! Revision 1.6  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2007/03/31 23:56:17  haselbac
! Removed superfluous close parentheses
!
! Revision 1.3  2006/09/18 20:37:02  fnajjar
! Activated tile datastructure for inflow bc
!
! Revision 1.2  2005/01/21 17:23:10  fnajjar
! Included vapor energy in IO capability
!
! Revision 1.1  2004/08/23 23:06:53  fnajjar
! Initial revision
!
! ******************************************************************************







