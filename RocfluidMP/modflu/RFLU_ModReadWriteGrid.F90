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
! Purpose: Suite of routines to read and write grid files.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModReadWriteGrid.F90,v 1.10 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModReadWriteGrid

  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global 
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region  
  USE ModGrid, ONLY: t_grid
  USE ModMPI
    
  IMPLICIT NONE
     
  PRIVATE
  PUBLIC :: RFLU_ReadGridWrapper, & 
            RFLU_WriteGridWrapper

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModReadWriteGrid.F90,v $ $Revision: 1.10 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS









! ******************************************************************************
!
! Purpose: Read grid in ASCII ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. For GENX runs, read file from time zero if restarting. This is for 
!      convenience and will have to be changed once grid adaptation is used.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadGridASCII(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic, &
                                 BuildFileNameUnsteady 

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,sectionString,timeString1,timeString2
    INTEGER :: errorFlag,i,iFile,iPatch,j,k,loopCounter,nBCellsVirt, &
               nBQuadsTot,nBTrisTot,nBVertTot,nHexsTot,nPatches,nPrisTot, &
               nPyrsTot,nTetsTot,nVertTot,p,r
    REAL(RFREAL) :: currentTime
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_patch), POINTER :: pPatch  
    TYPE(t_global), POINTER :: global  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion   
  
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadGridASCII',&
  'RFLU_ModReadWriteGrid.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII grid file...'
    END IF ! global%verbLevel

    iFile = IF_GRID

    IF ( global%flowType == FLOW_UNSTEADY .AND. &
         (pRegion%mixtInput%moveGrid .EQV. .TRUE.) ) THEN   
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.grda', & 
                                 pRegion%iRegionGlobal,global%currentTime, & 
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN    
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal       
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                            global%currentTime
      END IF ! global%verbLevel                                      
    ELSE 
      CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.grda', & 
                              pRegion%iRegionGlobal,iFileName)    

        IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal  
      END IF ! global%verbLevel
    END IF ! global

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Read header stuff
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel 

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU grid file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM

! ------------------------------------------------------------------------------
!   Precision and range
! ------------------------------------------------------------------------------

    READ(iFile,'(A)') sectionString  
    IF ( TRIM(sectionString) /= '# Precision and range' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
    END IF ! TRIM

    READ(iFile,'(2(I8))') p,r
    IF ( p < PRECISION(1.0_RFREAL) .OR. r < RANGE(1.0_RFREAL) ) THEN 
      CALL ErrorStop(global,ERR_PREC_RANGE,__LINE__)
    END IF ! p

! ------------------------------------------------------------------------------
!   Initial residual and physical time
! ------------------------------------------------------------------------------
  
    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Physical time' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
    END IF ! TRIM    
  
    READ(iFile,'(E23.16)') currentTime 

    IF ( global%flowType == FLOW_UNSTEADY .AND. & 
         (pRegion%mixtInput%moveGrid .EQV. .TRUE.) ) THEN
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

! ==============================================================================
!   Dimensions
! ==============================================================================

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM
  
    pGrid => pRegion%grid 

    READ(iFile,'(5(I8))') nVertTot,nTetsTot,nHexsTot,nPrisTot,nPyrsTot

! ==============================================================================
!   Check dimensions (against those read from dimensions file)
! ==============================================================================

    IF ( nVertTot /= pGrid%nVertTot ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nVertTot

    IF ( nTetsTot /= pGrid%nTetsTot ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nTetsTot

    IF ( nHexsTot /= pGrid%nHexsTot ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nHexsTot

    IF ( nPrisTot /= pGrid%nPrisTot ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nPrisTot

    IF ( nPyrsTot /= pGrid%nPyrsTot ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nPyrsTot

! ==============================================================================
!   Rest of file
! ==============================================================================

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1
    
      READ(iFile,'(A)') sectionString
    
      SELECT CASE ( TRIM(sectionString) ) 
   
! ------------------------------------------------------------------------------
!       Coordinates
! ------------------------------------------------------------------------------

        CASE ( '# Coordinates' )       
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Coordinates...'
          END IF ! global%verbLevel

          DO i = 1,3
            READ(iFile,'(5(E23.16))') (pGrid%xyz(i,j),j=1,pGrid%nVertTot)
          END DO ! i      

! ------------------------------------------------------------------------------
!       Tetrahedra
! ------------------------------------------------------------------------------
      
        CASE ( '# Tetrahedra' )
          IF ( pGrid%nTetsTot == 0 ) THEN 
            CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
          END IF ! pGrid%nTetsTot 
       
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Tetrahedra...'
          END IF ! global%verbLevel      

          DO i = 1,4
            READ(iFile,'(10(I8))') (pGrid%tet2v(i,j),j=1,pGrid%nTetsTot)
          END DO ! i        
      
! ------------------------------------------------------------------------------
!       Hexahedra
! ------------------------------------------------------------------------------ 
      
        CASE ( '# Hexahedra' ) 
          IF ( pGrid%nHexsTot == 0 ) THEN 
            CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
          END IF ! pGrid%nHexsTot   

          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Hexahedra...'
          END IF ! global%verbLevel        

          DO i = 1,8
            READ(iFile,'(10(I8))') (pGrid%hex2v(i,j),j=1,pGrid%nHexsTot)
          END DO ! i

! ------------------------------------------------------------------------------
!       Prisms
! ------------------------------------------------------------------------------ 
      
        CASE ( '# Prisms' ) 
          IF ( pGrid%nPrisTot == 0 ) THEN 
            CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
          END IF ! pGrid%nPrisTot  

          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Prisms...'
          END IF ! global%verbLevel  

          DO i = 1,6
            READ(iFile,'(10(I8))') (pGrid%pri2v(i,j),j=1,pGrid%nPrisTot)
          END DO ! i 
 
! ------------------------------------------------------------------------------
!       Pyramids
! ------------------------------------------------------------------------------ 
      
        CASE ( '# Pyramids' ) 
          IF ( pGrid%nPyrsTot == 0 ) THEN 
            CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
          END IF ! pGrid%nPyrsTot       

          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Pyramids...'
          END IF ! global%verbLevel  

          DO i = 1,5
            READ(iFile,'(10(I8))') (pGrid%pyr2v(i,j),j=1,pGrid%nPyrsTot)
          END DO ! i      
      
! ------------------------------------------------------------------------------
!       Boundaries (format v1)
! ------------------------------------------------------------------------------ 
      
        CASE ( '# Boundaries' ) 
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN   
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary information...'   
          END IF ! global%verbLevel

          READ(iFile,*) nPatches

          IF ( nPatches /= pGrid%nPatches ) THEN 
            CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
          END IF ! nPatches

! ------- Loop over patches ----------------------------------------------------

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)  

! --------- Read dimensions 

            READ(iFile,'(2(I8))') nBTrisTot,nBQuadsTot        

! --------- Check dimensions

            IF ( nBTrisTot /= pPatch%nBTrisTot ) THEN 
              CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
            END IF ! nBTrisTot

            IF ( nBQuadsTot /= pPatch%nBQuadsTot ) THEN 
              CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
            END IF ! nBQuadsTot

! --------- Read data

            IF ( pPatch%nBTrisTot > 0 ) THEN
              DO j = 1,3
                READ(iFile,'(10(I8))') (pPatch%bTri2v(j,k),k=1,pPatch%nBTrisTot)
              END DO ! j
            END IF ! pPatch%nBTrisTot

            IF ( pPatch%nBQuadsTot > 0 ) THEN
              DO j = 1,4
                READ(iFile,'(10(I8))') (pPatch%bQuad2v(j,k), & 
                                        k=1,pPatch%nBQuadsTot)
              END DO ! j
            END IF ! pPatch%nBQuadsTot      
          END DO ! iPatch           
              
! ------------------------------------------------------------------------------
!       Boundaries (format v2)
! ------------------------------------------------------------------------------ 
      
        CASE ( '# Boundaries (v2)' ) 
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN   
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary information...'   
          END IF ! global%verbLevel

          READ(iFile,*) nPatches

          IF ( nPatches /= pGrid%nPatches ) THEN 
            CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
          END IF ! nPatches

! ------- Loop over patches ----------------------------------------------------

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)  

! --------- Read dimensions 

            READ(iFile,'(3(I8))') nBTrisTot,nBQuadsTot,nBCellsVirt        

! --------- Check dimensions

            IF ( nBTrisTot /= pPatch%nBTrisTot ) THEN 
              CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
            END IF ! nBTrisTot

            IF ( nBQuadsTot /= pPatch%nBQuadsTot ) THEN 
              CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
            END IF ! nBQuadsTot

            IF ( nBCellsVirt /= pPatch%nBCellsVirt ) THEN 
              CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
            END IF ! nBCellsVirt

! --------- Read data

            IF ( pPatch%nBTrisTot > 0 ) THEN
              DO j = 1,3
                READ(iFile,'(10(I8))') (pPatch%bTri2v(j,k),k=1,pPatch%nBTrisTot)
              END DO ! j
            END IF ! pPatch%nBTrisTot

            IF ( pPatch%nBQuadsTot > 0 ) THEN
              DO j = 1,4
                READ(iFile,'(10(I8))') (pPatch%bQuad2v(j,k), & 
                                        k=1,pPatch%nBQuadsTot)
              END DO ! j
            END IF ! pPatch%nBQuadsTot
            
            IF ( pPatch%nBCellsVirt > 0 ) THEN
              READ(iFile,'(10(I8))') (pPatch%bvc(k),k=1,pPatch%nBCellsVirt)
            END IF ! pPatch%nBCellsVirt       
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
!       Invalid section string
! ------------------------------------------------------------------------------ 
      
        CASE DEFAULT
          IF ( global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(3X,A)') sectionString
          END IF ! global%verbLevel        

          CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)

      END SELECT ! TRIM
  
! ==============================================================================
!     Guard against infinite loop - might be unnecessary because of read errors?
! ==============================================================================  

      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter  

    END DO ! <empty>

#ifdef CHECK_DATASTRUCT
! ******************************************************************************
!   Write out connectivity so can check data structure
! ******************************************************************************

    WRITE(STDOUT,'(A)') SOLVER_NAME
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Cell connectivity'
    
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Number of tetrahedra:', & 
                                   pGrid%nTetsTot
    DO i = 1,pGrid%nTetsTot
      WRITE(STDOUT,'(A,5(1X,I6))') SOLVER_NAME,i,pGrid%tet2v(1:4,i)
    END DO ! i  
    
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Number of hexahedra:', & 
                                   pGrid%nHexsTot
    DO i = 1,pGrid%nHexsTot
      WRITE(STDOUT,'(A,9(1X,I6))') SOLVER_NAME,i,pGrid%hex2v(1:8,i)
    END DO ! i      
    
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Number of prisms:', & 
                                    pGrid%nPrisTot           
    DO i = 1,pGrid%nPrisTot
      WRITE(STDOUT,'(A,7(1X,I6))') SOLVER_NAME,i,pGrid%pri2v(1:6,i)
    END DO ! i      
    
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Number of pyramids:', & 
                                   pGrid%nPyrsTot            
    DO i = 1,pGrid%nPyrsTot
      WRITE(STDOUT,'(A,6(1X,I6))') SOLVER_NAME,i,pGrid%pyr2v(1:5,i)
    END DO ! i      
    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Coordinates'
    DO i = 1,pGrid%nVertTot
        WRITE(STDOUT,'(A,1X,I6,3(1X,E18.9))') SOLVER_NAME,i,pGrid%xyz(1:3,i)
    END DO ! i
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
    WRITE(STDOUT,'(A)') SOLVER_NAME               
#endif 

! ******************************************************************************
!   Close file
! ******************************************************************************

    CLOSE(iFile,IOSTAT=errorFlag)   
    global%error = errorFlag   
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)  

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN    
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading ASCII grid file done.'
    END IF ! global%verbLevel           
  
  END SUBROUTINE RFLU_ReadGridASCII








! ******************************************************************************
!
! Purpose: Read grid in binary ROCFLU format.
!
! Description: None.
!
! Input:
!   pRegion     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. For GENX runs, read file from time zero if restarting. This is for 
!      convenience and will have to be changed once grid adaptation is used.
!
! ******************************************************************************

  SUBROUTINE RFLU_ReadGridBinary(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic, &
                                 BuildFileNameUnsteady  

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,sectionString,timeString1,timeString2 
    INTEGER :: errorFlag,i,iFile,iPatch,j,k,loopCounter,nBCellsVirt, &
               nBQuadsTot,nBTrisTot,nBVertTot,nHexsTot,nPatches,nPrisTot, &
               nPyrsTot,nTetsTot,nVertTot,p,r
    REAL(RFREAL) :: currentTime
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_patch), POINTER :: pPatch  
    TYPE(t_global), POINTER :: global  
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion   
  
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadGridBinary',&
  'RFLU_ModReadWriteGrid.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary grid file...'
    END IF ! global%verbLevel

    iFile = IF_GRID
  
    IF ( global%flowType == FLOW_UNSTEADY .AND. &
         (pRegion%mixtInput%moveGrid .EQV. .TRUE.) ) THEN 
      CALL BuildFileNameUnsteady(global,FILEDEST_INDIR,'.grd', & 
                                 pRegion%iRegionGlobal,global%currentTime, &
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN    
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal       
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                            global%currentTime
      END IF ! global%verbLevel                                      
    ELSE 
      CALL BuildFileNameBasic(global,FILEDEST_INDIR,'.grd', & 
                              pRegion%iRegionGlobal,iFileName)    

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal  
      END IF ! global%verbLevel
    END IF ! global

    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD",IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Read header stuff
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel 

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU grid file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM

! ------------------------------------------------------------------------------
!   Precision and range
! ------------------------------------------------------------------------------

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Precision and range' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
    END IF ! TRIM

    READ(iFile) p,r
    IF ( p < PRECISION(1.0_RFREAL) .OR. r < RANGE(1.0_RFREAL) ) THEN 
      CALL ErrorStop(global,ERR_PREC_RANGE,__LINE__)
    END IF ! p

! ------------------------------------------------------------------------------
!   Initial residual and physical time
! ------------------------------------------------------------------------------

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Physical time' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,iFileName)
    END IF ! TRIM    

    READ(iFile) currentTime 

    IF ( global%flowType == FLOW_UNSTEADY .AND. & 
         (pRegion%mixtInput%moveGrid .EQV. .TRUE.) ) THEN
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

! ==============================================================================
!   Dimensions
! ==============================================================================

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# Dimensions' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
    END IF ! TRIM

    pGrid => pRegion%grid 

    READ(iFile) nVertTot,nTetsTot,nHexsTot,nPrisTot,nPyrsTot

! ==============================================================================
!   Check dimensions (against those read from dimensions file)
! ==============================================================================

    IF ( nVertTot /= pGrid%nVertTot ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nVertTot

    IF ( nTetsTot /= pGrid%nTetsTot ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nTetsTot

    IF ( nHexsTot /= pGrid%nHexsTot ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nHexsTot

    IF ( nPrisTot /= pGrid%nPrisTot ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nPrisTot

    IF ( nPyrsTot /= pGrid%nPyrsTot ) THEN 
      CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
    END IF ! nPyrsTot

! ==============================================================================
!   Rest of file
! ==============================================================================

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile) sectionString

      SELECT CASE ( TRIM(sectionString) ) 
   
! ------------------------------------------------------------------------------
!       Coordinates
! ------------------------------------------------------------------------------

        CASE ( '# Coordinates' ) 
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Coordinates...'
          END IF ! global%verbLevel

          DO i = 1,3
            READ(iFile) (pGrid%xyz(i,j),j=1,pGrid%nVertTot)
          END DO ! i      

! ------------------------------------------------------------------------------
!       Tetrahedra
! ------------------------------------------------------------------------------
      
        CASE ( '# Tetrahedra' )
          IF ( pGrid%nTetsTot == 0 ) THEN 
            CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
          END IF ! pGrid%nTetsTot 

          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Tetrahedra...'
          END IF ! global%verbLevel      

          DO i = 1,4
            READ(iFile) (pGrid%tet2v(i,j),j=1,pGrid%nTetsTot)
          END DO ! i        
      
! ------------------------------------------------------------------------------
!     Hexahedra
! ------------------------------------------------------------------------------ 
      
        CASE ( '# Hexahedra' ) 
          IF ( pGrid%nHexsTot == 0 ) THEN 
            CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
          END IF ! pGrid%nHexsTot 

          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Hexahedra...'
          END IF ! global%verbLevel        

          DO i = 1,8
            READ(iFile) (pGrid%hex2v(i,j),j=1,pGrid%nHexsTot)
          END DO ! i

! ------------------------------------------------------------------------------
!     Prisms
! ------------------------------------------------------------------------------ 
      
        CASE ( '# Prisms' ) 
          IF ( pGrid%nPrisTot == 0 ) THEN 
            CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
          END IF ! pGrid%nPrisTot  

          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Prisms...'
          END IF ! global%verbLevel  

          DO i = 1,6
            READ(iFile) (pGrid%pri2v(i,j),j=1,pGrid%nPrisTot)
          END DO ! i 
 
! ------------------------------------------------------------------------------
!       Pyramids
! ------------------------------------------------------------------------------ 

        CASE ( '# Pyramids' ) 
          IF ( pGrid%nPyrsTot == 0 ) THEN 
            CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)
          END IF ! pGrid%nPyrsTot       

          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Pyramids...'
          END IF ! global%verbLevel  

          DO i = 1,5
            READ(iFile) (pGrid%pyr2v(i,j),j=1,pGrid%nPyrsTot)
          END DO ! i      
      
! ------------------------------------------------------------------------------
!       Boundaries (format v1)
! ------------------------------------------------------------------------------ 
      
        CASE ( '# Boundaries' ) 
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN   
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary information...'   
          END IF ! global%verbLevel

          READ(iFile) nPatches

          IF ( nPatches /= pGrid%nPatches ) THEN 
            CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
          END IF ! nPatches
        
! ------- Loop over patches ----------------------------------------------------

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)  

! --------- Read dimensions 

            READ(iFile) nBTrisTot,nBQuadsTot

! --------- Check dimensions

            IF ( nBTrisTot /= pPatch%nBTrisTot ) THEN 
              CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
            END IF ! nBTris

            IF ( nBQuadsTot /= pPatch%nBQuadsTot ) THEN 
              CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
            END IF ! nBQuads

! --------- Read data

            IF ( pPatch%nBTrisTot > 0 ) THEN
              DO j = 1,3
                READ(iFile) (pPatch%bTri2v(j,k),k=1,pPatch%nBTrisTot)
              END DO ! j
            END IF ! pPatch%nBTrisTot

            IF ( pPatch%nBQuadsTot > 0 ) THEN
              DO j = 1,4
                READ(iFile) (pPatch%bQuad2v(j,k),k=1,pPatch%nBQuadsTot)
              END DO ! j
            END IF ! pPatch%nBQuadsTot 
          END DO ! iPatch          
      
! ------------------------------------------------------------------------------
!       Boundaries (format v2)
! ------------------------------------------------------------------------------ 
      
        CASE ( '# Boundaries (v2)' ) 
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN   
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary information...'   
          END IF ! global%verbLevel

          READ(iFile) nPatches

          IF ( nPatches /= pGrid%nPatches ) THEN 
            CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
          END IF ! nPatches
        
! ------- Loop over patches ----------------------------------------------------

          DO iPatch = 1,pGrid%nPatches
            pPatch => pRegion%patches(iPatch)  

! --------- Read dimensions 

            READ(iFile) nBTrisTot,nBQuadsTot,nBCellsVirt     

! --------- Check dimensions

            IF ( nBTrisTot /= pPatch%nBTrisTot ) THEN 
              CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
            END IF ! nBTris

            IF ( nBQuadsTot /= pPatch%nBQuadsTot ) THEN 
              CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
            END IF ! nBQuads
            
            IF ( nBCellsVirt /= pPatch%nBCellsVirt ) THEN 
              CALL ErrorStop(global,ERR_DIMENS_INVALID,__LINE__)
            END IF ! nBCellsVirt            

! --------- Read data

            IF ( pPatch%nBTrisTot > 0 ) THEN
              DO j = 1,3
                READ(iFile) (pPatch%bTri2v(j,k),k=1,pPatch%nBTrisTot)
              END DO ! j
            END IF ! pPatch%nBTrisTot

            IF ( pPatch%nBQuadsTot > 0 ) THEN
              DO j = 1,4
                READ(iFile) (pPatch%bQuad2v(j,k),k=1,pPatch%nBQuadsTot)
              END DO ! j
            END IF ! pPatch%nBQuadsTot 

            IF ( pPatch%nBCellsVirt > 0 ) THEN
              READ(iFile) (pPatch%bvc(k),k=1,pPatch%nBCellsVirt)
            END IF ! pPatch%nBCellsVirt
          END DO ! iPatch      
      
! ------------------------------------------------------------------------------
!       End marker
! ------------------------------------------------------------------------------ 

        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel       

          EXIT
      
! ------------------------------------------------------------------------------
!       Invalid section string
! ------------------------------------------------------------------------------ 
      
        CASE DEFAULT
          IF ( global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(3X,A)') sectionString
          END IF ! global%verbLevel        

          CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)

      END SELECT ! TRIM
  
! ==============================================================================
!     Guard against infinite loop - might be unnecessary because of read errors?
! ==============================================================================  
  
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter  

    END DO ! <empty>

#ifdef CHECK_DATASTRUCT
! ******************************************************************************
!   Write out connectivity so can check data structure
! ******************************************************************************

    WRITE(STDOUT,'(A)') SOLVER_NAME
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### START CHECK OUTPUT ###'
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Cell connectivity'
    
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Number of tetrahedra:', & 
                                   pGrid%nTetsTot
    DO i = 1,pGrid%nTetsTot
      WRITE(STDOUT,'(A,5(1X,I6))') SOLVER_NAME,i,pGrid%tet2v(1:4,i)
    END DO ! i  
    
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Number of hexahedra:', & 
                                   pGrid%nHexsTot
    DO i = 1,pGrid%nHexsTot
      WRITE(STDOUT,'(A,9(1X,I6))') SOLVER_NAME,i,pGrid%hex2v(1:8,i)
    END DO ! i      
    
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Number of prisms:', & 
                                    pGrid%nPrisTot           
    DO i = 1,pGrid%nPrisTot
      WRITE(STDOUT,'(A,7(1X,I6))') SOLVER_NAME,i,pGrid%pri2v(1:6,i)
    END DO ! i      
    
    WRITE(STDOUT,'(A,1X,A,1X,I6)') SOLVER_NAME,'Number of pyramids:', & 
                                   pGrid%nPyrsTot            
    DO i = 1,pGrid%nPyrsTot
      WRITE(STDOUT,'(A,6(1X,I6))') SOLVER_NAME,i,pGrid%pyr2v(1:5,i)
    END DO ! i      
    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Coordinates'
    DO i = 1,pGrid%nVertTot
        WRITE(STDOUT,'(A,1X,I6,3(1X,E18.9))') SOLVER_NAME,i,pGrid%xyz(1:3,i)
    END DO ! i
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'### END CHECK OUTPUT ###'   
    WRITE(STDOUT,'(A)') SOLVER_NAME               
#endif 

! ******************************************************************************
!   Close file
! ******************************************************************************

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error
    
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN    
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading binary grid file done.'
    END IF ! global%verbLevel           

    CALL DeregisterFunction(global)  
  
  END SUBROUTINE RFLU_ReadGridBinary






! ******************************************************************************
!
! Purpose: Wrapper for reading of grid files in ROCFLU format.
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

  SUBROUTINE RFLU_ReadGridWrapper(pRegion)
  
#ifdef GENX
    USE RFLU_ModRocstarIO, ONLY: RFLU_GENX_DecideReadFile, & 
                              RFLU_GENX_GetGrid
#endif  
  
    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion    
   
! ==============================================================================
!   Local variables
! ==============================================================================

    TYPE(t_global), POINTER :: global
  
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_ReadGridWrapper',&
  'RFLU_ModReadWriteGrid.F90')

! ******************************************************************************
!   Read grid files
! ******************************************************************************
  
#ifdef GENX
    IF ( RFLU_GENX_DecideReadFile(global) .EQV. .FALSE. ) THEN 
#endif
      IF ( global%gridFormat == FORMAT_ASCII ) THEN
        CALL RFLU_ReadGridASCII(pRegion)
      ELSE IF ( global%gridFormat == FORMAT_BINARY ) THEN
        CALL RFLU_ReadGridBinary(pRegion)
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! global%gridFormat
#ifdef GENX
    ELSE 
      CALL RFLU_GENX_GetGrid(pRegion) 
    END IF ! RFLU_GENX_DecideReadFile         
#endif
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global) 

  END SUBROUTINE RFLU_ReadGridWrapper






! ******************************************************************************
!
! Purpose: Write grid in ASCII ROCFLU format.
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

  SUBROUTINE RFLU_WriteGridASCII(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic, &
                                 BuildFileNameUnsteady

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,i,iFile,iPatch,j,k
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion  
  
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteGridASCII',&
  'RFLU_ModReadWriteGrid.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII grid file...'
    END IF ! global%verbLevel

    iFile = IF_GRID

    IF ( global%flowType == FLOW_UNSTEADY .AND. &
         (pRegion%mixtInput%moveGrid .EQV. .TRUE.) ) THEN 
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.grda', & 
                                 pRegion%iRegionGlobal,global%currentTime, & 
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN        
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal          
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                            global%currentTime   
      END IF ! global%verbLevel
    ELSE 
      CALL BuildFileNameBasic(global,FILEDEST_OUTDIR,'.grda', & 
                              pRegion%iRegionGlobal,iFileName)    

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal                                 
      END IF ! global%verbLevel
    END IF ! global

    OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", & 
         IOSTAT=errorFlag)   
    global%error = errorFlag        
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel 

    sectionString = '# ROCFLU grid file'  
    WRITE(iFile,'(A)') TRIM(sectionString)  

    sectionString = '# Precision and range'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(2(I8))') PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)

    sectionString = '# Physical time'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(E23.16)') global%currentTime 

! ==============================================================================
!   Dimensions
! ==============================================================================
  
    pGrid => pRegion%grid  

    sectionString = '# Dimensions'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(5(I8))') pGrid%nVertTot,pGrid%nTetsTot,pGrid%nHexsTot, &
                           pGrid%nPrisTot,pGrid%nPyrsTot

! ==============================================================================
!   Coordinates
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Coordinates...'
    END IF ! global%verbLevel

    sectionString = '# Coordinates'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    DO i = 1,3
      WRITE(iFile,'(5(E23.16))') (pGrid%xyz(i,j),j=1,pGrid%nVertTot)
    END DO ! i

! ==============================================================================
!   Connectivity
! ==============================================================================

    IF ( pGrid%nTetsTot > 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN  
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Tetrahedra...'
      END IF ! global%verbLevel   

      sectionString = '# Tetrahedra'
      WRITE(iFile,'(A)') TRIM(sectionString)    
      DO i = 1,4
        WRITE(iFile,'(10(I8))') (pGrid%tet2v(i,j),j=1,pGrid%nTetsTot)
      END DO ! i
    END IF ! pGrid%nTetsTot

    IF ( pGrid%nHexsTot > 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_LOW ) THEN  
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Hexahedra...'
      END IF ! global%verbLevel   

      sectionString = '# Hexahedra'
      WRITE(iFile,'(A)') TRIM(sectionString)  
      DO i = 1,8
        WRITE(iFile,'(10(I8))') (pGrid%hex2v(i,j),j=1,pGrid%nHexsTot)
      END DO ! i
    END IF ! pGrid%nHexsTot

    IF ( pGrid%nPrisTot > 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_LOW ) THEN  
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Prisms...'
      END IF ! global%verbLevel   

      sectionString = '# Prisms'
      WRITE(iFile,'(A)') TRIM(sectionString)   
      DO i = 1,6
        WRITE(iFile,'(10(I8))') (pGrid%pri2v(i,j),j=1,pGrid%nPrisTot)
      END DO ! i
    END IF ! pGrid%nPrisTot

    IF ( pGrid%nPyrsTot > 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_LOW ) THEN  
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Pyramids...'
      END IF ! global%verbLevel   

      sectionString = '# Pyramids'
      WRITE(iFile,'(A)') TRIM(sectionString) 
      DO i = 1,5
        WRITE(iFile,'(10(I8))') (pGrid%pyr2v(i,j),j=1,pGrid%nPyrsTot)
      END DO ! i
    END IF ! pGrid%nPyrsTot  

! ==============================================================================
!   Boundary information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary information...'   
    END IF ! global%verbLevel

    sectionString = '# Boundaries (v2)'
    WRITE(iFile,'(A)') TRIM(sectionString) 
    WRITE(iFile,'(I8)') pGrid%nPatches

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)  

      WRITE(iFile,'(3(I8))') pPatch%nBTrisTot,pPatch%nBQuadsTot, &
                             pPatch%nBCellsVirt

      IF ( pPatch%nBTrisTot > 0 ) THEN 
        DO j = 1,3
          WRITE(iFile,'(10(I8))') (pPatch%bTri2v(j,k),k=1,pPatch%nBTrisTot)
        END DO ! j      
      END IF ! pPatch%nBTrisTot

      IF ( pPatch%nBQuadsTot > 0 ) THEN 
        DO j = 1,4
          WRITE(iFile,'(10(I8))') (pPatch%bQuad2v(j,k),k=1,pPatch%nBQuadsTot)
        END DO ! j      
      END IF ! pPatch%nBQuadsTot
      
      IF ( pPatch%nBCellsVirt > 0 ) THEN 
        WRITE(iFile,'(10(I8))') (pPatch%bvc(k),k=1,pPatch%nBCellsVirt)
      END IF ! pPatch%nBCellsVirt       
    END DO ! iPatch

! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile,'(A)') TRIM(sectionString) 

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%myProcid == MASTERPROC .AND. &
         global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error

! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII grid file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)
  
  END SUBROUTINE RFLU_WriteGridASCII







! ******************************************************************************
!
! Purpose: Write grid in binary ROCFLU format.
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

  SUBROUTINE RFLU_WriteGridBinary(pRegion)

    USE ModBuildFileNames, ONLY: BuildFileNameBasic, &
                                 BuildFileNameUnsteady

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
!   Local variables
! ==============================================================================

    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,i,iFile,iPatch,j,k
    TYPE(t_grid), POINTER :: pGrid  
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_global), POINTER :: global
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion  
  
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteGridBinary',&
  'RFLU_ModReadWriteGrid.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary grid file...'
    END IF ! global%verbLevel

    iFile = IF_GRID

    IF ( global%flowType == FLOW_UNSTEADY .AND. &
         (pRegion%mixtInput%moveGrid .EQV. .TRUE.) ) THEN 
      CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.grd', & 
                                 pRegion%iRegionGlobal,global%currentTime, & 
                                 iFileName)

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN        
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal          
        WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                            global%currentTime   
      END IF ! global%verbLevel
    ELSE 
      CALL BuildFileNameBasic(global,FILEDEST_OUTDIR,'.grd', & 
                              pRegion%iRegionGlobal,iFileName)    

      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal                                 
      END IF ! global%verbLevel
    END IF ! global

    OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", & 
         IOSTAT=errorFlag)
    global%error = errorFlag           
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
    END IF ! global%error

! ==============================================================================
!   Header and general information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel 

    sectionString = '# ROCFLU grid file'  
    WRITE(iFile) sectionString  

    sectionString = '# Precision and range'
    WRITE(iFile) sectionString 
    WRITE(iFile) PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)

    sectionString = '# Physical time'
    WRITE(iFile) sectionString 
    WRITE(iFile) global%currentTime 
  
! ==============================================================================
!   Dimensions
! ==============================================================================

    pGrid => pRegion%grid  

    sectionString = '# Dimensions'
    WRITE(iFile) sectionString  
    WRITE(iFile) pGrid%nVertTot,pGrid%nTetsTot,pGrid%nHexsTot,pGrid%nPrisTot, & 
                 pGrid%nPyrsTot

! ==============================================================================
!   Coordinates
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Coordinates...'
    END IF ! global%verbLevel

    sectionString = '# Coordinates'
    WRITE(iFile) sectionString
    DO i = 1,3
      WRITE(iFile) (pGrid%xyz(i,j),j=1,pGrid%nVertTot)
    END DO ! i

! ==============================================================================
!   Connectivity
! ==============================================================================

    IF ( pGrid%nTetsTot > 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_NONE ) THEN  
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Tetrahedra...'
      END IF ! global%verbLevel   

      sectionString = '# Tetrahedra'
      WRITE(iFile) sectionString   
      DO i = 1,4
        WRITE(iFile) (pGrid%tet2v(i,j),j=1,pGrid%nTetsTot)
      END DO ! i
    END IF ! pGrid%nTetsTot

    IF ( pGrid%nHexsTot > 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_LOW ) THEN  
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Hexahedra...'
      END IF ! global%verbLevel   

      sectionString = '# Hexahedra'
      WRITE(iFile) sectionString  
      DO i = 1,8
        WRITE(iFile) (pGrid%hex2v(i,j),j=1,pGrid%nHexsTot)
      END DO ! i
    END IF ! pGrid%nHexsTot

    IF ( pGrid%nPrisTot > 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_LOW ) THEN  
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Prisms...'
      END IF ! global%verbLevel   

      sectionString = '# Prisms'
      WRITE(iFile) sectionString     
      DO i = 1,6
        WRITE(iFile) (pGrid%pri2v(i,j),j=1,pGrid%nPrisTot)
      END DO ! i
    END IF ! pGrid%nPrisTot

    IF ( pGrid%nPyrsTot > 0 ) THEN
      IF ( global%myProcid == MASTERPROC .AND. &
           global%verbLevel > VERBOSE_LOW ) THEN  
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Pyramids...'
      END IF ! global%verbLevel   

      sectionString = '# Pyramids'
      WRITE(iFile) sectionString  
      DO i = 1,5
        WRITE(iFile) (pGrid%pyr2v(i,j),j=1,pGrid%nPyrsTot)
      END DO ! i
    END IF ! pGrid%nPyrsTot  

! ==============================================================================
!   Boundary information
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN   
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Boundary information...'   
    END IF ! global%verbLevel

    sectionString = '# Boundaries (v2)'
    WRITE(iFile) sectionString  
    WRITE(iFile) pGrid%nPatches

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)  

      WRITE(iFile) pPatch%nBTrisTot,pPatch%nBQuadsTot,pPatch%nBCellsVirt

      IF ( pPatch%nBTrisTot > 0 ) THEN 
        DO j = 1,3
          WRITE(iFile) (pPatch%bTri2v(j,k),k=1,pPatch%nBTrisTot)
        END DO ! pPatch%nBTrisTot      
      END IF ! bound

      IF ( pPatch%nBQuadsTot > 0 ) THEN 
        DO j = 1,4
          WRITE(iFile) (pPatch%bQuad2v(j,k),k=1,pPatch%nBQuadsTot)
        END DO ! j      
      END IF ! pPatch%nBQuadsTot 
      
      IF ( pPatch%nBCellsVirt > 0 ) THEN 
        WRITE(iFile) (pPatch%bvc(k),k=1,pPatch%nBCellsVirt)
      END IF ! pPatch%nBCellsVirt     
    END DO ! iPatch
    
! ==============================================================================
!   End marker
! ==============================================================================

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile) sectionString 

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%myProcid == MASTERPROC .AND. &
         global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
    END IF ! global%error
   
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing binary grid file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_WriteGridBinary






! ******************************************************************************
!
! Purpose: Wrapper for writing of grid files in ROCFLU format.
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

  SUBROUTINE RFLU_WriteGridWrapper(pRegion)

#ifdef GENX
    USE RFLU_ModRocstarIO, ONLY: RFLU_GENX_DecideWriteFile, & 
                              RFLU_GENX_PutGrid
#endif

    IMPLICIT NONE
  
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  
   
! ==============================================================================
!   Local variables
! ==============================================================================

    TYPE(t_global), POINTER :: global
  
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion  
  
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_WriteGridWrapper',&
  'RFLU_ModReadWriteGrid.F90')

! ******************************************************************************
!   Read solution files
! ******************************************************************************
     
#ifdef GENX
    IF ( RFLU_GENX_DecideWriteFile(global) .EQV. .FALSE. ) THEN     
#endif
      IF ( global%gridFormat == FORMAT_ASCII ) THEN
        CALL RFLU_WriteGridASCII(pRegion)
      ELSE IF ( global%gridFormat == FORMAT_BINARY ) THEN
        CALL RFLU_WriteGridBinary(pRegion)
      ELSE
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END IF ! global%gridFormat
#ifdef GENX
    ELSE 
        CALL RFLU_GENX_PutGrid(pRegion) 
    END IF ! RFLU_GENX_DecideReadFile         
#endif    
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)
 
  END SUBROUTINE RFLU_WriteGridWrapper






! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModReadWriteGrid


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModReadWriteGrid.F90,v $
! Revision 1.10  2008/12/06 08:44:23  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:34  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.7  2006/03/30 20:50:14  haselbac
! Added CASEs for backward compatibility
!
! Revision 1.6  2006/03/25 21:55:41  haselbac
! Changes bcos of sype patches
!
! Revision 1.5  2005/09/14 15:53:03  haselbac
! Bug fix: Now get proper time printed when reading/writing mv grid
!
! Revision 1.4  2005/05/03 03:05:30  haselbac
! Bug fix in reading/writing of binary/ASCII files
!
! Revision 1.3  2004/11/03 17:04:06  haselbac
! Removed IO of vertex and cell flags, and code related to HACK_PERIODIC
!
! Revision 1.2  2004/10/19 19:28:27  haselbac
! Adapted to changes in GENX logic, rm r/w of bv and gs
!
! Revision 1.1  2004/07/06 15:14:31  haselbac
! Initial revision
!
! ******************************************************************************












