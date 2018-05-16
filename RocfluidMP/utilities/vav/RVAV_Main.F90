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
! Purpose: Verification And Validation (VAV) Tool to compare 2 DataStreams.
!
! Description: currently supported formats are:
!              - ASCII 
!              - Binary
!
! Input: case name from the list of arguments
!
! Output: to Standard Output.
!
! Notes: 
!
!******************************************************************************
!
! $Id: RVAV_Main.F90,v 1.3 2008/12/06 08:45:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

PROGRAM ROCVAV_Post

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE ModMixture, ONLY    : t_mixt
  USE ModInterfaces, ONLY : RFLO_ReadRegionTopology, RFLO_InitInputValues, &
        ReadInputFile, RFLO_DerivedInputValues, RFLO_GetDimensDummyNodes, &
        RFLO_GetDimensDummy, RFLO_GetNodeOffset, RFLO_GetCellOffset, &
        RFLO_ReadGridRegion, RFLO_ReadSolutionRegion
  USE ModMPI
  USE ModParameters
  USE RVAV_ModParameters
  USE RVAV_ModGlobal
  USE RVAV_ModInterfaces, ONLY : BuildVersionString, RVAV_ReadFileStream1, &
           RVAV_ReadFileStream2, RVAV_ExtractVariables, RVAV_ComputeError, &
           RVAV_PlotResults,     RVAV_ComputeSimilarField,                 &
           RVAV_ComputeAnalyticalSolution,                                 &
           RVAV_ReadInputFile
  USE RVAV_ModDataStruct
  IMPLICIT NONE

! ... loop variables
  INTEGER :: iReg, iLev, i, j, k

! ... local variables
  CHARACTER(CHRLEN) :: msg, verbosity, versionString, headerString

  TYPE(t_region) , POINTER :: regions(:)
  TYPE(t_grid)   , POINTER :: grid
  TYPE(t_mixt)   , POINTER :: mixt
  TYPE(t_region) , POINTER :: regionsS1(:), regionsS2(:) 
  TYPE(t_compare), POINTER :: RVAVcompare
  TYPE(t_global) , POINTER :: global
  
! ... we are setting the indCp and indMol values here
! ... if the Cp and Molecular weight changes from cell to cell you will
! ... need to use the computed values of indCp and indMol in this program

  INTEGER :: indCp = 0, indMol = 0
  INTEGER :: iCompare 
  INTEGER :: ibegS1,iendS1,ijumpS1
  INTEGER :: jbegS1,jendS1,jjumpS1
  INTEGER :: kbegS1,kendS1,kjumpS1
  INTEGER :: iCOffS1,ijCOffS1
  
  INTEGER :: ibegS2,iendS2,ijumpS2
  INTEGER :: jbegS2,jendS2,jjumpS2
  INTEGER :: kbegS2,kendS2,kjumpS2
  INTEGER :: iCOffS2,ijCOffS2
  
  INTEGER :: variableIndexS1, variableIndexS2
  INTEGER :: fileTypeS1, fileTypeS2
  INTEGER :: similarityTypeS2
  INTEGER :: iNodesS1, jNodesS1, kNodesS1
  INTEGER :: iNodesS2, jNodesS2, kNodesS2
  INTEGER :: iRegS1, iRegS2
  INTEGER :: margin, versionWidth, errorFlag
  INTEGER, PARAMETER :: headerWidth = 54

  REAL(RFREAL) :: evS1Min, evS2Min, evS1Max, evS2Max
  
!******************************************************************************

  ALLOCATE( global ) 
  
  global%nFunTree  = 0
  CALL RegisterFunction( global, 'ROCVAV_Post',&
  'RVAV_Main.F90' )

! initialize global parameters ------------------------------------------------

  global%verbLevel = VERBOSE_NONE

  global%flowType    = FLOW_STEADY  ! stationary flow
  global%currentTime = -1._RFREAL   ! no physical time set
  global%currentIter = -1           ! no iteration

  global%inDir  = './'              ! directory path
  global%outDir = './'

  global%nProcAlloc = 1
  global%myProcid   = MASTERPROC    ! default process number (not an MPI code)
  global%mpierr     = ERR_NONE
  global%error      = ERR_NONE

  global%pi  = 4._RFREAL*ATAN(1._RFREAL)
  global%rad = global%pi/180._RFREAL

! print header ----------------------------------------------------------------

#ifdef MPI
  CALL MPI_Init( global%mpierr )
  IF (global%mpierr /=0 ) CALL ErrorStop( global, ERR_MPI_TROUBLE,__LINE__ )
#endif

  CALL BuildVersionString( versionString )

  headerString = ' '
  versionWidth = LEN_TRIM(versionString)
  margin       = (headerWidth-versionWidth)/2
  headerString(margin+1:margin+versionWidth) = versionString(1:versionWidth)
  headerString(1:1) = '*'
  headerString(headerWidth:headerWidth) = '*'

  WRITE(STDOUT,'(/,A)') '  ******************************************************'
  WRITE(STDOUT,  '(A)') '  *                                                    *'
  WRITE(STDOUT,  '(A)') '  *      ROCVAV: Verification And Validation Tool      *'
  WRITE(STDOUT,  '(A)') '  *      ========================================      *'
  WRITE(STDOUT,  '(A)') '  *                                                    *'
  WRITE(STDOUT,  '(A)') '  '//TRIM(headerString)
  WRITE(STDOUT,  '(A)') '  *     Copyright (c) by the University of Illinois    *'
  WRITE(STDOUT,  '(A)') '  *                                                    *'
  WRITE(STDOUT,'(A,/)') '  ******************************************************'

! read argument list ----------------------------------------------------------

  CALL GETARG(1,global%casename)
  CALL GETARG(2,verbosity)

  IF (LEN_TRIM(global%casename)==0 .OR. LEN_TRIM(verbosity)==0) THEN
    WRITE(STDOUT,'(/,A,//,5(A,/))') &
      'Usage: rocvav <casename> <verbosity>', &
      '       verbosity = 0 - no output', &
      '                 = 1 - moderate output',&
      '                 = 2 - output all'
#ifdef MPI
    CALL MPI_Finalize( global%mpierr )
#endif
    STOP
  ENDIF
  
  READ(verbosity,*) global%verbLevel  
  
! set MultiGrid Start Level ---------------------------------------------------

  global%startLevel = 1

! read region topology --------------------------------------------------------

  IF (global%verbLevel /= VERBOSE_NONE) &
    WRITE(STDOUT,'(/,A)') 'RFLO Reading region topology ...'

  CALL RFLO_ReadRegionTopology( global,regionsS1 )

  DO iReg=1,global%nRegions
    regionsS1(iReg)%startLevel = global%startLevel
    regionsS1(iReg)%currLevel  = global%startLevel
    IF (regionsS1(iReg)%nGridLevels < regionsS1(iReg)%currLevel) THEN
      WRITE(msg,1000) iReg,global%startLevel
      CALL ErrorStop( global, ERR_GRID_LEVEL,__LINE__,msg )
    ENDIF
  ENDDO ! iReg

! get user parameters ---------------------------------------------------------

  IF (global%verbLevel /= VERBOSE_NONE) &
    WRITE(STDOUT,'(/,A)') 'RFLO Reading user input ...'

  CALL RFLO_InitInputValues( regionsS1 )
  CALL ReadInputFile( regionsS1 )
  CALL RFLO_DerivedInputValues( regionsS1 )

! reset Iteration Number or Timestamp to Zero

  IF (global%flowType == FLOW_STEADY) THEN
    global%currentIter = 0
  ELSE
    global%timeStamp   = 0.0_RFREAL
  ENDIF

! read Input File Pertinent to RocVAV

  IF (global%verbLevel /= VERBOSE_NONE) &
    WRITE(STDOUT,'(/,A)') 'RocVAV Reading user input ...'

  CALL RVAV_ReadInputFile( global )

! read File for Stream 1

  IF (global%verbLevel /= VERBOSE_NONE) &
    WRITE(STDOUT,'(/,A)') 'RocVAV Reading Stream1 Data ...'

  CALL RVAV_ReadFileStream1( regionsS1 )

! generate Analytical Solution File 

  IF (global%verbLevel /= VERBOSE_NONE) &
    WRITE(STDOUT,'(/,A)') 'RocVAV Compute Analytical Solution Data for Stream2...'

  CALL RVAV_ComputeAnalyticalSolution( globalRVAV%similarityTypeS2, regionsS1 )

! read File for Stream 2

  IF (global%verbLevel /= VERBOSE_NONE) &
    WRITE(STDOUT,'(/,A)') 'RocVAV Reading Stream2 Data ...'

  CALL RVAV_ReadFileStream2( regionsS1, regionsS2 )

! - extraction of variables to be compared and the compute errors

  DO iCompare=1,globalRVAV%nComparisons

    RVAVcompare => globalRVAV%RVAVcompare(iCompare)
    
    IF (global%verbLevel /= VERBOSE_NONE) THEN
      WRITE(STDOUT,'(/,A,I5.5,A)') 'In comparison = ',iCompare,'we are comparing:'
      WRITE(STDOUT,'(A,I5.5)') 'Stream1 block    = ' ,RVAVcompare%blockS1
      WRITE(STDOUT,'(A,I5.5)') 'Stream2 block    = ' ,RVAVcompare%blockS2
      WRITE(STDOUT,'(A,I5.5)') 'Stream1 variable = ' ,RVAVcompare%variableIndexS1
      WRITE(STDOUT,'(A,I5.5)') 'Stream2 variable = ' ,RVAVcompare%variableIndexS2
    END IF ! verbLevel
    
! - checking the contents of blockS1 and blockS2

    IF (globalRVAV%fileTypeS1 == FILE_COMPUTED .AND. &
        globalRVAV%fileTypeS2 == FILE_COMPUTED) THEN

      IF (RVAVcompare%blockS1 /= RVAVcompare%blockS2)THEN
        WRITE(STDOUT,'(/,A)')'Streams 1 and 2 are both Computed Results' 
        WRITE(STDOUT,'(A)')'Block numbers on Stream1 and Stream2 do not match'
        WRITE(STDOUT,'(A)')'RocVAV will abort'
        CALL ErrorStop( global, ERR_PREVIOUS_ERRORS,__LINE__ )
      ENDIF

    ENDIF ! fileTypeS1

    IF (RVAVcompare%operationS1 /= RVAVcompare%operationS2) THEN
      WRITE(STDOUT,'(/,A)')'Operations in Stream1 and Stream2 do not match'
      WRITE(STDOUT,'(A)')'RocVAV will abort'
      CALL ErrorStop( global, ERR_PREVIOUS_ERRORS,__LINE__ )
    ENDIF

    ibegS1  = RVAVcompare%ibegS1
    iendS1  = RVAVcompare%iendS1
    ijumpS1 = RVAVcompare%ijumpS1
    
    jbegS1  = RVAVcompare%jbegS1
    jendS1  = RVAVcompare%jendS1
    jjumpS1 = RVAVcompare%jjumpS1
    
    kbegS1  = RVAVcompare%kbegS1
    kendS1  = RVAVcompare%kendS1
    kjumpS1 = RVAVcompare%kjumpS1

    variableIndexS1 = RVAVcompare%variableIndexS1  
    fileTypeS1      = globalRVAV%fileTypeS1

    ibegS2  = RVAVcompare%ibegS2
    iendS2  = RVAVcompare%iendS2
    ijumpS2 = RVAVcompare%ijumpS2

    jbegS2  = RVAVcompare%jbegS2
    jendS2  = RVAVcompare%jendS2
    jjumpS2 = RVAVcompare%jjumpS2

    kbegS2  = RVAVcompare%kbegS2
    kendS2  = RVAVcompare%kendS2
    kjumpS2 = RVAVcompare%kjumpS2

    iCOffS1  = globalRVAV%iCOffS1
    ijCOffS1 = globalRVAV%iCOffS1
    
    iCOffS2  = globalRVAV%iCOffS2
    ijCOffS2 = globalRVAV%iCOffS2
    
    variableIndexS2 = RVAVcompare%variableIndexS2
    fileTypeS2      = globalRVAV%fileTypeS2
    similarityTypeS2 = globalRVAV%similarityTypeS2
    
    iNodesS1 = INT(REAL( (iendS1-ibegS1)/ijumpS1,KIND=RFREAL))+1
    jNodesS1 = INT(REAL( (jendS1-jbegS1)/jjumpS1,KIND=RFREAL))+1
    kNodesS1 = INT(REAL( (kendS1-kbegS1)/kjumpS1,KIND=RFREAL))+1

    iNodesS2 = INT(REAL( (iendS2-ibegS2)/ijumpS2,KIND=RFREAL))+1
    jNodesS2 = INT(REAL( (jendS2-jbegS2)/jjumpS2,KIND=RFREAL))+1
    kNodesS2 = INT(REAL( (kendS2-kbegS2)/kjumpS2,KIND=RFREAL))+1

    IF (global%verbLevel /= VERBOSE_NONE) THEN
      WRITE(STDOUT,'(/,A,3(I5,3X))') 'Stream1 iNodes,jNodes,kNodes = ', &
                                      iNodesS1,jNodesS1,kNodesS1
      WRITE(STDOUT,'(A,3(I5,3X))')  'Stream2 iNodes,jNodes,kNodes = ', &
                                     iNodesS2,jNodesS2,kNodesS2
    END IF ! verbLevel 

! - Allocate evS1 and evS2 

    IF (.NOT. ASSOCIATED(globalRVAV%evS1))THEN
      ALLOCATE( globalRVAV%evS1(iNodesS1,jNodesS1,kNodesS1),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE, __LINE__ )
    ENDIF  

    IF (.NOT. ASSOCIATED(globalRVAV%evS2))THEN
      ALLOCATE( globalRVAV%evS2(iNodesS2,jNodesS2,kNodesS2),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global, ERR_ALLOCATE, __LINE__ )
    ENDIF  

! - initialize evS1 and evS2

    globalRVAV%evS1 = 0.0_RFREAL
    globalRVAV%evS2 = 0.0_RFREAL

! - extract variables from Stream1

    iRegS1 = RVAVcompare%blockS1

    IF (global%verbLevel /= VERBOSE_NONE) &
      WRITE(STDOUT,'(/,A)') 'Entering RVAV_ExtractVariables-EVS1'

    CALL RVAV_ExtractVariables( global, regionsS1(iRegS1), &
                                ibegS1,iendS1,ijumpS1,     &
                                jbegS1,jendS1,jjumpS1,     &
                                kbegS1,kendS1,kjumpS1,     &
                                iCOffS1, ijCOffS1,         &
                                variableIndexS1,           &
                                fileTypeS1,                &
                                indCp,indMol,globalRVAV%evS1 )

    IF (global%verbLevel /= VERBOSE_NONE) &
      WRITE(STDOUT,'(/,A)') 'Exiting RVAV_ExtractVariables-EVS1'

! - extract variables from Stream2

    iRegS2 = RVAVcompare%blockS2

    IF (global%verbLevel /= VERBOSE_NONE) &
      WRITE(STDOUT,'(/,A)') 'Entering RVAV_ExtractVariables-EVS2'

    CALL RVAV_ExtractVariables( global, regionsS2(iRegS2), &
                                ibegS2,iendS2,ijumpS2,  &
                                jbegS2,jendS2,jjumpS2,  &
                                kbegS2,kendS2,kjumpS2,  &
                                iCOffS2, ijCOffS2,      &
                                variableIndexS2,        &
                                fileTypeS2,             &
                                indCp,indMol,globalRVAV%evS2 )
                                
    IF (global%verbLevel /= VERBOSE_NONE) &
      WRITE(STDOUT,'(/,A)') 'Exiting RVAV_ExtractVariables-EVS2'
      
! - apply similarity analysis on evS1 if needed

    IF ( fileTypeS2       == FILE_ANALYTICAL .AND. &
         similarityTypeS2 /= 0 ) THEN
      IF ( global%verbLevel/=VERBOSE_NONE ) &
        WRITE(STDOUT,'(/,A)') 'Entering RVAV_ComputeSimilarField-EVS2'

       CALL RVAV_ComputeSimilarField( global,                     &
                                      iNodesS1,jNodesS1,kNodesS1, &
                                      similarityTypeS2,           &
                                      variableIndexS1,            &
                                      globalRVAV%evS1)

      IF ( global%verbLevel/=VERBOSE_NONE ) &
        WRITE(STDOUT,'(/,A)') 'Entering RVAV_ComputeSimilarField-EVS2'

    END IF ! fileTypeS2
    
! - extract Min Max values from evS1, evS2

    IF (global%verbLevel /= VERBOSE_NONE) THEN
      evS1Min = +1.0E+30_RFREAL
      evS1Max = -1.0E+30_RFREAL
      evS2Min = +1.0E+30_RFREAL
      evS2Max = -1.0E+30_RFREAL

      DO k=1,kNodesS1
        DO j=1,jNodesS1
          DO i=1,iNodesS1
            evS1Min = MIN(globalRVAV%evS1(i,j,k), evS1Min)
            evS1Max = MAX(globalRVAV%evS1(i,j,k), evS1Max)
          ENDDO ! i
        ENDDO   ! j
      ENDDO     ! k

      DO k=1,kNodesS2
        DO j=1,jNodesS2
          DO i=1,iNodesS2
            evS2Min = MIN(globalRVAV%evS2(i,j,k), evS2Min)
            evS2Max = MAX(globalRVAV%evS2(i,j,k), evS2Max)
          ENDDO ! i
        ENDDO   ! j
      ENDDO     ! k
    
      WRITE(STDOUT,'(A,2E14.5)') 'MIN-MAX of EVS1',evS1min,evS1max
      WRITE(STDOUT,'(A,2E14.5)') 'MIN-MAX of EVS2',evS2min,evS2max
      
    END IF ! verbLevel
    
! - perform operations according to input flag

    iNodesS1 = INT(REAL((iendS1-ibegS1)/ijumpS1,KIND=RFREAL)) + 1
    jNodesS1 = INT(REAL((jendS1-jbegS1)/jjumpS1,KIND=RFREAL)) + 1
    kNodesS1 = INT(REAL((kendS1-kbegS1)/kjumpS1,KIND=RFREAL)) + 1

    iNodesS2 = INT(REAL((iendS2-ibegS2)/ijumpS2,KIND=RFREAL)) + 1
    jNodesS2 = INT(REAL((jendS2-jbegS2)/jjumpS2,KIND=RFREAL)) + 1
    kNodesS2 = INT(REAL((kendS2-kbegS2)/kjumpS2,KIND=RFREAL)) + 1

    IF (global%verbLevel /= VERBOSE_NONE) THEN
      WRITE(STDOUT,'(/,A,3(I5,3X))') 'Stream1 iNodes,jNodes,kNodes = ', &
                                      iNodesS1,jNodesS1,kNodesS1
      WRITE(STDOUT,'(/,A,3(I5,3X))') 'Stream2 iNodes,jNodes,kNodes = ', &
                                     iNodesS2,jNodesS2,kNodesS2
    END IF ! verbLevel
    
    IF (RVAVcompare%operationS1 ==  COMPUTE_ERRORS_ONLY) THEN 
      IF (iNodesS1 /= iNodesS2) THEN
        WRITE(STDOUT,'(/,A,2(I5,3X))') &
          'Non Matching Number of I-Nodes for Streams 1 and 2: ', &
          iNodesS1,iNodesS2
        CALL ErrorStop( global, ERR_PREVIOUS_ERRORS,__LINE__ )
      ENDIF ! iNodes
      
      IF (jNodesS1 /= jNodesS2) THEN
        WRITE(STDOUT,'(/,A,2(I5,3X))') &
          'Non Matching Number of J-Nodes for Streams 1 and 2: ', &
          jNodesS1,jNodesS2
        CALL ErrorStop( global, ERR_PREVIOUS_ERRORS,__LINE__ )
      ENDIF ! jNodes
      
      IF (kNodesS1 /= kNodesS2) THEN
        WRITE(STDOUT,'(/,A,2(I5,3X))') &
          'Non Matching Number of K-Nodes for Streams 1 and 2: ', &
          kNodesS1,kNodesS2
        CALL ErrorStop( global, ERR_PREVIOUS_ERRORS,__LINE__ )
      ENDIF ! kNodes

    ENDIF   ! operationS1   

    IF (RVAVcompare%operationS1 == COMPUTE_ERRORS_ONLY) THEN  
      CALL RVAV_ComputeError( global,iCompare,iNodesS1,jNodesS1,kNodesS1 )
    ENDIF   ! operationsS1

!    IF (RVAVcompare%operationS1 ==  PLOT_ERRORS_ONLY) THEN  
!      CALL RVAV_PlotResults( global, iCompare,iNodesS1,jNodesS1,kNodesS1 )
!    ENDIF   ! operationsS1

! - Deallocate evS1 and evS2 

    DEALLOCATE( globalRVAV%evS1,stat=errorFlag )
    global%error = errorFlag
    IF( global%error /= 0 ) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )

    DEALLOCATE( globalRVAV%evS2,stat=errorFlag )
    global%error = errorFlag
    IF( global%error /= 0 ) CALL ErrorStop( global,ERR_DEALLOCATE,__LINE__ )  

  ENDDO   ! iCompare

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

  WRITE(STDOUT,'(/,A)') 'RocVAV Finished.'

#ifdef MPI
  CALL MPI_Finalize( global%mpierr )
#endif

1000 FORMAT('Region ',I5,', grid level= ',I2,'.')

END PROGRAM ROCVAV_Post

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RVAV_Main.F90,v $
! Revision 1.3  2008/12/06 08:45:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:19  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 22:43:24  fnajjar
! Initial revision after changing case
!
! Revision 1.16  2003/11/20 16:40:41  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.12  2003/05/15 02:57:08  jblazek
! Inlined index function.
!
! Revision 1.11  2002/10/12 03:20:51  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.10  2002/09/20 22:22:37  jblazek
! Finalized integration into GenX.
!
! Revision 1.9  2002/09/11 16:20:21  jblazek
! Added directory path to input/output files (needed for GENX).
!
! Revision 1.8  2002/09/10 00:01:45  f-najjar
! Variable global moved into regions()
!
! Revision 1.7  2002/08/16 21:33:48  jblazek
! Changed interface to MixtureProperties.
!
! Revision 1.6  2002/07/12 21:50:08  jblazek
! Added tool to split single grid into multiple regions.
!
! Revision 1.5  2002/06/19 14:40:14  f-najjar
! Included verbLevel calls for cleanup
!
! Revision 1.4  2002/06/18 03:18:20  f-najjar
! Included RVAV_computeAnalyticalSolution
!
! Revision 1.3  2002/06/17 17:02:07  f-najjar
! Fix Calling sequence for ExtractVariables
!
! Revision 1.2  2002/06/14 17:00:52  jblazek
! Added version string.
!
! Revision 1.1.1.1  2002/06/03 21:41:29  f-najjar
! Initial Import of RocVaV
!
!******************************************************************************







