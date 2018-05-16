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
! Purpose: Collection of routines for particle surface statistics.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: PLAG_ModSurfStats.F90,v 1.13 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_ModSurfStats

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY: t_plag, t_surfstats_plag
  USE ModBndPatch, ONLY: t_patch  
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI
  USE PLAG_ModParameters

  USE ModBuildFileNames, ONLY: BuildFileNameSteady, &
                               BuildFileNameUnsteady

  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: PLAG_ModSurfStats.F90,v $ $Revision: 1.13 $'        

! ==============================================================================
! Public functions
! ==============================================================================
  
  PUBLIC :: PLAG_CreateSurfStats, & 
            PLAG_DestroySurfStats, & 
            PLAG_GatherSurfStats, &
            PLAG_ReadSurfStatsWrapper, &
            PLAG_WriteSurfStatsWrapper

! ==============================================================================
! Private functions
! ==============================================================================

  PRIVATE :: PLAG_CloseSurfStats, & 
             PLAG_CreateSurfStatsKernel, &
             PLAG_DestroySurfStatsKernel, &
             PLAG_InitSurfStats, &
             PLAG_NullifySurfStats, &
             PLAG_OpenSurfStatsASCII, &
             PLAG_OpenSurfStatsBinary, &
             PLAG_ReadSurfStatsASCII, & 
             PLAG_ReadSurfStatsBinary, & 
             PLAG_ReadSurfStatsKernelASCII, & 
             PLAG_ReadSurfStatsKernelBinary, & 
             PLAG_SetBinIndex, &
             PLAG_WriteSurfStatsASCII, & 
             PLAG_WriteSurfStatsBinary, &
             PLAG_WriteSurfStatsKernelASCII, & 
             PLAG_WriteSurfStatsKernelBinary

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Close surface statitstics file.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_CloseSurfStats(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,iFile
    TYPE(t_global), POINTER :: global    
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_CloseSurfStats',&
  'PLAG_ModSurfStats.F90')
             
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Closing PLAG surface statistics file...'
    END IF ! global%verbLevel             

    iFile = IF_PLAG_SURF_STATS
        
! ******************************************************************************
!   Close file
! ******************************************************************************

    CLOSE(iFile,IOSTAT=errorFlag)
    global%error = errorFlag      
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__)
    END IF ! global%error

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Closing PLAG surface statistics file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_CloseSurfStats  
  
  
  
  
  

! *******************************************************************************
!
! Purpose: Create surface statistics.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_CreateSurfStats(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,iPatch
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_CreateSurfStats',&
  'PLAG_ModSurfStats.F90')
    
    pGrid => pRegion%grid
        
! ******************************************************************************
!   Allocate memory
! ******************************************************************************
    
! ==============================================================================
!   Patch data
! ==============================================================================

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN        
        ALLOCATE(pPatch%statsPlag(pPatch%nBFaces),STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPatch%statsPlag')
        END IF ! global%error
     
        CALL PLAG_CreateSurfStatsKernel(pRegion,pPatch%statsPlag)                           
      END IF ! pPatch%plotStatsFlag
    END DO ! iPatch
    
! ==============================================================================
!   Plane data
! ==============================================================================

! ******************************************************************************
!   Initialize memory
! ******************************************************************************
    
! ==============================================================================
!   Patch data
! ==============================================================================

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)
      IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
        CALL PLAG_InitSurfStats(pRegion,pPatch%statsPlag)                           
      END IF ! pPatch%plotStatsFlag
    END DO ! iPatch
    
! ==============================================================================
!   Plane data
! ==============================================================================

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_CreateSurfStats
  
  
  
  
  

! *******************************************************************************
!
! Purpose: Kernel to create surface statistics.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pStatsPlag          Pointer to statsPlag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_CreateSurfStatsKernel(pRegion,pStatsPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
     
    TYPE(t_region), POINTER :: pRegion
    TYPE(t_surfstats_plag), DIMENSION(:), POINTER :: pStatsPlag
 
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ifl,nBFaces,nBins,nCont,nVars

    TYPE(t_global), POINTER :: global    

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'PLAG_CreateSurfStatsKernel',&
  'PLAG_ModSurfStats.F90')

    nBFaces = SIZE(pStatsPlag,DIM=1)   

    nCont = pRegion%plagInput%nCont
    nVars = PLAG_SURF_STATS_LAST +nCont
     
! TEMPORARY
    nBins = 20
!    nBins = pRegion%plagInput%nBins
! END TEMPORARY

! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    DO ifl = 1,nBFaces
      ALLOCATE(pStatsPlag(ifl)%nHits(nBins),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pStatsPlag%nHits')
      END IF ! global%error
        
      ALLOCATE(pStatsPlag(ifl)%vars(nVars,nBins),STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pStatsPlag%vars')
      END IF ! global%error     
    END DO ! ifl                           

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_CreateSurfStatsKernel 






! *******************************************************************************
!
! Purpose: Destroy surface statistics.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_DestroySurfStats(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,iPatch
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_DestroySurfStats',&
  'PLAG_ModSurfStats.F90')
    
    pGrid => pRegion%grid
        
! ******************************************************************************
!   Deallocate memory
! ******************************************************************************
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN       
        CALL PLAG_DestroySurfStatsKernel(global,pPatch%statsPlag)

        DEALLOCATE(pPatch%statsPlag,STAT=errorFlag)
        global%error = errorFlag
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%statsPlag')
        END IF ! global%error
      END IF ! pPatch%plotStatsFlag  
    END DO ! iPatch

! ******************************************************************************
!   Nullify memory
! ******************************************************************************
    
    CALL PLAG_NullifySurfStats(pRegion)    

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_DestroySurfStats
  





! *******************************************************************************
!
! Purpose: Kernel to destroy surface statistics.
!
! Description: None.
!
! Input:
!   global              Pointer to global
!   pStatsPlag          Pointer to statsPlag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_DestroySurfStatsKernel(global,pStatsPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
     
    TYPE(t_global), POINTER :: global  
    TYPE(t_surfstats_plag), DIMENSION(:), POINTER :: pStatsPlag
 
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,ifl,nBFaces

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    CALL RegisterFunction(global,'PLAG_DestroySurfStatsKernel',&
  'PLAG_ModSurfStats.F90')

    nBFaces = SIZE(pStatsPlag,DIM=1)        

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************
        
    DO ifl = 1,nBFaces
      DEALLOCATE(pStatsPlag(ifl)%nHits,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pStatsPlag%nHits')
      END IF ! global%error
        
      DEALLOCATE(pStatsPlag(ifl)%vars,STAT=errorFlag)
      global%error = errorFlag
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPatch%statsPlag%vars')
      END IF ! global%error     
    END DO ! ifl                           

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_DestroySurfStatsKernel
  
  
  
  
  

! *******************************************************************************
!
! Purpose: Gather surface statistics.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPlag               Pointer to plag
!   pStatsPlag          Pointer to statsPlag
!   ifl                 Local face index  
!   iPcl                Particle index
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_GatherSurfStats(pRegion,pPlag,pStatsPlag,ifl,iPcl,thetaAngle)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion   
    TYPE(t_plag), POINTER :: pPlag
    TYPE(t_surfstats_plag), DIMENSION(:), POINTER :: pStatsPlag

    INTEGER,INTENT(IN) :: ifl,iPcl
    REAL(KIND=RFREAL),INTENT(IN) :: thetaAngle
    
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: binMethod,errorFlag,iBin,iCont,nBins,nCont
    REAL(KIND=RFREAL) :: diamBin,diamMinMicron,diamMaxMicron,diamPlagMicron, &
                         enerPlag,massPlag,momeMagnPlag

    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_GatherSurfStats',&
  'PLAG_ModSurfStats.F90')

    nCont = pRegion%plagInput%nCont

! TEMPORARY
    nBins = 20
!    nBins = pRegion%plagInput%nBins
! END TEMPORARY

! TEMPORARY
  binMethod = BIN_METHOD_LINEAR
! binMethod = pRegion%plagInput%binMethod
! END TEMPORARY

! TEMPORARY
  diamMinMicron =  10.0_RFREAL
  diamMaxMicron = 210.0_RFREAL
  diamBin       = (diamMaxMicron-diamMinMicron)/REAL(nBins,KIND=RFREAL)
! END TEMPORARY
  
  diamPlagMicron = pPlag%dv(DV_PLAG_DIAM,iPcl)*1.0E+6_RFREAL       
  
  massPlag = SUM( pPlag%cv(pPlag%cvPlagMass(:),iPcl) )
  
  momeMagnPlag = SQRT( pPlag%cv(CV_PLAG_XMOM,iPcl)**2.0_RFREAL &
                      +pPlag%cv(CV_PLAG_YMOM,iPcl)**2.0_RFREAL &
                      +pPlag%cv(CV_PLAG_ZMOM,iPcl)**2.0_RFREAL )
  
  enerPlag = pPlag%cv(CV_PLAG_ENER,iPcl)

! ******************************************************************************
!   Gather statistics
! ******************************************************************************
    
! ==============================================================================
!   Determine bin
! ==============================================================================

    iBin = PLAG_SetBinIndex(global,binMethod,nBins,diamPlagMicron,diamBin, &
                            diamMinMicron,diamMaxMicron)

! ==============================================================================
!   Cummulate data on bins
! ==============================================================================

    pStatsPlag(ifl)%nHits(iBin) = pStatsPlag(ifl)%nHits(iBin) +1
 
    pStatsPlag(ifl)%vars(PLAG_SURF_STATS_DIAM3,iBin) = &
      pStatsPlag(ifl)%vars(PLAG_SURF_STATS_DIAM3,iBin) +diamPlagMicron**3.0_RFREAL 

    pStatsPlag(ifl)%vars(PLAG_SURF_STATS_DIAM4,iBin) = &
      pStatsPlag(ifl)%vars(PLAG_SURF_STATS_DIAM4,iBin) +diamPlagMicron**4.0_RFREAL 

    pStatsPlag(ifl)%vars(PLAG_SURF_STATS_THETA,iBin) = &
      pStatsPlag(ifl)%vars(PLAG_SURF_STATS_THETA,iBin) +thetaAngle 

    pStatsPlag(ifl)%vars(PLAG_SURF_STATS_MOME1,iBin) = &
      pStatsPlag(ifl)%vars(PLAG_SURF_STATS_MOME1,iBin) +momeMagnPlag

    pStatsPlag(ifl)%vars(PLAG_SURF_STATS_MOME2,iBin) = &
      pStatsPlag(ifl)%vars(PLAG_SURF_STATS_MOME2,iBin) +momeMagnPlag**2.0_RFREAL

    pStatsPlag(ifl)%vars(PLAG_SURF_STATS_ENER,iBin) = &
      pStatsPlag(ifl)%vars(PLAG_SURF_STATS_ENER,iBin) +enerPlag

    pStatsPlag(ifl)%vars(PLAG_SURF_STATS_MASS,iBin) = &
      pStatsPlag(ifl)%vars(PLAG_SURF_STATS_MASS,iBin) +massPlag

    DO iCont = 1, nCont
      pStatsPlag(ifl)%vars(PLAG_SURF_STATS_LAST+iCont,iBin) = &
       pStatsPlag(ifl)%vars(PLAG_SURF_STATS_LAST+iCont,iBin)  &
      +pPlag%cv(pPlag%cvPlagMass(iCont),iPcl)
    END DO ! iCont 

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_GatherSurfStats
  





! *******************************************************************************
!
! Purpose: Initialize surface statistics.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pStatsPlag          Pointer to statsPlag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_InitSurfStats(pRegion,pStatsPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
     
    TYPE(t_region), POINTER :: pRegion  
    TYPE(t_surfstats_plag), DIMENSION(:), POINTER :: pStatsPlag
  
! ==============================================================================
!   Locals
! ==============================================================================
 
    INTEGER :: errorFlag,iBin,ifl,iPatch,iVar,nBFaces,nBins,nCont,nVars

    TYPE(t_global), POINTER :: global  

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'PLAG_InitSurfStats',&
  'PLAG_ModSurfStats.F90')

    nBFaces = SIZE(pStatsPlag,DIM=1)          

    nCont = pRegion%plagInput%nCont
    nVars = PLAG_SURF_STATS_LAST +nCont

! TEMPORARY
    nBins = 20
!    nBins = pRegion%plagInput%nBins
! END TEMPORARY 
       
! ******************************************************************************
!   Initialize memory
! ******************************************************************************
        
    DO ifl = 1,nBFaces
      DO iBin = 1, nBins
        pStatsPlag(ifl)%nHits(iBin) = 0
      END DO ! iBin

      DO iBin = 1, nBins
      DO iVar = 1, nVars
        pStatsPlag(ifl)%vars(iVar,iBin) = 0.0_RFREAL
      END DO ! iVar
      END DO ! iBin                  
    END DO ! ifl                      

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_InitSurfStats






! *******************************************************************************
!
! Purpose: Nullify surface statistics.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_NullifySurfStats(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag,iPatch
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_DestroySurfStats',&
  'PLAG_ModSurfStats.F90')
    
    pGrid => pRegion%grid

! ******************************************************************************
!   Nullify pointer
! ******************************************************************************
     
! ==============================================================================
!     Patch data
! ==============================================================================

    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
        NULLIFY(pPatch%statsPlag)
      END IF ! pPatch%plotStatsFlag
    END DO ! iPatch
   
! ==============================================================================
!     Plane data
! ==============================================================================

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_NullifySurfStats






! *******************************************************************************
!
! Purpose: Open surface statistics file in ASCII format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   fileStatus          File status
!
! Output: 
!   fileExists          Flag indicating existence of file
!
! Notes:
!   1. File may not exist when it is opened in postprocessor before solver 
!      was run, so need to deal with this gracefully.
!
! ******************************************************************************

  SUBROUTINE PLAG_OpenSurfStatsASCII(pRegion,fileStatus,fileExists)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    INTEGER, INTENT(IN) :: fileStatus
    LOGICAL, INTENT(OUT), OPTIONAL :: fileExists
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: errorFlag,iFile
    TYPE(t_global), POINTER :: global 
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_OpenSurfStatsASCII',&
  'PLAG_ModSurfStats.F90')
         
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Opening ASCII PLAG surface statistics file...'
    END IF ! global%verbLevel     
     
    iFile = IF_PLAG_SURF_STATS
   
! ******************************************************************************
!   Build file name
! ******************************************************************************

    CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.plag_ssta', & 
                               pRegion%iRegionGlobal,global%currentTime, & 
                               iFileName)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN                                             
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
      WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                          global%currentTime  
    END IF ! global%verbLevel

! ******************************************************************************
!   Open file
! ******************************************************************************

    IF ( fileStatus == FILE_STATUS_OLD ) THEN 
      INQUIRE(FILE=iFileName,EXIST=fileExists)
    
      IF ( fileExists .EQV. .TRUE. ) THEN     
        OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD", &
             IOSTAT=errorFlag)
        global%error = errorFlag          
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
        END IF ! global%error
      END IF ! fileExists
    ELSE IF ( fileStatus == FILE_STATUS_UNKNOWN ) THEN 
      OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN", &
           IOSTAT=errorFlag)
      global%error = errorFlag          
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
      END IF ! global%error      
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)      
    END IF ! fileStatus       

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Opening ASCII PLAG surface statistics file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_OpenSurfStatsASCII






! *******************************************************************************
!
! Purpose: Open surface statistics file in binary format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   fileStatus          File status
!
! Output: 
!   fileExists          Flag indicating existence of file
!
! Notes: 
!   1. File may not exist when it is opened in postprocessor before solver 
!      was run, so need to deal with this gracefully.
!
! ******************************************************************************

  SUBROUTINE PLAG_OpenSurfStatsBinary(pRegion,fileStatus,fileExists)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    INTEGER, INTENT(IN) :: fileStatus
    LOGICAL, INTENT(OUT), OPTIONAL :: fileExists
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: errorFlag,iFile
    TYPE(t_global), POINTER :: global 
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_OpenSurfStatsBinary',&
  'PLAG_ModSurfStats.F90')
         
    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Opening binary PLAG surface statistics file...'
    END IF ! global%verbLevel     
     
    iFile = IF_PLAG_SURF_STATS     
     
! ******************************************************************************
!   Build file name
! ******************************************************************************

    CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.plag_sst', & 
                               pRegion%iRegionGlobal,global%currentTime, & 
                               iFileName)

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN                                             
      WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                       pRegion%iRegionGlobal
      WRITE(STDOUT,'(A,3X,A,1X,1PE11.5)') SOLVER_NAME,'Current time:', & 
                                          global%currentTime  
    END IF ! global%verbLevel

! ******************************************************************************
!   Open file
! ******************************************************************************

    IF ( fileStatus == FILE_STATUS_OLD ) THEN 
      INQUIRE(FILE=iFileName,EXIST=fileExists)
    
      IF ( fileExists .EQV. .TRUE. ) THEN     
        OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="OLD", &
             IOSTAT=errorFlag)
        global%error = errorFlag          
        IF ( global%error /= ERR_NONE ) THEN 
          CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
        END IF ! global%error
      END IF ! fileExists
    ELSE IF ( fileStatus == FILE_STATUS_UNKNOWN ) THEN 
      OPEN(iFile,FILE=iFileName,FORM="UNFORMATTED",STATUS="UNKNOWN", &
           IOSTAT=errorFlag)
      global%error = errorFlag          
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
      END IF ! global%error      
    ELSE 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)      
    END IF ! fileStatus         

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Opening binary PLAG surface statistics file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_OpenSurfStatsBinary
  
  
  
  
  
  
! *******************************************************************************
!
! Purpose: Read surface statistics in ASCII format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_ReadSurfStatsASCII(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,iFile,ifl,iPatch,loopCounter
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_ReadSurfStatsASCII',&
  'PLAG_ModSurfStats.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Reading ASCII PLAG surface statistics file...'
    END IF ! global%verbLevel
    
    pGrid => pRegion%grid

    iFile = IF_PLAG_SURF_STATS
         
! ******************************************************************************
!   Header and general information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    READ(iFile,'(A)') sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU surface statistics file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM

! ******************************************************************************
!   Read rest of file
! ******************************************************************************

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile,'(A)') sectionString

      SELECT CASE ( TRIM(sectionString) ) 

! ==============================================================================
!       Patch data
! ==============================================================================            

        CASE ( '# Patch' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Patch...'
          END IF ! global%verbLevel 

         DO iPatch = 1,pGrid%nPatches
           pPatch => pRegion%patches(iPatch)
           IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
             CALL PLAG_ReadSurfStatsKernelASCII(pRegion,pPatch%statsPlag)
           END IF ! pPatch%plotStatsFlag
         END DO ! iPatch

! ==============================================================================
!       Plane data
! ==============================================================================            

! ==============================================================================
!       End marker
! ==============================================================================
      
        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel           

          EXIT
      
! ==============================================================================
!       Invalid section string
! ==============================================================================
      
        CASE DEFAULT
          IF ( global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
          END IF ! verbosityLevel           

          CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)                               
      END SELECT ! TRIM
  
! ==============================================================================
!     Guard against infinite loop - unnecessary because of read errors?
! ==============================================================================
  
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter
    END DO ! <empty>       

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Reading ASCII PLAG surface statistics file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_ReadSurfStatsASCII






! *******************************************************************************
!
! Purpose: Read surface statistics in binary format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_ReadSurfStatsBinary(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,iFile,ifl,iPatch,loopCounter
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_ReadSurfStatsBinary',&
  'PLAG_ModSurfStats.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Reading ASCII PLAG surface statistics file...'
    END IF ! global%verbLevel
    
    pGrid => pRegion%grid

    iFile = IF_PLAG_SURF_STATS
         
! ******************************************************************************
!   Header and general information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    READ(iFile) sectionString
    IF ( TRIM(sectionString) /= '# ROCFLU surface statistics file' ) THEN 
      CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString) 
    END IF ! TRIM

! ******************************************************************************
!   Read rest of file
! ******************************************************************************

    loopCounter = 0

    DO ! set up infinite loop
      loopCounter = loopCounter + 1

      READ(iFile) sectionString

      SELECT CASE ( TRIM(sectionString) ) 

! ==============================================================================
!       Patch data
! ==============================================================================            

        CASE ( '# Patch' )
          IF ( global%myProcid == MASTERPROC .AND. &
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Patch...'
          END IF ! global%verbLevel 

         DO iPatch = 1,pGrid%nPatches
           pPatch => pRegion%patches(iPatch)
           IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
             CALL PLAG_ReadSurfStatsKernelBinary(pRegion,pPatch%statsPlag)
           END IF ! pPatch%plotStatsFlag
         END DO ! iPatch

! ==============================================================================
!       Plane data
! ==============================================================================            

! ==============================================================================
!       End marker
! ==============================================================================
      
        CASE ( '# End' ) 
          IF ( global%myProcid == MASTERPROC .AND. & 
               global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
          END IF ! global%verbLevel           

          EXIT
      
! ==============================================================================
!       Invalid section string
! ==============================================================================
      
        CASE DEFAULT
          IF ( global%verbLevel > VERBOSE_LOW ) THEN  
            WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,sectionString
          END IF ! verbosityLevel           

          CALL ErrorStop(global,ERR_INVALID_MARKER,__LINE__,sectionString)                               
      END SELECT ! TRIM
  
! ==============================================================================
!     Guard against infinite loop - unnecessary because of read errors?
! ==============================================================================
  
      IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! loopCounter
    END DO ! <empty>       

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Reading Binary PLAG surface statistics file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_ReadSurfStatsBinary






! *******************************************************************************
!
! Purpose: Kernel to read surface statistics in ASCII format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pStatsPlag          Pointer to statsPlag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_ReadSurfStatsKernelASCII(pRegion,pStatsPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    TYPE(t_surfstats_plag), DIMENSION(:), POINTER :: pStatsPlag
   
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,iBin,iFile,ifl,iPatch,iVar,nBFaces,nBins,nCont,nVars  

    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_ReadSurfStatsKernelASCII',&
  'PLAG_ModSurfStats.F90')

    iFile = IF_PLAG_SURF_STATS

    nBFaces = SIZE(pStatsPlag,DIM=1)

    nCont = pRegion%plagInput%nCont
    nVars = PLAG_SURF_STATS_LAST +nCont
        
! TEMPORARY
    nBins = 20
!    nBins = pRegion%plagInput%nBins
! END TEMPORARY

! ******************************************************************************
!   Read surface statistics to file
! ******************************************************************************
    
! ==============================================================================
!   Number of hits
! ==============================================================================

    READ(iFile,'(8(I10))') ((pStatsPlag(ifl)%nHits(iBin),iBin=1,nBins),&
                             ifl=1,nBFaces)        

! ==============================================================================
!   Surface variables
! ==============================================================================

    READ(iFile,'(5(E23.16))') (((pStatsPlag(ifl)%vars(iVar,iBin),         &
                                 iVar=1,nVars),iBin=1,nBins),ifl=1,nBFaces)
   
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_ReadSurfStatsKernelASCII






! *******************************************************************************
!
! Purpose: Kernel to read surface statistics in binary format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pStatsPlag          Pointer to statsPlag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_ReadSurfStatsKernelBinary(pRegion,pStatsPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    TYPE(t_surfstats_plag), DIMENSION(:), POINTER :: pStatsPlag
   
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,iBin,iFile,ifl,iPatch,iVar,nBFaces,nBins,nCont,nVars  

    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_ReadSurfStatsKernelBinary',&
  'PLAG_ModSurfStats.F90')

    iFile = IF_PLAG_SURF_STATS

    nBFaces = SIZE(pStatsPlag,DIM=1)

    nCont = pRegion%plagInput%nCont
    nVars = PLAG_SURF_STATS_LAST +nCont
        
! TEMPORARY
    nBins = 20
!    nBins = pRegion%plagInput%nBins
! END TEMPORARY

! ******************************************************************************
!   Read surface statistics to file
! ******************************************************************************
    
! ==============================================================================
!   Number of hits
! ==============================================================================

    READ(iFile) ((pStatsPlag(ifl)%nHits(iBin),iBin=1,nBins),ifl=1,nBFaces)        

! ==============================================================================
!   Surface variables
! ==============================================================================

    READ(iFile) (((pStatsPlag(ifl)%vars(iVar,iBin),         &
                   iVar=1,nVars),iBin=1,nBins),ifl=1,nBFaces)
   
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_ReadSurfStatsKernelBinary












! *******************************************************************************
!
! Purpose: Wrapper for reading surface statistics.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_ReadSurfStatsWrapper(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================

    LOGICAL :: fileExists    
    TYPE(t_global), POINTER :: global  
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_ReadSurfStatsWrapper',&
  'PLAG_ModSurfStats.F90')
    
    SELECT CASE(global%solutFormat)
      CASE( FORMAT_ASCII )
        CALL PLAG_OpenSurfStatsASCII(pRegion,FILE_STATUS_OLD,fileExists)
      
      IF ( fileExists .EQV. .TRUE. ) THEN
        CALL PLAG_ReadSurfStatsASCII(pRegion)
        CALL PLAG_CloseSurfStats(pRegion) 
      END IF ! fileExists

      CASE( FORMAT_BINARY )
        CALL PLAG_OpenSurfStatsBinary(pRegion,FILE_STATUS_OLD,fileExists)
      
         IF ( fileExists .EQV. .TRUE. ) THEN 
           CALL PLAG_ReadSurfStatsBinary(pRegion)
           CALL PLAG_CloseSurfStats(pRegion)      
         END IF ! fileExists
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! global%solutFormat   

! ******************************************************************************
!   Write warning if file is missing
! ******************************************************************************
    
    IF ( fileExists .EQV. .FALSE. ) THEN 
      global%warnCounter = global%warnCounter + 1
    
      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel > VERBOSE_NONE ) THEN 
        WRITE(STDOUT,'(A,2(1X,A))') SOLVER_NAME,'*** WARNING ***', &
                              'Plag surface statistics file missing, not read.'        
      END IF ! global%myProcid
    END IF ! fileExists     

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_ReadSurfStatsWrapper
  
  
  
  
  

! *******************************************************************************
!
! Purpose: Set bin index based on method selected..
!
! Description: None.
!
! Input:
!   global              Pointer to global
!   binMethod           Binning method selected
!   diamMicron          Particle diameter in microns
!   diamMin             Minimum diameter
!   diamMax             Maximum diameter
!
! Output: Bin index.
!
! Notes: None.
!
! ******************************************************************************

  INTEGER FUNCTION PLAG_SetBinIndex(global,binMethod,nBins,diamMicron,&
                                    diamBin,diamMin,diamMax)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_global), POINTER :: global   
    
    INTEGER, INTENT(IN) :: binMethod,nBins
    REAL(KIND=RFREAL), INTENT(IN) :: diamBin,diamMicron,diamMax,diamMin
    
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    CALL RegisterFunction(global,'PLAG_SetBinMethod',&
  'PLAG_ModSurfStats.F90')

! ******************************************************************************
!   Determine bin index
! ******************************************************************************

    SELECT CASE(binMethod)
      CASE( BIN_METHOD_LINEAR)
        PLAG_SetBinIndex = MIN(INT((diamMicron-diamMin)/diamBin)+1,nBins)

      CASE DEFAULT 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! binMethod

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END FUNCTION PLAG_SetBinIndex  
  
  
  
  

  
! *******************************************************************************
!
! Purpose: Write surface statistics in ASCII format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_WriteSurfStatsASCII(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,iFile,ifl,iPatch
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_WriteSurfStatsASCII',&
  'PLAG_ModSurfStats.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing ASCII PLAG surface statistics file...'
    END IF ! global%verbLevel
    
    pGrid => pRegion%grid

    iFile = IF_PLAG_SURF_STATS
         
! ******************************************************************************
!   Header and general information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    sectionString = '# ROCFLU surface statistics file'
    WRITE(iFile,'(A)') TRIM(sectionString)
                            
! ******************************************************************************
!   Write patch data
! ******************************************************************************

    sectionString = '# Patch'
    WRITE(iFile,'(A)') TRIM(sectionString)
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
        CALL PLAG_WriteSurfStatsKernelASCII(pRegion,pPatch%statsPlag)
      END IF ! pPatch%plotStatsFlag
    END DO ! iPatch
    
! ******************************************************************************
!   Write plane data
! ******************************************************************************

! ******************************************************************************
!   End marker
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile,'(A)') TRIM(sectionString)

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Writing ASCII PLAG surface statistics file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_WriteSurfStatsASCII





! *******************************************************************************
!
! Purpose: Write surface statistics in binary format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_WriteSurfStatsBinary(pRegion)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: sectionString
    INTEGER :: errorFlag,iFile,iPatch
    TYPE(t_global), POINTER :: global    
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_patch), POINTER :: pPatch
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_WriteSurfStatsBinary',&
  'PLAG_ModSurfStats.F90')

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN 
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Writing binary PLAG surface statistics file...'
    END IF ! global%verbLevel
    
    pGrid => pRegion%grid

    iFile = IF_PLAG_SURF_STATS

! ******************************************************************************
!   Header and general information
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Header information...'
    END IF ! global%verbLevel

    sectionString = '# ROCFLU surface statistics file'
    WRITE(iFile) sectionString
              
! ******************************************************************************
!   Write patch data
! ******************************************************************************

    sectionString = '# Patch'
    WRITE(iFile) sectionString
        
    DO iPatch = 1,pGrid%nPatches
      pPatch => pRegion%patches(iPatch)

      IF ( pPatch%plotStatsFlag .EQV. .TRUE. ) THEN
        CALL PLAG_WriteSurfStatsKernelBinary(pRegion,pPatch%statsPlag)
      END IF ! pPatch%plotStatsFlag
    END DO ! iPatch
    
! ******************************************************************************
!   Write plane data
! ******************************************************************************
    
          

! ******************************************************************************
!   End marker
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_LOW ) THEN  
      WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
    END IF ! global%verbLevel

    sectionString = '# End'
    WRITE(iFile) sectionString

! ******************************************************************************
!   End  
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN   
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                               'Writing binary PLAG surface statistics file done.'
    END IF ! global%verbLevel

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_WriteSurfStatsBinary






! *******************************************************************************
!
! Purpose: Kernel to write surface statistics in ASCII format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pStatsPlag          Pointer to statsPlag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_WriteSurfStatsKernelASCII(pRegion,pStatsPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    TYPE(t_surfstats_plag), DIMENSION(:), POINTER :: pStatsPlag
   
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName,sectionString
    INTEGER :: errorFlag,iBin,iFile,ifl,iPatch,iVar,nBFaces,nBins,nCont,nVars  

    TYPE(t_global), POINTER :: global     

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_WriteSurfStatsKernelASCII',&
  'PLAG_ModSurfStats.F90')

    iFile = IF_PLAG_SURF_STATS

    nBFaces = SIZE(pStatsPlag,DIM=1)

    nCont = pRegion%plagInput%nCont
    nVars = PLAG_SURF_STATS_LAST +nCont
      
! TEMPORARY
    nBins = 20
!    nBins = pRegion%plagInput%nBins
! END TEMPORARY

! ******************************************************************************
!   Write surface statistics to file
! ******************************************************************************
    
! ==============================================================================
!   Number of hits
! ==============================================================================

    WRITE(iFile,'(8(I10))') ((pStatsPlag(ifl)%nHits(iBin),iBin=1,nBins),&
                              ifl=1,nBFaces)        

! ==============================================================================
!   Surface variables
! ==============================================================================

    WRITE(iFile,'(5(E23.16))') (((pStatsPlag(ifl)%vars(iVar,iBin),         &
                                  iVar=1,nVars),iBin=1,nBins),ifl=1,nBFaces)
   
! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_WriteSurfStatsKernelASCII






! *******************************************************************************
!
! Purpose: Kernel to write surface statistics in binary format.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pStatsPlag          Pointer to statsPlag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_WriteSurfStatsKernelBinary(pRegion,pStatsPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================

    TYPE(t_region), POINTER :: pRegion   

    TYPE(t_surfstats_plag), DIMENSION(:), POINTER :: pStatsPlag
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    CHARACTER(CHRLEN) :: iFileName
    INTEGER :: errorFlag,iBin,iFile,ifl,iPatch,iVar,nBFaces,nBins,nCont,nVars

    TYPE(t_global), POINTER :: global  
 
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'PLAG_WriteSurfStatsKernelBinary',&
  'PLAG_ModSurfStats.F90')

    iFile = IF_PLAG_SURF_STATS 

    nBFaces = SIZE(pStatsPlag,DIM=1)

    nCont = pRegion%plagInput%nCont
    nVars = PLAG_SURF_STATS_LAST +nCont
          
! TEMPORARY
    nBins = 20
!    nBins = pRegion%plagInput%nBins
! END TEMPORARY

! ******************************************************************************
!   Write surface statistics to file
! ******************************************************************************
    
! ==============================================================================
!   Number of hits
! ==============================================================================

    WRITE(iFile) ((pStatsPlag(ifl)%nHits(iBin),iBin=1,nBins),ifl=1,nBFaces)       
    
! ==============================================================================
!   Surface variables
! ==============================================================================

    WRITE(iFile) (((pStatsPlag(ifl)%vars(iVar,iBin),         &
                    iVar=1,nVars),iBin=1,nBins),ifl=1,nBFaces)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_WriteSurfStatsKernelBinary





! *******************************************************************************
!
! Purpose: Wrapper for writing surface statistics.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_WriteSurfStatsWrapper(pRegion)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
    TYPE(t_region), POINTER :: pRegion
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    TYPE(t_global), POINTER :: global  
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_WriteSurfStatsWrapper',&
  'PLAG_ModSurfStats.F90')
    
    SELECT CASE( global%solutFormat ) 
      CASE( FORMAT_ASCII )
        CALL PLAG_OpenSurfStatsASCII(pRegion,FILE_STATUS_UNKNOWN)
        CALL PLAG_WriteSurfStatsASCII(pRegion)
        CALL PLAG_CloseSurfStats(pRegion)            
      CASE( FORMAT_BINARY )
        CALL PLAG_OpenSurfStatsBinary(pRegion,FILE_STATUS_UNKNOWN)
        CALL PLAG_WriteSurfStatsBinary(pRegion)
        CALL PLAG_CloseSurfStats(pRegion)      
      CASE DEFAULT 
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! global%solutFormat   
    
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_WriteSurfStatsWrapper






END MODULE PLAG_ModSurfStats

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModSurfStats.F90,v $
! Revision 1.13  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.12  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.11  2006/05/12 22:50:32  fnajjar
! Fixed binary write of strings by remove TRIM
!
! Revision 1.10  2006/05/11 18:17:44  fnajjar
! Fixed inconsistent naming for register functions
!
! Revision 1.9  2006/05/09 14:38:10  fnajjar
! Bug for a formatted read from a binary file
!
! Revision 1.8  2006/05/04 13:54:58  fnajjar
! Added and activated binary read for patch statistics file
!
! Revision 1.7  2006/05/02 17:46:23  fnajjar
! Allowed surface statistics to be gathered on active patches
!
! Revision 1.6  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.5  2005/12/30 16:27:18  fnajjar
! Changed diamMaxMicron to have size of 10 in bins
!
! Revision 1.4  2005/01/06 22:26:43  fnajjar
! Implemented clean version for PLAG_SetBinIndex
!
! Revision 1.3  2005/01/06 16:00:31  fnajjar
! Reformulated IDINT to INT with KIND in PLAG_SetBinIndex
!
! Revision 1.2  2005/01/05 19:26:53  fnajjar
! Fixed definition of diamBin to maintain bounds
!
! Revision 1.1  2004/12/21 15:07:02  fnajjar
! Initial import of surface statistics module
!
! ******************************************************************************



























