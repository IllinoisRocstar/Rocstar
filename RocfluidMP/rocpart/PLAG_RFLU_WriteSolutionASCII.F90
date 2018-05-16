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
! Purpose: Write flow file for particles in ASCII ROCFLU format.
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
! $Id: PLAG_RFLU_WriteSolutionASCII.F90,v 1.8 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_WriteSolutionASCII(pRegion)

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

  CHARACTER(CHRLEN) :: iFileName,sectionString,RCSIdentString
  INTEGER :: errorFlag,iCont,iFile,ifl,iMass,iPatch,iVar,j,nCont,nVars 
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv
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

  RCSIdentString = '$RCSfile: PLAG_RFLU_WriteSolutionASCII.F90,v $ $Revision: 1.8 $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_WriteSolutionASCII',&
  'PLAG_RFLU_WriteSolutionASCII.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII particle file...'
  END IF ! global%verbLevel

  CALL BuildFileNameUnsteady(global,FILEDEST_OUTDIR,'.plag_sola', & 
                             pRegion%iRegionGlobal,global%currentTime, & 
                             iFileName)

  iFile = IF_SOLUT
  OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="UNKNOWN",IOSTAT=errorFlag)
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

  sectionString = '# ROCFLU particle file'
  WRITE(iFile,'(A)') sectionString
  
  sectionString = '# Precision and range'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(2(I8))') PRECISION(1.0_RFREAL),RANGE(1.0_RFREAL)
  
  sectionString = '# Physical time'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(E23.16)') global%currentTime 

! ==============================================================================
! Dimensions
! ==============================================================================
  
  pGrid => pRegion%grid  
  pPlag => pRegion%plag
   
  nCont = pRegion%plagInput%nCont   
   
  nVars = 13 ! Hard-coded for now
   
  sectionString = '# Dimensions'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(2(I8))') pPlag%nPcls,nVars
  
! ==============================================================================
! State vector
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Conserved variables...'
  END IF ! global%verbLevel

  sectionString = '# Particle x-momentum'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(5(E23.16))') (pPlag%cv(CV_PLAG_XMOM,j),j=1,pPlag%nPcls)

  sectionString = '# Particle y-momentum'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(5(E23.16))') (pPlag%cv(CV_PLAG_YMOM,j),j=1,pPlag%nPcls)

  sectionString = '# Particle z-momentum'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(5(E23.16))') (pPlag%cv(CV_PLAG_ZMOM,j),j=1,pPlag%nPcls)
 
  sectionString = '# Particle energy'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(5(E23.16))') (pPlag%cv(CV_PLAG_ENER,j),j=1,pPlag%nPcls) 

  sectionString = '# Particle x-location'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(5(E23.16))') (pPlag%cv(CV_PLAG_XPOS,j),j=1,pPlag%nPcls)
  
  sectionString = '# Particle y-location'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(5(E23.16))') (pPlag%cv(CV_PLAG_YPOS,j),j=1,pPlag%nPcls)
  
  sectionString = '# Particle z-location'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(5(E23.16))') (pPlag%cv(CV_PLAG_ZPOS,j),j=1,pPlag%nPcls)    
  
  sectionString = '# Particle vapor energy'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(5(E23.16))') (pPlag%cv(CV_PLAG_ENERVAPOR,j),j=1,pPlag%nPcls)    
  
  DO iCont = 1,pRegion%plagInput%nCont
    iMass = pPlag%cvPlagMass(iCont)

    sectionString = '# Particle mass'       
    WRITE(iFile,'(A)') sectionString
    WRITE(iFile,'(5(E23.16))') (pPlag%cv(iMass,j),j=1,pPlag%nPcls)  
  END DO ! iCont
 
! ==============================================================================
! Additional variables
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Additional variables...'
  END IF ! global%verbLevel

  sectionString = '# Particle superloading'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(5(E23.16))') (pPlag%arv(ARV_PLAG_SPLOAD,j),j=1,pPlag%nPcls)
 
  sectionString = '# Particle initial identifier'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(5(I8))') (pPlag%aiv(AIV_PLAG_PIDINI,j),j=1,pPlag%nPcls)
 
  sectionString = '# Particle initial region'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(5(I8))') (pPlag%aiv(AIV_PLAG_REGINI,j),j=1,pPlag%nPcls) 

  sectionString = '# Particle cell'
  WRITE(iFile,'(A)') sectionString
  WRITE(iFile,'(5(I8))') (pPlag%aiv(AIV_PLAG_ICELLS,j),j=1,pPlag%nPcls) 
  
! ==============================================================================
! Patch data
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Patch data...'
  END IF ! global%verbLevel

  sectionString = '# Patch data'
  WRITE(iFile,'(A)') sectionString      
    
  DO iPatch = 1,pGrid%nPatches
    pPatch => pRegion%patches(iPatch)
    
    IF ( (pPatch%bcType >= BC_INJECTION .AND. pPatch%bcType <= BC_INJECTION + BC_RANGE) .OR. &
         (pPatch%bcType >= BC_INFLOW    .AND. pPatch%bcType <= BC_INFLOW    + BC_RANGE)      ) THEN 
      pTilePlag   => pPatch%tilePlag      
      
      DO ifl = 1,pPatch%nBFaces
        WRITE(iFile,'(2(E23.16))') pTilePlag%cv(CV_TILE_MOMNRM,ifl), & 
                                   pTilePlag%cv(CV_TILE_ENER  ,ifl)
      END DO ! ifl      
      
      DO iCont = 1,nCont
        iMass = pTilePlag%cvTileMass(iCont)
        WRITE(iFile,'(5(E23.16))') (pTilePlag%cv(iMass,ifl), &
                                    ifl=1,pPatch%nBFaces)
      END DO ! iCont
      
      DO ifl = 1,pPatch%nBFaces
        WRITE(iFile,'(3(E23.16))') pTilePlag%dv(DV_TILE_COUNTDOWN,ifl), & 
                                   pTilePlag%dv(DV_TILE_DIAM    ,ifl), & 
                                   pTilePlag%dv(DV_TILE_SPLOAD  ,ifl)
      END DO ! ifl                               
    END IF ! pPatch%bcType
  END DO ! iPatch

! ==============================================================================
! End marker
! ==============================================================================

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN  
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'End marker...'
  END IF ! global%verbLevel

  sectionString = '# End'
  WRITE(iFile,'(A)') sectionString  

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
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Writing ASCII particle file done.'
  END IF ! global%verbLevel

  CALL DeregisterFunction(global)
    
! ******************************************************************************
! End
! ******************************************************************************
 
END SUBROUTINE PLAG_RFLU_WriteSolutionASCII


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_WriteSolutionASCII.F90,v $
! Revision 1.8  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2006/09/18 20:37:02  fnajjar
! Activated tile datastructure for inflow bc
!
! Revision 1.5  2005/01/21 17:23:10  fnajjar
! Included vapor energy in IO capability
!
! Revision 1.4  2004/06/16 23:03:09  fnajjar
! Renamed DV_TILE_TIMEFCTR to DV_TILE_COUNTDOWN for CRE kernel
!
! Revision 1.3  2004/06/16 20:01:21  haselbac
! Added use of ModBuildFileNames, cosmetics
!                                   
! Revision 1.2  2004/03/05 23:19:33  haselbac                                  
! Added writing of additional patch data and modified writing of constituents  
!
! Revision 1.1  2004/02/26 21:00:52  haselbac                                  
! Initial revision                                                             
!
! ******************************************************************************







