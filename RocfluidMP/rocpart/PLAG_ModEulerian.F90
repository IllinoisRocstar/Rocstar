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
! Purpose: Collection of routines for Eulerian-based particle infrastructure.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: PLAG_ModEulerian.F90,v 1.8 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_ModEulerian

  USE ModParameters
  USE ModDataTypes  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY: t_plag 
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModMPI
  USE PLAG_ModParameters

#ifdef RFLO
#include "Indexing.h"  
  USE ModInterfaces, ONLY: RFLO_GetCellOffSet, &
                           RFLO_GetDimensDummy
#endif
 
  IMPLICIT NONE
  
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Private data
! ==============================================================================

  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: PLAG_ModEulerian.F90,v $ $Revision: 1.8 $'        

! ==============================================================================
! Public functions
! ==============================================================================
  
  PUBLIC :: PLAG_CalcEulerianField,    &
            PLAG_CreateEulerianField,  & 
            PLAG_DestroyEulerianField, & 
            PLAG_InitEulerianField

! ==============================================================================
! Private functions
! ==============================================================================

!  PRIVATE :: 

! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
  
  
  
  
  
  
  

  





! *******************************************************************************
!
! Purpose: Compute Eulerian-based field from Lagrangian datastructure.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!
! Output: None.
!
! Notes: Routine relevant to RFLO.
!
! ******************************************************************************

  SUBROUTINE PLAG_CalcEulerianField(pRegion)

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

#ifdef RFLO
    INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend  
    INTEGER :: iLev,iCellOffset, ijCellOffset   
#endif

    INTEGER :: errorFlag,ibc,ic,iCont,iec,iPcl,iVar,nCont,nPcls,nPlagCell,nVars
    INTEGER, POINTER, DIMENSION(:)   :: cvMass    

    REAL(RFREAL) :: diamMicron,massL,massSqrL,nPlagCellInv
    REAL(RFREAL), POINTER, DIMENSION(:,:) :: cv,dv,ev

    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'PLAG_CalcEulerianField',&
  'PLAG_ModEulerian.F90')

#ifdef RFLO
    iLev  = pRegion%currLevel
    pGrid => pRegion%levels(iLev)%grid
    pPlag => pRegion%levels(iLev)%plag

    CALL RFLO_GetDimensDummy( pRegion,iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( pRegion,iLev,iCellOffset,ijCellOffset )
    ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCellOffset,ijCellOffset)
    iec = IndIJK(idcend,jdcend,kdcend,iCellOffset,ijCellOffset)
#endif

#ifdef RFLU
    pGrid => pRegion%grid
    pPlag => pRegion%plag
    
    ibc = 1
    iec = pGrid%nCells
#endif
   
    nCont = pRegion%plagInput%nCont
    nPcls = pPlag%nPcls
    nVars = pPlag%nEv

    cvMass => pPlag%cvPlagMass
    cv => pPlag%cv
    dv => pPlag%dv
    ev => pPlag%ev

    IF ( pPlag%nPcls > 0 ) THEN

! ******************************************************************************
!     Reinitialize Eulerian field to compute instantaneous variables
! ******************************************************************************

      DO ic = ibc,iec
      DO iVar = 1, nVars
        ev(iVar,ic) = 0.0_RFREAL
      END DO ! iVar
      END DO ! ic

! ******************************************************************************
!     Compute instantaneous Eulerian field from Lagrangian datastructure
! ******************************************************************************

      DO iPcl = 1,nPcls
        ic = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)

        diamMicron = dv(DV_PLAG_DIAM,iPcl)*1.0E+06_RFREAL      
        massL = SUM(cv(cvMass(:),iPcl))

        ev(EV_PLAG_DIAM3,ic) = ev(EV_PLAG_DIAM3,ic) +diamMicron**3.0_RFREAL
        ev(EV_PLAG_DIAM4,ic) = ev(EV_PLAG_DIAM4,ic) +diamMicron**4.0_RFREAL

        ev(EV_PLAG_NUMDENS,ic) = ev(EV_PLAG_NUMDENS,ic) +1.0_RFREAL
 
        ev(EV_PLAG_UVEL,ic) = ev(EV_PLAG_UVEL,ic) +dv(DV_PLAG_UVEL,iPcl)
        ev(EV_PLAG_VVEL,ic) = ev(EV_PLAG_VVEL,ic) +dv(DV_PLAG_VVEL,iPcl)
        ev(EV_PLAG_WVEL,ic) = ev(EV_PLAG_WVEL,ic) +dv(DV_PLAG_WVEL,iPcl)

        ev(EV_PLAG_MASS,ic) = ev(EV_PLAG_MASS,ic) +massL

        DO iCont = 1, nCont
          ev(EV_PLAG_LAST+iCont,ic) = ev(EV_PLAG_LAST+iCont,ic) &
                                    + cv(cvMass(iCont),iPcl)
        END DO ! iCont

      END DO ! iPcl

! ******************************************************************************
!     Scale field by cell-based number of particles
!       Do not apply scaling to number density
! ******************************************************************************

      DO ic = ibc,iec
        nPlagCell = ev(EV_PLAG_NUMDENS,ic)
   
        IF ( nPlagCell > 0 ) THEN
          nPlagCellInv = 1.0/REAL(nPlagCell,KIND=RFREAL)
      
          DO iVar = 1, EV_PLAG_NUMDENS-1
            ev(iVar,ic) = ev(iVar,ic) *nPlagCellInv
          END DO ! iVar
      
          DO iVar = EV_PLAG_NUMDENS+1, nVars
            ev(iVar,ic) = ev(iVar,ic) *nPlagCellInv
          END DO ! iVar

         ENDIF ! nPlagCell 
       END DO ! ic  

    END IF ! pPlag%nPcls 

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_CalcEulerianField

  
  
  
  
  

! *******************************************************************************
!
! Purpose: Create eulerian field.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPlag               Pointer to plag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_CreateEulerianField(pRegion,pPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
#ifdef RFLO
    TYPE(t_region) :: pRegion 
#endif

#ifdef RFLU
    TYPE(t_region), POINTER :: pRegion
#endif

    TYPE(t_plag), POINTER :: pPlag

! ==============================================================================
!   Locals
! ==============================================================================

#ifdef RFLO
    INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend  
    INTEGER :: iLev,iCellOffset, ijCellOffset   
#endif

    INTEGER :: errorFlag,ibc,iec,nEv
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    
! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************
    
    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_CreateEulerian',&
  'PLAG_ModEulerian.F90')

    nEv = pPlag%nEv        

#ifdef RFLO
    iLev  = pRegion%currLevel
    pGrid => pRegion%levels(iLev)%grid

    CALL RFLO_GetDimensDummy( pRegion,iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( pRegion,iLev,iCellOffset,ijCellOffset )
    ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCellOffset,ijCellOffset)
    iec = IndIJK(idcend,jdcend,kdcend,iCellOffset,ijCellOffset)
#endif

#ifdef RFLU
    pGrid => pRegion%grid
    ibc = 1
    iec = pGrid%nCells
#endif
  
! ******************************************************************************
!   Allocate memory
! ******************************************************************************

    ALLOCATE(pPlag%ev(nEv,ibc:iec),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'pPlag%ev')
    END IF ! global%error                       

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_CreateEulerianField






! *******************************************************************************
!
! Purpose: Destroy Eulerian field.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPlag               Pointer to plag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_DestroyEulerianField(pRegion,pPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
    
#ifdef RFLO
    TYPE(t_region) :: pRegion 
#endif

#ifdef RFLU
    TYPE(t_region), POINTER :: pRegion
#endif

    TYPE(t_plag), POINTER :: pPlag
    
! ==============================================================================
!   Locals
! ==============================================================================
    
    INTEGER :: errorFlag
    TYPE(t_global), POINTER :: global    

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
    
    CALL RegisterFunction(global,'PLAG_DestroyEulerian',&
  'PLAG_ModEulerian.F90')

! ******************************************************************************
!   Deallocate memory
! ******************************************************************************
        
    DEALLOCATE(pPlag%ev,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'pPlag%ev')
    END IF ! global%error

! ******************************************************************************
!   Nullify memory
! ******************************************************************************
 
    NULLIFY(pPlag%ev)

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_DestroyEulerianField

  



! *******************************************************************************
!
! Purpose: Initialize Eulerian field.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   pPlag               Pointer to plag
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE PLAG_InitEulerianField(pRegion,pPlag)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================
     
#ifdef RFLO
    TYPE(t_region) :: pRegion
#endif

#ifdef RFLU
    TYPE(t_region), POINTER :: pRegion
#endif

    TYPE(t_plag), POINTER :: pPlag
  
! ==============================================================================
!   Locals
! ==============================================================================

#ifdef RFLO
    INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend  
    INTEGER :: iLev,iCellOffset, ijCellOffset   
#endif

    INTEGER :: errorFlag,ibc,iec,iCell,iVar,nVars
    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'PLAG_InitEulerianField',&
  'PLAG_ModEulerian.F90')

    nVars = pPlag%nEv       

#ifdef RFLO
    iLev  = pRegion%currLevel
    pGrid => pRegion%levels(iLev)%grid

    CALL RFLO_GetDimensDummy( pRegion,iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( pRegion,iLev,iCellOffset,ijCellOffset )
    ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCellOffset,ijCellOffset)
    iec = IndIJK(idcend,jdcend,kdcend,iCellOffset,ijCellOffset)
#endif

#ifdef RFLU
    pGrid => pRegion%grid
    ibc = 1
    iec = pGrid%nCells
#endif

! ******************************************************************************
!   Initialize memory
! ******************************************************************************
        
    DO iCell = ibc,iec
      DO iVar = 1, nVars
        pPlag%ev(iVar,iCell)  = 0.0_RFREAL
      END DO ! iVar                  
    END DO ! iCell                      

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_InitEulerianField







#ifdef RFLU

! *******************************************************************************
!
! Purpose: Compute Eulerian-based field from Lagrangian datastructure.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region
!   ev                  array containing Eulerian field
!
! Output: None.
!
! Notes: Routine relevant to write TECPLOT files in RFLU.
!
! ******************************************************************************

  SUBROUTINE PLAG_RFLU_CalcEulerianField(pRegion,ev)

    IMPLICIT NONE
        
! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************
    
! ==============================================================================
!   Arguments
! ==============================================================================

    REAL(RFREAL), DIMENSION(:,:) :: ev     
    TYPE(t_region), POINTER :: pRegion
  
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: errorFlag,ibc,icg,iCont,iec,iPcl,iVar,nCells1D,nCont,nPcls, & 
               nPlagCell,nVars,nVarsPlot
    INTEGER, POINTER, DIMENSION(:)   :: cvMass    

    REAL(RFREAL) :: diamMicron,massL,massSqrL,nPlagCellInv
    REAL(RFREAL) :: massFracL,densL,densMixt,spLoadL,volL,volMixt
    REAL(RFREAL), POINTER, DIMENSION(:,:) :: cv,dv

    TYPE(t_global), POINTER :: global
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_plag), POINTER :: pPlag

! ******************************************************************************
!   Start, set pointers and variables
! ******************************************************************************

    global => pRegion%global
 
    CALL RegisterFunction(global,'PLAG_RFLU_CalcEulerianField',&
  'PLAG_ModEulerian.F90')

    pGrid => pRegion%grid
    pPlag => pRegion%plag

    ibc = 1
    iec = pGrid%nCells
   
    nCont = pRegion%plagInput%nCont
    nPcls = pPlag%nPcls
    nVars = pPlag%nEv

    cvMass => pPlag%cvPlagMass
    cv => pPlag%cv
    dv => pPlag%dv

    IF ( pPlag%nPcls > 0 ) THEN

! ******************************************************************************
!     Check if size of array sent is aligned with the expected size
! ******************************************************************************

! CHECK
    nVarsPlot = pPlag%nEv -nCont
    
    IF ( nVarsPlot /= SIZE(ev,1) ) THEN 
      WRITE(*,*) 'nVarsPlot   = ',nVars
      WRITE(*,*) 'SIZE(ev,1) = ',SIZE(ev,1)
      CALL ErrorStop(global,ERR_DATADIM_MISMATCH,__LINE__)
    END IF ! nVars 

    IF ( iec /= SIZE(ev,2) ) THEN 
      WRITE(*,*) 'iec        = ',iec
      WRITE(*,*) 'SIZE(ev,2) = ',SIZE(ev,2)
      CALL ErrorStop(global,ERR_DATADIM_MISMATCH,__LINE__)
    END IF ! nVars 
        
! END CHECK

! ******************************************************************************
!     Reinitialize Eulerian field to compute instantaneous variables
! ******************************************************************************

      DO icg = ibc,iec
      DO iVar = 1, nVarsPlot
        ev(iVar,icg) = 0.0_RFREAL
      END DO ! iVar
      END DO ! icg

! ******************************************************************************
!     Compute instantaneous Eulerian field from Lagrangian datastructure
! ******************************************************************************

      DO iPcl = 1,nPcls
        icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcl)

        diamMicron = dv(DV_PLAG_DIAM,iPcl)*1.0E+06_RFREAL      
        massL      = SUM(cv(cvMass(:),iPcl))
        spLoadL    = pPlag%arv(ARV_PLAG_SPLOAD,iPcl)
        densL      = dv(DV_PLAG_DENS,iPcl)
        
        volL       = dv(DV_PLAG_VOLU,iPcl)
        massFracL  = densL*volL*spLoadL

        ev(EV_PLAG_DIAM3,icg) = ev(EV_PLAG_DIAM3,icg) +diamMicron**3.0_RFREAL
        ev(EV_PLAG_DIAM4,icg) = ev(EV_PLAG_DIAM4,icg) +diamMicron**4.0_RFREAL

        ev(EV_PLAG_NUMDENS,icg) = ev(EV_PLAG_NUMDENS,icg) +1.0_RFREAL
 
        ev(EV_PLAG_UVEL,icg) = ev(EV_PLAG_UVEL,icg) +dv(DV_PLAG_UVEL,iPcl)
        ev(EV_PLAG_VVEL,icg) = ev(EV_PLAG_VVEL,icg) +dv(DV_PLAG_VVEL,iPcl)
        ev(EV_PLAG_WVEL,icg) = ev(EV_PLAG_WVEL,icg) +dv(DV_PLAG_WVEL,iPcl)
        
        ev(EV_PLAG_TEMP,icg) = ev(EV_PLAG_TEMP,icg) +dv(DV_PLAG_TEMP,iPcl)

        ev(EV_PLAG_MASS,icg) = ev(EV_PLAG_MASS,icg) +massFracL

! TEMPORARY
!    Not used for RFLU in tecplot files
!        DO iCont = 1, nCont
!          ev(EV_PLAG_LAST+iCont,icg) = ev(EV_PLAG_LAST+iCont,icg) &
!                                    + cv(cvMass(iCont),iPcl)
!        END DO ! iCont
! END TEMPORARY

      END DO ! iPcl

! ******************************************************************************
!     Scale field by cell-based number of particles
!       Do not apply scaling to number density and mass fraction
! ******************************************************************************

      DO icg = ibc,iec
        nPlagCell = ev(EV_PLAG_NUMDENS,icg)
   
        IF ( nPlagCell > 0 ) THEN
          nPlagCellInv = 1.0/REAL(nPlagCell,KIND=RFREAL)
      
          DO iVar = 1, EV_PLAG_NUMDENS-1
            ev(iVar,icg) = ev(iVar,icg) *nPlagCellInv
          END DO ! iVar

          DO iVar = EV_PLAG_NUMDENS+1, EV_PLAG_TEMP     
! TEMPORARY
!          DO iVar = EV_PLAG_NUMDENS+1, nVars
! END TEMPORARY

            ev(iVar,icg) = ev(iVar,icg) *nPlagCellInv
          END DO ! iVar

! Compute Mass Fraction
          
          densMixt = pRegion%mixt%cv(CV_MIXT_DENS,icg)
          volMixt  = pRegion%grid%vol(icg)
          ev(EV_PLAG_MASS,icg) = ev(EV_PLAG_MASS,icg)/ &
                                (ev(EV_PLAG_MASS,icg) +densMixt *volMixt)
         ENDIF ! nPlagCell 
       END DO ! icg  
    END IF ! pPlag%nPcls 

! ******************************************************************************
!   End  
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE PLAG_RFLU_CalcEulerianField
#endif





END MODULE PLAG_ModEulerian

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModEulerian.F90,v $
! Revision 1.8  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.5  2006/02/16 14:51:58  fnajjar
! Bug fix for missing ev pointer array in PLAG_CalcEulerian
!
! Revision 1.4  2006/01/31 17:20:21  fnajjar
! Bug fix to embed PLAG_RFLU_CalcEulerian in ifdef RFLU
!
! Revision 1.3  2005/11/30 22:20:21  fnajjar
! Added PLAG_RFLU_CalcEulerianField
!
! Revision 1.2  2005/02/16 14:46:48  fnajjar
! Bug fix for w-velocity component of ev and cosmetics cleanup
!
! Revision 1.1  2005/01/08 20:44:32  fnajjar
! Initial import for PLAG statistics
!
! ******************************************************************************











