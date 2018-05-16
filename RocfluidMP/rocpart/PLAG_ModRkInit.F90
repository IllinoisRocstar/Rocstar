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
!*******************************************************************************
!
! Purpose: Suite of routines to initialize Runge-Kutta schemes.
!
! Description: None.
!
! Notes: None.
!
!*******************************************************************************
!
! $Id: PLAG_ModRkInit.F90,v 1.8 2009/03/02 00:19:36 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
!*******************************************************************************

MODULE PLAG_ModRkInit
  
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModDataTypes
  USE ModGlobal,     ONLY: t_global
  USE ModGrid,       ONLY: t_grid
  USE ModPartLag,    ONLY: t_plag, t_plag_input, t_tile_plag
  USE ModBndPatch,   ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModError

  USE PLAG_ModParameters
  USE INRT_ModParameters

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: PLAG_InjcTileRkInit, &
            PLAG_RkInitPrimary
        
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: PLAG_ModRkInit.F90,v $ $Revision: 1.8 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS




!******************************************************************************
!
! Purpose: Initialize Runge-Kutta scheme for Integer variables.
!
! Description: None.
!
! Input: 
!   region      Region data
!   iStage      Runge-Kutta stage
!   icBeg       Beginning index for cell update
!   icEnd       Ending index for cell update
!   ivBeg       Beginning index for variable update
!   ivEnd       Ending index for variable update
!   aiv         Integer variables
!   aivOld      Old integer variables
!
! Output: 
!   ivOld       Old integer variables
!
! Notes: None.
!
!******************************************************************************

SUBROUTINE PLAG_RkInitGenericInt(region,iStage,icBeg,icEnd,ivBeg,ivEnd,&
                                 aiv,aivOld)

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: icBeg,icEnd,iStage,ivBeg,ivEnd
  INTEGER, DIMENSION(:,:), POINTER :: aiv,aivOld
  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: ic,iv
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ModRkInit.F90,v $ $Revision: 1.8 $'

  global => region%global

  CALL RegisterFunction(global,'PLAG_RkInitGenericInt',&
  'PLAG_ModRkInit.F90')

! *****************************************************************************
! Initialize Runge-Kutta scheme
! *****************************************************************************

  SELECT CASE ( global%rkScheme ) 
    CASE ( RK_SCHEME_4_CLASSICAL ) 
      IF ( iStage == 1 ) THEN
        DO ic = icBeg,icEnd
          DO iv = ivBeg,ivEnd
            aivOld(iv,ic) = aiv(iv,ic)
          END DO ! iv
        END DO ! ic
      END IF ! iStage
    CASE ( RK_SCHEME_3_WRAY ) 
      DO ic = icBeg,icEnd
        DO iv = ivBeg,ivEnd
          aivOld(iv,ic) = aiv(iv,ic)
        END DO ! iv
      END DO ! ic   
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%rkScheme

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RkInitGenericInt





!******************************************************************************
!
! Purpose: Initialize Runge-Kutta scheme for Real variables.
!
! Description: None.
!
! Input: 
!   region      Region data
!   iStage      Runge-Kutta stage
!   icBeg       Beginning index for cell update
!   icEnd       Ending index for cell update
!   ivBeg       Beginning index for variable update
!   ivEnd       Ending index for variable update
!   rv          Real variables
!   rvOld       Old real variables
!
! Output: 
!   rvOld       Old real variables
!
! Notes: None.
!
!******************************************************************************

SUBROUTINE PLAG_RkInitGenericReal(region,iStage,icBeg,icEnd,ivBeg,ivEnd,&
                                  rv,rvOld)

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: icBeg,icEnd,iStage,ivBeg,ivEnd
  REAL(RFREAL), DIMENSION(:,:), POINTER :: rv,rvOld
  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: ic,iv
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ModRkInit.F90,v $ $Revision: 1.8 $'

  global => region%global

  CALL RegisterFunction(global,'PLAG_RkInitGenericReal',&
  'PLAG_ModRkInit.F90')

! *****************************************************************************
! Initialize Runge-Kutta scheme
! *****************************************************************************

  SELECT CASE ( global%rkScheme ) 
    CASE ( RK_SCHEME_4_CLASSICAL ) 
      IF ( iStage == 1 ) THEN
        DO ic = icBeg,icEnd
          DO iv = ivBeg,ivEnd
            rvOld(iv,ic) = rv(iv,ic)
          END DO ! iv
        END DO ! ic
      END IF ! iStage
    CASE ( RK_SCHEME_3_WRAY ) 
      DO ic = icBeg,icEnd
        DO iv = ivBeg,ivEnd
          rvOld(iv,ic) = rv(iv,ic)
        END DO ! iv
      END DO ! ic   
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! global%rkScheme

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RkInitGenericReal





!******************************************************************************
!
! Purpose: Driver that initializes primary variables.
!
! Description: none.
!
! Input: istage = RK stage
!        region = data of current region.
!
! Output: region%levels%plag%cvOld 
!         region%levels%plag%aivOld
!         region%levels%plag%arvOld
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE PLAG_RkInitPrimary( region, iStage )

! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: iStage
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

#ifdef RFLO
  INTEGER :: iLev
#endif
  INTEGER :: nAiv,nArv,nCv,nPcls
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv, pAivOld

  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv, pArvOld, pCv, pCvOld
  
  TYPE(t_plag),   POINTER :: pPlag  
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ModRkInit.F90,v $ $Revision: 1.8 $'

  global => region%global

  CALL RegisterFunction(global,'PLAG_RkInitPrimary',&
  'PLAG_ModRkInit.F90')

!******************************************************************************
! Set pointers
!******************************************************************************

#ifdef RFLO   
  iLev  = region%currLevel
  pPlag   => region%levels(iLev)%plag
#endif
#ifdef RFLU
  pPlag   => region%plag
#endif

  pAiv    => pPlag%aiv
  pArv    => pPlag%arv   
  pCv     => pPlag%cv
  
  pAivOld => pPlag%aivOld
  pArvOld => pPlag%arvOld
  pCvOld  => pPlag%cvOld
  
!******************************************************************************
! Get dimensions
!******************************************************************************

#ifdef RFLO
  nPcls = region%levels(iLev)%plag%nPcls 
#endif
#ifdef RFLU
  nPcls = region%plag%nPcls
#endif
  
  nAiv = pPlag%nAiv
  nArv = pPlag%nArv           
  nCv  = pPlag%nCv

!******************************************************************************
! Initialize previous solution
!******************************************************************************
  
  CALL PLAG_RkInitGenericInt(  region,iStage,1,nPcls,1,nAiv,pAiv,pAivOld )
  CALL PLAG_RkInitGenericReal( region,iStage,1,nPcls,1,nArv,pArv,pArvOld )
  CALL PLAG_RkInitGenericReal( region,iStage,1,nPcls,1,nCv ,pCv ,pCvOld  )

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RkInitPrimary







!******************************************************************************
!
! Purpose: Driver that initializes variables for tiles.
!
! Description: none.
!
! Input: region = current region.
!        iStage = current RK stage.
!
! Output: region%levels%tilePlag%cvOld 
!
! Notes: none.
!
!******************************************************************************

SUBROUTINE PLAG_InjcTileRkInit( region, iStage )


! *****************************************************************************
! Definitions and declarations
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: iStage
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: bcType,iPatch,nCv,nPatches,nTiles
#ifdef RFLO
  INTEGER :: iLev,n1,n2
#endif
  
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv, pCvOld

  TYPE(t_patch),     POINTER :: pPatch
  TYPE(t_tile_plag), POINTER :: pTilePlag  
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ModRkInit.F90,v $ $Revision: 1.8 $'

  global => region%global

  CALL RegisterFunction(global,'PLAG_InjcTileRkInit',&
  'PLAG_ModRkInit.F90')
  
!****************************************************************************** 
! Get dimensions
!****************************************************************************** 

#ifdef RFLO
  iLev     = region%currLevel
  nPatches = region%nPatches
#endif
#ifdef RFLU
  nPatches = region%grid%nPatches
#endif

!******************************************************************************     
! Loop over patches
!****************************************************************************** 

  DO iPatch=1,nPatches

#ifdef RFLO
    pPatch => region%levels(iLev)%patches(iPatch)
#endif
#ifdef RFLU
    pPatch => region%patches(iPatch)
#endif

    bcType = pPatch%bcType

! =============================================================================
!   Select injection boundary condition
! =============================================================================

#ifdef RFLU
    IF ( (bcType >= BC_INJECTION .AND. bcType <= BC_INJECTION + BC_RANGE) .OR. &
         (bcType >= BC_INFLOW    .AND. bcType <= BC_INFLOW    + BC_RANGE)      ) THEN 
#else
    IF ( (bcType >= BC_INJECTION .AND. bcType <= BC_INJECTION + BC_RANGE)) THEN
#endif

!  ----------------------------------------------------------------------------   
!     Get tile dimensions and set pointers
!  ----------------------------------------------------------------------------   

#ifdef RFLO
      n1      = ABS(pPatch%l1end   -pPatch%l1beg   ) + 1    
      n2      = ABS(pPatch%l2end   -pPatch%l2beg   ) + 1
      nTiles  = n1*n2
#endif
#ifdef RFLU
      nTiles  = pPatch%nBFaces
#endif
      
      pTilePlag => pPatch%tilePlag
      
      nCv  = pTilePlag%nCv

      pCv     => pTilePlag%cv
      pCvOld  => pTilePlag%cvOld

!  ----------------------------------------------------------------------------  
!     Initialize previous solution
!  ---------------------------------------------------------------------------- 
      
      CALL PLAG_RkInitGenericReal( region,iStage,1,nTiles,1,nCv,pCv,pCvOld  )

    END IF !bcType
          
  END DO ! iPatch

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_InjcTileRkInit

! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE PLAG_ModRkInit

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModRkInit.F90,v $
! Revision 1.8  2009/03/02 00:19:36  mtcampbe
! Added some ifdefs around Rocflo to disable particle injection on INFLOW
! boundaries and added some checks around MPI tags utilizing a new global
! data item, global%mpiTagMax.
!
! Revision 1.7  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/09/18 20:30:18  fnajjar
! Activated tile datastructure for inflow bc
!
! Revision 1.4  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.3  2005/05/26 23:02:05  fnajjar
! Bug fix in setting iLev before defining pPlag pointer
!
! Revision 1.2  2005/05/23 18:42:15  fnajjar
! Bug fix to define pointers before setting dimensions
!
! Revision 1.1  2005/05/19 16:02:43  fnajjar
! Initial import
!
!******************************************************************************










