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
! Purpose: initialize the tiles for the multiphase injection algorithm.
!
! Description: none.
!
! Input: region    = current region
!
! Output: regions(iReg)%levels%patch%tile = initial tile values
!                                           of current region.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_InjcTileInitialize.F90,v 1.4 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_InjcTileInitialize( region )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag_input, t_tile_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModRandom, ONLY     : Rand1Uniform
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetPatchIndices, RFLO_GetPatchDirection, &
                            RFLO_GetCellOffset, RFLO_GetNodeOffset
#endif
  USE PLAG_ModInterfaces, ONLY : PLAG_injcMakeParticle
  USE ModError
  USE ModParameters
  USE ModMPI
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iPatch, iTile

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: bcType, nPatches, nTiles
#ifdef RFLO
  INTEGER :: iLev, n1, n2
#endif
  INTEGER :: injcDiamDist
  INTEGER :: nCv, nDv
  INTEGER :: ejecModel,iCont,nCont
  
  REAL(RFREAL) :: randUnif
  REAL(RFREAL) :: poolVolumeInit, ratioPhiDensInv   
  REAL(RFREAL), POINTER, DIMENSION(:) :: injcMassFluxRatio, densPlag

  TYPE(t_patch),     POINTER :: pPatch
  TYPE(t_tile_plag), POINTER :: pTilePlag
  TYPE(t_global),    POINTER :: global
    
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InjcTileInitialize.F90,v $ $Revision: 1.4 $'

  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_InjcTileInitialize',&
  'PLAG_InjcTileInitialize.F90' )

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Initializing tile memory for PLAG...'
  END IF ! global%verbLevel
  
! Get dimensions --------------------------------------------------------------

#ifdef RFLO
  iLev     = region%currLevel
  nPatches = region%nPatches
#endif
#ifdef RFLU
  nPatches = region%grid%nPatches
#endif

  nCont        = region%plagInput%nCont  
  injcDiamDist = region%plagInput%injcDiamDist
  ejecModel    = region%plagInput%ejecModel                

  injcMassFluxRatio => region%plagInput%injcMassFluxRatio
  densPlag          => region%plagInput%dens

! Loop over patches -----------------------------------------------------------

  DO iPatch=1,nPatches

#ifdef RFLO
    pPatch  => region%levels(iLev)%patches(iPatch)
#endif
#ifdef RFLU
    pPatch  => region%patches(iPatch)
#endif

    bcType = pPatch%bcType

! - Select injection boundary condition ---------------------------------------

    IF (bcType>=BC_INJECTION .AND. bcType<=BC_INJECTION+BC_RANGE) THEN

! - Get dimensions -----------------------------------------------------------

#ifdef RFLO
      n1      = ABS(pPatch%l1end   -pPatch%l1beg   ) + 1    
      n2      = ABS(pPatch%l2end   -pPatch%l2beg   ) + 1
      nTiles  = n1*n2
#endif
#ifdef RFLU
      nTiles = pPatch%nBFaces
#endif

      pTilePlag => pPatch%tilePlag

      nCv = pTilePlag%nCv
      nDv = pTilePlag%nDv

      pTilePlag%nPclsInjc(:nTiles)   = 0
      
      pTilePlag%cv(:nCv,:nTiles)     = 0.0_RFREAL
      pTilePlag%dv(:nDv,:nTiles)     = 0.0_RFREAL
        
      pTilePlag%rhs(:nCv,:nTiles)    = 0.0_RFREAL
        
      pTilePlag%cvOld(:nCv,:nTiles)  = 0.0_RFREAL
      pTilePlag%rhsSum(:nCv,:nTiles) = 0.0_RFREAL
                        
! --- Loop over tile patch ----------------------------------------------------

      DO iTile = 1, nTiles

        CALL PLAG_injcMakeParticle(region, injcDiamDist,         &
                                   pTilePlag%dv(DV_TILE_DIAM,iTile), &
                                   pTilePlag%dv(DV_TILE_SPLOAD,iTile) )
 
        randUnif = Rand1Uniform(region%randData) 

! ---- Avoid error if randUnif is 0 -------------------------------------------
! ---- and set timefactor such that EXP(-50) = 1.9E-22 ------------------------

        IF ( randUnif <= 0.0_RFREAL) THEN          
          pTilePlag%dv(DV_TILE_COUNTDOWN,iTile) = 50.0_RFREAL 
        ELSE
          pTilePlag%dv(DV_TILE_COUNTDOWN,iTile) = -LOG(randUnif)
        END IF ! randUnif

! ---- Set initial pool volume

       IF ( ejecModel == PLAG_EJEC_CRE ) THEN
         poolVolumeInit  = 0.0_RFREAL
         ratioPhiDensInv = poolVolumeInit  &
                         *1.0_RFREAL/SUM(injcMassFluxRatio(:)/densPlag(:))

         DO iCont = 1, nCont
           pTilePlag%cv(pTilePlag%cvTileMass(iCont),iTile) = poolVolumeInit    &
                                                           * injcMassFluxRatio(iCont)
         END DO ! iCont
       ENDIF ! ejecModel       
      END DO ! iTile
                                       
    ENDIF  ! bcType

  ENDDO    ! iPatch
  
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_injcTileInitialize

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InjcTileInitialize.F90,v $
! Revision 1.4  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 20:57:45  fnajjar
! Initial revision after changing case
!
! Revision 1.8  2004/06/30 15:43:27  fnajjar
! Added initialization step for CRE model
!
! Revision 1.7  2004/06/16 23:06:33  fnajjar
! Renamed variabled for CRE kernel
!
! Revision 1.6  2004/03/03 00:30:52  fnajjar
! Incorrect ifdef construct for RFLU pPatch pointer
!
! Revision 1.5  2003/11/26 22:00:28  fnajjar
! Removed improper comment symbol after USE ModRandom
!
! Revision 1.4  2003/11/21 22:35:51  fnajjar
! Update Random Number Generator
!
! Revision 1.3  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.2  2003/01/16 20:15:11  f-najjar
! Removed iRegionGlobal
!
! Revision 1.1  2002/10/25 14:16:31  f-najjar
! Initial Import of Rocpart
!
!******************************************************************************







