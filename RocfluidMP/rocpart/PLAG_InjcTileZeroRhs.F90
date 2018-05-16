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
! Purpose: Zero out residuals for tile data infrastructure.
!
! Description: none.
!
! Input: region = current region.
!
! Output: region%levels%tilePlag%rhs
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_InjcTileZeroRhs.F90,v 1.4 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_injcTileZeroRhs( region )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag_input, t_tile_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt_input
  USE ModError  
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  
  INTEGER        :: iStage

! ... loop variables
  INTEGER :: iCont, iPatch, iTile

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER     :: bcType, nCont, nPatches, nTiles 
#ifdef RFLO
  INTEGER     :: iLev, n1, n2
#endif
  INTEGER, POINTER, DIMENSION(:)        :: pCvTileMass
  
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pRhs
  
  TYPE(t_patch),     POINTER :: pPatch
  TYPE(t_tile_plag), POINTER :: pTilePlag    
  TYPE(t_global),    POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InjcTileZeroRhs.F90,v $ $Revision: 1.4 $'

  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_InjcTileZeroRhs',&
  'PLAG_InjcTileZeroRhs.F90' )
  
! Get dimensions --------------------------------------------------------------

#ifdef RFLO
  iLev     = region%currLevel  
  nPatches = region%nPatches
#endif
#ifdef RFLU
  nPatches = region%grid%nPatches
#endif
  nCont    = region%plagInput%nCont 
 
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
    
! -- Get tile dimensions and pointers -----------------------------------------

#ifdef RFLO
      n1      = ABS(pPatch%l1end   -pPatch%l1beg   ) + 1    
      n2      = ABS(pPatch%l2end   -pPatch%l2beg   ) + 1
      nTiles  = n1*n2
#endif
#ifdef RFLU
      nTiles  = pPatch%nBFaces
#endif
      
      pTilePlag   => pPatch%tilePlag
        
      pCvTileMass => pTilePlag%cvTileMass
      pRhs        => pTilePlag%rhs

! -- Zero out residuals -------------------------------------------------------

      DO iTile = 1, nTiles 
        pRhs(CV_TILE_MOMNRM,iTile) = 0.0_RFREAL
        pRhs(CV_TILE_ENER  ,iTile) = 0.0_RFREAL

        DO iCont = 1, nCont
          pRhs( pCvTileMass(iCont),iTile) = 0.0_RFREAL                               
        END DO ! iCont
 
      END DO !iTile
       
    END IF !bcType
    
  END DO ! iPatch
  
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_InjcTileZeroRhs

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InjcTileZeroRhs.F90,v $
! Revision 1.4  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 20:57:49  fnajjar
! Initial revision after changing case
!
! Revision 1.4  2004/03/08 22:25:13  fnajjar
! Modified routine to be RFLU-aware
!
! Revision 1.3  2004/02/25 21:56:15  fnajjar
! Moved tile pointers outside do-loop
!
! Revision 1.2  2003/01/16 20:15:11  f-najjar
! Removed iRegionGlobal
!
! Revision 1.1  2002/10/25 14:16:32  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************







