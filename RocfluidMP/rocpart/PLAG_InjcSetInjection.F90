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
! Purpose: set the next particle diameter and superparticle loading
!          for the multiphase injection algorithm.
!
! Description: none.
!
! Input: injcDiamDist  = injection model type
!        diam       = particle diameter
!        spLoad     = superparticle loading
!
! Output: diam and spLoad 
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_InjcSetInjection.F90,v 1.3 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_InjcSetInjection( region, pTilePlag, iTile, tCoeff,  tSum, &
                                  poolVol,  injectQ, ratio )

  USE ModDataTypes  
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModRandom, ONLY     : Rand1Uniform
  USE ModPartLag, ONLY    : t_tile_plag
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region)             :: region
  TYPE(t_tile_plag), POINTER :: pTilePlag
  
  INTEGER      :: iTile
  LOGICAL      :: injectQ 
    
  REAL(RFREAL) :: poolVol, ratio, tCoeff, tSum

! ... loop variables

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  
  REAL(RFREAL) :: coefRatio, pi, poolVolNew, randUnif, volPart
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pDv
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InjcSetInjection.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global, 'PLAG_InjcSetInjection',&
  'PLAG_InjcSetInjection.F90' )

! Set pointer -----------------------------------------------------------------

  pDv => pTilePlag%dv
  
! Set particle volume and new pool volume -------------------------------------

  pi = global%pi
  
  volPart    = pi/6.0_RFREAL* pDv(DV_TILE_DIAM,iTile)**3  
  poolVolNew = poolVol - pDv(DV_TILE_SPLOAD,iTile)*volPart
  
  ratio   = 0.0_RFREAL  
  injectQ = .FALSE.
  
  IF ( poolVolNew > 0.0_RFREAL ) THEN
  
    coefRatio = tCoeff/poolVolNew
    tSum      = tSum + coefRatio*pDv(DV_TILE_COUNTDOWN,iTile) 

    IF ( tSum <= 1.0_RFREAL ) THEN
      injectQ = .TRUE.
      ratio   = volPart/poolVol
      
      randUnif = Rand1Uniform(region%randData) 
      IF ( randUnif <= 0.0_RFREAL) THEN         
        pDv(DV_TILE_COUNTDOWN,iTile) = 50.0_RFREAL 
      ELSE
        pDv(DV_TILE_COUNTDOWN,iTile) = -LOG(randUnif)
      END IF ! randUnif
      
    ELSE 
    
! - Decrease time factor to account for time waited in this step --------------

      pDv(DV_TILE_COUNTDOWN,iTile) = (tSum-1.0_RFREAL)/coefRatio 
    
    END IF ! tSum
    
  END IF !poolVolNew  
            
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_InjcSetInjection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InjcSetInjection.F90,v $
! Revision 1.3  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:42  fnajjar
! Initial revision after changing case
!
! Revision 1.3  2004/06/16 23:06:33  fnajjar
! Renamed variabled for CRE kernel
!
! Revision 1.2  2003/11/21 22:35:51  fnajjar
! Update Random Number Generator
!
! Revision 1.1  2002/10/25 14:16:31  f-najjar
! Initial Import of Rocpart
!
!******************************************************************************







