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
! Purpose: update step for injection tiles.
!
! Description: none.
!
! Input: region = current region
!        iReg   = current region number
!        iStage = current RK stage.
!
! Output: region%plag = plag variables
!
! Notes: This corresponds to Part I Step 2 in RocfluidMP framework.
!
!******************************************************************************
!
! $Id: PLAG_InjcTileUpdate.F90,v 1.3 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_injcTileUpdate( region, iReg, iStage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt_input
  USE ModPartLag, ONLY    : t_plag
  USE ModError
  USE ModParameters

  USE PLAG_ModParameters 
  
  USE PLAG_ModInterfaces, ONLY: PLAG_InjcTileRKUpdate        
#ifdef RFLO
  USE PLAG_ModInterfaces, ONLY: PLAG_InjcTileCalcRhs 
#endif
#ifdef RFLU
  USE PLAG_ModInterfaces, ONLY: PLAG_RFLU_InjcTileCalcRhs
#endif
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  
  INTEGER :: iReg, iStage
  
! ... loop variables

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_global), POINTER :: global
#ifdef RFLU
  TYPE(t_region), POINTER :: pRegion
#endif
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InjcTileUpdate.F90,v $ $Revision: 1.3 $'
  
  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_InjcTileUpdate',&
  'PLAG_InjcTileUpdate.F90' )

! - Calculate rhs for tile ----------------------------------------------------

#ifdef RFLO
!  WRITE(STDOUT,'(A)') '    Entering PLAG_InjcTileCalcRhs'
  CALL PLAG_InjcTileCalcRhs( region )
#endif

#ifdef RFLU
  pRegion => region%pRegion
!  WRITE(STDOUT,'(A)') '    Entering PLAG_RFLU_InjcTileCalcRhs'
  CALL PLAG_RFLU_InjcTileCalcRhs( pRegion )
#endif

! - Invoke RK update ----------------------------------------------------------

!  WRITE(STDOUT,'(A)') '    Entering PLAG_InjcTileRKUpdate'
  CALL PLAG_InjcTileRKUpdate( region, iStage )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_InjcTileUpdate

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InjcTileUpdate.F90,v $
! Revision 1.3  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:48  fnajjar
! Initial revision after changing case
!
! Revision 1.6  2004/03/08 22:24:39  fnajjar
! Modified routine to be RFLU-aware and added PLAG_RFLU_InjcTileCalcRhs call
!
! Revision 1.5  2004/02/25 21:56:15  fnajjar
! Moved tile pointers outside do-loop
!
! Revision 1.4  2003/05/01 22:51:12  fnajjar
! Removed PLAG_CalcFaceCentroids in PLAG_ModInterfaces list
!
! Revision 1.3  2003/04/15 23:02:38  fnajjar
! Removed dead code section
!
! Revision 1.2  2003/03/28 19:52:22  fnajjar
! Removed initialization step for tiles
!
! Revision 1.1  2003/02/04 19:09:46  f-najjar
! Initial Import
!
!******************************************************************************







