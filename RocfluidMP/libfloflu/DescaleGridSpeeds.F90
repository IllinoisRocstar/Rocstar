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
! Purpose: Descale grid speeds.
!
! Description: None. 
!
! Input: 
!   region     Pointer to region
!
! Output: None.
!
! Notes: 
!   1. Grid speeds need to be descaled because they are scaled during the 
!      Runge-Kutta stages, but for consistent restarts, they need to be 
!      written out consistent with their values when they were used for the 
!      time-step computation.
!
!******************************************************************************
!
! $Id: DescaleGridSpeeds.F90,v 1.3 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE DescaleGridSpeeds( region )

  USE ModDataTypes
  USE ModGrid, ONLY      : t_grid
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY  : t_patch
  USE ModGlobal, ONLY    : t_global
  USE ModError
  USE ModParameters
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset
#include "Indexing.h"
#endif
#ifdef RFLU
  USE RFLU_ModGrid
#endif

  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  INTEGER :: ifc, iPatch, irk, irkStep
  REAL(RFREAL) :: scaleFactor, term
  REAL(RFREAL) :: ark(5), grk(5)
#ifdef RFLO
  INTEGER :: iLev, idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ibn, ien, iNOff, ijNOff, in
  REAL(RFREAL), POINTER :: siVel(:), sjVel(:), skVel(:)
#endif
#ifdef RFLU
  REAL(RFREAL), DIMENSION(:), POINTER :: gs
  TYPE(t_patch), POINTER :: pPatch
#endif
  TYPE(t_global), POINTER :: global

! *****************************************************************************
! Start
! *****************************************************************************

  global => region%global

  CALL RegisterFunction(global,'DescaleGridSpeeds',&
  'DescaleGridSpeeds.F90')
  
! *****************************************************************************
! Set pointers and variables
! *****************************************************************************

  irkStep = region%irkStep

  ark(:) = region%mixtInput%ark(:)
  grk(:) = region%mixtInput%grk(:)

! *****************************************************************************
! Descale grid speeds
! *****************************************************************************

! =============================================================================
! Determine scaling factor (compare with ScaleGridSpeeds.F90)
! ============================================================================= 

  term = 0.0_RFREAL

  DO irk = 1,global%nrkSteps-1
    term = term + grk(irk)/ark(irk)
  END DO ! irk

  scaleFactor = 1.0_RFREAL/(1.0_RFREAL/ark(global%nrkSteps) - term)

! =============================================================================
! Interior faces
! ============================================================================= 

#ifdef RFLO
  iLev = region%currLevel
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

  siVel => region%levels(iLev)%grid%siVel
  sjVel => region%levels(iLev)%grid%sjVel
  skVel => region%levels(iLev)%grid%skVel

  DO in=ibn,ien
    siVel(in) = scaleFactor*siVel(in)
    sjVel(in) = scaleFactor*sjVel(in)
    skVel(in) = scaleFactor*skVel(in)
  ENDDO
#endif

#ifdef RFLU
  gs => region%grid%gs

  DO ifc = 1,region%grid%nFaces
    gs(ifc) = scaleFactor*gs(ifc)
  END DO ! ifc

! =============================================================================
! Patch faces
! ============================================================================= 

  DO iPatch = 1,region%grid%nPatches
    pPatch => region%patches(iPatch)
    gs     => pPatch%gs

    DO ifc = 1,pPatch%nBFaces
      gs(ifc) = scaleFactor*gs(ifc)      
    END DO ! ifc
  END DO ! iPatch
#endif

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE DescaleGridSpeeds

!******************************************************************************
!
! RCS Revision history:
!
! $Log: DescaleGridSpeeds.F90,v $
! Revision 1.3  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:48:29  haselbac
! Initial revision after changing case
!
! Revision 1.7  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.4  2003/10/01 23:52:09  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.3  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.2  2003/03/15 16:14:30  haselbac
! Changed loop limit
!
! Revision 1.1  2003/02/25 21:42:24  haselbac
! Initial revision
!
!******************************************************************************







