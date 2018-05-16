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
! Purpose: initializes Runge-Kutta data for particles
!
! Description: none.
!
! Input: none.
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_RkInit.F90,v 1.5 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RkInit( region,iStage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartLag,    ONLY : t_plag
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  USE PLAG_ModInterfaces, ONLY: PLAG_InjcTileZeroRhs, &
                                PLAG_ZeroRhs
  USE PLAG_ModRkInit
  
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

  INTEGER, INTENT(IN) :: iStage

! ... loop variables

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

#ifdef RFLO
  INTEGER :: iLev
#endif

  TYPE(t_plag) ,  POINTER :: pPlag
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RkInit.F90,v $ $Revision: 1.5 $'

  global => region%global

  CALL RegisterFunction( global,'PLAG_RkInit',&
  'PLAG_RkInit.F90' )

  IF (global%plagUsed) THEN

#ifdef RFLO
    iLev  =  region%currLevel
    pPlag => region%levels(iLev)%plag
#endif
#ifdef RFLU
    pPlag => region%plag
#endif

    IF ( pPlag%nPcls > 0 ) THEN
      CALL PLAG_ZeroRhs( region )
      CALL PLAG_RkInitPrimary( region, iStage )
    END IF ! nPcls

  END IF ! plagUsed

! Main algorithm for Tile evolution -------------------------------------------

! - Zero out residuals --------------------------------------------------------

  CALL PLAG_InjcTileZeroRhs( region )

! - Load old values -----------------------------------------------------------

  CALL PLAG_InjcTileRkInit( region, iStage )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RkInit

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RkInit.F90,v $
! Revision 1.5  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/05/19 16:02:20  fnajjar
! Added calls to generic routines for initialization
!
! Revision 1.2  2005/04/22 18:34:32  fnajjar
! Removed IF statements to properly handle RK3 scheme
!
! Revision 1.1  2004/12/01 20:58:16  fnajjar
! Initial revision after changing case
!
! Revision 1.8  2004/03/22 23:49:28  fnajjar
! Removed RFLO-based ifdef for tile kernel
!
! Revision 1.7  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2004/02/26 21:02:21  haselbac
! Commented out RFLO-specific tile init routines
!
! Revision 1.5  2004/02/06 21:18:20  fnajjar
! Initial Integration of Rocpart with Rocflu
!
! Revision 1.4  2003/11/12 21:33:07  fnajjar
! Moved FC computations to PLAG_RFLO_SetMetrics
!
! Revision 1.3  2003/11/03 22:39:25  fnajjar
! Added call to copy face vectors
!
! Revision 1.2  2003/03/28 19:51:50  fnajjar
! Included initialization step for tiles
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!
!******************************************************************************







