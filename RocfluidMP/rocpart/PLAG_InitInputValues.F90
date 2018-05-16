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
! Purpose: initialize user input parameters for Lagrangians particles
!          to default values.
!
! Description: none.
!
! Input: none.
!
! Output: regions = initial input values.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_InitInputValues.F90,v 1.5 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_InitInputValues( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartLag,    ONLY : t_plag_input
  USE ModError
  USE ModParameters
  USE ModMaterials, ONLY  : t_material
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_plag_input), POINTER :: pPlagInput
  TYPE(t_global),     POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InitInputValues.F90,v $ $Revision: 1.5 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PLAG_InitInputValues',&
  'PLAG_InitInputValues.F90' )

! global values ---------------------------------------------------------------

! none currently

! region related values -------------------------------------------------------

  DO iReg=LBOUND(regions,1),UBOUND(regions,1)

    pPlagInput => regions(iReg)%plagInput

! plagInput quantities

    pPlagInput%nPclsMax          = 1000
    pPlagInput%nPclsBuffTot      = 100
    pPlagInput%nPclsBuffCECellsMax = 10
    pPlagInput%ejecModel         = PLAG_EJEC_MODEL1
    pPlagInput%injcDiamDist      = PLAG_INJC_LOGNORM
    pPlagInput%injcVelRatio      = 0.0_RFREAL
    pPlagInput%spLoad            = 1.0_RFREAL
    pPlagInput%injcDiamMean      = 10.0E-06_RFREAL
    pPlagInput%injcDiamMin       = 1.0E-06_RFREAL
    pPlagInput%injcDiamMax       = 100.0E-06_RFREAL
    pPlagInput%injcStdDev        = 0.0_RFREAL
    pPlagInput%injcBeta          = 1.0_RFREAL
    pPlagInput%intrplMixtModel   = ZEROTH_ORDER
    pPlagInput%nCont             = 1
    pPlagInput%breakupModel      = PLAG_BREAKUP_NOMODEL
    pPlagInput%breakupFac        = 1.0_RFREAL
    pPlagInput%breakupWebSwi     = PLAG_BREAKUP_NOWEBSWI
    pPlagInput%readStatus        = -1 ! not read
    pPlagInput%findPclMethod     = FIND_PCL_METHOD_TRAJ_SAFE

! initial quantities

    pPlagInput%nPclsIni          = 0

  ENDDO  ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_InitInputValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InitInputValues.F90,v $
! Revision 1.5  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/03/06 23:13:13  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.2  2005/03/11 02:22:33  haselbac
! Changed default for tracking to safe method
!
! Revision 1.1  2004/12/01 20:57:36  fnajjar
! Initial revision after changing case
!
! Revision 1.14  2004/10/11 22:12:23  haselbac
! Bug fix
!
! Revision 1.13  2004/10/09 16:37:37  fnajjar
! Removed initialization of initFlag
!
! Revision 1.12  2004/10/08 22:12:43  haselbac
! Added initialization of findPclMethod
!
! Revision 1.11  2004/08/20 23:27:13  fnajjar
! Added Infrastructure for Plag prep tool
!
! Revision 1.10  2004/06/17 15:19:03  fnajjar
! Added infrastructure for ejection model
!
! Revision 1.9  2004/06/16 23:03:49  fnajjar
! Renamed variabled for CRE kernel
!
! Revision 1.8  2004/03/10 23:09:50  fnajjar
! Added maximum buffer size for corner-edge cells
!
! Revision 1.7  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.6  2004/02/26 21:02:14  haselbac
! Changed loop limits for generality
!
! Revision 1.5  2004/02/16 23:31:11  fnajjar
! Included default values for injcDiamMin and injcDiamMax
!
! Revision 1.4  2003/11/21 22:42:16  fnajjar
! Added plagActive
!
! Revision 1.3  2003/09/13 20:14:21  fnajjar
! Added infrastructure for Breakup model
!
! Revision 1.2  2003/09/10 23:35:50  fnajjar
! Removed flags that are subsumed with Rocinteract
!
! Revision 1.1  2003/04/14 14:33:15  fnajjar
! Initial import for proper input initialization
!
!******************************************************************************







