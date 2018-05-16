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
! Purpose: Wrapper routine for moving grid.
!
! Description: None. 
!
! Input: 
!   regions     Region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_MoveGridWrapper.F90,v 1.5 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_MoveGridWrapper(regions)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModDataStruct, ONLY: t_region  
  USE ModMPI
  USE ModParameters

#ifdef GENX
  USE RFLU_ModRocstarTools
#endif

  USE ModInterfaces, ONLY: RFLU_MoveGridDisp, & 
                           RFLU_MoveGridXyz

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************
  INTEGER :: iReg
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_MoveGridWrapper.F90,v $ $Revision: 1.5 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_MoveGridWrapper',&
  'RFLU_MoveGridWrapper.F90')

! ******************************************************************************
! Call appropriate grid motion routine
! ******************************************************************************
#ifdef GENX
  DO iReg = 1,global%nRegionsLocal
     pRegion => regions(iReg)
     IF(global%cnstrCaseRad > 0.0) THEN
        CALL RFLU_GENX_ConstrainDisp(pRegion)
     ENDIF
  END DO ! iReg
#endif
  IF ( regions(1)%mixtInput%moveGridType == MOVEGRID_TYPE_DISP ) THEN 
    CALL RFLU_MoveGridDisp(regions)  
  ELSE IF ( regions(1)%mixtInput%moveGridType == MOVEGRID_TYPE_XYZ ) THEN 
    CALL RFLU_MoveGridXyz(regions,MOVEGRID_CONTEXT_MOVESMOOTH)
#ifdef GENX    
  ELSE IF ( regions(1)%mixtInput%moveGridType == MOVEGRID_TYPE_GENX ) THEN 
    CALL RFLU_GENX_MoveGrid(regions)
#endif    
  ELSE 
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END IF ! regions(1)%mixtInput%moveGridType

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_MoveGridWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_MoveGridWrapper.F90,v $
! Revision 1.5  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:42  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/04/14 14:25:40  mtcampbe
! Updated for rocket case constraints
!
! Revision 1.2  2004/10/19 19:29:22  haselbac
! Added GENX grid motion option, cosmetics
!
! Revision 1.1  2003/03/31 16:05:39  haselbac
! Initial revision
!
! ******************************************************************************







