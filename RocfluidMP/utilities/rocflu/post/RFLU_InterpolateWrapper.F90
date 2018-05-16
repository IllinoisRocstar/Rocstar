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
! Purpose: Compute variables other than conserved variables at vertices
!
! Description: None.
!
! Input:
!   pRegion        Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_InterpolateWrapper.F90,v 1.15 2008/12/06 08:45:05 mtcampbe Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_InterpolateWrapper(pRegion)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid

  USE RFLU_ModInterpolation, ONLY: RFLU_InterpCells2Verts, &
                                   RFLU_InterpSimpleCells2Verts

  USE ModInterfaces, ONLY: RFLU_ComputeVertexVariables

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Local variables
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InterpolateWrapper.F90,v $ $Revision: 1.15 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_InterpolateWrapper', &
                        'RFLU_InterpolateWrapper.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pGrid => pRegion%grid

! ******************************************************************************
! Mixture
! ******************************************************************************

! ==============================================================================
! Conserved variables
! ==============================================================================

  IF ( global%postInterpType == INTERP_TYPE_PROPER ) THEN
    CALL RFLU_InterpCells2Verts(pRegion,global%postInterpOrder, &
                                pRegion%mixtInput%nCv,pRegion%mixt%cv, &
                                pRegion%mixt%cvVert)
  ELSE IF ( global%postInterpType == INTERP_TYPE_SIMPLE ) THEN
    CALL RFLU_InterpSimpleCells2Verts(pRegion,pRegion%mixtInput%nCv, & 
                                      pRegion%mixt%cv,pRegion%mixt%cvVert)
  ELSE
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END IF ! global%postInterpType

! ==============================================================================
! Gas variables
! ==============================================================================

  SELECT CASE ( pRegion%mixtInput%fluidModel )
    CASE ( FLUID_MODEL_COMP ) 
      IF ( pRegion%mixtInput%nGvAct == 0 ) THEN
        pRegion%mixt%gvVert(GV_MIXT_CP ,:) = pRegion%mixt%gv(GV_MIXT_CP ,:)
        pRegion%mixt%gvVert(GV_MIXT_MOL,:) = pRegion%mixt%gv(GV_MIXT_MOL,:)
      ELSE 
        IF ( global%postInterpType == INTERP_TYPE_PROPER ) THEN
          CALL RFLU_InterpCells2Verts(pRegion,global%postInterpOrder, &
                                      pRegion%mixtInput%nGv,pRegion%mixt%gv, &
                                      pRegion%mixt%gvVert)
        ELSE IF ( global%postInterpType == INTERP_TYPE_SIMPLE ) THEN
          CALL RFLU_InterpSimpleCells2Verts(pRegion,pRegion%mixtInput%nGv, & 
                                            pRegion%mixt%gv,pRegion%mixt%gvVert)
        ELSE
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
        END IF ! global%postInterpType
      END IF ! pRegion%mixtInput%nGvAct
    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%fluidModel 

! ******************************************************************************
! Physical modules
! ******************************************************************************

#ifdef SPEC
  IF ( global%specUsed .EQV. .TRUE. ) THEN

! ==============================================================================
!   Conserved variables
! ==============================================================================

    IF ( global%postInterpType == INTERP_TYPE_PROPER ) THEN
      CALL RFLU_InterpCells2Verts(pRegion,global%postInterpOrder, &
                                  pRegion%specInput%nSpecies, &
                                  pRegion%spec%cv,pRegion%spec%cvVert)
    ELSE IF ( global%postInterpType == INTERP_TYPE_SIMPLE ) THEN
      CALL RFLU_InterpSimpleCells2Verts(pRegion,pRegion%specInput%nSpecies, &
                                        pRegion%spec%cv,pRegion%spec%cvVert)
    ELSE
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! global%postInterpType
  END IF ! global%specUsed
#endif

! ******************************************************************************
! Compute remaining vertex variables
! ******************************************************************************

  CALL RFLU_ComputeVertexVariables(pRegion)

! ******************************************************************************
! Interpolate plotting variables from cell centroids to vertices
! ******************************************************************************

  IF ( pRegion%mixtInput%fluidModel == FLUID_MODEL_COMP ) THEN 
    IF ( global%postInterpType == INTERP_TYPE_PROPER ) THEN
      CALL RFLU_InterpCells2Verts(pRegion,global%postInterpOrder, &
                                  pRegion%plot%nPv,pRegion%plot%pv, &
                                  pRegion%plot%pvVert)
    ELSE IF ( global%postInterpType == INTERP_TYPE_SIMPLE ) THEN
      CALL RFLU_InterpSimpleCells2Verts(pRegion,pRegion%plot%nPv, & 
                                        pRegion%plot%pv,pRegion%plot%pvVert)
    ELSE
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! global%postInterpType
  END IF ! pRegion%mixtInput%fluidModel

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InterpolateWrapper

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InterpolateWrapper.F90,v $
! Revision 1.15  2008/12/06 08:45:05  mtcampbe
! Updated license.
!
! Revision 1.14  2008/11/19 22:18:16  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.13  2007/03/19 21:41:44  haselbac
! Adapted to changes related to plotting variables
!
! Revision 1.12  2006/03/26 20:22:33  haselbac
! Added support for GL model, changed order
!
! Revision 1.11  2005/11/14 17:05:01  haselbac
! Generalized to support pseudo-gas model
!
! Revision 1.10  2005/11/10 02:47:52  haselbac
! Added support for gas models
!
! Revision 1.9  2005/10/31 21:09:39  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.8  2005/05/28 18:02:33  haselbac
! Replaced CV_MIXT_ENER by pRegion%mixtInput%nCv
!
! Revision 1.7  2005/05/01 14:22:19  haselbac
! Added interp of plotting vars
!
! Revision 1.6  2004/11/14 19:57:14  haselbac
! Added IF for fluid model
!
! Revision 1.5  2004/07/28 15:29:21  jferry
! created global variable for spec use
!
! Revision 1.4  2004/07/21 14:56:59  haselbac
! Added capability of using simple interpolation
!
! Revision 1.3  2004/03/18 03:34:11  haselbac
! Changed call to renamed interpolation function
!
! Revision 1.2  2003/12/04 03:29:32  haselbac
! Adapted to changes in RFLU_ModInterpolation
!
! Revision 1.1  2003/11/25 21:04:05  haselbac
! Initial revision
!
! ******************************************************************************







