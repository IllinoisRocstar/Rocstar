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
! Purpose: Set variable info for mixture.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_SetVarInfo.F90,v 1.6 2008/12/06 08:44:13 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_SetVarInfo(pRegion)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
      
  USE ModInterfaces, ONLY: RFLU_GetCvLoc    
      
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: cvMixtPres,cvMixtXVel,cvMixtYVel,cvMixtZVel
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_SetVarInfo.F90,v $ $Revision: 1.6 $'

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_SetVarInfo',&
  'RFLU_SetVarInfo.F90')

! ******************************************************************************
! Set variable info depending on fluid model
! ******************************************************************************

  SELECT CASE ( pRegion%mixtInput%fluidModel ) 

! ==============================================================================
!   Incompressible fluid
! ==============================================================================    

    CASE ( FLUID_MODEL_INCOMP ) 
      cvMixtXVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_XVEL)
      cvMixtYVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_YVEL)
      cvMixtZVel = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_ZVEL)
      cvMixtPres = RFLU_GetCvLoc(global,FLUID_MODEL_INCOMP,CV_MIXT_PRES)                  
      
      pRegion%mixt%cvInfo(cvMixtXVel) = VAR_INFO_POSNEG
      pRegion%mixt%cvInfo(cvMixtYVel) = VAR_INFO_POSNEG
      pRegion%mixt%cvInfo(cvMixtZVel) = VAR_INFO_POSNEG
      pRegion%mixt%cvInfo(cvMixtPres) = VAR_INFO_POS           
      
! ==============================================================================
!   Compressible fluid
! ==============================================================================    
    
    CASE ( FLUID_MODEL_COMP ) 
      pRegion%mixt%cvInfo(CV_MIXT_DENS) = VAR_INFO_POS
      pRegion%mixt%cvInfo(CV_MIXT_XVEL) = VAR_INFO_POSNEG
      pRegion%mixt%cvInfo(CV_MIXT_YVEL) = VAR_INFO_POSNEG
      pRegion%mixt%cvInfo(CV_MIXT_ZVEL) = VAR_INFO_POSNEG
      pRegion%mixt%cvInfo(CV_MIXT_PRES) = VAR_INFO_POS      

! ==============================================================================
!   Default
! ==============================================================================    

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! pRegion%mixtInput%fluidModel

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_SetVarInfo

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_SetVarInfo.F90,v $
! Revision 1.6  2008/12/06 08:44:13  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:26  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/03/26 20:21:47  haselbac
! Removed error trap for GL model
!
! Revision 1.3  2004/11/14 19:44:00  haselbac
! Generalized setting of variable info
!
! Revision 1.2  2004/11/06 03:17:34  haselbac
! Cosmetics only
!
! Revision 1.1  2004/11/02 02:26:50  haselbac
! Initial revision
!
! ******************************************************************************







