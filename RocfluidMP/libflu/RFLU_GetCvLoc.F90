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
! Purpose: Get variable location for mixture conserved variables.
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
! $Id: RFLU_GetCvLoc.F90,v 1.4 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2004-2006 by the University of Illinois
!
! ******************************************************************************

INTEGER FUNCTION RFLU_GetCvLoc(global,fluidModel,var)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
      
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: fluidModel,var
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_GetCvLoc.F90,v $ $Revision: 1.4 $'

  CALL RegisterFunction(global,'RFLU_GetCvLoc',&
  'RFLU_GetCvLoc.F90')

! ******************************************************************************
! Set variable info depending on fluid model
! ******************************************************************************

  SELECT CASE ( fluidModel ) 

! ==============================================================================
!   Incompressible fluid
! ==============================================================================    

    CASE ( FLUID_MODEL_INCOMP ) 
      SELECT CASE ( var )
        CASE ( CV_MIXT_XVEL )  
          RFLU_GetCvLoc = 1
        CASE ( CV_MIXT_YVEL ) 
          RFLU_GetCvLoc = 2
        CASE ( CV_MIXT_ZVEL )
          RFLU_GetCvLoc = 3
        CASE ( CV_MIXT_PRES ) 
          RFLU_GetCvLoc = 4
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! var
      
! ==============================================================================
!   Compressible fluid
! ==============================================================================    
    
    CASE ( FLUID_MODEL_COMP ) 
      SELECT CASE ( var )
        CASE ( CV_MIXT_DENS ) 
          RFLU_GetCvLoc = 1
        CASE ( CV_MIXT_XMOM:CV_MIXT_XVEL )  
          RFLU_GetCvLoc = 2        
        CASE ( CV_MIXT_YMOM:CV_MIXT_YVEL ) 
          RFLU_GetCvLoc = 3     
        CASE ( CV_MIXT_ZMOM:CV_MIXT_ZVEL )
          RFLU_GetCvLoc = 4       
        CASE ( CV_MIXT_ENER:CV_MIXT_PRES ) 
          RFLU_GetCvLoc = 5    
        CASE DEFAULT 
          CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
      END SELECT ! var   

! ==============================================================================
!   Default
! ==============================================================================    

    CASE DEFAULT
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END SELECT ! fluidModel

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_GetCvLoc

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_GetCvLoc.F90,v $
! Revision 1.4  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/03/26 20:21:33  haselbac
! Removed error trap for GL model
!
! Revision 1.1  2004/11/06 03:16:43  haselbac
! Initial revision
!
! ******************************************************************************







