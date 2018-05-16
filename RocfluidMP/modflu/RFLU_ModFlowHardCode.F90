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
! Purpose: Suite of routines to support hard-coded flow solutions.
!
! Description: None.
!
! Notes: 
!   1. Collected routines in single module because setting of parameters is
!      needed in at least three places: initialization of solution, setting 
!      of boundary profiles, and computation of errors. 
!
! ******************************************************************************
!
! $Id: RFLU_ModFlowHardCode.F90,v 1.6 2008/12/06 08:44:21 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModFlowHardCode

  USE ModParameters
  USE ModDataTypes
      
  IMPLICIT NONE
     
  PRIVATE
  PUBLIC :: RFLU_GetParamsHardCodePAcoust, & 
            RFLU_GetParamsHardCodeProudman, &
            RFLU_GetParamsHardCodeRingleb, & 
            RFLU_GetParamsHardCodeSsVortex

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
   
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModFlowHardCode.F90,v $ $Revision: 1.6 $'        
             
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  


! ******************************************************************************
!
! Purpose: Get parameters for pipe acoustics.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   pTot        Total pressure
!   aTot        Total speed of sound
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_GetParamsHardCodePAcoust(pTot,aTot)

    REAL(RFREAL), INTENT(OUT) :: aTot,pTot
    
    aTot = 340.0_RFREAL
    pTot = 1.0E+5_RFREAL
    
  END SUBROUTINE RFLU_GetParamsHardCodePAcoust
  
  
  
! ******************************************************************************
!
! Purpose: Get parameters for Proudman-Culick flow.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   dInc        Density
!   mInj        Injection mass flux
!   vInj        Injection velocity
!   pTot        Total pressure
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_GetParamsHardCodeProudman(dInc,mInj,vInj,pTot)

    REAL(RFREAL), INTENT(OUT) :: dInc,mInj,pTot,vInj
    
    dInc = 1.0_RFREAL
    mInj = 2.42_RFREAL
    vInj = mInj/dInc
    pTot = 1.0E+5_RFREAL
    
  END SUBROUTINE RFLU_GetParamsHardCodeProudman



! ******************************************************************************
!
! Purpose: Get parameters for Ringleb flow.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   pTot        Total pressure
!   tTot        Total temperature
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_GetParamsHardCodeRingleb(pTot,tTot)

    REAL(RFREAL), INTENT(OUT) :: pTot,tTot
    
    pTot = 1.0E+5_RFREAL
    tTot = 288.15_RFREAL
    
  END SUBROUTINE RFLU_GetParamsHardCodeRingleb



! ******************************************************************************
!
! Purpose: Get parameters for supersonic vortex flow.
!
! Description: None.
!
! Input: None.
!
! Output: 
!   ri          Inner radius
!   Mi          Mach number at inner radius    
!   pTot        Total pressure
!   tTot        Total temperature
!
! Notes: None.
!
! ******************************************************************************
  
  SUBROUTINE RFLU_GetParamsHardCodeSsVortex(ri,Mi,pTot,tTot)

    REAL(RFREAL), INTENT(OUT) :: Mi,pTot,ri,tTot
    
    ri   = 1.0_RFREAL 
    Mi   = 2.25_RFREAL
    pTot = 1.0E+5_RFREAL
    tTot = 288.15_RFREAL
    
  END SUBROUTINE RFLU_GetParamsHardCodeSsVortex




! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModFlowHardCode


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModFlowHardCode.F90,v $
! Revision 1.6  2008/12/06 08:44:21  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:32  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.3  2005/03/15 20:45:06  haselbac
! Added routine to get parameters for pipe acoustics
!
! Revision 1.2  2004/07/06 15:14:39  haselbac
! Cosmetics only
!
! Revision 1.1  2004/02/23 23:01:49  haselbac
! Initial revision
!
! ******************************************************************************






