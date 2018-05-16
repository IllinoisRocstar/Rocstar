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
! Purpose: Check input values for species.
!
! Description: None.
!
! Input: 
!   regions             Data associated with regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: SPEC_CheckUserInput.F90,v 1.6 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_CheckUserInput(regions)

  USE ModGlobal, ONLY: t_global
  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModError
  USE ModParameters
          
  USE ModSpecies, ONLY: t_spec_type          
            
  USE ModInterfaces, ONLY: MixtPerf_R_M            
              
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:), POINTER :: regions

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: iReg,iSpec
  REAL(RFREAL) :: cp,g,gc
  TYPE(t_global), POINTER :: global
  TYPE(t_spec_type), POINTER :: pSpecType  

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: SPEC_CheckUserInput.F90,v $ $Revision: 1.6 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'SPEC_CheckUserInput',&
  'SPEC_CheckUserInput.F90')

! ******************************************************************************
! Check input values
! ******************************************************************************

  DO iReg = LBOUND(regions,1),UBOUND(regions,1)

! ==============================================================================
!   Gas properties
! ==============================================================================

    DO iSpec = 1,regions(iReg)%specInput%nSpecies
      pSpecType => regions(iReg)%specInput%specType(iSpec)
      
      gc = MixtPerf_R_M(pSpecType%pMaterial%molw)
      cp = pSpecType%pMaterial%spht
      g  = cp/(cp-gc)
      
      IF ( (g < 1.0_RFREAL) .OR. (g > 2.0_RFREAL) ) THEN 
        CALL ErrorStop(global,ERR_SPEC_PROPS_INVALID,__LINE__)
      END IF ! pSpecType%pMaterial%spht
    END DO ! iSpec

! ==============================================================================
!   EEM velocity
! ==============================================================================

    SELECT CASE ( regions(iReg)%mixtInput%indSd )
      CASE ( 0 )
      CASE ( 1 )
        IF ( regions(iReg)%mixtInput%computeTv .EQV. .FALSE. ) THEN
          CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__, &
                         'Attempting to use EEM without viscosity')
        END IF ! regions(iReg)%mixtInput%computeTv

        IF ( regions(iReg)%mixtInput%moveGrid .EQV. .TRUE. ) THEN
          CALL ErrorStop(global,ERR_UNKNOWN_OPTION,__LINE__, &
                         'EEM not yet implemented with moving grid')
        END IF ! regions(iReg)%mixtInput%moveGrid
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! regions(iReg)%mixtInput%indSd
  END DO ! iReg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_CheckUserInput

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_CheckUserInput.F90,v $
! Revision 1.6  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:52  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.3  2005/11/10 02:33:38  haselbac
! Clean-up, added checks for gas properties
!
! Revision 1.2  2004/07/30 22:47:37  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.1  2003/11/25 21:08:37  haselbac
! Initial revision
!
! ******************************************************************************







