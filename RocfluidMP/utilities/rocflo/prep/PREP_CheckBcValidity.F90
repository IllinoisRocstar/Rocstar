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
! Purpose: check validity of BC input.
!
! Description: comparison made between .bc and .top files to check whether
!              all bc required in the .top file are defined in the .bc file.
!
! Input: regions = data of all regions.
!
! Output: error msg for inconsistency.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PREP_CheckBcValidity.F90,v 1.7 2008/12/06 08:44:50 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE CheckBcValidity( regions )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: regions(:)
  TYPE(t_patch), POINTER  :: patch

! ... loop variables
  INTEGER :: iReg, iPatch

! ... local variables
  INTEGER :: iLev, bcType

!******************************************************************************

  global => regions(1)%global
  CALL RegisterFunction( global,'CheckBcValidity',&
  'PREP_CheckBcValidity.F90' )

! start -----------------------------------------------------------------------

  iLev = global%startLevel
  
  DO iReg = 1,global%nRegions
    DO iPatch=1,regions(iReg)%nPatches
      
      patch  => regions(iReg)%levels(iLev)%patches(iPatch)
      bcType =  patch%bcType

      SELECT CASE(bcType)

      CASE (BC_SLIPWALL)
        IF (.NOT. global%prepBcDefined(BC_SLIPWALL)) &
          CALL ErrorStop( global,ERR_NO_BCSPECIFIED,__LINE__, &
            'BC_SLIPW is used in .top but not defined in .bc file' )

      CASE (BC_NOSLIPWALL)
        IF (.NOT. global%prepBcDefined(BC_NOSLIPWALL)) &
          CALL ErrorStop( global,ERR_NO_BCSPECIFIED,__LINE__, &
            'BC_NOSLIP is used in .top but not defined in .bc file' )

      CASE (BC_INFLOW)
        IF (.NOT. global%prepBcDefined(BC_INFLOW) .AND. &
            .NOT. global%prepBcDefined(BC_INFLOW_TOTANG) .AND. &
            .NOT. global%prepBcDefined(BC_INFLOW_VELTEMP) .AND. &
            .NOT. global%prepBcDefined(BC_INFLOW_VELPRESS)) &
          CALL ErrorStop( global,ERR_NO_BCSPECIFIED,__LINE__, &
            'BC_INFLOW_... is used in .top but not defined in .bc file' )

      CASE (BC_OUTFLOW)
        IF (.NOT. global%prepBcDefined(BC_OUTFLOW)) &
          CALL ErrorStop( global,ERR_NO_BCSPECIFIED,__LINE__, &
            'BC_OUTFLOW is used in .top but not defined in .bc file' )
    
      CASE (BC_FARFIELD)
        IF (.NOT. global%prepBcDefined(BC_FARFIELD)) &
          CALL ErrorStop( global,ERR_NO_BCSPECIFIED,__LINE__, &
            'BC_FARF is used in .top but not defined in .bc file' )
    
      CASE (BC_INJECTION)
        IF (.NOT. global%prepBcDefined(BC_INJECTION) .AND. &
            .NOT. global%prepBcDefined(BC_INJECTION_APN)) &
          CALL ErrorStop( global,ERR_NO_BCSPECIFIED,__LINE__, &
            'BC_INJECT is used in .top but not defined in .bc file' )

        IF (global%prepBcDefined(BC_INJECTION_APN) .AND. &
            patch%bcCoupled == BC_EXTERNAL) &
          CALL ErrorStop( global,ERR_EXTERNAL_FUNCT,__LINE__, &
            'BC_INJECT_APN used in .bc with interaction flag in .top file' )

      END SELECT

    ENDDO  ! iReg
  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE CheckBcValidity

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PREP_CheckBcValidity.F90,v $
! Revision 1.7  2008/12/06 08:44:50  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:18:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2006/01/30 20:16:47  wasistho
! added safety for using injectionAPN
!
! Revision 1.4  2006/01/29 09:35:00  wasistho
! added injection_apn
!
! Revision 1.3  2005/04/29 20:20:16  wasistho
! bug fixed, removed second bc_inflow in select case
!
! Revision 1.2  2005/04/29 03:31:33  wasistho
! added distribution bc file generator
!
! Revision 1.1  2004/12/03 02:20:08  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:40:49  wasistho
! lower to upper case
!
! Revision 1.1  2004/07/27 20:29:47  wasistho
! added readBcInputFile and checkBcValidity
!
!
!******************************************************************************











