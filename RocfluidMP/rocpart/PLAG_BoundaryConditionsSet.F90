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
! Purpose: set boundary conditions or exchange data between adjacent
!          regions being on same processor.
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: regions(iReg)%levels%plag = updated values (aiv,arv,cv,dv,tv)
!                                     in destination region.
!
! Notes: Current support for conforming region interface only.
!
!******************************************************************************
!
! $Id: PLAG_BoundaryConditionsSet.F90,v 1.3 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_BoundaryConditionsSet( regions,iReg )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE PLAG_ModInterfaces, ONLY : PLAG_PatchExchangeConf
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER :: iReg

! ... loop variables
  INTEGER :: iPatch

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iLev, nPatches, bcType, iRegDes, iPatchDes

  TYPE(t_patch), POINTER  :: patch, patchDes
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_BoundaryConditionsSet.F90,v $ $Revision: 1.3 $'

  global => regions(iReg)%global

  CALL RegisterFunction( global,'PLAG_BoundaryConditionsSet',&
  'PLAG_BoundaryConditionsSet.F90' )

! get dimensions --------------------------------------------------------------

  iLev     = regions(iReg)%currLevel
  nPatches = regions(iReg)%nPatches

! loop over patches -----------------------------------------------------------

  DO iPatch=1,nPatches

    patch => regions(iReg)%levels(iLev)%patches(iPatch)

    bcType    = patch%bcType
    iRegDes   = patch%srcRegion
    iPatchDes = patch%srcPatch

! - conforming region interface

    IF (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) THEN
      patchDes => regions(iRegDes)%levels(iLev)%patches(iPatchDes)

      IF (regions(iRegDes)%procid == global%myProcid) THEN
        CALL PLAG_PatchExchangeConf( regions(iReg),regions(iRegDes), &
                                     patch,patchDes, iReg, iRegDes )
      ENDIF

! - non-conforming region interface (integer)

    ELSE IF ( (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT+BC_RANGE)  .OR. &
              (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
              (bcType>=BC_TRA_PERI .AND. bcType<=BC_TRA_PERI+BC_RANGE)     .OR. &
              (bcType>=BC_ROT_PERI .AND. bcType<=BC_ROT_PERI+BC_RANGE)   ) THEN
      patchDes => regions(iRegDes)%levels(iLev)%patches(iPatchDes)
      
      CALL ErrorStop( global,ERR_UNKNOWN_BC,__LINE__ )

    ENDIF  ! bcType

  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_BoundaryConditionsSet

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_BoundaryConditionsSet.F90,v $
! Revision 1.3  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:56:54  fnajjar
! Initial revision after changing case
!
! Revision 1.2  2003/01/13 19:29:35  f-najjar
! Included iReg and iRegDes in calling sequence for PLAG_patchExchangeConf
!
! Revision 1.1  2003/01/13 19:14:50  f-najjar
! Initial Import of Multiblock capability
!
!******************************************************************************







