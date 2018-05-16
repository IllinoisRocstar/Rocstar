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
! Purpose: Compute the L2-norm of the density changes.
!
! Description: None.
!
! Input: 
!   regions     Array of regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ResidualNorm.F90,v 1.8 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_ResidualNorm(regions)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,icg,iReg
  REAL(RFREAL) :: dr,drSum,drSumTot
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCv,pCvOld
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ResidualNorm.F90,v $ $Revision: 1.8 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_ResidualNorm',&
  'RFLU_ResidualNorm.F90')

! ******************************************************************************
! Sum density changes over local regions 
! ******************************************************************************

  drSum = 0.0_RFREAL

  DO iReg = 1,global%nRegionsLocal
    pCv    => regions(iReg)%mixt%cv
    pCvOld => regions(iReg)%mixt%cvOld

    DO icg = 1,regions(iReg)%grid%nCells
      dr    = pCv(CV_MIXT_DENS,icg) - pCvOld(CV_MIXT_DENS,icg)
      drSum = drSum + dr*dr
    END DO ! icg
  END DO ! iReg

! ******************************************************************************
! Reduce across processes 
! ******************************************************************************

  IF ( global%nRegions > 1 ) THEN 
    CALL MPI_AllReduce(drSum,drSumTot,1,MPI_RFREAL,MPI_SUM,global%mpiComm, &
                       errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN 
      CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
    END IF ! global%error
  ELSE
    drSumTot = drSum
  END IF ! global%nRegions

! ******************************************************************************
! Finalize
! ******************************************************************************

  global%residual = SQRT(drSumTot)

  IF ( global%currentIter == 1 ) THEN
    global%resInit = global%residual
  END IF ! global

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_ResidualNorm

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ResidualNorm.F90,v $
! Revision 1.8  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:43  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.5  2005/04/15 15:07:24  haselbac
! Converted to MPI
!
! Revision 1.4  2002/09/09 15:51:56  haselbac
! global now under region
!
! Revision 1.3  2002/07/25 14:30:40  haselbac
! Added FEM call to find proper residual norm for parallel runs
!
! Revision 1.2  2002/06/14 20:19:59  haselbac
! Deleted ModLocal, changed local%nRegions to global%nRegionsLocal
!
! Revision 1.1  2002/05/04 17:01:59  haselbac
! Initial revision
!
! ******************************************************************************







