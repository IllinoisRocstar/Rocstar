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
! Purpose: Check that have coupled boundaries defined for coupled simulations.
!
! Description: None.
!
! Input: 
!   regions     Regions data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_CheckCouplingInput.F90,v 1.4 2008/12/06 08:44:01 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_CheckCouplingInput(regions)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModError
  USE ModParameters
  USE ModMPI

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), DIMENSION(:) :: regions

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN)   :: RCSIdentString
  INTEGER :: cntrGlob,cntrLoc,errorFlag,iPatch,iReg
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_CheckCouplingInput.F90,v $ $Revision: 1.4 $'

  global => regions(1)%global

  CALL RegisterFunction(global,'RFLU_CheckCouplingInput',&
  'RFLU_CheckCouplingInput.F90')

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking coupling input...' 
  END IF ! global%verbLevel

  cntrGlob = 0
  cntrLoc  = 0

! ******************************************************************************
! Count number of interacting patches on this processor
! ******************************************************************************
  
  DO iReg = LBOUND(regions,1),UBOUND(regions,1)
    DO iPatch = 1,regions(iReg)%grid%nPatches
      pPatch => regions(iReg)%patches(iPatch)

      IF ( pPatch%bcCoupled /= BC_NOT_COUPLED ) THEN 
        cntrLoc = cntrLoc + 1
      END IF ! pPatch%bcCoupled
    END DO ! iPatch  
  END DO ! iReg
  
! ******************************************************************************
! Reduce to total number of interacting patches, exit if none
! ******************************************************************************

  CALL MPI_Reduce(cntrLoc,cntrGlob,1,MPI_INTEGER,MPI_SUM,MASTERPROC, & 
                  global%mpiComm,errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_TROUBLE,__LINE__)
  END IF ! global%errorFlag
  
  IF ( global%myProcid == MASTERPROC .AND. cntrGlob == 0 ) THEN 
    CALL ErrorStop(global,ERR_BCCOUPLED_NONE,__LINE__)
  END IF ! global%myProcid

! ******************************************************************************
! End
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Checking coupling input done.' 
  END IF ! global%verbLevel   

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_CheckCouplingInput

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_CheckCouplingInput.F90,v $
! Revision 1.4  2008/12/06 08:44:01  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/04/15 15:05:58  haselbac
! Converted to MPI, cosmetics
!
! Revision 1.1  2003/05/01 14:04:23  haselbac
! Initial revision
!
! ******************************************************************************







