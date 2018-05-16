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
! Purpose: wrapper for update step in PLAG module.
!
! Description: none.
!
! Input: pRegion = data of current region,
!        iReg    = current region,
!        istage  = current RK stage
!
! Output: regions%levels%mixt = new solution after one time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_RkUpdateWrapper.F90,v 1.3 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RkUpdateWrapper( region, iReg, istage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PLAG_ModInterfaces, ONLY: PLAG_InjcTileUpdate

#ifdef RFLO  
  USE PLAG_ModInterfaces, ONLY: PLAG_Update
#endif

#ifdef RFLU
  USE PLAG_ModInterfaces, ONLY: PLAG_RFLU_Update
#endif    

  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  
  INTEGER, INTENT(IN) :: iReg, istage

! ... loop variables

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_global), POINTER :: global

#ifdef RFLU
  TYPE(t_region), POINTER :: pRegion
#endif 

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RkUpdateWrapper.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'PLAG_rkUpdateWrapper',&
  'PLAG_RkUpdateWrapper.F90' )

! update tiles ----------------------------------------------------------------

  CALL PLAG_InjcTileUpdate( region, iReg, iStage )

! update evolution equations --------------------------------------------------

#ifdef RFLO
!  WRITE(*,*) 'Entering PLAG_Update: iReg, iStage = ',iReg, iStage
  CALL PLAG_Update( region, iReg, iStage )
#endif

#ifdef RFLU
!  WRITE(*,*) 'Entering PLAG_RFLU_Update: iReg, iStage = ',iReg, iStage
  pRegion => region%pRegion
  CALL PLAG_RFLU_Update( pRegion, iStage )
#endif

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RkUpdateWrapper

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RkUpdateWrapper.F90,v $
! Revision 1.3  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:58:17  fnajjar
! Initial revision after changing case
!
! Revision 1.4  2004/03/26 21:32:36  fnajjar
! Cleaned up routine for separate RFLO and RFLU calls within ifdef constructs and added PLAG_RFLU_Update call
!
! Revision 1.3  2004/03/08 22:26:29  fnajjar
! Removed ifdef RFLO around PLAG_InjcTileUpdate since injection is active with RFLU
!
! Revision 1.2  2004/02/26 21:02:22  haselbac
! Commented out RFLO-specific tile update routines
!
! Revision 1.1  2003/03/28 19:53:10  fnajjar
! Initial import of wrapper routines for RocfluidMP
!
!******************************************************************************







