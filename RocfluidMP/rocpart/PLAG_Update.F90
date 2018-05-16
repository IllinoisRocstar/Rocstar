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
! Purpose: update step for PLAG module.
!
! Description: none.
!
! Input: region = current region
!        iReg   = current region number
!        iStage = current RK stage.
!
! Output: region%plag = plag variables
!
! Notes: This corresponds to Part IV Steps 10d and Part V 
!        in RocfluidMP framework.
!
!******************************************************************************
!
! $Id: PLAG_Update.F90,v 1.4 2009/10/26 00:19:32 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_Update( region, iReg, iStage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt_input
  USE ModPartLag, ONLY    : t_plag
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters         
  
  USE ModInterfaces, ONLY: RkUpdateGeneric
  
  USE PLAG_ModInterfaces, ONLY: PLAG_CalcBreakup,            &
                                PLAG_CalcRhsPosition,        &
                                PLAG_GetCellIndices,         &
                                PLAG_GetCellIndicesOutflow,  &
                                PLAG_InjcEjectParticle,      &
                                PLAG_PatchRemoveDataOutflow, &
                                PLAG_UpdateDataStruct,       &
                                PLAG_WallBounce

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region
  
  INTEGER, INTENT(IN) :: iReg, iStage
  
! ... loop variables

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: iLev

  TYPE(t_plag),   POINTER :: pPlag  
  TYPE(t_global), POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_Update.F90,v $ $Revision: 1.4 $'
  
  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_Update',&
  'PLAG_Update.F90' )
  
! get dimensions and set pointer ----------------------------------------------

  iLev  =  region%currLevel
  pPlag => region%levels(iLev)%plag

! Main algorithm for discrete particle evolution ------------------------------
!  Active for non-zero number of particles in region --------------------------

  IF ( pPlag%nPcls > 0 ) THEN

! - Calculate RHS for position vector -----------------------------------------

!    WRITE(STDOUT,'(A)') '    Entering PLAG_CalcRhsPosition', iReg, iStage
    CALL PLAG_CalcRhsPosition( region )

! - Invoke RK Update ----------------------------------------------------------

!    WRITE(STDOUT,'(A)') '    Entering PLAG_RKUpdate', iReg, iStage
    CALL RkUpdateGeneric(region,VAR_TYPE_POINT,iStage,1,pPlag%nPcls,1, & 
                         pPlag%nCv,pPlag%cv,pPlag%cvOld,pPlag%rhs,pPlag%rhsSum)

! - Apply wall bouncing algorithm ---------------------------------------------

!    WRITE(STDOUT,'(A)') '    Entering PLAG_WallBounce'
    CALL PLAG_WallBounce( region )

! - Invoke particle location search for outflow -------------------------------

!    WRITE(STDOUT,'(A)') '    Entering PLAG_GetCellIndicesOutflow'
    CALL PLAG_GetCellIndicesOutflow( region )

! - Remove exiting particle from datastructure --------------------------------

!    WRITE(STDOUT,'(A)') '    Entering PLAG_PatchRemoveDataOutflow'
    CALL PLAG_PatchRemoveDataOutflow( region, iReg )

! - Invoke particle location search -------------------------------------------

!    WRITE(STDOUT,'(A)') '    Entering PLAG_GetCellIndices'
    CALL PLAG_GetCellIndices( region, iReg )

! - Remove particle with a delete status from datastructure -------------------
    CALL PLAG_UpdateDataStruct( region,iReg )

  END IF ! nPcls

! Invoke injection algorithm at last RK stage ----------------------------------

  IF ( iStage == global%nrkSteps ) THEN 

!    WRITE(STDOUT,'(A)') '  Entering PLAG_InjcEjectParticle'
    CALL PLAG_InjcEjectParticle( region, iReg )

!    WRITE(STDOUT,'(A,I4,I8)') '     iReg: nPcls = ',iReg, pPlag%nPcls

  END IF ! iStage

! Invoke breakup algorithm at last RK stage -----------------------------------

  IF ( region%plagInput%breakupModel > PLAG_BREAKUP_NOMODEL .AND. &
       iStage == global%nrkSteps                                  ) THEN 
    CALL PLAG_CalcBreakup( region, iReg )
  END IF ! iStage

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_Update

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_Update.F90,v $
! Revision 1.4  2009/10/26 00:19:32  mtcampbe
! Updates for completion of NATIVE_MP_IO
!
! Revision 1.3  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:58:18  fnajjar
! Initial revision after changing case
!
! Revision 1.20  2004/11/17 16:44:11  haselbac
! Replaced RKUpdate call with call to RkUpdateGeneric
!
! Revision 1.19  2004/08/05 20:09:54  fnajjar
! Reactivate PLAG_UpdateDataStruct
!
! Revision 1.18  2004/04/09 23:16:36  fnajjar
! Added call to PLAG_UpdateDataStruct for appropriate deletion kernel
!
! Revision 1.17  2004/03/26 21:30:54  fnajjar
! Removed RFLU-specific calls since PLAG_RFLU_Update is created
!
! Revision 1.16  2004/03/22 23:48:30  fnajjar
! Activated RFLU-aware trajectory-based search algorithm
!
! Revision 1.15  2004/03/15 21:09:04  haselbac
! Changed name of particle location routine
!
! Revision 1.14  2004/03/08 22:34:36  fnajjar
! Moved pRegion section and added PLAG_RFLU_InjectionDriver
!
! Revision 1.13  2004/03/05 23:21:00  haselbac
! Added (as yet commented-out) call to PLAG_RFLU_FindParticleCellsTraj
!
! Revision 1.12  2004/02/26 21:02:23  haselbac
! Added initial RFLU support
!
! Revision 1.11  2004/02/06 21:18:20  fnajjar
! Initial Integration of Rocpart with Rocflu
!
! Revision 1.10  2003/09/13 20:14:22  fnajjar
! Added infrastructure for Breakup model
!
! Revision 1.9  2003/05/28 15:16:11  fnajjar
! Removed obsolete PLAG_mixt calls as embedded in Rocinteract
!
! Revision 1.8  2003/05/07 00:02:07  fnajjar
! Included I/O for nPcls
!
! Revision 1.7  2003/04/15 23:02:38  fnajjar
! Removed dead code section
!
! Revision 1.6  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
! Revision 1.5  2003/02/04 19:07:33  f-najjar
! Aligned calls in routine with RocfluidMP framework
!
! Revision 1.4  2003/01/17 19:31:15  f-najjar
! Included iReg in calling sequence for PLAG_PatchRemoveDataOutflow
!
! Revision 1.3  2003/01/16 20:44:21  f-najjar
! Include iReg in calling sequence of PLAG_getCellIndices
!
! Revision 1.2  2003/01/10 19:09:10  f-najjar
! Included iReg in calling sequence
!
! Revision 1.1  2002/10/25 14:20:32  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************







