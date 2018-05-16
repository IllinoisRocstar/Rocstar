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
! Purpose: wrapper for corner and edge Cells update step in PLAG module.
!
! Description: none.
!
! Input: regions = data of all regions,
!
! Output: regions%levels%plag = solution in corner and edge cells
!                               after one RK stage.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_CECellsWrapper.F90,v 1.3 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CECellsWrapper( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PLAG_ModInterfaces, ONLY : PLAG_CECellsExchange,        &
                                 PLAG_CECellsGetBufferSize,   &
                                 PLAG_CECellsLoadDataWrapper, &
                                 PLAG_CECellsSendRecvWrapper 

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CECellsWrapper.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PLAG_CECellsWrapper',&
  'PLAG_CECellsWrapper.F90' )

! check if module is active in any region======================================

  IF (.NOT. global%plagUsed) GOTO 999

! Load buffer data for on-processor regions -----------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE ) THEN             ! on my processor

#ifdef PLAG_CECELLS_DEBUG
      WRITE(*,*) '----------------------------------------------------'
      WRITE(*,*) 'Entering PLAG_CECellsGetBufferSize  : iReg=',iReg
#endif
      CALL PLAG_CECellsGetBufferSize(   regions(iReg), iReg )

#ifdef PLAG_CECELLS_DEBUG
      WRITE(*,*) '----------------------------------------------------'
      WRITE(*,*) 'Entering PLAG_CECellsLoadDataWrapper: iReg=',iReg
#endif
      CALL PLAG_CECellsLoadDataWrapper( regions      , iReg )
    ENDIF ! regions

  ENDDO ! iReg

! Exchange buffer for on-processor regions ------------------------------------

  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE ) THEN             ! on my processor

#ifdef PLAG_CECELLS_DEBUG
      WRITE(*,*) '----------------------------------------------------'
      WRITE(*,*) 'Entering PLAG_CECellsExchange: iReg=',iReg
#endif
      CALL PLAG_CECellsExchange( regions, iReg )

    ENDIF ! regions

  ENDDO ! iReg

! ******************************************************************************
! Synchronize through an MPI barrier 
! ******************************************************************************

#ifdef MPI
  CALL MPI_Barrier( global%mpiComm,global%mpierr )
  IF ( global%mpierr /= ERR_NONE ) &
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
#endif

! Communicate buffer data for off-processor regions ---------------------------

  CALL PLAG_CECellsSendRecvWrapper( regions )

! finalize ====================================================================

999 CONTINUE
  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CECellsWrapper

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CECellsWrapper.F90,v $
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:21  fnajjar
! Initial revision after changing case
!
! Revision 1.4  2004/03/10 23:14:04  fnajjar
! Cleaned dead sections and included call for MPI-based corner-edge cells wrapper
!
! Revision 1.3  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.2  2004/01/26 22:54:43  fnajjar
! Included new wrapper for load data, PLAG_cECellsLoadDataWrapper
!
! Revision 1.1  2003/11/12 21:37:59  fnajjar
! Initial import of Corner-Edge cells Infrastructure
!
!******************************************************************************







