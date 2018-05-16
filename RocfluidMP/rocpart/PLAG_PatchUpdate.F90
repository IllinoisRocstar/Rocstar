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
! Purpose: update step for patch data buffers.
!
! Description: none.
!
! Input: region = current region
!        iReg   = current region number
!        
!
! Output: region%plag = plag variables
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_PatchUpdate.F90,v 1.3 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_patchUpdate( regions, iReg )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE PLAG_ModInterfaces, ONLY: PLAG_BoundaryConditionsSet,  &
                                PLAG_PatchGetBufferSize,     &
                                PLAG_PatchLoadDataBuffers                     
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters         
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  
  INTEGER :: iReg
  
! ... loop variables
  
! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
 
  TYPE(t_global), POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_PatchUpdate.F90,v $ $Revision: 1.3 $'
  
  global => regions(iReg)%global
  
  CALL RegisterFunction( global, 'PLAG_patchUpdate',&
  'PLAG_PatchUpdate.F90' )

! Obtain patch buffer size ----------------------------------------------------

!  WRITE(STDOUT,'(A,I2)') '    Entering PLAG_PatchGetBufferSize: iReg', iReg

  CALL PLAG_PatchGetBufferSize( regions(iReg) )
  
! Load data buffer and reshuffle particle datastructure -----------------------
  
!  WRITE(STDOUT,'(A,I2)') '    Entering PLAG_PatchLoadDataBuffers: iReg', iReg

  CALL PLAG_PatchLoadDataBuffers( regions(iReg), iReg )

! Impose boundary conditions --------------------------------------------------

!  WRITE(STDOUT,'(A,I2)') '    Entering PLAG_BoundaryConditionsSet: iReg', iReg

  CALL PLAG_BoundaryConditionsSet( regions, iReg )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_patchUpdate

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_PatchUpdate.F90,v $
! Revision 1.3  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:59  fnajjar
! Initial revision after changing case
!
! Revision 1.6  2004/02/26 21:02:17  haselbac
! Removed iStage argument
!
! Revision 1.5  2003/01/23 17:03:22  f-najjar
! Moved calls to PLAG_AppendDataFromBuffers into PLAG_PatchBufferSendRecv and removed call to PLAG_PatchBufferSendRecv
!
! Revision 1.4  2003/01/17 19:33:28  f-najjar
! Included iReg in call sequence for PLAG_PatchLoadDataBuffers
!
! Revision 1.3  2003/01/16 22:33:00  f-najjar
! Added Calling for PLAG_AppendDatafromBuffers and PLAG_PatchBufferSendRecv
!
! Revision 1.2  2003/01/13 19:03:24  f-najjar
! Removed PLAG_allocateDataBuffers and PLAG_deallocateDataBuffers from calling sequence
!
! Revision 1.1  2002/10/25 14:18:35  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************







