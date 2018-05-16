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
! Purpose: driver to communicate buffer data for adjacent regions 
!          on different processors.
!
! Description: none.
!
! Input: regions = data of all regions.
!
! Output: buffer sizes, data and appended PLAG data.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_PatchBufferSendRecv.F90,v 1.3 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_PatchBufferSendRecv( regions )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_buffer_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag
  USE PLAG_ModInterfaces, ONLY: PLAG_AppendDatafromBuffers, &
                                PLAG_BufferDataRecv,        &
                                PLAG_BufferDataSend,        &                                
                                PLAG_BufferSizeRecv,        &
                                PLAG_BufferSizeSend,        &
                                PLAG_ClearSizeSendRequests, &
                                PLAG_ClearDataSendRequests
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_PatchBufferSendRecv.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PLAG_PatchBufferSendRecv',&
  'PLAG_PatchBufferSendRecv.F90' )

! communicate buffer size -----------------------------------------------------
                           
! send buffer size ------------------------------------------------------------  

  DO iReg = 1, global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor        
      CALL PLAG_BufferSizeSend(regions, iReg)    
    ENDIF       ! regions    
  ENDDO         ! iReg
                     
! receive buffer size ---------------------------------------------------------

  DO iReg = 1, global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor        
      CALL PLAG_BufferSizeRecv(regions, iReg)
    ENDIF ! regions
  ENDDO ! iReg 

! wait for data being received by other processors ----------------------------

  DO iReg = 1, global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor        
      CALL PLAG_ClearSizeSendRequests(regions, iReg)
    ENDIF ! regions  
  ENDDO ! iReg

! communicate buffer arrays ---------------------------------------------------
                           
! send buffer data ------------------------------------------------------------  

  DO iReg = 1, global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor        
      CALL PLAG_BufferDataSend(regions, iReg)    
    ENDIF       ! regions    
  ENDDO         ! iReg
                     
! receive buffer data ---------------------------------------------------------

  DO iReg = 1, global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor        
      CALL PLAG_BufferDataRecv(regions, iReg)
    ENDIF ! regions
  ENDDO ! iReg 

! wait for data being received by other processors ----------------------------

  DO iReg = 1, global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor        
      CALL PLAG_ClearDataSendRequests(regions, iReg)
    ENDIF ! regions  
  ENDDO ! iReg
  
! Append data from buffers ----------------------------------------------------

  DO iReg = 1, global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor
  
!     WRITE(STDOUT,'(A,I2)') '    Entering PLAG_AppendDatafromBuffers: iReg', iReg

      CALL PLAG_AppendDatafromBuffers( regions(iReg), iReg )

    ENDIF ! regions  
  ENDDO ! iReg
  
#ifdef MPI
  CALL MPI_Barrier( global%mpiComm,global%mpierr )
  IF ( global%mpierr /= ERR_NONE ) &
    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

!  CALL MPI_Finalize(global%mpierr)
!  IF ( global%mpierr /= ERR_NONE ) &
!    CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
!
!  STOP
#endif
  
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_PatchBufferSendRecv

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_PatchBufferSendRecv.F90,v $
! Revision 1.3  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:54  fnajjar
! Initial revision after changing case
!
! Revision 1.4  2003/02/21 17:07:55  fnajjar
! Included calls to Data Send and Recv
!
! Revision 1.3  2003/01/24 22:40:27  f-najjar
! Added MPI_Barrier and MPI_Finalize for testing
!
! Revision 1.2  2003/01/24 22:08:31  f-najjar
! Include function call to PLAG_ClearRequests
!
! Revision 1.1  2003/01/23 18:51:49  f-najjar
! Initial Import for MPI
!
!******************************************************************************







