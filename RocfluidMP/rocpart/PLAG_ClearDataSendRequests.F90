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
! Purpose: wait until data is received by other processors
!          communicating with the current region.
!
! Description: none.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_ClearDataSendRequests.F90,v 1.3 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_ClearDataSendRequests( regions, iReg )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_buffer_plag
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag
  USE ModError
  USE ModMPI
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  
  INTEGER :: iReg

! ... loop variables
  INTEGER :: iPatch

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

#ifdef MPI
  INTEGER :: statusPlag(MPI_STATUS_SIZE)
#endif

  INTEGER :: bcType, iLev, iRegDes, iRequestPlag, nBuffSizeSrc, nPatches

  LOGICAL :: doWait

  TYPE(t_patch), POINTER  :: patch
  TYPE(t_plag),   POINTER :: pPlag
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ClearDataSendRequests.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PLAG_ClearDataSendRequests',&
  'PLAG_ClearDataSendRequests.F90' )

#ifdef MPI

! get dimensions --------------------------------------------------------------

  iLev     = regions(iReg)%currLevel
  nPatches = regions(iReg)%nPatches

  pPlag => regions(iReg)%levels(iLev)%plag 

  DO iPatch = 1, nPatches
    patch => regions(iReg)%levels(iLev)%patches(iPatch)

    bcType  = patch%bcType
    iRegDes = patch%srcRegion
    iRequestPlag = patch%bufferPlag%iRequest

! - extract buffer size

    IF ( (bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
         (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
         (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
         (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
         (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE) ) THEN

        nBuffSizeSrc = patch%bufferPlag%nBuffSize
 
    ENDIF ! bcType
       
! - region interface, periodic boundary ---------------------------------------

    doWait = ((bcType>=BC_REGIONCONF .AND. bcType<=BC_REGIONCONF+BC_RANGE) .OR. &
              (bcType>=BC_REGIONINT  .AND. bcType<=BC_REGIONINT +BC_RANGE) .OR. &
              (bcType>=BC_REGNONCONF .AND. bcType<=BC_REGNONCONF+BC_RANGE) .OR. &
              (bcType>=BC_TRA_PERI   .AND. bcType<=BC_TRA_PERI  +BC_RANGE) .OR. &
              (bcType>=BC_ROT_PERI   .AND. bcType<=BC_ROT_PERI  +BC_RANGE))
              
    IF (iRegDes > 0) THEN
      IF ( doWait .AND. (regions(iRegDes)%procid /= global%myProcid) .AND. &
          (nBuffSizeSrc /= 0) ) THEN
        CALL MPI_Wait( pPlag%requestsI(iRequestPlag), statusPlag, global%mpierr )
        IF (global%mpierr /= ERR_NONE) & 
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
        
        CALL MPI_Wait( pPlag%requestsR(iRequestPlag), statusPlag, global%mpierr )
        IF (global%mpierr /= ERR_NONE) & 
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
      ENDIF ! doWait
    ENDIF ! iRegDes

  ENDDO ! iPatch

#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_ClearDataSendRequests

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ClearDataSendRequests.F90,v $
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:23  fnajjar
! Initial revision after changing case
!
! Revision 1.1  2003/02/21 17:09:54  fnajjar
! Initial import
!
!******************************************************************************







