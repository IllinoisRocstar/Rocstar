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
!          communicating with the current region adjacent to regions with
!          corner or edge cells.
!
! Description: none.
!
! Input: 
!   regions = data of all regions
!   iReg    = index of current region.
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_CECellsClearRequestsData.F90,v 1.3 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CECellsClearRequestsData( regions, iReg )

  USE ModDataTypes
  USE ModError
  USE ModMPI
  USE ModParameters
  USE ModGlobal,     ONLY : t_global
  USE ModDataStruct, ONLY : t_region, t_level, t_dCellTransf

  USE ModPartLag,    ONLY : t_plag

  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)
  
  INTEGER, INTENT(IN)     :: iReg

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
 
#ifdef MPI
  INTEGER :: statusPlag(MPI_STATUS_SIZE)
#endif
  INTEGER :: iLev,ir,iRequestPlag,nBuffSizePlag

  LOGICAL :: doWait

  TYPE(t_global),      POINTER :: global
  TYPE(t_level),       POINTER :: pLevel
  TYPE(t_plag),        POINTER :: pPlag
  TYPE(t_region),      POINTER :: pRegion 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CECellsClearRequestsData.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PLAG_CECellsClearRequestsData',&
  'PLAG_CECellsClearRequestsData.F90' )

#ifdef MPI
  
! ******************************************************************************
! Get dimensions 
! ******************************************************************************

  iLev = regions(iReg)%currLevel

! ******************************************************************************
! Set pointers 
! ******************************************************************************
 
  pRegion => regions(iReg)
  pLevel  => regions(iReg)%levels(iLev)
  pPlag   => pLevel%plag

! ******************************************************************************
!   Wait for edges & corners being received by other processors 
! ******************************************************************************

  DO ir=1,global%nRegions
    IF (regions(iReg)%levels(iLev)%sendEcCells(ir)%nCells > 0) THEN
      iRequestPlag = pLevel%sendEcCells(ir)%iRequestPlag
      nBuffSizePlag = pLevel%sendEcCells(ir)%nBuffSizePlag

      IF ( (regions(ir)%procid /= global%myProcid) .AND. &
           (nBuffSizePlag /= 0) ) THEN
        CALL MPI_Wait( pPlag%requestsCECellsI(iRequestPlag), statusPlag, global%mpierr )
        IF (global%mpierr /= ERR_NONE) & 
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
        
        CALL MPI_Wait( pPlag%requestsCECellsR(iRequestPlag), statusPlag, global%mpierr )
        IF (global%mpierr /= ERR_NONE) & 
          CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )

      ENDIF ! nBuffSizePlag
    ENDIF ! nCells
  ENDDO ! ir
#endif

! ******************************************************************************
! Finalize 
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CECellsClearRequestsData

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CECellsClearRequestsData.F90,v $
! Revision 1.3  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:06  fnajjar
! Initial revision after changing case
!
! Revision 1.1  2004/03/18 21:43:27  fnajjar
! Initial import for MPI-based data buffer communication
!
! Revision 1.1  2004/03/10 23:16:09  fnajjar
! Initial import of routines to MPI-communicate buffer sizes
!
!******************************************************************************







