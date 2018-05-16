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
! Purpose: wrapper routine to load data buffer size for corner and edge cells
!          and shrinks particle datastructure. 
!
! Description: none.
!
! Input: regions = data of all regions,
!        iReg    = current region number.
!
! Output: region%level%edgeCells%buffExchPlag%aiv,arv,cv,dv,tv   = buffer data.
!         region%level%cornerCells%buffExchPlag%aiv,arv,cv,dv,tv = buffer data.
!         region%level%plag%aiv,arv,cv,dv,tv = Plag data.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: PLAG_CECellsLoadDataWrapper.F90,v 1.3 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CECellsLoadDataWrapper( regions, iReg )

  USE ModDataTypes
  USE ModPartLag, ONLY    : t_plag
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global  
  USE ModError
  USE ModParameters
  USE PLAG_ModInterfaces, ONLY : PLAG_CornCellsLoadData, &
                                 PLAG_EdgeCellsLoadData
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

  INTEGER, INTENT(IN)     :: iReg

! ... loop variables

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  
  INTEGER :: iLev, nPcls

  TYPE(t_global),      POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PLAG_CECellsLoadDataWrapper.F90,v $ $Revision: 1.3 $'

  global => regions(iReg)%global
    
  CALL RegisterFunction( global, 'PLAG_CECellsLoadDataWrapper',&
  'PLAG_CECellsLoadDataWrapper.F90' )

! Get dimensions --------------------------------------------------------------

  iLev  = regions(iReg)%currLevel
  nPcls = regions(iReg)%levels(iLev)%plag%nPcls 

! Exit for null number of particles -------------------------------------------

  IF (nPcls == 0) GOTO 999 

! Load data buffers for corner cells ------------------------------------------

#ifdef PLAG_CECELLS_DEBUG
  WRITE(*,*) ' Entering PLAG_CornCellsLoadData: pid, iReg = ',global%myProcId, iReg
#endif
!
  CALL PLAG_CornCellsLoadData( regions, iReg )
  
! Load data buffers for edge cells ------------------------------------------

#ifdef PLAG_CECELLS_DEBUG
  WRITE(*,*) ' Entering PLAG_EdgeCellsLoadData: pid, iReg = ',global%myProcId, iReg
#endif

  CALL PLAG_EdgeCellsLoadData( regions, iReg )

! finalize --------------------------------------------------------------------

999  CONTINUE
  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CECellsLoadDataWrapper

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CECellsLoadDataWrapper.F90,v $
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:14  fnajjar
! Initial revision after changing case
!
! Revision 1.1  2004/01/26 22:56:28  fnajjar
! Initial import for corner-edge cells to load buffer data
!
!******************************************************************************







