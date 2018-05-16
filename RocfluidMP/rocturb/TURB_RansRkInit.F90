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
! Purpose: set initial solution field to turbulence transport equations, 
!          if any (selected model is of RaNS class), for first RK stage.
!
! Description: none.
!
! Input: region = data of current region,
!        iStage = current RK stage.
!
! Output: region%levels%turb = new solution after one time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_RansRkInit.F90,v 1.6 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_RansRkInit( region, istage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  USE TURB_ModParameters
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset

#include "Indexing.h"
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  
  INTEGER, INTENT(IN) :: istage

! ... loop variables
  INTEGER :: ic, idx

! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: ibc, iec, idxbeg, idxend
  REAL(RFREAL), POINTER :: tcv(:,:), tcvOld(:,:), tdiss(:,:)
  
#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff
#endif

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'TURB_RansRkInit',&
  'TURB_RansRkInit.F90' )

  IF (region%turbInput%modelClass /= MODEL_RANS) GOTO 999
  
! get dimensions and pointers -------------------------------------------------

#ifdef RFLO
  iLev = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  tcv     => region%levels(iLev)%turb%cv
  tcvOld  => region%levels(iLev)%turb%cvOld
  tdiss   => region%levels(iLev)%turb%diss
#endif

#ifdef RFLU
  ibc = 1
  iec = region%grid%nCellsTot        

  tcv     => region%turb%cv
  tcvOld  => region%turb%cvOld
  tdiss   => region%turb%diss
#endif

! select start and end index of 1st dimension depending on RaNS model selected

  IF ((region%mixtInput%turbModel == TURB_MODEL_SA) .OR. &
      (region%mixtInput%turbModel == TURB_MODEL_DESSA) .OR. &
      (region%mixtInput%turbModel == TURB_MODEL_HDESSA)) THEN
    idxbeg = CV_SA_NUTIL
    idxend = CV_SA_NUTIL
  ENDIF

! store previous solution; set dissipation to zero ----------------------------
  
  IF (iStage == 1) THEN
    DO ic=ibc,iec
      DO idx=idxbeg, idxend
        tcvOld(idx,ic) = tcv(idx,ic)
        tdiss(idx,ic)  = 0._RFREAL
      ENDDO
    ENDDO
  ELSE
    DO ic=ibc,iec
      DO idx=idxbeg, idxend
        tdiss(idx,ic) = 0._RFREAL
      ENDDO
    ENDDO
  ENDIF ! iStage

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_RansRkInit

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_RansRkInit.F90,v $
! Revision 1.6  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/03/09 06:35:47  wasistho
! incorporated HDESSA
!
! Revision 1.3  2004/07/03 02:07:35  wasistho
! delete ! public
!
! Revision 1.2  2004/03/20 00:29:17  wasistho
! set turb_rflo_ransNumericalDiss to turb_ransNumerical..
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.1  2003/10/07 02:17:03  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







