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
! Purpose: set initial solution field to FLD transport equation, 
!          if any (selected model is FLDTRAN), for first RK stage.
!
! Description: none.
!
! Input: region = data of current region,
!        iStage = current RK stage.
!
! Output: region%levels%radi = new solution after one time step.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_FlimRkInit.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FlimRkInit( region, istage )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModMPI
  USE ModParameters
  USE RADI_ModParameters
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
  REAL(RFREAL), POINTER :: cv(:,:), cvOld(:,:), diss(:,:)
  
#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff
#endif

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RADI_FlimRkInit',&
  'RADI_FlimRkInit.F90' )

  IF (region%radiInput%radiModel /= RADI_MODEL_FLDTRAN) GOTO 999
  
! get dimensions and pointers -------------------------------------------------

#ifdef RFLO
  iLev = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  cv     => region%levels(iLev)%radi%cv
  cvOld  => region%levels(iLev)%radi%cvOld
  diss   => region%levels(iLev)%radi%diss
#endif

#ifdef RFLU
  ibc = 1
  iec = region%grid%nCellsTot        

  cv     => region%radi%cv
  cvOld  => region%radi%cvOld
  diss   => region%radi%diss
#endif

! select start and end index of 1st dimension depending on RADI model selected

  IF (region%radiInput%radiModel == RADI_MODEL_FLDTRAN) THEN
    idxbeg = CV_RADI_ENER
    idxend = CV_RADI_ENER
  ENDIF

! store previous solution; set dissipation to zero ----------------------------
  
  IF (iStage == 1) THEN
    DO ic=ibc,iec
      DO idx=idxbeg, idxend
        cvOld(idx,ic) = cv(idx,ic)
        diss(idx,ic)  = 0._RFREAL
      ENDDO
    ENDDO
  ELSE
    DO ic=ibc,iec
      DO idx=idxbeg, idxend
        diss(idx,ic) = 0._RFREAL
      ENDDO
    ENDDO
  ENDIF ! iStage

! finalize --------------------------------------------------------------------

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_FlimRkInit

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_FlimRkInit.F90,v $
! Revision 1.3  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/09/30 17:49:10  wasistho
! prepared for full FLD radiation model
!
!
!
!******************************************************************************







