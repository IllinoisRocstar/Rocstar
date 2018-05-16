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
! Purpose: set initial solution field to radiant Energy transport equation, 
!          if any (selected model is FLDTRAN), for first EMS stage,
!          and initiate dissipation at stages > 1.
!
! Description: none.
!
! Input: region = data of current region,
!        iStage = current EMS stage.
!
! Output: region%levels%radi = initial values
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_FlimEmsInit.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FlimEmsInit( region, istage )

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

  INTEGER :: ibc, iec, idxbeg, idxend, ldiss(5)
  REAL(RFREAL) :: blend1, betrk(5)
  REAL(RFREAL), POINTER :: rcv(:,:), rcvOld(:,:), rdiss(:,:)
  
#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff
#endif

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RADI_FlimEmsInit',&
  'RADI_FlimEmsInit.F90' )
  
! get parameters, dimensions and pointers -------------------------------------

  ldiss(:) = region%mixtInput%ldiss(:)
  betrk(:) = region%mixtInput%betrk(:)

#ifdef RFLO
  iLev = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  rcv     => region%levels(iLev)%radi%cv
  rcvOld  => region%levels(iLev)%radi%cvOld
  rdiss   => region%levels(iLev)%radi%diss
#endif

#ifdef RFLU
  ibc = 1
  iec = region%grid%nCellsTot        

  rcv     => region%radi%cv
  rcvOld  => region%radi%cvOld
  rdiss   => region%radi%diss
#endif

! select start and end index of 1st dimension depending on RADI model selected

  IF (region%radiInput%radiModel == RADI_MODEL_FLDTRAN) THEN
    idxbeg = CV_RADI_ENER
    idxend = CV_RADI_ENER
  ENDIF

! store previous solution and initialize dissipation --------------------------
  
  IF (iStage == 1) THEN
    DO ic=ibc,iec
      DO idx=idxbeg, idxend
        rcvOld(idx,ic) = rcv(idx,ic)
        rdiss(idx,ic)  = 0._RFREAL
      ENDDO
    ENDDO
  ENDIF
  IF (iStage>1 .AND. ldiss(iStage)/=0) THEN
    blend1 = 1._RFREAL - betrk(iStage)
    DO ic=ibc,iec
      DO idx=idxbeg, idxend
        rdiss(idx,ic) = blend1*rdiss(idx,ic)
      ENDDO
    ENDDO
  ENDIF ! iStage

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_FlimEmsInit

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_FlimEmsInit.F90,v $
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







