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
! Purpose: update FLD radiation solution for unsteady flow using RK scheme
!
! Description: none.
!
! Input: region = data of current region,
!
! Output: region%levels%radi%cv = new solution after each stage.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RADI_FlimRkUpdate.F90,v 1.3 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_FlimRkUpdate( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset, &
                            RFLO_ExtrapIntCellVec

#include "Indexing.h"
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: ic, idx
  
! ... local variables
#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev, iCOff, ijCOff
#endif
  INTEGER :: iStage, ibc, iec, idxb, idxe

  LOGICAL :: moveGrid

  REAL(RFREAL), POINTER :: cv(:,:), cvOld(:,:), rhs(:,:), rhsSum(:,:)
  REAL(RFREAL), POINTER :: vol(:), volOld(:)
  REAL(RFREAL)          :: ark(5), grk(5), fac, adtv, volRat
  
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RADI_FlimRkUpdate',&
  'RADI_FlimRkUpdate.F90' )

! get stage number ------------------------------------------------------------

  iStage = region%irkStep

! get dimensions and pointers -------------------------------------------------

  ark(:)   = region%mixtInput%ark(:)
  grk(:)   = region%mixtInput%grk(:)

  moveGrid = region%mixtInput%moveGrid

  IF (region%radiInput%radiModel == RADI_MODEL_FLDTRAN) THEN
    idxb = CV_RADI_ENER
    idxe = CV_RADI_ENER
  ENDIF

#ifdef RFLO
  iLev  = region%currLevel

  CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
  ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
  iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

  cv     => region%levels(iLev)%radi%cv
  cvOld  => region%levels(iLev)%radi%cvOld
  rhs    => region%levels(iLev)%radi%rhs
  rhsSum => region%levels(iLev)%radi%rhsSum
  vol    => region%levels(iLev)%grid%vol
  IF (moveGrid) volOld => region%levels(iLev)%gridOld%vol
#endif

#ifdef RFLU
  ibc = 1
  iec = region%grid%nCellsTot        

  cv     => region%radi%cv
  cvOld  => region%radi%cvOld
  rhs    => region%radi%rhs
  rhsSum => region%radi%rhsSum
  vol    => region%grid%vol
  IF (moveGrid) volOld => region%gridOld%vol
#endif

! update FLD radiation solution, sum up residuals -----------------------------

  fac = ark(istage)*global%dtMin
    
  IF (moveGrid) THEN

! - grid is moving (but only at every time step, NOT stage)

    IF (istage == 1) THEN
      DO ic=ibc,iec
        adtv   = fac/vol(ic)
        volRat = volOld(ic)/vol(ic)
        DO idx = idxb,idxe
          cv(idx,ic)     = volRat*cvOld(idx,ic) - adtv*rhs(idx,ic)
          rhsSum(idx,ic) = rhs(idx,ic)
        ENDDO
      ENDDO ! ic
      
    ELSE IF (istage == global%nrkSteps) THEN
      DO ic=ibc,iec
        adtv   = fac/vol(ic)
        volRat = volOld(ic)/vol(ic)                                        
        DO idx = idxb,idxe
          cv(idx,ic) = volRat*cvOld(idx,ic) - adtv*(rhs(idx,ic)+ &
                                                    rhsSum(idx,ic))
        ENDDO
      ENDDO ! ic
      
    ELSE
      DO ic=ibc,iec
        adtv   = fac/vol(ic)
        volRat = volOld(ic)/vol(ic)
        DO idx = idxb,idxe
          cv(idx,ic)     = volRat*cvOld(idx,ic) - adtv*rhs(idx,ic)
          rhsSum(idx,ic) = rhsSum(idx,ic) + grk(istage)*rhs(idx,ic)
        ENDDO
      ENDDO ! ic
    ENDIF   ! istage

  ELSE      ! moveGrid

! - grid is fixed

    IF (istage == 1) THEN
      DO ic=ibc,iec
        adtv = fac/vol(ic)                                
        DO idx = idxb,idxe
          cv(idx,ic)     = cvOld(idx,ic) - adtv*rhs(idx,ic)
          rhsSum(idx,ic) = rhs(idx,ic)
        ENDDO
      ENDDO ! ic
      
    ELSE IF (istage == global%nrkSteps) THEN
      DO ic=ibc,iec
        adtv = fac/vol(ic)
        DO idx = idxb,idxe
          cv(idx,ic) = cvOld(idx,ic) - adtv*(rhs(idx,ic)+ &
                                             rhsSum(idx,ic))
        ENDDO
      ENDDO !ic
      
    ELSE
      DO ic=ibc,iec
        adtv = fac/vol(ic)              
        DO idx = idxb,idxe
          cv(idx,ic)     = cvOld(idx,ic) - adtv*rhs(idx,ic)
          rhsSum(idx,ic) = rhsSum(idx,ic) + grk(istage)*rhs(idx,ic)
        ENDDO
      ENDDO !ic
      
    ENDIF   ! istage
  ENDIF     ! moveGrid

! extrapolate solution to dummy cells

#ifdef RFLO
  CALL RFLO_ExtrapIntCellVec( region,idxb,idxe,cv )
#endif
#ifdef RFLU
!  CALL RFLU_ExtrapIntCellVec( region,idxb,idxe,cv )
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RADI_FlimRkUpdate

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_FlimRkUpdate.F90,v $
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







