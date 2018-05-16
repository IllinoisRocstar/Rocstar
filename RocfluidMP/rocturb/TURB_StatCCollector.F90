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
! Purpose: Collect cell variables of interest to be time averaged for check or
!          comparison
!
! Description: As a generic routine, number and kind of collected variables
!              are not specified. The cell values are accumulated in time.
!
! Input: region = data of current region
!        iBegSv = begin index of averaged field sv
!        iEndSv = end index of averaged field sv
!        colVar = array of collected variables
!
! Output: turb%st(iBegSv:iEndSv,:) at cell centers including dummy cells.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_StatCCollector.F90,v 1.3 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_StatCCollector( region,iBegSt,iEndSt,colVar )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, &
                            RFLO_GetCellOffset

#include "Indexing.h"
#endif
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  INTEGER        :: iBegSt, iEndSt
  REAL(RFREAL), POINTER :: colVar(:,:)

! ... loop variables
  INTEGER :: i, j, k, l, m

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: ijkC
  REAL(RFREAL), POINTER :: st(:,:)

#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: iLev,iCOff,ijCOff
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_StatCCollector.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'Turb_StatCCollector',&
  'TURB_StatCCollector.F90' )

! check some input arguments --------------------------------------------------

  IF (iEndSt > region%turbInput%nSt) THEN
    CALL ErrorStop( global,ERR_TURB_STATSINPUT,__LINE__, &
        'index of collected vars larger than nSt (allocated nmbr stats vars)' )
  ENDIF
  IF (iBegSt > iEndSt) THEN
    CALL ErrorStop( global,ERR_TURB_STATSINPUT,__LINE__, &
        'begin index of collected vars larger than end index' )
  ENDIF

#ifdef RFLO
! get parameters --------------------------------------------------------

  iLev  = region%currLevel

! get dimensions and pointers

  CALL RFLO_GetDimensDummy( region,ilev,idcbeg,idcend, &
                            jdcbeg,jdcend,kdcbeg,kdcend )
  CALL RFLO_GetCellOffset( region,ilev,iCOff,ijCOff )

  st => region%levels(ilev)%turb%st

! interpolate from face to cell

  DO k=kdcbeg,kdcend
    DO j=jdcbeg,jdcend
      DO i=idcbeg,idcend

        ijkC  = IndIJK(i     ,j     ,k     ,iCOff,ijCOff)

        DO l=iBegSt,iEndSt   
          m = l - iBegSt + 1  
          st(l,ijkC) = st(l,ijkC)+colVar(m,ijkC)
        ENDDO

      ENDDO ! i
    ENDDO   ! j
  ENDDO     ! k
#endif
#ifdef RFLU
! get dimensions and pointers ------------------------------------------

  st => region%turb%st
  st =  0._RFREAL

  DO ijkC=1,region%grid%nCells

! - region (here denoted as 'global') to local (type) mapping

    DO l=iBegSt,iEndSt   
      m = l - iBegSt + 1  
      st(l,ijkC) = st(l,ijkC)+colVar(m,ijkC)
    ENDDO   ! l
  ENDDO     ! ijkC
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_StatCCollector

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_StatCCollector.F90,v $
! Revision 1.3  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/10/22 23:22:18  wasistho
! added statistics collector based on cell centered variables
!
!
!******************************************************************************







