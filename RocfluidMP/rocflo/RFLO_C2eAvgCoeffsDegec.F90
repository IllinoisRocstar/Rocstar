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
! ******************************************************************************
!
! Purpose: correct cell2edge averaging coeffs. in degenerated edges/corners.
!
! Description: as the e/c is practically no exist, no real coordinates are 
!              available, hence the 4 point coefficients are set to 0.25.
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: RFLO_C2eAvgCoeffsDegec.F90,v 1.3 2008/12/06 08:44:25 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLO_C2eAvgCoeffsDegec( region )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global 
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetEdgeCellsIndices, &
                            RFLO_GetCornerCellsIndices
  USE ModError
  USE ModParameters
  USE ModMPI
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iedge, icorner, i, j, k
   
! ... local variables
  TYPE(t_global), POINTER :: global
  TYPE(t_level), POINTER  :: level

  INTEGER :: iLev, ijkN, iNOff, ijNOff
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend
  REAL(RFREAL), POINTER :: avgCoI(:,:), avgCoJ(:,:), avgCoK(:,:)
  
! ******************************************************************************

  global => region%global

  CALL RegisterFunction(global,'RFLO_C2eAvgCoeffsDegec',&
  'RFLO_C2eAvgCoeffsDegec.F90')

! get parameters and pointers --------------------------------------------------

  iLev   =  region%currLevel
  level  => region%levels(iLev)
  avgCoI => region%levels(iLev)%grid%c2eCoI
  avgCoJ => region%levels(iLev)%grid%c2eCoJ
  avgCoK => region%levels(iLev)%grid%c2eCoK

  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! treat degenerat edges      

  DO iedge = 1,12
    IF (level%edgeCells(iedge)%degenrt /= DEGENERAT_NONE) THEN
      CALL RFLO_GetEdgeCellsIndices( region,iLev,iedge, &
                                     ibeg,iend,jbeg,jend,kbeg,kend )
      CALL CalcAvgCo
    ENDIF  ! degenrt
  ENDDO    ! iedge          
      
! treat degenerat corners

  DO icorner = 1,8
    IF (level%cornerCells(icorner)%degenrt /= DEGENERAT_NONE) THEN
      CALL RFLO_GetCornerCellsIndices( region,iLev,icorner, &
                                       ibeg,iend,jbeg,jend,kbeg,kend )
      CALL CalcAvgCo
    ENDIF  ! degenrt
  ENDDO    ! icorner

! finalize ---------------------------------------------------------------------

  CALL DeregisterFunction( global )

! =============================================================================
!   Averaging coefficients assignment
! =============================================================================

CONTAINS

  SUBROUTINE CalcAvgCo

! - I-edge
    DO k=kbeg,kend+1
      DO j=jbeg,jend+1
        DO i=ibeg,iend
          ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff)        
          avgCoI(:,ijkN) = 0.25_RFREAL
        ENDDO
      ENDDO
    ENDDO

! - J-face
    DO k=kbeg,kend+1
      DO j=jbeg,jend
        DO i=ibeg,iend+1
          ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff)        
          avgCoJ(:,ijkN) = 0.25_RFREAL
        ENDDO
      ENDDO
    ENDDO

! - K-face
    DO k=kbeg,kend
      DO j=jbeg,jend+1
        DO i=ibeg,iend+1
          ijkN = IndIJK(i ,j ,k ,iNOff,ijNOff)        
          avgCoK(:,ijkN) = 0.25_RFREAL
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE CalcAvgCo
  
END SUBROUTINE RFLO_C2eAvgCoeffsDegec

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_C2eAvgCoeffsDegec.F90,v $
! Revision 1.3  2008/12/06 08:44:25  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:36  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:38  wasistho
! lower to upper case
!
! Revision 1.1  2004/08/25 07:49:15  wasistho
! added RFLO_C2f/eAvgCoeffsDegec and RFLO_InitAvgCoeffs
!
!
! ******************************************************************************







