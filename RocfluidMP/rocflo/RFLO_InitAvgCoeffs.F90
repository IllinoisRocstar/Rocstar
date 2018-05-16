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
! Purpose: initialize cell2face and cell2edge averaging coefficients
!
! Description: 2 point c2f is set to 0.5, while 4 point c2e to 0.25. 
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: RFLO_InitAvgCoeffs.F90,v 1.3 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLO_InitAvgCoeffs( region )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global 
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : RFLO_GetNodeOffset, RFLO_GetDimensDummyNodes
  USE ModError
  USE ModParameters
  USE ModMPI
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: ijkN
   
! ... local variables
  TYPE(t_global), POINTER :: global

  INTEGER :: iLev, ibn, ien, iNOff, ijNOff
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  REAL(RFREAL), POINTER :: avgCoI(:,:), avgCoJ(:,:), avgCoK(:,:)
  
! ******************************************************************************

  global => region%global

  CALL RegisterFunction(global,'RFLO_InitAvgCoeffs',&
  'RFLO_InitAvgCoeffs.F90')

! get dimensions ---------------------------------------------------------------

  iLev = region%currLevel
  CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                 jdnbeg,jdnend,kdnbeg,kdnend )
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
  ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
  ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

! initialize cell2face coeffs to 0.5
  avgCoI => region%levels(iLev)%grid%c2fCoI
  avgCoJ => region%levels(iLev)%grid%c2fCoJ
  avgCoK => region%levels(iLev)%grid%c2fCoK

  DO ijkN=ibn,ien
    avgCoI(:,ijkN) = 0.5_RFREAL
    avgCoJ(:,ijkN) = 0.5_RFREAL
    avgCoK(:,ijkN) = 0.5_RFREAL
  ENDDO

! initialize cell2edge coeffs to 0.25
  avgCoI => region%levels(iLev)%grid%c2eCoI
  avgCoJ => region%levels(iLev)%grid%c2eCoJ
  avgCoK => region%levels(iLev)%grid%c2eCoK

  DO ijkN=ibn,ien
    avgCoI(:,ijkN) = 0.25_RFREAL
    avgCoJ(:,ijkN) = 0.25_RFREAL
    avgCoK(:,ijkN) = 0.25_RFREAL
  ENDDO

! finalize ---------------------------------------------------------------------

  CALL DeregisterFunction( global )
  
END SUBROUTINE RFLO_InitAvgCoeffs

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_InitAvgCoeffs.F90,v $
! Revision 1.3  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.1  2004/08/25 07:49:15  wasistho
! added RFLO_C2f/eAvgCoeffsDegec and RFLO_InitAvgCoeffs
!
!
! ******************************************************************************







