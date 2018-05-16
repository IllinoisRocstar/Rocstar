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
! Purpose: copy face vectors from RFLO datastructure.
!
! Description: none.
!
! Input: region = current region.
!
! Output: plag%si, sj, sk = face vectors.
!
! Notes: face vectors are not properly defined in dummy corner and edge cells.
!
!******************************************************************************
!
! $Id: PLAG_CopyFaceVectors.F90,v 1.3 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CopyFaceVectors( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, &
                            RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: ijk, iLev

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: ibn, idnbeg, idnend, ien, iNOff, ijNOff, &
             jdnbeg, jdnend, kdnbeg, kdnend

  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pSi, pSj, pSk
  
  TYPE(t_global), POINTER :: global
  TYPE(t_plag),   POINTER :: pPlag
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CopyFaceVectors.F90,v $ $Revision: 1.3 $'

  global => region%global
  
  CALL RegisterFunction( global,'PLAG_CopyFaceVectors',&
  'PLAG_CopyFaceVectors.F90' )
    
! Loop over all grid levels ---------------------------------------------------

  DO iLev=1,region%nGridLevels
  
! - Get dimensions ------------------------------------------------------------

    CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

    ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
    ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)

! - Set pointers --------------------------------------------------------------

    pPlag => region%levels(iLev)%plag
    pSi   => pPlag%si
    pSj   => pPlag%sj
    pSk   => pPlag%sk

! - i-direction ---------------------------------------------------------------

    DO ijk = ibn, ien
      pSi(XCOORD:ZCOORD,ijk) = region%levels(iLev)%grid%si(XCOORD:ZCOORD,ijk) 
    ENDDO

! - j-direction ---------------------------------------------------------------

    DO ijk = ibn, ien
      pSj(XCOORD:ZCOORD,ijk) = region%levels(iLev)%grid%sj(XCOORD:ZCOORD,ijk) 
    ENDDO

! - k-direction ---------------------------------------------------------------

    DO ijk = ibn, ien
      pSk(XCOORD:ZCOORD,ijk) = region%levels(iLev)%grid%sk(XCOORD:ZCOORD,ijk) 
    ENDDO

  ENDDO   ! iLev

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CopyFaceVectors

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CopyFaceVectors.F90,v $
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:26  fnajjar
! Initial revision after changing case
!
! Revision 1.1  2003/11/03 21:23:16  fnajjar
! Initial import of face vectors for PLAG datastructure
!
!******************************************************************************







