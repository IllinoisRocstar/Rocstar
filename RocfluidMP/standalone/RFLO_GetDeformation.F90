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
! Purpose: obtain deformation to boundary nodes from an external source.
!
! Description: none.
!
! Input: region     = grid dimensions and topology
!        boundMoved = flag for boundaries of region which have moved.
!
! Output: dNode = deformations at all six boundaries.
!
! Notes: variable dNode contains the whole 3-D field.
!
!******************************************************************************
!
! $Id: RFLO_GetDeformation.F90,v 1.6 2008/12/06 08:44:45 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_GetDeformation( region,boundMoved,dNode )

  USE ModDataTypes
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  LOGICAL :: boundMoved(6)

  REAL(RFREAL), POINTER :: dNode(:,:)

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iPatch, i, j, k, ijkN

! ... local variables
  INTEGER :: iLev, lbound, ibeg, iend, jbeg, jend, kbeg, kend, iNOff, ijNOff

  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER  :: patch

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'RFLO_GetDeformation',&
  'RFLO_GetDeformation.F90' )

! initialize variables --------------------------------------------------------

  iLev = 1     ! finest grid
  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

! zero out displacements

  dNode(:,:) = 0._RFREAL

! reset boundary movement flags

  boundMoved(:)  = .false.

! obtain displacements --------------------------------------------------------

  DO iPatch=1,region%nPatches
    patch  => region%levels(iLev)%patches(iPatch)
    lbound =  patch%lbound

    IF (patch%mixt%setMotion) THEN
      CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                      ibeg,iend,jbeg,jend,kbeg,kend )
      boundMoved(lbound) = .true.
      
      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkN = IndIJK(i,j,k,iNOff,ijNOff)
            dNode(XCOORD,ijkN) = patch%mixt%bndVel(XCOORD)*global%dtMin
            dNode(YCOORD,ijkN) = patch%mixt%bndVel(YCOORD)*global%dtMin
            dNode(ZCOORD,ijkN) = patch%mixt%bndVel(ZCOORD)*global%dtMin
          ENDDO
        ENDDO
      ENDDO

    ENDIF    ! motion prescribed internally
  ENDDO      ! iPatch

  IF (global%internDeform==0) THEN

! - give error message

    CALL ErrorStop( region%global,ERR_EXTERNAL_FUNCT,__LINE__, &
                   'No internal-motion at boundaries' )
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_GetDeformation

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_GetDeformation.F90,v $
! Revision 1.6  2008/12/06 08:44:45  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:57  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/08/19 15:41:05  mparmar
! Renamed patch variables
!
! Revision 1.3  2005/09/30 01:00:41  wasistho
! modified safety, replace used by global%internDeform
!
! Revision 1.2  2005/09/09 03:27:18  wasistho
! added internal boundary motion kernel
!
! Revision 1.1  2004/12/01 21:29:36  haselbac
! Initial revision after changing case
!
! Revision 1.4  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.3  2002/09/20 22:22:37  jblazek
! Finalized integration into GenX.
!
! Revision 1.2  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/08/15 19:43:11  jblazek
! Prepared solver to be optionally linked with an external driver.
!
!******************************************************************************







