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
! Purpose: output surface grid for all patches interacting with GenX.
!
! Description: none.
!
! Input: iReg   = global number of current region
!        region = dimensions of patches, types of BC`s, grid.
!
! Output: to file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: SURF_WriteSurfaceGrid.F90,v 1.4 2008/12/06 08:44:52 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE WriteSurfaceGrid( iReg,region )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModBndPatch, ONLY   : t_patch
  USE ModDataStruct, ONLY : t_region
  USE SURF_ModInterfaces, ONLY : RFLO_GetPatchIndicesNodes, RFLO_GetNodeOffset
  USE ModError
  USE ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  INTEGER :: iReg

  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iPatch, i, j, k, ijkN, n1, n2, ng1, ng2

! ... local variables
  INTEGER :: iLev, bcType, lbound, iNOff, ijNOff, dims(2)
  INTEGER :: ibeg, iend, jbeg, jend, kbeg, kend, icount, pid, errorFlag

  REAL(RFREAL), POINTER :: xyz(:,:), surfCoord(:,:,:)

  TYPE(t_global), POINTER :: global
  TYPE(t_patch) , POINTER :: patch

!******************************************************************************

  global => region%global

  CALL RegisterFunction( global,'WriteSurfaceGrid',&
  'SURF_WriteSurfaceGrid.F90' )

! store pointer to coordinates ------------------------------------------------

  iLev   = region%currLevel
  icount = 0

  CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )

  xyz => region%levels(iLev)%grid%xyz

! loop over all cells of the patch (if an interface) --------------------------

  DO iPatch=1,region%nPatches

    patch   => region%levels(iLev)%patches(iPatch)
    bcType  =  patch%bcType
    lbound  =  patch%lbound

    IF (patch%bcCoupled == BC_EXTERNAL) THEN        ! interacting
      icount  = icount + 1
      pid     = iReg*REGOFF + icount

! --- get dimensions, allocate memory

      dims(1) = ABS(patch%l1end-patch%l1beg) + 2    ! nodal values
      dims(2) = ABS(patch%l2end-patch%l2beg) + 2
      ALLOCATE( surfCoord(3,dims(1),dims(2)),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

      CALL RFLO_GetPatchIndicesNodes( region,patch,iLev, &
                                      ibeg,iend,jbeg,jend,kbeg,kend )

! --- copy coordinates to temporary array

      DO k=kbeg,kend
        DO j=jbeg,jend
          DO i=ibeg,iend
            ijkN = IndIJK(i,j,k,iNOff,ijNOff)
            IF      (lbound==1 .OR. lbound==2) THEN
              IF (lbound == 1) THEN
                ng1 = j - jbeg + 1
              ELSE
                ng1 = jend - j + 1
              ENDIF
              ng2 = k - kbeg + 1
            ELSE IF (lbound==3 .OR. lbound==4) THEN
              ng1 = k - kbeg + 1
              IF (lbound == 3) THEN
                ng2 = i - ibeg + 1
              ELSE
                ng2 = iend - i + 1
              ENDIF
            ELSE IF (lbound==5 .OR. lbound==6) THEN
              IF (lbound == 5) THEN
                ng1 = i - ibeg + 1
              ELSE
                ng1 = iend - i + 1
              ENDIF
              ng2 = j - jbeg + 1
            ENDIF
            surfCoord(1,ng1,ng2) = xyz(XCOORD,ijkN)
            surfCoord(2,ng1,ng2) = xyz(YCOORD,ijkN)
            surfCoord(3,ng1,ng2) = xyz(ZCOORD,ijkN)
          ENDDO
        ENDDO
      ENDDO

! --- write to file

      WRITE(IF_PLOT,*) pid
      WRITE(IF_PLOT,*) dims(1),dims(2)
      DO ng2=1,dims(2)
        DO ng1=1,dims(1)
!RAF          WRITE(IF_PLOT,*) surfCoord(1,ng1,ng2), &
          WRITE(IF_PLOT,'(3E24.14)') surfCoord(1,ng1,ng2), &
                           surfCoord(2,ng1,ng2), &
                           surfCoord(3,ng1,ng2)
        ENDDO
      ENDDO

! --- release temporary storage

      DEALLOCATE( surfCoord,stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
      NULLIFY( surfCoord )

    ENDIF  ! external BC
  ENDDO    ! iPatch

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE WriteSurfaceGrid

!******************************************************************************
!
! RCS Revision history:
!
! $Log: SURF_WriteSurfaceGrid.F90,v $
! Revision 1.4  2008/12/06 08:44:52  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:18:02  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/12/03 03:35:56  wasistho
! rflo_modinterfacessurf to surf_modinterfaces
!
! Revision 1.1  2004/12/03 02:47:00  wasistho
! added prefix
!
! Revision 1.1  2004/12/03 00:49:09  wasistho
! lower to upper case
!
! Revision 1.6  2004/06/30 04:08:08  wasistho
! moved Genx related parameter REGOFF to ModParameters
!
! Revision 1.5  2004/04/14 20:28:38  rfiedler
! Use formatted output to avoid exponentials with D, which C hates.
!
! Revision 1.4  2003/05/15 02:57:07  jblazek
! Inlined index function.
!
! Revision 1.3  2003/03/20 22:35:02  haselbac
! Renamed ModInterfaces
!
! Revision 1.2  2003/03/20 19:48:09  haselbac
! Corrected mistake in phased check-in
!
! Revision 1.1  2002/10/19 00:40:31  jblazek
! Added utility (rflosurf) to write out surface grids for GenX.
!
!******************************************************************************







