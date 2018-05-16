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
! Purpose: searches list of materials for one with the given name
!
! Description: none.
!
! Input:  name = input name of material
!
! Output: material points to the correct element of global%materials(:)
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_SetMaterial.F90,v 1.3 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_SetMaterial(global,material,name)

  USE ModDataTypes
  USE ModGlobal,    ONLY : t_global
  USE ModMaterials, ONLY : t_material
  USE ModError
  IMPLICIT NONE

! ... parameters
  TYPE(t_global),   POINTER    :: global
  TYPE(t_material), POINTER    :: material
  CHARACTER(*),     INTENT(in) :: name

! ... loop variables
  INTEGER :: iMat

! ... local variables
  LOGICAL :: found

  CHARACTER(CHRLEN) :: RCSIdentString

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_SetMaterial.F90,v $ $Revision: 1.3 $'

  CALL RegisterFunction( global,'INRT_SetMaterial',&
  'INRT_SetMaterial.F90' )

! begin -----------------------------------------------------------------------

  found = .FALSE.

  DO iMat = 1,global%nMaterials

    IF (TRIM(name) == TRIM(global%materials(iMat)%name)) THEN
      material => global%materials(iMat)
      found = .TRUE.
      EXIT
    END IF ! name

  END DO ! iMat

  IF (.NOT.found) CALL ErrorStop( global,ERR_INRT_BADMAT,__LINE__,TRIM(name))

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_SetMaterial

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_SetMaterial.F90,v $
! Revision 1.3  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:56:42  fnajjar
! Initial revision after changing case
!
! Revision 1.1  2003/03/24 23:23:25  jferry
! converted from libfloflu routine to rocinteract routine
!
!******************************************************************************







