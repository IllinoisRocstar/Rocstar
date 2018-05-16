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
! Purpose: close window to fluids code.
!
! Description: none.
!
! Input: winName = name of fluids window.
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: Rocflo_unload_module.F90,v 1.3 2008/12/06 08:44:01 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE Rocflo_unload_module( winName )

  USE ModRocstar, ONLY : t_globalGenx, associate_pointer
  IMPLICIT NONE
  INCLUDE 'comf90.h'

  INTERFACE 
    SUBROUTINE COM_get_pointer( attr,ptr,asso )
      USE ModRocstar, ONLY : t_globalGenx
      CHARACTER(*), INTENT(IN) :: attr
      TYPE(t_globalGenx), POINTER :: ptr
      EXTERNAL asso
    END SUBROUTINE COM_get_pointer
  END INTERFACE

! ... parameters
  CHARACTER(*) :: winName

! ... local variables
  TYPE(t_globalGenx), POINTER  :: glb

!******************************************************************************

  CALL COM_get_pointer( winName//'.global',glb,associate_pointer )
  DEALLOCATE( glb )

  CALL COM_delete_window( winName )
  
END SUBROUTINE Rocflo_unload_module

!******************************************************************************
!
! RCS Revision history:
!
! $Log: Rocflo_unload_module.F90,v $
! Revision 1.3  2008/12/06 08:44:01  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:23:57  haselbac
! Initial revision after changing case
!
! Revision 1.2  2002/09/25 19:01:00  jblazek
! Removed :: in front of asso (SUN does not like it).
!
! Revision 1.1  2002/09/20 22:22:34  jblazek
! Finalized integration into GenX.
!
!******************************************************************************






