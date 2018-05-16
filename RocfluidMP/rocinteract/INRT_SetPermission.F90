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
! Purpose: converts a real number input to a Permission level
!
! Description: none.
!
! Input: val = real number from a ReadSection call
!
! Output: perm = integer used to designate Permission level
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_SetPermission.F90,v 1.3 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************
SUBROUTINE INRT_SetPermission( global,val,perm )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  USE INRT_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER     :: global
  REAL(RFREAL),   INTENT(IN)  :: val
  INTEGER,        INTENT(OUT) :: perm

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_SetPermission.F90,v $ $Revision: 1.3 $'

  CALL RegisterFunction( global,'INRT_SetPermission',&
  'INRT_SetPermission.F90' )

! begin -----------------------------------------------------------------------

  SELECT CASE (NINT(val))

  CASE (0)                 ! value used in input deck
    perm = INRT_PERM_BLOCK ! value used inside the code

  CASE (1)                 ! value used in input deck
    perm = INRT_PERM_PMASS ! value used inside the code

  CASE (2)                 ! value used in input deck
    perm = INRT_PERM_PMOME ! value used inside the code

  CASE (3)                 ! value used in input deck
    perm = INRT_PERM_PALL  ! value used inside the code

  CASE DEFAULT
    CALL ErrorStop( global,ERR_INRT_BADPERM,__LINE__ )

  END SELECT ! val

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_SetPermission

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_SetPermission.F90,v $
! Revision 1.3  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:56:45  fnajjar
! Initial revision after changing case
!
! Revision 1.1  2003/04/02 22:32:04  jferry
! codified Activeness and Permission structures for rocinteract
!
!******************************************************************************







