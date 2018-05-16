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
! Purpose: converts a real number input to an Activeness
!
! Description: none.
!
! Input: val = real number from a ReadSection call
!
! Output: actv = integer used to designate Activeness
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_SetActiveness.F90,v 1.3 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_SetActiveness( global,val,actv )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  USE INRT_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER     :: global
  REAL(RFREAL),   INTENT(IN)  :: val
  INTEGER,        INTENT(OUT) :: actv

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_SetActiveness.F90,v $ $Revision: 1.3 $'

  CALL RegisterFunction( global,'INRT_SetActiveness',&
  'INRT_SetActiveness.F90' )

! begin -----------------------------------------------------------------------

! Convert input deck convention:  Active means Activeness = 1
!            to code convention:  Active means Activeness = INRT_ACT_ACTIVE

  actv = NINT(val) - 1 + INRT_ACT_ACTIVE

  IF (actv > INRT_ACT_ACTIVE) &
    CALL ErrorStop( global,ERR_INRT_BADACTV,__LINE__ )

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_SetActiveness

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_SetActiveness.F90,v $
! Revision 1.3  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:56:41  fnajjar
! Initial revision after changing case
!
! Revision 1.2  2004/04/15 16:04:21  jferry
! minor formatting (removed trailing spaces)
!
! Revision 1.1  2003/04/02 22:32:04  jferry
! codified Activeness and Permission structures for rocinteract
!
!******************************************************************************







