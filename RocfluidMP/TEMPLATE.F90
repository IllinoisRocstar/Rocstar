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
! Purpose:
!
! Description:
!
! Input:
!
! Output:
!
! Notes:
!
!******************************************************************************
!
! $Id: TEMPLATE.F90,v 1.5 2008/12/06 08:44:00 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE 
FUNCTION 

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : 
  IMPLICIT NONE

! ... parameters


! ... loop variables


! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

!******************************************************************************

  RCSIdentString = '$RCSfile: TEMPLATE.F90,v $ $Revision: 1.5 $'

  CALL RegisterFunction( global,'function name',&
  'TEMPLATE.F90' )

! comment ---------------------------------------------------------------------


! comment



! - comment



! --- comment



  CALL DeregisterFunction( global )

END SUBROUTINE 
END FUNCTION 

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TEMPLATE.F90,v $
! Revision 1.5  2008/12/06 08:44:00  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:14  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2002/09/05 18:29:43  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/07/05 23:20:45  jblazek
! Corrected bug in perfgasDependentVars.F90; did some cosmetics.
!
! Revision 1.1.1.1  2001/12/03 21:44:04  jblazek
! Import of RocfluidMP
!
!******************************************************************************







