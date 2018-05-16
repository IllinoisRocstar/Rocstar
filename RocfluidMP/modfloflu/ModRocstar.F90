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
! Purpose: define global variable for coupling with GenX.
!
! Description: none
!
! Notes: none
!
!******************************************************************************
!
! $Id: ModGenx.F90,v 1.6 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

MODULE ModRocstar

  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModDataStruct, ONLY : t_region
#endif
#ifdef RFLU
  USE ModDataStruct, ONLY: t_level
#endif
  IMPLICIT NONE

#ifdef GENX
  INCLUDE 'rocmanf90.h'
#endif

  TYPE t_globalGenx
    LOGICAL :: isDummy

    TYPE(t_global), POINTER :: global
#ifdef RFLO
    TYPE(t_region), POINTER :: regions(:)
#endif 
#ifdef RFLU
    TYPE(t_level), POINTER :: levels(:)
#endif
  END TYPE t_globalGenx

CONTAINS  
  SUBROUTINE associate_pointer( attr, ptr)
    TYPE(t_globalGenx), POINTER :: attr, ptr
    ptr => attr
  END SUBROUTINE ASSOCIATE_POINTER

END MODULE ModRocstar

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModGenx.F90,v $
! Revision 1.6  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2002/11/15 22:11:31  haselbac
! Enclosed INCLUDE within ifdef
!
! Revision 1.3  2002/11/15 21:24:08  haselbac
! Added INCLUDE rocmanf90.h for use in interfaces
!
! Revision 1.2  2002/10/05 18:59:11  haselbac
! Integration of RFLU into GENX
!
! Revision 1.1  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
!******************************************************************************






