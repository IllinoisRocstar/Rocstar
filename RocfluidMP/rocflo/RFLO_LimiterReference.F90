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
! Purpose: calculate reference values for the limiter (upwind scheme).
!
! Description: none.
!
! Input: regions = all regions
!
! Output: regions = limiter reference values.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_LimiterReference.F90,v 1.3 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_LimiterReference( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... local variables
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RFLO_LimiterReference',&
  'RFLO_LimiterReference.F90' )

  global%limRef(1) = global%refDensity
  global%limRef(2) = global%refVelocity
  global%limRef(3) = global%refPressure
  global%limVolRef = global%refLength**3

  CALL DeregisterFunction( global )

END SUBROUTINE RFLO_LimiterReference

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_LimiterReference.F90,v $
! Revision 1.3  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.9  2003/05/29 17:28:43  jblazek
! Implemented Roe scheme.
!
! Revision 1.8  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.7  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.6  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.5  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.3  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.2  2002/01/08 22:09:17  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.1  2001/12/19 23:09:22  jblazek
! Added routines to read grid and solution.
!
!******************************************************************************







