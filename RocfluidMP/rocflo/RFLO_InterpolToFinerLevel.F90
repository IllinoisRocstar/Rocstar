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
! Purpose: interpolate solution from coarser to finer grid.
!
! Description: none.
!
! Input: region = conservative variables of the finer grid.
!
! Output: region = conservative variables of the coarser grid.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_InterpolToFinerLevel.F90,v 1.3 2008/12/06 08:44:27 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_InterpolToFinerLevel( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModInterfaces, ONLY : 
  USE ModError
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables


! ... local variables

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_InterpolToFinerLevel',&
  'RFLO_InterpolToFinerLevel.F90' )



  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_InterpolToFinerLevel

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_InterpolToFinerLevel.F90,v $
! Revision 1.3  2008/12/06 08:44:27  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:38  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/11/29 20:51:39  wasistho
! lower to upper case
!
! Revision 1.9  2003/11/20 16:40:39  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/05/15 02:57:04  jblazek
! Inlined index function.
!
! Revision 1.4  2002/09/05 17:40:21  jblazek
! Variable global moved into regions().
!
! Revision 1.3  2002/02/21 23:25:06  jblazek
! Blocks renamed as regions.
!
! Revision 1.2  2002/01/23 03:51:25  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.1  2002/01/16 22:03:35  jblazek
! Added time-stepping routines.
!
!******************************************************************************







