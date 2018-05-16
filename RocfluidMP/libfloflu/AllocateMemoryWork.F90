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
! Purpose: allocate memory for all variables associated with the mixture
!          for all active regions on current processor.
!
! Description: none.
!
! Input: region = current region
!
! Output: region%work1D = 1D work array
!         region%work2D = 2D work array.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: AllocateMemoryWork.F90,v 1.3 2008/12/06 08:44:08 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE AllocateMemoryWork( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... local variables
  INTEGER :: errorFlag

!******************************************************************************

  CALL RegisterFunction( region%global,'AllocateMemoryWork',&
  'AllocateMemoryWork.F90' )

! allocate work space

#ifdef RFLO
  ALLOCATE( region%work1D(region%dimWork1D),stat=errorFlag )
  region%global%error = errorFlag
  IF (region%global%error /= 0) &
    CALL ErrorStop( region%global,ERR_ALLOCATE,__LINE__ )

  ALLOCATE( region%work2D(region%dimWork2D(1),region%dimWork2D(2)), &
            stat=errorFlag )
  region%global%error = errorFlag
  IF (region%global%error /= 0) &
    CALL ErrorStop( region%global,ERR_ALLOCATE,__LINE__ )
#endif

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE AllocateMemoryWork

!******************************************************************************
!
! RCS Revision history:
!
! $Log: AllocateMemoryWork.F90,v $
! Revision 1.3  2008/12/06 08:44:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:47:48  haselbac
! Initial revision after changing case
!
! Revision 1.9  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.8  2002/10/08 15:48:35  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.7  2002/09/05 17:40:19  jblazek
! Variable global moved into regions().
!
! Revision 1.6  2002/08/30 19:08:58  jblazek
! Dimensions of work arrays now set in derivedInputValues.
!
! Revision 1.5  2002/08/24 03:12:54  wasistho
! put safety within #ifdef TURB
!
! Revision 1.4  2002/08/07 20:40:07  jblazek
! Added RFLO preprocessor directive around memory allocation.
!
! Revision 1.3  2002/07/15 22:03:34  jblazek
! Extended allocation of work space to physical modules.
!
! Revision 1.2  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.1  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
!******************************************************************************







