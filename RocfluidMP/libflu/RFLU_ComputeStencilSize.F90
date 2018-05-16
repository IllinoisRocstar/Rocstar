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
! Purpose: Compute stencil size for given order.
!
! Description: None.
!
! Input:
!   global              Global pointer
!   factor              Safety factor
!   order               Polynomial order
!
! Output: 
!   RFLU_ComputeStencilSize     (Minimum) number of points in stencil
!
! Notes: 
!   1. Include additional support (+1) to allow for interpolation.
!   2. Make larger than necessary to allow for least-squares if order > 1.
!   3. If order = 0, do not include factor. 
!
!******************************************************************************
!
! $Id: RFLU_ComputeStencilSize.F90,v 1.4 2008/12/06 08:44:12 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

INTEGER FUNCTION RFLU_ComputeStencilSize(global,factor,order)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  USE ModGlobal, ONLY: t_global
      
  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  INTEGER, INTENT(IN) :: factor,order
  TYPE(t_global), POINTER :: global

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_ComputeStencilSize.F90,v $ $Revision: 1.4 $'

  CALL RegisterFunction(global,'RFLU_ComputeStencilSize',&
  'RFLU_ComputeStencilSize.F90')

! *****************************************************************************
! Compute stencil size
! *****************************************************************************

  IF ( order > 0 ) THEN 
    RFLU_ComputeStencilSize = factor*(order + 1)*(order + 2)*(order + 3)/6
  ELSE 
    RFLU_ComputeStencilSize = 1
  END IF ! order

! *****************************************************************************
! End
! *****************************************************************************

  CALL DeregisterFunction(global)

END FUNCTION RFLU_ComputeStencilSize

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ComputeStencilSize.F90,v $
! Revision 1.4  2008/12/06 08:44:12  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:25  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:16  haselbac
! Removed tabs
!
! Revision 1.1  2003/12/04 03:23:34  haselbac
! Initial revision
!
!******************************************************************************







