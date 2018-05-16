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
! ******************************************************************************
!
! Purpose: Reflect vector about plane.
!
! Description: None.
!
! Input: 
!   nx          x-component of normal vector of plane
!   ny          y-component of normal vector of plane
!   nz          z-component of normal vector of plane
!   xComp       x-component of vector
!   yComp       y-component of vector
!   zComp       z-component of vector
!
! Output: 
!   xComp       x-component of reflected vector
!   yComp       y-component of reflected vector
!   zComp       z-component of reflected vector
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReflectVector.F90,v 1.5 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReflectVector(nx,ny,nz,xComp,yComp,zComp)

  USE ModDataTypes
  USE ModParameters
   
  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  REAL(RFREAL), INTENT(IN) :: nx,ny,nz
  REAL(RFREAL), INTENT(INOUT) :: xComp,yComp,zComp

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  REAL(RFREAL) :: term 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: ReflectVector.F90,v $ $Revision: 1.5 $'

! ******************************************************************************
! Reflect vector
! ******************************************************************************

  term = 2.0_RFREAL*(xComp*nx + yComp*ny + zComp*nz)
  
  xComp = xComp - term*nx
  yComp = yComp - term*ny  
  zComp = zComp - term*nz
  
! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE ReflectVector

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReflectVector.F90,v $
! Revision 1.5  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.2  2004/05/05 20:37:03  fnajjar
! Fixed bug: Vector reflection does not need rc
!
! Revision 1.1  2004/04/08 01:32:31  haselbac
! Initial revision
!
! ******************************************************************************






