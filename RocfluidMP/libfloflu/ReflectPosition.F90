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
! Purpose: Reflect position about plane.
!
! Description: None.
!
! Input: 
!   nx          x-component of normal vector of plane
!   ny          y-component of normal vector of plane
!   nz          z-component of normal vector of plane
!   xc          x-component of centroid of plane
!   yc          y-component of centroid of plane
!   zc          z-component of centroid of plane
!   xComp       x-component of position
!   yComp       y-component of position
!   zComp       z-component of position
!
! Output: 
!   xComp       x-component of reflected position
!   yComp       y-component of reflected position
!   zComp       z-component of reflected position
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReflectPosition.F90,v 1.4 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReflectPosition(nx,ny,nz,xc,yc,zc,xComp,yComp,zComp)

  USE ModDataTypes
  USE ModParameters
   
  IMPLICIT NONE

! ******************************************************************************
! Declarations
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  REAL(RFREAL), INTENT(IN) :: nx,ny,nz,xc,yc,zc
  REAL(RFREAL), INTENT(INOUT) :: xComp,yComp,zComp

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  REAL(RFREAL) :: term 

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: ReflectPosition.F90,v $ $Revision: 1.4 $'

! ******************************************************************************
! Reflect position
! ******************************************************************************

  term = 2.0_RFREAL*((xComp-xc)*nx + (yComp-yc)*ny + (zComp-zc)*nz)
  
  xComp = xComp - term*nx
  yComp = yComp - term*ny  
  zComp = zComp - term*nz
  
! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE ReflectPosition

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReflectPosition.F90,v $
! Revision 1.4  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.1  2004/05/05 20:34:05  fnajjar
! Initial revision
!
! ******************************************************************************






