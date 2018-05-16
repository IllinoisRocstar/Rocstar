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
! Purpose: Scale and rotate vector.
!
! Description: None.
!
! Input: vect = vector to scale and rotate.
!
! Output: vect = scaled and rotated.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: ScaleRotateVector.F90,v 1.3 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ScaleRotateVector( global,vect )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModError
  USE ModParameters

  IMPLICIT NONE
    
! ... parameters
  REAL(RFREAL), POINTER :: vect(:,:)

  TYPE(t_global), POINTER :: global
  
! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: ic, nComp

  REAL(RFREAL) :: angleX, angleY, angleZ, cx, cy, cz, sx, sy, sz, x, y, z

!******************************************************************************

  RCSIdentString = '$RCSfile: ScaleRotateVector.F90,v $ $Revision: 1.3 $'

  IF ( global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Scaling and rotating vector...'
    WRITE(STDOUT,'(A,3X,A,1X,E15.6)') SOLVER_NAME,'X-scale:',global%scaleX
    WRITE(STDOUT,'(A,3X,A,1X,E15.6)') SOLVER_NAME,'Y-scale:',global%scaleY
    WRITE(STDOUT,'(A,3X,A,1X,E15.6)') SOLVER_NAME,'Z-scale:',global%scaleZ
    WRITE(STDOUT,'(A,3X,A,1X,E15.6)') SOLVER_NAME,'X-angle:',global%angleX
    WRITE(STDOUT,'(A,3X,A,1X,E15.6)') SOLVER_NAME,'Y-angle:',global%angleY
    WRITE(STDOUT,'(A,3X,A,1X,E15.6)') SOLVER_NAME,'Z-angle:',global%angleZ
  END IF ! global%verbLevel

  nComp = SIZE(vect,DIM=2)

! scale vector

  DO ic = 1,nComp
    vect(1,ic) = global%scaleX*vect(1,ic)
    vect(2,ic) = global%scaleY*vect(2,ic)
    vect(3,ic) = global%scaleZ*vect(3,ic)
  END DO ! ic

! rotate vector

  angleX = global%angleX*global%rad
  angleY = global%angleY*global%rad
  angleZ = global%angleZ*global%rad   

  cx = COS(angleX)
  sx = SIN(angleX)
  cy = COS(angleY)
  sy = SIN(angleY)
  cz = COS(angleZ)
  sz = SIN(angleZ)
  
  DO ic = 1,nComp
    x = vect(1,ic)
    y = vect(2,ic)
    z = vect(3,ic)

    vect(1,ic) = cy*cz*x - (sx*sy*cz + cx*sz)*y - (cx*sy*cz - sx*sz)*z
    vect(2,ic) = cy*sz*x - (sx*sy*sz - cx*cz)*y - (cx*sy*sz + sx*cz)*z
    vect(3,ic) =    sy*x + (           sx*cy)*y + (           cx*cy)*z
  END DO ! ic
  
  IF ( global%verbLevel > VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Scaling and rotating vector done.'
  END IF ! global%verbLevel
  
END SUBROUTINE ScaleRotateVector

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ScaleRotateVector.F90,v $
! Revision 1.3  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:51:21  haselbac
! Initial revision after changing case
!
! Revision 1.5  2004/01/29 22:52:48  haselbac
! Clean up
!
! Revision 1.4  2003/02/01 00:28:20  haselbac
! Some clean up, added diagnostic output
!
! Revision 1.3  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.2  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.1  2002/05/07 18:51:27  haselbac
! Initial revision
!
!******************************************************************************






