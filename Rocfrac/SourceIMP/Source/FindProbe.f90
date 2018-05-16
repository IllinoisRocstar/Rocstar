!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
SUBROUTINE FindProbe(glb, myid)

  USE ROCSTAR_RocFrac

!!****f* Rocfrac/Source/feminp.f90
!!
!!  NAME
!!     feminp
!!
!!  FUNCTION
!!
!!     READ INPUT INFORMATION (i.e. Analysis Deck File)
!!
!!  INPUTS
!!     glb -- global array
!!     myid -- processor id (starting at 0)
!!
!!****

  IMPLICIT NONE

  TYPE(ROCFRAC_GLOBAL) :: glb

  INTEGER :: myid

  REAL*8 :: xx,yy,zz,size1,size2,size3,size4,size5,size6,size7,size8,size9,size10,size11,size12
  REAL*8 :: LongestEdge

  CHARACTER*4 :: ichr1, ichr2
  
  INTEGER :: i, j, ierr
  LOGICAL :: inside
  LOGICAL :: PointOnProc

  LongestEdge = 0.d0

  IF(glb%iElType.EQ.8)THEN
     DO i = 1, glb%NumElVol
        CALL SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(1,i), glb%ElConnVol(4,i), size1)
        CALL SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(5,i), glb%ElConnVol(8,i), size2)
        CALL SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(6,i), glb%ElConnVol(7,i), size3)
        CALL SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(2,i), glb%ElConnVol(3,i), size4)
        
        CALL SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(4,i), glb%ElConnVol(8,i), size5)
        CALL SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(8,i), glb%ElConnVol(7,i), size6)
        CALL SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(7,i), glb%ElConnVol(3,i), size7)
        CALL SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(3,i), glb%ElConnVol(4,i), size8)
        
        CALL SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(1,i), glb%ElConnVol(5,i), size9)
        CALL SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(5,i), glb%ElConnVol(6,i), size10)
        CALL SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(6,i), glb%ElConnVol(2,i), size11)
        CALL SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(2,i), glb%ElConnVol(1,i), size12)
        
        LongestEdge = MAX(size1, size2, size3, size4, size5, size6, size7, size8, size9, size10, size11, size12, LongestEdge)
     ENDDO
     
  ELSE
     
     DO i = 1, glb%NumElVol
!
! -- Find the size of the smallest element
!
        xx = glb%MeshCoor(1,glb%ElConnVol(1,i)) - glb%MeshCoor(1,glb%ElConnVol(2,i))
        yy = glb%MeshCoor(2,glb%ElConnVol(1,i)) - glb%MeshCoor(2,glb%ElConnVol(2,i))
        zz = glb%MeshCoor(3,glb%ElConnVol(1,i)) - glb%MeshCoor(3,glb%ElConnVol(2,i))
        size1 = SQRT(xx*xx+yy*yy+zz*zz)
        xx = glb%MeshCoor(1,glb%ElConnVol(2,i)) - glb%MeshCoor(1,glb%ElConnVol(3,i))
        yy = glb%MeshCoor(2,glb%ElConnVol(2,i)) - glb%MeshCoor(2,glb%ElConnVol(3,i))
        zz = glb%MeshCoor(3,glb%ElConnVol(2,i)) - glb%MeshCoor(3,glb%ElConnVol(3,i))
        size2 = SQRT(xx*xx+yy*yy+zz*zz)
        xx = glb%MeshCoor(1,glb%ElConnVol(3,i)) - glb%MeshCoor(1,glb%ElConnVol(1,i))
        yy = glb%MeshCoor(2,glb%ElConnVol(3,i)) - glb%MeshCoor(2,glb%ElConnVol(1,i))
        zz = glb%MeshCoor(3,glb%ElConnVol(3,i)) - glb%MeshCoor(3,glb%ElConnVol(1,i))
        size3 = SQRT(xx*xx+yy*yy+zz*zz)
        xx = glb%MeshCoor(1,glb%ElConnVol(4,i)) - glb%MeshCoor(1,glb%ElConnVol(1,i))
        yy = glb%MeshCoor(2,glb%ElConnVol(4,i)) - glb%MeshCoor(2,glb%ElConnVol(1,i))
        zz = glb%MeshCoor(3,glb%ElConnVol(4,i)) - glb%MeshCoor(3,glb%ElConnVol(1,i))
        size4 = SQRT(xx*xx+yy*yy+zz*zz)
        xx = glb%MeshCoor(1,glb%ElConnVol(4,i)) - glb%MeshCoor(1,glb%ElConnVol(2,i))
        yy = glb%MeshCoor(2,glb%ElConnVol(4,i)) - glb%MeshCoor(2,glb%ElConnVol(2,i))
        zz = glb%MeshCoor(3,glb%ElConnVol(4,i)) - glb%MeshCoor(3,glb%ElConnVol(2,i))
        size5 = SQRT(xx*xx+yy*yy+zz*zz)
        xx = glb%MeshCoor(1,glb%ElConnVol(4,i)) - glb%MeshCoor(1,glb%ElConnVol(3,i))
        yy = glb%MeshCoor(2,glb%ElConnVol(4,i)) - glb%MeshCoor(2,glb%ElConnVol(3,i))
        zz = glb%MeshCoor(3,glb%ElConnVol(4,i)) - glb%MeshCoor(3,glb%ElConnVol(3,i))
        size6 = SQRT(xx*xx+yy*yy+zz*zz)
        LongestEdge = MAX(size1,size2,size3,size4,size5,size6,LongestEdge)
     
     ENDDO
  ENDIF

  ALLOCATE(glb%PointOnProc(1:glb%NumProbesNd))

  DO i = 1, glb%NumProbesNd

     glb%PointOnProc(i) = .FALSE.

     DO j = 1, glb%NumNP

        CALL sphere_imp_contains_point_3d ( LongestEdge, glb%ProbeCoorNd(:,i), glb%MeshCoor(:,j), inside )

        IF(inside)THEN

           glb%PointOnProc(i) = .TRUE.
           glb%ProbeNd(i) = j
           LongestEdge = SQRT ( SUM ( ( glb%ProbeCoorNd(1:3,i) - glb%MeshCoor(1:3,j)  )**2 ) )
        ENDIF

     ENDDO

     IF(glb%PointOnProc(i))THEN
        WRITE(ichr1,'(i4.4)') i
        WRITE(ichr2,'(I4.4)') myid

        OPEN(440+i,FILE='Rocfrac/Rocout/Probe.'//ichr1//'.'//ichr2,POSITION='REWIND')
        WRITE(440+i,*) '# Probe Data Coordinate', glb%MeshCoor(1:3,glb%ProbeNd(i))
        CLOSE(440+i)
     ENDIF

  ENDDO

END SUBROUTINE FindProbe


SUBROUTINE sphere_imp_contains_point_3d ( r, center, p, inside )

!*******************************************************************************
!
!! SPHERE_IMP_CONTAINS_POINT_3D: point in implicit sphere in 3D?
!
!  Discussion:
!
!    An implicit sphere in 3D satisfies the equation:
!
!      sum ( ( P(1:NDIM) - CENTER(1:NDIM) )**2 ) = R**2
!
!  Modified:
!
!    05 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the sphere.
!
!    Input, real ( kind = 8 ) CENTER(3), the center of the sphere.
!
!    Input, real ( kind = 8 ) P(3), the point to be checked.
!
!    Output, logical INSIDE, is TRUE if the point is inside the sphere.
!
    IMPLICIT NONE

    INTEGER, PARAMETER :: ndim = 3

    REAL ( kind = 8 ) center(ndim)
    LOGICAL inside
    REAL ( kind = 8 ) p(ndim)
    REAL ( kind = 8 ) r

    IF ( SUM ( ( p(1:ndim) - center(1:ndim) )**2 ) <= r * r ) THEN
       inside = .TRUE.
    ELSE
       inside = .FALSE.
    END IF
    
    RETURN
  END SUBROUTINE sphere_imp_contains_point_3d

