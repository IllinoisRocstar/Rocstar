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
SUBROUTINE max_dt_solid(dt_courant,glb)

!!****f* Rocfrac/Rocfrac/Source/max_dt_solid.f90
!!
!!  NAME
!!     max_dt_solid
!!
!!  FUNCTION
!!     Determines the maximum time step
!!
!!  INPUTS
!!     glb -- global array
!!
!!  OUTPUT
!!     dt_courant -- maximum time step
!!****

  USE ROCSTAR_RocFrac
  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  TYPE(ROCFRAC_GLOBAL) :: glb
  
  REAL*8 :: xx,yy,zz,size1,size2,size3,size4,size5,size6,size7,size8,size9,size10,size11,size12
  REAL*8 :: dhmin,dhminp,dt_courant, dt_courantp,dt_courantHT
  
  INTEGER :: i, ierr
  
  dhminp = 1000000000.d0

  IF(glb%iElType.EQ.8)THEN
     DO i = 1, glb%NumElVol
        call SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(1,i), glb%ElConnVol(4,i), size1)
        call SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(5,i), glb%ElConnVol(8,i), size2)
        call SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(6,i), glb%ElConnVol(7,i), size3)
        call SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(2,i), glb%ElConnVol(3,i), size4)
        
        call SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(4,i), glb%ElConnVol(8,i), size5)
        call SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(8,i), glb%ElConnVol(7,i), size6)
        call SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(7,i), glb%ElConnVol(3,i), size7)
        call SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(3,i), glb%ElConnVol(4,i), size8)
        
        call SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(1,i), glb%ElConnVol(5,i), size9)
        call SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(5,i), glb%ElConnVol(6,i), size10)
        call SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(6,i), glb%ElConnVol(2,i), size11)
        call SmallestElement(glb%NumNP, glb%MeshCoor, glb%ElConnVol(2,i), glb%ElConnVol(1,i), size12)
     
        dhminp = MIN(size1, size2, size3, size4, size5, size6, size7, size8, size9, size10, size11, size12, dhminp)
     enddo
     
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
     dhminp = MIN(size1,size2,size3,size4,size5,size6,dhminp)
     
  ENDDO
  ENDIF

!fix parallel, need to min across all processors

  
  CALL MPI_ALLREDUCE(dhminp,dhmin,1,MPI_DOUBLE_PRECISION, &
       MPI_MIN,glb%MPI_COMM_ROCFRAC,ierr)
  
  dt_courant = dhmin/glb%cd_fastest

  IF(glb%HeatTransSoln)THEN

     dt_courantHT = dhmin**2/glb%ThermalDiffusivity/6. ! is this 6., 3. or 2. ?

     IF(glb%ThermalExpansion)THEN
        
        dt_courant = MIN( dt_courant, dt_courantHT)

     ELSE
        dt_courant = dt_courantHT
     ENDIF

  Endif

  dt_courant = dt_courant*glb%CourantRatio

  
  RETURN
END SUBROUTINE max_dt_solid

