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
MODULE IntegralCheck

CONTAINS
  
  SUBROUTINE CheckIntegral(glb, IntegralArray)

!!****f* Rocfrac/Rocfrac/Source/IntegralCheck.f90
!!
!!  NAME
!!    CheckIntegral
!!
!!  FUNCTION
!!    Sums the conservation quantities over all the
!!    processors. MPI_REDUCE
!!
!!  INPUTS
!!    glb -- global array
!!
!!  OUTPUT
!!    IntegralArray -- Conservation term:
!!                1. Volume - Deformed
!!                2. Mass
!!                3. x-momentum
!!                4. y-momentum
!!                5. z-momentum
!!                6. energy
!!                7. burning area
!!                8. non-burning area
!!                9. Volume - Undeformed
!!***
    USE ROCSTAR_RocFrac
  
    IMPLICIT NONE
  
    INCLUDE 'mpif.h'
    INCLUDE 'rocmanf90.h'
  
    TYPE(ROCFRAC_GLOBAL), POINTER :: glb
    INTEGER :: ierr
    REAL*8  :: IntegralArray(MAN_INTEG_SIZE)
    
    REAL*8, dimension(1:3) :: tmpSnd, tmpRcv

    tmpSnd(1)= glb%TotalGeomVolp
    tmpSnd(2)= glb%TotalGeomUndefVolp
    tmpSnd(3)= glb%TotalMassSolidp

    CALL MPI_REDUCE(tmpSnd,tmpRcv,3,MPI_DOUBLE_PRECISION, &
         MPI_SUM,0,glb%MPI_COMM_ROCFRAC,ierr)
    
    IntegralArray(MAN_INTEG_VOL) = tmpRcv(1)
    IntegralArray(MAN_INTEG_MASS) = tmpRcv(3) 
    IntegralArray(MAN_INTEG_XMOM) = 0.d0
    IntegralArray(MAN_INTEG_YMOM) = 0.d0
    IntegralArray(MAN_INTEG_ZMOM) = 0.d0
    IntegralArray(MAN_INTEG_ENER) = 0.d0
    IntegralArray(MAN_INTEG_IBAREA) = 0.d0
    IntegralArray(MAN_INTEG_INBAREA) = 0.d0
    IntegralArray(MAN_INTEG_VOL_UND) = tmpRcv(2)
    RETURN
    
  END SUBROUTINE  CheckIntegral
END MODULE IntegralCheck


