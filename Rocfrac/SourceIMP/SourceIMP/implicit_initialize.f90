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

!!****
!!
!!  NAME
!!     implicit_initialize
!!
!!  FUNCTION
!!     Creates all of the data structures and matrices needed
!!     by the implicit solids solver
!!
!!  INPUTS
!!     none
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     MPI
!!
!!****

SUBROUTINE implicit_initialize(global)

  USE implicit_global
  USE ROCSTAR_RocFracComm

  IMPLICIT NONE

  TYPE(ROCFRAC_GLOBAL) :: global
  INTEGER :: i, j, m


!
! Let the code know it will need to initialize the acceleration array
!
  initAccel = .TRUE.

!
! Set up the global to local node mapping
!
  ALLOCATE(Global2Local(1:GNumNP))
  Global2Local(:) = -1
  DO i = 1, global%NumNP
     Global2Local(Local2Global(i)) = i
  ENDDO


!
! Initialize MPI communications
!
  IF(myid==0) PRINT*,'INITIALIZING COMMUNICATIONS'
  IF (nprocs > 1) THEN
     CALL InitComm1(global)
     CALL InitComm2(global)
  ENDIF
  

!
! Set up some constants
!  
  nstart_km = GNumNp
  DO m = 1, GNumNp
     DO i = 1, global%NumNp
        IF (Local2Global(i) == m) THEN
           IF (NodeProc(i) == myid) THEN
              nstart_km = MIN(nstart_km,Local2Global(i))
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  nstart_km = 3 * (nstart_km - 1) + 1
  nrows_km = 3*LNumNp


!
! Set up the nodal DOF vectors
!
  ALLOCATE(node_flag(1:global%NumNp,1:3))
  ALLOCATE(boundary_value(1:global%NumNp,1:3))
  node_flag(1:global%NumNp,1:3) = 0
  DO i = 1, global%NumNdsBC
     DO j = 1, 3
        IF ( global%BCFlag(j+1,i) == 0 ) THEN
           node_flag(global%BCFlag(1,i),j) = 8
           boundary_value(global%BCFlag(1,i),j) = global%BCvalue(j,i)
        ENDIF
        IF ( global%BCFlag(j+1,i) == 1 ) THEN
           node_flag(global%BCFlag(1,i),j) = 7
           boundary_value(global%BCFlag(1,i),j) = global%BCvalue(j,i)
        ENDIF
     ENDDO
  ENDDO


!
! Create the mass matrix
!
  CALL createM(global)


!
! Create the stiffness matrix
!
  CALL createK(global) 
  

END SUBROUTINE implicit_initialize

