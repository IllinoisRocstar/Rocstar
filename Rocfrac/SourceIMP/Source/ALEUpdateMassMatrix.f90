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
SUBROUTINE ALEUpdateMassMatrix(glb)

!!****f* Rocfrac/Rocfrac/Source/ALEUpdateMassMatrix.f90
!!
!!  NAME
!!    ALEUpdateMassMatrix
!!
!!  FUNCTION
!!    Updates the nodal coordinates as a result
!!    of the regressing boundaries. Calculates
!!    the new lumped inverse mass matrix from the new
!!    nodal coordinates, MPI calls handle 
!!    communication between partition boundaries
!!
!!  USED BY
!!    RocfracMain
!!  
!!  USES
!!    ROCSTAR_RocFrac -- Global variables
!!    V3D4_MASS, V3D4N_MASS, V3D10_MASS
!!
!!  INPUTS
!!    glb -- global array
!!
!!  OUTPUTS
!!    xmass -- Inverse lumped mass matrix
!!
!!***

  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER :: j, j1, k , k1
  INTEGER :: ierr
  INTEGER :: TotNumNdComm3
  
!-- UPDATE MESH POSITION
          
  DO j = 1, glb%NumNP
     j1 = j*3
     glb%MeshCoor(1,j) = glb%coor(1,j) + glb%DispBar(j1 - 2)
     glb%MeshCoor(2,j) = glb%coor(2,j) + glb%DispBar(j1 - 1)
     glb%MeshCoor(3,j) = glb%coor(3,j) + glb%DispBar(j1    )
  ENDDO

!-- UPDATE MASS MATRIX due to change in undeformed configuration
!        ***Warning ! initialize xm(j) first

  glb%xmass(:) = 0.d0
  glb%TotalMassSolidp = 0.d0
  IF(glb%iElType.EQ.4)THEN
     IF( .NOT.(glb%NdBasedEl) )THEN
        CALL V3D4_MASS(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho,glb%xmass, &
             glb%NumNP,glb%NumElVol,glb%NumMatVol,1,glb%NumElPartBndry,glb%TotalMassSolidp)
     ELSE IF(glb%NdBasedEl)THEN
        CALL V3D4N_MASS(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho,glb%xmass, &
             glb%NumNP,glb%NumElVol,glb%NumMatVol,1,glb%NumElVol,glb%TotalMassSolidp, &
             glb%NumElNeigh,glb%ElConnNd,glb%AlphaR,glb%VolUndfmd)
     ENDIF
  ELSE
     CALL V3D10_MASS(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho,glb%xmass, &
          glb%NumNP,glb%NumElVol,glb%NumMatVol,1,glb%NumElPartBndry,glb%TotalMassSolidp)
  ENDIF
!
!----- FORM THE BUFFER CONTAINING COMMUNICATED MASS MATRIX NODAL VALUES
!
  TotNumNdComm3 = glb%TotNumNdComm/3
  ALLOCATE(buf(1:glb%TotNumNdComm3))
  k1 = 1
  DO j1 = 1, glb%TotNumNeighProcs
     k = glb%NeighProcList(j1)
     DO j = 1, glb%NumNdComm(j1)
        buf(k1) = glb%xmass( glb%NdCommList(j1)%NdId(j) )
        k1 = k1 + 1
     ENDDO
  ENDDO
!     
!-MPI- RECEIVE THE RECIPRICAL MASS MATRIX DIAGONAL FROM THE NEIGHBORS
!
  DO j1 = 1, glb%TotNumNeighProcs
     k = glb%NeighProcList(j1)+1
     CALL MPI_IRECV(glb%RecvDataFrm(k)%rcvbuf(1),glb%NumNdComm(j1), &
          MPI_DOUBLE_PRECISION,k-1,10,glb%MPI_COMM_ROCFRAC,glb%ReqRcv(j1),ierr)
  ENDDO
!     
!-MPI- SEND THE RECIPRICAL MASS MATRIX DIAGONAL TO THE NEIGHBORS
!
  k1 = 1
  DO j1 = 1, glb%TotNumNeighProcs
     k = glb%NeighProcList(j1)
     CALL MPI_ISEND(buf(k1),glb%NumNdComm(j1), &
          MPI_DOUBLE_PRECISION,k,10,glb%MPI_COMM_ROCFRAC,glb%ReqSnd(j1),ierr)
     k1 = k1 + glb%NumNdComm(j1)
  ENDDO
!     
!----- CALCULATE THE INTERIOR SUBMESH'S RECIPRICAL MASS MATRIX DIAGONAL
!
  IF(glb%iElType.EQ.4)THEN
     IF(glb%iSolnType.NE.10)THEN
        CALL V3D4_MASS(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho,glb%xmass, &
             glb%NumNP,glb%NumElVol,glb%NumMatVol,glb%NumElPartBndry+1,glb%NumElVol,glb%TotalMassSolidp)
     ENDIF
  ELSE
     CALL V3D10_MASS(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho,glb%xmass, &
          glb%NumNP,glb%NumElVol,glb%NumMatVol,glb%NumElPartBndry+1,glb%NumElVol,glb%TotalMassSolidp)
  ENDIF
!
!-MPI- WAIT FOR THE RECIPRICAL MASS MATRIX DIAGONAL COMMUNICATION 
!      TO COMPLETE
!
  IF(glb%TotNumNeighProcs.GT.0)THEN
     CALL MPI_WAITALL(glb%TotNumNeighProcs,glb%ReqRcv,Glb%StatRcv,ierr)
     CALL MPI_WAITALL(glb%TotNumNeighProcs,glb%ReqSnd,glb%StatSnd,ierr)
  ENDIF
  DEALLOCATE(buf)
!
!----- ADD NEIGHBOR'S CONTRIBUTION TO THE RECIPRICAL MASS 
!      MATRIX DIAGONAL
!     
  DO j1 = 1, glb%TotNumNeighProcs
     k = glb%NeighProcList(j1)+1
     DO j = 1, glb%NumNdComm(j1)
        k1 = glb%NdCommList(j1)%NdId(j)
        glb%xmass(k1) = glb%xmass(k1) + glb%RecvDataFrm(k)%rcvbuf(j)
     ENDDO
  ENDDO
          
  glb%xmass(:) = 1.d0/glb%xmass(:)

  RETURN 
  END

