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
SUBROUTINE UpdateRbar(glb,Rnet)
  
  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(ROCFRAC_GLOBAL) :: glb

  REAL*8, ALLOCATABLE, DIMENSION(:) :: buf
  REAL*8, DIMENSION(1:3*glb%NumNP) :: Rnet

  INTEGER :: j, j1, k , k1, k2
  INTEGER :: ierr
            
  IF(glb%iElType.EQ.4)THEN 
     CALL v3d4_r_bar(glb%DispBar,Rnet, &
          glb%NumNP,glb%NumElVol,glb%ElConnVol,glb%MeshCoor,1,glb%NumElPartBndry)
  ELSE
     CALL V3D10_R_BAR(glb%DispBar,Rnet, &
          glb%NumNP,glb%NumElVol,glb%ElConnVol,glb%MeshCoor,1,glb%NumElPartBndry)
  ENDIF
!
!----- FORM THE BUFFER CONTAINING COMMUNICATED NODAL VALUES
!
  ALLOCATE(buf(1:glb%TotNumNdComm))
  k1 = 1
  DO j1 = 1, glb%TotNumNeighProcs
     k = glb%NeighProcList(j1)
     DO j = 1, glb%NumNdComm(j1)
        k2 = 3*glb%NdCommList(j1)%NdId(j)
        buf(k1)   = Rnet( k2 - 2 )
        buf(k1+1) = Rnet( k2 - 1 )
        buf(k1+2) = Rnet( k2 )
        k1 = k1 + 3
     ENDDO
  ENDDO
      
!     
!-MPI- RECEIVE THE INTERNAL FORCE VECTOR FROM THE NEIGHBORS
!
  DO j1 = 1, glb%TotNumNeighProcs
     k = glb%NeighProcList(j1)+1
     CALL MPI_IRECV(glb%RecvDataFrm(k)%rcvbuf(1), &
          glb%NumNdComm(j1)*3,MPI_DOUBLE_PRECISION,k-1,10,glb%MPI_COMM_ROCFRAC, &
          glb%ReqRcv(j1),ierr)
  ENDDO
!     
!-MPI- SEND THE INTERNAL FORCE VECTOR TO THE NEIGHBORS
!     
  k1 = 1
  DO j1 = 1, glb%TotNumNeighProcs
     k = glb%NeighProcList(j1)
     CALL MPI_ISEND(buf(k1),glb%NumNdComm(j1)*3,MPI_DOUBLE_PRECISION, &
          k,10,glb%MPI_COMM_ROCFRAC,glb%ReqSnd(j1),ierr)
     k1 = k1 + glb%NumNdComm(j1)*3
  ENDDO
!     
!----- CALCULATE THE INTERIOR SUBMESH'S INTERNAL FORCE VECTOR
!
  IF(glb%iElType.EQ.4)THEN
     CALL v3d4_r_bar(glb%DispBar,Rnet, &
          glb%NumNP,glb%NumElVol,glb%ElConnVol,glb%MeshCoor,glb%NumElPartBndry+1,glb%NumElVol)
  ELSE
     CALL V3D10_R_BAR(glb%DispBar,Rnet, &
          glb%NumNP,glb%NumElVol,glb%ElConnVol,glb%MeshCoor,glb%NumElPartBndry+1,glb%NumElVol)
  ENDIF
!
!-MPI- WAIT FOR INTERNAL FORCE VECTOR COMMUNICATION TO COMPLETE
!
  IF(glb%TotNumNeighProcs.GT.0)THEN
     CALL MPI_WAITALL(glb%TotNumNeighProcs,glb%ReqRcv,Glb%StatRcv,ierr)
     CALL MPI_WAITALL(glb%TotNumNeighProcs,glb%ReqSnd,glb%StatSnd,ierr)
  ENDIF
  DEALLOCATE(buf)
!
!----- ADD NEIGHBOR'S CONTRIBUTION TO THE INTERNAL FORCE VECTOR
!
  DO j1 = 1, glb%TotNumNeighProcs
     k = glb%NeighProcList(j1)+1
     k1 = 1
     DO j = 1, glb%NumNdComm(j1)
        k2 = ( glb%NdCommList(j1)%NdId(j) )*3
        Rnet(k2-2)= Rnet(k2-2) + glb%RecvDataFrm(k)%rcvbuf(k1)
        Rnet(k2-1)= Rnet(k2-1) + glb%RecvDataFrm(k)%rcvbuf(k1+1)
        Rnet(k2)  = Rnet(k2)   + glb%RecvDataFrm(k)%rcvbuf(k1+2)
        k1 = k1 + 3
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE UpdateRbar

