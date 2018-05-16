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
SUBROUTINE UpdateMassMatrix(glb)

!!****f* Rocfrac/Rocfrac/Source/UpdateMassMatrix.f90
!!
!!  NAME
!!     UpdateMassMatrix
!!
!!  FUNCTION
!!     Updates the mass matrix due to change in undeformed
!!     configuration. Calls the appropriate mass matrix subroutine
!!     and handles the parallel communication.
!!
!!  INPUTS
!!     glb -- global variables
!!
!!  USES
!!    ROCSTAR_RocFrac, V3D4_MASS, V3D4N_MASS, V3D10_MASS, V3D8_MASS
!!
!!****
 
  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  TYPE(ROCFRAC_GLOBAL) :: glb

  REAL*8, ALLOCATABLE, DIMENSION(:) :: buf

  INTEGER :: j, j1, k , k1
  INTEGER :: ierr
  INTEGER :: TotNumNdComm3
  INTEGER :: ElemStart, ElemEnd


  glb%xmass(:) = 0.d0

  IF(glb%HeatTransSoln) glb%CapctInv(:) = 0.d0
!  glb%TotalMassSolidp = 0.d0
!  glb%TotalGeomVolp = 0.d0

  IF(glb%iElType.EQ.4)THEN
        
     IF(.NOT.(glb%NdBasedEl))THEN

        ElemStart = 1
        DO j = 1, glb%NumMatVol
           ElemEnd = glb%NumElPartBndryMat(j) + ElemStart - 1
           CALL V3D4_MASS(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho,glb%xmass, &
                glb%NumNP,glb%NumElVol,glb%NumMatVol,ElemStart,ElemEnd)

           IF(glb%HeatTransSoln) CALL V3D4_capacitance(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho,glb%Cp,glb%CapctInv, &
                glb%NumNP,glb%NumElVol,glb%NumMatVol,ElemStart,ElemEnd)

           ElemStart = glb%NumElVolMat(j) + ElemStart
        ENDDO

     ELSE IF(glb%NdBasedEl)THEN
!        CALL VolNodal(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol, &
!             glb%NumNP,glb%NumElVol,glb%NumMatVol,1,glb%NumElVol,glb%TotalMassSolidp, &
!             glb%NumElNeigh,glb%ElConnNd,glb%AlphaR,glb%VolUndfmd)

        CALL V3D4N_MASS( &
             glb%MeshCoor,&
             glb%ElConnVol,&
             glb%MatIdVol,&
             glb%rho, &
             glb%xmass, &
             glb%NumNP,glb%NumElVol,glb%NumMatVol,1,glb%NumElVol,glb%TotalMassSolidp, &
             glb%NumElNeigh,glb%ElConnNd,glb%AlphaR,glb%VolUndfmd)

     ENDIF

  ELSE IF(glb%iElType.EQ.10)THEN
     ElemStart = 1
     DO j = 1, glb%NumMatVol
        ElemEnd = glb%NumElPartBndryMat(j) + ElemStart - 1
        CALL V3D10_MASS(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho,glb%xmass, &
             glb%NumNP,glb%NumElVol,glb%NumMatVol,ElemStart,ElemEnd,glb%TotalMassSolidp)
        IF(glb%HeatTransSoln) CALL V3D10_capacitance(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho,glb%Cp,glb%CapctInv, &
             glb%NumNP,glb%NumElVol,glb%NumMatVol,ElemStart,ElemEnd)

        ElemStart = glb%NumElVolMat(j) + ElemStart
     ENDDO

  ELSE IF(glb%iElType.EQ.8)THEN
     ElemStart = 1
     DO j = 1, glb%NumMatVol
        ElemEnd = glb%NumElPartBndryMat(j) + ElemStart - 1
        CALL V3D8_MASS(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho,glb%xmass, &
             glb%NumNP,glb%NumElVol,glb%NumMatVol,ElemStart,ElemEnd)
        ElemStart = glb%NumElVolMat(j) + ElemStart
     ENDDO
  ENDIF

  IF(.NOT.(glb%HeatTransSoln))THEN
!
!----- FORM THE BUFFER CONTAINING COMMUNICATED MASS MATRIX NODAL VALUES
!
     TotNumNdComm3 = glb%TotNumNdComm/3
     ALLOCATE(buf(1:TotNumNdComm3))
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

  ELSE  ! Pass the heat transfer capacitance to neighboring processors

!
!----- FORM THE BUFFER CONTAINING COMMUNICATED MASS MATRIX NODAL VALUES
!
     TotNumNdComm3 = ( glb%TotNumNdComm/3 )*2
     ALLOCATE(buf(1:TotNumNdComm3))
     k1 = 1
     DO j1 = 1, glb%TotNumNeighProcs
        k = glb%NeighProcList(j1)
        DO j = 1, glb%NumNdComm(j1)
           buf(k1) = glb%xmass( glb%NdCommList(j1)%NdId(j) )
           buf(k1+1) = glb%CapctInv( glb%NdCommList(j1)%NdId(j) )
           k1 = k1 + 2
        ENDDO
     ENDDO


!     
!-MPI- RECEIVE THE RECIPRICAL MASS MATRIX DIAGONAL FROM THE NEIGHBORS
!
     DO j1 = 1, glb%TotNumNeighProcs
        k = glb%NeighProcList(j1)+1
        CALL MPI_IRECV(glb%RecvDataFrm(k)%rcvbuf(1),glb%NumNdComm(j1)*2, &
             MPI_DOUBLE_PRECISION,k-1,10,glb%MPI_COMM_ROCFRAC,glb%ReqRcv(j1),ierr)
     ENDDO
!     
!-MPI- SEND THE RECIPRICAL MASS MATRIX DIAGONAL TO THE NEIGHBORS
!
     k1 = 1
     DO j1 = 1, glb%TotNumNeighProcs
        k = glb%NeighProcList(j1)
        CALL MPI_ISEND(buf(k1),glb%NumNdComm(j1)*2, &
             MPI_DOUBLE_PRECISION,k,10,glb%MPI_COMM_ROCFRAC,glb%ReqSnd(j1),ierr)
        k1 = k1 + glb%NumNdComm(j1)*2
     ENDDO

  ENDIF

!     
!----- CALCULATE THE INTERIOR SUBMESH'S RECIPRICAL MASS MATRIX DIAGONAL
!

  IF(glb%iElType.EQ.4)THEN
     IF(glb%iSolnType(1).LT.10)THEN
        ElemEnd = 0
        DO j = 1, glb%NumMatVol
           ElemStart =  ElemEnd + glb%NumElPartBndryMat(j) + 1
           ElemEnd = glb%NumElVolMat(j) + ElemEnd
           CALL V3D4_MASS(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho,glb%xmass, &
                glb%NumNP,glb%NumElVol,glb%NumMatVol,ElemStart,ElemEnd)

           IF(glb%HeatTransSoln) CALL V3D4_capacitance(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho,glb%Cp,glb%CapctInv, &
                glb%NumNP,glb%NumElVol,glb%NumMatVol,ElemStart,ElemEnd)
        ENDDO
     ENDIF
  ELSE IF(glb%iElType.EQ.10)THEN
     ElemEnd = 0
     DO j = 1, glb%NumMatVol
        ElemStart =  ElemEnd + glb%NumElPartBndryMat(j) + 1
        ElemEnd = glb%NumElVolMat(j) + ElemEnd
        CALL V3D10_MASS(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho,glb%xmass, &
             glb%NumNP,glb%NumElVol,glb%NumMatVol,ElemStart,ElemEnd,glb%TotalMassSolidp)

        IF(glb%HeatTransSoln) CALL V3D10_capacitance(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho,glb%Cp,glb%CapctInv, &
             glb%NumNP,glb%NumElVol,glb%NumMatVol,ElemStart,ElemEnd)

     ENDDO

  ELSE IF(glb%iElType.EQ.8)THEN
     ElemEnd = 0
     DO j = 1, glb%NumMatVol
        ElemStart =  ElemEnd + glb%NumElPartBndryMat(j) + 1
        ElemEnd = glb%NumElVolMat(j) + ElemEnd
        CALL V3D8_MASS(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho,glb%xmass, &
             glb%NumNP,glb%NumElVol,glb%NumMatVol,ElemStart,ElemEnd)
     ENDDO


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


  IF(.NOT.(glb%HeatTransSoln))THEN
     DO j1 = 1, glb%TotNumNeighProcs
        k = glb%NeighProcList(j1)+1
        DO j = 1, glb%NumNdComm(j1)
           k1 = glb%NdCommList(j1)%NdId(j)
           glb%xmass(k1) = glb%xmass(k1) + glb%RecvDataFrm(k)%rcvbuf(j)
        ENDDO
     ENDDO

  ELSE

     DO j1 = 1, glb%TotNumNeighProcs
        k = glb%NeighProcList(j1)+1
        DO j = 1, glb%NumNdComm(j1)
           k1 = glb%NdCommList(j1)%NdId(j)
           glb%xmass(k1) = glb%xmass(k1) + glb%RecvDataFrm(k)%rcvbuf(j*2-1)
           glb%CapctInv(k1) = glb%CapctInv(k1) + glb%RecvDataFrm(k)%rcvbuf(j*2)
        ENDDO
     ENDDO

  ENDIF
  glb%xmass(:) = 1.d0/glb%xmass(:)

  IF(glb%HeatTransSoln) glb%CapctInv(:) = 1./glb%CapctInv(:)


  IF(glb%iSolnType(1).GE.10)THEN
     
     ALLOCATE(buf(1:TotNumNdComm3))
     k1 = 1
     DO j1 = 1, glb%TotNumNeighProcs
        k = glb%NeighProcList(j1)
        DO j = 1, glb%NumNdComm(j1)
           buf(k1) = glb%VolUndfmd( glb%NdCommList(j1)%NdId(j) )
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
           glb%VolUndfmd(k1) = glb%VolUndfmd(k1) + glb%RecvDataFrm(k)%rcvbuf(j)
        ENDDO
     ENDDO
  ENDIF

  RETURN 
END SUBROUTINE UpdateMassMatrix

