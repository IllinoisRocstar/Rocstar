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
MODULE RocFracMain

!
!!****h* Rocfrac/Rocfrac/Source/RocFracMain
!!
!!  NAME
!!    RocFracMain.f90
!!
!!  FUNCTION
!!    3D Dynamic  Explicit Code  with  ALE  formulation for regressing 
!!    boundaries Finite Element Analysis Code with additional fracture
!!    simulation using cohesive elements.
!!
!!  USAGE
!!    Finite  Element code to solve the 3-Dimensional TRANSIENT structural problem
!!
!!  USES
!!    RocFracSubInterface, UpdateStructuralSoln, feminp, VolRatio, vol_elem_mat,
!!    RocFracInterfaceInitial, RocFracInterfaceBuff,  UpdateMassMatrix, V3D4_volume,
!!    max_dt_solid, UpdateRbar,v3d4_ale,V3D10_ALE,FluidPressLoad, TractPressLoad,
!!    UpdateStructural,principal_stress,bc_enforce
!!    
!!    Global variables stored in modules :  ROCSTAR_RocFrac,ROCSTAR_RocFracComm,ROCSTAR_RocFracInterp
!!
!!  COPYRIGHT
!!
!!    University of Illinois, Urbana-Champaign, (C) 2003
!!
!!  AUTHOR
!!
!!  principal : M.S. Breitenfeld, P.H. Geubelle
!!    - email : brtnfld@uiuc.edu, geubelle@uiuc.edu
!!
!!  contributing : Changyu Huang, Amit Acharya
!!
!!
!!  CREATION DATE
!!       2001
!!
!!
!!*****

  
  USE ROCSTAR_RocFrac 
  USE ROCSTAR_RocFracComm
  USE ROCSTAR_RocFracInterp
  USE RocFracSubInterface
  USE UpdateStructuralSoln
  USE implicit_global

  IMPLICIT NONE 

  PUBLIC :: RocFracFinalize, RocFracSoln, RocFracInitialize

CONTAINS
!
! ------------------------------------------------------------------ RocFracInitialize

  SUBROUTINE RocFracInitialize( glb, InitialTime, MPI_COMM_ROCSTAR, MAN_init, surfIn, volIn, obtain_attr)

    USE implicit_global

    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INCLUDE 'roccomf90.h'
    
    TYPE(ROCFRAC_GLOBAL),POINTER  :: glb
    REAL*8, INTENT(IN)            :: InitialTime
    INTEGER, INTENT(IN)           :: MPI_COMM_ROCSTAR
    INTEGER, INTENT(IN)           :: MAN_init, obtain_attr
    CHARACTER(*), INTENT(IN)      :: surfIn, volIn
    
    INTEGER :: i,j,j1,jj,k,k1,k2,idum,iaux,iaux1
    REAL*8 :: aux1,aux
    INTEGER :: NdBCflag

    INTEGER :: MyId, NumProcs, ierr
    
    INTEGER :: id2d
    INTEGER :: iunit2
    INTEGER :: kk
    INTEGER, DIMENSION(3) :: ndsurf
    REAL*8, ALLOCATABLE, DIMENSION(:) :: buf
    INTEGER :: IdPacket
    REAL*8 :: mag

    REAL*8 :: IOversion

    INTEGER :: ios

    CHARACTER(len=2) :: chr2

    CHARACTER*120 :: fracHDFFname, meshFile

    character(LEN=1), POINTER, DIMENSION(:) :: names
    character(LEN=3) :: ichr03
    character(LEN=4) :: ichr04

    integer :: NumElTypes2D 

    character(LEN=4) :: ChrElType
    integer :: endPt, startPt, chrlngth

 
    INTEGER :: NumParComm
    INTEGER, pointer, dimension(:) :: ParComm

    INTEGER, pointer :: ArrayTmp
    INTEGER, pointer, dimension(:) :: ArrayTmp1
    integer :: iptr
    integer :: icnt1, icnt2, icnt3, icnt4

    integer :: iprocs

    ! ... For RocTherm validation, triangular Temp distribution, COstoich 03/18/09
    REAL*8 :: SliceTemp, zcoord
    INTEGER :: dir

    glb%MPI_COMM_ROCFRAC = MPI_COMM_ROCSTAR
    rocstar_communicator = MPI_COMM_ROCSTAR

    CALL MPI_COMM_RANK(glb%MPI_COMM_ROCFRAC,MyId,ierr)
    CALL MPI_COMM_SIZE(glb%MPI_COMM_ROCFRAC,NumProcs,ierr)
    WRITE(glb%MyIdChr,'(i4.4)') MyId



!
! -- Read the control deck file

    IF(MyId.EQ.0) PRINT*,'RocFrac :: Reading input deck'

    glb%io_input = 10
    glb%io_sum   = 11
    
    CALL feminp(glb, MyId)

    glb%InterfaceSFNumNodes = 0
    glb%InterfaceSFnbNumNodes = 0

    PRINT*,'RocFrac :: ....Done Reading input deck',myid
    CALL MPI_BARRIER(glb%MPI_COMM_ROCFRAC,ierr)


! -- Read in roc_face data structure for extracted 2D mesh.
!    mapnode maps the surface nodes to the volume mesh's nodes

    IF(glb%iElType.EQ.4)THEN
       glb%iStrGss = 1
       glb%LwrBnd = 1
       glb%UppBnd = 3
    ELSE IF(glb%iElType.EQ.10)THEN
       glb%iStrGss = 4
       glb%LwrBnd = 4
       glb%UppBnd = 6
    ELSE IF(glb%iElType.EQ.8)THEN
       glb%iStrGss = 1
       glb%LwrBnd = 1
       glb%UppBnd = 4
    ENDIF

!   glb%InterfaceSVbar = 0.d0
    
   IF(myid.EQ.0) PRINT*,'RocFrac :: Reading Mesh'

! - Start 
! - Process the 3D mesh

! Number of Nodes
   CALL COM_new_window( volWin)
 
   CALL COM_get_size( VolIn//".nc", MyId+1, glb%NumNP)
   CALL COM_set_size( volWin//'.nc', MyId+1, glb%NumNP   )

!   ALLOCATE(glb%coor(1:3,1:glb%NumNP))
   ALLOCATE(glb%MeshCoor(1:3,1:glb%NumNP))
   ALLOCATE(glb%xmass(1:glb%NumNP))
   IF(glb%HeatTransSoln) ALLOCATE(glb%CapctInv(1:glb%NumNP))

!!$          
!!$          ! Nodal Coordinates
!!$          DO i=1,glb%NumNP
!!$             ! NodeID, X, Y, Z ! , DummyFlag1
!!$             READ(14,*) j,glb%coor(1,i),glb%coor(2,i),glb%coor(3,i) ! ,iaux
!!$             
!!$             glb%Coor(1,j) =  glb%Coor(1,j)+ mag
!!$             glb%Coor(3,j) =  glb%Coor(3,j)+ mag
!!$             glb%MeshCoor(1:3,i) = glb%coor(1:3,i)
!!$          ENDDO
!!$

   CALL COM_get_size( VolIn//".bcnode", MyId+1,glb%NumNdsBCcrypt)

   ALLOCATE(glb%BCFlagCrypt(1:2,1:glb%NumNdsBCcrypt)) ! should this be half NumNdsBCdrypt

   ALLOCATE(glb%BCValueGlb(1:glb%NumNdsBCcrypt*6))



!!$          ! _________________________________________
!!$          ! ___ Read Structural Nodal Boundary Flags
!!$          ! _________________________________________
!!$
!!$          READ(14,*) glb%NumNdsBC ! , iaux
!!$          ALLOCATE(glb%BCFlag(1:4,1:glb%NumNdsBC),glb%BCvalue(1:3,1:glb%NumNdsBC))
!!$          DO i = 1, glb%NumNdsBC
!!$             READ(14,*) glb%BCFlag(1,i),NdBCflag !, iaux
!!$             glb%BCFlag(2,i) = glb%bcCond(NdBCflag)%BCtypeX
!!$             glb%BCFlag(3,i) = glb%bcCond(NdBCflag)%BCtypeY
!!$             glb%BCFlag(4,i) = glb%bcCond(NdBCflag)%BCtypeZ
!!$             glb%BCvalue(1,i) = glb%bcCond(NdBCflag)%BCvalueX
!!$             glb%BCvalue(2,i) = glb%bcCond(NdBCflag)%BCvalueY
!!$             glb%BCvalue(3,i) = glb%bcCond(NdBCflag)%BCvalueZ            
!!$          ENDDO
!!$
!!$!          DEALLOCATE(glb%bcCond)
!!$
!!$       CASE(4)
!!$
!!$          ! Packet 04
!!$          ! __________________________________________
!!$          ! ___ Read Nodal Mesh Motion Boundary Flags
!!$          ! __________________________________________
!!$
!!$
!!$          READ(14,*) glb%NumNdsBCmm !, iaux
!!$          IF(glb%NumNdsBCmm.NE.0) ALLOCATE(glb%BCFlagmm(1:4,1:glb%NumNdsBCmm),glb%BCvaluemm(1:3,1:glb%NumNdsBCmm))
!!$
!!$          IF(glb%ALEenabled)THEN
!!$
!!$             DO i=1,glb%NumNdsBCmm
!!$                READ(14,*) glb%BCFlagmm(1,i) ,NdBCflag !, iaux
!!$                glb%BCFlagmm(2,i) = glb%bcCondmm(NdBCflag)%BCtypeX
!!$                glb%BCFlagmm(3,i) = glb%bcCondmm(NdBCflag)%BCtypeY
!!$                glb%BCFlagmm(4,i) = glb%bcCondmm(NdBCflag)%BCtypeZ
!!$                glb%BCvaluemm(1,i) = glb%bcCondmm(NdBCflag)%BCvalueX
!!$                glb%BCvaluemm(2,i) = glb%bcCondmm(NdBCflag)%BCvalueY
!!$                glb%BCvaluemm(3,i) = glb%bcCondmm(NdBCflag)%BCvalueZ 
!!$             ENDDO
!!$          ELSE
!!$
!!$             DO i=1,glb%NumNdsBCmm
!!$                READ(14,'()')
!!$             ENDDO
!!$
!!$          ENDIF
!!$             
!!$          
!!$!          IF(glb%NumNdsBCmm.NE.0) DEALLOCATE(glb%bcCondmm)
!!$ 
         
!!$          ! __________________________________________
!!$          ! ___ Read Element Data
!!$          ! __________________________________________

    CALL COM_get_connectivities(VolIn,MyId+1,NumElTypes2D,names)

    startPt = 1
    DO i = 1, NumElTypes2D
      ! Search for the next attribute name
      endPt = startPt
      chrlngth = 0
      DO WHILE (endPt .LE. UBOUND(names,1))
         IF (names(endPt) .NE. ' ') THEN
            chrlngth = chrlngth + 1
            ChrElType(chrlngth:chrlngth) = names(endPt)
            endPt = endPt + 1
         ELSE
            EXIT
         END IF
      END DO 

      startPt = endPt + 1

      IF(ChrElType(1:chrlngth).EQ.':T10')THEN
         CALL COM_get_size( VolIn//".:T10", MyId+1, glb%NumElVol)
         CALL COM_set_size( volWin//'.:T10', MyId+1, glb%NumElVol)
      ELSE IF(ChrElType(1:chrlngth).EQ.':T4')THEN
         CALL COM_get_size( VolIn//".:T4", MyId+1, glb%NumElVol)
         CALL COM_set_size( volWin//'.:T4', MyId+1, glb%NumElVol)
      ELSE IF(ChrElType(1:chrlngth).EQ.':H8')THEN
         CALL COM_get_size( VolIn//".:H8", MyId+1, glb%NumElVol)
         CALL COM_set_size( volWin//'.:H8', MyId+1, glb%NumElVol)
      ELSE
         PRINT*,'ROCFRAC ERROR: Volume mesh type element not supported'
         PRINT*,'Read in Element Type :: ', ChrElType(1:chrlngth)
         CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)
      ENDIF

   END DO

   CALL COM_free_buffer(names)

   CALL COM_get_array_const(VolIn//".IMP",MyId+1,ArrayTmp)   
   IF ( ArrayTmp == 1 ) THEN
      glb%IMP = .TRUE.
   ELSE 
      glb%IMP = .FALSE.
   END IF
   
   IF ( glb%IMP ) THEN
      
      ALLOCATE(NodeProc(1:glb%NumNP))
      ALLOCATE(Local2Global(1:glb%NumNp))
      
      CALL COM_get_array_const(VolIn//".NumNPLocal",MyId+1,ArrayTmp)
      LNumNP = ArrayTmp
      
      CALL COM_get_array_const(VolIn//".NumNPGlobal",MyId+1,ArrayTmp) 
      GNumNP = ArrayTmp
      
      CALL COM_get_array_const(VolIn//".NodeNumGlobal",MyId+1,ArrayTmp1)
      Local2Global = ArrayTmp1
      
      CALL COM_get_array_const(VolIn//".NodeProc",MyId+1,ArrayTmp1)
      NodeProc = ArrayTmp1
      
      myid = MyId
      nprocs = NumProcs
      
   END IF

!!$   CALL COM_call_function( obtain_attr, 2, &
!!$        COM_get_attribute_handle_const( VolIn//".NumElPartBndry"), &
!!$        COM_get_attribute_handle( VolIn//".NumElPartBndry"))
!   stop
   CALL COM_get_array_const(VolIn//".NumElPartBndry",MyId+1,ArrayTmp)

   glb%NumElPartBndry = ArrayTmp


!!$   CALL COM_call_function( obtain_attr, 2, &
!!$        COM_get_attribute_handle_const( VolIn//".NumElVolMat"), &
!!$        COM_get_attribute_handle( VolIn//".NumElVolMat"))
   CALL COM_get_array_const(VolIn//".NumElVolMat",MyId+1,ArrayTmp1)

   ALLOCATE(glb%NumElVolMat(1:glb%NumMatVol))

   DO i = 1, glb%NumMatVol
      glb%NumElVolMat(i) = ArrayTmp1(i)
   ENDDO
   CALL COM_get_array_const(VolIn//".NumElPartBndryMat",MyId+1,ArrayTmp1)

   ALLOCATE(glb%NumElPartBndryMat(1:glb%NumMatVol))
   DO i = 1, glb%NumMatVol
      glb%NumElPartBndryMat(i) = ArrayTmp1(i)
   ENDDO

   IF(.NOT.(glb%DebondPart).AND..NOT.(glb%DebondPart_Matous)) THEN
      ALLOCATE(glb%ci(1:9,1:glb%NumMatVol))
      ALLOCATE(glb%cj(1:9,1:glb%NumMatVol))
      ALLOCATE(glb%ci_full(1:6,1:6,1:glb%NumMatVol))
   ENDIF

   i = 0
          
   DO kk = 1, glb%NumMatVol
      IF(kk.EQ.1)THEN
         ALLOCATE(glb%MatIdVol(1:glb%NumElVol))
         ALLOCATE(glb%ElConnVol(1:glb%iElType,1:glb%NumElVol))
      ENDIF
!
! -- Read element data.
       
      IF(glb%NdBasedEl)THEN
         IF(kk.EQ.1)THEN
                   
            ALLOCATE(glb%NumElNeigh(1:glb%NumNP))
            ALLOCATE(glb%ElConnNd(1:glb%NumNP,1:40))
            ALLOCATE(glb%AlphaR(1:4,1:glb%NumElVol))
            ALLOCATE(glb%VolUndfmd(1:glb%NumNP))
            
            glb%NumElNeigh(:) = 0
            glb%ElConnNd(:,:) = 0 ! fix not work for more then one material ??
         ENDIF
      ENDIF
      DO jj = 1,glb%NumElVolMat(kk)
         i = i + 1
!!$         IF ( glb%iElType==4) THEN
!!$            READ(14,*) glb%MatIdVol(i),glb%ElConnVol(1,i), &
!!$                 glb%ElConnVol(2,i), glb%ElConnVol(3,i), &
!!$                 glb%ElConnVol(4,i) !, iaux, iaux
!!$         ELSE IF( glb%iElType==10) THEN
!!$            READ(14,*) glb%MatIdVol(i),glb%ElConnVol(1,i), &
!!$                 glb%ElConnVol(2,i), glb%ElConnVol(3,i), &
!!$                 glb%ElConnVol(4,i), glb%ElConnVol(5,i), &
!!$                 glb%ElConnVol(6,i), glb%ElConnVol(7,i), &
!!$                 glb%ElConnVol(8,i), glb%ElConnVol(9,i), &
!!$                 glb%ElConnVol(10,i) !, iaux, iaux
!!$         ELSE IF( glb%iElType==8) THEN
!!$            READ(14,*) glb%MatIdVol(i),glb%ElConnVol(1,i), &
!!$                 glb%ElConnVol(2,i), glb%ElConnVol(3,i), &
!!$                 glb%ElConnVol(4,i), glb%ElConnVol(5,i), &
!!$                 glb%ElConnVol(6,i), glb%ElConnVol(7,i), &
!!$                 glb%ElConnVol(8,i) !, iaux, iaux
!!$         END IF
                
         IF(glb%NdBasedEl)THEN
                   
            CALL VolRatio(glb%ElConnVol(1,i),glb%ElConnVol(2,i),glb%ElConnVol(3,i),glb%ElConnVol(4,i),glb%AlphaR(1:4,i),&
                 glb%MeshCoor,glb%NumNp,glb%NdMassLump)
                   
            DO j = 1, 4
               glb%NumElNeigh(glb%ElConnVol(j,i)) = glb%NumElNeigh(glb%ElConnVol(j,i)) + 1
               glb%ElConnNd(glb%ElConnVol(j,i),glb%NumElNeigh(glb%ElConnVol(j,i))) = i
            ENDDO
         ENDIF
      ENDDO
   ENDDO
! __________________________________________
! _______  Parallel Communication
! __________________________________________

     IF(NumProcs.NE.0)THEN



!        CALL COM_set_size( volWin//'.pconn', ip, NumNeighProcs*2+MaxNumNodesComm)    
!        CALL COM_set_array(volWin//'.pconn', ip, Pconn_Comm, 1)


!!$        CALL COM_new_attribute( volWin//'.NodesToCommunicate', 'p', COM_INTEGER, 1, '')
!!$        CALL COM_set_size( volWin//'.NodesToCommunicate', ip, MaxNumNodesComm)    
!!$        CALL COM_set_array(volWin//'.NodesToCommunicate', ip, NodesToCommunicate, 1)
!!$        
!!$
!!$        CALL COM_new_attribute( volWin//'.ID_sendto_List', 'p', COM_INTEGER, 1, '')
!!$        CALL COM_set_size( volWin//'.ID_sendto_List', ip, NumNeighProcs)    
!!$        CALL COM_set_array(volWin//'.ID_sendto_List', ip, ID_sendto_List, 1)
!!$        
!!$        CALL COM_new_attribute( volWin//'.NumNeighProcs_List', 'p', COM_INTEGER, 1, '')
!!$        CALL COM_set_size( volWin//'.NumNeighProcs_List', ip, NumNeighProcs)    
!!$        CALL COM_set_array(volWin//'.NumNeighProcs_List', ip, NumNeighProcs_List, 1)

 
     ELSE
        glb%TotNumNeighProcs = 1
     ENDIF

     

     CALL COM_get_size( VolIn//".pconn", MyId+1, NumParComm)
     CALL COM_set_size( volWin//'.pconn', MyId+1, NumParComm)    
     CALL COM_allocate_array(volWin//'.pconn', MyId+1, ParComm, 1)

!!$!
!!$! -- Read lst nodes that need to be sent to other processors
!!$!
!!$          READ(14,*) glb%TotNumNeighProcs
!!$          
!!$!          IF(glb%TotNumNeighProcs.NE.0)THEN
!!$             ALLOCATE(glb%NeighProcList(1:glb%TotNumNeighProcs))
!!$!          ENDIF
!!$          
!!$          ALLOCATE(glb%NumNdComm(1:glb%TotNumNeighProcs))
!!$          ALLOCATE(glb%NdCommList(1:glb%TotNumNeighProcs))
!!$          glb%TotNumNdComm = 0
!!$
!!$! NOTE: NeighProcList is the actual processor starting at id = 0, not 1
!!$
!!$          DO i = 1, glb%TotNumNeighProcs
!!$             READ(14,*) glb%NeighProcList(i),glb%NumNdComm(i)
!!$             !       print*,glb%NeighProcList(i),glb%NumNdComm(i)
!!$             
!!$             glb%TotNumNdComm = glb%TotNumNdComm + glb%NumNdComm(i)
!!$             
!!$             ALLOCATE(glb%NdCommList(i)%NdId(1:glb%NumNdComm(i)))
!!$             DO j=1,glb%NumNdComm(i)
!!$                READ(14,*) glb%NdCommList(i)%NdId(j)
!!$             ENDDO
!!$          ENDDO
!!$          glb%TotNumNdComm = glb%TotNumNdComm*3 ! x3 (x,y,z R_in)
!!$
!!$       CASE(7)
!!$
!!$          ! Packet 07
!!$          ! __________________________________________
!!$          ! ___ Nodal History 
!!$          ! __________________________________________
!!$
!!$          READ(14,*) glb%NumNodeIO
!!$          ALLOCATE(glb%NodeIO(1:glb%NumNodeIO))
!!$          DO i = 1, glb%NumNodeIO
!!$             READ(14,*) glb%NodeIO(i)
!!$             WRITE(chr2,'(i2.2)') i
!!$             OPEN(33+i,FILE='Rocfrac/NdHistory.'//glb%MyIdChr//'.'//chr2)
!!$             WRITE(33+i,*) '#', glb%coor(1:3,glb%NodeIO(i))
!!$          ENDDO
!!$
!!$          CASE(8)
!!$
!!$          ! Packet 08
!!$          ! __________________________________________
!!$          ! ___ Nodal History 
!!$          ! __________________________________________  
!!$
!!$             READ(14,*) glb%NumNdsBCHT !, iaux
!!$             IF(glb%NumNdsBCHT.NE.0) ALLOCATE(glb%BCFlagHT(1:2,1:glb%NumNdsBCHT),glb%BCvalueHT(1,1:glb%NumNdsBCHT))
!!$             DO i=1,glb%NumNdsBCHT
!!$                READ(14,*) glb%BCFlagHT(1,i) ,NdBCflag !, iaux
!!$                glb%BCFlagHT(2,i) = glb%bcCondHT(NdBCflag)%BCtypeX
!!$                glb%BCvalueHT(1,i) = glb%bcCondHT(NdBCflag)%BCvalueX 
!!$             ENDDO
!!$       CASE(99)
!!$
!!$          ! Packet 99
!!$          ! _____________________________
!!$          ! ___ End of Input File
!!$          ! _____________________________
!!$
!!$          EXIT
!!$          
!!$       CASE DEFAULT
!!$          PRINT*, 'ERROR: ROCFRAC'
!!$          PRINT*, 'Input Packet = ',IdPacket, 'Is Not A Valid Option'
!!$          PRINT*, 'Check {prefix}.####.inp files and try again'
!!$
!!$          CALL MPI_BARRIER(glb%MPI_COMM_ROCFRAC,ierr)
!!$          CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)
!!$          STOP
!!$       END SELECT
!!$
!!$    ENDDO
!!$
!!$    CLOSE(14)
    
    IF(MyId.EQ.0) PRINT*,'RocFrac :: Finished Reading Solids Mesh'

    IF(.NOT.(glb%DebondPart).AND. .NOT.(glb%DebondPart_Matous) )THEN

       IF ( glb%NumMatOrtho == 0 ) THEN
          ALLOCATE(glb%E11o(1:0))
          ALLOCATE(glb%E22o(1:0))
          ALLOCATE(glb%E33o(1:0))
          ALLOCATE(glb%xnu12o(1:0))
          ALLOCATE(glb%xnu13o(1:0))
          ALLOCATE(glb%xnu23o(1:0))
          ALLOCATE(glb%G12o(1:0))
          ALLOCATE(glb%G13o(1:0))
          ALLOCATE(glb%G23o(1:0))
          ALLOCATE(glb%vx1o(1:0))
          ALLOCATE(glb%vy1o(1:0))
          ALLOCATE(glb%vz1o(1:0))
          ALLOCATE(glb%vx2o(1:0))
          ALLOCATE(glb%vy2o(1:0))
          ALLOCATE(glb%vz2o(1:0))
          ALLOCATE(glb%vx3o(1:0))
          ALLOCATE(glb%vy3o(1:0))
          ALLOCATE(glb%vz3o(1:0))
       END IF

       CALL VOL_ELEM_MAT(glb%E,glb%xnu,glb%ci,glb%cj,glb%NumMatVol,glb%iElIntgratn,glb%MatOrtho)
       CALL VOL_ELEM_MAT_ORTHO( glb%ci_full, glb%ci, glb%NumMatVol, glb%NumMatOrtho, glb%MatOrtho, &
             glb%E11o, glb%E22o, glb%E33o, glb%xnu12o, glb%xnu13o, glb%xnu23o, &
             glb%G12o, glb%G13o, glb%G23o, glb%vx1o, glb%vy1o, glb%vz1o, &
             glb%vx2o, glb%vy2o, glb%vz2o, glb%vx3o, glb%vy3o, glb%vz3o )

    ELSE IF(glb%DebondPart_Matous)THEN

       ALLOCATE( glb%StrainOld(1:4,1:glb%NumElVol*6) ) 
       ALLOCATE( glb%SoftParam(1:4,1:glb%NumElVol) ) 
       ALLOCATE( glb%cd(1:4,1:glb%NumElVol) )

       CALL VOL_ELEM_MAT_MATOUS(glb)

    ELSE IF(glb%DebondPart)THEN

       ALLOCATE(glb%STATEV_Part1(1:glb%NumElVol))
       ALLOCATE(glb%STATEV_Part2(1:glb%NumElVol))

       glb%STATEV_Part1(:) = 1
       glb%STATEV_Part2(:) = 1

       ALLOCATE(glb%StrainTrace(1:glb%NumElVol))
    ENDIF

    
    IF(glb%iElType.EQ.8)THEN

       allocate(glb%mixed_map(1:8,1:9,1:12))
       allocate(glb%enhanced_map(1:8,1:9,1:9))
       allocate(glb%Aenh(1:9,1:glb%NumElVol))

       glb%mixed_map = 0.d0
       glb%enhanced_map = 0.d0
       glb%Aenh(:,:) = 0.d0
            
       CALL enhanced_elem_maps_hex(glb%mixed_map,glb%enhanced_map)

       ! ... If thermal solver, dmat is NumMatVol x 3 x 3.
       IF (glb%HeatTransSoln) THEN
          allocate(glb%dmat(1:glb%NumMatVol,1:3,1:3))
          CALL ConductivityTensor(glb%NumMatVol,glb%KappaHT,glb%dmat)
       ELSE
          allocate(glb%dmat(1:glb%NumMatVol,1:9,1:9))
          
          DO i = 1, glb%NumMatVol
             CALL get_mat_stiffness(glb%E(i),glb%xnu(i),glb%dmat(i,:,:))
          ENDDO
       ENDIF
    ENDIF

!-----Setting initial conditions
    ALLOCATE(glb%Disp(1:3*glb%NumNP))
    glb%Disp(:)  = 0.d0

    ALLOCATE(glb%Accel(1:3*glb%NumNP))
    glb%Accel(:) = 0.d0

    ALLOCATE(glb%DispOld(1:3*glb%NumNP))
    ALLOCATE(glb%S11(1:glb%iStrGss,1:glb%NumElVol),glb%S22(1:glb%iStrGss,1:glb%NumElVol))
    ALLOCATE(glb%S33(1:glb%iStrGss,1:glb%NumElVol) )
    ALLOCATE(glb%S12(1:glb%iStrGss,1:glb%NumElVol),glb%S23(1:glb%iStrGss,1:glb%NumElVol))
    ALLOCATE(glb%S13(1:glb%iStrGss,1:glb%NumElVol) )
    ALLOCATE(glb%SVonMises(1:glb%NumElVol))

    IF(glb%HeatTransSoln)THEN
       ALLOCATE(glb%Temperature(1:glb%NumNP))
       glb%Temperature(1:glb%NumNP) = glb%Temperature0
    ENDIF

    IF(glb%ArtificialDamping)THEN
       ALLOCATE(glb%DetF_Old(1:glb%iStrGss,1:glb%NumElVol))
       glb%DetF_Old(:,:) = 1.d0
    endif

! new: ale formulation
    ALLOCATE(glb%DispBar(1:3*glb%NumNP),glb%VeloBar(1:3*glb%NumNP),glb%AccelBar(1:3*glb%NumNP))
    glb%DispBar(:) = 0.d0
    glb%VeloBar(:) = 0.d0
    glb%AccelBar(:) = 0.d0
    

    ALLOCATE(glb%VeloHalf(1:3*glb%NumNP),glb%VeloBarOld(1:3*glb%NumNP))
    glb%VeloHalf(:) = 0.d0
    glb%VeloBarOld(:) = 0.d0

    ALLOCATE(glb%DispTotal(1:3*glb%NumNP))

    glb%S11(1:glb%iStrGss,1:glb%NumElVol) = 0.d0
    glb%S22(1:glb%iStrGss,1:glb%NumElVol) = 0.d0
    glb%S33(1:glb%iStrGss,1:glb%NumElVol) = 0.d0
    glb%S12(1:glb%iStrGss,1:glb%NumElVol) = 0.d0
    glb%S23(1:glb%iStrGss,1:glb%NumElVol) = 0.d0
    glb%S13(1:glb%iStrGss,1:glb%NumElVol) = 0.d0
    glb%SVonMises(1:glb%NumElVol) = 0.d0
    glb%Accel = 0.d0
    glb%VeloBar = 0.d0
!    glb%AccelBndry = 0.d0 ! not zero for restart.
!    glb%VeloBndry = 0.d0 ! not zero for restart.
    
    glb%DispBar = 0.d0
    glb%DispTotal = 0.d0

! Create window for HDF output

    CALL COM_new_attribute( volWin//'.disp', 'n', COM_DOUBLE, 3, 'm')
    CALL COM_new_attribute( volWin//'.disp_burn', 'n',COM_DOUBLE, 3, 'm')
    CALL COM_new_attribute( volWin//'.velo', 'n', COM_DOUBLE, 3, 'm/s')
    CALL COM_new_attribute( volWin//'.stress', 'e', COM_DOUBLE, 1, 'Pa')
    CALL COM_new_attribute( volWin//'.accel', 'n', COM_DOUBLE, 3, 'm/s^2')
    CALL COM_new_attribute( volWin//'.vbar', 'n', COM_DOUBLE, 3, 'm/s')
    CALL COM_new_attribute( volWin//'.S11', 'e', COM_DOUBLE, glb%iStrGss, 'Pa')
    CALL COM_new_attribute( volWin//'.S22', 'e', COM_DOUBLE, glb%iStrGss, 'Pa')
    CALL COM_new_attribute( volWin//'.S33', 'e', COM_DOUBLE, glb%iStrGss, 'Pa')
    CALL COM_new_attribute( volWin//'.S12', 'e', COM_DOUBLE, glb%iStrGss, 'Pa')
    CALL COM_new_attribute( volWin//'.S23', 'e', COM_DOUBLE, glb%iStrGss, 'Pa')
    CALL COM_new_attribute( volWin//'.S13', 'e', COM_DOUBLE, glb%iStrGss, 'Pa')
    
    CALL COM_new_attribute( volWin//'.NumElPartBndry', 'p', COM_INTEGER, 1, '')
    CALL COM_new_attribute( volWin//'.NumElVolMat', 'p', COM_INTEGER, 1, '')
    CALL COM_new_attribute( volWin//'.NumElPartBndryMat', 'p', COM_INTEGER, 1, '')

    IF(glb%DebondPart)THEN
       CALL COM_new_attribute( volWin//'.StrainTrace', 'e', COM_DOUBLE, 1, ' ')
       CALL COM_new_attribute( volWin//'.DebondLg', 'e', COM_DOUBLE, 1, ' ')
       CALL COM_new_attribute( volWin//'.DebondSm', 'e', COM_DOUBLE, 1, ' ')
    ENDIF

    IF(glb%DebondPart_Matous)THEN
       CALL COM_new_attribute( volWin//'.StrainOld', 'e', COM_DOUBLE, 4, ' ')
       CALL COM_new_attribute( volWin//'.SoftParam', 'e', COM_DOUBLE, 4, ' ')
    ENDIF

    IF(glb%HeatTransSoln) CALL COM_new_attribute( volWin//'.Temp', 'n', COM_DOUBLE, 1, 'K')

    IF ( glb%NumNP > 0) THEN
!
! Registering Coordinates
!

!!$    CALL COM_init_mesh( volWin//'.nc', MyId+1, glb%MeshCoor, glb%NumNP)

       CALL COM_set_array(volWin//'.nc', MyId+1, glb%MeshCoor,3)
!
! Registering Element Connectivity
!
       IF(glb%iElType.EQ.4)THEN

!!$       CALL COM_init_mesh( volWin//'.T4', MyId+1, glb%ElConnVol, glb%NumElVol)

          CALL COM_set_array( volWin//'.:T4', MyId+1,  glb%ElConnVol,4)

       ELSE IF(glb%iElType.EQ.10)THEN

!!$       CALL COM_init_mesh( volWin//'.T10', MyId+1, glb%ElConnVol, glb%NumElVol)

          CALL COM_set_array( volWin//'.:T10', MyId+1,  glb%ElConnVol,10)

       ELSE IF(glb%iElType.EQ.8)THEN
!!$       CALL COM_init_mesh( volWin//'.H8', MyId+1, glb%ElConnVol, glb%NumElVol)

          CALL COM_set_array( volWin//'.:H8', MyId+1,  glb%ElConnVol,8)

       ENDIF
       
       CALL COM_set_array( volWin//'.disp', MyId+1, glb%Disp,3) ! was DispTotal
       CALL COM_set_array( volWin//'.disp_burn', MyId+1, glb%DispBar,3)
       CALL COM_set_array( volWin//'.velo', MyId+1, glb%VeloHalf,3)
       CALL COM_set_array( volWin//'.stress', MyId+1, glb%SVonMises,1)
       CALL COM_set_array( volWin//'.accel', MyId+1, glb%Accel,3)
       CALL COM_set_array( volWin//'.vbar', MyId+1, glb%VeloBar,3)
       CALL COM_set_array( volWin//'.S11', MyId+1, glb%S11)
       CALL COM_set_array( volWin//'.S22', MyId+1, glb%S22)
       CALL COM_set_array( volWin//'.S33', MyId+1, glb%S33)
       CALL COM_set_array( volWin//'.S12', MyId+1, glb%S12)
       CALL COM_set_array( volWin//'.S23', MyId+1, glb%S23)
       CALL COM_set_array( volWin//'.S13', MyId+1, glb%S13)

       CALL COM_set_size( volWin//'.NumElPartBndry', MyId+1, 1)
       CALL COM_set_array( volWin//'.NumElPartBndry', MyId+1,glb%NumElPartBndry,1 )

       CALL COM_set_size( volWin//'.NumElVolMat', MyId+1, glb%NumMatVol)
       CALL COM_set_array( volWin//'.NumElVolMat', MyId+1, glb%NumElVolMat, 1)
       
       CALL COM_set_size( volWin//'.NumElPartBndryMat', MyId+1, glb%NumMatVol)
       CALL COM_set_array( volWin//'.NumElPartBndryMat',MyId+1, glb%NumElPartBndryMat, 1 )

       IF(glb%HeatTransSoln) CALL COM_set_array( volWin//'.Temp', MyId+1, glb%Temperature,1)
       IF(glb%DebondPart) CALL COM_set_array( volWin//'.DebondLg', MyId+1, glb%STATEV_Part1,1)
       IF(glb%DebondPart) CALL COM_set_array( volWin//'.DebondSm', MyId+1, glb%STATEV_Part2,1)
       IF(glb%DebondPart) CALL COM_set_array( volWin//'.StrainTrace', MyId+1, glb%StrainTrace,1)
       IF(glb%DebondPart_Matous) THEN
          CALL COM_set_array( volWin//'.StrainOld', MyId+1, glb%StrainOld,4)
          CALL COM_set_array( volWin//'.SoftParam', MyId+1, glb%SoftParam,4)
       ENDIF

    ENDIF


    IF(InitialTime.NE.0.d0)THEN
       PRINT*,'RocFrac :: RESTARTING, SOLIDS'
       glb%ReStart = .TRUE.
    ENDIF
   CALL COM_new_attribute(volwin//'.BCValue','p',COM_DOUBLE, 1, '')
   CALL COM_set_size( volwin//'.BCValue', MyId+1, glb%NumNdsBCcrypt*6)
   CALL COM_set_array( volwin//'.BCValue', MyId+1, glb%BCValueGlb, 1)

   

   CALL COM_new_attribute(volwin//'.bcnode','p',COM_INTEGER, 2, '')
   CALL COM_set_size( volwin//'.bcnode', MyId+1, glb%NumNdsBCcrypt)
   CALL COM_set_array( volwin//'.bcnode', MyId+1, glb%BCFlagCrypt, 2)
   
   CALL COM_new_attribute(volwin//'.MatType','e',COM_INTEGER, 1, '')
   CALL COM_set_array( volwin//'.MatType', MyId+1, glb%MatIdVol, 1)
   
   glb%NumNdsBC   = 0
   glb%NumNdsBCmm = 0
   glb%NumNdsBCHT = 0   

   CALL COM_window_init_done( volWin)

   CALL COM_call_function( obtain_attr, 2, &
        COM_get_attribute_handle_const( VolIn//".all"), &
        COM_get_attribute_handle( VolWin//".all"))





   glb%TotNumNeighProcs = 0
 
!
! -- Read lst nodes that need to be sent to other processors
!
!!$   iptr = 1
!!$   DO
!!$      IF(iptr.GT.NumParComm)EXIT
!!$      glb%TotNumNeighProcs = glb%TotNumNeighProcs + 1
!!$      iptr = iptr + ParComm( iptr + 1) + 2
!!$   ENDDO

    ! crashes on copper if no if statement
    IF(NumProcs.GT.1) glb%TotNumNeighProcs = ParComm(1)

          
!          IF(glb%TotNumNeighProcs.NE.0)THEN
   ALLOCATE(glb%NeighProcList(1:glb%TotNumNeighProcs))
!          ENDIF
 
     
   ALLOCATE(glb%NumNdComm(1:glb%TotNumNeighProcs))
   ALLOCATE(glb%NdCommList(1:glb%TotNumNeighProcs))
   glb%TotNumNdComm = 0


! NOTE: NeighProcList is the actual processor starting at id = 0, not 1
  
   iptr = 1
   DO i = 1, glb%TotNumNeighProcs
      iptr = iptr + 1
      glb%NeighProcList(i) = ParComm(iptr) -1
      iptr = iptr + 1
      glb%NumNdComm(i) = ParComm(iptr)
      glb%TotNumNdComm = glb%TotNumNdComm + glb%NumNdComm(i)
             
      ALLOCATE(glb%NdCommList(i)%NdId(1:glb%NumNdComm(i)))

      DO j=1,glb%NumNdComm(i)
         iptr = iptr + 1
         glb%NdCommList(i)%NdId(j) = ParComm(iptr)
      ENDDO

   ENDDO


   glb%TotNumNdComm = glb%TotNumNdComm*3 ! x3 (x,y,z R_in)

    IF(glb%TotNumNeighProcs.NE.0)THEN
       ALLOCATE(glb%RecvDataFrm(1:NumProcs) )
       ALLOCATE(glb%ReqRcv     (1:glb%TotNumNeighProcs)  )
       ALLOCATE(glb%ReqSnd     (1:glb%TotNumNeighProcs)  )
       
       ALLOCATE(glb%StatSnd  (1:MPI_STATUS_SIZE,1:glb%TotNumNeighProcs)  )
       ALLOCATE(glb%StatRcv  (1:MPI_STATUS_SIZE,1:glb%TotNumNeighProcs)  )
    ENDIF

    IF(glb%HeatTransSoln)THEN
       DO j = 1, glb%TotNumNeighProcs
          k = glb%NeighProcList(j) + 1
          ALLOCATE(glb%RecvDataFrm(k)%rcvbuf(1:glb%NumNdComm(j)*4))
       ENDDO
    ELSE
       DO j = 1, glb%TotNumNeighProcs
          k = glb%NeighProcList(j) + 1
          ALLOCATE(glb%RecvDataFrm(k)%rcvbuf(1:glb%NumNdComm(j)*3))
       ENDDO
    ENDIF

   
   DO i = 1, glb%NumNdsBCcrypt
      NdBCflag = MOD(glb%BCFlagCrypt(2,i),100)
! structural boundary conditions
      IF(NdBCflag.GT.0) glb%NumNdsBC = glb%NumNdsBC + 1

      NdBCflag = glb%BCFlagCrypt(2,i)/10000
! mesh motion boundary conditions
      IF(NdBCflag.GT.0) glb%NumNdsBCmm = glb%NumNdsBCmm  + 1

      NdBCflag = MOD(glb%BCFlagCrypt(2,i),10000)/100
! thermal boundary conditions
      IF(NdBCflag.GT.0) glb%NumNdsBCHT = glb%NumNdsBCHT + 1
   ENDDO

   IF(glb%NumNdsBC.NE.0)THEN
      ALLOCATE(glb%BCFlag(1:4,1:glb%NumNdsBC),glb%BCvalue(1:3,1:glb%NumNdsBC))
      ALLOCATE(glb%AccelBndry(1:3*glb%NumNdsBC),glb%VeloBndry(1:3*glb%NumNdsBC))
   ENDIF
   IF(glb%NumNdsBCmm.NE.0) ALLOCATE(glb%BCFlagmm(1:4,1:glb%NumNdsBCmm),glb%BCvaluemm(1:3,1:glb%NumNdsBCmm))
   IF(glb%NumNdsBCHT.NE.0) ALLOCATE(glb%BCFlagHT(1:2,1:glb%NumNdsBCHT),glb%BCvalueHT(1,1:glb%NumNdsBCHT))

   icnt1 = 0
   icnt2 = 0
   icnt3 = 0
   icnt4 = 1
   
!!$   IF ( glb%NumNdsBC > 0) THEN ! boundary conditions
!!$
!!$      CALL COM_new_attribute( volWin//'.velobndry', 'p', COM_DOUBLE, 3, 'm/s')
!!$      CALL COM_new_attribute( volWin//'.accbndry', 'p', COM_DOUBLE, 3, 'm/s^2')
!!$      
!!$      CALL COM_set_size( volWin//'.velobndry', MyId+1, glb%NumNdsBC  )
!!$      CALL COM_set_array( volWin//'.velobndry', MyId+1, glb%VeloBndry, 3)
!!$      
!!$      CALL COM_set_size( volWin//'.accbndry', MyId+1, glb%NumNdsBC  )
!!$      CALL COM_set_array( volWin//'.accbndry', MyId+1, glb%AccelBndry, 3)
!!$      
!!$   ENDIF   



   DO i = 1, glb%NumNdsBCcrypt


      NdBCflag = MOD(glb%BCFlagCrypt(2,i),100)



! structural boundary conditions
      IF(NdBCflag.GT.0) THEN
         icnt1 = icnt1 + 1
         glb%BCFlag(1,icnt1) = glb%BCFlagCrypt(1,i)

         glb%BCFlag(2,icnt1) = glb%bcCond(NdBCflag)%BCtypeX
         glb%BCFlag(3,icnt1) = glb%bcCond(NdBCflag)%BCtypeY
         glb%BCFlag(4,icnt1) = glb%bcCond(NdBCflag)%BCtypeZ
         glb%BCvalue(1,icnt1) = glb%bcCond(NdBCflag)%BCvalueX
         glb%BCvalue(2,icnt1) = glb%bcCond(NdBCflag)%BCvalueY
         glb%BCvalue(3,icnt1) = glb%bcCond(NdBCflag)%BCvalueZ
         
         glb%AccelBndry(icnt1*3-2:icnt1*3) = glb%BCValueGlb(icnt4:icnt4+2)
         icnt4 = icnt4 + 3
         glb%VeloBndry(icnt1*3-2:icnt1*3) = glb%BCValueGlb(icnt4:icnt4+2)
         icnt4 = icnt4 + 3
         
      ENDIF

      NdBCflag = glb%BCFlagCrypt(2,i)/10000



! mesh motion boundary conditions
      IF(NdBCflag.GT.0) THEN

         icnt2 = icnt2 + 1

         glb%BCFlagmm(1,icnt2) = glb%BCFlagCrypt(1,i)
      

         glb%BCFlagmm(2,icnt2) = glb%bcCondmm(NdBCflag)%BCtypeX
         glb%BCFlagmm(3,icnt2) = glb%bcCondmm(NdBCflag)%BCtypeY
         glb%BCFlagmm(4,icnt2) = glb%bcCondmm(NdBCflag)%BCtypeZ
         glb%BCvaluemm(1,icnt2) = glb%bcCondmm(NdBCflag)%BCvalueX
         glb%BCvaluemm(2,icnt2) = glb%bcCondmm(NdBCflag)%BCvalueY
         glb%BCvaluemm(3,icnt2) = glb%bcCondmm(NdBCflag)%BCvalueZ
      ENDIF

      NdBCflag = MOD(glb%BCFlagCrypt(2,i),10000)/100
! thermal boundary conditions

      IF(NdBCflag.GT.0) THEN

         icnt3 = icnt3 + 1

         glb%BCFlagHT(1,icnt3) = glb%BCFlagCrypt(1,i)
         
         glb%BCFlagHT(2,icnt3) = glb%bcCondHT(NdBCflag)%BCtypeX
         glb%BCvalueHT(1,icnt3) = glb%bcCondHT(NdBCflag)%BCvalueX
      ENDIF
   ENDDO


   CALL RocFracInterfaceInitial( glb,obtain_attr,surfIn)



!!    ! ... Triangular initial temperature distribution in 
!!    ! ... 1-D bar validation problem. COstoich 03/18/09
!!    ! ... Change the direction (MeshCoor(#,i)) depending on orientation of bar

!!    ! ... Initialize zcoord ... 
!!    zcoord = -1.0
!!    ! ... Bar orientation (ex: 1 is in the x-dir)
!!    dir = 2

!!   IF(glb%HeatTransSoln)THEN
!!      DO i = 1, glb%NumNP
!!         if (glb%MeshCoor(dir,i) > zcoord) zcoord = glb%MeshCoor(dir,i)
!!         print*,'zcoord,i',zcoord,i,maxval(glb%MeshCoor(dir,:)),minval(glb%MeshCoor(dir,:))
!!         if ((zcoord - minval(glb%MeshCoor(dir,:)))<=((maxval(glb%MeshCoor(dir,:))-minval(glb%MeshCoor(dir,:)))/2)) THEN
!!            SliceTemp =  0 + 2 * glb%Temperature0 / (maxval(glb%MeshCoor(dir,:))-minval(glb%MeshCoor(dir,:))) * (zcoord - minval(glb%MeshCoor(dir,:)))
!!         else
!!            SliceTemp =  glb%Temperature0 - 2 * glb%Temperature0 / (maxval(glb%MeshCoor(dir,:))-minval(glb%MeshCoor(dir,:))) * (zcoord - maxval(glb%MeshCoor(dir,:))/2 - minval(glb%MeshCoor(dir,:))/2)
           
!!         end if
!!         print*,'SliceTemp'
!!         do j=1,glb%NumNP
!!            print*,'z',glb%MeshCoor(dir,j),'zcoord',zcoord
!!         end do
!!         do j = 1, glb%NumNP
!!            if ((glb%MeshCoor(dir,j)<(zcoord+.001)).and.(glb%MeshCoor(dir,j)>(zcoord-.001))) THEN
!!               glb%Temperature(j) = Slicetemp 
              
!!            end if
!!         end do
!!         print*,myid,glb%Temperature
!!         read*,j
!!      end DO
!!   ENDIF


!!   ! ... Smooth initial temperature distribution in 
!!   ! ... 1-D coupled problem. COstoich 04/02/09
  
!!    zcoord = -1.0
!!    dir=3

!!   IF(glb%HeatTransSoln)THEN
!!      DO i = 1, glb%NumNP
!!         if (glb%MeshCoor(dir,i) > zcoord) zcoord = glb%MeshCoor(dir,i)
!! !        print*,'zcoord,i',zcoord,i,maxval(glb%MeshCoor(dir,:)),minval(glb%MeshCoor(dir,:))
!! !        if ((zcoord - minval(glb%MeshCoor(2,:)))<=((maxval(glb%MeshCoor(dir,:))-minval(glb%MeshCoor(dir,:)))/2)) THEN
!!            SliceTemp =  glb%Temperature0 +  (330.0-glb%Temperature0) / (maxval(glb%MeshCoor(dir,:))-minval(glb%MeshCoor(dir,:))) * (zcoord - minval(glb%MeshCoor(dir,:)))
!! !        else
!! !           SliceTemp =  glb%Temperature0 - 2 * glb%Temperature0 / (maxval(glb%MeshCoor(2,:))-minval(glb%MeshCoor(2,:))) * (zcoord - maxval(glb%MeshCoor(2,:))/2 - minval(glb%MeshCoor(2,:))/2)
           
!! !        end if
!! !        print*,'SliceTemp'
!!         do j=1,glb%NumNP
!! !           print*,'z',glb%MeshCoor(dir,j),'zcoord',zcoord
!!         end do
!!         do j = 1, glb%NumNP
!!            if ((glb%MeshCoor(dir,j)<(zcoord+.001)).and.(glb%MeshCoor(dir,j)>(zcoord-.001))) THEN
!!               glb%Temperature(j) = Slicetemp 
              
!!            end if
!!         end do
!! !        print*,glb%Temperature
!! !        read*,j
!!      end DO
!!   ENDIF


!!$   DO j = 1, glb%NumNP
!!$      j1 = j*3
!!$      glb%MeshCoor(1,j) = glb%coor(1,j) + glb%DispBar(j1 - 2)
!!$      glb%MeshCoor(2,j) = glb%coor(2,j) + glb%DispBar(j1 - 1)
!!$      glb%MeshCoor(3,j) = glb%coor(3,j) + glb%DispBar(j1    )
!!$   ENDDO

   CALL RocFracInterfaceBuff(glb)

   CALL COM_call_function( MAN_init, 2, surWin, volWin)
   
   IF(.NOT.(glb%ALEenabled))THEN
      CALL UpdateMassMatrix(glb)
      
   ENDIF
11 FORMAT(A,'_',A,A1)

!!$    IF(myid.EQ.0) PRINT*,'ROCFRAC :: Finished RocFracInitialize'
!!$   CALL MPI_BARRIER(glb%MPI_COMM_ROCFRAC,ierr)
!!$   CALL MPI_FINALIZE(glb%MPI_COMM_ROCFRAC,ierr)

! -- Mass/Volume Conservation

   glb%TotalMassSolidp = 0.d0
   glb%TotalGeomVolp = 0.d0
   
   IF(glb%iElType.EQ.4.OR.glb%iElType.EQ.10 .AND.(.NOT.(glb%NdBasedEl)))THEN
      CALL V3D4_volume(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho, &
           glb%NumNP,glb%NumElVol,glb%NumMatVol,glb%Disp,1,glb%NumElVol,&
           glb%TotalMassSolidp,glb%TotalGeomVolp,glb%TotalGeomUndefVolp, &
           glb%iElType)
   ENDIF
   IF(glb%NumProbesNd.NE.0)THEN
      CALL FindProbe(glb,myid)
   ENDIF
   glb%iAmpCnt = 1
   
   glb%prop = 0.d0
   glb%slope = 0.d0



!!$  PRINT*,'sdfsdflkjsdflkj',myid
!!$  PRINT*,'finished readsdv',myid
!!$  CALL MPI_BARRIER(glb%MPI_COMM_ROCFRAC,i)
!!$  stop

   inquire(FILE='Rocfrac/Rocin/OverlayMappings.txt',EXIST=glb%OverlayExist)

   IF(glb%OverlayExist)THEN
      CALL readsdv(glb,myid)

      OPEN(456,FILE ='Rocfrac/Rocin/OverlayMappings.txt')
      DO

         READ(456,*,IOSTAT=ios) iprocs 
         IF(ios.LT.0) EXIT
         
         IF(iprocs.EQ.myid)THEN
            
!!$        READ(456,*) glb%nf1
!!$        ALLOCATE(MapFaceEl2Vol1a(1:600),FaceOfVolEL1a(1:600))
!!$        
!!$        DO i = 1, glb%nf1
!!$           READ(456,'()') 
!!$        ENDDO
!!$        
!!$        READ(456,*) iaux
!!$        
!!$        DO i = 1, iaux
!!$           READ(456,'()')
!!$        ENDDO
!!$        PRINT*,';;',myid
            
!!$        READ(456,*) glb%nf1
!!$        PRINT*,glb%nf1,myid
            
            READ(456,*) glb%nf1
            ALLOCATE(glb%MapFaceEl2Vol1(1:glb%nf1),glb%FaceOfVolEL1(1:glb%nf1))
            
            DO i = 1, glb%nf1
               READ(456,*) glb%MapFaceEl2Vol1(i), glb%FaceOfVolEL1(i)
            ENDDO
            
            READ(456,*) glb%nf2
            
            ALLOCATE(glb%MapFaceEl2Vol2(1:glb%nf2),glb%FaceOfVolEL2(1:glb%nf2))
            
            DO i = 1, glb%nf2
               READ(456,*) glb%MapFaceEl2Vol2(i), glb%FaceOfVolEL2(i)
            ENDDO
!!$
!!$        READ(456,*) iaux
!!$        
!!$        DO i = 1, iaux
!!$           READ(456,'()')
!!$        ENDDO
!!$
!!$
!!$        DO i = 1, glb%nf1
!!$           READ(456,'()')
!!$        ENDDO
        EXIT
        
     ELSE
        READ(456,*) iaux
        
        DO i = 1, iaux
           READ(456,'()') 
        ENDDO
        
        READ(456,*) iaux
        
        DO i = 1, iaux
           READ(456,'()')
        ENDDO
     ENDIF
  ENDDO
  close(456)

endif

!!$  IF(myid.EQ.1)THEN
!!$     DO i = 1, glb%nf1
!!$        PRINT*, i,glb%MapFaceEl2Vol1(i), glb%FaceOfVolEL1(i)
!!$     ENDDO
!!$     
!!$     DO i = 1, glb%nf2
!!$        PRINT*,i, glb%MapFaceEl2Vol2(i), glb%FaceOfVolEL2(i)
!!$     ENDDO
!!$  ENDIF

!!$   do i = 1, glb%NumNp
!!$      IF( glb%MeshCoor(1,i).LE.1.0001.AND.glb%MeshCoor(2,i).GT.0.02499.AND.glb%MeshCoor(3,i).GT.0.02499)THEN
!!$         glb%NumNodeIO = i
!!$      ENDIF
!!$   enddo



    IF ( glb%IMP .eqv. .true.) THEN
       
       ! ... For development of thermal code, must change in the future
       ! ... at present, if the thermal code is run, then the structural 
       ! ... code is not.
       IF ( glb%HeatTransSoln .eqv. .true.) THEN
          CALL thermal_initialize(glb)
       ELSE
          CALL implicit_initialize(glb)
       END IF
    END IF

   
 END SUBROUTINE RocFracInitialize


!
! ----------------------------------------------------------------------- RocFracFinalize

  SUBROUTINE RocFracFinalize( glb)

    USE implicit_global

    IMPLICIT NONE
    INCLUDE 'roccomf90.h'
     
    TYPE(ROCFRAC_GLOBAL), POINTER :: glb

    CALL COM_delete_window( surWin)
    CALL COM_delete_window( volWin)

    IF ( glb%IMP ) THEN
       IF(glb%HeatTransSoln .EQV. .TRUE.) THEN

          CALL thermal_finalize(glb)

       ELSE

          CALL implicit_finalize(glb)

       END IF

    END IF
    
  END SUBROUTINE RocFracFinalize
!
! ----------------------------------------------------------------------- RocFracSoln

  SUBROUTINE RocFracSoln( glb, CurrentTime, CurrentTimeStep, MAN_update_inbuff)
    
    IMPLICIT NONE
    INCLUDE 'mpif.h'

    TYPE(ROCFRAC_GLOBAL), POINTER :: glb
    REAL*8, INTENT(IN)            :: CurrentTime, CurrentTimeStep
    INTEGER, INTENT(IN)           :: MAN_update_inbuff

    INTEGER, SAVE:: istep

    
    INTEGER j,jj,jjj,k,k1,k2,j1,i,j2,j3,k3,k4 ! loop counters & dummy variables
    
    INTEGER :: ntime_out   ! processor's id ( 0-(nprocs-1) )
!--   First section Rnet = R_bar
!--   Force Sum: cohesive traction + (-internal) + (-damping)+ external 
!--   i.e. R_co - R_in - R_damp + R_ex
    REAL*8, ALLOCATABLE, DIMENSION(:) :: Rnet
    REAL*8, ALLOCATABLE, DIMENSION(:) :: RnetHT
!--   reciprical of mass matrix diagonal
    REAL*8 :: tvol, tvol_com
    REAL*8, DIMENSION(2) :: trry_dt,trry_vol,trry_calc,trry_dv_com
    REAL*8, DIMENSION(2) :: trry_energy,trry_io,trry_coh,trry_vol_com
    REAL*8, ALLOCATABLE, DIMENSION(:) :: buf


! Dummy variable
    REAL*8 :: a1
    CHARACTER*8 :: ichr8   ! output file number (a character)
    CHARACTER*3 :: ai,ai1  ! output file number (a character)
    INTEGER :: MyId, ierr  ! mpi error variable
    
    INTEGER :: isubstep,nsubstep
    REAL*8 :: dt_solid_sub, alpha
    REAL*8 :: CurrentTimestepSolid
    
    REAL*8 :: prin1, prin2, prin3
    
    INTEGER, SAVE :: ifirststep

    REAL*8 :: maxvalue
    INTEGER :: maxlocat
    INTEGER :: jaux
    INTEGER :: ElemStart, ElemEnd, NumProcs
    REAL*8 :: DispPrev
    
    INTEGER :: NdBCflag, icnt1, icnt5


    logical :: debug
    CHARACTER*4 :: ichr1, ichr2

    REAL*8 :: tempf1, tempf2
! _____________________________
! _____ SOLUTION LOOP
! _____________________________
         
    CALL MPI_COMM_RANK(glb%MPI_COMM_ROCFRAC,MyId,ierr)    
    CALL MPI_COMM_SIZE(glb%MPI_COMM_ROCFRAC,NumProcs,ierr)

    debug = glb%debug_state

    tvol = 0.d0; tvol_com = 0.d0
    nsubstep = 1
    
    CurrentTimestepSolid = CurrentTimestep
    ifirststep = 0
    
    glb%DispOld(:) = glb%Disp(:)
    ALLOCATE(Rnet(1:3*glb%NumNP))
    IF(glb%HeatTransSoln) ALLOCATE(RnetHT(1:glb%NumNP))

    IF(MyId.EQ.0)THEN
       PRINT*,'RocFrac ::  Time Step             Dt'
       PRINT*,'RocFrac :: -------------------------'
    ENDIF

    glb%CurrTime = CurrentTime

!    open(3001,file='ndhistory')
   
    DO WHILE (nsubstep.EQ.1)

       CALL max_dt_solid(dt_solid_sub,glb)

       IF(dt_solid_sub.GE.CurrentTimestepSolid)THEN
          nsubstep = 0
          glb%DT = CurrentTimestepSolid
          CurrentTimestepSolid = 0
       ELSE
          nsubstep = 1
          glb%DT = dt_solid_sub
          CurrentTimestepSolid = CurrentTimestepSolid - glb%DT
       ENDIF

       glb%DTInv = 1.d0 / glb%DT

       istep = istep + 1

       ! Compute alpha here
       alpha = 1.d0 - CurrentTimestepSolid / CurrentTimeStep
       CALL COM_call_function( MAN_update_inbuff, 1, alpha)

       glb%CurrTime = CurrentTime + (CurrentTimestep-CurrentTimestepSolid)

       IF(MyId.EQ.0) WRITE(*,'(a10,i10,4e12.4)')'RocFrac ::', istep, &
            (CurrentTimestep-CurrentTimestepSolid), glb%DT,CurrentTimestep,glb%CurrTime


       !-- (0) INITIALIZE


       ntime_out = 0
       Rnet(:) = 0.d0
       IF(glb%HeatTransSoln) RnetHT(:) = 0.d0


       IF( glb%AmplitudeTable )THEN

          IF(glb%iAmpCnt.LE.glb%NumEntries-1)THEN
             IF(glb%CurrTime.GT.glb%AmpTable(1,glb%iAmpCnt+1))THEN
                glb%iAmpCnt = glb%iAmpCnt + 1
             ENDIF
          ENDIF
          glb%slope = glb%AmpTable(2,glb%iAmpCnt)
          glb%prop = glb%CurrTime*glb%AmpTable(2,glb%iAmpCnt) + glb%AmpTable(3,glb%iAmpCnt)

       ENDIF


       IF(glb%ALEenabled)THEN  

          ! -- Enforce Mesh Velocity BC's of Mesh Motion Given by the fluids.

          DO j = 1, glb%InterfaceSFNumNodes
             k3 = 3*glb%MapNodeSF(j)
             k2 = k3 - 1
             k1 = k3 - 2
             glb%VeloBar(k1) = glb%InterfaceSFVbar(1,j)
             glb%VeloBar(k2) = glb%InterfaceSFVbar(2,j)
             glb%VeloBar(k3) = glb%InterfaceSFVbar(3,j)
          ENDDO


          !-- (2) update mesh displacement vector (ale)

          glb%DispBar(:) = glb%DT*glb%VeloBar(:)   + glb%DispBar(:)

          !-- (3) UPDATE MESH POSITION



          DO j = 1, glb%NumNP
             j1 = j*3
             glb%MeshCoor(1,j) = glb%Meshcoor(1,j) + glb%DT*glb%VeloBar(j1 - 2)
             glb%MeshCoor(2,j) = glb%Meshcoor(2,j) + glb%DT*glb%VeloBar(j1 - 1)
             glb%MeshCoor(3,j) = glb%Meshcoor(3,j) + glb%DT*glb%VeloBar(j1    )

          ENDDO



          CALL UpdateMassMatrix(glb)

          !-- (5) CALCULATE R_BAR

          CALL UpdateRbar(glb,Rnet)

          !-- (6) CALCULATE MESH VELOCITY VECTOR
          DO j = 1, glb%NumNP
             jaux = 3*j
             DO k = jaux-2, jaux
                glb%VeloBarOld(k) = glb%VeloBar(k)
                glb%VeloBar(k) = glb%xmass(j)*glb%rho(1)*(-Rnet(k))*glb%kappa ! move - to subroutine
             ENDDO
          ENDDO

          !---(7) CALCULATE R_CO

          !     ! fix for parallel to do boundary then inside
          !      if(numco .gt. 0) then     ! change idpressload bounds to 1:2
          !         call CST_COH(glb%MeshCoor,lmc,matc,Rnet,d,deltan,deltat,
          !     $        sigmax,taumax,glb%Sinit,Sthresh_cst,istep,glb%NumNP,numco,
          !     $        numat_coh,delta,numploadelem,idpressload, pressload,
          !     $        glb%NumNdsBCmesh, idmesh, rmesh, nboundtype, num_fail_coh, 
          !     $        index_fail_coh, npress, nmesh, ielem_coh, iface_coh)
          !      endif

          !-- (8) Enforce Fixed BC's of Mesh for other boundaries not in contact with fluids

!!$         DO j = 1, glb%NumNdsBCmm
!!$            k3 = glb%BCFlagmm(1,j)*3
!!$            k2 = k3-1 
!!$            k1 = k3-2
!!$            IF(glb%BCFlagmm(2,j).eq.0) glb%VeloBar(k1) = 0.d0 ! rmesh(1,j)*meshvelo
!!$            IF(glb%BCFlagmm(3,j).eq.0) glb%VeloBar(k2) = 0.d0 ! rmesh(2,j)*meshvelo
!!$            IF(glb%BCFlagmm(4,j).eq.0) glb%VeloBar(k3) = 0.d0 ! rmesh(3,j)*meshvelo
!!$         ENDDO

          ! -- Enforce Mesh Velocity BC's of Mesh Motion (Not on F/S interface)

          DO j = 1, glb%InterfaceSNumNodes
             k3 = 3*glb%MapNodeS(j)
             k2 = k3 - 1
             k1 = k3 - 2
             glb%VeloBar(k1) = glb%InterfaceSVbar(1,j)
             glb%VeloBar(k2) = glb%InterfaceSVbar(2,j)
             glb%VeloBar(k3) = glb%InterfaceSVbar(3,j)
             IF(glb%InterfaceSVbar(1,j).NE.0.d0) glb%VeloBar(k1) = glb%InterfaceSVbar(1,j)
             IF(glb%InterfaceSVbar(2,j).NE.0.d0) glb%VeloBar(k2) = glb%InterfaceSVbar(2,j)
             IF(glb%InterfaceSVbar(3,j).NE.0.d0) glb%VeloBar(k3) = glb%InterfaceSVbar(3,j) 
          ENDDO



          ! -- Enforce Mesh Velocity BC's of Mesh Motion Given by the fluids.

          DO j = 1, glb%InterfaceSFNumNodes
             k3 = 3*glb%MapNodeSF(j)
             k2 = k3 - 1
             k1 = k3 - 2
             glb%VeloBar(k1) = glb%InterfaceSFVbar(1,j)
             glb%VeloBar(k2) = glb%InterfaceSFVbar(2,j)
             glb%VeloBar(k3) = glb%InterfaceSFVbar(3,j) 
             !             IF(glb%InterfaceSFVbar(1,j).NE.0.d0) glb%VeloBar(k1) = glb%InterfaceSFVbar(1,j)
             !             IF(glb%InterfaceSFVbar(2,j).NE.0.d0) glb%VeloBar(k2) = glb%InterfaceSFVbar(2,j)
             !             IF(glb%InterfaceSFVbar(3,j).NE.0.d0) glb%VeloBar(k3) = glb%InterfaceSFVbar(3,j) 
          ENDDO

          ! -- Enforce Mesh Velocity BC's of Mesh Motion (Not on F/S interface)

!!$          DO j = 1, glb%InterfaceSNumNodes
!!$             k3 = 3*glb%MapNodeS(j)
!!$             k2 = k3 - 1
!!$             k1 = k3 - 2
!!$             glb%VeloBar(k1) = glb%InterfaceSVbar(1,j)
!!$             glb%VeloBar(k2) = glb%InterfaceSVbar(2,j)
!!$             glb%VeloBar(k3) = glb%InterfaceSVbar(3,j)
!!$!             IF(glb%InterfaceSVbar(1,j).NE.0.d0) glb%VeloBar(k1) = glb%InterfaceSVbar(1,j)
!!$!             IF(glb%InterfaceSVbar(2,j).NE.0.d0) glb%VeloBar(k2) = glb%InterfaceSVbar(2,j)
!!$!             IF(glb%InterfaceSVbar(3,j).NE.0.d0) glb%VeloBar(k3) = glb%InterfaceSVbar(3,j)  
!!$          ENDDO
          ! Used for mass conservation
!!$          DO j = 1, glb%InterfaceSFNumNodes
!!$             k3 = 3*glb%MapNodeSF(j)
!!$             k2 = k3 - 1
!!$             k1 = k3 - 2
!!$             glb%VeloBar(k1) = 0.d0
!!$             glb%VeloBar(k2) = 0.d0
!!$             glb%VeloBar(k3) = 1.d0
!!$          ENDDO

          DO j = 1, glb%NumNdsBCmm
             k3 = glb%BCFlagmm(1,j)*3
             k2 = k3-1 
             k1 = k3-2
             IF(glb%BCFlagmm(2,j).EQ.0) glb%VeloBar(k1) = 0.d0 ! rmesh(1,j)*meshvelo
             IF(glb%BCFlagmm(3,j).EQ.0) glb%VeloBar(k2) = 0.d0 ! rmesh(2,j)*meshvelo
             IF(glb%BCFlagmm(4,j).EQ.0) glb%VeloBar(k3) = 0.d0 ! rmesh(3,j)*meshvelo
          ENDDO

          Rnet(:) = 0.d0

          ! -- Calculate mesh acceleration

          glb%AccelBar(:) = ( glb%VeloBar(:) - glb%VeloBarOld(:) ) * glb%DTInv
          IF(glb%iElType.EQ.4)THEN

             CALL v3d4_ale(glb%VeloBar,glb%AccelBar,glb%Disp,glb%VeloHalf,Rnet, &
                  glb%E,glb%xnu,glb%rho,glb%NumNP,glb%NumMatVol, &
                  glb%NumElVol,glb%MatIdVol,glb%ElConnVol,glb%MeshCoor, & !Fix not work for more then one material
                  1,glb%NumElPartBndry)

          ELSE IF(glb%iElType.EQ.10)THEN
             CALL V3D10_ALE(glb%VeloBar,glb%AccelBar,glb%Disp,glb%VeloHalf,Rnet, &
                  glb%E,glb%xnu,glb%rho,glb%NumNP,glb%NumMatVol,  &
                  glb%NumElVol,glb%MatIdVol,glb%ElConnVol,glb%MeshCoor, & !Fix not work for more then one material
                  1,glb%NumElPartBndry)

          ENDIF

       ENDIF ! END ALE OPTION

       IF(DEBUG) print*,'finished ale'




       ! -- Transfer the Tractions Due to Fluid Pressure to the 3D Solid


!!$       IF(glb%DefConfig)THEN

!!$          CALL TractLoadDef(Rnet,glb%NumNP, &
!!$            glb%InterfaceSFElemTract, &
!!$            glb%InterfaceSFNumElems, glb%InterfaceSFNumNodes, &
!!$            glb%InterfaceSFElemConn, &
!!$            glb%MapNodeSF,glb%LwrBnd,glb%UppBnd,glb%Meshcoor,glb%Disp)



!!$       ELSE
!!$
       !          CALL TractLoad(Rnet,glb%NumNP, &
       !               glb%InterfaceSFElemTract, &
       !               glb%InterfaceSFNumElems, glb%InterfaceSFNumNodes, &
       !               glb%InterfaceSFElemConn, &
       !               glb%MapNodeSF,glb%LwrBnd,glb%UppBnd,glb%Meshcoor)
!!$
!!$       ENDIF

       IF(DEBUG) PRINT*,'Applying Pressure Loading'

       IF(glb%iElType.EQ.8)THEN

          IF(glb%InterfaceSFnbNumElems.GT.0)THEN

             ! ... For temporary thermal solver development, must change!! COstoich, 02/23/09
             IF ( glb%HeatTransSoln ) THEN

                CALL HeatLoad_Hex(RnetHT,glb%NumNP,glb%InterfaceSFnbHeatFlux, &
                     glb%InterfaceSFnbNumElems, glb%InterfaceSFnbNumNodes, &
                     glb%InterfaceSFnbElemConn, &
                     glb%MapNodeSFnb,glb%LwrBnd,glb%UppBnd,glb%Meshcoor)

             ELSE

                CALL TractLoad_Hex(Rnet,glb%NumNP, &  ! Don't imposed pressure on deformed coordinates
                     glb%InterfaceSFnbElemTract, &
                     glb%InterfaceSFnbNumElems, glb%InterfaceSFnbNumNodes, &
                     glb%InterfaceSFnbElemConn, &
                     glb%MapNodeSFnb,glb%LwrBnd,glb%UppBnd,glb%Meshcoor)

             ENDIF

          ENDIF
       ELSE

          CALL FluidPressLoad(glb%NumNP,Rnet, &
               glb%InterfaceSFNumElems, glb%InterfaceSFNumNodes, &
               glb%InterfaceSFElemConn, &
               glb%MapNodeSF,glb%LwrBnd,glb%UppBnd,glb%Meshcoor,glb%Disp,glb%MapSFElVolEl,&
               glb%ElConnVol,glb%iElType,glb%NumElVol,glb%InterfaceSFElemTract)

          CALL FluidPressLoad(glb%NumNP,Rnet, &
               glb%InterfaceSFnbNumElems, glb%InterfaceSFnbNumNodes, &
               glb%InterfaceSFnbElemConn, &
               glb%MapNodeSFnb,glb%LwrBnd,glb%UppBnd,glb%Meshcoor,glb%Disp,glb%MapSFnbElVolEl,&
               glb%ElConnVol,glb%iElType,glb%NumElVol,glb%InterfaceSFnbElemTract)
       ENDIF

       !!        ! ... This was in the Rocstar code, I replaced it with what is above COstoich, 5-14-09
       !!        IF(glb%iElType.EQ.8)THEN

       !!           CALL TractLoad_Hex(Rnet,glb%NumNP, &  ! Don't imposed pressure on deformed coordinates
       !!                glb%InterfaceSFElemTract, &
       !!                glb%InterfaceSFNumElems, glb%InterfaceSFNumNodes, &
       !!                glb%InterfaceSFElemConn, &
       !!                glb%MapNodeSF,glb%LwrBnd,glb%UppBnd,glb%Meshcoor)

       !!        ELSE

       !!           CALL FluidPressLoad(glb%NumNP,Rnet, &
       !!                glb%InterfaceSFNumElems, glb%InterfaceSFNumNodes, &
       !!                glb%InterfaceSFElemConn, &
       !!                glb%MapNodeSF,glb%LwrBnd,glb%UppBnd,glb%Meshcoor,glb%Disp,glb%MapSFElVolEl,&
       !!                glb%ElConnVol,glb%iElType,glb%NumElVol,glb%InterfaceSFElemTract)

       !!           CALL FluidPressLoad(glb%NumNP,Rnet, &
       !!                glb%InterfaceSFnbNumElems, glb%InterfaceSFnbNumNodes, &
       !!                glb%InterfaceSFnbElemConn, &
       !!                glb%MapNodeSFnb,glb%LwrBnd,glb%UppBnd,glb%Meshcoor,glb%Disp,glb%MapSFnbElVolEl,&
       !!                glb%ElConnVol,glb%iElType,glb%NumElVol,glb%InterfaceSFnbElemTract)
       !!        ENDIF

       IF(DEBUG) PRINT*,'End Pressure Loading'

!!$       CALL TractLoad(Rnet,glb%NumNP, &
!!$            glb%InterfaceSFElemTract, &
!!$            glb%InterfaceSFNumElems, glb%InterfaceSFNumNodes, &
!!$            glb%InterfaceSFElemConn, &
!!$            glb%MapNodeSF,glb%LwrBnd,glb%UppBnd,glb%Meshcoor,glb%Disp,glb%MapSFElVolEl,&
!!$            glb%ElConnVol,glb%iElType,glb%NumElVol)

       ! -- Enforce tractions (if any) where there are no fluid's domains.
       !    Uses the non-solid/fluid surface mesh.


       IF(DEBUG) PRINT*,'Start traction loading'

       IF(glb%EnforceTractionS.OR.glb%EnforceTractionSF)THEN    
          !         IF(glb%UnDefConfig)THEN
          ! if using cauchy stress
          !             CALL TractPressLoadDef(Rnet,glb%NumNP, &
          !                  glb%InterfaceSNumElems, glb%InterfaceSNumNodes, &
          !                  glb%InterfaceSElemConn, &
          !                  glb%MapNodeS,glb%LwrBnd,glb%UppBnd,glb%Meshcoor,glb%Disp,glb%MapSElVolEl,&
          !                  glb%ElConnVol,glb%iElType,glb%NumElVol,glb%DummyTractVal*glb%prop)
          !             
          !             IF(myid.EQ.0) PRINT*,'Pressure =', glb%DummyTractVal*glb%prop

          !          ELSE

          IF(glb%iElType.EQ.8)THEN

             ! fix need other option
             IF(glb%EnforceTractionS)THEN
                CALL TractLoadPress_Hex(Rnet,glb%NumNP, &
                     glb%InterfaceSNumElems, glb%InterfaceSNumNodes, &
                     glb%InterfaceSElemConn, &
                     glb%MapNodeS,glb%LwrBnd,glb%UppBnd,glb%Meshcoor,glb%DummyTractVal*glb%prop)
             ENDIF

          ELSE

             IF(glb%EnforceTractionS)THEN


                CALL TractPressLoad(Rnet,glb%NumNP, &
                     glb%InterfaceSNumElems, glb%InterfaceSNumNodes, &
                     glb%InterfaceSElemConn, &
                     glb%MapNodeS,glb%LwrBnd,glb%UppBnd,glb%Meshcoor,glb%Disp,glb%MapSElVolEl,&
                     glb%ElConnVol,glb%iElType,glb%NumElVol,glb%DummyTractVal*glb%prop)

                IF(myid.EQ.0) PRINT*,'Pressure Solid =', glb%DummyTractVal*glb%prop

             ENDIF

             IF(glb%EnforceTractionSF)THEN

                CALL TractPressLoad(Rnet,glb%NumNP, &
                     glb%InterfaceSFNumElems, glb%InterfaceSFNumNodes, &
                     glb%InterfaceSFElemConn, &
                     glb%MapNodeSF,glb%LwrBnd,glb%UppBnd,glb%Meshcoor,glb%Disp,glb%MapSFElVolEl,&
                     glb%ElConnVol,glb%iElType,glb%NumElVol,glb%DummyTractVal*glb%prop)

                IF(myid.EQ.0) PRINT*,'Pressure Solid/Fluid =', glb%DummyTractVal*glb%prop


                CALL TractPressLoad(Rnet,glb%NumNP, &
                     glb%InterfaceSFnbNumElems, glb%InterfaceSFnbNumNodes, &
                     glb%InterfaceSFnbElemConn, &
                     glb%MapNodeSFnb,glb%LwrBnd,glb%UppBnd,glb%Meshcoor,glb%Disp,glb%MapSFnbElVolEl,&
                     glb%ElConnVol,glb%iElType,glb%NumElVol,glb%DummyTractVal*glb%prop)

             ENDIF

          ENDIF
          !          ENDIF
       ENDIF


       IF ( .NOT.( glb%IMP ) ) THEN ! Explicit Solver


          !       Rnet(:) = 0.d0

          ! -- Calculate R_in, R_damp

          IF(DEBUG) print*,'start UpdateStructural'

          IF(glb%HeatTransSoln)THEN

             CALL UpdateStructuralHT(glb,NumProcs,Rnet,RnetHT)

             IF(DEBUG) PRINT*,'MaxTemperature',MAXVAL(glb%Temperature)

             ! should not you have updatestrural here
          ELSE
             CALL UpdateStructural(glb,NumProcs,Rnet)
          ENDIF


          IF(DEBUG) print*,'finished UpdateStructural'

          CALL principal_stress(glb%S11,glb%S22,glb%S33, &
               glb%S12,glb%S23,glb%S13, &
               glb%iStrGss,glb%NumElVol,glb%SVonMises)


          ! -- (12) UPDATE THE ACCELERATION AND VELOCITY


          IF(glb%ALEenabled)THEN
             DO j = 1, glb%NumNP
                jaux = 3*j
                DO k = jaux-2, jaux

                   !
                   ! Accelerations are found by
                   !       ..           -1        ext        in
                   !     { D }     = [M]   * ( { F    } - { F   } )
                   !          i+1                           
                   !
                   a1 = Rnet(k) * glb%xmass(j)

                   ! Store Accelerations

                   glb%Accel(k) = a1
                ENDDO
             ENDDO

             !-- (2) update mesh displacement vector (ale)

             !       glb%DispBar(:) = glb%DT*glb%VeloBar(:)   + glb%DispBar(:)

          ELSE IF(.NOT.(glb%DampEnabled))THEN


             DO j = 1, glb%NumNP
                jaux = 3*j
                DO k = jaux-2, jaux
                   !
                   ! Accelerations are found by
                   !       ..           -1        ext        in
                   !     { D }     = [M]   * ( { F    } - { F   } )
                   !          i+1                           
                   !
                   a1 = Rnet(k) * glb%xmass(j)


                   !
                   ! Velocities are found by
                   !       .          .                 ..        ..
                   !     { D }    = { D }  + 0.5*Dt*( { D }  +  { D }    )
                   !          i+1        i                 i         i+1

                   glb%VeloHalf(k) = glb%VeloHalf(k) + glb%DT * ( glb%Accel(k)+a1 ) * 0.5d0

                   ! Store Accelerations

                   glb%Accel(k) = a1
                ENDDO
             ENDDO

          ELSE

             ! Damping Enabled

             DO j = 1, glb%NumNP
                jaux = 3*j
                DO k = jaux-2, jaux
                   DispPrev = glb%Disp(k)

                   ! Displacement is found by
                   !
                   !                         -1   ext    int    damp                 .
                   ! { D }    = Dt**2 * [ M ]  ( F    - F    - F    ) + { D } + Dt*{ D } 
                   !      i+1                     i      i                 i            i-1/2

                   glb%Disp(k)= glb%DT*glb%DT*Rnet(k)*glb%xmass(j) + glb%Disp(k) + glb%DT*glb%VeloHalf(k)

                   ! Update the velocity [ at time (i+1/2)Dt ]
                   !
                   !     .            -1     
                   !   { D }      = Dt   ( { D }      - { D }  )
                   !        i+1/2               n + 1        n


                   glb%VeloHalf(k) = ( glb%Disp(k)-DispPrev )* glb%DTInv

                ENDDO
             ENDDO
          ENDIF

          !     
          ! -- (13) APPLY BOUNDARY CONDITIONS 
          !

          IF(DEBUG) print*,'start bc'

          IF(glb%NumNdsBC.NE.0)THEN
             CALL bc_enforce(glb%NumNdsBC,glb%NumNP,glb%BCFlag,glb%BCvalue,glb%slope,glb%prop,&
                  glb%VeloBndry,glb%AccelBndry,glb%VeloHalf,glb%Accel,glb%Disp,glb%DT,Rnet,&
                  glb%xmass,glb%DampEnabled,glb%CurrTime)
          ENDIF


          IF(DEBUG) print*,'finished bc'


          !-- (1) UPDATE VELOCITY AND DISPLACEMENT VECTORS (structural)

          IF(glb%ALEenabled)THEN
             glb%VeloHalf(:) = glb%DT*glb%Accel(:)     + glb%VeloHalf(:)
             glb%Disp(:)     = glb%DT*glb%VeloHalf(:) + glb%Disp(:)
          ELSE IF(.NOT.(glb%DampEnabled))THEN
             !          glb%VeloHalf(k) = glb%VeloHalf(k) + glb%DT * ( glb%Accel(k)+a1 ) * 0.5d0

             ! Displacement is found by
             !
             !                         -1   ext    int    damp                 .
             ! { D }    = Dt**2 * [ M ]  ( F    - F    - F    ) + { D } + Dt*{ D } 
             !      i+1                     i      i                 i            i-1/2


             !          PRINT*,glb%DT*glb%DT*glb%Accel(:)*0.5d0,glb%DT*glb%VeloHalf(:)

             glb%Disp(:) = glb%DT*glb%DT*glb%Accel(:)*0.5d0 + glb%DT*glb%VeloHalf(:) + glb%Disp(:)
          ENDIF

       ELSE ! Implicit Solver

          ! ... Temporary, for thermal solver development, CHANGE!! COstoich 02/23/09
          IF ( glb%HeatTransSoln ) THEN

             CALL thermal_soln(CurrentTimeStep,CurrentTime,RnetHT,glb,istep)
          ELSE

             CALL implicit_soln(CurrentTimeStep,CurrentTime,Rnet,glb)
          END IF

          nsubstep = 0  ! stop trying to sub-step

          ! DEBUG
          tempf1 = 0.0
          DO i = 1, glb%NumNP
             tempf2 = glb%Disp(3*i-2)*glb%Disp(3*i-2) + glb%Disp(3*i-1)*glb%Disp(3*i-1) + glb%Disp(3*i)*glb%Disp(3*i)
             tempf1 = MAX(tempf1,tempf2)
          ENDDO
          tempf1 = SQRT(tempf1)
          print*,myid,CurrentTime+CurrentTimeStep,'999 999',tempf1
          IF ( ABS(MAXVAL(glb%Disp(:))) > ABS(MINVAL(glb%Disp(:))) ) THEN
             print*,myid,CurrentTime+CurrentTimeStep,'888 888',MAXVAL(glb%Disp(:))
          ELSE
             print*,myid,CurrentTime+CurrentTimeStep,'888 888',MINVAL(glb%Disp(:))
          END IF
          ! END DEBUG

       END IF ! IMP

       IF(DEBUG) PRINT*,'MAX DISPLACEMENT =', MAXVAL(glb%Disp(:))

!!$       DO i = 1, 3*glb%numnp
!!$          IF(glb%Disp(i).NE.0.d0) PRINT*,i,glb%Disp(i)
!!$       ENDDO

       glb%DispTotal(:) = glb%Disp(:) + glb%DispBar(:)

!!$       DO i = 1, glb%NumNodeIO
!!$          WRITE(33+i,*) glb%CurrTime, glb%DispTotal(3*glb%NodeIO(i)-2:3*glb%NodeIO(i))
!!$       ENDDO

       !       IF(glb%NumNodeIO.NE.0)THEN
       !          write(400,*) glb%CurrTime, glb%Disp(glb%NumNodeIO*3-2)
       ! print*,myid,glb%Meshcoor(1:3,glb%NumNodeIO)
       !       endif

       !       PRINT*,'Max Displacement', MAXval(glb%Disp)
       !        PRINT*,'Min Displacement', minval(glb%Disp)

       ! -- Mass/Volume Conservation


       IF(DEBUG) print*,'finished mass volume conservation'

       glb%TotalMassSolidp = 0.d0
       glb%TotalGeomVolp = 0.d0

       IF(glb%iElType.EQ.4 .OR. glb%iElType.EQ.10 .AND.(.NOT.(glb%NdBasedEl)) )THEN
          CALL V3D4_volume(glb%MeshCoor,glb%ElConnVol,glb%MatIdVol,glb%rho, &
               glb%NumNP,glb%NumElVol,glb%NumMatVol,glb%Disp,1,glb%NumElVol,&
               glb%TotalMassSolidp,glb%TotalGeomVolp,glb%TotalGeomUndefVolp, &
               glb%iElType)
       ENDIF

       IF (glb%HeatTransSoln)THEN
          DO i = 1, glb%NumProbesNd
             IF(glb%PointOnProc(i))THEN
                WRITE(ichr1,'(i4.4)') i
                WRITE(ichr2,'(I4.4)') myid

                OPEN(440+i,FILE='Rocfrac/Rocout/Probe.'//ichr1//'.'//ichr2,POSITION='APPEND')
                WRITE(440+i,*) glb%CurrTime, glb%Temperature(glb%ProbeNd(i))
                close(440+i)
             ENDIF
          ENDDO
       ELSE
          DO i = 1, glb%NumProbesNd
             IF(glb%PointOnProc(i))THEN
                WRITE(ichr1,'(i4.4)') i
                WRITE(ichr2,'(I4.4)') myid

                OPEN(440+i,FILE='Rocfrac/Rocout/Probe.'//ichr1//'.'//ichr2,POSITION='APPEND')
                WRITE(440+i,*) glb%CurrTime, glb%Disp( glb%ProbeNd(i)*3-2), glb%Disp( glb%ProbeNd(i)*3-1) , glb%Disp( glb%ProbeNd(i)*3)  ! glb%Temperature(glb%ProbeNd(i))
                close(440+i)
             ENDIF
          ENDDO
       ENDIF


    ENDDO  ! substepping
 
    DEALLOCATE(Rnet)
 

    IF(DEBUG) print*,'structural boundary conditions start'

    icnt1 = 1
    icnt5 = 0
    DO i = 1, glb%NumNdsBCcrypt

       NdBCflag = MOD(glb%BCFlagCrypt(2,i),100)
       ! structural boundary conditions
       IF(NdBCflag.GT.0) THEN

          icnt5 = icnt5 + 1

          glb%BCValueGlb(icnt1:icnt1+2) = glb%AccelBndry(icnt5*3-2:icnt5*3)
          icnt1 = icnt1 + 3
          glb%BCValueGlb(icnt1:icnt1+2) = glb%VeloBndry(icnt5*3-2:icnt5*3)
          icnt1 = icnt1 + 3

       ENDIF
    ENDDO

    IF(myid.EQ.0) PRINT*,'RocFrac :: END SOLID STEP'

    CALL RocFracInterfaceBuff(glb)

    RETURN
  END SUBROUTINE RocFracSoln

!
! ----------------------------------------------------------------------- RocFracInterfaceBuff

  SUBROUTINE RocFracInterfaceBuff(glb)

!!****f* Rocfrac/Rocfrac/Source/RocFracMain
!!
!!  NAME
!!    RocFracInterfaceBuff
!!
!!  FUNCTION
!!  Passes the variables to the fluid code. It transfers the new mesh 
!!  velocity and displacement to the interface mesh  arrays that are
!!  registered with RocCom 
!!
!!
!!***
    
    IMPLICIT NONE

    TYPE(ROCFRAC_GLOBAL) :: glb

    INTEGER :: iInterfaceNode
    INTEGER :: k,k1,k2,k3
    
    INTEGER :: myid,ierr

    DO iInterfaceNode = 1, glb%InterfaceSFNumNodes
       k = glb%MapNodeSF(iInterfaceNode)
       k3 = 3*k
       k2 = k3 - 1
       k1 = k3 - 2

! Move the interface mesh to match the volume mesh
       glb%InterfaceSFNodalCoors(1:3,iInterfaceNode) = glb%MeshCoor(1:3,k)

! Total Displacement due to structural tractions (not due to mesh motion)

       glb%InterfaceSFTotalNodalDisps(1,iInterfaceNode) = glb%Disp(k1)
       glb%InterfaceSFTotalNodalDisps(2,iInterfaceNode) = glb%Disp(k2)
       glb%InterfaceSFTotalNodalDisps(3,iInterfaceNode) = glb%Disp(k3)

! 
! Incremental Total Displacement
!       
       glb%InterfaceSFNodalDisps(1,iInterfaceNode) = glb%Disp(k1) + glb%DispBar(k1) 
       glb%InterfaceSFNodalDisps(2,iInterfaceNode) = glb%Disp(k2) + glb%DispBar(k2) 
       glb%InterfaceSFNodalDisps(3,iInterfaceNode) = glb%Disp(k3) + glb%DispBar(k3)

! 
! Acceleration
!       
       glb%InterfaceSFNodalAccel(1,iInterfaceNode) = glb%Accel(k1) 
       glb%InterfaceSFNodalAccel(2,iInterfaceNode) = glb%Accel(k2)
       glb%InterfaceSFNodalAccel(3,iInterfaceNode) = glb%Accel(k3) 

! Passing the velocity, Vs
! Structural Velocity ! same
            
       glb%InterfaceSFNodalVels(1,iInterfaceNode) = glb%VeloHalf(k1) ! fix: need more terms
       glb%InterfaceSFNodalVels(2,iInterfaceNode) = glb%VeloHalf(k2) ! fix: need more terms
       glb%InterfaceSFNodalVels(3,iInterfaceNode) = glb%VeloHalf(k3) ! fix: need more terms

    ENDDO

    DO iInterfaceNode = 1, glb%InterfaceSFnbNumNodes
       k = glb%MapNodeSFnb(iInterfaceNode)
       k3 = 3*k
       k2 = k3 - 1
       k1 = k3 - 2

! Move the interface mesh to match the volume mesh
       glb%InterfaceSFnbNodalCoors(1:3,iInterfaceNode) = glb%MeshCoor(1:3,k)

! Total Displacement due to structural tractions (not due to mesh motion)

       glb%InterfaceSFnbTotalNodalDisps(1,iInterfaceNode) = glb%Disp(k1)
       glb%InterfaceSFnbTotalNodalDisps(2,iInterfaceNode) = glb%Disp(k2)
       glb%InterfaceSFnbTotalNodalDisps(3,iInterfaceNode) = glb%Disp(k3)

! 
! Incremental Total Displacement
!       
       glb%InterfaceSFnbNodalDisps(1,iInterfaceNode) = glb%Disp(k1) + glb%DispBar(k1) 
       glb%InterfaceSFnbNodalDisps(2,iInterfaceNode) = glb%Disp(k2) + glb%DispBar(k2) 
       glb%InterfaceSFnbNodalDisps(3,iInterfaceNode) = glb%Disp(k3) + glb%DispBar(k3)

! Passing the velocity, Vs
! Structural Velocity ! same
            
       glb%InterfaceSFnbNodalVels(1,iInterfaceNode) = glb%VeloHalf(k1) ! fix: need more terms
       glb%InterfaceSFnbNodalVels(2,iInterfaceNode) = glb%VeloHalf(k2) ! fix: need more terms
       glb%InterfaceSFnbNodalVels(3,iInterfaceNode) = glb%VeloHalf(k3) ! fix: need more terms

! Acceleration
            
       glb%InterfaceSFnbNodalAccel(1,iInterfaceNode) = glb%Accel(k1)
       glb%InterfaceSFnbNodalAccel(2,iInterfaceNode) = glb%Accel(k2)
       glb%InterfaceSFnbNodalAccel(3,iInterfaceNode) = glb%Accel(k3)


    ENDDO

    IF(glb%HeatTransSoln)THEN
       ! ... if implicit thermal solver (Roctherm), then only non-burning
       ! ... interaction is supported
       IF (glb%IMP .EQV. .TRUE.) THEN
          DO iInterfaceNode = 1, glb%InterfaceSFnbNumNodes
             k = glb%MapNodeSFnb(iInterfaceNode)
             
             glb%InterfaceSFnbNodalTemp(iInterfaceNode) = glb%Temperature(k)
          ENDDO
       ELSE
          DO iInterfaceNode = 1, glb%InterfaceSFNumNodes
             k = glb%MapNodeSF(iInterfaceNode)
             
             glb%InterfaceSFNodalTemp(iInterfaceNode) = glb%Temperature(k)
             
          ENDDO
       ENDIF
    ENDIF


    
! Is this needed still?
!!$    DO iInterfaceNode = 1, glb%InterfaceSNumNodes
!!$
!!$       k = glb%MapNodeS(iInterfaceNode)
!!$
!!$! Move the interface mesh to match the volume mesh
!!$       glb%InterfaceSNodalCoors(1:3,iInterfaceNode) = glb%MeshCoor(1:3,k)
!!$    ENDDO

    RETURN
         
  END SUBROUTINE RocFracInterfaceBuff

  SUBROUTINE RocFracUpdateInbuff( glb, alpha)
    TYPE(ROCFRAC_GLOBAL), POINTER :: glb
    REAL*8, INTENT(IN)            :: alpha
    INTEGER :: TriConn(1:3)
    INTEGER :: i
    REAL*8 :: radius, signx,signy,signz
    REAL*8 :: prop 
    INTEGER, SAVE :: icnt
    REAL*8 :: TractVal
    REAL*8 :: slope,mag

!!$    icnt = icnt + 1
!!$
!!$    IF(mod(icnt,26).eq.0) THEN
!!$       prop = float(icnt/26)*glb%DummyTractVal
!!$       print*,'**PRESSURE=', prop
!!$       write(44,*) glb%CurrTime,prop
!!$    endif

! Wind Block 

!!$    DO i = 1, glb%InterfaceSFNumElems
!!$     
!!$       TriConn(1:3) = glb%InterfaceSFElemConn(glb%LwrBnd:glb%UppBnd,i)
!!$
!!$       IF(glb%coor(1,glb%MapNodeSF(TriConn(1))) .LT. 0.276d0 .AND. &
!!$            glb%coor(1,glb%MapNodeSF(TriConn(2))) .LT. 0.276d0 .AND. &
!!$            glb%coor(1,glb%MapNodeSF(TriConn(3))) .LT. 0.276d0 ) THEN
!!$
!!$          glb%InterfaceSFElemTract(1,i) = glb%DummyTractVal
!!$          glb%InterfaceSFElemTract(2:3,i) = 0.d0
!!$       ELSE
!!$          glb%InterfaceSFElemTract(:,i) = 0.d0
!!$       ENDIF
!!$       
!!$    ENDDO
!!$
!!$! -- Enforce Mesh Velocity BC's of Mesh Motion (Not on F/S interface)
!!$
!!$!    glb%InterfaceSVbar(:,1:glb%InterfaceSNumNodes) = 0.d0
!!$
!!$! -- Enforce Mesh Velocity BC's of Mesh Motion Given by the fluids.
!!$
!!$    DO i = 1, glb%InterfaceSFNumNodes
!!$
!!$       IF(glb%coor(1,glb%MapNodeSF(i)) .LT. 0.276d0) THEN
!!$          glb%InterfaceSFVbar(1,i) = glb%DummyBurnRate
!!$          glb%InterfaceSFVbar(2:3,i) = 0.d0
!!$       ELSE
!!$          glb%InterfaceSFVbar(:,i) = 0.d0
!!$       ENDIF
!!$       
!!$    ENDDO

! Scaleability

!!$     DO i = 1, glb%InterfaceSFNumElems
!!$     
!!$       TriConn(1:3) = glb%InterfaceSFElemConn(glb%LwrBnd:glb%UppBnd,i)
!!$
!!$       signx= SUM( glb%coor(1,glb%MapNodeSF(TriConn(1:3))) ) / 3.0d0
!!$       signy= SUM( glb%coor(2,glb%MapNodeSF(TriConn(1:3))) ) / 3.0d0
!!$
!!$       radius = SQRT( signx**2 + signy**2)
!!$
!!$       glb%InterfaceSFElemTract(1,i) = signx/radius*glb%DummyTractVal
!!$       glb%InterfaceSFElemTract(2,i) = signy/radius*glb%DummyTractVal
!!$       glb%InterfaceSFElemTract(3,i) = 0.d0
!!$       
!!$    ENDDO
!!$
!!$! -- Enforce Mesh Velocity BC's of Mesh Motion (Not on F/S interface)
!!$
!!$!    glb%InterfaceSVbar(:,1:glb%InterfaceSNumNodes) = 0.d0
!!$
!!$! -- Enforce Mesh Velocity BC's of Mesh Motion Given by the fluids.
!!$
!!$    DO i = 1, glb%InterfaceSFNumNodes
!!$
!!$
!!$       radius = SQRT( glb%coor(1,glb%MapNodeSF(i))**2 + &
!!$            glb%coor(2,glb%MapNodeSF(i))**2)
!!$ 
!!$       glb%InterfaceSFVbar(1,i) = glb%coor(1,glb%MapNodeSF(i))/radius*glb%DummyBurnRate
!!$       glb%InterfaceSFVbar(2,i) = glb%coor(2,glb%MapNodeSF(i))/radius*glb%DummyBurnRate
!!$       glb%InterfaceSFVbar(3,i) = 0.d0
!!$       
!!$    ENDDO      

! Mass Conservation
!!$
!!$     DO i = 1, glb%InterfaceSFNumElems
!!$     
!!$       TriConn(1:3) = glb%InterfaceSFElemConn(glb%LwrBnd:glb%UppBnd,i)
!!$
!!$       glb%InterfaceSFElemTract(1:2,i) = 0.d0
!!$       glb%InterfaceSFElemTract(3,i) = glb%DummyTractVal
!!$       
!!$    ENDDO
!!$
!!$! -- Enforce Mesh Velocity BC's of Mesh Motion (Not on F/S interface)
!!$
!!$!    glb%InterfaceSVbar(:,1:glb%InterfaceSNumNodes) = 0.d0
!!$
!!$! -- Enforce Mesh Velocity BC's of Mesh Motion Given by the fluids.
!!$
!!$    DO i = 1, glb%InterfaceSFNumNodes
!!$
!!$       glb%InterfaceSFVbar(3,i) = glb%DummyBurnRate
!!$       glb%InterfaceSFVbar(1:2,i) = 0.d0
!!$       
!!$    ENDDO 

 
! Hollow Sphere

!!$     DO i = 1, glb%InterfaceSFNumElems
!!$     
!!$       TriConn(1:3) = glb%InterfaceSFElemConn(glb%LwrBnd:glb%UppBnd,i)
!!$
!!$       signx= SUM( glb%coor(1,glb%MapNodeSF(TriConn(1:3))) ) / 3.0d0
!!$       signy= SUM( glb%coor(2,glb%MapNodeSF(TriConn(1:3))) ) / 3.0d0
!!$       signz= SUM( glb%coor(3,glb%MapNodeSF(TriConn(1:3))) ) / 3.0d0
!!$
!!$       radius = SQRT( signx**2 + signy**2 + signz**2)
!!$
!!$       glb%InterfaceSFElemTract(1,i) = signx/radius*glb%DummyTractVal
!!$       glb%InterfaceSFElemTract(2,i) = signy/radius*glb%DummyTractVal
!!$       glb%InterfaceSFElemTract(3,i) = signz/radius*glb%DummyTractVal
!!$       
!!$    ENDDO
!!$
!!$! -- Enforce Mesh Velocity BC's of Mesh Motion (Not on F/S interface)
!!$
!!$!    glb%InterfaceSVbar(:,1:glb%InterfaceSNumNodes) = 0.d0
!!$
!!$! -- Enforce Mesh Velocity BC's of Mesh Motion Given by the fluids.
!!$
!!$    DO i = 1, glb%InterfaceSFNumNodes
!!$
!!$       radius = SQRT( glb%coor(1,glb%MapNodeSF(i))**2 + &
!!$            glb%coor(2,glb%MapNodeSF(i))**2 + glb%coor(3,glb%MapNodeSF(i))**2)
!!$ 
!!$       glb%InterfaceSFVbar(1,i) = glb%coor(1,glb%MapNodeSF(i))/radius*glb%DummyBurnRate
!!$       glb%InterfaceSFVbar(2,i) = glb%coor(2,glb%MapNodeSF(i))/radius*glb%DummyBurnRate
!!$       glb%InterfaceSFVbar(3,i) = glb%coor(3,glb%MapNodeSF(i))/radius*glb%DummyBurnRate
!!$       
!!$    ENDDO  
! -- Enforce Mesh Velocity BC's of Mesh Motion (Not on F/S interface)

!    glb%InterfaceSVbar(:,1:glb%InterfaceSNumNodes) = 0.d0

! -- Enforce Mesh Velocity BC's of Mesh Motion Given by the fluids.

    DO i = 1, glb%InterfaceSFNumNodes

       glb%InterfaceSFVbar(1,i) = 0.00075 ! glb%DummyBurnRate
       glb%InterfaceSFVbar(2,i) = 0.
       glb%InterfaceSFVbar(3,i) = 0.
       
    ENDDO  


! 2 Material Beam


!!$     DO i = 1, glb%InterfaceSFNumElems
!!$     
!!$       TriConn(1:3) = glb%InterfaceSFElemConn(glb%LwrBnd:glb%UppBnd,i)
!!$
!!$       glb%InterfaceSFElemTract(1,i) = 0.d0
!!$       glb%InterfaceSFElemTract(2,i) = -glb%DummyTractVal
!!$       glb%InterfaceSFElemTract(3,i) = 0.d0
!!$       
!!$    ENDDO
!!$

! 5 sided cube

!!$    glb%InterfaceSFElemTract(:,:) = 0.d0
!!$    IF(glb%CurrTime.LE.0.1d0)THEN
!!$       slope = glb%DummyTractVal/0.1d0
!!$       
!!$       prop = slope*glb%CurrTime
!!$       print*,' PRESSURE =', prop
!!$    else
!!$       prop = 0.d0
!!$    endif
!!$
!!$    mag = 0.d0
!!$
!!$    DO i = 1, glb%InterfaceSFNumElems
!!$     
!!$       TriConn(1:3) = glb%InterfaceSFElemConn(1:3,i)
!!$
!!$! y = 0.9 face
!!$
!!$       IF(glb%coor(2,glb%MapNodeSF(TriConn(1))) .GT. 0.89999d0 .AND. &
!!$            glb%coor(2,glb%MapNodeSF(TriConn(2))) .GT. 0.89999d0 .AND. &
!!$            glb%coor(2,glb%MapNodeSF(TriConn(3))) .GT. 0.89999d0 ) THEN
!!$          
!!$          glb%InterfaceSFElemTract(1,i) = 0.d0
!!$          glb%InterfaceSFElemTract(2,i) = -prop
!!$          glb%InterfaceSFElemTract(3,i) = 0.d0
!!$
!!$! z = 0.9 face
!!$       ELSE IF(glb%coor(3,glb%MapNodeSF(TriConn(1))) .GT. 0.89999d0+mag .AND. &
!!$            glb%coor(3,glb%MapNodeSF(TriConn(2))) .GT. 0.89999d0+mag .AND. &
!!$            glb%coor(3,glb%MapNodeSF(TriConn(3))) .GT. 0.89999d0+mag ) THEN
!!$          
!!$          glb%InterfaceSFElemTract(1,i) = 0.d0
!!$          glb%InterfaceSFElemTract(2,i) = 0.d0
!!$          glb%InterfaceSFElemTract(3,i) = -prop
!!$
!!$! z = 0.0 face
!!$       ELSE IF(glb%coor(3,glb%MapNodeSF(TriConn(1))) .LT. 0.00001d0+mag .AND. &
!!$            glb%coor(3,glb%MapNodeSF(TriConn(2))) .LT. 0.00001d0 +mag.AND. &
!!$            glb%coor(3,glb%MapNodeSF(TriConn(3))) .LT. 0.00001d0+mag ) THEN
!!$          
!!$          glb%InterfaceSFElemTract(1,i) = 0.d0
!!$          glb%InterfaceSFElemTract(2,i) = 0.d0
!!$          glb%InterfaceSFElemTract(3,i) = prop
!!$
!!$! x = 0.9 face
!!$       ELSE IF(glb%coor(1,glb%MapNodeSF(TriConn(1))) .GT. 0.89999d0+mag .AND. &
!!$            glb%coor(1,glb%MapNodeSF(TriConn(2))) .GT. 0.89999d0+mag .AND. &
!!$            glb%coor(1,glb%MapNodeSF(TriConn(3))) .GT. 0.89999d0+mag ) THEN
!!$          
!!$          glb%InterfaceSFElemTract(1,i) = -prop
!!$          glb%InterfaceSFElemTract(2,i) = 0.d0
!!$          glb%InterfaceSFElemTract(3,i) = 0.d0
!!$
!!$! x = 0.0 face
!!$       ELSE IF(glb%coor(1,glb%MapNodeSF(TriConn(1))) .LT. 0.00001d0+mag .AND. &
!!$            glb%coor(1,glb%MapNodeSF(TriConn(2))) .LT. 0.00001d0+mag .AND. &
!!$            glb%coor(1,glb%MapNodeSF(TriConn(3))) .LT. 0.00001d0+mag ) THEN
!!$          
!!$          glb%InterfaceSFElemTract(1,i) = prop 
!!$          glb%InterfaceSFElemTract(2,i) = 0.d0
!!$          glb%InterfaceSFElemTract(3,i) = 0.d0
!!$
!!$       endif
!!$    
!!$    ENDDO

! Cantilever beam

!!$    glb%InterfaceSFElemTract(:,:) = 0.d0
!!$
!!$    IF(glb%CurrTime.LE.0.05d0)THEN
!!$
!!$       slope = glb%DummyTractVal/0.05d0
!!$
!!$       prop = slope*glb%CurrTime
!!$       print*,' PRESSURE =', prop
!!$    ELSE
!!$       prop = 0.d0
!!$    ENDIF
!!$
!!$    DO i = 1, glb%InterfaceSFNumElems
!!$     
!!$       TriConn(1:3) = glb%InterfaceSFElemConn(1:3,i)
!!$
!!$! y = 0.25 face
!!$
!!$       IF(glb%coor(2,glb%MapNodeSF(TriConn(1))) .GT.   0.2499d0 .AND. &
!!$            glb%coor(2,glb%MapNodeSF(TriConn(2))) .GT. 0.2499d0 .AND. &
!!$            glb%coor(2,glb%MapNodeSF(TriConn(3))) .GT. 0.2499d0 ) THEN
!!$          
!!$          glb%InterfaceSFElemTract(1,i) = 0.d0
!!$          glb%InterfaceSFElemTract(2,i) = -prop
!!$          glb%InterfaceSFElemTract(3,i) = 0.d0
!!$
!!$       ENDIF    
!!$    ENDDO

! Propellent/Case

!!$     DO i = 1, glb%InterfaceSFNumElems
!!$     
!!$       TriConn(1:3) = glb%InterfaceSFElemConn(glb%LwrBnd:glb%UppBnd,i)
!!$
!!$       signx= SUM( glb%coor(1,glb%MapNodeSF(TriConn(1:3))) ) / 3.0d0
!!$       signy= SUM( glb%coor(2,glb%MapNodeSF(TriConn(1:3))) ) / 3.0d0
!!$
!!$       radius = SQRT( signx**2 + signy**2)
!!$
!!$       TractVal = MIN(1.d0, Prop/25.d0  )*glb%DummyTractVal
!!$
!!$       glb%InterfaceSFElemTract(1,i) = signx/radius*TractVal
!!$       glb%InterfaceSFElemTract(2,i) = signy/radius*TractVal
!!$       glb%InterfaceSFElemTract(3,i) = 0.d0
!!$       
!!$    ENDDO
      
  END SUBROUTINE RocFracUpdateInbuff
   
END MODULE RocFracMain

