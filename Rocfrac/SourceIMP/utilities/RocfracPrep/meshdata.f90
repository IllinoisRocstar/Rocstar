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
MODULE meshdata


  USE Generic_List, ONLY : Link_Ptr_Type,Link_Type,List_Type
  USE Generic_List, ONLY : LI_Init_List,LI_Add_To_Head,LI_Get_Head, &
       LI_Remove_Head,LI_Get_Next,LI_Associated,LI_Get_Len

  IMPLICIT NONE
  
  INTEGER :: io_input       ! unit id of control deck
  PARAMETER(io_input = 10)

!  - DECK FILE PRAMETERS
  CHARACTER*20 :: prefx   ! I/O file prefx
  INTEGER :: prefx_lngth  ! length of prefx
!  - Fastest Dilatitial wave speed
  REAL*8 :: cd_fastest
!  - Derive Type for the boundary conditions
  TYPE bcvalues
     INTEGER :: b1, b2, b3
     REAL*8 :: bc1, bc2, bc3
  END TYPE bcvalues
!  - Stores the boundary conditions data values
  TYPE(bcvalues), DIMENSION(32) :: bc_conditions

! PRIMARY MESH DATA (i.e. mesh created with ansys, patran, or mesh3d):

!   - Number of nodes 
  INTEGER :: numnp_prmry
!   - Number of volumetric elements  
  INTEGER :: numelv_prmry
!   - meshing software flag
  INTEGER :: iansys         ! 0- no, 1-yes
  INTEGER :: ipatran        ! 0- no, 1-yes
  INTEGER :: itetmesh       ! 0- no, 1-yes
  INTEGER :: ipatcohin      ! 0- no, 1-yes
  INTEGER :: itetcohin      ! 0- no, 1-yes
  
  INTEGER :: numvertx
  
  REAL*8 :: ConvertUnit
      

!   - Coordinate array
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: coor
!   - Coordinate array
  REAL*8, ALLOCATABLE, DIMENSION(:) :: press_nodal
!   - Connectivity array for volumetric elements
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lmelv_prmry


!   - Nodal corner boundary flag
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ibcaxi
!   - Partioned element's processor (from METIS)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: epart

!   - Number of surface triangles with applied tractions
!  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: numel_2d
!   - Number of surface nodes with applied mesh motion velocity
  INTEGER, ALLOCATABLE, DIMENSION(:) :: numnp_2d
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NumElPerProc
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NumNdPerProc
! sub
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: lmtri
  INTEGER, DIMENSION(:), ALLOCATABLE :: elm_2D
  INTEGER, DIMENSION(:), ALLOCATABLE :: elm_2D_flag
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: neigh_2d
  INTEGER, DIMENSION(:), ALLOCATABLE :: epart_2d
! sub

!   - Node tracking
  INTEGER :: NumNodeIO
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NodeIO

!   - Node 

  INTEGER :: IOformat     ! 0 = unformatted (default), 1 = formatted 
        ! 0 = no output, 1 = output
  INTEGER :: iopmvis  ! PMVIS software to view the partition 0=no, 1 = yes
  INTEGER :: ipress  ! flag for if we have pressure loading 0=no, 1 = yes

  INTEGER :: IntFaceFlag

!
! -- derived surface type
!
!$$$        TYPE request
!$$$            INTEGER, DIMENSION(1:6) :: Conn ! list of nodes
!$$$            INTEGER :: NumEl
!$$$            INTEGER :: NumNp
!$$$            TYPE(request), POINTER :: next
!$$$        END TYPE request
!$$$      
!$$$        TYPE(request), POINTER :: head, tail
!$$$        TYPE(request), POINTER :: item, first

  INTEGER, ALLOCATABLE, DIMENSION(:) :: iNdsBurnFlg
  
  INTEGER :: NumMat
  
  INTEGER, DIMENSION(:), ALLOCATABLE :: MatId
  INTEGER, DIMENSION(:), ALLOCATABLE :: ElTypeId

      ! limit of 10 nodes to monitor history

  INTEGER, DIMENSION(1:10) :: NdHistory
  INTEGER :: NumNdHistory

  INTEGER :: numbc_prmry
  INTEGER :: numbc_prmry_mm
  INTEGER :: numbc_prmry_ht
!----- dimension:  <1>  old node number   <2> processor
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ik1

! User-defined list element
! The Link_Type field MUST be the FIRST in the user-defined list element 
! Note pointer to data so as to easily create sublists

  TYPE User_Type_BC
     TYPE(Link_Type) :: Link
     TYPE(User_Data_Type_BC), POINTER :: DATA  
  END TYPE User_Type_BC
  
  TYPE User_Data_Type_BC
     INTEGER :: BC_nodeGlb
     INTEGER :: BC_flagGlb
  END TYPE User_Data_Type_BC
  
! Auxilliary data type required for the transfer function
  TYPE User_Ptr_Type_BC
     TYPE(User_Type_BC), POINTER :: P
  END TYPE User_Ptr_Type_BC
  
  TYPE(User_Ptr_Type_BC)  :: User_BC

  TYPE(List_Type) :: BC_structural
  TYPE(List_Type) :: BC_meshmotion
  TYPE(List_Type) :: BC_thermal
  
  TYPE User_Type_MapNd
     TYPE(Link_Type) :: Link
     TYPE(User_Data_Type_MapNd), POINTER :: DATA  
  END TYPE User_Type_MapNd
  
  TYPE User_Data_Type_MapNd
     INTEGER :: LocNd
     INTEGER :: Proc
  END TYPE User_Data_Type_MapNd

! Auxilliary data type required for the transfer function
  TYPE User_Ptr_Type_MapNd
     TYPE(User_Type_MapNd), POINTER :: P
  END TYPE User_Ptr_Type_MapNd
  
  TYPE(User_Ptr_Type_MapNd)  :: User_MapNd
  
  TYPE(List_Type), ALLOCATABLE, DIMENSION(:) :: MapNd_Glb2LocProc
  
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NumBC_structural, NumBC_meshmotion, NumBC_thermal
  
  
  INTEGER, PARAMETER :: MaxNumberOfProcsToShareNode = 8
  
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ProcNdList
  INTEGER, ALLOCATABLE, DIMENSION(:) ::  NumProcPerNd
  
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NodeFlag

! User-defined list element
! The Link_Type field MUST be the FIRST in the user-defined list element 
! Note pointer to data so as to easily create sublists

  TYPE User_Type_SurfMesh_tri3
     TYPE(Link_Type) :: Link
     TYPE(User_Data_Type_SurfMesh_tri3), POINTER :: DATA  
  END TYPE User_Type_SurfMesh_tri3

  TYPE User_Data_Type_SurfMesh_tri3
     INTEGER, DIMENSION(1:4) :: ElemData
  END TYPE User_Data_Type_SurfMesh_tri3
  
! Auxilliary data type required for the transfer function
  TYPE User_Ptr_Type_SurfMesh_tri3
     TYPE(User_Type_SurfMesh_tri3), POINTER :: P
  END TYPE User_Ptr_Type_SurfMesh_tri3
  
  TYPE(User_Ptr_Type_SurfMesh_tri3)  :: User_SurfMesh_tri3
  
  TYPE(List_Type) :: SurfMesh_tri3_S
  TYPE(List_Type) :: SurfMesh_tri3_SF
  
! User-defined list element
! The Link_Type field MUST be the FIRST in the user-defined list element 
! Note pointer to data so as to easily create sublists

  TYPE User_Type_SurfMesh_tri6
     TYPE(Link_Type) :: Link
     TYPE(User_Data_Type_SurfMesh_tri6), POINTER :: DATA  
  END TYPE User_Type_SurfMesh_tri6
  
  TYPE User_Data_Type_SurfMesh_tri6
     INTEGER, DIMENSION(1:7) :: ElemData
  END TYPE User_Data_Type_SurfMesh_tri6

! Auxilliary data type required for the transfer function
  TYPE User_Ptr_Type_SurfMesh_tri6
     TYPE(User_Type_SurfMesh_tri6), POINTER :: P
  END TYPE User_Ptr_Type_SurfMesh_tri6
  
  TYPE(User_Ptr_Type_SurfMesh_tri6)  :: User_SurfMesh_tri6
  
  TYPE(List_Type) :: SurfMesh_tri6_S
  TYPE(List_Type) :: SurfMesh_tri6_SF
  
! User-defined list element
! The Link_Type field MUST be the FIRST in the user-defined list element 
! Note pointer to data so as to easily create sublists

  TYPE User_Type_SurfMesh_hex8
     TYPE(Link_Type) :: Link
     TYPE(User_Data_Type_SurfMesh_hex8), POINTER :: DATA  
  END TYPE User_Type_SurfMesh_hex8
  
  TYPE User_Data_Type_SurfMesh_hex8
     INTEGER, DIMENSION(1:5) :: ElemData
  END TYPE User_Data_Type_SurfMesh_hex8

! Auxilliary data type required for the transfer function
  TYPE User_Ptr_Type_SurfMesh_hex8
     TYPE(User_Type_SurfMesh_hex8), POINTER :: P
  END TYPE User_Ptr_Type_SurfMesh_hex8
  
  TYPE(User_Ptr_Type_SurfMesh_hex8)  :: User_SurfMesh_hex8

  TYPE(List_Type) :: SurfMesh_hex8_S
  TYPE(List_Type) :: SurfMesh_hex8_SF

! User-defined list element
! The Link_Type field MUST be the FIRST in the user-defined list element 
! Note pointer to data so as to easily create sublists

  TYPE User_Type_ProcNodeList
     TYPE(Link_Type) :: Link
     TYPE(User_Data_Type_ProcNodeList), POINTER :: DATA  
  END TYPE User_Type_ProcNodeList

  TYPE User_Data_Type_ProcNodeList
     INTEGER :: LocNdNum
  END TYPE User_Data_Type_ProcNodeList

! Auxilliary data type required for the transfer function
  TYPE User_Ptr_Type_ProcNodeList
     TYPE(User_Type_ProcNodeList), POINTER :: P
  END TYPE User_Ptr_Type_ProcNodeList
  
  TYPE(User_Ptr_Type_ProcNodeList)  :: User_ProcNodeList
  
  TYPE(List_Type), POINTER, DIMENSION(:) :: ProcNodeList

! type of surface mesh
      ! = 3 for all tri
      ! = 4 for all quad
      ! = 5 for mixed

  INTEGER :: MeshType2D
  
  INTEGER, DIMENSION(:), POINTER :: ElFlag
  
  INTEGER, DIMENSION(:,:),POINTER :: NumElHex2D
  INTEGER, DIMENSION(:,:),POINTER :: NumEltet2D

  INTEGER, DIMENSION(:,:), POINTER :: ElConnTable
  INTEGER, DIMENSION(:), POINTER :: matType

  
  INTEGER, DIMENSION(:), POINTER :: NumElVolMat
  INTEGER, DIMENSION(:), POINTER :: NumElPartBndryMat


  INTEGER, DIMENSION(:,:),POINTER :: NodeFlag_str
  INTEGER, DIMENSION(:,:),POINTER :: NodeFlag_mm
  INTEGER, DIMENSION(:,:),POINTER :: NodeFlag_th

  INTEGER, DIMENSION(:), POINTER :: NumNeighProcs_List

  INTEGER, DIMENSION(:), POINTER :: NodesToCommunicate,ID_sendto_List
  INTEGER :: NodesToCommunicate_cnt

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: BC_Flag
  
  INTEGER, DIMENSION(:), ALLOCATABLE :: NumBC_Flag
  INTEGER, POINTER, DIMENSION(:) :: Pconn_Comm

  INTEGER :: MaxNumBC_str,MaxNumBC_mm,MaxNumBC_th, NumSerBC

  INTEGER, POINTER, DIMENSION(:,:) :: BC_values_mm, BC_values_str
  INTEGER, POINTER, DIMENSION(:) :: BC_values_th
  REAL*8, POINTER, DIMENSION(:) :: BCValue

  LOGICAL :: InteractMesh
  LOGICAL :: OverlayMesh

! Node numbering for implicit solver
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NumNP_loc_implicit, StartNumNP_loc_implicit
  INTEGER, ALLOCATABLE, DIMENSION(:) :: MapNodeImp, NodeProcImpGlobal
  INTEGER, POINTER, DIMENSION(:) :: NodeNumGlobalImp
  INTEGER, POINTER, DIMENSION(:) :: NodeProcImp
  LOGICAL :: IMP
! end


END MODULE meshdata

