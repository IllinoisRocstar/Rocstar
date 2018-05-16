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
PROGRAM RocfracPrep
  
  USE CommGlobal
  USE meshdata
  USE Linked_List2
  USE linked_list
  
  IMPLICIT NONE

  INCLUDE 'comf90.h'

!     Type definition 

  TYPE nodebc
     INTEGER :: node,bc
  END TYPE nodebc
  
  TYPE nodebc_lst
     INTEGER :: node_lst,bc_lst
  END TYPE nodebc_lst

  TYPE sndrcv_buf
     INTEGER :: sndn,rcvn,nodes
  END TYPE sndrcv_buf

  TYPE(sndrcv_buf), ALLOCATABLE, DIMENSION(:) :: sndrcvnod
  TYPE(sndrcv_buf), ALLOCATABLE, DIMENSION(:) :: sndrcvnod_lst

!-- number of nodes for triangle
  INTEGER :: numvertx2d
!-- dummy variables
  INTEGER :: i,ii,j,jj,k,kk
  INTEGER :: iaux,iaux1,iaux2,iaux3,n,mm
  REAL*8 :: aux
  INTEGER :: ntime
  
  INTEGER :: nk
!-- Stores the material number associated with lst
  INTEGER, ALLOCATABLE, DIMENSION(:) :: iaux89
!-- Stores the cohesive material number associated with cohesive
  INTEGER :: flag
  
  INTEGER :: numbcss
  INTEGER :: indx,jndx
  INTEGER :: edgecut,nn,nprocs,iunit
  CHARACTER*3 :: ai4
  CHARACTER*4 :: ichr4,ai1
  INTEGER :: i1,i2,i3
  INTEGER :: NumNeighProcs
  

  INTEGER, ALLOCATABLE, DIMENSION(:) :: ncoor

!-- Element Neighbor array from triangle
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lmtri_2d
!-- Connectivity array (OUTPUT) for cohesive elements
!----- dimension:  <1>  6 nodes   <2> local element id  <3> processor
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: lmcoh
!-- Tempory holding array of cohesive and no-cohesive (input to METIS)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: elmnts
!-- Tempory holding array of partitioned 'elmnts' array
  INTEGER, ALLOCATABLE, DIMENSION(:) :: epart_p
!-- Number of nodal points & elements on each processor
  INTEGER, ALLOCATABLE, DIMENSION(:) :: numnp
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: numel
!-- Number of cohesive elements on each processor
  INTEGER, ALLOCATABLE, DIMENSION(:) :: numclst
!-- For the lst: relates the old node numbering to the new nodes

  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jk1
!-- Number of nodal b.c. for each processor 
  INTEGER, ALLOCATABLE, DIMENSION(:) :: numbc
!-- Number of nodal mesh motionb.c. for each processor 
  INTEGER, ALLOCATABLE, DIMENSION(:) :: numbc_mm
!-- Number of nodal mesh motionb.c. for each processor 
  INTEGER, ALLOCATABLE, DIMENSION(:) :: numbc_ht
!-- Number of cohesive elements on processor boundary for each processor
  INTEGER, ALLOCATABLE, DIMENSION(:) :: num_border_coh
!-- Number of neighboring processors for R_co
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nproc_neigh
!-- Number of neighboring processors for R_in
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nproc_neigh_lst
!-- Number of neighboring processors for R_co on No/cohesive boundary
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nproc_neigh_rco
!-- Stores the neighboring processors that I'm communcating with.
  INTEGER, ALLOCATABLE, DIMENSION(:) :: my_neigh
!-- Stores if node is on the processor boundary
  INTEGER, ALLOCATABLE, DIMENSION(:) :: bord_node_org
!--
  INTEGER :: err
  INTEGER, ALLOCATABLE, DIMENSION(:) :: matclst,lmtemp
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nnntemp_lst
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nnntemp_rco
  
  TYPE proclist
     INTEGER, DIMENSION(:), POINTER :: proc_list
  END TYPE proclist
  TYPE proclistrco
     INTEGER, DIMENSION(:), POINTER :: proc_list_rco
  END TYPE proclistrco
!-- wave properties and stability varibles
  REAL*8 :: xx,yy,zz,size1,size2,size3,size4,size5,size6
  REAL*8 :: dhmin,dt_courant
  REAL*8 :: dh_courant,cd_courant
  INTEGER :: maxdim
!-

  INTEGER :: numnp_total,nface

  INTEGER, ALLOCATABLE,DIMENSION(:,:) :: NumNp2D
  INTEGER, ALLOCATABLE,DIMENSION(:,:,:) :: jk1_2D
  INTEGER, ALLOCATABLE,DIMENSION(:,:,:) :: ik1_2D
  
  CHARACTER*3 :: ai3
  CHARACTER*1 ichr1
  CHARACTER*12 ichr12
  INTEGER :: num_zones,numcoh_zone,numelv_prmry_zone
  INTEGER, ALLOCATABLE, DIMENSION(:) :: numelv_prmry_zones, ik1_z
  INTEGER, ALLOCATABLE, DIMENSION(:) :: MapGlbEl2LocEl
  
  INTEGER :: numel_z,numnp_z,p2max_lst,num_rco_border,ip,icount

  INTEGER iaux23
  INTEGER :: num_vol,num_coh
  INTEGER :: iflag
  INTEGER, DIMENSION(1:4) :: iflag89
  INTEGER :: ntri
  INTEGER :: ios
  INTEGER :: jk1_size
  
  TYPE(coh_info_type), POINTER :: coh_item
  TYPE(coh_list_type), TARGET, ALLOCATABLE, DIMENSION(:) :: coh_list
  TYPE(coh_comm_info_type), POINTER :: coh_comm_item
  TYPE(coh_comm_list_type), TARGET, ALLOCATABLE, DIMENSION(:,:) :: coh_comm_list
  TYPE(vol_info_type), POINTER :: vol_item
  TYPE(vol_list_type), TARGET, ALLOCATABLE, DIMENSION(:) :: vol_list
  TYPE(vol_box_point), ALLOCATABLE, DIMENSION(:) :: ik1_c4
  
  INTEGER ni,node1,node2,node3,node4,gnode1,gnode2,gnode3,gnode4
  
  TYPE(bcvalues), DIMENSION(16) :: bc_mshmtn
  
  
  INTEGER :: numploadelem
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: idpressload
  INTEGER, ALLOCATABLE, DIMENSION(:) :: ipressflag
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: pressload
  
  LOGICAL :: ElOnPartBndry
  
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: icnt
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NumBorderVol
  
  INTEGER :: myid
  REAL*8 :: shift
  INTEGER :: iargc
  CHARACTER*3 chr_arg
  CHARACTER*4 chr_procs
  CHARACTER*20 chr_units
  
  INTEGER :: imat, itmp1, itmp2
  INTEGER, POINTER :: tmpptr
  
  INTEGER, DIMENSION(1:10) :: NdHistoryFlag
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NumNdHistoryP
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NdHistoryP
  INTEGER :: GlbNd
  INTEGER :: ielem, ElemCount
  
  TYPE(Link_Ptr_Type)  :: Link
  
  TYPE(ProcElemList_data_ptr), POINTER :: ptr2
  TYPE(BC_ptr), POINTER :: ptr_BC
  
  CHARACTER(*), PARAMETER :: volWin = "vfrac"

  REAL*8, POINTER, DIMENSION(:,:) :: MeshCoor

  INTEGER :: MaxNumNodesComm, icounter

  INTEGER :: write_attr, set_option, vol_all, errFlg, comp_pconn
  
  ConvertUnit = 1.d0
 
  n = iargc()
  IF(n.LE.1) GOTO 554
  
  j = 1
  DO i = 1, n/2
     CALL getarg( j, chr_arg)
     IF(chr_arg.EQ.'-np')THEN
        CALL getarg( j+1, chr_procs)
        READ(chr_procs,*) nprocs
        j = j +2
     ELSE IF(chr_arg.EQ.'-un')THEN ! units
        CALL getarg( j+1, chr_units)
        READ(chr_units,*) ConvertUnit
        j = j + 2
        PRINT*,'UNITS CONVERSION'
        PRINT*,' multiply by ', ConvertUnit
        
     ELSE
        GOTO 554
     ENDIF
  ENDDO
  
  GOTO 555

554 CONTINUE
  PRINT*,'Usage:'
  PRINT*,'   MeshTran -np #'
  PRINT*,'         - where # is the number of processors'
  STOP

555 CONTINUE

  myid = 0
  
  CALL readinp(ntime)
  
  PRINT*,'numvertx=',numvertx
  
  !      CALL system('\rm -f -r '//prefx(1:prefx_lngth))


! -- mesh motion boundary conditions

            
  dhmin = 1000000000.d0

  IF(iansys.EQ.1)THEN
     PRINT*,'ANSYS NO LONGER SUPPORTED'
     PRINT*,'STOPPING'
     STOP
!        CALL read_ansys(numbcss, dhmin)
  ELSE IF(ipatran.EQ.1)THEN
     CALL read_patran(numvertx2d,dhmin,nprocs)
     numbcss = numnp_prmry
!$$$      ELSE IF(ipatcohin.EQ.1)THEN
!$$$         CALL read_patran_cohin(numbcss,numvertx2d,dhmin,nprocs)
!$$$      ELSE IF(itetcohin.EQ.1)THEN
!$$$         CALL read_tetmesh_cohin(numbcss,numvertx2d,dhmin,nprocs)
!$$$         numbcss = numnp_prmry
!$$$      ELSE IF(itetmesh.EQ.1)THEN
!NOTE:
!IF both ascii and binary temesh input files exist,then
! the binary file will be the one that gets read.
!$$$         CALL read_tetmesh(numbcss,numvertx2d,dhmin,nprocs)
!$$$         numbcss = numnp_prmry
  ENDIF
!
! -- Checking Courant condition for time step
!
  dt_courant = dhmin/cd_fastest
  cd_courant = cd_fastest
  dh_courant = dhmin
  
  PRINT*,' --------------------------------------'
  PRINT*,' ---- COURANT STABILITY CONDITION -----'
  PRINT*,' --------------------------------------'
  PRINT*,'   dt_courant =',dt_courant
  PRINT*,'   cd_courant =',cd_courant
  PRINT*,'    h_courant =',dh_courant
  PRINT*,' --------------------------------------'

  iaux = 1
  PRINT*,'Allocate(numel)'
  ALLOCATE(numel(1:NumMat,1:nprocs))
  
  PRINT*,'Allocate(numel).finish'
  numel(:,:) = 0
  
  DO n = 1, numelv_prmry    ! Loop over the tetrahedra elements
     imat = MatId(n)
     numel(imat,epart(n)) = numel(imat,epart(n)) + 1 
  ENDDO
  
  ALLOCATE(numnp(1:nprocs))
  
  ii  = 0
  
  numnp(1:nprocs) = 0
  numel(1:NumMat,1:nprocs) = 0
  ALLOCATE(vol_item)
  ALLOCATE(vol_list(1:NumMat))


!
! -- RENUMBER THE VOLUMETRIC ELEMENTS USING LOCAL NUMBERING
  PRINT*,'RENUMBERING VOLUMETRIC ELEMENTS'

!      allocate(MapGlbEl2LocEl(1:numelv_prmry))
!      allocate(icnt(1:NumMat,1:2))
!      icnt(:,:) = 0
!      allocate(ElOnPartBndry(1:numelv_prmry))
      
  ALLOCATE(ProcNodeList(1:nprocs))
  
  DO i = 1, nprocs
     CALL LI_Init_List(ProcNodeList(i))
  ENDDO
  
  ALLOCATE(NodeFlag(1:numnp_prmry))
  ALLOCATE(ElFlag(1:numelv_prmry))
  
!
!    Write each processor's input file
!
  ALLOCATE(nproc_neigh_lst(1:nprocs)) ! does not need to keep track of processors
  nproc_neigh_lst(:) = 0



  CALL COM_Init

! - load Rocout module
  CALL COM_set_verbose( 10)
  CALL SimOUT_load_module( 'OUT')

! Surface boundary meshes

!!$  OPEN(4001,FILE=prefx(1:prefx_lngth)//'/fracSF.im',STATUS='replace',FORM='formatted')
!!$  OPEN(4002,FILE=prefx(1:prefx_lngth)//'/fracS.im',STATUS='replace',FORM='formatted')
!!$  
!!$  WRITE(4001,*) nprocs,MeshType2D
!!$  WRITE(4002,*) nprocs,MeshType2D


  OPEN(4005,FILE='Rocin/isolid_in_00.000000.txt',STATUS='replace',FORM='formatted')
  WRITE(4005,*) '@Proc: *'
  !WRITE(4005,*) '@Files: Rocfrac/Rocin/SurfMesh.%4p.hdf'
  WRITE(4005,*) '@Files: Rocfrac/Rocin/SurfMesh.%4p.cgns'
  WRITE(4005,*) '@Panes: @BlockCyclic 100 100'
  close(4005)
  
  OPEN(4005,FILE='Rocin/solid_in_00.000000.txt',STATUS='replace',FORM='formatted')
  WRITE(4005,*) '@Proc: *'
  !WRITE(4005,*) '@Files: Rocfrac/Rocin/'//prefx(1:prefx_lngth)//'.%4p.hdf'
  WRITE(4005,*) '@Files: Rocfrac/Rocin/'//prefx(1:prefx_lngth)//'.%4p.cgns'
  WRITE(4005,*) '@Panes: @Cyclic 1'
  close(4005)




! -------------------------------------
! MAIN LOOP OVER PROCESSORS
! ------------------------------------------------

  DO ip = 1, nprocs

     CALL COM_new_window( volWin )

! -- Initialize link list variables

     vol_list(1:NumMat)%num_border_vol = 0
     DO i = 1, NumMat
        NULLIFY(vol_list(i)%vol_head)
        NULLIFY(vol_list(i)%vol_tail)
     ENDDO

     PRINT*,'Processor id =',ip
     
     WRITE(ichr4,'(i4.4)') ip - 1


! --------------------------
! Nodes 
!-----------------------------------

! Renumber nodes locally, keep track of which node already renumbered
     NodeFlag(:) = 0
     ElFlag(:) = 0

! go through the Processor's element link list


     ptr2 => ProcElemList(ip)%head

     
     iaux = 0
     
     DO WHILE(ASSOCIATED(ptr2))
        
        ielem = ptr2%GlbElNum
        
        imat = MatId(ielem)
        numel(imat,ip) = numel(imat,ip) + 1
        ElOnPartBndry  = .FALSE.
        DO k = 1, ElTypeId(ielem)
           nk = lmelv_prmry(k,ielem)
           
           IF(NumProcPerNd(nk).GT.1) ElOnPartBndry =.TRUE.
           
           IF(NodeFlag(nk).EQ.0)THEN
              numnp(ip) = numnp(ip) + 1
              iaux = iaux + 1
              NodeFlag(nk) = numnp(ip)
           ENDIF
           
           vol_item%mat_vol = imat
           vol_item%lmvol(k) = NodeFlag(nk)
           
        ENDDO

        
        ptr2 => ptr2%next
        
        vol_item%iface = 0
        vol_item%press = ielem ! 0.d0

!     Add item to volumetric element list
!
            
        IF(ElOnPartBndry)THEN ! Element contains a node that is on the partition boundary
           CALL vol_insert_head(vol_list(imat),vol_item)
!               icnt(imat,1) = icnt(imat,1) + 1
!               ElFlag(ielem) = icnt(imat,1)
!               MapGlbEl2LocEl(ip) = icnt(imat,1)
        ELSE
           CALL vol_insert_tail(vol_list(imat),vol_item)
!               icnt(imat,2) = icnt(imat,2) + 1
!               ElFlag(ielem) = icnt(imat,2)
!               MapGlbEl2LocEl(ip) = icnt(imat,2)
        ENDIF

     ENDDO


     ALLOCATE(MeshCoor(1:3,1:NumNP(ip))) !**

     DO i = 1, numnp_prmry
        IF(NodeFlag(i).NE.0) MeshCoor(1:3,NodeFlag(i)) = coor(1:3,i)
     END DO

     CALL COM_new_dataitem( volWin//'.nc', 'n', COM_DOUBLE, 3, 'm')
     CALL COM_set_size( volWin//'.nc', ip, NumNdPerProc(ip) )
     CALL COM_set_array(volWin//'.nc', ip, MeshCoor, 3)

     PRINT*,'registerd Coordinates'

     IF(iaux.NE.NumNdPerProc(ip))THEN
        PRINT*,'Test failed when renumbering'
        PRINT*,'iaux=',iaux
        PRINT*,'NumNdPerProc',NumNdPerProc(ip)
        STOP
     ENDIF

! BOUNDARY CONDITIONS
! 10000*MeshMotionBC + 100*ThermalBC + StructuralBC

     PRINT*,'Boundary Conditions'

     IF(NumBC_Flag(ip).NE.0)THEN

        ALLOCATE( NodeFlag_str(1:2,1:NumBC_Flag(ip)) )


        iaux = 0

        DO i = 1, numnp_prmry
           IF(NodeFlag(i).NE.0.AND.BC_Flag(1,i).NE.0)THEN
              iaux = iaux + 1
              NodeFlag_str(1,iaux) = NodeFlag(i)
              NodeFlag_str(2,iaux) = BC_Flag(1,i)
           ENDIF
           IF(NodeFlag(i).NE.0.AND.BC_Flag(2,i).NE.0)THEN
              iaux = iaux + 1
              NodeFlag_str(1,iaux) = NodeFlag(i)
              NodeFlag_str(2,iaux) = BC_Flag(2,i)
           ENDIF
           IF(NodeFlag(i).NE.0.AND.BC_Flag(3,i).NE.0)THEN
              iaux = iaux + 1
              NodeFlag_str(1,iaux) = NodeFlag(i)
              NodeFlag_str(2,iaux) = BC_Flag(3,i)
           ENDIF
        enddo

        IF(iaux.NE.NumBC_Flag(ip))THEN
           PRINT*,'ERROR, number of BCs in partitioned mesh',iaux
           PRINT*,'    inconsistant with that of serial mesh',NumBC_Flag(ip)
           PRINT*, 'stopping'
           stop
        ENDIF

        CALL COM_new_dataitem( volWin//'.bcnode', 'p', COM_INTEGER, 2, '')
        CALL COM_set_size( volWin//'.bcnode', ip, NumBC_Flag(ip) )    
        CALL COM_set_array(volWin//'.bcnode', ip, NodeFlag_str, 2)

        
        allocate(BCValue(1:NumBC_Flag(ip)*6))
        BCValue(:) = 0.d0
        CALL COM_new_dataitem( volWin//'.BCValue', 'p', COM_DOUBLE, 1, '')
        CALL COM_set_size( volWin//'.BCValue', ip,NumBC_Flag(ip)*6 )
        CALL COM_set_array(volWin//'.BCValue', ip, BCValue, 1)
     ENDIF

!!$     CALL COM_new_dataitem( volWin//'.NumBC_str', 'p', COM_INTEGER, 1, '')
!!$     CALL COM_set_size( volWin//'.NumBC_str', ip, 1 )    
!!$     CALL COM_set_array(volWin//'.NumBC_str', ip, NumBC_structural(ip), 1)
!!$
!!$     CALL COM_new_dataitem( volWin//'.NumBC_mm', 'p', COM_INTEGER, 1, '')
!!$     CALL COM_set_size( volWin//'.NumBC_mm', ip, 1 )    
!!$     CALL COM_set_array(volWin//'.NumBC_mm', ip, NumBC_meshmotion(ip), 1)
!!$
!!$     CALL COM_new_dataitem( volWin//'.NumBC_th', 'p', COM_INTEGER, 1, '')
!!$     CALL COM_set_size( volWin//'.NumBC_th', ip, 1 )    
!!$     CALL COM_set_array(volWin//'.NumBC_th', ip, NumBC_thermal(ip), 1)

!!$! ------------------------------
!!$!
!!$! STRUCTURAL BOUNDARY CONDITIONS
!!$! --------------------------------------------
!!$
!!$!     WRITE(4000,*) 3
!!$!     WRITE(4000,*) NumBC_structural(ip),0
!!$
!!$     IF(NumBC_structural(ip).NE.0)THEN
!!$        
!!$        
!!$        ALLOCATE( NodeFlag_str(1:2,1:NumBC_structural(ip)) )
!!$        
!!$        ptr_BC => BC_structural_head
!!$        iaux = 0
!!$        DO WHILE(ASSOCIATED(ptr_BC))
!!$           
!!$           GlbNd = ptr_BC%BC_nodeGlb
!!$           iflag = ptr_BC%BC_flagGlb
!!$           IF(NodeFlag(GlbNd).NE.0)THEN
!!$              iaux = iaux + 1
!!$              NodeFlag_str(1,iaux) = NodeFlag(GlbNd)
!!$              NodeFlag_str(2,iaux) = iflag
!!$           ENDIF
!!$           ptr_BC => ptr_BC%next
!!$        ENDDO
!!$        
!!$        IF(iaux.NE.NumBC_structural(ip))THEN
!!$           PRINT*,'ERROR, number of structural BCs in linked list'
!!$           PRINT*,'    inconsistant with that of serierial mesh'
!!$           PRINT*, 'stopping'
!!$        ENDIF
!!$
!!$        CALL COM_new_dataitem( volWin//'.NodeFlag_str', 'p', COM_INTEGER, 2, '')
!!$        CALL COM_set_size( volWin//'.NodeFlag_str', ip, NumBC_structural(ip) )    
!!$        CALL COM_set_array(volWin//'.NodeFlag_str', ip, NodeFlag_str, 2)
!!$     ENDIF
!!$
!!$         
!!$! ------------------------------
!!$!
!!$! MESH MOTION BOUNDARY CONDITIONS
!!$! --------------------------------------------
!!$
!!$!     WRITE(4000,*) 4
!!$!     WRITE(4000,*) NumBC_meshmotion(ip),0
!!$
!!$     IF(NumBC_meshmotion(ip).NE.0)THEN
!!$        
!!$        ALLOCATE( NodeFlag_mm(1:2,1:NumBC_meshmotion(ip)) )
!!$        
!!$        
!!$        ptr_BC => BC_meshmotion_head
!!$        iaux = 0
!!$        DO WHILE(ASSOCIATED(ptr_BC))
!!$           
!!$           GlbNd = ptr_BC%BC_nodeGlb
!!$           iflag = ptr_BC%BC_flagGlb
!!$           IF(NodeFlag(GlbNd).NE.0)THEN
!!$              iaux = iaux + 1
!!$              NodeFlag_mm(1,iaux) = NodeFlag(GlbNd)
!!$              NodeFlag_mm(2,iaux) = iflag
!!$           ENDIF
!!$           ptr_BC => ptr_BC%next
!!$        ENDDO
!!$        
!!$        IF(iaux.NE.NumBC_meshmotion(ip))THEN
!!$           PRINT*,'ERROR, number of mesh motion BCs in linked list'
!!$           PRINT*,'    inconsistant with that of serierial mesh'
!!$           PRINT*, 'stopping'
!!$        ENDIF
!!$        
!!$        CALL COM_new_dataitem( volWin//'.NodeFlag_mm', 'p', COM_INTEGER, 2, '')
!!$        CALL COM_set_size( volWin//'.NodeFlag_mm', ip, NumBC_meshmotion(ip) )    
!!$        CALL COM_set_array(volWin//'.NodeFlag_mm', ip, NodeFlag_mm, 2)
!!$
!!$     ENDIF
!!$
!!$! ------------------------------
!!$!
!!$! THERMAL BOUNDARY CONDITIONS
!!$! --------------------------------------------
!!$
!!$!     WRITE(4000,*) 8
!!$!     WRITE(4000,*) NumBC_thermal(ip),0
!!$
!!$     IF(NumBC_thermal(ip).NE.0)THEN
!!$        
!!$        ALLOCATE( NodeFlag_th(1:2, 1:NumBC_thermal(ip) ) )
!!$        
!!$        ptr_BC => BC_thermal_head
!!$        iaux = 0
!!$        DO WHILE(ASSOCIATED(ptr_BC))
!!$           
!!$           GlbNd = ptr_BC%BC_nodeGlb
!!$           iflag = ptr_BC%BC_flagGlb
!!$           IF(NodeFlag(GlbNd).NE.0)THEN
!!$              iaux = iaux + 1
!!$              NodeFlag_th(1,iaux) = NodeFlag(GlbNd)
!!$              NodeFlag_th(2,iaux) = iflag
!!$           ENDIF
!!$           ptr_BC => ptr_BC%next
!!$        ENDDO
!!$        
!!$        IF(iaux.NE.NumBC_thermal(ip))THEN
!!$           PRINT*,'ERROR, number of thermal BCs in linked list'
!!$           PRINT*,'    inconsistant with that of serierial mesh'
!!$           PRINT*, 'stopping'
!!$        ENDIF
!!$        
!!$        CALL COM_new_dataitem( volWin//'.NodeFlag_th', 'p', COM_INTEGER, 2, '')
!!$        CALL COM_set_size( volWin//'.NodeFlag_th', ip, NumBC_thermal(ip) )    
!!$        CALL COM_set_array(volWin//'.NodeFlag_th', ip, NodeFlag_th, 2)
!!$
!!$     ENDIF

        
! ----------------------------------------------------
! --  WRITE VOLUMETRIC ELEMENT CONNECTIVITY ARRAY 
! ------------------------------------------------
!
! No. of 4-node tetrahedral
! No. of 10-node tetrahedral
! No. of lst on the partioned mesh boundary


     itmp1 = SUM(numel(1:NumMat,ip))
     itmp2 = SUM(vol_list(1:NumMat)%num_border_vol)
     

     ALLOCATE(NumElVolMat(1:NumMat),NumElPartBndryMat(1:NumMat))

     ALLOCATE(matType(1:itmp1))
     ALLOCATE(ElConnTable(1:numvertx,itmp1))

     ! ... ElemCount added to keep track of elements put into this processors 
     ! ... connectivity table (ElConnTable in Print_vol_list) regardless of the 
     ! ... material.  COstoich 10/27/09
     ElemCount = 0
     DO ii = 1, NumMat
        
!        WRITE(4000,'(6i9)') itmp1,itmp2,numel(ii,ip),vol_list(ii)%num_border_vol,numvertx,0
        
        NumElVolMat(ii) = numel(ii,ip)
        NumElPartBndryMat(ii) = vol_list(ii)%num_border_vol

        CALL Print_vol_list(vol_list(ii),ii,NumElVolMat(ii),ElemCount)
        ElemCount = ElemCount + numel(ii,ip)
     ENDDO

     CALL COM_new_dataitem( volWin//'.NumElPartBndry', 'p', COM_INTEGER, 1, '')
     CALL COM_set_size( volWin//'.NumElPartBndry', ip, 1)
     CALL COM_resize_array(volWin//'.NumElPartBndry', ip, tmpptr)
     tmpptr = itmp2

     CALL COM_new_dataitem( volWin//'.NumElVolMat', 'p', COM_INTEGER, 1, '')
     CALL COM_set_size( volWin//'.NumElVolMat', ip, NumMat)    
     CALL COM_set_array(volWin//'.NumElVolMat', ip, NumElVolMat, 1)

     CALL COM_new_dataitem( volWin//'.NumElPartBndryMat', 'p', COM_INTEGER, 1, '')
     CALL COM_set_size( volWin//'.NumElPartBndryMat', ip, NumMat)    
     CALL COM_set_array(volWin//'.NumElPartBndryMat', ip, NumElPartBndryMat, 1)
!
! Registering Element Connectivity
!
       IF(numvertx.EQ.4)THEN

!!$       CALL COM_init_mesh( volWin//'.T4', MyId+1, glb%ElConnVol, glb%NumElVol)

          CALL COM_set_size( volWin//'.:T4', ip, itmp1)
          CALL COM_set_array( volWin//'.:T4', ip,  ElConnTable, 4)
     

       ELSE IF(numvertx.EQ.10)THEN

!!$       CALL COM_init_mesh( volWin//'.T10', MyId+1, glb%ElConnVol, glb%NumElVol)

          CALL COM_set_size( volWin//'.:T10', ip, itmp1)
          CALL COM_set_array( volWin//'.:T10', ip, ElConnTable, 10)

       ELSE IF(numvertx.EQ.8)THEN
!!$       CALL COM_init_mesh( volWin//'.H8', MyId+1, glb%ElConnVol, glb%NumElVol)

          CALL COM_set_size( volWin//'.:H8', ip, itmp1)
          CALL COM_set_array( volWin//'.:H8', ip,  ElConnTable, 8)

       ENDIF

! Element Material Type Flag
       CALL COM_new_dataitem( volWin//'.MatType', 'e', COM_INTEGER, 1, '')
       CALL COM_set_array(volWin//'.MatType', ip, matType, 1)

       


! --------------------------------------------
! --  WRITE MPI COMMUNICATION INFORMATION
! --------------------------------------------------


!     Nodeal Force calculaton communciation

!     WRITE(4000,*) 6

!     Determine the neighbor of processors 'i' is communicating with.
     
       NumNeighProcs = 0
       MaxNumNodesComm = 0
       DO j = 1,nprocs
          IF(ID_sendto(ip,j)%num_border_comm.NE.0) THEN
             NumNeighProcs = NumNeighProcs + 1
             MaxNumNodesComm = MaxNumNodesComm + ID_sendto(ip,j)%num_border_comm
          ENDIF
       ENDDO

!     Number of neighboring proc. involved in R_in calculation

!     WRITE(4000,*) NumNeighProcs

! 1D array format: 
!
!         (1) number of communicating panes
!         (2) <communicating pane id> 
!         (3) # shared node between panes (3) List of nodes .. repeat

!     ALLOCATE(NumNeighProcs_List(1:NumNeighProcs))
!     ALLOCATE(NodesToCommunicate(1:MaxNumNodesComm))
!     ALLOCATE(ID_sendto_List(1:NumNeighProcs))

!     List these neighboring processors

!should this be moved, to you always need to register
       ALLOCATE(Pconn_Comm(1:NumNeighProcs*2+MaxNumNodesComm+1))
    !   print*,'NumNeighProcs*2+MaxNumNodesComm',NumNeighProcs*2+MaxNumNodesComm
    !   print*,'MaxNumNodesComm',MaxNumNodesComm

       icounter = 1
       Pconn_Comm(icounter) = NumNeighProcs

     !NodesToCommunicate_cnt = 0
       DO j=1,nprocs          ! receiving processor
          IF(ID_sendto(j,ip)%num_border_comm.NE.0)THEN
           !  print*,'ID_sendto(j,ip)%num_border_comm',ID_sendto(j,ip)%num_border_comm
!     Number of nodes that need to be communicated for R_in calculation


           

             icounter = icounter+1
             Pconn_Comm(icounter) = j

             icounter = icounter+1
             Pconn_Comm(icounter) = ID_sendto(j,ip)%num_border_comm ! why plus  + 1
             
!     List of nodes that need to be communicated for R_in calculation
             CALL print_comm_list(ID_sendto(j,ip),ip,icounter)
             
          ENDIF
       ENDDO

            
       IF(icounter.NE.NumNeighProcs*2+MaxNumNodesComm+1)THEN
          PRINT*,'ERROR in Communication Pack array'
          STOP
       ENDIF


       CALL COM_new_dataitem( volWin//'.pconn', 'p', COM_INTEGER, 1, '')
       CALL COM_set_size( volWin//'.pconn', ip, NumNeighProcs*2+MaxNumNodesComm+1)    
       CALL COM_set_array(volWin//'.pconn', ip, Pconn_Comm, 1)
        
!!$        CALL COM_new_dataitem( volWin//'.ID_sendto_List', 'p', COM_INTEGER, 1, '')
!!$        CALL COM_set_size( volWin//'.ID_sendto_List', ip, NumNeighProcs)    
!!$        CALL COM_set_array(volWin//'.ID_sendto_List', ip, ID_sendto_List, 1)
!!$        
!!$        CALL COM_new_dataitem( volWin//'.NumNeighProcs_List', 'p', COM_INTEGER, 1, '')
!!$        CALL COM_set_size( volWin//'.NumNeighProcs_List', ip, NumNeighProcs)    
!!$        CALL COM_set_array(volWin//'.NumNeighProcs_List', ip, NumNeighProcs_List, 1)
   
   
       CALL COM_window_init_done( volWin)

     ! Load Rocmap using a name "MyPrivateMAP"
!!$     CALL COM_load_module( "Rocmap", "MyPrivateMAP")
!!$
!!$     ! Call compute_pconn
!!$     comp_pconn = COM_get_function_handle( "MyPrivateMAP.compute_pconn")
!!$     CALL COM_call_function( comp_pconn, 2, &
!!$          COM_get_dataitem_handle_const( volWin//'.mesh'), &
!!$          COM_get_dataitem_handle(volWin//'.pconn'))
!!$
!!$! Unload Rocmap.
!!$     CALL COM_unload_module( "Rocmap", "MyPrivateMAP") 

! obtain function handle ------------------------------------------------------

     write_attr = COM_get_function_handle( 'OUT.write_dataitem')
     set_option = COM_get_function_handle( 'OUT.set_option')

     CALL COM_call_function( set_option, 2, 'mode', 'w')
     ! Masoud: switching to HDF4 for this module
     !CALL COM_call_function( set_option, 2, 'format', 'HDF4')
     ! End

! do not append process rank -----------------

     CALL COM_call_function( set_option, 2, 'rankwidth', '0')
! write volume window ------------------------
     vol_all = Com_get_dataitem_handle( volWin//'.all')

     CALL COM_call_function( write_attr, 4, 'Rocin/'//prefx(1:prefx_lngth)//'.'//ichr4, vol_all,&
          "solid","00.000000")

! delete volume  window ----------

     CALL COM_delete_window( volWin)


     DEALLOCATE(MeshCoor)
  
     IF(ASSOCIATED(NodeFlag_str)) DEALLOCATE(NodeFlag_str)
     IF(ASSOCIATED(NodeFlag_mm)) DEALLOCATE(NodeFlag_mm)
     IF(ASSOCIATED(NodeFlag_th)) DEALLOCATE(NodeFlag_th)

     IF(ASSOCIATED(BCValue)) DEALLOCATE(BCValue)

     DEALLOCATE(NumElVolMat,NumElPartBndryMat,matType,ElConnTable)
     !DEALLOCATE(NodesToCommunicate)
     !deallocate(NumNeighProcs_List,ID_sendto_List)

     IF(ASSOCIATED(Pconn_Comm)) deallocate(Pconn_Comm)

     IF(InteractMesh) CALL mesh2d(nprocs,ip,ichr4)

     PRINT*,'mesh2dOverlay'
     IF(OverlayMesh) THEN
        IF(ip.EQ.1) THEN
           OPEN(456,FILE ='Rocin/OverlayMappings.txt')
        ENDIF
        CALL mesh2dOverlay(nprocs,ip,ichr4)
        PRINT*,'Finsihed mesh2dOverlay', ip
     ENDIF
  ENDDO

  CLOSE(4001)
  CLOSE(4002)


!--START

!!$  OPEN(4001,FILE=prefx(1:prefx_lngth)//'/fracSF.im',STATUS='replace',FORM='formatted')
!!$  OPEN(4002,FILE=prefx(1:prefx_lngth)//'/fracS.im',STATUS='replace',FORM='formatted')
!!$  
!!$  WRITE(4001,*) nprocs,MeshType2D
!!$  WRITE(4002,*) nprocs,MeshType2D
!!$      
!!$  DO ip = 1, nprocs
!!$
!!$
!!$     CALL COM_new_window( volWin )
!!$
!!$
!!$  !   CALL RocstarInitSolution( gridLevel,iReg,regions,wins,winv )
!!$
!!$   !  CALL RocstarWriteSolution( gridLevel,iReg,regions(iReg),wins,winv )
!!$
!!$
!!$
!!$! -- Initialize link list variables
!!$
!!$     vol_list(1:NumMat)%num_border_vol = 0
!!$     DO i = 1, NumMat
!!$        NULLIFY(vol_list(i)%vol_head)
!!$        NULLIFY(vol_list(i)%vol_tail)
!!$     ENDDO
!!$     
!!$     PRINT*,'Processor id =',ip
!!$     
!!$     WRITE(ichr4,'(i4.4)') ip - 1
!!$
!!$! Output To Each Processors Files
!!$
!!$     OPEN(4000,FILE= &
!!$          prefx(1:prefx_lngth)//'/'//prefx(1:prefx_lngth)//'.'//ichr4//'.inp', &
!!$          STATUS='replace',FORM='formatted')
!!$
!!$! --------------------------
!!$! Version of ROCSTAR_DATA
!!$! -----------------------------------
!!$!
!!$     WRITE(4000,*) 1
!!$!     WRITE(4000,*) 2.5
!!$! --------------------------
!!$! Nodes 
!!$!-----------------------------------
!!$     WRITE(4000,*) 2
!!$     WRITE(4000,*) NumNdPerProc(ip),0,0,0,0
!!$
!!$! Renumber nodes locally, keep track of which node already renumbered
!!$     NodeFlag(:) = 0
!!$     ElFlag(:) = 0
!!$
!!$! go through the Processor's element link list
!!$
!!$     ptr2 => ProcElemList(ip)%head
!!$     
!!$     iaux = 0
!!$     
!!$     DO WHILE(ASSOCIATED(ptr2))
!!$        
!!$        ielem = ptr2%GlbElNum
!!$        
!!$        imat = MatId(ielem)
!!$        numel(imat,ip) = numel(imat,ip) + 1
!!$        ElOnPartBndry  = .FALSE.
!!$        DO k = 1, ElTypeId(ielem)
!!$           nk = lmelv_prmry(k,ielem)
!!$           
!!$           IF(NumProcPerNd(nk).GT.1) ElOnPartBndry =.TRUE.
!!$           
!!$           IF(NodeFlag(nk).EQ.0)THEN
!!$              numnp(ip) = numnp(ip) + 1
!!$              WRITE(4000,'(i9,3(1x,e16.9),2i9)') numnp(ip), coor(1:3,nk),0
!!$              iaux = iaux + 1
!!$              NodeFlag(nk) = numnp(ip)
!!$           ENDIF
!!$           
!!$           
!!$           vol_item%mat_vol = imat
!!$           vol_item%lmvol(k) = NodeFlag(nk)
!!$           
!!$        ENDDO
!!$        
!!$        ptr2 => ptr2%next
!!$        
!!$        vol_item%iface = 0
!!$        vol_item%press = ielem ! 0.d0
!!$
!!$!     Add item to volumetric element list
!!$!
!!$            
!!$        IF(ElOnPartBndry)THEN ! Element contains a node that is on the partition boundary
!!$           CALL vol_insert_head(vol_list(imat),vol_item)
!!$!               icnt(imat,1) = icnt(imat,1) + 1
!!$!               ElFlag(ielem) = icnt(imat,1)
!!$!               MapGlbEl2LocEl(ip) = icnt(imat,1)
!!$        ELSE
!!$           CALL vol_insert_tail(vol_list(imat),vol_item)
!!$!               icnt(imat,2) = icnt(imat,2) + 1
!!$!               ElFlag(ielem) = icnt(imat,2)
!!$!               MapGlbEl2LocEl(ip) = icnt(imat,2)
!!$        ENDIF
!!$
!!$     ENDDO
!!$
!!$     IF(iaux.NE.NumNdPerProc(ip))THEN
!!$        PRINT*,'Test failed when renumbering'
!!$        PRINT*,'iaux=',iaux
!!$        PRINT*,'NumNdPerProc',NumNdPerProc(ip)
!!$        STOP
!!$     ENDIF
!!$
!!$! ------------------------------
!!$!
!!$! STRUCTURAL BOUNDARY CONDITIONS
!!$! --------------------------------------------
!!$
!!$     WRITE(4000,*) 3
!!$     WRITE(4000,*) NumBC_structural(ip),0
!!$
!!$         
!!$     ptr_BC => BC_structural_head
!!$     iaux = 0
!!$     DO WHILE(ASSOCIATED(ptr_BC))
!!$        
!!$        GlbNd = ptr_BC%BC_nodeGlb
!!$        iflag = ptr_BC%BC_flagGlb
!!$        IF(NodeFlag(GlbNd).NE.0)THEN
!!$           WRITE(4000,'(4i9)') NodeFlag(GlbNd), iflag, 0
!!$           iaux = iaux + 1
!!$        ENDIF
!!$        ptr_BC => ptr_BC%next
!!$     ENDDO
!!$     
!!$     IF(iaux.NE.NumBC_structural(ip))THEN
!!$        PRINT*,'ERROR, number of structural BCs in linked list'
!!$        PRINT*,'    inconsistant with that of serierial mesh'
!!$        PRINT*, 'stopping'
!!$     ENDIF
!!$         
!!$! ------------------------------
!!$!
!!$! MESH MOTION BOUNDARY CONDITIONS
!!$! --------------------------------------------
!!$     WRITE(4000,*) 4
!!$     WRITE(4000,*) NumBC_meshmotion(ip),0
!!$     
!!$     ptr_BC => BC_meshmotion_head
!!$     iaux = 0
!!$     DO WHILE(ASSOCIATED(ptr_BC))
!!$        
!!$        GlbNd = ptr_BC%BC_nodeGlb
!!$        iflag = ptr_BC%BC_flagGlb
!!$        IF(NodeFlag(GlbNd).NE.0)THEN
!!$           WRITE(4000,'(4i9)') NodeFlag(GlbNd), iflag, 0
!!$           iaux = iaux + 1
!!$        ENDIF
!!$        ptr_BC => ptr_BC%next
!!$     ENDDO
!!$     
!!$     IF(iaux.NE.NumBC_meshmotion(ip))THEN
!!$        PRINT*,'ERROR, number of mesh motion BCs in linked list'
!!$        PRINT*,'    inconsistant with that of serierial mesh'
!!$        PRINT*, 'stopping'
!!$     ENDIF
!!$     
!!$
!!$! ------------------------------
!!$!
!!$! THERMAL BOUNDARY CONDITIONS
!!$! --------------------------------------------
!!$
!!$     WRITE(4000,*) 8
!!$     WRITE(4000,*) NumBC_thermal(ip),0
!!$     
!!$! cycle through list, no longer is nice order as before
!!$
!!$     ptr_BC => BC_thermal_head
!!$     iaux = 0
!!$     DO WHILE(ASSOCIATED(ptr_BC))
!!$        
!!$        GlbNd = ptr_BC%BC_nodeGlb
!!$        iflag = ptr_BC%BC_flagGlb
!!$        IF(NodeFlag(GlbNd).NE.0)THEN
!!$           WRITE(4000,'(4i9)') NodeFlag(GlbNd), iflag, 0
!!$           iaux = iaux + 1
!!$        ENDIF
!!$        ptr_BC => ptr_BC%next
!!$     ENDDO
!!$     
!!$     IF(iaux.NE.NumBC_thermal(ip))THEN
!!$        PRINT*,'ERROR, number of thermal BCs in linked list'
!!$        PRINT*,'    inconsistant with that of serierial mesh'
!!$        PRINT*, 'stopping'
!!$     ENDIF
!!$
!!$
!!$! ----------------------------------------------------
!!$! --  WRITE VOLUMETRIC ELEMENT CONNECTIVITY ARRAY 
!!$! ------------------------------------------------
!!$!
!!$! No. of 4-node tetrahedral
!!$! No. of 10-node tetrahedral
!!$! No. of lst on the partioned mesh boundary
!!$
!!$
!!$     itmp1 = SUM(numel(1:NumMat,ip))
!!$     itmp2 = SUM(vol_list(1:NumMat)%num_border_vol)
!!$     
!!$     WRITE(4000,*) 5
!!$     
!!$     DO ii = 1, NumMat
!!$        
!!$        WRITE(4000,'(6i9)') itmp1,itmp2,numel(ii,ip),vol_list(ii)%num_border_vol,numvertx,0
!!$        
!!$        CALL print_vol_list(vol_list(ii),ii)
!!$     ENDDO
!!$
!!$
!!$! --------------------------------------------
!!$! --  WRITE MPI COMMUNICATION INFORMATION
!!$! --------------------------------------------------
!!$
!!$!     Nodeal Force calculaton communciation
!!$
!!$     WRITE(4000,*) 6
!!$
!!$!     Determine the neighbor of processors 'i' is communicating with.
!!$     
!!$     NumNeighProcs = 0
!!$     DO j = 1,nprocs
!!$        IF(ID_sendto(ip,j)%num_border_comm.NE.0) &
!!$             NumNeighProcs = NumNeighProcs + 1
!!$     ENDDO
!!$
!!$!     Number of neighboring proc. involved in R_in calculation
!!$
!!$     WRITE(4000,*) NumNeighProcs
!!$
!!$!     List these neighboring processors
!!$
!!$     DO j=1,nprocs          ! receiving processor
!!$        IF(ID_sendto(j,ip)%num_border_comm.NE.0)THEN
!!$!     Number of nodes that need to be communicated for R_in calculation
!!$           WRITE(4000,*) j-1,ID_sendto(j,ip)%num_border_comm ! common
!!$!     List of nodes that need to be communicated for R_in calculation
!!$           CALL print_comm_list(ID_sendto(j,ip),ip)
!!$        ENDIF
!!$     ENDDO
!!$     
!!$     WRITE(4000,*) 99
!!$     CLOSE(4000)
!!$     
!!$
!!$     CALL mesh2d(nprocs,ip)
!!$       
!!$
!!$  ENDDO
!!$
!!$  CLOSE(4001)
!!$  CLOSE(4002)
      
END PROGRAM RocfracPrep

