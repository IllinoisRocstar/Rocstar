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
!!     implicit_global
!!
!!  FUNCTION
!!     The variables within this module are used throughout the program.
!!     Many of them are communications arrays initialized with InitComm1
!!     and InitComm2.  Others are local and global mappings and numbers 
!!     of nodes, which are both initialized when the mesh is read in.
!!
!!  INPUTS
!!     none
!!
!!  OUTPUTS
!!     none
!!
!!  USES
!!     none
!!
!!****

MODULE implicit_global

  USE Precision
  USE ROCSTAR_RocFrac
  
  IMPLICIT NONE

!
! Node Numbering
!
     INTEGER :: LNumNp  ! Number of local nodes
     INTEGER :: GNumNp  ! Number of global nodes
     INTEGER, ALLOCATABLE, DIMENSION(:) :: Local2Global  ! Local to global mapping of node numbers
     INTEGER, ALLOCATABLE, DIMENSION(:) :: Global2Local  ! Global to local mapping of node numbers

!
! Communication Constructs
!
     INTEGER, ALLOCATABLE, DIMENSION(:) :: NodeProc  ! Assignment of nodes to a processor
     INTEGER :: NumCommProcs1, NumCommProcs2  ! Number of processors that must be communicated with
     INTEGER, ALLOCATABLE, DIMENSION(:) :: CommProcs1, CommProcs2  ! Processors that must be communicated with
     INTEGER, ALLOCATABLE, DIMENSION(:) :: NumCommNodes1, NumCommNodes2  ! Number of nodes that need to be communicated to a processor
     INTEGER, ALLOCATABLE, DIMENSION(:,:) :: CommNodes1, CommNodes2  ! Nodes to be communicated to a processor
     INTEGER :: MaxNumCommNodes1, MaxNumCommNodes2  ! Maximum number of nodes to be sent to a processor
     INTEGER :: NumCommProcsFrom1, NumCommProcsFrom2  ! Number of procs to recieve data from
     INTEGER, ALLOCATABLE, DIMENSION(:) :: CommProcsFrom1, CommProcsFrom2  ! Procs to recieve data from
     INTEGER, ALLOCATABLE, DIMENSION(:) :: NumCommNodesFrom1, NumCommNodesFrom2  ! Number of nodes to be received from processors

!
! MPI Communication Variables
!
     INTEGER :: myid ! this processor's ID
     INTEGER :: nprocs ! number of procs available
     INTEGER :: ierr ! MPI error code
     TYPE(rcv_buf), ALLOCATABLE, DIMENSION(:) :: frmproc ! multi-proc receive buffer
     INTEGER,ALLOCATABLE,DIMENSION(:) :: req_rcv, req_snd  ! send and receive requests
     INTEGER,ALLOCATABLE,DIMENSION(:,:) :: stat_rcv, stat_snd  ! send and receive statuses

!
! Nodal BC Flags for each DOF
!
     INTEGER,ALLOCATABLE,DIMENSION(:,:) :: node_flag
     REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: boundary_value


!!!!!!!!!!!!!!!!!!!!!!! ... Structural Solver Data Structure ... !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! External forces from input file
!
     REAL(wp),ALLOCATABLE,DIMENSION(:) :: fext_imp

!
! Storage for matrices and associated variables
!

     INTEGER :: nstart_km, nrows_km  ! Dimensions of parts of M  and K matrices

     INTEGER :: nnz_m  ! Number of nonzeros in the mass matrix
     REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: aval_m  ! M matrix in compressed row storage
     INTEGER,ALLOCATABLE,DIMENSION(:) :: rp_m, cval_m  ! M matrix in compressed row storage

     INTEGER :: nnz_k  ! number of nonzeros in the K matrix
     REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: aval_k  ! K matrix in compressed row storage
     INTEGER,ALLOCATABLE,DIMENSION(:) :: rp_k, cval_k  ! K matrix in compressed row storage
     
     REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: aval_meff  ! Effective M matrix in compressed row storage
     INTEGER,ALLOCATABLE,DIMENSION(:) :: cval_meff, rp_meff  ! Effective M matrix in compressed row storage
     INTEGER :: nnz_meff, nstart_meff, nrows_meff  ! Dimensions of parts of effective M matrix

!
! Keep track of whether the initial acceleration has been solved for
!
     
     LOGICAL :: initAccel

!!!!!!!!!!!!!!!!!!!!!!! ... Thermal Solver Data Structure ... !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! ... External thermal loads from input file
     REAL(wp),ALLOCATABLE,DIMENSION(:) :: rext_imp

     ! ... Thermal load vector from previous timestep
     ! ... for use in the beta method in thermal_soln.f90
     REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: r_old

     ! ... Storage for matrices and associated variables

     ! ... Dimensions of parts of the thermal stiffness [Kt] and thermal capacitance [C] matrices 
     INTEGER :: nstart_ktc, nrows_ktc

     ! ... Thermal Capacitance Matrix [C]
     ! ... Number of nonzeros in the thermal capacitance matrix 
     INTEGER :: nnz_c
     ! ... Capacitance matix in compressed row storage
     REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: aval_c 
     INTEGER,ALLOCATABLE,DIMENSION(:) :: rp_c, cval_c 

     ! ... Thermal Stiffness Matrix [Kt]
     ! ... Number of nonzeros in the thermal stiffness matrix 
     INTEGER :: nnz_kt
     ! ... Thermal stiffness matix in compressed row storage
     REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: aval_kt
     INTEGER,ALLOCATABLE,DIMENSION(:) :: rp_kt, cval_kt
     
     ! ... Effective Thermal Capacitance Matrix [C_eff]
     ! ... Number of nonzeros in the effective capacitance matrix
     REAL(kind=wp),ALLOCATABLE,DIMENSION(:) :: aval_ceff
     ! ... Effective capacitance matrix in compressed row storage 
     INTEGER,ALLOCATABLE,DIMENSION(:) :: cval_ceff, rp_ceff
     INTEGER :: nnz_ceff, nstart_ceff, nrows_ceff

     REAL, ALLOCATABLE, DIMENSION(:) :: Ceff_fpTP

END MODULE implicit_global
