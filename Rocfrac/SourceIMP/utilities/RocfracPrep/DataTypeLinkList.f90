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
MODULE DataTypeLinkList

  IMPLICIT NONE

! -- Start -- A Processor's Element List

  TYPE ProcElemList_data_ptr
     INTEGER :: GlbElNum
     TYPE(ProcElemList_data_ptr), POINTER :: next
  END TYPE ProcElemList_data_ptr
  
  TYPE :: ProcElemList_data_type
     TYPE(ProcElemList_data_ptr), POINTER :: head
     TYPE(ProcElemList_data_ptr), POINTER :: tail
  endtype ProcElemList_data_type
  
  TYPE(ProcElemList_data_type), TARGET, ALLOCATABLE, DIMENSION(:) :: ProcElemList
  
  TYPE(ProcElemList_data_ptr), POINTER :: ProcElem_Item
  
! -- End

! -- Start -- Boundary Conditions

  TYPE BC_ptr
     INTEGER :: BC_nodeGlb
     INTEGER :: BC_flagGlb
     TYPE(BC_ptr), POINTER :: next
  END TYPE BC_ptr

  TYPE(BC_ptr), POINTER :: BC_structural_head
  TYPE(BC_ptr), POINTER :: BC_structural_tail
  
  TYPE(BC_ptr), POINTER :: BC_meshmotion_head
  TYPE(BC_ptr), POINTER :: BC_meshmotion_tail
  
  TYPE(BC_ptr), POINTER :: BC_thermal_head
  TYPE(BC_ptr), POINTER :: BC_thermal_tail
  
  TYPE(BC_ptr), POINTER :: BC_structural_item
  TYPE(BC_ptr), POINTER :: BC_meshmotion_item
  TYPE(BC_ptr), POINTER :: BC_thermal_item
  
! -- End

! -- Start Surface Mesh

 ! -- 3 node triangle ( 4 node tet)
  TYPE SurfMesh_tri3_ptr
     INTEGER, DIMENSION(1:5) :: ElemData
     TYPE(SurfMesh_tri3_ptr), POINTER :: next
  END TYPE SurfMesh_tri3_ptr

  TYPE(SurfMesh_tri3_ptr), POINTER :: SurfMesh_tri3_S_head
  TYPE(SurfMesh_tri3_ptr), POINTER :: SurfMesh_tri3_S_tail
  TYPE(SurfMesh_tri3_ptr), POINTER :: SurfMesh_tri3_SF_head
  TYPE(SurfMesh_tri3_ptr), POINTER :: SurfMesh_tri3_SF_tail
  TYPE(SurfMesh_tri3_ptr), POINTER :: SurfMesh_tri3_SF_NonIgnt_head
  TYPE(SurfMesh_tri3_ptr), POINTER :: SurfMesh_tri3_SF_NonIgnt_tail
  TYPE(SurfMesh_tri3_ptr), POINTER, DIMENSION(:) :: SurfMesh_tri3_Ov_head 
!  TYPE(SurfMesh_tri3_ptr), POINTER, DIMENSION(:) :: SurfMesh_tri3_Ov_tail
  TYPE(SurfMesh_tri3_ptr), POINTER :: SurfMesh_tri3_Ov1_head, SurfMesh_tri3_Ov1_tail
  TYPE(SurfMesh_tri3_ptr), POINTER :: SurfMesh_tri3_Ov2_head, SurfMesh_tri3_Ov2_tail

 ! -- 6 node triangle ( 10 node tet)
  TYPE SurfMesh_tri6_ptr
     INTEGER, DIMENSION(1:8) :: ElemData
     TYPE(SurfMesh_tri6_ptr), POINTER :: next
  END TYPE SurfMesh_tri6_ptr

  TYPE(SurfMesh_tri6_ptr), POINTER :: SurfMesh_tri6_S_head
  TYPE(SurfMesh_tri6_ptr), POINTER :: SurfMesh_tri6_S_tail
  TYPE(SurfMesh_tri6_ptr), POINTER :: SurfMesh_tri6_SF_head
  TYPE(SurfMesh_tri6_ptr), POINTER :: SurfMesh_tri6_SF_tail
  TYPE(SurfMesh_tri6_ptr), POINTER :: SurfMesh_tri6_SF_NonIgnt_head
  TYPE(SurfMesh_tri6_ptr), POINTER :: SurfMesh_tri6_SF_NonIgnt_tail

! -- 4 node quad ( 8 node hex)
  TYPE SurfMesh_hex8_ptr
     INTEGER, DIMENSION(1:6) :: ElemData
     TYPE(SurfMesh_hex8_ptr), POINTER :: next
  END TYPE SurfMesh_hex8_ptr
  
  TYPE(SurfMesh_hex8_ptr), POINTER :: SurfMesh_hex8_S_head
  TYPE(SurfMesh_hex8_ptr), POINTER :: SurfMesh_hex8_S_tail
  TYPE(SurfMesh_hex8_ptr), POINTER :: SurfMesh_hex8_SF_head
  TYPE(SurfMesh_hex8_ptr), POINTER :: SurfMesh_hex8_SF_tail
  TYPE(SurfMesh_hex8_ptr), POINTER :: SurfMesh_hex8_SF_NonIgnt_head
  TYPE(SurfMesh_hex8_ptr), POINTER :: SurfMesh_hex8_SF_NonIgnt_tail
  
! add item

  TYPE(SurfMesh_tri3_ptr), POINTER :: SurfMesh_tri3_item
  TYPE(SurfMesh_tri6_ptr), POINTER :: SurfMesh_tri6_item
  TYPE(SurfMesh_hex8_ptr), POINTER :: SurfMesh_hex8_item

  
END MODULE DataTypeLinkList

