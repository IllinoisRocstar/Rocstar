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
MODULE Linked_List2

  USE DataTypeLinkList

  IMPLICIT NONE

contains

   SUBROUTINE init_ProcElemList(head, tail)

      TYPE(ProcElemList_data_ptr), POINTER :: head, tail

      NULLIFY(head,tail)

    END SUBROUTINE init_ProcElemList

    SUBROUTINE add_ProcElemList(new, head, tail)

      TYPE(ProcElemList_data_ptr), POINTER :: new, head, tail

      IF(ASSOCIATED(head))THEN
         tail%next => new
         nullify(new%next)
         tail => new
      ELSE
         head => new
         tail => new
         nullify(tail%next)
      ENDIF
    END SUBROUTINE add_ProcElemList

SUBROUTINE list_ProcElemList(head)

  TYPE(ProcElemList_data_ptr), pointer :: head
  TYPE(ProcElemList_data_ptr), pointer :: ptr

  IF(.NOT. ASSOCIATED(head) ) THEN
     PRINT*,'list empty'
  ELSE
     ptr => head
     DO WHILE (ASSOCIATED(ptr))
        ptr => ptr%next
     ENDDO
  ENDIF
END SUBROUTINE list_ProcElemList

!! Boundary Conditions

SUBROUTINE add_BC(new, head, tail)

  TYPE(BC_ptr), POINTER :: new, head, tail

  IF(ASSOCIATED(head))THEN
     tail%next => new
     nullify(new%next)
     tail => new
  ELSE
     head => new
     tail => new
     nullify(tail%next)
  ENDIF
END SUBROUTINE add_BC

!! Surface Mesh

!! Boundary Conditions

SUBROUTINE add_SurfMesh_tri3(new, head, tail)

  TYPE(SurfMesh_tri3_ptr), POINTER :: new, head, tail

  IF(ASSOCIATED(head))THEN
     tail%next => new
     nullify(new%next)
     tail => new
  ELSE
     head => new
     tail => new
     nullify(tail%next)
  ENDIF
END SUBROUTINE add_SurfMesh_tri3

SUBROUTINE add_SurfMesh_tri6(new, head, tail)

  TYPE(SurfMesh_tri6_ptr), POINTER :: new, head, tail

  IF(ASSOCIATED(head))THEN
     tail%next => new
     nullify(new%next)
     tail => new
  ELSE
     head => new
     tail => new
     nullify(tail%next)
  ENDIF
END SUBROUTINE add_SurfMesh_tri6


SUBROUTINE add_SurfMesh_hex8(new, head, tail)

  TYPE(SurfMesh_hex8_ptr), POINTER :: new, head, tail

  IF(ASSOCIATED(head))THEN
     tail%next => new
     nullify(new%next)
     tail => new
  ELSE
     head => new
     tail => new
     nullify(tail%next)
  ENDIF
END SUBROUTINE add_SurfMesh_hex8


!!$SUBROUTINE list_ProcElemList(head)
!!$
!!$  TYPE(ProcElemList_data_ptr), pointer :: head
!!$  TYPE(ProcElemList_data_ptr), pointer :: ptr
!!$
!!$  IF(.NOT. ASSOCIATED(head) ) THEN
!!$     PRINT*,'list empty'
!!$  ELSE
!!$     ptr => head
!!$     DO WHILE (ASSOCIATED(ptr))
!!$        ptr => ptr%next
!!$     ENDDO
!!$  ENDIF
!!$END SUBROUTINE list_ProcElemList

!!$INTEGER FUNCTION Get_Len_ProcElemList(List)
!!$  IMPLICIT NONE
!!$  TYPE(ProcElemList_data_type), INTENT(IN), TARGET :: List
!!$  TYPE(ProcElemList_data_type), INTENT(IN), TARGET :: ptr
!!$  INTEGER N
!!$  
!!$  ptr => List
!!$  N = 0
!!$  DO WHILE(ASSOCIATED(ptr))
!!$     ptr => ptr%next
!!$     N = N + 1
!!$  ENDDO
!!$  Get_Len_ProcElemList = N
!!$
!!$  RETURN
!!$END FUNCTION Get_Len_ProcElemList

END MODULE Linked_List2

