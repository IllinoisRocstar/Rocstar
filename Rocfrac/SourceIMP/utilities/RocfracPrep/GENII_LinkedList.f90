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
! The problem is how to treat generic lists in Fortran 90. 
! Lists can be (1) homogeneous (elements all of the same type) or 
! (2) heterogeneous.
! 
! (1) Often there is a need to work with many (homogeneous) lists which
! may be of different types. Patrice Lignelet has shown how generic list
! properties may be logically treated in Fortran 90. However it still
! appears that the list operations (initialization, addition/removal of
! elements etc) must be separately defined for each list type, which
! leads to considerable duplication of code.
! 
! (2) Jean Vezina has shown how to handle a heterogeneous list by employing the 
! F90 TRANSFER() function. 
! 
! Peter McGavin at Industrial Research Limited (p.mcgavin@irl.cri.nz)
! has constructed a species of generic list for Fortran 90. The method
! is based on 2 ideas: the properties of the TRANSFER() function, and
! the fact that a pointer to a derived data type also points to the
! *first field* within the data type (and conversely). 
! 
! Since the TRANSFER() function does not accept pointer arguments the
! method requires the introduction of 2 auxilliary data types which
! contain the pointers, one in the generic list module and one for each
! list type in the calling program (the same construction is adopted
! when defining "arrays of pointers"). In spite of this complication the
! method represents a big saving, both conceptually and practically, when
! many lists of different types are involved.
! 
! To make clear the method we present a simple generic list module
! together with a calling program. The list module defines a
! uni-directional linked list with a few sample operations, but
! obviously more complicated generic lists could be substituted in its
! place (eg include back pointers).
! 
! Roger Young
! Peter McGavin
! 
! .........................................................................

MODULE Generic_List
! Defines a generic uni-directional linked list with a few sample operations

IMPLICIT NONE

PRIVATE

PUBLIC :: &
     Link_Type,        &! Put a Link_Type field first in your structure
     Link_Ptr_Type,    &! Mold this to and from your type ptr with TRANSFER
     List_Type          ! You should declare a List_Type variable

PUBLIC :: & 
     LI_Init_List,        &! Initialise the List_Type variable before use
     LI_Get_Head,         &! Returns the first Link in the list
     LI_Get_Next,         &! Return the next Link after a given one
     LI_Add_To_Head,      &! Add a Link to the head of the list
     LI_Remove_Head,      &! Remove the first Link and return it
     LI_Get_Len,          &! Compute list length
     LI_Associated,       &! Check if list member is associated
     LI_Check_List         ! Aborts program if list is invalid or corrupt

TYPE Link_Type
  PRIVATE
  TYPE(Link_Type), POINTER :: Next
END TYPE Link_Type

! Auxilliary data type required for the transfer function 
TYPE Link_Ptr_Type       ! Use TRANSFER() function to mold Link_Ptr_Type
  PRIVATE                ! to your pointer type and vice versa
  TYPE(Link_Type), POINTER :: P
END TYPE Link_Ptr_Type

TYPE List_Type             
  PRIVATE
  TYPE(Link_Type) :: Head   ! Dummy Link always at head of list
END TYPE List_Type

CONTAINS

!-----------------------------------------------------------------------
SUBROUTINE Abort(Message)
IMPLICIT NONE
CHARACTER *(*) Message

WRITE(6,*) Message
WRITE(6,*) 'Program aborted'
STOP

END SUBROUTINE Abort

!-----------------------------------------------------------------------
SUBROUTINE LI_Check_List(List,Message)
IMPLICIT NONE
TYPE(List_Type) List
CHARACTER *(*) Message

IF(.NOT.ASSOCIATED(List%Head%Next))THEN
   WRITE(6,*) Message
   CALL Abort('List is not initialised in call to LI_Check_List()')
ENDIF

END SUBROUTINE LI_Check_List

!-----------------------------------------------------------------------
SUBROUTINE LI_Init_List(List)
  IMPLICIT NONE
  integer :: NumProcs
  TYPE(List_Type),INTENT(INOUT),TARGET :: List

  NULLIFY(List%Head%Next)

  RETURN
END SUBROUTINE LI_Init_List

!-----------------------------------------------------------------------
SUBROUTINE LI_Add_To_Head(Link,List)
  IMPLICIT NONE
  TYPE(List_Type),INTENT(INOUT)     :: List
  TYPE(Link_Ptr_Type),INTENT(INOUT) :: Link

  Link%P%Next => List%Head%Next
  List%Head%Next => Link%P

  RETURN
END SUBROUTINE LI_Add_To_Head

!-----------------------------------------------------------------------
INTEGER FUNCTION LI_Get_Len(List)
  IMPLICIT NONE
  TYPE(List_Type), INTENT(IN),TARGET :: LIST
  TYPE(Link_Ptr_Type) :: Link
  INTEGER N

  Link%P => List%Head
  N = 0
  DO WHILE(ASSOCIATED(Link%P%Next))
     Link%P => Link%P%Next
     N = N+1
  ENDDO
  LI_Get_Len = N

  RETURN
END FUNCTION LI_Get_Len

!-----------------------------------------------------------------------
FUNCTION LI_Associated(Link)
  IMPLICIT NONE
  LOGICAL :: LI_Associated
  TYPE(Link_Ptr_Type),INTENT(IN) :: Link

  LI_Associated = .FALSE.
  IF(ASSOCIATED(Link%P))LI_Associated=.TRUE.

  RETURN
END FUNCTION LI_Associated

!-----------------------------------------------------------------------
FUNCTION LI_Get_Next(Link)
  IMPLICIT NONE
  Type(Link_Ptr_Type)           :: LI_Get_Next
  TYPE(Link_Ptr_Type),INTENT(IN) :: Link

  IF(.NOT.ASSOCIATED(Link%P%Next))THEN
     NULLIFY(LI_Get_Next%P)
  ELSE   
     LI_Get_Next%P => Link%P%Next
  ENDIF

  RETURN
END FUNCTION LI_Get_Next

!-----------------------------------------------------------------------
FUNCTION LI_Get_Head(List)
  IMPLICIT NONE
  TYPE(Link_Ptr_Type)               :: LI_Get_Head
  TYPE(List_Type),INTENT(IN),TARGET :: List

  LI_Get_Head%P => List%Head%Next

  RETURN
END FUNCTION LI_Get_Head

!-----------------------------------------------------------------------
FUNCTION LI_Remove_Head(List)
  IMPLICIT NONE
  TYPE(Link_Ptr_Type)                  :: LI_Remove_Head
  TYPE(List_Type),INTENT(INOUT),TARGET :: List
  TYPE(Link_Ptr_Type) :: Link

  Link%P => List%Head%Next
  IF(ASSOCIATED(Link%P))THEN
     List%Head%Next => Link%P%Next
     NULLIFY(Link%P%Next)
  ENDIF
     LI_Remove_Head%P => Link%P

  RETURN
END FUNCTION LI_Remove_Head

!-----------------------------------------------------------------------
END MODULE Generic_List


