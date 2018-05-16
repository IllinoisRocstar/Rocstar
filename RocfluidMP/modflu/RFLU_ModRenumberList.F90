! *********************************************************************
! * Rocstar Simulation Suite                                          *
! * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
! *                                                                   *
! * Illinois Rocstar LLC                                              *
! * Champaign, IL                                                     *
! * www.illinoisrocstar.com                                           *
! * sales@illinoisrocstar.com                                         *
! *                                                                   *
! * License: See LICENSE file in top level of distribution package or *
! * http://opensource.org/licenses/NCSA                               *
! *********************************************************************
! *********************************************************************
! * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
! * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
! * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
! * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
! * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
! * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
! * Arising FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
! * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
! *********************************************************************
! ******************************************************************************
!
! Purpose: Suite of routines to renumber and denumber lists.
!
! Description: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_ModRenumberList.F90,v 1.6 2008/12/06 08:44:23 mtcampbe Exp $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModRenumberList

  USE ModDataTypes
  USE ModParameters
  USE ModError
  USE ModGlobal, ONLY: t_global

  USE ModSortSearch

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: RFLU_DenumberList, & 
            RFLU_RenumberList, & 
            RFLU_RenumberList2       
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
        
  CHARACTER(CHRLEN) :: & 
    RCSIdentString = '$RCSfile: RFLU_ModRenumberList.F90,v $ $Revision: 1.6 $' 
              
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS





! ******************************************************************************
!
! Purpose: Denumber list.
!
! Description: None.
!
! Input:
!   global      Global pointer
!   listDim1    Leading dimension of list
!   listDim2    Trailing dimension of list
!   list        List to be denumbered
!   keyDim      Dimension of key
!   key         Key for renumbering
!
! Output: 
!   list        Denumbered list
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_DenumberList(global,listDim1,listDim2,list,keyDim,key)    

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================
  
    INTEGER, INTENT(IN) :: listDim1,listDim2,keyDim
    INTEGER, INTENT(INOUT) :: list(listDim1,listDim2)
    INTEGER, INTENT(IN) :: key(keyDim)
    TYPE(t_global), POINTER :: global
  
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: i,j,k

! ******************************************************************************
!   Start
! ******************************************************************************  

    CALL RegisterFunction(global,'RFLU_DenumberList',&
  'RFLU_ModRenumberList.F90')

! ******************************************************************************
!   Denumber list
! ******************************************************************************  

    DO j = 1,listDim2
      DO i = 1,listDim1
        k = list(i,j)

        IF ( k /= VERT_NONE ) THEN         
          IF ( k > keyDim ) THEN 
            CALL ErrorStop(global,ERR_DENUMBER_LIST,__LINE__)
          ELSE         
            list(i,j) = key(k)
          END IF ! k
        END IF ! k
      END DO ! i 
    END DO ! j

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_DenumberList




! ******************************************************************************
!
! Purpose: Renumber list.
!
! Description: None.
!
! Input:
!   global      Global pointer
!   listDim1    Leading dimension of list
!   listDim2    Trailing dimension of list
!   list        List to be renumbered
!   keyDim      Dimension of key
!   key         Key for renumbering
!
! Output: 
!   list        Renumbered list
!
! Notes: 
!   1. IMPORTANT: key is assumed to be in ascending order - otherwise binary
!      search will fail.
!
! ******************************************************************************

  SUBROUTINE RFLU_RenumberList(global,listDim1,listDim2,list,keyDim,key)    

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================
  
    INTEGER, INTENT(IN) :: listDim1,listDim2,keyDim
    INTEGER, INTENT(INOUT) :: list(listDim1,listDim2)
    INTEGER, INTENT(IN) :: key(keyDim)
    TYPE(t_global), POINTER :: global
  
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: i,iLoc,j,k

! ******************************************************************************
!   Start
! ******************************************************************************  

    CALL RegisterFunction(global,'RFLU_RenumberList',&
  'RFLU_ModRenumberList.F90')

! ******************************************************************************
!   Renumber list
! ******************************************************************************  

    DO j = 1,listDim2
      DO i = 1,listDim1
        k = list(i,j)

        IF ( k /= VERT_NONE ) THEN 
          CALL BinarySearchInteger(key,keyDim,k,iLoc)

          IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
            list(i,j) = iLoc
          ELSE 
            CALL ErrorStop(global,ERR_BINARY_SEARCH,__LINE__)
          END IF ! iLoc  
        END IF ! k
      END DO ! i 
    END DO ! j

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RenumberList








! ******************************************************************************
!
! Purpose: Renumber list based on two keys.
!
! Description: None.
!
! Input:
!   global      Global pointer
!   listDim1    Leading dimension of list
!   listDim2    Trailing dimension of list
!   list        List to be renumbered
!   keyDim      Dimension of key
!   key1        Key for renumbering
!   key2        Second key for renumbering
!
! Output: 
!   list        Renumbered list
!
! Notes: 
!   1. IMPORTANT: key1 is assumed to be in ascending order - otherwise binary
!      search will fail.
!
! ******************************************************************************

  SUBROUTINE RFLU_RenumberList2(global,listDim1,listDim2,list,keyDim,key1, & 
                                key2)    

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************  

! ==============================================================================
!   Arguments
! ==============================================================================
  
    INTEGER, INTENT(IN) :: listDim1,listDim2,keyDim
    INTEGER, INTENT(INOUT) :: list(listDim1,listDim2)
    INTEGER, INTENT(IN) :: key1(keyDim),key2(keyDim)
    TYPE(t_global), POINTER :: global
  
! ==============================================================================
!   Locals
! ==============================================================================

    INTEGER :: i,iLoc,j,k

! ******************************************************************************
!   Start
! ******************************************************************************  

    CALL RegisterFunction(global,'RFLU_RenumberList',&
  'RFLU_ModRenumberList.F90')

! ******************************************************************************
!   Renumber list
! ******************************************************************************  

    DO j = 1,listDim2
      DO i = 1,listDim1
        k = list(i,j)

        IF ( k /= VERT_NONE ) THEN                 
          CALL BinarySearchInteger(key1,keyDim,k,iLoc)

          IF ( iLoc /= ELEMENT_NOT_FOUND ) THEN 
            list(i,j) = key2(iLoc)
          ELSE 
            CALL ErrorStop(global,ERR_BINARY_SEARCH,__LINE__)
          END IF ! iLoc  
        END IF ! k
      END DO ! i 
    END DO ! j

! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_RenumberList2








! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModRenumberList


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModRenumberList.F90,v $
! Revision 1.6  2008/12/06 08:44:23  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:34  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2006/04/07 15:19:20  haselbac
! Removed tabs
!
! Revision 1.3  2006/03/25 21:56:15  haselbac
! Cosmetics only
!
! Revision 1.2  2004/12/04 03:32:20  haselbac
! Added new renumbering routine, used for partitioning
!
! Revision 1.1  2004/10/19 19:27:02  haselbac
! Initial revision
!
! ******************************************************************************









