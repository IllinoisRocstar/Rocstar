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
! Purpose: Set explicit interfaces to common subroutines and functions for
!   Rocflu.
!
! Description: None
!
! Notes: 
!   1. The subroutines contained in this interface file are common to at least
!      two Rocflu modules in the sense that these modules contain subroutines
!      with these names, but the actual source is NOT common. 
!
! ******************************************************************************
!
! $Id: RFLU_ModInterfacesCommon.F90,v 1.6 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************
  
MODULE RFLU_ModInterfacesCommon

  IMPLICIT NONE

  INTERFACE

  SUBROUTINE RFLU_AllocateMemory(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_AllocateMemory

  SUBROUTINE RFLU_AllocateMemoryWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_AllocateMemoryWrapper
  
  SUBROUTINE RFLU_DeallocateMemory(pRegion)
    USE ModDataStruct, ONLY: t_region 
    TYPE(t_region), POINTER :: pRegion    
  END SUBROUTINE RFLU_DeallocateMemory

  SUBROUTINE RFLU_DeallocateMemoryWrapper(pRegion)
    USE ModDataStruct, ONLY: t_region
    TYPE(t_region), POINTER :: pRegion
  END SUBROUTINE RFLU_DeallocateMemoryWrapper
  
  SUBROUTINE RFLU_PrintHeader(global)
    USE ModGlobal, ONLY: t_global
    TYPE(t_global), POINTER :: global
  END SUBROUTINE RFLU_PrintHeader

  END INTERFACE

END MODULE RFLU_ModInterfacesCommon

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModInterfacesCommon.F90,v $
! Revision 1.6  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2004/10/19 19:28:04  haselbac
! Adapted interface of RFLU_AllocateMemoryWrapper
!
! Revision 1.3  2004/03/19 21:19:21  haselbac
! Removed interfaces for alloc/dealloc routines
!
! Revision 1.2  2004/02/26 21:02:01  haselbac
! Added/deleted entries for memory allocation due to PLAG integration
!
! Revision 1.1  2003/04/10 14:37:10  haselbac
! Initial revision
!
! ******************************************************************************






