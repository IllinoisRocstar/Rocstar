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
!******************************************************************************
!
! Purpose: place holder for MP non-dummy global communication if any.
!
! Description: none.
!
! Input: regions = data of all regions
!
! Output: global% = new global parameters after global communication
!
! Notes: none.
!
!******************************************************************************
!
! $Id: GlobalCommunicationMP.F90,v 1.3 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE GlobalCommunicationMP( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
#ifdef TURB
  USE ModInterfacesTurbulence, ONLY : TURB_LesCommunication
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  
! ... local variables
  TYPE(t_global), POINTER :: global

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'GlobalCommunicationMP',&
  'GlobalCommunicationMP.F90' )
  
! perform global communication for MP -----------------------------------------

#ifdef TURB
  IF (global%turbActive) THEN
    CALL TURB_LesCommunication( regions )
  ENDIF
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE GlobalCommunicationMP

!******************************************************************************
!
! RCS Revision history:
!
! $Log: GlobalCommunicationMP.F90,v $
! Revision 1.3  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:48:37  haselbac
! Initial revision after changing case
!
! Revision 1.1  2004/02/26 21:29:51  wasistho
! install globalCommunicationMP
!
!
!******************************************************************************







