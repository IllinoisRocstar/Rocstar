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
! Purpose: collect source terms pertinent to periodic flows rocperi
!
! Description: Depending on flow kind, different source terms are added to
!              the right hand side.
!
! Input: regions = data of all regions
!        iReg    = index of current region.
!
! Output: rhs updated by PERI source terms
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PERI_SourceTerms.F90,v 1.5 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_SourceTerms( region ) ! PUBLIC

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE PERI_ModInterfaces, ONLY : PERI_CoCprSlowTerms, PERI_CnlForceTerm
  USE ModError
  USE PERI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_SourceTerms.F90,v $ $Revision: 1.5 $'

  global => region%global
  CALL RegisterFunction( global,'PERI_SourceTerms',&
  'PERI_SourceTerms.F90' )

! collect PERI source terms -------- -----------------------------------------

  IF (region%periInput%flowKind == PERI_FLOW_CPR) THEN
    CALL PERI_CoCprSlowTerms( region )
  ELSEIF (region%periInput%flowKind == PERI_FLOW_CHANNEL) THEN
    CALL PERI_CnlForceTerm( region )
  ENDIF

! finalize -------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PERI_SourceTerms

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_SourceTerms.F90,v $
! Revision 1.5  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/06/11 21:49:40  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/06/09 01:21:41  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/06/08 23:56:56  wasistho
! changed nomenclature
!
! Revision 1.1.1.1  2003/03/29 03:36:30  wasistho
! install ROCPERI
!
!
!******************************************************************************







