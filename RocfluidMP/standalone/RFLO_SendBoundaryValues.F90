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
! Purpose: send boundary values to an external source.
!
! Description: none.
!
! Input: region     = dimensions and topology
!        initialize = initial data to external program (true/false).
!
! Output: whatever is needed.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: RFLO_SendBoundaryValues.F90,v 1.3 2008/12/06 08:44:45 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLO_SendBoundaryValues( region,initialize )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

  LOGICAL :: initialize

!******************************************************************************

  CALL RegisterFunction( region%global,'RFLO_SendBoundaryValues',&
  'RFLO_SendBoundaryValues.F90' )

! give error message

  IF (region%mixtInput%externalBc) THEN
    CALL ErrorStop( region%global,ERR_EXTERNAL_FUNCT,__LINE__ )
  ENDIF

! finalize

  CALL DeregisterFunction( region%global )

END SUBROUTINE RFLO_SendBoundaryValues

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLO_SendBoundaryValues.F90,v $
! Revision 1.3  2008/12/06 08:44:45  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:57  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:29:37  haselbac
! Initial revision after changing case
!
! Revision 1.5  2003/07/03 21:48:45  jblazek
! Implemented dual-time stepping.
!
! Revision 1.4  2003/05/15 02:57:06  jblazek
! Inlined index function.
!
! Revision 1.3  2002/09/20 22:22:37  jblazek
! Finalized integration into GenX.
!
! Revision 1.2  2002/09/05 17:40:22  jblazek
! Variable global moved into regions().
!
! Revision 1.1  2002/08/15 19:43:11  jblazek
! Prepared solver to be optionally linked with an external driver.
!
!******************************************************************************







