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
! Purpose: Build unique patch identifier for GENX runs.
!
! Description: None.
!
! Input:
!   iRegion     Region index
!   iPatch      Patch index
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: BuildPatchIdentifier.F90,v 1.4 2008/12/06 08:44:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

INTEGER FUNCTION BuildPatchIdentifier(iRegion,iPatch)

  USE ModDataTypes
  USE ModError
  USE ModParameters
  
  IMPLICIT NONE

! ... parameters
  INTEGER, INTENT(IN) :: iPatch,iRegion

! ... loop variables


! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

!******************************************************************************

  RCSIdentString = '$RCSfile: BuildPatchIdentifier.F90,v $ $Revision: 1.4 $'

! start -----------------------------------------------------------------------

  BuildPatchIdentifier = iRegion*REGION_INDEX_OFFSET + iPatch

! end -------------------------------------------------------------------------

END FUNCTION BuildPatchIdentifier

!******************************************************************************
!
! RCS Revision history:
!
! $Log: BuildPatchIdentifier.F90,v $
! Revision 1.4  2008/12/06 08:44:08  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 16:48:16  haselbac
! Initial revision after changing case
!
! Revision 1.1  2002/10/27 18:46:19  haselbac
! Initial revision
!
!******************************************************************************






