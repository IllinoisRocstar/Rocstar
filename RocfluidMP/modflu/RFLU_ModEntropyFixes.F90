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
! Purpose: Collection definitions of entropy fixes.
!
! Description: None
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_ModEntropyFixes.F90,v 1.4 2008/12/06 08:44:21 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

MODULE RFLU_ModEntropyFixes

  USE ModDataTypes

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: EntropyFixHartenHyman

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
  
! ==============================================================================
! Local variables
! ==============================================================================
 
  CHARACTER(CHRLEN), PARAMETER :: & 
    RCSIdentString = '$RCSfile: RFLU_ModEntropyFixes.F90,v $ $Revision: 1.4 $'

! ******************************************************************************
! Module subroutines
! ******************************************************************************

  CONTAINS
  
    FUNCTION EntropyFixHartenHyman(l,d)
      
      REAL(RFREAL), INTENT(IN) :: l,d
      REAL(RFREAL) :: EntropyFixHartenHyman
      
      IF ( l > d ) THEN 
        EntropyFixHartenHyman = l
      ELSE  
        EntropyFixHartenHyman = (l*l + d*d)/(2.0_RFREAL*d)
      END IF ! l
      
    END FUNCTION EntropyFixHartenHyman

! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLU_ModEntropyFixes


! ******************************************************************************
!
! RCS Revision history:
!
!   $Log: RFLU_ModEntropyFixes.F90,v $
!   Revision 1.4  2008/12/06 08:44:21  mtcampbe
!   Updated license.
!
!   Revision 1.3  2008/11/19 22:17:32  mtcampbe
!   Added Illinois Open Source License/Copyright
!
!   Revision 1.2  2004/01/22 16:03:59  haselbac
!   Made contents of modules PRIVATE, only procs PUBLIC, to avoid errors on ALC and titan
!
!   Revision 1.1  2003/11/25 21:03:30  haselbac
!   Initial revision
!
! ******************************************************************************






