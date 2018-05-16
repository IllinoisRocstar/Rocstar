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
! Purpose: Read in user input related to preprocessor.
!
! Description: None.
!
! Input: 
!   global      Pointer to global type
!
! Output: None.
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: ReadPrepSection.F90,v 1.4 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReadPrepSection(global)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModInterfaces, ONLY: ReadSection
  USE ModError
  USE ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER, PARAMETER :: NVALS_MAX = 1  

  LOGICAL :: defined(NVALS_MAX)
  CHARACTER(10) :: keys(NVALS_MAX)
  CHARACTER(CHRLEN) :: RCSIdentString    
  INTEGER :: nVals  
  REAL(RFREAL) :: vals(NVALS_MAX)

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: ReadPrepSection.F90,v $ $Revision: 1.4 $'

  CALL RegisterFunction(global,'ReadPrepSection',&
  'ReadPrepSection.F90')

! ******************************************************************************
! Specify keywords and search for them
! ******************************************************************************

#ifdef RFLU
  nVals = NVALS_MAX

  keys(1) = 'PARTMODE'

  CALL ReadSection(global,IF_INPUT,nVals,keys,vals,defined)

  IF ( defined(1) .EQV. .TRUE. ) THEN
    IF ( NINT(vals(1)) == PARTITION_MODE_PROPER ) THEN 
      global%prepPartMode = PARTITION_MODE_PROPER
    ELSE 
      global%prepPartMode = PARTITION_MODE_IMPOSED
    END IF ! NINT(vals)
  END IF ! defined
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE ReadPrepSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadPrepSection.F90,v $
! Revision 1.4  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.1  2004/12/01 16:50:41  haselbac
! Initial revision after changing case
!
! Revision 1.7  2004/10/19 19:25:47  haselbac
! Removed SURFFLAG, cosmetics
!
! Revision 1.6  2003/08/07 15:29:01  haselbac
! Changed var names
!
! Revision 1.5  2003/05/16 02:27:43  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.4  2003/05/07 00:21:40  haselbac
! Changed logic for check
!
! Revision 1.3  2003/05/01 14:06:40  haselbac
! Fixed bug: Did not update properly
!
! Revision 1.2  2003/04/29 21:47:37  haselbac
! Added SURFFLAG
!
! Revision 1.1  2003/04/29 14:55:48  haselbac
! Initial revision
!
! ******************************************************************************







