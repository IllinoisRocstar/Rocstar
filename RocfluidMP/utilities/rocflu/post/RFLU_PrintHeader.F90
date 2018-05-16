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
! Purpose: Print header including version number and date.
!
! Description: None.
!
! Input: 
!   global      Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_PrintHeader.F90,v 1.11 2008/12/06 08:45:06 mtcampbe Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_PrintHeader(global)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters

  USE ModInterfaces, ONLY: BuildVersionString

#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY: PLAG_BuildVersionString
#endif
  
#ifdef SPEC
  USE ModInterfacesSpecies, ONLY: SPEC_BuildVersionString
#endif
  
  IMPLICIT NONE

! ******************************************************************************
! Definitions and declarations
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: headerString,RCSIdentString,versionString
  INTEGER, PARAMETER :: headerWidth = 38
  INTEGER :: margin,versionWidth

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PrintHeader.F90,v $ $Revision: 1.11 $'

  CALL RegisterFunction(global,'RFLU_PrintHeader', & 
                        'RFLU_PrintHeader.F90')

! ==============================================================================
! Build version string
! ==============================================================================

  CALL BuildVersionString(versionString)

! ==============================================================================
! Build header string
! ==============================================================================

  headerString = ' '
  
  versionWidth = LEN_TRIM(versionString)
  margin       = (headerWidth - versionWidth)/2 ! Note integer division
  
  headerString(margin+1:margin+versionWidth) = versionString(1:versionWidth)

! ==============================================================================
! Print header 
! ==============================================================================

  WRITE(STDOUT,'(A)')      SOLVER_NAME
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'======================================'
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'                                      '
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'              rflupost                '
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'                                      '
  
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,TRIM(headerString)
 
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,' Copyright (c) University of Illinois '
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'                                      '
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'======================================'
  WRITE(STDOUT,'(A)')      SOLVER_NAME

! ==============================================================================
! Write out messages from conditional compilation 
! ==============================================================================
  
#ifdef PLAG
  WRITE(STDOUT,'(A)') SOLVER_NAME
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Compiled with particle module'
  CALL PLAG_BuildVersionString(headerString)
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,TRIM(headerString) 
#endif

#ifdef SPEC
  WRITE(STDOUT,'(A)') SOLVER_NAME
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Compiled with species module'
  CALL SPEC_BuildVersionString(headerString)
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,TRIM(headerString) 
#endif

  WRITE(STDOUT,'(A)') SOLVER_NAME
    
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PrintHeader

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PrintHeader.F90,v $
! Revision 1.11  2008/12/06 08:45:06  mtcampbe
! Updated license.
!
! Revision 1.10  2008/11/19 22:18:17  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.9  2006/04/07 15:19:26  haselbac
! Removed tabs
!
! Revision 1.8  2005/05/03 03:12:11  haselbac
! Changed header
!
! Revision 1.7  2004/12/02 15:28:39  haselbac
! Added printing of module version strings, cosmetics
!
! Revision 1.6  2003/11/25 21:03:56  haselbac
! Cosmetic changes only
!
! Revision 1.5  2003/03/20 20:07:19  haselbac
! Modified RegFun call to avoid probs with
! long 'RFLU_PrintHeader.F90' names
!
! Revision 1.4  2002/09/09 16:42:28  haselbac
! global and mixtInput now under regions
!
! Revision 1.3  2002/07/25 15:23:32  haselbac
! Added output for CHECK_DATASTRUCT
!
! Revision 1.2  2002/06/17 13:37:35  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.1  2002/06/10 21:38:42  haselbac
! Initial revision
!
! ******************************************************************************








