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
!   global              Pointer to global type
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_PrintHeader.F90,v 1.19 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
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
  
#ifdef TURB
  USE ModInterfacesTurbulence, ONLY: TURB_BuildVersionString
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

  RCSIdentString = '$RCSfile: RFLU_PrintHeader.F90,v $ $Revision: 1.19 $'

  CALL RegisterFunction(global,'RFLU_PrintHeader',&
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
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'======================================='
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'                                       '
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'               RocfluidMP              '
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'                                       '
  
!  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,TRIM(headerString)
 
  WRITE(STDOUT,'(A,1X,A,A)') SOLVER_NAME,'Copyright (C) 2015 Illinois Rocstar',&
                                       ' LLC.'
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'                                       '
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'======================================='
  WRITE(STDOUT,'(A)')      SOLVER_NAME
  
! ==============================================================================
! Write out messages from conditional compilation 
! ==============================================================================
  
#ifdef STATS
  WRITE(STDOUT,'(A)') SOLVER_NAME
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'###############################'
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Compiled with statistics module'
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'###############################'
#endif

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

#ifdef TURB
  WRITE(STDOUT,'(A)') SOLVER_NAME
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Compiled with turbulence module'
  CALL TURB_BuildVersionString(headerString)
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
! Revision 1.19  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.18  2008/11/19 22:17:43  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.17  2005/05/03 03:06:11  haselbac
! Changed header
!
! Revision 1.16  2005/04/20 14:43:34  haselbac
! Removed CHECK_UNIFLOW code section
!
! Revision 1.15  2004/12/02 15:26:38  haselbac
! Added printing of module version strings, removed some stuff
!
! Revision 1.14  2004/11/03 17:05:34  haselbac
! Removed HACK_PERIODIC ifdef, cosmetics
!
! Revision 1.13  2004/02/02 22:51:44  haselbac
! Added writing of header for particles
!
! Revision 1.12  2003/11/25 21:04:43  haselbac
! Added support for rocspecies, cosmetic changes
!
! Revision 1.11  2003/03/15 18:52:45  haselbac
! Cosmetics only
!
! Revision 1.10  2002/10/27 19:16:30  haselbac
! Removed tabs
!
! Revision 1.9  2002/09/09 15:51:56  haselbac
! global now under region
!
! Revision 1.8  2002/07/25 14:33:20  haselbac
! Added messages for checking of gradients and periodic hack
!
! Revision 1.7  2002/07/24 17:31:08  wasistho
! Added Rocturb header
!
! Revision 1.6  2002/06/27 16:30:32  haselbac
! Added statement for CHECK_DATASTRUCT
!
! Revision 1.5  2002/06/17 16:09:45  wasistho
! Added STATS compilation message
!
! Revision 1.4  2002/06/17 13:54:24  haselbac
! Cosmetic change, added SOLVER_NAME, added message for STATS
!
! Revision 1.3  2002/06/17 13:34:12  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.2  2002/06/14 21:54:35  wasistho
! Added time avg statistics
!
! Revision 1.1  2002/06/10 21:28:09  haselbac
! Initial revision
!
! ******************************************************************************







