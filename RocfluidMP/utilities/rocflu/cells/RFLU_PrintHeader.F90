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
! Purpose: Print header including version number and date.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: RFLU_PrintHeader.F90,v 1.3 2008/12/06 08:44:54 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RFLU_PrintHeader(global)

  USE ModDataTypes
  USE ModError
  USE ModGlobal, ONLY: t_global
  USE ModParameters
  
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... loop variables


! ... local variables
  CHARACTER(CHRLEN) :: headerString,RCSIdentString,versionString
  INTEGER, PARAMETER :: headerWidth = 53
  INTEGER :: margin,versionWidth

!******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_PrintHeader.F90,v $'

  CALL RegisterFunction(global,'RFLU_PrintHeader', &
                        'RFLU_PrintHeader.F90')

! Build version string --------------------------------------------------------

  CALL BuildVersionString(versionString)

! Build header string ---------------------------------------------------------

  headerString = ' '
  
  versionWidth = LEN_TRIM(versionString)
  margin       = (headerWidth - versionWidth)/2 ! Note integer division
  
  headerString(margin+1:margin+versionWidth) = versionString(1:versionWidth)
    
  headerString(1:1) = '*'
  headerString(headerWidth:headerWidth) = '*'

! Print header ----------------------------------------------------------------

  WRITE(STDOUT,'(A)')      SOLVER_NAME
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'*****************************************************'
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'*                                                   *'
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'*         ROCFLU-MP: Special cell converter         *'
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'*         =================================         *'
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'*                                                   *'
  
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,TRIM(headerString)
 
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'*    Copyright (c) by the University of Illinois    *'
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'*                                                   *'
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'*****************************************************'
  WRITE(STDOUT,'(A)')      SOLVER_NAME
    
! Done ------------------------------------------------------------------------

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_PrintHeader

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_PrintHeader.F90,v $
! Revision 1.3  2008/12/06 08:44:54  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:18:04  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2003/04/01 17:02:46  haselbac
! Initial revision
!
! Revision 1.2  2003/03/20 19:53:36  haselbac
! Modified RegFun call to avoid probs with
! long 'RFLU_PrintHeader.F90' names
!
! Revision 1.1.1.1  2003/03/19 17:22:22  haselbac
! Initial revision
!
!******************************************************************************








