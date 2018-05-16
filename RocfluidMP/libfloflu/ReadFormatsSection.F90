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
! Purpose: read in user input related to formats of grid and solution file.
!
! Description: none.
!
! Input: user input file.
!
! Output: global = format of grid and solution file.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ReadFormatsSection.F90,v 1.4 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadFormatsSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  INTEGER :: nVals
  INTEGER, PARAMETER :: NVALS_MAX = 3

  CHARACTER(10) :: keys(NVALS_MAX)

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL) :: vals(NVALS_MAX)

!******************************************************************************

  CALL RegisterFunction( global,'ReadFormatsSection',&
  'ReadFormatsSection.F90' )

! specify keywords and search for them

  nVals = NVALS_MAX
  
  keys(1) = 'GRID'
  keys(2) = 'SOLUTION'

#ifdef RFLU
  keys(3) = 'GRIDSRC'
#endif   
  
  CALL ReadSection( global,IF_INPUT,nVals,keys(1:nVals),vals(1:nVals), & 
                    defined(1:nVals) )

  IF (defined(1).eqv. .true.) THEN
                                       global%gridFormat = FORMAT_ASCII
    IF (vals(1)>0.9 .AND. vals(1)<1.1) global%gridFormat = FORMAT_BINARY
    IF (vals(1) > 1.9)                 global%gridFormat = FORMAT_HDF
  ENDIF
  IF (defined(2).eqv. .true.) THEN
                                       global%solutFormat = FORMAT_ASCII
    IF (vals(2)>0.9 .AND. vals(2)<1.1) global%solutFormat = FORMAT_BINARY
    IF (vals(2) > 1.9)                 global%solutFormat = FORMAT_HDF
  ENDIF

#ifdef RFLU
  IF (defined(3).eqv. .true.) THEN
    global%gridSource = NINT(vals(3))
  ENDIF ! defined  
#endif 

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadFormatsSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadFormatsSection.F90,v $
! Revision 1.4  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2008/10/23 18:20:55  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.1  2004/12/01 16:50:18  haselbac
! Initial revision after changing case
!
! Revision 1.11  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.8  2003/05/16 02:27:43  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.7  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.6  2003/03/15 16:23:51  haselbac
! Added KIND qualifyer
!
! Revision 1.5  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/03/26 19:02:34  haselbac
! Added ROCFLU functionality
!
! Revision 1.3  2002/01/11 17:18:31  jblazek
! Updated description of I/O variables.
!
! Revision 1.2  2001/12/22 00:09:38  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1  2001/12/07 16:54:31  jblazek
! Added files to read user input.
!
!******************************************************************************







