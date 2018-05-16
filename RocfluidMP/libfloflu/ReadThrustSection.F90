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
! Purpose: read in user input related to calculation of thrust.
!
! Description: none.
!
! Input: user input file.
!
! Output: global = parameters for calculation of thrust.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ReadThrustSection.F90,v 1.4 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadThrustSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  CHARACTER(10) :: keys(7)

  LOGICAL :: defined(7)

  REAL(RFREAL) :: vals(7)

!******************************************************************************

  CALL RegisterFunction( global,'ReadThrustSection',&
  'ReadThrustSection.F90' )

! specify keywords and search for them

  keys(1) = 'TYPE'
  keys(2) = 'PLANE'
  keys(3) = 'COORD'
  keys(4) = 'WRITIME'
  keys(5) = 'WRIITER'
  keys(6) = 'OPENCLOSE'
  keys(7) = 'PAMB'
  
  CALL ReadSection( global,IF_INPUT,7,keys,vals,defined )

  IF (defined(1).eqv..true.) THEN
                                       global%thrustType = THRUST_NONE
    IF (vals(1)>0.9 .AND. vals(1)<1.1) global%thrustType = THRUST_MOM
    IF (vals(1) > 1.9)                 global%thrustType = THRUST_MOMP
  ENDIF
  IF (defined(2).eqv..true.) THEN
    IF (vals(2) < 1.1)                  global%thrustPlane = XCOORD
    IF (vals(2)>=1.1 .AND. vals(2)<2.1) global%thrustPlane = YCOORD
    IF (vals(2) >= 2.1)                 global%thrustPlane = ZCOORD
  ENDIF
  IF (defined(3).eqv..true.) THEN
    global%thrustCoord = vals(3)
  ENDIF
  IF (defined(4).eqv..true.) global%thrustSaveTime = ABS(vals(4))
  IF (defined(5).eqv..true.) THEN
    global%thrustSaveIter = INT(ABS(vals(5))+0.5_RFREAL)
    global%thrustSaveIter = MAX(1,global%thrustSaveIter)
  ENDIF
  IF (defined(6).eqv..true.) THEN
    IF (vals(6) < 0.5_RFREAL) THEN
      global%thrustOpenClose = .false.
    ELSE
      global%thrustOpenClose = .true.
    ENDIF
  ENDIF
  IF (defined(7).eqv..true.) THEN
    global%thrustPamb = ABS(vals(7))
  ENDIF

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadThrustSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadThrustSection.F90,v $
! Revision 1.4  2008/12/06 08:44:10  mtcampbe
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
! Revision 1.1  2004/12/01 16:50:52  haselbac
! Initial revision after changing case
!
! Revision 1.4  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.1  2003/06/02 17:12:00  jblazek
! Added computation of thrust.
!
!******************************************************************************







