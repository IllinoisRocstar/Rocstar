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
! Purpose: read in user input related to reference values.
!
! Description: none.
!
! Input: user input file.
!
! Output: global = reference variables (viscosity set in DerivedInputValues).
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ReadReferenceSection.F90,v 1.4 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadReferenceSection( global )

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
  INTEGER, PARAMETER :: NVALS_MAX = 11

  CHARACTER(10) :: keys(NVALS_MAX)

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL) :: vals(NVALS_MAX)

!******************************************************************************

  CALL RegisterFunction( global,'ReadReferenceSection',&
  'ReadReferenceSection.F90' )

! specify keywords and search for them

  nVals = NVALS_MAX

  keys( 1) = 'ABSVEL'
  keys( 2) = 'PRESS'
  keys( 3) = 'DENS'
  keys( 4) = 'CP'
  keys( 5) = 'GAMMA'
  keys( 6) = 'LENGTH'
  keys( 7) = 'RENUM'
  keys( 8) = 'PRLAM'
  keys( 9) = 'PRTURB'
  keys(10) = 'SCNLAM'
  keys(11) = 'SCNTURB'
  
  CALL ReadSection( global,IF_INPUT,nVals,keys,vals,defined )

  IF (defined( 1).eqv..true.) global%refVelocity = ABS(vals( 1))
  IF (defined( 2).eqv..true.) global%refPressure = ABS(vals( 2))
  IF (defined( 3).eqv..true.) global%refDensity  = ABS(vals( 3))
  IF (defined( 4).eqv..true.) global%refCp       = ABS(vals( 4))
  IF (defined( 5).eqv..true.) global%refGamma    = ABS(vals( 5))
  IF (defined( 6).eqv..true.) global%refLength   = ABS(vals( 6))
  IF (defined( 7).eqv..true.) global%refREnum    = ABS(vals( 7))
  IF (defined( 8).eqv..true.) global%prLam       = ABS(vals( 8))
  IF (defined( 9).eqv..true.) global%prTurb      = ABS(vals( 9))
  IF (defined(10).eqv..true.) global%scnLam      = ABS(vals(10))
  IF (defined(11).eqv..true.) global%scnTurb     = ABS(vals(11))

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadReferenceSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadReferenceSection.F90,v $
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
! Revision 1.1  2004/12/01 16:50:48  haselbac
! Initial revision after changing case
!
! Revision 1.9  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.6  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.5  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/04/11 18:46:16  haselbac
! Use parameter to specify size of keys,vals,defined
!
! Revision 1.3  2002/01/11 17:18:31  jblazek
! Updated description of I/O variables.
!
! Revision 1.2  2001/12/22 00:09:38  jblazek
! Added routines to store grid and solution.
!
! Revision 1.1  2001/12/07 16:54:32  jblazek
! Added files to read user input.
!
!******************************************************************************







