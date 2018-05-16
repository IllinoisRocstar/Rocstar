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
! Purpose: read in user input related to time stepping.
!
! Description: none.
!
! Input: user input file.
!
! Output: global = flow type and global time-stepping parameters.
!
! Notes: 
!   1. RFLU: Must not overwrite global%currentTime if running within GENX 
!      because solver gets actual time from Roccom. In preprocessor, however,
!      do not get time from Roccom, so need to read timeStamp into separate
!      variable. This will be used in RFLU_GetUserInput to set the variable
!      global%currentTime iff RFLU_GetUserInput is called from within the
!      preprocessing module. 
!
!******************************************************************************
!
! $Id: ReadTimestepSection.F90,v 1.6 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadTimestepSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  INTEGER, PARAMETER :: NVALS_MAX = 19

  CHARACTER(10) :: keys(NVALS_MAX)

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL) :: vals(NVALS_MAX)

!******************************************************************************

  CALL RegisterFunction( global,'ReadTimestepSection',&
  'ReadTimestepSection.F90' )

! specify keywords and search for them

  keys( 1) = 'FLOWTYPE'
  keys( 2) = 'TIMESTEP'
  keys( 3) = 'MAXTIME'
  keys( 4) = 'WRITIME'
  keys( 5) = 'PRNTIME'
  keys( 6) = 'MAXITER'
  keys( 7) = 'RESTOL'
  keys( 8) = 'WRIITER'
  keys( 9) = 'PRNITER'
  keys(10) = 'STARTTIME'
  keys(11) = 'STARTITER'
  keys(12) = 'SOLVERTYPE'
  keys(13) = 'ORDER'
  keys(14) = 'MAXSUBITER'
  keys(15) = 'TOLSUBITER'
  keys(16) = 'PREDICTSOL'
  keys(17) = 'DTMINLIMIT'
  keys(18) = 'RKSCHEME'
  keys(19) = 'DTFIXED'
  
  CALL ReadSection( global,IF_INPUT,NVALS_MAX,keys,vals,defined )

  IF (defined(1).eqv..true.) THEN
                              global%flowType = FLOW_STEADY
    IF (vals(1) > 0.9_RFREAL) global%flowType = FLOW_UNSTEADY
  ENDIF
  IF (defined( 2).eqv..true.) global%dtImposed     = ABS(vals(2))
  IF (defined( 3).eqv..true.) global%maxTime       = ABS(vals(3))
  IF (defined( 4).eqv..true.) global%writeTime     = ABS(vals(4))
  IF (defined( 5).eqv..true.) global%printTime     = ABS(vals(5))
  IF (defined( 6).eqv..true.) global%maxIter       = INT(ABS(vals(6))+0.5_RFREAL)
  IF (defined( 7).eqv..true.) global%resTol        = ABS(vals(7))
  IF (defined( 8).eqv..true.) global%writeIter     = MAX(1,INT(ABS(vals(8))+0.5_RFREAL))
  IF (defined( 9).eqv..true.) global%printIter     = MAX(1,INT(ABS(vals(9))+0.5_RFREAL))
#ifndef GENX
  IF (defined(10).eqv..true.) global%timeStamp     = ABS(vals(10))
#else
  IF (defined(10).eqv..true.) global%timeStampPrep = ABS(vals(10))
#endif
  IF (defined(11).eqv..true.) global%currentIter   = INT(ABS(vals(11))+0.5_RFREAL)
#ifdef RFLO
  IF (defined(12).eqv..true.) THEN
                               global%solverType = SOLV_EXPLICIT
    IF (vals(12) > 0.9_RFREAL) global%solverType = SOLV_IMPLICIT
  ENDIF
#endif
  IF (defined(13).eqv..true.) global%tstepOrder = MAX(2,INT(ABS(vals(13))+0.5_RFREAL))
  IF (defined(14).eqv..true.) global%maxSubIter = MAX(1,INT(ABS(vals(14))+0.5_RFREAL))
  IF (defined(15).eqv..true.) global%tolSubIter = ABS(vals(15))
  IF (defined(16).eqv..true.) THEN
    IF (vals(16) < 0.9_RFREAL) THEN
      global%predictSol = .false.
    ELSE
      global%predictSol = .true.
    ENDIF
  ENDIF
  IF ( defined(18) .eqv..true.) THEN 
    global%rkScheme = INT(ABS(vals(18)) + 0.5_RFREAL)
  END IF ! defined
#ifdef RFLO
  IF (defined(19).eqv..true.) THEN
    IF (vals(19) < 0.9_RFREAL) THEN
      global%dtFixed = .false.
    ELSE
      global%dtFixed = .true.
    ENDIF
  ENDIF
#endif

#ifdef RFLU
#ifndef GENX
  IF ( defined(10) .EQV. .TRUE. ) THEN
    global%currentTime = global%timeStamp
  ELSE 
    global%currentTime = 0.0_RFREAL
  END IF ! defined
#else
  IF ( defined(10) .EQV. .FALSE. ) THEN 
    global%timeStampPrep = 0.0_RFREAL
  END IF ! defined
#endif
  IF ( defined(12) .EQV. .TRUE. ) THEN 
    IF ( NINT(vals(12)) == SOLV_EXPLICIT ) THEN
      global%solverType = SOLV_EXPLICIT
    ELSE 
      global%solverType = SOLV_IMPLICIT_NK
    END IF ! NINT
  ENDIF
  IF ( defined(17) .EQV. .TRUE. ) THEN 
    global%dtMinLimit = ABS(vals(17))
  END IF ! defined 
#endif

! finalize

  CALL DeregisterFunction( global )

END SUBROUTINE ReadTimestepSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadTimestepSection.F90,v $
! Revision 1.6  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2008/10/23 18:20:55  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.3  2006/05/09 23:36:17  wasistho
! added DTFIXED for implicit
!
! Revision 1.2  2005/08/03 18:14:07  hdewey2
! Added reading of solverType
!
! Revision 1.1  2004/12/01 16:50:55  haselbac
! Initial revision after changing case
!
! Revision 1.19  2004/11/17 16:23:21  haselbac
! Added RKSCHEME as input parameter
!
! Revision 1.18  2004/08/10 00:22:24  wasistho
! added RFREAL to real number 0.9 in IF statements
!
! Revision 1.17  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.14  2003/10/15 02:38:58  haselbac
! Added new key and parameter NVALS_MAX
!
! Revision 1.13  2003/07/03 21:48:44  jblazek
! Implemented dual-time stepping.
!
! Revision 1.12  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.11  2003/03/25 19:21:00  haselbac
! Removed old DEBUG variable
!
! Revision 1.10  2003/01/30 19:06:40  haselbac
! Added timeStampPrep variable, see note
!
! Revision 1.9  2002/10/16 21:09:59  haselbac
! Fixed bug in RFLU code segment
!
! Revision 1.8  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.7  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.6  2002/05/28 13:46:55  haselbac
! Set currentTime to timeStamp for RFLU
!
! Revision 1.5  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/02/01 00:00:24  jblazek
! Edge and corner cells defined for each level.
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







