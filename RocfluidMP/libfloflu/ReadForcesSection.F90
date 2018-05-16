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
! Purpose: Read in user input related to calculation of forces.
!
! Description: None.
!
! Input: 
!    global     Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReadForcesSection.F90,v 1.8 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReadForcesSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadSection
  USE ModError
  USE ModParameters
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

#ifdef RFLO
  INTEGER, PARAMETER :: NVALS_MAX = 13
#endif
#ifdef RFLU
  INTEGER, PARAMETER :: NVALS_MAX = 7  
#endif

  LOGICAL :: defined(NVALS_MAX)
  CHARACTER(10) :: keys(NVALS_MAX)
  INTEGER :: nVals   
  REAL(RFREAL) :: vals(NVALS_MAX)

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction( global,'ReadForcesSection',&
  'ReadForcesSection.F90' )

! ******************************************************************************
! Specify keywords and search for them
! ******************************************************************************

#ifdef RFLO
  nVals = NVALS_MAX

  keys(1)  = 'TYPE'
  keys(2)  = 'AEROCOEFFS'
  keys(3)  = 'REFLENGTH'
  keys(4)  = 'REFAREA'
  keys(5)  = 'REFXCOORD'
  keys(6)  = 'REFYCOORD'
  keys(7)  = 'REFZCOORD'
  keys(8)  = 'BNDBOXXMIN'
  keys(9)  = 'BNDBOXXMAX'
  keys(10) = 'BNDBOXYMIN'
  keys(11) = 'BNDBOXYMAX'
  keys(12) = 'BNDBOXZMIN'
  keys(13) = 'BNDBOXZMAX'
#endif

#ifdef RFLU
  nVals = NVALS_MAX

  keys(1) = 'FLAG'
  keys(2) = 'REFLENGTH'
  keys(3) = 'REFAREA'
  keys(4) = 'REFXCOORD'
  keys(5) = 'REFYCOORD'
  keys(6) = 'REFZCOORD'
  keys(7) = 'PATCHFLAG'    
#endif  
  
  CALL ReadSection(global,IF_INPUT,nVals,keys,vals,defined)

! ******************************************************************************
! Set variables
! ******************************************************************************

#ifdef RFLO
  IF (defined(1).eqv..true.) THEN
                                       global%forcesOn = FORCES_NONE
    IF (vals(1)>0.9 .AND. vals(1)<1.1) global%forcesOn = FORCES_PRESS
    IF (vals(1) > 1.9)                 global%forcesOn = FORCES_VISC
  ENDIF

  IF (defined(2).eqv..true.) THEN
                     global%aeroCoeffs = OFF
    IF (vals(1)>0.9) global%aeroCoeffs = ACTIVE
  ENDIF

  IF (defined(3).eqv..true.) THEN 
    global%forceRefLength = vals(3)
  END IF ! defined
  
  IF (defined(4).eqv..true.) THEN 
    global%forceRefArea   = vals(4)
  END IF ! defined
  
  IF (defined(5).eqv..true.) THEN 
    global%forceRefXCoord = vals(5)
  END IF ! defined 
  
  IF (defined(6).eqv..true.) THEN 
    global%forceRefYCoord = vals(6)
  END IF ! defined 
  
  IF (defined(7).eqv..true.) THEN 
    global%forceRefZCoord = vals(7)
  END IF ! defined  
  
  IF (defined(8).eqv..true.) THEN 
    global%acBndBoxXmin   = vals(8)
  END IF ! defined  
  
  IF (defined(9).eqv..true.) THEN 
    global%acBndBoxXmax   = vals(9)
  END IF ! defined  
  
  IF (defined(10).eqv..true.) THEN 
    global%acBndBoxYmin   = vals(10)
  END IF ! defined  
  
  IF (defined(11).eqv..true.) THEN 
    global%acBndBoxYmax   = vals(11)
  END IF ! defined  
  
  IF (defined(12).eqv..true.) THEN 
    global%acBndBoxZmin   = vals(12)
  END IF ! defined  
  
  IF (defined(13).eqv..true.) THEN 
    global%acBndBoxZmax   = vals(13)
  END IF ! defined  
#endif

#ifdef RFLU
  IF ( defined(1) .EQV. .TRUE. ) THEN 
    IF ( NINT(vals(1)) == 1 ) THEN 
      global%forceFlag = .TRUE.
    ELSE 
      global%forceFlag = .FALSE.
    END IF ! NINT
  END IF ! defined

  IF ( defined(2) .EQV. .TRUE. ) THEN 
    global%forceRefLength = vals(2)
  END IF ! defined
  
  IF ( defined(3) .EQV. .TRUE. ) THEN 
    global%forceRefArea = vals(3)
  END IF ! defined
  
  IF ( defined(4) .EQV. .TRUE. ) THEN 
    global%forceRefXCoord = vals(4)
  END IF ! defined 
  
  IF ( defined(5) .EQV. .TRUE. ) THEN 
    global%forceRefYCoord = vals(5)
  END IF ! defined 
  
  IF ( defined(6) .EQV. .TRUE. ) THEN 
    global%forceRefZCoord = vals(6)
  END IF ! defined  
  
  IF ( defined(7) .EQV. .TRUE. ) THEN
    IF ( NINT(vals(1)) == 1 ) THEN 
      global%patchCoeffFlag = .TRUE.
    ELSE 
      global%patchCoeffFlag = .FALSE.
    END IF ! NINT
  END IF ! defined        
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE ReadForcesSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadForcesSection.F90,v $
! Revision 1.8  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2008/10/23 18:20:55  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.5  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.4  2006/03/10 01:46:33  wasistho
! read acBndBox coordinates for Rocflo
!
! Revision 1.3  2006/03/09 20:47:48  wasistho
! prepared for aerodyn.coeffs calc in RFLO
!
! Revision 1.2  2005/08/09 00:52:39  haselbac
! Added reading of PATCHFLAG
!
! Revision 1.1  2004/12/01 16:50:15  haselbac
! Initial revision after changing case
!
! Revision 1.9  2004/06/16 20:00:13  haselbac
! Added RFLU code
!
! Revision 1.8  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.4  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
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
! ******************************************************************************







