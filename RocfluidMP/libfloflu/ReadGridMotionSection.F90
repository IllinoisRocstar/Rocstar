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
! Purpose: Read in user input related to the grid motion scheme.
!
! Description: None.
!
! Input: from file.
!
! Output:
!   global  = global parameters related to grid motion (RFLO)
!   regions = region data (RFLU).
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReadGridMotionSection.F90,v 1.21 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2003-2004 by the University of Illinois
!
! ******************************************************************************

#ifdef RFLO
SUBROUTINE ReadGridMotionSection( global )
#endif
#ifdef RFLU
SUBROUTINE ReadGridMotionSection(regions)
#endif

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModDataStruct, ONLY: t_region
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

#ifdef RFLO
  TYPE(t_global), POINTER :: global
#endif
#ifdef RFLU
  TYPE(t_region), POINTER :: regions(:)
#endif

! ==============================================================================
! Locals
! ==============================================================================

#ifdef RFLU
  INTEGER :: iReg
#endif
  INTEGER, PARAMETER :: NVALS_MAX = 16
  INTEGER :: nVals

  CHARACTER(10) :: keys(NVALS_MAX)

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL) :: vals(NVALS_MAX)

! ******************************************************************************
! Start, specify keywords and search for them
! ******************************************************************************

  keys(1) = 'TYPE'
  keys(2) = 'NITER'
#ifdef RFLO
  keys(3) = 'VITER'
  keys(4) = 'SITER'
  keys(5) = 'WEIGHT'
  keys(6) = 'AMPLIFX'
  keys(7) = 'AMPLIFY'
  keys(8) = 'AMPLIFZ'
  keys(9) = 'POWER'
  keys(10)= 'NEIGHBOR'
  keys(11)= 'NSURFMATCH'
  keys(12)= 'ORTHODIR'
  keys(13)= 'ORTHOWGHTX'
  keys(14)= 'ORTHOWGHTY'
  keys(15)= 'ORTHOWGHTZ'
  keys(16)= 'ORTHOCELL'

  nVals = NVALS_MAX

  CALL RegisterFunction( global,'ReadGridMotionSection',&
  'ReadGridMotionSection.F90' )

  CALL ReadSection( global,IF_INPUT,nVals,keys,vals,defined )

  IF (defined(1).eqv..true.) THEN
                       global%moveGridScheme = MOVEGRID_BLOCKS
    IF (vals(1) > 0.9 .AND. vals(1) < 1.9) &
                       global%moveGridScheme = MOVEGRID_GLOBAL
    IF (vals(1) > 1.9 .AND. vals(1) < 2.9) &
                       global%moveGridScheme = MOVEGRID_FRAME
    IF (vals(1) > 2.9 .AND. vals(1) < 3.9) &
                       global%moveGridScheme = MOVEGRID_FOMS
    IF (vals(1) > 3.9 .AND. vals(1) < 4.9) &
                       global%moveGridScheme = MOVEGRID_ELGLOBAL
    IF (vals(1) > 4.9 .AND. vals(1) < 5.9) &
                       global%moveGridScheme = MOVEGRID_ELFRAME
    IF (vals(1) > 5.9) &
                       global%moveGridScheme = MOVEGRID_VMS
  ENDIF
  IF (defined(2).eqv..true.) THEN
    global%moveGridNiter   = INT(ABS(vals(2))+0.5_RFREAL)
  ENDIF
  IF (defined(3).eqv..true.) THEN
    global%moveGridViter   = INT(ABS(vals(3))+0.5_RFREAL)
  ENDIF
  IF (defined(4).eqv..true.) THEN
    global%moveGridSiter   = INT(ABS(vals(4))+0.5_RFREAL)
  ENDIF
  IF (defined(5).eqv..true.) THEN
    global%moveGridWeight  = vals(5)
  ENDIF
  IF (defined(6).eqv..true.) THEN
    global%moveGridAmplifX = vals(6)
  ENDIF
  IF (defined(7).eqv..true.) THEN
    global%moveGridAmplifY = vals(7)
  ENDIF
  IF (defined(8).eqv..true.) THEN
    global%moveGridAmplifZ = vals(8)
  ENDIF
  IF (defined(9).eqv..true.) THEN
    global%moveGridPower   = vals(9)
  ENDIF
  IF (defined(10).eqv..true.) THEN
    global%moveGridNbour   = INT(ABS(vals(10))+0.5_RFREAL)
  ENDIF
  IF (defined(11).eqv..true.) THEN
    global%moveGridNsmatch = INT(ABS(vals(11))+0.5_RFREAL)
  ENDIF
  IF (defined(12).eqv..true.) THEN
    global%moveGridOrthDir = INT(ABS(vals(12))+0.5_RFREAL)
  ENDIF
  IF (defined(13).eqv..true.) THEN
    global%moveGridOrthWghtX = vals(13)
  ENDIF
  IF (defined(14).eqv..true.) THEN
    global%moveGridOrthWghtY = vals(14)
  ENDIF
  IF (defined(15).eqv..true.) THEN
    global%moveGridOrthWghtZ = vals(15)
  ENDIF
  IF (defined(16).eqv..true.) THEN
    global%moveGridOrthCell  = vals(16)
  ENDIF

  CALL DeregisterFunction( global )
#endif

#ifdef RFLU
  keys(3) = 'SFACT'

  nVals = 3

  CALL RegisterFunction(regions(1)%global,'ReadGridMotionSection',&
  'ReadGridMotionSection.F90')

  CALL ReadSection(regions(1)%global,IF_INPUT,nVals,keys,vals,defined) 
  
  IF ( defined(1) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%moveGridType = NINT(vals(1))
    END DO ! iReg    
  END IF ! defined  
  
  IF ( defined(2) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)
      regions(iReg)%mixtInput%moveGridNIter = NINT(vals(2))
    END DO ! iReg
  END IF ! defined
  
  IF ( defined(3) .EQV. .TRUE. ) THEN
    DO iReg = LBOUND(regions,1),UBOUND(regions,1)      
      regions(iReg)%mixtInput%moveGridSFact = vals(3)
    END DO ! iReg
  END IF ! defined 

  CALL DeregisterFunction( regions(1)%global )
#endif 

! ******************************************************************************
! End
! ******************************************************************************

END SUBROUTINE ReadGridMotionSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadGridMotionSection.F90,v $
! Revision 1.21  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.20  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.19  2008/10/23 18:20:55  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.18  2006/03/18 13:28:03  wasistho
! added orthDir and orthWghtX,Y,Z
!
! Revision 1.17  2006/03/08 06:38:41  wasistho
! added moveGridSiter and Viter
!
! Revision 1.16  2006/03/02 03:52:51  wasistho
! bug fixed, defined(11) was defined(10)
!
! Revision 1.15  2006/03/02 01:26:12  wasistho
! split movegrid_epde to elglobal and elframe
!
! Revision 1.14  2006/02/08 04:02:35  wasistho
! added movegrid_epde
!
! Revision 1.13  2005/10/28 22:46:34  wasistho
! added orthocell
!
! Revision 1.12  2005/10/28 05:41:48  wasistho
! read FOMS scheme
!
! Revision 1.11  2005/08/28 23:47:55  wasistho
! added orthoWght for block orthogonality of RFLO global-gridmotion
!
! Revision 1.10  2005/08/18 19:47:46  wasistho
! added moveGridNsmatch
!
! Revision 1.9  2005/06/23 05:50:21  wasistho
! changed NEIGHBOUR to NEIGHBOR
!
! Revision 1.8  2005/06/23 03:32:57  wasistho
! fixed bug in reading NEIGHBOUR
!
! Revision 1.7  2005/06/23 01:35:49  wasistho
! added input parameter NEIGHBOUR for rocflo
!
! Revision 1.6  2005/06/04 01:02:40  wasistho
! distinguished to AMPLIFX,Y,Z
!
! Revision 1.5  2005/06/02 22:57:53  wasistho
! added moveGridAmplif and moveGridPower
!
! Revision 1.4  2005/06/02 03:21:54  wasistho
! shuffle MoveGridVms with MoveGridFrame
!
! Revision 1.3  2005/05/28 21:24:16  wasistho
! added moveGridFrame
!
! Revision 1.2  2005/05/21 04:16:21  wasistho
! added MOVEGRID_VMS in Rocflo part
!
! Revision 1.1  2004/12/01 16:50:23  haselbac
! Initial revision after changing case
!
! Revision 1.13  2004/10/19 19:25:41  haselbac
! Cosmetics only
!
! Revision 1.12  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.9  2003/08/27 21:33:58  haselbac
! Changed logic so vars not overwritten if not present
!
! Revision 1.8  2003/08/25 21:51:24  jblazek
! Full version of global grid motion scheme.
!
! Revision 1.7  2003/08/11 21:51:17  jblazek
! Added basic global grid smoothing scheme.
!
! Revision 1.6  2003/05/16 02:27:43  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.5  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.4  2003/03/31 16:31:56  haselbac
! Added reading of new argument
!
! Revision 1.3  2003/03/15 16:26:10  haselbac
! Added KIND qualifyer
!
! Revision 1.2  2003/02/07 23:10:43  haselbac
! Fixed stupid mistake in setting of sFact
!
! Revision 1.1  2003/01/28 16:13:38  haselbac
! Initial revision
!
! ******************************************************************************








