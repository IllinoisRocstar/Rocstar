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
! Purpose: Read in user input related to postprocessor.
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
!******************************************************************************
!
! $Id: ReadPostSection.F90,v 1.15 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadPostSection(global)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModInterfaces, ONLY: ReadSection
  USE ModError
  USE ModParameters
  
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  INTEGER :: nVals  

#ifdef RFLO
  INTEGER, PARAMETER :: NVALS_MAX = 9 
#endif
#ifdef RFLU
  INTEGER, PARAMETER :: NVALS_MAX = 20
#endif

  CHARACTER(CHRLEN) :: RCSIdentString
  CHARACTER(10) :: keys(NVALS_MAX)  

  LOGICAL :: defined(NVALS_MAX)

  REAL(RFREAL) :: vals(NVALS_MAX)

!******************************************************************************

  RCSIdentString = '$RCSfile: ReadPostSection.F90,v $ $Revision: 1.15 $'

  CALL RegisterFunction(global,'ReadPostSection',&
  'ReadPostSection.F90')

! specify keywords and search for them

  nVals = NVALS_MAX

#ifdef RFLO
  keys( 1) = 'PLTTYPE'
  keys( 2) = 'TIME'
  keys( 3) = 'ITER'
  keys( 4) = 'OUTFORMAT'
  keys( 5) = 'STATSFLAG'
  keys( 6) = 'TURBFLAG'
  keys( 7) = 'PLAGFLAG'
  keys( 8) = 'RADIFLAG'
  keys( 9) = 'SPECFLAG'

  CALL ReadSection( global,IF_INPUT,nVals,keys(1:nVals),vals(1:nVals), &
                    defined(1:nVals) )

  global%postStatsFlag = .FALSE.
  global%postTurbFlag  = .FALSE.
  global%postPlagFlag  = .FALSE.
  global%postRadiFlag  = .FALSE.
  global%postSpecFlag  = .FALSE.

  IF (defined(1).eqv..true.) THEN
    IF (INT(vals(1)+0.5_RFREAL) == 1) THEN 
      global%postPlotType = PLOT_GRID_ONLY
    ELSE 
      global%postPlotType = PLOT_GRID_FLOW
    END IF
  END IF

  IF (defined(2).eqv..true.) THEN
    global%postTime   = ABS(vals(2))
  END IF

  IF (defined(3).eqv..true.) THEN
    global%postIter   = INT(vals(3)+0.5_RFREAL)
  END IF

  IF (defined(4).eqv..true.) THEN
    global%postOutFmt = INT(vals(4)+0.5_RFREAL)
  END IF
  
  IF (defined(5).eqv..true.) THEN
    IF (INT(vals(5)+0.5_RFREAL) == 1) THEN 
      global%postStatsFlag = .TRUE.
    END IF
  END IF
  
  IF (defined(6).eqv..true.) THEN
    IF (INT(vals(6)+0.5_RFREAL) == 1) THEN 
      global%postTurbFlag = .TRUE.
    END IF
  END IF
  
  IF (defined(7).eqv..true.) THEN
    IF (INT(vals(7)+0.5_RFREAL) == 1) THEN 
      global%postPlagFlag = .TRUE.
    END IF
  END IF
  
  IF (defined(8).eqv..true.) THEN
    IF (INT(vals(8)+0.5_RFREAL) == 1) THEN 
      global%postRadiFlag = .TRUE.
    END IF
  END IF
  
  IF (defined(9).eqv..true.) THEN
    IF (INT(vals(9)+0.5_RFREAL) == 1) THEN 
      global%postSpecFlag = .TRUE.
    END IF
  END IF

  IF (global%postOutFmt <= 1) THEN
    global%postOutFmt = PLOT_FMT_GENERIC
  ELSE IF (global%postOutFmt == 2) THEN
    global%postOutFmt = PLOT_FMT_TECPLOT
  ELSE
    global%postOutFmt = PLOT_FMT_TECASCII
  ENDIF
#endif

#ifdef RFLU
  keys( 1) = 'PLTTYPE'
  keys( 2) = 'MERGEFLAG'
  keys( 3) = 'PLTVOLFLAG'
  keys( 4) = 'INTERTYPE'
  keys( 5) = 'INTERORDER'
  keys( 6) = 'SPECFLAG'
  keys( 7) = 'EXTRFLAG'    
  keys( 8) = 'DISCFLAG'
  keys( 9) = 'SCHTYPE'
  keys(10) = 'SCHEXP'
  keys(11) = 'NFRINGES'
  keys(12) = 'VORTFLAG'
  keys(13) = 'COREFLAG'
  keys(14) = 'WRIMERGE'
  keys(15) = 'ERRFLAG'
  keys(16) = 'OUTFORMAT'
  keys(17) = 'NSERVERS'
  keys(18) = 'PLTPATFLAG'
  keys(19) = 'PEULFLAG'
  keys(20) = 'GRADFLAG'

  CALL ReadSection(global,IF_INPUT,nVals,keys,vals,defined)

  IF ( defined(1) .EQV. .TRUE. ) THEN
    IF ( NINT(vals(1)) == 1 ) THEN 
      global%postPlotType = PLOT_GRID_ONLY
    ELSE 
      global%postPlotType = PLOT_GRID_FLOW
    END IF ! NINT(vals)
  END IF ! defined

  IF ( defined(2) .EQV. .TRUE. ) THEN
    IF ( NINT(vals(2)) == 1 ) THEN 
      global%postMergeFlag = .TRUE.
    ELSE 
      global%postMergeFlag = .FALSE.
    END IF ! NINT(vals)
  END IF ! defined
  
  IF ( defined(3) .EQV. .TRUE. ) THEN
    IF ( NINT(vals(3)) == 1 ) THEN 
      global%postPlotVolFlag = .TRUE.
    ELSE 
      global%postPlotVolFlag = .FALSE.
    END IF ! NINT(vals)
  END IF ! defined 
  
  IF ( defined(4) .EQV. .TRUE. ) THEN 
    global%postInterpType = NINT(vals(4))
  END IF ! defined(4)
  
  IF ( defined(5) .EQV. .TRUE. ) THEN 
    global%postInterpOrder = NINT(vals(5))
  END IF ! defined
  
  IF ( defined(6) .EQV. .TRUE. ) THEN 
    IF ( NINT(vals(6)) == 1 ) THEN 
      global%postSpecFlag = .TRUE. 
    ELSE 
      global%postSpecFlag = .FALSE.
    END IF ! NINT(vals)
  END IF ! defined 
  
  IF ( defined(7) .EQV. .TRUE. ) THEN 
    IF ( NINT(vals(7)) == 1 ) THEN 
      global%postExtractFlag = .TRUE. 
    ELSE 
      global%postExtractFlag = .FALSE.
    END IF ! NINT(vals)
  END IF ! defined    

  IF ( defined(8) .EQV. .TRUE. ) THEN 
    IF ( NINT(vals(8)) == 1 ) THEN 
      global%postDiscFlag = .TRUE. 
    ELSE 
      global%postDiscFlag = .FALSE.
    END IF ! NINT(vals)
  END IF ! defined    

  IF ( defined(9) .EQV. .TRUE. ) THEN 
    global%postSchType = NINT(vals(9))
  END IF ! defined

  IF ( defined(10) .EQV. .TRUE. ) THEN 
    global%postSchExp = vals(10)
  END IF ! defined

  IF ( defined(11) .EQV. .TRUE. ) THEN 
    global%postNFringes = NINT(vals(11))
  END IF ! defined
  
  IF ( defined(12) .EQV. .TRUE. ) THEN 
    IF ( NINT(vals(12)) == 1 ) THEN 
      global%postVortFlag = .TRUE. 
    ELSE 
      global%postVortFlag = .FALSE.
    END IF ! NINT(vals)
  END IF ! defined  

  IF ( defined(13) .EQV. .TRUE. ) THEN 
    IF ( NINT(vals(13)) == 1 ) THEN 
      global%postVortCoreFlag = .TRUE. 
    ELSE 
      global%postVortCoreFlag = .FALSE.
    END IF ! NINT(vals)
  END IF ! defined  

  IF ( defined(14) .EQV. .TRUE. ) THEN 
    IF ( NINT(vals(14)) == 1 ) THEN 
      global%postWriteMergeFlag = .TRUE. 
    ELSE 
      global%postWriteMergeFlag = .FALSE.
    END IF ! NINT(vals)
  END IF ! defined  
  
  IF ( defined(15) .EQV. .TRUE. ) THEN 
    IF ( NINT(vals(15)) == 1 ) THEN 
      global%postCompErrFlag = .TRUE. 
    ELSE 
      global%postCompErrFlag = .FALSE.
    END IF ! NINT(vals)
  END IF ! defined 
  
  IF ( defined(16) .EQV. .TRUE. ) THEN 
    IF ( NINT(vals(16)) == 1 ) THEN 
      global%postOutputFormat = POST_OUTPUT_FORMAT_TECPLOT 
    ELSE IF ( NINT(vals(16)) == 2 ) THEN 
      global%postOutputFormat = POST_OUTPUT_FORMAT_ENSIGHT
    END IF ! NINT(vals)
  END IF ! defined  
  
  IF ( defined(17) .EQV. .TRUE. ) THEN 
    global%postNServers = NINT(vals(17))
  END IF ! defined          

  IF ( defined(18) .EQV. .TRUE. ) THEN
    IF ( NINT(vals(18)) == 1 ) THEN 
      global%postPlotPatchFlag = .TRUE.
    ELSE 
      global%postPlotPatchFlag = .FALSE.
    END IF ! NINT(vals)
  END IF ! defined

  IF ( defined(19) .EQV. .TRUE. ) THEN 
    IF ( NINT(vals(19)) == 1 ) THEN 
      global%postLag2EulFlag = .TRUE.
    ELSE 
      global%postLag2EulFlag = .FALSE.
    END IF ! NINT(vals)
  END IF ! defined
  
  IF ( defined(20) .EQV. .TRUE. ) THEN 
    IF ( NINT(vals(20)) == 1 ) THEN 
      global%postGradFlag = .TRUE.
    ELSE 
      global%postGradFlag = .FALSE.
    END IF ! NINT(vals)
  END IF ! defined
#endif

! finalize

  CALL DeregisterFunction(global)

END SUBROUTINE ReadPostSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadPostSection.F90,v $
! Revision 1.15  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.14  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.13  2008/10/23 18:20:55  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.12  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.11  2006/01/06 22:04:37  haselbac
! Added GRADFLAG
!
! Revision 1.10  2005/12/10 23:27:33  haselbac
! Removed pick variables
!
! Revision 1.9  2005/12/10 16:53:36  haselbac
! Added PEULFLAG
!
! Revision 1.8  2005/10/28 19:16:05  haselbac
! Added PLTPATFLAG
!
! Revision 1.7  2005/10/05 20:01:06  haselbac
! Modified for ENSIGHT
!
! Revision 1.6  2005/08/10 00:29:36  haselbac
! Added reading of COREFLAG
!
! Revision 1.5  2005/08/09 00:53:01  haselbac
! Added reading of WRIMERGE and ERRFLAG
!
! Revision 1.4  2005/07/25 12:19:25  haselbac
! Added vorticity flag, removed NINT
!
! Revision 1.3  2005/07/05 19:25:12  haselbac
! Added SCHTYPE and SCHEXP variables
!
! Revision 1.2  2005/05/01 14:17:48  haselbac
! Added capability of visualizing discontinuities
!
! Revision 1.1  2004/12/01 16:50:39  haselbac
! Initial revision after changing case
!
! Revision 1.9  2004/10/27 18:04:15  haselbac
! Changed flag name for extraction, must be no longer than 10 chars
!
! Revision 1.8  2004/10/26 15:15:56  haselbac
! Added setting of postExtractFlag
!
! Revision 1.7  2004/07/24 03:43:57  wasistho
! add rocflo options in readPostSection
!
! Revision 1.6  2004/07/21 14:54:06  haselbac
! Added postInterpType
!
! Revision 1.5  2003/08/07 15:29:20  haselbac
! Added and changed var names
!
! Revision 1.4  2003/07/22 01:53:06  haselbac
! Added postInterpOrder variable
!
! Revision 1.3  2003/05/16 02:27:43  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.2  2003/05/05 18:38:41  haselbac
! Added plotType, changed pltVolFlag to plotVolFlag under global
!
! Revision 1.1  2003/04/29 14:55:48  haselbac
! Initial revision
!
!******************************************************************************







