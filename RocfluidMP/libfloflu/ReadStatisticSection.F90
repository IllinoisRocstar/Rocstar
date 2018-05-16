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
! Purpose: read in user input related to time averaging statistics
!
! Description: read-in procedure makes use of utility string reader 
!              ReadStringSection for flexible number of statistics variables
!
! Input: user input file
!
! Output: global = statistics parameters and variable indices:
!                  global%doStat    = do-statistics switch
!                  global%reStat    = statistics restart switch
!                  global%mixtNStat = number of mixture statistics variable
!                  global%turbNStat = number of TURB statistics variable
!                  global%mixtStatId= mixture statistics variable ID
!                  global%turbStatId= TURB statistics variable ID
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ReadStatisticSection.F90,v 1.9 2008/12/06 08:44:10 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE ReadStatisticSection( global )

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : ReadStringSection
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... local variables
  INTEGER            :: errorFlag,l, ival, nVals
  INTEGER, PARAMETER :: NVALS_MAX = 10
  CHARACTER(10)      :: keys(NVALS_MAX)
  LOGICAL            :: defined(NVALS_MAX)
  CHARACTER(256)     :: line(NVALS_MAX)
  INTEGER            :: vals(NVALS_MAX)

!******************************************************************************

  CALL RegisterFunction( global,'ReadStatisticSection',&
  'ReadStatisticSection.F90' )

! specify keywords and search for them --------------------------------------

  nVals  = NVALS_MAX

  keys( 1) = 'DOSTAT'
  keys( 2) = 'RESTART'
  keys( 3) = 'MIXTNSTAT'
  keys( 4) = 'MIXTSTATID'
  keys( 5) = 'TURBNSTAT'
  keys( 6) = 'TURBSTATID'
  keys( 7) = 'PLAGNSTAT'
  keys( 8) = 'PLAGSTATID'
  keys( 9) = 'PEULNSTAT'
  keys(10) = 'PEULSTATID'

  CALL ReadStringSection( global,IF_INPUT,nVals,keys,line,defined )

  IF (defined( 1).eqv..true.) READ (line(1),*) global%doStat
  IF (defined( 2).eqv..true.) READ (line(2),*) global%reStat
  IF (defined( 3).eqv..true.) READ (line(3),*) global%mixtNStat
#ifdef TURB
  IF (defined( 5).eqv..true.) READ (line(5),*) global%turbNStat
#endif
#ifdef PLAG
  IF (defined( 7).eqv..true.) READ (line(7),*) global%plagNStat
#endif
#ifdef PEUL
  IF (defined( 9).eqv..true.) READ (line(9),*) global%peulNStat
#endif

! specific mixture

  IF ((defined( 1).eqv..true.).AND.(global%doStat==ACTIVE).AND. &
      (defined( 3).eqv..true.).AND.(global%mixtNStat > 0) .AND. &
      (defined( 4).eqv..true.)) THEN
    ALLOCATE( global%mixtStatId(2,global%mixtNStat),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    READ(line(4),*) (global%mixtStatId(1,l),l=1,global%mixtNStat)
  ELSE
    NULLIFY( global%mixtStatId )
  ENDIF     ! doStat

  IF ((defined(4).eqv..true.).AND.(global%doStat==ACTIVE).AND.(global%mixtNStat > 0)) THEN
    global%mixtStatId(2,:)= MOD(global%mixtStatId(1,:),10)
    global%mixtStatId(1,:)= (global%mixtStatId(1,:)-global%mixtStatId(2,:))/10
  ENDIF

! specific turbulence

#ifdef TURB
  IF ((defined( 1).eqv..true.).AND.(global%doStat==ACTIVE).AND. &
      (defined( 5).eqv..true.).AND.(global%turbNStat > 0) .AND. &
      (defined( 6).eqv..true.)) THEN
    ALLOCATE( global%turbStatId(2,global%turbNStat),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    READ(line(6),*) (global%turbStatId(1,l),l=1,global%turbNStat)
  ELSE
    NULLIFY( global%turbStatId )
  ENDIF     ! doStat

  IF ((defined(6).eqv..true.).AND.(global%doStat==ACTIVE).AND.(global%turbNStat > 0)) THEN
    global%turbStatId(2,:)= MOD(global%turbStatId(1,:),10)
    global%turbStatId(1,:)= (global%turbStatId(1,:)-global%turbStatId(2,:))/10
  ENDIF
#endif

! specific PLAG

#ifdef PLAG
  IF ((defined( 1).eqv..true.).AND.(global%doStat==ACTIVE).AND. &
      (defined( 7).eqv..true.).AND.(global%plagNStat > 0) .AND. &
      (defined( 8).eqv..true.)) THEN
    ALLOCATE( global%plagStatId(2,global%plagNStat),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    READ(line(8),*) (global%plagStatId(1,l),l=1,global%plagNStat)
  ELSE
    NULLIFY( global%plagStatId )
  ENDIF     ! doStat

  IF ((defined(8).eqv..true.).AND.(global%doStat==ACTIVE).AND.(global%plagNStat > 0)) THEN
    global%plagStatId(2,:)= MOD(global%plagStatId(1,:),10)
    global%plagStatId(1,:)= (global%plagStatId(1,:)-global%plagStatId(2,:))/10
  ENDIF
#endif

! specific PEUL

#ifdef PEUL
  IF ((defined( 1).eqv..true.).AND.(global%doStat==ACTIVE).AND. &
      (defined( 9).eqv..true.).AND.(global%peulNStat > 0) .AND. &
      (defined(10).eqv..true.)) THEN
    ALLOCATE( global%peulStatId(2,global%peulNStat),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )
    READ(line(10),*) (global%peulStatId(1,l),l=1,global%peulNStat)
  ELSE
    NULLIFY( global%peulStatId )
  ENDIF     ! doStat

  IF ((defined(10).eqv..true.).AND.(global%doStat==ACTIVE).AND.(global%peulNStat > 0)) THEN
    global%peulStatId(2,:)= MOD(global%peulStatId(1,:),10)
    global%peulStatId(1,:)= (global%peulStatId(1,:)-global%peulStatId(2,:))/10
  ENDIF
#endif

! check input parameters ------------------------------------------------------ 

  IF (global%doStat == ACTIVE) THEN

! - general

    ival=2
    IF (.NOT. (defined(ival).eqv..true.)) THEN
      CALL ErrorStop( global,ERR_STATS_INPUT,__LINE__,keys(ival) ) 
    ENDIF
    IF ((.NOT. (defined(3).eqv..true.)) .AND. &
        (.NOT. (defined(5).eqv..true.))) THEN
      CALL ErrorStop( global,ERR_STATS_INPUT,__LINE__, &
           'mixture or module NSTAT not defined' ) 
    ENDIF

#ifndef STATS
      CALL ErrorStop( global,ERR_STATS_INPUT,__LINE__, &
           'DOSTAT=1 but executable is not compiled with STATS=1' ) 
#endif    

! - mixture

    ival=3
    IF (defined(ival).eqv..true.) THEN
      IF (global%mixtNStat > 0) THEN
        ival=4
        IF (.NOT. (defined(ival).eqv..true.)) THEN
          CALL ErrorStop( global,ERR_STATS_INPUT,__LINE__,keys(ival) ) 
        ENDIF
      ELSEIF (global%mixtNStat < 0) THEN
        CALL ErrorStop( global,ERR_STATS_INPUT,__LINE__,'mixtNSTAT < 0' ) 
      ENDIF
    ENDIF

! - turbulence

#ifdef TURB
    ival=5
    IF (defined(ival).eqv..true.) THEN
      IF (global%turbNStat > 0) THEN
        ival=6
        IF (.NOT. (defined(ival).eqv..true.)) THEN
          CALL ErrorStop( global,ERR_STATS_INPUT,__LINE__,keys(ival) ) 
        ENDIF
      ELSEIF (global%turbNStat < 0) THEN
        CALL ErrorStop( global,ERR_STATS_INPUT,__LINE__,'turbNSTAT < 0' ) 
      ENDIF
    ENDIF
#endif

! - lagrangian particles

#ifdef PLAG
    ival=7
    IF (defined(ival).eqv..true.) THEN
      IF (global%plagNStat > 0) THEN
        ival=8
        IF (.NOT. (defined(ival).eqv..true.)) THEN
          CALL ErrorStop( global,ERR_STATS_INPUT,__LINE__,keys(ival) ) 
        ENDIF
      ELSEIF (global%plagNStat < 0) THEN
        CALL ErrorStop( global,ERR_STATS_INPUT,__LINE__,'plagNSTAT < 0' ) 
      ENDIF
    ENDIF
#endif

! - eulerian particles 

#ifdef PEUL
    ival=9
    IF (defined(ival).eqv..true.) THEN
      IF (global%peulNStat > 0) THEN
        ival=10
        IF (.NOT. (defined(ival).eqv..true.)) THEN
          CALL ErrorStop( global,ERR_STATS_INPUT,__LINE__,keys(ival) ) 
        ENDIF
      ELSEIF (global%peulNStat < 0) THEN
        CALL ErrorStop( global,ERR_STATS_INPUT,__LINE__,'peulNSTAT < 0' ) 
      ENDIF
    ENDIF
#endif
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE ReadStatisticSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadStatisticSection.F90,v $
! Revision 1.9  2008/12/06 08:44:10  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2008/10/23 18:20:55  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.6  2006/01/14 23:00:22  wasistho
! added safety if read but dostat=0
!
! Revision 1.5  2005/03/07 18:45:58  wasistho
! check if stats module is compiled
!
! Revision 1.4  2005/01/08 20:34:05  fnajjar
! Bug fix of line for PLAG and PEUL
!
! Revision 1.3  2004/12/30 19:35:46  wasistho
! set NVALS_MAX to 10
!
! Revision 1.2  2004/12/29 23:28:38  wasistho
! prepared statistics for PLAG and PEUL
!
! Revision 1.1  2004/12/01 16:50:50  haselbac
! Initial revision after changing case
!
! Revision 1.12  2004/11/29 17:14:06  wasistho
! use ModInterfacesStatistics
!
! Revision 1.11  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.8  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.7  2002/11/02 01:50:25  wasistho
! Added TURB statistics
!
! Revision 1.6  2002/10/08 15:48:35  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.5  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.4  2002/07/22 22:59:10  jblazek
! Some more clean up.
!
! Revision 1.3  2002/07/22 15:47:53  wasistho
! Cleaned-up conforming Coding Rule
!
! Revision 1.2  2002/06/14 23:33:08  wasistho
! update error msg
!
! Revision 1.1  2002/06/14 21:17:01  wasistho
! Added time avg statistics
!
!
!******************************************************************************







