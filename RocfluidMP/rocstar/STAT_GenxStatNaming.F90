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
! Purpose: provide name and units of statistics variables
!
! Description: strings are defined for the name and units of the time/ensemble
!              averaged variables based on the user input statistics
!
! Input: global%mixtStatId : mixture statistics ID from user input
!        global%turbStatId : TURB statistics ID from user input
!
! Output: global%mixtStatNm : mixture statistics names
!         global%turbStatNm : TURB statistics names
!
! Notes: none.
!
!******************************************************************************
!
! $Id: STAT_GenxStatNaming.F90,v 1.3 2008/12/06 08:44:01 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE GenxStatNaming( global, fluidType )

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global
  INTEGER :: fluidType

! ... loop variables
  INTEGER :: l

! ... local variables
  CHARACTER(CHRLEN), POINTER :: statName(:,:,:)

  INTEGER :: nStat
  INTEGER, POINTER :: statCode(:,:,:)

!******************************************************************************

  CALL RegisterFunction( global,'GenxStatNaming',&
  'STAT_GenxStatNaming.F90' )

! get dimensions and pointers -------------------------------------------------

  IF (fluidType == FTYPE_MIXT) THEN
    nStat    =  global%mixtNStat
    statCode => global%mixtStatCode
    statName => global%mixtStatNm
  ELSEIF (fluidType == FTYPE_TURB) THEN
#ifdef TURB
    nStat    =  global%turbNStat
    statCode => global%turbStatCode
    statName => global%turbStatNm
#endif
  ENDIF

! Quantities to be time-averaged is determined from the index selected by user.
! Data accumulation proceeds afterwards for each quantity.

  DO l=1,nStat

    IF ((statCode(1,1,l)==STAT_NONE).AND.(statCode(2,1,l)==STAT_NONE)) &
    GOTO 999

    IF (statCode(1,1,l)==STAT_NONE) THEN
      statName(1,1,l) = '<'//TRIM(statName(2,1,l))//'>'
      statName(1,2,l) =      TRIM(statName(2,2,l))
    ELSEIF (statCode(2,1,l)==STAT_NONE) THEN
      statName(1,1,l) = '<'//TRIM(statName(1,1,l))//'>'
    ELSE
      IF (statName(2,1,l)==statName(1,1,l)) THEN
        statName(1,1,l) = '<'//TRIM(statName(2,1,l))//'^2>'
        IF (statName(2,2,l)/=' ') THEN
          statName(1,2,l) = '('//TRIM(statName(2,2,l))//')^2'
        ENDIF
      ELSE
        IF (LEN_TRIM(statName(1,1,l))>1 .OR. LEN_TRIM(statName(2,1,l))>1) THEN
          statName(1,1,l) = '<'//TRIM(statName(1,1,l))//' '// &
                                 TRIM(statName(2,1,l))//'>'
        ELSE
          statName(1,1,l) = '<'//TRIM(statName(1,1,l))// &
                                 TRIM(statName(2,1,l))//'>'
        ENDIF
        statName(1,2,l) = TRIM(statName(1,2,l))//' '//TRIM(statName(2,2,l))
      ENDIF
    ENDIF

  ENDDO

! finalize --------------------------------------------------------------------

999 CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE GenxStatNaming

!******************************************************************************
!
! RCS Revision history:
!
! $Log: STAT_GenxStatNaming.F90,v $
! Revision 1.3  2008/12/06 08:44:01  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:15  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2006/01/03 06:34:41  wasistho
! initial import
!
! Revision 1.1  2004/12/01 21:23:42  haselbac
! Initial revision after changing case
!
! Revision 1.1  2004/06/07 23:05:56  wasistho
! provide Genx statistics names, units, and anytime-activation
!
!
!******************************************************************************







