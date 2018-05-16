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
! Purpose: mapping from TURB ID to statistics ID
!
! Description: mapping based on parameter turbStatId(:,:) input by user
!
! Input: global%turbStatId : TURB statistics ID from user input
!
! Output: global%turbStatCode : mapped TURB statistics ID
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_StatMapping.F90,v 1.5 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_StatMapping( global ) ! PUBLIC

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
#ifdef GENX
  USE ModInterfacesStatistics, ONLY : GenxStatNaming
#endif
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_global), POINTER :: global

! ... loop variables
  INTEGER  :: l, n

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
#ifdef GENX
  CHARACTER(CHRLEN), POINTER :: statName(:,:,:)
#endif
  INTEGER :: errorFlag
  INTEGER, POINTER  :: statId(:,:), statCode(:,:,:)

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_StatMapping.F90,v $'

  CALL RegisterFunction( global,'TURB_StatMapping',&
  'TURB_StatMapping.F90' )

! allocate TURB variables and set pointers ---------------------------------

  ALLOCATE( global%turbStatCode(2,2,global%turbNStat),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  statId   => global%turbStatId
  statCode => global%turbStatCode

#ifdef GENX
  ALLOCATE( global%turbStatNm(2,2,global%turbNStat),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  statName => global%turbStatNm
#endif

! set TURB mapping ---------------------------------------------------------

  DO l=1,global%turbNStat
    DO n=1,2
      IF (statId(n,l)==0) THEN 
        statCode(n,:,l) = STAT_NONE
      ELSE IF (statId(n,l)==1) THEN 
        statCode(n,1,l) = STAT_TV
        statCode(n,2,l) = TV_MIXT_MUET
#ifdef GENX
        statName(n,1,l) = 'mut'
        statName(n,2,l) = 'kg/ms'
#endif
      ELSE IF (statId(n,l)==2) THEN 
        statCode(n,1,l) = STAT_TV
        statCode(n,2,l) = TV_MIXT_TCOT
#ifdef GENX
        statName(n,1,l) = 'tcot'
        statName(n,2,l) = 'kg m/Ks^3'
#endif
      ELSE IF (statId(n,l)==3) THEN 
        statCode(n,1,l) = STAT_DV
        statCode(n,2,l) = DV_TURB_CDYN
#ifdef GENX
        statName(n,1,l) = 'Cd'
        statName(n,2,l) = ' '
#endif
      ELSE IF (statId(n,l)==4) THEN 
        statCode(n,1,l) = STAT_ST
        statCode(n,2,l) = ST_TURB_VAR1
#ifdef GENX
        statName(n,1,l) = 'alpha_1'
        statName(n,2,l) = 'Nm/s'
#endif
      ELSE IF (statId(n,l)==5) THEN 
        statCode(n,1,l) = STAT_ST
        statCode(n,2,l) = ST_TURB_VAR2
#ifdef GENX
        statName(n,1,l) = 'alpha_2+3'
        statName(n,2,l) = 'Nm/s'
#endif
      ELSE IF (statId(n,l)==6) THEN 
        statCode(n,1,l) = STAT_ST
        statCode(n,2,l) = ST_TURB_VAR3
#ifdef GENX
        statName(n,1,l) = 'alpha_4'
        statName(n,2,l) = 'Nm/s'
#endif
!      ELSE IF (statId(n,l)==4) THEN 
!        statCode(n,1,l) = STAT_SV
!        statCode(n,2,l) = E11
!#ifdef GENX
!        statName(n,1,l) = 'tau11'
!        statName(n,2,l) = 'N/m^2'
!#endif
!      ELSE IF (statId(n,l)==5) THEN 
!        statCode(n,1,l) = STAT_SV
!        statCode(n,2,l) = E22
!#ifdef GENX
!        statName(n,1,l) = 'tau22'
!        statName(n,2,l) = 'N/m^2'
!#endif
!      ELSE IF (statId(n,l)==6) THEN 
!        statCode(n,1,l) = STAT_SV
!        statCode(n,2,l) = E33
!#ifdef GENX
!        statName(n,1,l) = 'tau33'
!        statName(n,2,l) = 'N/m^2'
!#endif
      ELSE IF (statId(n,l)==7) THEN 
        statCode(n,1,l) = STAT_SV
        statCode(n,2,l) = E12
#ifdef GENX
        statName(n,1,l) = 'tau12'
        statName(n,2,l) = 'N/m^2'
#endif
      ELSE IF (statId(n,l)==8) THEN 
        statCode(n,1,l) = STAT_ST
        statCode(n,2,l) = ST_TURB_VAR4
#ifdef GENX
        statName(n,1,l) = 'MMij'
        statName(n,2,l) = 'N^2/m^4'
#endif
      ELSE IF (statId(n,l)==9) THEN 
        statCode(n,1,l) = STAT_ST
        statCode(n,2,l) = ST_TURB_VAR5
#ifdef GENX
        statName(n,1,l) = 'MLij'
        statName(n,2,l) = 'N^2/m^4'
#endif
      ELSE
        CALL ErrorStop( global,ERR_STATS_INDEXING,__LINE__, &
                        'TURB index is not defined.' ) 
      ENDIF
    ENDDO
  ENDDO

#ifdef GENX
! defined names and units of turbulence statistics

  CALL GenxStatNaming( global, FTYPE_TURB )
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_StatMapping

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_StatMapping.F90,v $
! Revision 1.5  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2004/10/22 23:33:29  wasistho
! added new turb.statistics ID, alpha1 to alpha4
!
! Revision 1.2  2004/06/07 23:10:21  wasistho
! provide Genx statistics names, units, and anytime-activation
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.3  2003/10/09 23:49:07  wasistho
! added flag ! PUBLIC
!
! Revision 1.2  2003/05/24 02:07:49  wasistho
! turbulence statistics expanded
!
! Revision 1.1  2002/11/02 02:06:29  wasistho
! Added TURB statistics
!
!
!******************************************************************************







