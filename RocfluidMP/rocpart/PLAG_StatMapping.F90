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
! Purpose: mapping from PLAG ID to statistics ID
!
! Description: mapping based on parameter plagStatId(:,:) input by user
!
! Input: global%plagStatId : PLAG statistics ID from user input
!
! Output: global%plagStatCode : mapped PLAG statistics ID
!
! Notes: 
!   1. To access composition add to EV_PLAG_MASS values of 1 and 2.
!      Statistics can only performed on a droplet with two consistuents.
!
!******************************************************************************
!
! $Id: PLAG_StatMapping.F90,v 1.4 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_StatMapping( global ) ! PUBLIC

  USE ModDataTypes
  USE ModGlobal, ONLY : t_global
#ifdef GENX
  USE ModInterfacesStatistics, ONLY : GenxStatNaming
#endif
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
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

  RCSIdentString = '$RCSfile: PLAG_StatMapping.F90,v $'

  CALL RegisterFunction( global,'PLAG_StatMapping',&
  'PLAG_StatMapping.F90' )

! allocate PLAG variables and set pointers ---------------------------------

  ALLOCATE( global%plagStatCode(2,2,global%plagNStat),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  statId   => global%plagStatId
  statCode => global%plagStatCode

#ifdef GENX
  ALLOCATE( global%plagStatNm(2,2,global%plagNStat),stat=errorFlag )
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  statName => global%plagStatNm
#endif

! set PLAG mapping ---------------------------------------------------------

  DO l=1,global%plagNStat
    DO n=1,2
      IF (statId(n,l)==0) THEN 
        statCode(n,:,l) = STAT_NONE
      ELSE IF (statId(n,l)==1) THEN 
        statCode(n,1,l) = STAT_PLAGEV
        statCode(n,2,l) = EV_PLAG_DIAM3
#ifdef GENX
        statName(n,1,l) = ' '
        statName(n,2,l) = ' '
#endif
      ELSE IF (statId(n,l)==2) THEN 
        statCode(n,1,l) = STAT_PLAGEV
        statCode(n,2,l) = EV_PLAG_DIAM4
#ifdef GENX
        statName(n,1,l) = ' '
        statName(n,2,l) = ' '
#endif
      ELSE IF (statId(n,l)==3) THEN 
        statCode(n,1,l) = STAT_PLAGEV
        statCode(n,2,l) = EV_PLAG_NUMDENS
#ifdef GENX
        statName(n,1,l) = ' '
        statName(n,2,l) = ' '
#endif
      ELSE IF (statId(n,l)==4) THEN 
        statCode(n,1,l) = STAT_PLAGEV
        statCode(n,2,l) = EV_PLAG_UVEL
#ifdef GENX
        statName(n,1,l) = ' '
        statName(n,2,l) = ' '
#endif
      ELSE IF (statId(n,l)==5) THEN 
        statCode(n,1,l) = STAT_PLAGEV
        statCode(n,2,l) = EV_PLAG_VVEL
#ifdef GENX
        statName(n,1,l) = ' '
        statName(n,2,l) = ' '
#endif
      ELSE IF (statId(n,l)==6) THEN 
        statCode(n,1,l) = STAT_PLAGEV
        statCode(n,2,l) = EV_PLAG_WVEL
#ifdef GENX
        statName(n,1,l) = ' '
        statName(n,2,l) = ' '
#endif
      ELSE IF (statId(n,l)==7) THEN 
        statCode(n,1,l) = STAT_PLAGEV
        statCode(n,2,l) = EV_PLAG_MASS
#ifdef GENX
        statName(n,1,l) = ' '
        statName(n,2,l) = ' '
#endif
      ELSE IF (statId(n,l)==8) THEN 
        statCode(n,1,l) = STAT_PLAGEV
        statCode(n,2,l) = EV_PLAG_MASS+1
#ifdef GENX
        statName(n,1,l) = ' '
        statName(n,2,l) = ' '
#endif
       ELSE IF (statId(n,l)==9) THEN 
        statCode(n,1,l) = STAT_PLAGEV
        statCode(n,2,l) = EV_PLAG_MASS+2
#ifdef GENX
        statName(n,1,l) = ' '
        statName(n,2,l) = ' '
#endif
     ELSE
        CALL ErrorStop( global,ERR_STATS_INDEXING,__LINE__, &
                        'PLAG index is not defined.' ) 
      ENDIF
    ENDDO
  ENDDO

#ifdef GENX
! defined names and units of turbulence statistics

  CALL GenxStatNaming( global, FTYPE_PLAG )
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_StatMapping

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_StatMapping.F90,v $
! Revision 1.4  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2005/01/08 20:43:21  fnajjar
! Assigned properly statCode based on PLAG statistics infrastructure
!
! Revision 1.1  2004/12/29 23:30:42  wasistho
! prepared statistics for PLAG
!
!
!
!******************************************************************************







