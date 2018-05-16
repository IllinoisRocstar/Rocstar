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
! Purpose: generate metrics specific for turbulence model selected, if needed
!
! Description: none.
!
! Input: regions = data of all grid regions.
!
! Output: regions = relevant metrics computed.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_CalcMetrics.F90,v 1.11 2009/08/26 12:28:52 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_CalcMetrics( regions, isInit ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE TURB_ModInterfaces, ONLY : TURB_CoRansWallDistOV
  USE ModError
  USE ModMPI
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  INTEGER :: isInit

! ... loop variables
  INTEGER :: iReg

! ... local variables
  TYPE(t_global), POINTER :: global
  INTEGER :: nRegions, globalWDistMethod, sndInteg, rcvInteg
  LOGICAL :: compWDist

!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'TURB_CalcMetrics',&
  'TURB_CalcMetrics.F90' )

! get general parameters ------------------------------------------------------

#ifdef RFLO
  nRegions = global%nRegions
#endif
#ifdef RFLU
  nRegions = global%nRegionsLocal
#endif

! compute RaNS/DES metrics ----------------------------------------------------

! first get total number of no-slip and injection wall faces

  IF (isInit == 1) THEN
#ifdef MPI
    sndInteg = global%turbWallDim
    CALL MPI_ALLREDUCE( sndInteg,rcvInteg,1,MPI_INTEGER, MPI_SUM, &
                        global%mpiComm, global%mpierr )
    IF (global%mpierr /=0 ) CALL ErrorStop( global,ERR_MPI_TROUBLE,__LINE__ )
    global%turbWallDim = rcvInteg
#endif
  ENDIF

! RaNS wall distances or DES length scales ------------------------------------

! global method of wall distance computation

  globalWDistMethod = 0

  DO iReg=1,nRegions
    globalWDistMethod = &
            MAX( globalWDistMethod , regions(iReg)%turbInput%wDistMethod )
  ENDDO

#ifdef RFLU
#ifdef CHARM
    sndInteg = globalWDistMethod
    CALL FEM_REDUCE( global%fieldFlagTurbInt,sndInteg,rcvInteg,FEM_MAX )
    globalWDistMethod = rcvInteg
#endif
#endif

! frequency of wall distance computation

  compWDist = .FALSE.
  IF (global%turbCalcWDistFreq == CALCWDIST_INI .AND. isInit == 1) THEN
    compWDist = .TRUE.
  ENDIF
  IF (global%turbCalcWDistFreq == CALCWDIST_FDT .AND. &
      regions(1)%irkStep == global%nrkSteps) THEN
    compWDist = .TRUE.
  ENDIF
  IF (global%turbCalcWDistFreq == CALCWDIST_SDT .AND. &
      regions(1)%irkStep == global%nrkSteps .AND. &
      (global%currentTime+global%dtMin) >= global%dTimeSystem) THEN
    compWDist = .TRUE.
  ENDIF
  IF (global%turbCalcWDistFreq == CALCWDIST_REM) THEN
    CALL ErrorStop( global,ERR_TURB_RANSINPUT,__LINE__, &
               'remesh flag is required but still missing')
  ENDIF

! compute wall distance if needed by selected model, at the right frequency

  IF ((global%turbCalcWDist .eqv. .true.) .AND. (compWDist .eqv. .true.)) THEN
    IF (globalWDistMethod == WDIST_DIRECT) THEN
      CALL TURB_CoRansWallDistOV( regions )
    ELSE
      CALL ErrorStop( global,ERR_TURB_RANSINPUT,__LINE__, &
                 'only direct method of wall distance is currently possible')
    ENDIF
  ENDIF

! compute LES and WLM metrics -------------------------------------------------
! (to be moved from TURB_InitSolution and TURB_CoViscousFluxes)

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_CalcMetrics

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_CalcMetrics.F90,v $
! Revision 1.11  2009/08/26 12:28:52  mtcampbe
! Ported to Hera.   Fixed logical expression syntax errors.  Replaced all
! IF (logical_variable)  with IF (logical_variable .eqv. .true.) as
! consistent with the specification.  Also changed: IF( ASSOCIATED(expr) )
! to IF ( ASSOCIATED(expr) .eqv. .true. ).   Intel compilers produce code
! which silently fails for some mal-formed expressions, so these changes
! are a net which should ensure that they are evaluated as intended.
!
! Revision 1.10  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.9  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.8  2005/12/29 19:46:05  wasistho
! removed CHARM stuff
!
! Revision 1.7  2005/04/15 15:07:34  haselbac
! Removed Charm/FEM stuff
!
! Revision 1.6  2004/04/20 20:45:46  wasistho
! added user option for frequency of computing wall distance
!
! Revision 1.5  2004/03/19 02:44:47  wasistho
! prepared for RFLU
!
! Revision 1.4  2004/03/13 04:28:32  wasistho
! prepared for RFLU
!
! Revision 1.3  2004/03/13 03:14:30  wasistho
! get rid of flo/flu identifier in TURB_Co.. routines
!
! Revision 1.2  2004/03/08 23:29:47  wasistho
! changed turb nomenclature
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.2  2004/02/04 22:30:30  wasistho
! move MPIsum of global%turbWallDim from allocateMemory to calcMetrics
!
! Revision 1.1  2003/10/07 02:17:02  wasistho
! initial installation of RaNS-SA and DES
!
!
!******************************************************************************







