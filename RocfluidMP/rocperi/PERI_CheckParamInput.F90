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
! Purpose: Check PERI parameters either specified by user or set in the code.
!
! Description: The checking includes the existency and order of parameters.
!
! Input: regions = input parameters contained in periInput of all regions.
!
! Output: Error msg.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PERI_CheckParamInput.F90,v 1.5 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_CheckParamInput( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPeriodic, ONLY   : t_peri_input
  USE ModBndPatch, ONLY   : t_patch
  USE ModError
  USE ModMPI
  USE ModParameters
  USE PERI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER      :: global
  TYPE(t_peri_input), POINTER  :: input  
  TYPE(t_patch), POINTER       :: patches(:)

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_CheckParamInput.F90,v $'

  global => regions(1)%global
  CALL RegisterFunction( global,'PERI_CheckParamInput',&
  'PERI_CheckParamInput.F90' )

! check fixed parameters setting ---------------------------------------------

  IF (CPR_RHO /= 1 .OR. &
      CPR_RUC /= 2 .OR. &
      CPR_RVC /= 3 .OR. &
      CPR_UVE /= 4 .OR. &
      CPR_VVE /= 5 .OR. &
      CPR_TMP /= 6 .OR. &
      CPR_PRS /= 7) THEN
    CALL ErrorStop( global,ERR_PERI_FIXPARAM,__LINE__,'cprVar id inconsistent' )
  ENDIF   

  IF (global%nProcAlloc /= global%nRegions) THEN
    CALL ErrorStop( global,ERR_PERI_MPI,__LINE__, &
        'nProcAlloc should be = nRegions due to PERI_pgradUpdate procedure' )
  ENDIF
 
! check PERI parameter selection regionwise

#ifdef RFLO
  DO iReg=1,global%nRegions
    IF (regions(iReg)%procid==global%myProcid .AND. &   ! region active and
        regions(iReg)%active==ACTIVE) THEN              ! on my processor

      patches => regions(iReg)%levels(1)%patches(:)
#endif
#ifdef RFLU
  DO iReg = LBOUND(regions,1),UBOUND(regions,1)

      patches => regions(iReg)%patches(:)
#endif

! --- check if injection (cpr) or noslip wall (channel) bc on correct patch

      IF ((patches(1)%bcType >= BC_INJECTION  .AND. &
         patches(1)%bcType <= BC_INJECTION+BC_RANGE) .OR. &
        (patches(2)%bcType >= BC_INJECTION  .AND. &
         patches(2)%bcType <= BC_INJECTION+BC_RANGE) .OR. &
        (patches(5)%bcType >= BC_INJECTION  .AND. &
         patches(5)%bcType <= BC_INJECTION+BC_RANGE) .OR. &
        (patches(6)%bcType >= BC_INJECTION  .AND. &
         patches(6)%bcType <= BC_INJECTION+BC_RANGE)) THEN
        CALL ErrorStop( regions(iReg)%global,ERR_PERI_CPRBC,__LINE__, &
                      'cpr injection bc should be on patch 3 and 4' )
      ENDIF

      IF ((patches(1)%bcType >= BC_NOSLIPWALL  .AND. &
         patches(1)%bcType <= BC_NOSLIPWALL+BC_RANGE) .OR. &
        (patches(2)%bcType >= BC_NOSLIPWALL  .AND. &
         patches(2)%bcType <= BC_NOSLIPWALL+BC_RANGE) .OR. &
        (patches(5)%bcType >= BC_NOSLIPWALL  .AND. &
         patches(5)%bcType <= BC_NOSLIPWALL+BC_RANGE) .OR. &
        (patches(6)%bcType >= BC_NOSLIPWALL  .AND. &
         patches(6)%bcType <= BC_NOSLIPWALL+BC_RANGE)) THEN
        CALL ErrorStop( regions(iReg)%global,ERR_PERI_INPUT,__LINE__, &
                      'channel noslip wall bc should be on patch 3 and 4' )
      ENDIF

! --- peri input check
    
      input => regions(iReg)%periInput

      IF ((input%flowKind < PERI_FLOW_NONE).AND. &
          (input%flowKind > PERI_FLOW_CHANNEL)) THEN
        CALL ErrorStop( global,ERR_PERI_INPUT,__LINE__, &
        'flowkind out of range' )
      ENDIF

      IF (input%flowKind == PERI_FLOW_BOLA) THEN
        CALL ErrorStop( global,ERR_PERI_INPUT,__LINE__, &
        'boundary layer flow is not ready yet' )
      ENDIF

#ifdef RFLO
      IF (input%split(JCOORD) /= OFF) THEN
        IF (input%split(ICOORD)/=OFF .OR. input%split(KCOORD)/=OFF) THEN
          CALL ErrorStop( global,ERR_PERI_INPUT,__LINE__, &
          'cannot combine splitting domain in wall normal with other directions' )
        ENDIF
      ENDIF
#endif

      IF (input%flowKind == PERI_FLOW_CPR) THEN
        IF (input%nVar /= CPR_NVAR) THEN
          CALL ErrorStop( global,ERR_PERI_INPUT,__LINE__,'cpr nVar /= 7' )
        ENDIF
        IF (input%minjRate /= &
          patches(3)%mixt%vals(BCDAT_INJECT_MFRATE,0)) THEN
          CALL ErrorStop( global,ERR_PERI_CPRBC,__LINE__, &
          'mass injection rate in inj. patches different' )
        ENDIF
        IF (input%bulkmFlux < 0._RFREAL) THEN
          CALL ErrorStop( global,ERR_PERI_INPUT,__LINE__,'cpr bulk mflux < 0' )
        ENDIF
        IF (input%cprEpsilon < 0._RFREAL) THEN
          CALL ErrorStop( global,ERR_PERI_INPUT,__LINE__,'cpr epsilon < 0' )
        ENDIF
        IF (input%headPres < 0._RFREAL) THEN
          CALL ErrorStop( global,ERR_PERI_INPUT,__LINE__,'cpr head pres. < 0' )
        ENDIF
        IF (input%headTemp < 0._RFREAL) THEN
          CALL ErrorStop( global,ERR_PERI_INPUT,__LINE__,'cpr head temp. < 0' )
        ENDIF
      ELSEIF (input%flowKind == PERI_FLOW_CHANNEL) THEN
        IF ((input%pgradType < 0) .OR. (input%pgradType > 1)) THEN
          CALL ErrorStop( global,ERR_PERI_INPUT,__LINE__, &
                          'PGRADTYPE should be 0 or 1' )
        ENDIF
      ENDIF

#ifdef RFLO
    ENDIF ! region active
#endif
  ENDDO   ! iReg

! finalize ---------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PERI_CheckParamInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_CheckParamInput.F90,v $
! Revision 1.5  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/08/19 15:40:08  mparmar
! Renamed patch variables
!
! Revision 1.2  2004/06/11 21:49:40  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/06/08 23:56:56  wasistho
! changed nomenclature
!
! Revision 1.4  2003/10/17 20:25:39  wasistho
! generalize check of nregions = nprocs
!
! Revision 1.3  2003/09/18 01:56:17  wasistho
! added ijksplit and pgradType in PERI_PgradUpdate
!
! Revision 1.2  2003/04/03 00:31:07  wasistho
! enable channel flow
!
! Revision 1.1.1.1  2003/03/29 03:36:30  wasistho
! install ROCPERI
!
!
!
!******************************************************************************







