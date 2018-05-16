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
! Purpose: Allocate memory for all variables associated with rocperi
!          for all active regions on current processor.
!
! Description: none.
!
! Input: region = info of current region data
!
! Output: region%peri = info of current rocperi data
!
! Notes: none
!
!******************************************************************************
!
! $Id: PERI_AllocateMemory.F90,v 1.6 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_AllocateMemory( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
#ifdef RFLO
  USE ModInterfaces, ONLY : RFLO_GetDimensDummy, RFLO_GetCellOffset

#include "Indexing.h"
#endif
  USE ModPeriodic, ONLY   : t_peri
  USE ModError
  USE ModParameters
  USE PERI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), TARGET :: region

! ... loop variables
  INTEGER :: iLev

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_peri), POINTER   :: peri

  INTEGER :: nSndVar, ibc, iec, ncv, errFl
#ifdef RFLO
  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: jdc, iCOff, ijCOff
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_AllocateMemory.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'PERI_AllocateMemory',&
  'PERI_AllocateMemory.F90' )

#ifdef RFLO

! get variables and pointers --------------------------------------

! loop over all grid levels
  
  DO iLev=1,region%nGridLevels

    peri      => region%levels(iLev)%peri

! - get cell and node dimensions ----------------------------------

    CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
    ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
    iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)
#endif
#ifdef RFLU
! - get variables and pointers ----------------------------------------

! - get cell and node dimensions --------------------------------------
    ibc = 1
    iec = region%grid%nCellsTot
#endif

! - PERI derived variables

    IF (region%periInput%flowKind /= PERI_FLOW_NONE) THEN
      nSndVar = region%periInput%nVar + GAS_NVAR
      ncv     = CV_MIXT_NEQS
#ifdef RFLO
      jdc = jdcend-jdcbeg+1
      ALLOCATE( peri%varSend(jdc,nSndVar),stat=errFl ); IF (errFl>0) GOTO 88
      ALLOCATE( peri%varRecv(jdc,nSndVar),stat=errFl ); IF (errFl>0) GOTO 88
      ALLOCATE( peri%cvMean(jdc,ncv),stat=errFl ); IF (errFl>0) GOTO 88
#endif
#ifdef RFLU
      ALLOCATE( peri%varSend(iec,nSndVar),stat=errFl ); IF (errFl>0) GOTO 88
      ALLOCATE( peri%varRecv(iec,nSndVar),stat=errFl ); IF (errFl>0) GOTO 88
      ALLOCATE( peri%cvMean(iec,ncv),stat=errFl ); IF (errFl>0) GOTO 88
#endif
    ELSE
      NULLIFY( peri%varSend, peri%varRecv )
    ENDIF

    IF (region%periInput%flowKind == PERI_FLOW_CPR) THEN
#ifdef RFLO
      ALLOCATE( peri%cprVar(CPR_NCOMP,jdcbeg:jdcend),stat=errFl )
                                                     IF (errFl>0) GOTO 88
#endif
#ifdef RFLU
      ALLOCATE( peri%cprVar(CPR_NCOMP,iec),stat=errFl )
                                           IF (errFl>0) GOTO 88
#endif
    ELSE
      NULLIFY( peri%cprVar )
    ENDIF

#ifdef RFLO
  ENDDO  !iLev
#endif

  GOTO 999

! finalize ----------------------------------------------------------

88   CONTINUE

  global%error = errFl
  CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

999  CONTINUE

  CALL DeregisterFunction( global )

END SUBROUTINE PERI_AllocateMemory

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_AllocateMemory.F90,v $
! Revision 1.6  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.5  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.4  2005/03/07 05:06:42  wasistho
! install hybrid DESSA turbulence model
!
! Revision 1.3  2004/06/17 20:02:54  wasistho
! compiled with RFLU
!
! Revision 1.2  2004/06/11 21:49:40  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/06/08 23:56:56  wasistho
! changed nomenclature
!
! Revision 1.3  2003/08/29 01:40:05  wasistho
! Added TARGET attribute to region variable, since pointers are cached into it
!
! Revision 1.2  2003/05/15 02:57:05  jblazek
! Inlined index function.
!
! Revision 1.1.1.1  2003/03/29 03:36:30  wasistho
! install ROCPERI
!
!
!******************************************************************************












