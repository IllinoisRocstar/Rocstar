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
! Purpose: Allocate memory for all variables associated with rocrad
!          for all active regions on current processor.
!
! Description: none.
!
! Input: region = info of current region data
!
! Output: region%radi = info of current rocrad data
!
! Notes: none
!
!******************************************************************************
!
! $Id: RADI_AllocateMemory.F90,v 1.4 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE RADI_AllocateMemory( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModInterfaces, ONLY : RFLO_GetDimensDummyNodes, RFLO_GetNodeOffset, &
                            RFLO_GetDimensDummy, RFLO_GetCellOffset
  USE ModRadiation, ONLY  : t_radi
  USE ModError
  USE ModParameters
  USE RADI_ModParameters
  IMPLICIT NONE

#include "Indexing.h"

! ... parameters
  TYPE(t_region), TARGET :: region

! ... loop variables
  INTEGER :: iLev

! ... local variables
  TYPE(t_global), POINTER :: global
  TYPE(t_radi), POINTER   :: radi

  INTEGER :: idcbeg, idcend, jdcbeg, jdcend, kdcbeg, kdcend
  INTEGER :: idnbeg, idnend, jdnbeg, jdnend, kdnbeg, kdnend
  INTEGER :: ibc,iec, ibn,ien, iCOff,ijCOff, iNOff,ijNOff, errFl
  INTEGER :: radiModel, media, fluxLim, solMethod
  INTEGER :: nOrdin, nPol, nAzi, nAng, nCv, nDv, nGrad

!******************************************************************************

  global => region%global
  CALL RegisterFunction( global,'RADI_AllocateMemory',&
  'RADI_AllocateMemory.F90' )

#ifdef RFLO

! get variables ---------------------------------------------------

  radiModel = region%radiInput%radiModel
  media     = region%radiInput%media
  fluxLim   = region%radiInput%fluxLim
  solMethod = region%radiInput%solMethod
  nOrdin    = region%radiInput%nOrdin
  nPol      = region%radiInput%nPol
  nAzi      = region%radiInput%nAzi
  nAng      = region%radiInput%nAng

  nCv       = region%radiInput%nCv
  nDv       = region%radiInput%nDv
  nGrad     = region%radiInput%nGrad

! loop over all grid levels
  
  DO iLev=1,region%nGridLevels

    radi => region%levels(iLev)%radi

! - get cell and node dimensions ----------------------------------

    CALL RFLO_GetDimensDummy( region,iLev,idcbeg,idcend, &
                              jdcbeg,jdcend,kdcbeg,kdcend )
    CALL RFLO_GetCellOffset( region,iLev,iCOff,ijCOff )
    ibc = IndIJK(idcbeg,jdcbeg,kdcbeg,iCOff,ijCOff)
    iec = IndIJK(idcend,jdcend,kdcend,iCOff,ijCOff)

    CALL RFLO_GetDimensDummyNodes( region,iLev,idnbeg,idnend, &
                                   jdnbeg,jdnend,kdnbeg,kdnend )
    CALL RFLO_GetNodeOffset( region,iLev,iNOff,ijNOff )
    ibn = IndIJK(idnbeg,jdnbeg,kdnbeg,iNOff,ijNOff)
    ien = IndIJK(idnend,jdnend,kdnend,iNOff,ijNOff)
#endif
#ifdef RFLU
! - get variables -----------------------------------------------------

    radiModel = radiInput%radiModel
    media     = radiInput%media
    fluxLim   = radiInput%fluxLim
    solMethod = radiInput%solMethod
    nOrdin    = radiInput%nOrdin
    nPol      = radiInput%nPol
    nAzi      = radiInput%nAzi
    nAng      = radiInput%nAng

    nCv       = radiInput%nCv
    nDv       = radiInput%nDv
    nGrad     = radiInput%nGrad

    ALLOCATE( radiInput%angles(nAng,RADI_ANGLE_NCOMP),stat=errFl )
    IF (errFl>0) GOTO 88

! - get cell and node dimensions --------------------------------------
    ibc =
    iec =
    ibn =
    ien =
#endif

! - general radiation variables

    ALLOCATE( radi%qri(ibn:ien)     ,stat=errFl ); IF (errFl>0) GOTO 88
    ALLOCATE( radi%qrj(ibn:ien)     ,stat=errFl ); IF (errFl>0) GOTO 88
    ALLOCATE( radi%qrk(ibn:ien)     ,stat=errFl ); IF (errFl>0) GOTO 88
    ALLOCATE( radi%goFact(ibc:iec)  ,stat=errFl ); IF (errFl>0) GOTO 88

    ALLOCATE( radi%wvInt(XCOORD:ZCOORD,ibc:iec),stat=errFl )
    IF (errFl>0) GOTO 88
    ALLOCATE( radi%radInt(nAng        ,ibc:iec),stat=errFl )
    IF (errFl>0) GOTO 88
    ALLOCATE( radi%radCoef(ibc:iec,RADI_COEFF_NCOMP),stat=errFl )
    IF (errFl>0) GOTO 88

    IF (nCv > 0) THEN
      ALLOCATE( radi%cv(    nCv,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
      ALLOCATE( radi%cvOld( nCv,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
      ALLOCATE( radi%rhs(   nCv,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
      ALLOCATE( radi%rhsSum(nCv,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
      ALLOCATE( radi%diss(  nCv,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
    ELSE
      NULLIFY( radi%cv,radi%cvOld,radi%rhs,radi%rhsSum,radi%diss )
    ENDIF

#ifdef RFLO
    IF (nCv > 0) THEN
      ALLOCATE( radi%srad(ICOORD:KCOORD,ibc:iec),stat=errFl )
                                                 IF (errFl>0) GOTO 88
      IF ((global%flowType==FLOW_STEADY) .AND. &
          (region%radiInput%smoocf > 0._RFREAL)) THEN
        ALLOCATE( radi%epsIrs(ICOORD:KCOORD,ibc:iec), stat=errFl )
                                                      IF (errFl>0) GOTO 88
      ELSE
        NULLIFY( radi%epsIrs ) 
      ENDIF
    ELSE
      NULLIFY( radi%srad,radi%epsIrs )
    ENDIF
#endif

    IF (nDv > 0) THEN
      ALLOCATE( radi%dv(nDv,ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
    ELSE
      NULLIFY( radi%dv )
    ENDIF

    IF (nGrad > 0) THEN
      ALLOCATE( radi%gradi(nGrad,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88  
      ALLOCATE( radi%gradj(nGrad,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88  
      ALLOCATE( radi%gradk(nGrad,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
    ELSE
      NULLIFY( radi%gradi,radi%gradj,radi%gradk )
    ENDIF

    IF (radiModel == RADI_MODEL_ROSS    .OR. &
        radiModel == RADI_MODEL_FLDSRC  .OR. &
        radiModel == RADI_MODEL_FLDTRAN) THEN
      ALLOCATE( radi%fluxLim(ibc:iec),stat=errFl ); IF (errFl>0) GOTO 88
    ELSE
      NULLIFY( radi%fluxLim )
    ENDIF

    IF (radiModel == RADI_MODEL_FLDTRAN) THEN
      ALLOCATE( radi%ptens(TENSOR_SYMM_NELM,ibc:iec),stat=errFl )
      IF (errFl>0) GOTO 88
      ALLOCATE( radi%eddFact(ibc:iec),stat=errFl )
      IF (errFl>0) GOTO 88
    ELSE
      NULLIFY( radi%cv,radi%ptens,radi%eddFact )
    ENDIF

    IF (radiModel == RADI_MODEL_RTEGRAY .OR. &
        radiModel == RADI_MODEL_RTEBAND) THEN
      ALLOCATE( radi%dWghti(nAng,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
      ALLOCATE( radi%dWghtj(nAng,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
      ALLOCATE( radi%dWghtk(nAng,ibn:ien),stat=errFl ); IF (errFl>0) GOTO 88
    ELSE
      NULLIFY( radi%dWghti, radi%dWghtj, radi%dWghtk )
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

END SUBROUTINE RADI_AllocateMemory

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_AllocateMemory.F90,v $
! Revision 1.4  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/09/30 17:09:58  wasistho
! prepared for full FLD radiation model
!
! Revision 1.1  2004/09/22 02:35:49  wasistho
! changed file nomenclature from lower to upper case
!
! Revision 1.6  2004/09/22 01:31:13  wasistho
! switch LFD to FLD for flux limited diffusion
!
! Revision 1.5  2004/09/18 17:41:11  wasistho
! install Limited Flux Diffusion radiation
!
! Revision 1.4  2003/08/29 01:39:49  wasistho
! Added TARGET attribute to region variable, since pointers are cached into it
!
! Revision 1.3  2003/07/31 02:54:01  wasistho
! enter part and smoke data into radiation
!
! Revision 1.2  2003/07/18 01:38:41  wasistho
! removed bcModel from input data structure
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!
!******************************************************************************












