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
! Purpose: read in user input related to Eulerian particle module.
!
! Description: none.
!
! Input: user input file.
!
! Output: regions = PEUL information
!
! Notes:
!
!   Reads in Sc and vis2, even though these do not do anything
!
!******************************************************************************
!
! $Id: PEUL_ReadConPartSection.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_ReadConPartSection( regions,nPtypes,brbeg,brend )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartEul,    ONLY : t_peul_input
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters

  USE ModInterfaces, ONLY : ReadBothRegionSection
  USE ModInterfacesInteract, ONLY : INRT_SetMaterial
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  INTEGER, INTENT(IN)     :: nPtypes
  INTEGER, INTENT(OUT)    :: brbeg,brend

! ... loop  variables
  INTEGER :: iReg,iPtype

! ... local variables
  INTEGER, PARAMETER :: NSTRKEYS_MAX = 5
  INTEGER, PARAMETER :: NKEYS_MAX = 20

  CHARACTER(CHRLEN) :: RCSIdentString
  CHARACTER(20)     :: strKeys(NSTRKEYS_MAX),keys(NKEYS_MAX)
  CHARACTER(CHRLEN) :: strVals(NSTRKEYS_MAX)

  INTEGER :: nKeys,nStrKeys,nPtypesUsed,errorFlag,readStatus
  INTEGER :: iStrKeyMaterial
  INTEGER :: iKeyUsed,iKeyDiam,iKeyPuff,iKeyInitc
  INTEGER :: iKeySchm,iKeyK2,iKeyInvK4,iKeySmooCf,iKeyConstInit
  INTEGER :: iKeyNegReport,iKeyClipModel,iKeyMethodV

  LOGICAL :: strDefined(NSTRKEYS_MAX),defined(NKEYS_MAX)

  REAL(RFREAL) :: vals(NKEYS_MAX)

  TYPE(t_peul_input), POINTER :: input
  TYPE(t_global),     POINTER :: global

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PEUL_ReadConPartSection.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PEUL_ReadConPartSection',&
  'PEUL_ReadConPartSection.F90' )

! begin -----------------------------------------------------------------------

! define string keys

  iStrKeyMaterial = 1
  nStrKeys = 1

  IF (nStrKeys > NSTRKEYS_MAX) &
    CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

  strKeys(iStrKeyMaterial) = 'MATERIAL'

! define real keys

  iKeyUsed      =  1
  iKeyDiam      =  2
  iKeyPuff      =  3
  iKeyInitc     =  4
  iKeySchm      =  5
  iKeyK2        =  6
  iKeyInvK4     =  7
  iKeySmooCf    =  8
  iKeyConstInit =  9
  iKeyNegReport = 10
  iKeyClipModel = 11
  iKeyMethodV   = 12
  nKeys         = 12

  IF (nKeys > NKEYS_MAX) CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

  keys(iKeyUsed)      = 'USED'
  keys(iKeyDiam)      = 'DIAM'
  keys(iKeyPuff)      = 'PUFF'
  keys(iKeyInitc)     = 'INITC'
  keys(iKeySchm)      = 'SCHM'
  keys(iKeyK2)        = 'K2'
  keys(iKeyInvK4)     = '1/K4'
  keys(iKeySmooCf)    = 'SMOOCF'
  keys(iKeyConstInit) = 'CONSTINIT'
  keys(iKeyNegReport) = 'NEGREPORT'
  keys(iKeyClipModel) = 'CLIPMODEL'
  keys(iKeyMethodV)   = 'METH_VEL'

! specify default values for keys

  vals                =  -1.0_RFREAL   ! default value for undefined quantity
  vals(iKeyUsed)      =   1.0_RFREAL   ! used by default CONPART section exists
  vals(iKeyPuff)      =   1.0_RFREAL   ! default for puff factor
  vals(iKeyInitc)     =   1.E-9_RFREAL ! default for initial concentration
  vals(iKeySchm)      =   1.0_RFREAL   ! default for Schmidt number
  vals(iKeyK2)        =   0.0_RFREAL   ! default for k2
  vals(iKeyInvK4)     = 128.0_RFREAL   ! default for 1/k4
  vals(iKeySmooCf)    =  -1.0_RFREAL   ! default for residual smoothing (none)
  vals(iKeyConstInit) =   0.0_RFREAL   ! default for constant init (.FALSE.)
  vals(iKeyNegReport) =   0.0_RFREAL   ! default for negativity report
  vals(iKeyClipModel) =   0.0_RFREAL   ! default for clipping model
  vals(iKeyMethodV)   =   0.0_RFREAL   ! default method (0: smoke = fluid vel)
  nPtypesUsed = MAX(nPtypes,1)         ! must be at least one particle type

! read smoke section from input file

  CALL ReadBothRegionSection( global,IF_INPUT,nKeys,nStrKeys,keys,strKeys, &
                              vals,strVals,brbeg,brend,defined,strDefined )

  IF (NINT(vals(iKeyUsed)) == 1) THEN
    readStatus = 1 ! read and used
  ELSE
    readStatus = 0 ! read and not used
  END IF ! iKeyUsed

  DO iReg = brbeg,brend

    input => regions(iReg)%peulInput

    IF (input%readStatus == -1) THEN
      input%readStatus = readStatus
    ELSE
      CALL ErrorStop( global,ERR_SEC_READ_TWICE,__LINE__ )
    ENDIF

    IF (input%readStatus == 0) CYCLE

! - allocate ptypes

    ALLOCATE( input%ptypes(nPtypesUsed),stat=errorFlag )
    global%error = errorFlag
    IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

! - fill in input data

    input%nPtypes = nPtypesUsed

    DO iPtype=1,nPtypesUsed
      NULLIFY(input%ptypes(iPtype)%material)
      IF (strDefined(iStrKeyMaterial)) THEN
        CALL INRT_SetMaterial(global,input%ptypes(iPtype)%material, &
                              strVals(iStrKeyMaterial))
      END IF ! strDefined
    END DO ! iPtype

    input%ptypes(:)%diam    = vals(iKeyDiam)
    input%ptypes(:)%puff    = vals(iKeyPuff)
    input%ptypes(:)%initc   = vals(iKeyInitc)
    input%ptypes(:)%Sc      = vals(iKeySchm)
    input%ptypes(:)%vis2    = vals(iKeyK2)
    input%ptypes(:)%vis4    = vals(iKeyInvK4) ! reciprocal taken in
                                              ! PEUL_DerivedInputValues
    input%smoocf            = vals(iKeySmooCf)
    input%constInit         = (NINT(vals(iKeyConstInit)) == 1)
    input%ptypes(:)%negReport = NINT(vals(iKeyNegReport))
    input%ptypes(:)%clipModel = NINT(vals(iKeyClipModel))

    SELECT CASE (NINT(vals(iKeyNegReport)))
    CASE (0)
      input%ptypes(:)%negReport = PEUL_NEG_REPORT_NONE
    CASE (1)
      input%ptypes(:)%negReport = PEUL_NEG_REPORT_USED
    CASE DEFAULT
      CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )
    END SELECT ! vals(iKeyNegReport)

    SELECT CASE (NINT(vals(iKeyClipModel)))
    CASE (0)
      input%ptypes(:)%clipModel = PEUL_CLIP_MODEL_NONE
    CASE (1)
      input%ptypes(:)%clipModel = PEUL_CLIP_MODEL_USED
    CASE DEFAULT
      CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )
    END SELECT ! vals(iKeyClipModel)

    SELECT CASE (NINT(vals(iKeyMethodV)))
    CASE (0)
      input%ptypes(:)%methodV = PEUL_METHV_FLUIDVEL
    CASE (1)
      input%ptypes(:)%methodV = PEUL_METHV_EQEUL
    CASE DEFAULT
      CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )
    END SELECT ! vals(iKeyMethodV)

  END DO ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_ReadConPartSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_ReadConPartSection.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:45  haselbac
! Initial revision after changing case
!
! Revision 1.11  2004/05/03 15:09:42  jferry
! added equilibrium Eulerian capability for smoke
!
! Revision 1.10  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.9  2004/03/02 21:44:52  jferry
! Added clipping options
!
! Revision 1.8  2003/09/26 21:49:06  fnajjar
! Changed interface call for INRT_SetMaterial to ModInterfacesInteract
!
! Revision 1.7  2003/04/14 16:33:52  jferry
! added option to initialize to constant for t > 0
!
! Revision 1.6  2003/04/07 18:29:01  jferry
! added inflow boundary condition and initialization to a constant
!
! Revision 1.5  2003/03/24 23:30:53  jferry
! overhauled rocinteract to allow interaction design to use user input
!
! Revision 1.4  2003/03/11 16:04:57  jferry
! Created data type for material properties
!
! Revision 1.3  2003/03/04 19:26:47  jferry
! Cleaned up routines that read sections of input files
!
! Revision 1.2  2003/02/12 23:34:48  jferry
! Replaced [io]stat=global%error with local errorFlag
!
! Revision 1.1  2003/02/11 22:52:51  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







