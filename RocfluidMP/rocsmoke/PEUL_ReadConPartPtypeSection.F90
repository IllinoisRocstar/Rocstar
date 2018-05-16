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
! Purpose: read in user input related to individual Eulerian particle type
!
! Description: none.
!
! Input: user input file.
!
! Output: regions = PEUL information
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PEUL_ReadConPartPtypeSection.F90,v 1.3 2008/12/06 08:44:39 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PEUL_ReadConPartPtypeSection( regions,brbeg,brend,iPtype )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartEul,    ONLY : t_peul_input, t_peul_ptype
  USE ModError
  USE ModParameters
  USE PEUL_ModParameters

  USE ModInterfaces, ONLY : ReadBothSection
  USE ModInterfacesInteract, ONLY : INRT_SetMaterial
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)
  INTEGER, INTENT(IN)     :: brbeg,brend,iPtype

! ... loop  variables
  INTEGER :: iReg

! ... local variables
  INTEGER, PARAMETER :: NSTRKEYS_MAX = 5
  INTEGER, PARAMETER :: NKEYS_MAX = 20

  CHARACTER(CHRLEN) :: RCSIdentString
  CHARACTER(20)     :: strKeys(NSTRKEYS_MAX),keys(NKEYS_MAX)
  CHARACTER(CHRLEN) :: strVals(NSTRKEYS_MAX)

  INTEGER :: nKeys,nStrKeys,nPtypesUsed,errorFlag
  INTEGER :: iStrKeyMaterial
  INTEGER :: iKeyUsed,iKeyDiam,iKeyPuff,iKeyInitc
  INTEGER :: iKeySchm,iKeyK2,iKeyInvK4,iKeySmooCf
  INTEGER :: iKeyNegReport,iKeyClipModel,iKeyMethodV

  LOGICAL :: strDefined(NSTRKEYS_MAX),defined(NKEYS_MAX)

  REAL(RFREAL) :: vals(NKEYS_MAX)

  TYPE(t_global),     POINTER :: global
  TYPE(t_peul_ptype), POINTER :: ptype

!******************************************************************************

  RCSIdentString = &
    '$RCSfile: PEUL_ReadConPartPtypeSection.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'PEUL_ReadConPartPtypeSection',&
  'PEUL_ReadConPartPtypeSection.F90' )

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
  iKeyNegReport =  9
  iKeyClipModel = 10
  iKeyMethodV   = 11
  nKeys         = 11

  IF (nKeys > NKEYS_MAX) CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

  keys(iKeyUsed)   = 'USED'
  keys(iKeyDiam)   = 'DIAM'
  keys(iKeyPuff)   = 'PUFF'
  keys(iKeyInitc)  = 'INITC'
  keys(iKeySchm)   = 'SCHM'
  keys(iKeyK2)     = 'K2'
  keys(iKeyInvK4)  = '1/K4'
  keys(iKeySmooCf) = 'SMOOCF'
  keys(iKeyNegReport) = 'NEGREPORT'
  keys(iKeyClipModel) = 'CLIPMODEL'
  keys(iKeyMethodV)   = 'METH_VEL'

  IF (iPtype > 0) THEN

! - read smoke particle type section from input file

    CALL ReadBothSection( global,IF_INPUT,nKeys,nStrKeys,keys,strKeys, &
                          vals,strVals,defined,strDefined )

! - print warnings for keys that have no meaning

    IF (defined(iKeyUsed)) &
      WRITE(STDOUT,*) '### WARNING: key meaningful only for CONPART input, ', &
        'not CONPART_PTYPE: ', keys(iKeyUsed)

    IF (defined(iKeySmooCf)) &
      WRITE(STDOUT,*) '### WARNING: key meaningful only for CONPART input, ', &
        'not CONPART_PTYPE: ', keys(iKeySmooCf)

    DO iReg = brbeg,brend

      IF (regions(iReg)%peulInput%readStatus == -1) THEN
        CALL ErrorStop( global,ERR_PEUL_PTYPE,__LINE__ )
      ENDIF ! readStatus

      IF (regions(iReg)%peulInput%readStatus == 0) CYCLE

! --- fill in input data

      ptype => regions(iReg)%peulInput%ptypes(iPtype)

      IF (strDefined(iStrKeyMaterial)) &
        CALL INRT_SetMaterial(global,ptype%material,strVals(iStrKeyMaterial))

      IF (defined(iKeyDiam))  ptype%diam    = vals(iKeyDiam)
      IF (defined(iKeyPuff))  ptype%puff    = vals(iKeyPuff)
      IF (defined(iKeyInitc)) ptype%initc   = vals(iKeyInitc)
      IF (defined(iKeySchm))  ptype%Sc      = vals(iKeySchm)
      IF (defined(iKeyK2))    ptype%vis2    = vals(iKeyK2)
      IF (defined(iKeyInvK4)) ptype%vis4    = vals(iKeyInvK4) ! reciprocal
                                        ! taken in PEUL_DerivedInputValues
      IF (defined(iKeyNegReport)) ptype%negReport = NINT(vals(iKeyNegReport))
      IF (defined(iKeyClipModel)) ptype%clipModel = NINT(vals(iKeyClipModel))

      IF (defined(iKeyNegReport)) THEN
        SELECT CASE (NINT(vals(iKeyNegReport)))
        CASE (0)
          ptype%negReport = PEUL_NEG_REPORT_NONE
        CASE (1)
          ptype%negReport = PEUL_NEG_REPORT_USED
        CASE DEFAULT
          CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )
        END SELECT ! vals(iKeyNegReport)
      END IF ! defined(iKeyNegReport)

      IF (defined(iKeyClipModel)) THEN
        SELECT CASE (NINT(vals(iKeyClipModel)))
        CASE (0)
          ptype%clipModel = PEUL_CLIP_MODEL_NONE
        CASE (1)
          ptype%clipModel = PEUL_CLIP_MODEL_USED
        CASE DEFAULT
          CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )
        END SELECT ! vals(iKeyClipModel)
      END IF ! defined(iKeyClipModel)

      IF (defined(iKeyMethodV)) THEN
        SELECT CASE (NINT(vals(iKeyMethodV)))
        CASE (0)
          ptype%methodV = PEUL_METHV_FLUIDVEL
        CASE (1)
          ptype%methodV = PEUL_METHV_EQEUL
        CASE DEFAULT
          CALL ErrorStop( global,ERR_PEUL_BADVAL,__LINE__ )
        END SELECT ! vals(iKeyMethodV)
      END IF ! defined(iKeyMethodV)

    END DO ! iReg

  END IF ! iPtype

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PEUL_ReadConPartPtypeSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PEUL_ReadConPartPtypeSection.F90,v $
! Revision 1.3  2008/12/06 08:44:39  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:51  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:09:44  haselbac
! Initial revision after changing case
!
! Revision 1.10  2004/05/03 15:09:42  jferry
! added equilibrium Eulerian capability for smoke
!
! Revision 1.9  2004/03/05 22:09:04  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.8  2004/03/02 21:44:52  jferry
! Added clipping options
!
! Revision 1.7  2003/09/26 21:49:06  fnajjar
! Changed interface call for INRT_SetMaterial to ModInterfacesInteract
!
! Revision 1.6  2003/04/07 18:29:01  jferry
! added inflow boundary condition and initialization to a constant
!
! Revision 1.5  2003/03/24 23:30:53  jferry
! overhauled rocinteract to allow interaction design to use user input
!
! Revision 1.4  2003/03/12 17:00:35  jferry
! Added missing USE statement
!
! Revision 1.3  2003/03/11 16:04:57  jferry
! Created data type for material properties
!
! Revision 1.2  2003/03/04 19:26:47  jferry
! Cleaned up routines that read sections of input files
!
! Revision 1.1  2003/02/11 22:52:51  jferry
! Initial import of Rocsmoke
!
!******************************************************************************







