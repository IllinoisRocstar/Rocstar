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
! Purpose: Reads in information related to the interaction Drag
!
! Description: none.
!
! Input: regions = data of all regions
!
! Output: fills user data into region%inrtInput%inrts
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_ReadDrag.F90,v 1.4 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_ReadDrag( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInteract,   ONLY : t_inrt_input,t_inrt_interact
  USE ModError
  USE ModParameters
  USE INRT_ModParameters

#ifdef RFLO
  USE ModInterfaces,      ONLY : ReadRegionSection
#endif
#ifdef RFLU
  USE ModInterfaces,      ONLY : ReadSection
#endif
  USE INRT_ModInterfaces, ONLY : INRT_SetActiveness,INRT_SetPermission, &
                                 INRT_DetermineTokens,INRT_DefineDrag
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg,iPlag

! ... local variables
  INTEGER, PARAMETER :: NKEYS_MAX = 20

  CHARACTER(CHRLEN)  :: RCSIdentString
  CHARACTER(CHRLEN)  :: keys(NKEYS_MAX)

  INTEGER :: brbeg,brend,nEdges
  INTEGER :: nPlag
  INTEGER :: nImplKeys,nNodeKeys,nKeys
  INTEGER :: ind,indMixt,indPlag0
  INTEGER :: iKeyUsed,iKeyModel,iKeyNode0
  INTEGER :: iKeyMixtActv,iKeyPlagActv,iKeyMixtPerm,iKeyPlagPerm

  LOGICAL :: defined(NKEYS_MAX)

  REAL(RFREAL) :: vals(NKEYS_MAX)

  TYPE(t_inrt_input),    POINTER :: input
  TYPE(t_inrt_interact), POINTER :: inrt
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_ReadDrag.F90,v $ $Revision: 1.4 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'INRT_ReadDrag',&
  'INRT_ReadDrag.F90' )

! begin -----------------------------------------------------------------------

! define implementation-dependent keys

  iKeyUsed  = 1
  iKeyModel = 2
  nImplKeys = 2

  keys(iKeyUsed)  = 'USED'
  keys(iKeyModel) = 'MODEL'

! define Node keys

  iKeyNode0 = nImplKeys
  iKeyMixtActv = iKeyNode0 + 1
  iKeyPlagActv = iKeyNode0 + 2
  iKeyMixtPerm = iKeyNode0 + 3
  iKeyPlagPerm = iKeyNode0 + 4
  nNodeKeys = 4

  keys(iKeyMixtActv) = 'MIXT_ACTV'
  keys(iKeyPlagActv) = 'PLAG_ACTV'
  keys(iKeyMixtPerm) = 'MIXT_PERM'
  keys(iKeyPlagPerm) = 'PLAG_PERM'

  nKeys = iKeyNode0 + nNodeKeys

  IF (nKeys > NKEYS_MAX) CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

! Read interaction section from input file

#ifdef RFLO
  CALL ReadRegionSection( global,IF_INPUT,nKeys,keys,vals,brbeg,brend,defined )
#endif
#ifdef RFLU
  CALL ReadSection( global,IF_INPUT,nKeys,keys,vals,defined )
  brbeg = LBOUND(regions,1)
  brend = UBOUND(regions,1)
#endif

  DO iReg=brbeg,brend

    input => regions(iReg)%inrtInput
    inrt  => input%inrts(INRT_TYPE_DRAG)

! - Check that INRT_DEFAULT section has been read, and that interaction has not

    IF (.NOT. input%defaultRead) &
      CALL ErrorStop( global,ERR_INRT_DEFUNREAD,__LINE__ )

    IF (inrt%used) CALL ErrorStop( global,ERR_INRT_READ,__LINE__ )

! - Use local variables for some useful quantities

    nPlag = input%nPlag

    indMixt  = input%indMixt
    indPlag0 = input%indPlag0

! - Check if interaction is used

    inrt%used = .TRUE. ! used by default when its section appears

    IF (defined(iKeyUsed)) THEN
      IF (NINT(vals(iKeyUsed)) == 0) inrt%used = .FALSE.
    END IF ! defined(iKeyUsed)

    IF (nPlag < 1) inrt%used = .FALSE. ! cannot occur without particles

    IF (.NOT. inrt%used) CYCLE ! do not bother with unused interactions

! - Define interaction (using any relevant information from input deck)

    CALL INRT_DefineDrag(regions(iReg))

! - Check for switches

    inrt%switches(INRT_SWI_DRAG_MODEL) = INRT_DRAG_MODEL_DEFAULT

    IF (defined(iKeyModel)) THEN

      SELECT CASE (NINT(vals(iKeyModel)))

      CASE (1)
        inrt%switches(INRT_SWI_DRAG_MODEL) = INRT_DRAG_MODEL_STOKES

      CASE (2)
        inrt%switches(INRT_SWI_DRAG_MODEL) = INRT_DRAG_MODEL_SN

      CASE (3)
        inrt%switches(INRT_SWI_DRAG_MODEL) = INRT_DRAG_MODEL_SMRFLD

      CASE DEFAULT
        CALL ErrorStop( global,ERR_INRT_BADSWITCH,__LINE__ )

      END SELECT ! vals(iKeyModel)

    END IF ! defined(iKeyModel)

! - Check for Mixture controls

    IF (defined(iKeyMixtActv)) &
      CALL INRT_SetActiveness(global,vals(iKeyMixtActv), &
                              inrt%activeness(indMixt))

    IF (defined(iKeyMixtPerm)) &
      CALL INRT_SetPermission(global,vals(iKeyMixtPerm), &
                              inrt%permission(indMixt))

! - Check for Lagrangian particle controls

    DO iPlag=1,nPlag+1

      ind = indPlag0 + iPlag

      IF (defined(iKeyPlagActv)) &
        CALL INRT_SetActiveness(global,vals(iKeyPlagActv), &
                                inrt%activeness(ind))

      IF (defined(iKeyPlagPerm)) &
        CALL INRT_SetPermission(global,vals(iKeyPlagPerm), &
                                inrt%permission(ind))

    END DO ! iPlag

! - Determine permission Tokens

    CALL INRT_DetermineTokens(regions(iReg),inrt)

  END DO ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_ReadDrag

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_ReadDrag.F90,v $
! Revision 1.4  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2007/03/07 22:18:05  fnajjar
! Included Sommerfeld drag law
!
! Revision 1.1  2004/12/01 21:56:35  fnajjar
! Initial revision after changing case
!
! Revision 1.7  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.6  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.5  2004/03/02 21:49:23  jferry
! Added inrtUsed flag to mixture data structure
!
! Revision 1.4  2003/04/02 22:32:04  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.3  2003/03/24 23:30:52  jferry
! overhauled rocinteract to allow interaction design to use user input
!
! Revision 1.2  2003/03/11 16:09:39  jferry
! Added comments
!
! Revision 1.1  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
!******************************************************************************







