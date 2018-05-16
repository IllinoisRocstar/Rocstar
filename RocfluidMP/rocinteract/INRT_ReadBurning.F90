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
! Purpose: Reads in information related to the interaction Burning
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
! $Id: INRT_ReadBurning.F90,v 1.3 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_ReadBurning( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModMaterials,  ONLY : t_material
  USE ModInteract,   ONLY : t_inrt_input,t_inrt_interact
  USE ModError
  USE ModParameters
  USE INRT_ModParameters

#ifdef RFLO
  USE ModInterfaces,         ONLY : ReadBothRegionSection
#endif
#ifdef RFLU
  USE ModInterfaces,         ONLY : ReadBothSection
#endif
  USE ModInterfaces,         ONLY : MakeNumberedKeys

  USE ModInterfacesInteract, ONLY : INRT_SetMaterial
  USE INRT_ModInterfaces,    ONLY : INRT_SetActiveness,INRT_SetPermission, &
                                    INRT_DetermineTokens,INRT_DefineBurning
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg,iPlag,iPeul,iPeulOutEdge

! ... local variables
  INTEGER, PARAMETER :: NSTRKEYS_MAX = 5
  INTEGER, PARAMETER :: NKEYS_MAX    = 50
  INTEGER, PARAMETER :: NPEUL_MAX    = 10

  CHARACTER(CHRLEN)  :: RCSIdentString
  CHARACTER(CHRLEN)  :: strKeys(NSTRKEYS_MAX),keys(NKEYS_MAX)
  CHARACTER(CHRLEN)  :: strVals(NSTRKEYS_MAX)

  INTEGER :: brbeg,brend
  INTEGER :: nPlag,nPeul
  INTEGER :: iEdge,nEdges,nPeulOutEdges,nPeulOxEdges
  INTEGER :: nStrKeys,nFixedImplKeys,nFixedNodeKeys,nKeys
  INTEGER :: ind,indMixt,indPlag0,indPeul0,indIntl
  INTEGER :: iStrKeyMaterialIn,iStrKeyMaterialOut,iStrKeyMaterialOx
  INTEGER :: iKey,iKeyUsed,iKeyModel,iKeyOxUsed,iKeyHeatCoef
  INTEGER :: iKeyVaporMeth,iKeyVaporTemp
  INTEGER :: iKeyMfrcPlag,iKeyMfrcPeul0
  INTEGER :: iKeyNode0,iKeyMixtActv,iKeyMixtPerm,iKeyPlagActv,iKeyPlagPerm
  INTEGER :: iKeyPeulActv0,iKeyPeulPerm0,iKeyActv,iKeyPerm

  LOGICAL :: plagOutExists,oxUsed
  LOGICAL :: defined(NKEYS_MAX),strDefined(NSTRKEYS_MAX)

  REAL(RFREAL) :: coef
  REAL(RFREAL) :: vals(NKEYS_MAX)

  TYPE(t_material),      POINTER :: matIn,matOut,matOx
  TYPE(t_inrt_input),    POINTER :: input
  TYPE(t_inrt_interact), POINTER :: inrt
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_ReadBurning.F90,v $ $Revision: 1.3 $'

  global => regions(1)%global

  CALL RegisterFunction( global,'INRT_ReadBurning',&
  'INRT_ReadBurning.F90' )

! begin -----------------------------------------------------------------------

! define string keys

  iStrKeyMaterialIn  = 1
  iStrKeyMaterialOut = 2
  iStrKeyMaterialOx  = 3
  nStrKeys = 3

  strKeys(iStrKeyMaterialIn)  = 'MATERIAL_IN'
  strKeys(iStrKeyMaterialOut) = 'MATERIAL_OUT'
  strKeys(iStrKeyMaterialOx)  = 'MATERIAL_OX'

  IF (nStrKeys > NSTRKEYS_MAX) &
    CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

! define implementation-dependent keys

  iKeyUsed      = 1
  iKeyModel     = 2
  iKeyOxUsed    = 3
  iKeyVaporMeth = 4
  iKeyVaporTemp = 5
  iKeyHeatCoef  = 6
  iKeyMfrcPlag  = 7
  nFixedImplKeys = 8

  keys(iKeyUsed)      = 'USED'
  keys(iKeyModel)     = 'MODEL'
  keys(iKeyOxUsed)    = 'OX_USED'
  keys(iKeyVaporMeth) = 'VAPOR_METH'
  keys(iKeyVaporTemp) = 'VAPOR_TEMP'
  keys(iKeyHeatCoef)  = 'HEAT_COEF'
  keys(iKeyMfrcPlag)  = 'MFRC_PLAG'

  iKeyMfrcPeul0 = nFixedImplKeys

#ifdef RFLO
  CALL MakeNumberedKeys(keys,iKeyMfrcPeul0+1,'MFRC_PEUL',1,NPEUL_MAX,1)
#endif
#ifdef RFLU
  CALL MakeNumberedKeys(keys,iKeyMfrcPeul0+1,'MFRC_SPEC',1,NPEUL_MAX,1)
#endif

! define Node keys

  iKeyNode0 = iKeyMfrcPeul0 + NPEUL_MAX
  iKeyMixtActv = iKeyNode0 + 1
  iKeyMixtPerm = iKeyNode0 + 2
  iKeyPlagActv = iKeyNode0 + 3
  iKeyPlagPerm = iKeyNode0 + 4
  nFixedNodeKeys = 4

  keys(iKeyMixtActv) = 'MIXT_ACTV'
  keys(iKeyMixtPerm) = 'MIXT_PERM'
  keys(iKeyPlagActv) = 'PLAG_ACTV'
  keys(iKeyPlagPerm) = 'PLAG_PERM'

  iKeyPeulActv0 = iKeyNode0 + nFixedNodeKeys
  iKeyPeulPerm0 = iKeyPeulActv0 + NPEUL_MAX

#ifdef RFLO
  CALL MakeNumberedKeys(keys,iKeyPeulActv0+1,'PEUL',1,NPEUL_MAX,1)
  CALL MakeNumberedKeys(keys,iKeyPeulPerm0+1,'PEUL',1,NPEUL_MAX,1)
#endif
#ifdef RFLU
  CALL MakeNumberedKeys(keys,iKeyPeulActv0+1,'SPEC',1,NPEUL_MAX,1)
  CALL MakeNumberedKeys(keys,iKeyPeulPerm0+1,'SPEC',1,NPEUL_MAX,1)
#endif

  DO iPeul=1,NPEUL_MAX
    keys(iKeyPeulActv0+iPeul) = TRIM(keys(iKeyPeulActv0+iPeul))//'_ACTV'
    keys(iKeyPeulPerm0+iPeul) = TRIM(keys(iKeyPeulPerm0+iPeul))//'_PERM'
  END DO ! iPeul

  nKeys = iKeyPeulPerm0 + NPEUL_MAX

  IF (nKeys > NKEYS_MAX) CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

! Read interaction section from input file

#ifdef RFLO
  CALL ReadBothRegionSection( global,IF_INPUT,nKeys,nStrKeys,keys,strKeys, &
                              vals,strVals,brbeg,brend,defined,strDefined )
#endif
#ifdef RFLU
  CALL ReadBothSection( global,IF_INPUT,nKeys,nStrKeys,keys,strKeys, &
                        vals,strVals,defined,strDefined )
  brbeg = LBOUND(regions,1)
  brend = UBOUND(regions,1)
#endif

  DO iReg=brbeg,brend

    input => regions(iReg)%inrtInput
    inrt  => input%inrts(INRT_TYPE_BURNING)

! - Check that INRT_DEFAULT section has been read, and that interaction has not

    IF (.NOT. input%defaultRead) &
      CALL ErrorStop( global,ERR_INRT_DEFUNREAD,__LINE__ )

    IF (inrt%used) CALL ErrorStop( global,ERR_INRT_READ,__LINE__ )

! - Use local variables for some useful quantities

    nPlag = input%nPlag
    nPeul = input%nPeul

    IF (nPeul > NPEUL_MAX) &
      CALL ErrorStop( global,ERR_EXCEEDS_DECL_MEM,__LINE__ )

    indMixt  = input%indMixt
    indPlag0 = input%indPlag0
    indPeul0 = input%indPeul0
    indIntl  = input%indIntl

! - Check if interaction is used

    inrt%used = .TRUE. ! used by default when its section appears

    IF (defined(iKeyUsed)) THEN
      IF (NINT(vals(iKeyUsed)) == 0) inrt%used = .FALSE.
    END IF ! defined(iKeyUsed)

    IF (nPlag < 1) inrt%used = .FALSE. ! cannot occur without particles

    IF (.NOT. inrt%used) CYCLE ! do not bother with unused interactions

! - Check if oxidizer is used

    oxUsed = .FALSE.

    IF (defined(iKeyOxUsed) .AND. nPeul > 0) THEN

      SELECT CASE (NINT(vals(iKeyOxUsed)))

      CASE (0)
        oxUsed = .FALSE.

      CASE (1)
        oxUsed = .TRUE.

      CASE DEFAULT
        CALL ErrorStop( global,ERR_INRT_BADSWITCH,__LINE__ )

      END SELECT ! vals(iKeyOxUsed)

    END IF ! defined(iKeyOxUsed)

! - Set pointers to in, out, and ox materials

    IF (strDefined(iStrKeyMaterialIn)) THEN
      CALL INRT_SetMaterial(global,matIn,strVals(iStrKeyMaterialIn))
    ELSE
      CALL ErrorStop( global,ERR_INRT_MISSINGMAT,__LINE__ )
    END IF ! strDefined(iStrKeyMaterialIn)

    IF (strDefined(iStrKeyMaterialOut)) THEN
      CALL INRT_SetMaterial(global,matOut,strVals(iStrKeyMaterialOut))
    ELSE
      CALL ErrorStop( global,ERR_INRT_MISSINGMAT,__LINE__ )
    END IF ! strDefined(iStrKeyMaterialOut)

    IF (strDefined(iStrKeyMaterialOx)) THEN
      CALL INRT_SetMaterial(global,matOx,strVals(iStrKeyMaterialOx))
    ELSE
! --- Give missing material error if oxidizer used but material undefined
      IF (oxUsed) CALL ErrorStop( global,ERR_INRT_MISSINGMAT,__LINE__ )
! --- set matOx to point to something even if it is not defined
      CALL INRT_SetMaterial(global,matOx,strVals(iStrKeyMaterialIn))
    END IF ! strDefined(iStrKeyMaterialOx)

! - Define interaction (using any relevant information from input deck)

    CALL INRT_DefineBurning(regions(iReg),matIn%index,matOut%index, &
      matOx%index,oxUsed,plagOutExists)

! - Check for switches

! - Which burning model is used

    inrt%switches(INRT_SWI_BURNING_MODEL) = INRT_BURNING_MODEL_DEFAULT

    IF (defined(iKeyModel)) THEN

      SELECT CASE (NINT(vals(iKeyModel)))

      CASE (1)
        inrt%switches(INRT_SWI_BURNING_MODEL) = INRT_BURNING_MODEL_BECKSTEAD

      CASE DEFAULT
        CALL ErrorStop( global,ERR_INRT_BADSWITCH,__LINE__ )

      END SELECT ! vals(iKeyModel)

    END IF ! defined(iKeyModel)

! - Whether oxidizer is used

    IF (oxUsed) THEN
      inrt%switches(INRT_SWI_BURNING_OXUSED) = 1
      nPeulOxEdges = 1
    ELSE
      inrt%switches(INRT_SWI_BURNING_OXUSED) = 0
      nPeulOxEdges = 0
    END IF ! oxUsed

! - Which Vapor Energy method is used

    inrt%switches(INRT_SWI_BURNING_VAPOR_METH) = INRT_BURNING_VAPOR_METH_USED

    IF (defined(iKeyVaporMeth)) THEN

      SELECT CASE (NINT(vals(iKeyVaporMeth)))

      CASE (0)
        inrt%switches(INRT_SWI_BURNING_VAPOR_METH) = &
          INRT_BURNING_VAPOR_METH_NONE

      CASE (1)
        inrt%switches(INRT_SWI_BURNING_VAPOR_METH) = &
          INRT_BURNING_VAPOR_METH_USED

      CASE DEFAULT
        CALL ErrorStop( global,ERR_INRT_BADSWITCH,__LINE__ )

      END SELECT ! vals(iKeyVaporMeth)

    END IF ! defined(iKeyVaporMeth)

! - Check for data

! - temperature above which to shunt burning energy to vapor

    coef = matOut%Tboil

    IF (defined(iKeyVaporTemp)) coef = vals(iKeyVaporTemp)

    IF (coef < 0._RFREAL) CALL ErrorStop( global,ERR_INRT_BADVAL,__LINE__ )

    inrt%data(INRT_DAT_BURNING_VAPOR_TEMP) = coef

! - fraction of heat actually released

    coef = 1._RFREAL ! default heat release coefficient

    IF (defined(iKeyHeatCoef)) coef = vals(iKeyHeatCoef)

    IF (coef < 0._RFREAL) CALL ErrorStop( global,ERR_INRT_BADVAL,__LINE__ )

    inrt%data(INRT_DAT_BURNING_HEAT_COEF) = coef

! - fraction of substance produced that goes immediately back to the particle

    coef = 0._RFREAL ! default mass fraction that goes back to particle

    IF (defined(iKeyMfrcPlag)) coef = vals(iKeyMfrcPlag)

    IF (.NOT. plagOutExists) coef = 0._RFREAL ! no transfer to non-existent Node

    IF (coef < 0._RFREAL) CALL ErrorStop( global,ERR_INRT_BADVAL,__LINE__ )

    inrt%data(INRT_DAT_BURNING_MFRC_PLAG) = coef

! - fraction of substance going to smoke that goes to each smoke type

    nPeulOutEdges = inrt%nEdges - nPeulOxEdges - INRT_BURNING_NEDGES0

    DO iPeulOutEdge = 1,nPeulOutEdges ! loop over edges that output smoke

      iEdge = INRT_BURNING_X_MASS_S0 + nPeulOxEdges + iPeulOutEdge

      iPeul = inrt%edges(iEdge)%iNode(2) - indPeul0

      IF (iPeul < 1 .OR. iPeul > nPeul) &
        CALL ErrorStop( global,ERR_INRT_INDEXRANGE,__LINE__ )

      iKey = iKeyMfrcPeul0 + iPeul
      ind  = INRT_DAT_BURNING_MFRC_PEUL0 + iPeulOutEdge

      coef = 0._RFREAL ! default mass fraction that goes to a smoke type

      IF (defined(iKey)) coef = vals(iKey)

      IF (coef < 0._RFREAL) CALL ErrorStop( global,ERR_INRT_BADVAL,__LINE__ )

      inrt%data(ind) = coef

    END DO ! iPeulOutEdge

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

! - Check for Smoke controls

    DO iPeul = 1,nPeul

      iKeyActv = iKeyPeulActv0 + iPeul
      iKeyPerm = iKeyPeulPerm0 + iPeul
      ind      = indPeul0      + iPeul

      IF (defined(iKeyActv)) &
        CALL INRT_SetActiveness(global,vals(iKeyActv), &
                                inrt%activeness(ind))

      IF (defined(iKeyPerm)) &
        CALL INRT_SetPermission(global,vals(iKeyPerm), &
                                inrt%permission(ind))

    END DO ! iPeul

! - Define Activeness of Internal Node

! - For the Burning interaction, this is defined as the minimum of the
! - Activeness of the Gas and the Lagrangian particles

    inrt%activeness(indIntl) = MIN(inrt%activeness(indMixt), &
                                   inrt%activeness(indPlag0+1))

! - Determine permission Tokens

    CALL INRT_DetermineTokens(regions(iReg),inrt)

  END DO ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE INRT_ReadBurning

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_ReadBurning.F90,v $
! Revision 1.3  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:56:33  fnajjar
! Initial revision after changing case
!
! Revision 1.11  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.10  2004/03/08 21:57:36  jferry
! better error checking for burning without smoke case
!
! Revision 1.9  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.8  2004/03/02 21:47:29  jferry
! Added After Update interactions
!
! Revision 1.7  2003/09/26 21:46:54  fnajjar
! Modified ModInterfaces call to ModInterfacesInteract
!
! Revision 1.6  2003/09/25 15:46:31  jferry
! removed temporary comments
!
! Revision 1.5  2003/09/19 20:35:26  jferry
! Implemented oxidizer species for burning interaction
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







