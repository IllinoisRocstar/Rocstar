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
! Purpose: check interaction information provided by the user
!
! Description: none.
!
! Input: region = data of current region.
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_CheckUserInput.F90,v 1.5 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_CheckUserInput( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInteract
  USE ModError
  USE ModTools,      ONLY : FloatEqual
  USE ModParameters
  USE ModMPI
  USE INRT_ModParameters

  IMPLICIT NONE

! ... parameters
  TYPE (t_region), INTENT(INOUT) :: region

! ... loop variables
  INTEGER :: iInrt,iEdge,iNod,iPlag,iPeul,iPeulOutEdge

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  LOGICAL :: errorFlag,maxConEdgesFound,maxDisEdgesFound

  INTEGER :: nPlag,nPeul,nNodes,nUsedNodes,nPeulOutEdges,nPeulOxEdges
  INTEGER :: ind,indMixt,indPlag0,indPeul0,indIntl,indPlagJoint,indPeulOx
  INTEGER :: indPlagVapor,loclActDiff,maxConEdges,maxDisEdges,globActDiff
  INTEGER :: iwrite

  REAL(RFREAL) :: outmass,coef

  TYPE(t_inrt_input),    POINTER :: input
  TYPE(t_inrt_interact), POINTER :: inrt
  TYPE(t_inrt_edge),     POINTER :: edge
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_CheckUserInput.F90,v $ $Revision: 1.5 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_CheckUserInput',&
  'INRT_CheckUserInput.F90' )

! begin -----------------------------------------------------------------------

  IF (INRT_PERM_PMASS - INRT_PERM_BLOCK /= 1 .OR. &
      INRT_PERM_PMOME - INRT_PERM_PMASS /= 1 .OR. &
      INRT_PERM_PALL  - INRT_PERM_PMOME /= 1)     &
    CALL ErrorStop( global,ERR_INRT_PARAMETER,__LINE__ )

  input => region%inrtInput

  IF (.NOT. input%defaultRead) &
    CALL ErrorStop( global,ERR_INRT_DEFUNREAD,__LINE__ )

  nPlag  = input%nPlag
  nPeul  = input%nPeul
  nNodes = input%nNodes

  indMixt  = input%indMixt
  indPlag0 = input%indPlag0
  indPeul0 = input%indPeul0
  indIntl  = input%indIntl

  indPlagJoint = input%indPlagJoint
  indPlagVapor = input%indPlagVapor

  maxConEdges = input%maxConEdges
  maxDisEdges = input%maxDisEdges
  maxConEdgesFound = .FALSE.
  maxDisEdgesFound = .FALSE.

  iwrite = 1   ! option to write (1) or not (0) warnings

! Check that arrays are allocated an have the correct size

  errorFlag = .FALSE.

  IF (ASSOCIATED(input%globActiveness)) THEN
    IF (UBOUND(input%globActiveness,1) /= nNodes) errorFlag = .TRUE.
  ELSE
    IF (nNodes /= 0) errorFlag = .TRUE.
  END IF ! input%globActiveness

  IF (ASSOCIATED(input%globPermission)) THEN
    IF (UBOUND(input%globPermission,1) /= nNodes) errorFlag = .TRUE.
  ELSE
    IF (nNodes /= 0) errorFlag = .TRUE.
  END IF ! input%globPermission

  IF (ASSOCIATED(input%inrts)) THEN
    IF (UBOUND(input%inrts,1) /= INRT_TYPE_TOTAL) errorFlag = .TRUE.
  ELSE
    IF (INRT_TYPE_TOTAL /= 0) errorFlag = .TRUE.
  END IF ! input%inrts

  IF (errorFlag) CALL ErrorStop( global,ERR_INRT_ALLOCRANGE,__LINE__ )

! Active smoke not implemented for Rocflo
#ifdef RFLO
  DO ind = indPeul0+1,indPeul0+nPeul
    IF (input%globActiveness(ind) == INRT_ACT_ACTIVE) THEN
      CALL ErrorStop( global,ERR_INRT_BADACTV,__LINE__, &
        'Active smoke not implemented for Rocflo' )
    END IF
  END DO ! ind
#endif

! Active species not implemented for Rocflu (but soon?)
#ifdef RFLU
  DO ind = indPeul0+1,indPeul0+nPeul
    IF (input%globActiveness(ind) == INRT_ACT_ACTIVE) THEN
      CALL ErrorStop( global,ERR_INRT_BADACTV,__LINE__, &
        'Active species not implemented for Rocflu' )
    END IF
  END DO ! ind
#endif

! Check that indPlagJoint refers to a particle index

  IF (nPlag > 0) THEN

    IF (indPlagJoint < indPlag0+1 .OR. indPlagJoint > indPlag0+nPlag) &
      CALL ErrorStop( global,ERR_INRT_INDEXRANGE,__LINE__ )

  END IF ! nPlag

! Check that indPlagVapor is assigned correctly

  IF (indPlagVapor /= indPlag0+nPlag+1) &
      CALL ErrorStop( global,ERR_INRT_INDEXRANGE,__LINE__ )

  DO iInrt = 1,INRT_TYPE_TOTAL

    inrt => input%inrts(iInrt)

    IF (ASSOCIATED(inrt%switches)) THEN
      IF (UBOUND(inrt%switches,1) /= inrt%nSwitches) errorFlag = .TRUE.
    ELSE
      IF (inrt%nSwitches /= 0) errorFlag = .TRUE.
    END IF ! switches

    IF (ASSOCIATED(inrt%data)) THEN
      IF (UBOUND(inrt%data,1) /= inrt%nData) errorFlag = .TRUE.
    ELSE
      IF (inrt%nData /= 0) errorFlag = .TRUE.
    END IF ! data

    IF (ASSOCIATED(inrt%activeness)) THEN
      IF (UBOUND(inrt%activeness,1) /= nNodes) errorFlag = .TRUE.
    ELSE
      IF (nNodes /= 0) errorFlag = .TRUE.
    END IF ! inrt%activeness

    IF (ASSOCIATED(inrt%permission)) THEN
      IF (UBOUND(inrt%permission,1) /= nNodes) errorFlag = .TRUE.
    ELSE
      IF (nNodes /= 0) errorFlag = .TRUE.
    END IF ! inrt%permission

    IF (ASSOCIATED(inrt%edges)) THEN
      IF (UBOUND(inrt%edges,1) /= inrt%nEdges) errorFlag = .TRUE.
    ELSE
      IF (inrt%nEdges /= 0) errorFlag = .TRUE.
    END IF ! edges

    IF (errorFlag) CALL ErrorStop( global,ERR_INRT_ALLOCRANGE,__LINE__ )

    IF (.NOT. inrt%used) CYCLE ! do not worry about undefined interactions

! - Check that number of Edges is within proper range

    IF (inrt%pclsUsed) THEN
      IF (inrt%nEdges > maxDisEdges) THEN
        CALL ErrorStop( global,ERR_ILLEGAL_VALUE,__LINE__, &
          'Number of Edges more than maximal number of Edges')
      ELSE IF (inrt%nEdges == maxDisEdges) THEN
        maxDisEdgesFound = .TRUE.
      END IF
    ELSE
      IF (inrt%nEdges > maxConEdges) THEN
        CALL ErrorStop( global,ERR_ILLEGAL_VALUE,__LINE__, &
          'Number of Edges more than maximal number of Edges')
      ELSE IF (inrt%nEdges == maxConEdges) THEN
        maxConEdgesFound = .TRUE.
      END IF
    END IF

! - Check that number of internal Nodes is valid

    IF (inrt%nIntl < 0 .OR. inrt%nIntl > 1) &
      CALL ErrorStop( global,ERR_INRT_NINTL,__LINE__ )

! - Check various things if there is an Internal Node

    IF (inrt%nIntl == 1) THEN

! --- Check that Permission level is 3 (Permit All)

      IF (inrt%permission(indIntl) /= INRT_PERM_PALL) &
        CALL ErrorStop( global,ERR_INRT_PERMLEVINTL,__LINE__ )

! --- Check that there is at least one input and one output Edge

      IF (inrt%nInputEdges <= 0 .OR. inrt%nInputEdges >= inrt%nEdges) &
        CALL ErrorStop( global,ERR_INRT_NINPUTEDGES,__LINE__ )

      DO iEdge = 1,inrt%nInputEdges

        edge => inrt%edges(iEdge)

! ----- Check that input Edge is correctly connected to Internal Node

        IF (edge%iNode(1) == indIntl .OR. edge%iNode(2) /= indIntl) &
          CALL ErrorStop( global,ERR_INRT_CONNECTINTL,__LINE__ )

! ----- Check that input Edge has correct permission Token at Internal Node

        IF (edge%tEdge == INRT_EDGE_MOME_DUM .OR. &
            edge%tEdge == INRT_EDGE_MASS_GHO) THEN

          IF (edge%token(2) /= INRT_PERM_BLOCK) &
            CALL ErrorStop( global,ERR_INRT_PERMINTL,__LINE__ )

        ELSE

          IF (edge%token(2) /= INRT_PERM_PALL) &
            CALL ErrorStop( global,ERR_INRT_PERMINTL,__LINE__ )

        END IF ! edge%tEdge

      END DO ! iEdge

      DO iEdge = inrt%nInputEdges+1,inrt%nEdges

        edge => inrt%edges(iEdge)

! ----- Check that output Edge is correctly connected to Internal Node

        IF (edge%iNode(1) /= indIntl .OR. edge%iNode(2) == indIntl) &
          CALL ErrorStop( global,ERR_INRT_CONNECTINTL,__LINE__ )

! ----- Check that output Edge has correct permission Token at Internal Node

        IF (edge%token(1) /= INRT_PERM_BLOCK) &
          CALL ErrorStop( global,ERR_INRT_PERMINTL,__LINE__ )

      END DO ! iEdge

    END IF ! inrt%nIntl

    nUsedNodes = nNodes - 1 + inrt%nIntl

    DO iNod = 1,nUsedNodes

! --- Check that Activeness of Node has valid value

      IF (inrt%activeness(iNod) > INRT_ACT_ACTIVE) &
        CALL ErrorStop( global,ERR_INRT_BADACTV,__LINE__ )

! --- Warn if Node has been changed from active to passive, or vice-versa
! --- (unless Node is Internal or not the first Lagranigan particle type)

      IF ( iNod /= indIntl .AND. &
          (iNod <= indPlag0+1 .OR. iNod > indPlag0+nPlag+1) ) THEN

        IF (inrt%activeness(iNod)      == INRT_ACT_ACTIVE .NEQV. &
            input%globActiveness(iNod) == INRT_ACT_ACTIVE) THEN

          IF (input%globActiveness(iNod) == INRT_ACT_ACTIVE) THEN
            IF (global%myProcid==MASTERPROC .AND. iwrite==1) &
            WRITE(STDOUT,1030) SOLVER_NAME//'### INRT_WARNING: Node has '// &
              'been changed from active to passive for '//TRIM(inrt%name)
          ELSE
            IF (global%myProcid==MASTERPROC .AND. iwrite==1) &
            WRITE(STDOUT,1030) SOLVER_NAME//'### INRT_WARNING: Node has '// &
              'been changed from passive to active for '//TRIM(inrt%name)
          END IF ! input%globActiveness(iNod)

          IF (global%myProcid==MASTERPROC .AND. iwrite==1) &
          WRITE(STDOUT,1040)   SOLVER_NAME//'### INRT_WARNING:   '// &
            'Node number',iNod

          IF (input%consistent) THEN
            IF (global%myProcid==MASTERPROC .AND. iwrite==1)  &
            WRITE(STDOUT,1030) &
            SOLVER_NAME//'### INRT_WARNING: *** Consistency ruined! ***'
          ENDIF
          input%consistent = .FALSE.
    

        END IF ! inrt%activeness(iNod)

      END IF ! iNod

! --- Check that Permission level of Node has valid value

      IF (inrt%permission(iNod) < INRT_PERM_BLOCK .OR. &
          inrt%permission(iNod) > INRT_PERM_PALL)      &
        CALL ErrorStop( global,ERR_INRT_BADPERM,__LINE__ )

! --- Warn if Permission level less than 3 (Permit All) is on active Node

      IF (inrt%activeness(iNod) == INRT_ACT_ACTIVE .AND. &
          inrt%permission(iNod) /= INRT_PERM_PALL) THEN

        IF (global%myProcid==MASTERPROC .AND. iwrite==1 ) &
        WRITE(STDOUT,1030) SOLVER_NAME//'### INRT_WARNING: Permission '// &
          'restricted for '//TRIM(inrt%name)
        IF (global%myProcid==MASTERPROC .AND. iwrite==1 ) &
        WRITE(STDOUT,1040) SOLVER_NAME//'### INRT_WARNING:   on (active) '// &
          'Node number',iNod

        IF (input%consistent) THEN
          IF (global%myProcid==MASTERPROC .AND. iwrite==1 ) &
          WRITE(STDOUT,1030) &
          SOLVER_NAME//'### INRT_WARNING: *** Consistency ruined! ***'
        ENDIF
        input%consistent = .FALSE.

      END IF ! inrt%activeness(iNod)

    END DO ! iNod

! - Check that all Lagrangian particle constituents have the same Activeness
! - (Do not require Vapor Energy to have same Activeness as others, however)

    DO iPlag=1,nPlag

      IF (inrt%activeness(indPlag0+iPlag) /= inrt%activeness(indPlag0+1)) &
        CALL ErrorStop( global,ERR_INRT_ACTVPLAG,__LINE__ )

    END DO ! iPlag

! - Interaction-specific checks

    SELECT CASE (iInrt)

    CASE (INRT_TYPE_BURNING)

      nPeulOxEdges = 0
      IF (inrt%switches(INRT_SWI_BURNING_OXUSED) /= 0) THEN

        nPeulOxEdges = 1

! ----- Extract index of oxidizer smoke type

        indPeulOx = inrt%edges(INRT_BURNING_S_MASS_X0 + nPeulOxEdges)%iNode(1)

! ----- Check that index is in the smoke range

        IF (indPeulOx < indPeul0 + 1 .OR. indPeulOx > indPeul0 + nPeul) &
          CALL ErrorStop( global,ERR_INRT_INDEXRANGE,__LINE__ )

! ----- Check that if oxidizer is used, it is less active than gas

        IF (input%globActiveness(indPeulOx) >= input%globActiveness(indMixt)) &
          CALL ErrorStop( global,ERR_INRT_OX_ACTV,__LINE__ )

      END IF ! INRT_SWI_BURNING_OXUSED

      IF (inrt%activeness(indIntl) == INRT_ACT_ACTIVE) THEN

! ----- Check that Gas and Lagrangian particles are active

        IF (inrt%activeness(indMixt)    /= INRT_ACT_ACTIVE .OR. &
            inrt%activeness(indPlag0+1) /= INRT_ACT_ACTIVE)     &
          CALL ErrorStop( global,ERR_INRT_BURNING1,__LINE__ )

! ----- Warn if active smoke output has incorrect mass (within a tolerance)

        outmass = 0._RFREAL
        nPeulOutEdges = inrt%nEdges - nPeulOxEdges - INRT_BURNING_NEDGES0

        DO iPeulOutEdge = 1,nPeulOutEdges ! loop over edges that output smoke

          iEdge = INRT_BURNING_X_MASS_S0 + nPeulOxEdges + iPeulOutEdge

          iPeul = inrt%edges(iEdge)%iNode(2) - indPeul0

          IF (iPeul < 1 .OR. iPeul > nPeul) &
            CALL ErrorStop( global,ERR_INRT_INDEXRANGE,__LINE__ )

          IF (inrt%activeness(indPeul0+iPeul) == INRT_ACT_ACTIVE) THEN

            ind = INRT_DAT_BURNING_MFRC_PEUL0 + iPeulOutEdge

            outmass = outmass + inrt%data(ind)

          END IF ! inrt%activeness

        END DO ! iPeulOutEdge

        coef = inrt%data(INRT_DAT_BURNING_MFRC_PLAG)
        outmass = coef + (1._RFREAL - coef)*outmass

        IF (.NOT.FloatEqual(outmass,1._RFREAL)) THEN

          IF (global%myProcid==MASTERPROC .AND. iwrite==1) &
          WRITE(STDOUT,1030) SOLVER_NAME//'### INRT_WARNING: active output '// &
            'mass for '//TRIM(inrt%name)
          IF (global%myProcid==MASTERPROC .AND. iwrite==1) &
          WRITE(STDOUT,1050) SOLVER_NAME//'### INRT_WARNING:   sums not '// &
            'to 1, but to',outmass

          IF (input%consistent) THEN
            IF (global%myProcid==MASTERPROC .AND. iwrite==1) &
            WRITE(STDOUT,1030) &
            SOLVER_NAME//'### INRT_WARNING: *** Consistency ruined! ***'
          ENDIF
          input%consistent = .FALSE.

        END IF ! outmass

      END IF ! inrt%activeness(indIntl)

    END SELECT ! iInrt

  END DO ! iInrt

! Check that maximal number of Edges is correct

  IF (maxConEdges > 0 .AND. .NOT.maxConEdgesFound) THEN
    CALL ErrorStop( global,ERR_ILLEGAL_VALUE,__LINE__, &
      'Inconsistency in maximal number of Edges')
  ENDIF

  IF (maxDisEdges > 0 .AND. .NOT.maxDisEdgesFound) THEN
    CALL ErrorStop( global,ERR_ILLEGAL_VALUE,__LINE__, &
      'Inconsistency in maximal number of Edges')
  ENDIF

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1030 FORMAT(A)
1040 FORMAT(A,I3)
1050 FORMAT(A,ES14.6)

END SUBROUTINE INRT_CheckUserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_CheckUserInput.F90,v $
! Revision 1.5  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/03/06 18:13:44  wasistho
! added parameter iwrite to control printing warnings to stdout
!
! Revision 1.2  2005/03/06 00:18:38  wasistho
! refrain from writing WARNING from each processor
!
! Revision 1.1  2004/12/01 21:56:18  fnajjar
! Initial revision after changing case
!
! Revision 1.9  2004/07/28 15:42:12  jferry
! deleted defunct constructs: useDetangle, useSmokeDrag, useSmokeHeatTransfer
!
! Revision 1.8  2004/07/27 21:30:00  jferry
! integrated maxConEdges and maxDisEdges variables more fully
!
! Revision 1.7  2004/03/03 23:55:40  jferry
! Allowed particles to be run with Euler case
!
! Revision 1.6  2004/03/02 21:47:29  jferry
! Added After Update interactions
!
! Revision 1.5  2003/09/25 15:48:43  jferry
! implemented Boiling Regulation interaction
!
! Revision 1.4  2003/09/19 20:35:26  jferry
! Implemented oxidizer species for burning interaction
!
! Revision 1.3  2003/04/03 21:10:18  jferry
! implemented additional safety checks for rocinteract
!
! Revision 1.2  2003/04/02 22:32:03  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.1  2003/03/11 15:55:02  jferry
! Implemented routine to check user input for Rocinteract
!
!******************************************************************************







