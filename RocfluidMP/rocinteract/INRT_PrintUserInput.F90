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
! Purpose: write out user input for interactions for checking purposes.
!
! Description: none.
!
! Input: regions = user input.
!
! Output: to standard output.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_PrintUserInput.F90,v 1.3 2008/12/06 08:44:31 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_PrintUserInput( region ) ! PUBLIC

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModInteract
  USE ModError
  USE ModParameters
  USE INRT_ModParameters

  USE ModInterfaces, ONLY : MakeNumberedKeys
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(IN) :: region

! ... loop variables
  INTEGER :: iInrt,iNod,iPeul,iSwi,iDat,iEdge

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  CHARACTER(10)     :: peulStr(100),plagStr(100),tempStr
  CHARACTER(256)    :: edgeStr

  INTEGER :: nNodes,nPlag,nPeul,nInrt,indMixt,indPlag0,indPeul0,indIntl
  INTEGER :: indPlagVapor,maxConEdges,maxDisEdges

  TYPE(t_inrt_input),    POINTER :: input
  TYPE(t_inrt_interact), POINTER :: inrt
  TYPE(t_inrt_edge),     POINTER :: edge
  TYPE(t_global),        POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_PrintUserInput.F90,v $ $Revision: 1.3 $'

  global => region%global

  CALL RegisterFunction( global,'INRT_PrintUserInput',&
  'INRT_PrintUserInput.F90' )

! begin -----------------------------------------------------------------------

  input => region%inrtInput

  nNodes   = input%nNodes
  nPlag    = input%nPlag
  nPeul    = input%nPeul
  indMixt  = input%indMixt
  indPlag0 = input%indPlag0
  indPeul0 = input%indPeul0
  indIntl  = input%indIntl
  indPlagVapor = input%indPlagVapor
  maxConEdges  = input%maxConEdges
  maxDisEdges  = input%maxDisEdges

! Determine total number of interactions used

  nInrt = 0
  DO iInrt = 1,INRT_TYPE_TOTAL
    IF (input%inrts(iInrt)%used) nInrt = nInrt + 1
  END DO ! iInrt

! Make auxillary strings (used in the internal function nodeStr)

  CALL MakeNumberedKeys(plagStr,1,'L',1,nPlag,1)
  CALL MakeNumberedKeys(peulStr,1,'S',1,nPeul,1)

  IF (nInrt > 0) WRITE(STDOUT,*)
  WRITE(STDOUT,1010) SOLVER_NAME//'     Number of interactions',nInrt

  IF (nInrt > 0) THEN

    WRITE(STDOUT,*)
    WRITE(STDOUT,1030) SOLVER_NAME//'       *** Indexing ***'
    WRITE(STDOUT,1010) SOLVER_NAME//'         nNodes      ', nNodes
    WRITE(STDOUT,1010) SOLVER_NAME//'         nPlag       ', nPlag
    WRITE(STDOUT,1010) SOLVER_NAME//'         nPeul       ', nPeul
    WRITE(STDOUT,1010) SOLVER_NAME//'         indMixt     ', indMixt
    WRITE(STDOUT,1010) SOLVER_NAME//'         indPlag0    ', indPlag0
    WRITE(STDOUT,1010) SOLVER_NAME//'         indPeul0    ', indPeul0
    WRITE(STDOUT,1010) SOLVER_NAME//'         indIntl     ', indIntl
    WRITE(STDOUT,1010) SOLVER_NAME//'         indPlagJoint', input%indPlagJoint
    WRITE(STDOUT,1010) SOLVER_NAME//'         indPlagVapor', indPlagVapor
    WRITE(STDOUT,1010) SOLVER_NAME//'         maxConEdges ', maxConEdges
    WRITE(STDOUT,1010) SOLVER_NAME//'         maxDisEdges ', maxDisEdges
    WRITE(STDOUT,*)

    IF (input%consistent) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//'       *** Active Phases '// &
        'Consistent ***'
    ELSE
      WRITE(STDOUT,1030) SOLVER_NAME//'       *** WARNING: Active Phases '// &
        'Inconsistent ***'
    END IF ! input%consistent
    WRITE(STDOUT,*)

    WRITE(STDOUT,1030) SOLVER_NAME//'       *** Global Activeness ***'

    DO iNod = 1,nNodes-1
      WRITE(STDOUT,3010) SOLVER_NAME//'         '//TRIM(nodeStr(iNod)), &
        input%globActiveness(iNod),'  '// &
        TRIM(fullActvStr(input%globActiveness(iNod)))
    END DO ! iNod
    WRITE(STDOUT,*)

    WRITE(STDOUT,1030) SOLVER_NAME//'       *** Global Permission levels ***'

    DO iNod = 1,nNodes-1
      WRITE(STDOUT,3010) SOLVER_NAME//'         '//TRIM(nodeStr(iNod)), &
        input%globPermission(iNod),'  '// &
        TRIM(fullPermStr(input%globPermission(iNod)))
    END DO ! iNod
    WRITE(STDOUT,*)

  END IF ! nInrt

  DO iInrt = 1,INRT_TYPE_TOTAL

    inrt => input%inrts(iInrt)
    IF (.NOT. inrt%used) CYCLE

    WRITE(STDOUT,1010) SOLVER_NAME//'       *** Interaction',iInrt

    WRITE(STDOUT,1030) SOLVER_NAME//'         name          = '//TRIM(inrt%name)
    WRITE(STDOUT,1025) SOLVER_NAME//'         pclsUsed     ',inrt%pclsUsed
    WRITE(STDOUT,1015) SOLVER_NAME//'         order        ',inrt%order
    WRITE(STDOUT,1010) SOLVER_NAME//'         nIntl        ',inrt%nIntl
    WRITE(STDOUT,1010) SOLVER_NAME//'         nInputEdges  ',inrt%nInputEdges
    WRITE(STDOUT,1010) SOLVER_NAME//'         nSwitches    ',inrt%nSwitches
    WRITE(STDOUT,1010) SOLVER_NAME//'         nData        ',inrt%nData
    WRITE(STDOUT,1010) SOLVER_NAME//'         nEdges       ',inrt%nEdges

    DO iSwi = 1,inrt%nSwitches
      WRITE(STDOUT,2010) SOLVER_NAME//'         Switch',iSwi, &
        inrt%switches(iSwi)
    END DO ! iSwi

    DO iDat = 1,inrt%nData
      WRITE(STDOUT,2020) SOLVER_NAME//'         Data  ',iDat,inrt%data(iDat)
    END DO ! iDat

    DO iNod = 1,nNodes-1

      IF (inrt%activeness(iNod) /= input%globActiveness(iNod)) THEN
        WRITE(STDOUT,3010) SOLVER_NAME//'         Activeness override: '// &
          TRIM(nodeStr(iNod)), inrt%activeness(iNod),'  '// &
          TRIM(fullActvStr(inrt%activeness(iNod)))
      END IF ! inrt%activeness(iNod)

      IF (inrt%permission(iNod) /= input%globPermission(iNod)) THEN
        WRITE(STDOUT,3010) SOLVER_NAME//'         Permission override: '// &
          TRIM(nodeStr(iNod)), inrt%permission(iNod),'  '// &
          TRIM(fullPermStr(inrt%permission(iNod)))
      END IF ! inrt%permission(iNod)

    END DO ! iNod

    IF (inrt%nIntl > 0) THEN

      WRITE(STDOUT,3010) SOLVER_NAME//'         Internal Activeness: '// &
        TRIM(nodeStr(indIntl)), inrt%activeness(indIntl),'  '// &
        TRIM(fullActvStr(inrt%activeness(indIntl)))

    END IF ! inrt%nIntl

    DO iEdge = 1,inrt%nEdges

      edge => inrt%edges(iEdge)

      edgeStr = TRIM(nodeStr(edge%iNode(1)))//' '//permStr(edge%token(1))// &
        '-----'//tEdgeStr(edge%tEdge)//'----'//permStr(edge%token(2))

      tempStr = TRIM(nodeStr(edge%iNode(2)))

      IF (tempStr(1:1) == ' ') THEN
        edgeStr = TRIM(edgeStr)//TRIM(tempStr)
      ELSE
        edgeStr = TRIM(edgeStr)//' '//TRIM(tempStr)
      END IF

      WRITE(STDOUT,1030) SOLVER_NAME//'         '//TRIM(edgeStr)

    END DO ! iEdge

    WRITE(STDOUT,*)

  END DO ! iInrt

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1010 FORMAT(A,' =',I3)
1015 FORMAT(A,' =',I9)
1025 FORMAT(A,' = ',L1)
1030 FORMAT(A)
2010 FORMAT(A,I2,' =',I3)
2020 FORMAT(A,I2,' =',ES13.5)
3010 FORMAT(A,':',I3,A)

CONTAINS

  CHARACTER(13) FUNCTION defStr(eq)

    LOGICAL, INTENT(IN) :: eq

    IF (eq) THEN
      defStr = '(default)    '
    ELSE
      defStr = '(non-default)'
    END IF ! eq

  END FUNCTION defStr

  CHARACTER(10) FUNCTION nodeStr(ind)

    INTEGER, INTENT(IN) :: ind

    IF (ind == indMixt) THEN
      nodeStr = ' G'

    ELSE IF (indPlag0 + 1 <= ind .AND. ind <= indPlag0 + nPlag) THEN
      nodeStr = TRIM(plagStr(ind - indPlag0))

    ELSE IF (indPlagVapor == ind) THEN
      nodeStr = 'LV'

    ELSE IF (indPeul0 + 1 <= ind .AND. ind <= indPeul0 + nPeul) THEN
      nodeStr = TRIM(peulStr(ind - indPeul0))

    ELSE IF (ind == indIntl) THEN
      nodeStr = ' X'

    ELSE
      CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

    END IF ! ind

  END FUNCTION nodeStr

  CHARACTER(1) FUNCTION permStr(perm)

    INTEGER, INTENT(IN) :: perm

    SELECT CASE (perm)

    CASE (INRT_PERM_BLOCK)
      permStr = 'x'

    CASE (INRT_PERM_PMASS)
      permStr = 'm'

    CASE (INRT_PERM_PMOME)
      permStr = 'o'

    CASE (INRT_PERM_PALL)
      permStr = '-'

    CASE DEFAULT
      CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

    END SELECT ! perm

  END FUNCTION permStr

  CHARACTER(5) FUNCTION tEdgeStr(tEdge)

    INTEGER, INTENT(IN) :: tEdge

    SELECT CASE (tEdge)

    CASE (INRT_EDGE_BAD)
      tEdgeStr = 'ERROR'

    CASE (INRT_EDGE_MASS)
      tEdgeStr = 'MASS-'

    CASE (INRT_EDGE_MOME_DUM)
      tEdgeStr = 'DUMMY'

    CASE (INRT_EDGE_MOME)
      tEdgeStr = 'MOME-'

    CASE (INRT_EDGE_ENER)
      tEdgeStr = 'ENER-'

    CASE (INRT_EDGE_MASS_GHO)
      tEdgeStr = 'GMASS'

    CASE DEFAULT
      CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

    END SELECT ! tEdge

  END FUNCTION tEdgeStr

  CHARACTER(20) FUNCTION fullActvStr(actv)

    INTEGER, INTENT(IN) :: actv

    SELECT CASE (actv)

    CASE (INRT_ACT_ACTIVE)
      fullActvStr = '(active)'

    CASE (INRT_ACT_ACTIVE-1)
      fullActvStr = '(passive)'

    CASE (:INRT_ACT_ACTIVE-2)
      fullActvStr = '(very passive)'

    CASE DEFAULT
      fullActvStr = '(Error!)'

    END SELECT ! actv

  END FUNCTION fullActvStr

  CHARACTER(40) FUNCTION fullPermStr(perm)

    INTEGER, INTENT(IN) :: perm

    SELECT CASE (perm)

    CASE (INRT_PERM_BLOCK)
      fullPermStr = '(permit nothing)'

    CASE (INRT_PERM_PMASS)
      fullPermStr = '(permit mass only)'

    CASE (INRT_PERM_PMOME)
      fullPermStr = '(permit mass and momentum only)'

    CASE (INRT_PERM_PALL)
      fullPermStr = '(permit all)'

    CASE DEFAULT
      fullPermStr = '(Error!)'

    END SELECT ! perm

  END FUNCTION fullPermStr

END SUBROUTINE INRT_PrintUserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_PrintUserInput.F90,v $
! Revision 1.3  2008/12/06 08:44:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 21:56:30  fnajjar
! Initial revision after changing case
!
! Revision 1.11  2004/07/28 15:42:12  jferry
! deleted defunct constructs: useDetangle, useSmokeDrag, useSmokeHeatTransfer
!
! Revision 1.10  2004/07/27 21:30:00  jferry
! integrated maxConEdges and maxDisEdges variables more fully
!
! Revision 1.9  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.8  2004/03/02 21:48:09  jferry
! First phase of replacing Detangle interaction
!
! Revision 1.7  2003/09/19 20:35:26  jferry
! Implemented oxidizer species for burning interaction
!
! Revision 1.6  2003/04/07 18:26:16  jferry
! minor correction
!
! Revision 1.5  2003/04/03 21:10:18  jferry
! implemented additional safety checks for rocinteract
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







