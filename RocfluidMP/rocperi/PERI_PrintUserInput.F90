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
! Purpose: Write out PERI user input for checking purposes.
!
! Description: none.
!
! Input: region = user input in current region.
!
! Output: to standard output (monitor).
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PERI_PrintUserInput.F90,v 1.4 2008/12/06 08:44:36 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PERI_PrintUserInput( region ) ! PUBLIC

  USE ModDataTypes
  USE ModGlobal, ONLY     : t_global
  USE ModDataStruct, ONLY : t_region
  USE ModError
  USE ModParameters
  USE PERI_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global

  INTEGER :: flowKind, split(3)
  REAL(RFREAL) :: rgas

! ... cpr and channel local variables
  REAL(RFREAL) :: meanPgrad, bulkmFlux, delta

! ... cpr local variables
  REAL(RFREAL) :: minjRate, cprEpsilon, headPres, headTemp

! ... channel local variables
  INTEGER      :: pgradType
  REAL(RFREAL) :: reTau, utau, ucenter, tauWall, skinFriction

!******************************************************************************

  RCSIdentString = '$RCSfile: PERI_PrintUserInput.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'PERI_PrintUserInput',&
  'PERI_PrintUserInput.F90' )

! get pointers and parameters -------------------------------------------------

! general:
  rgas        = global%refCp*(1._RFREAL - 1._RFREAL/global%refGamma)
  flowKind    = region%periInput%flowKind
  split(1:3)  = region%periInput%split(ICOORD:KCOORD)

! CPR and channel:
  delta       = global%refLength    
  meanPgrad   = region%periInput%meanPgrad
  bulkmFlux   = region%periInput%bulkmFlux   

! CPR:
  minjRate    = region%periInput%minjRate    
  cprEpsilon  = region%periInput%cprEpsilon  
  headPres    = region%periInput%headPres    
  headTemp    = region%periInput%headTemp    

! channel:
  pgradType   = region%periInput%pgradType
  reTau       = region%periInput%cnlRetau    
  utau        = region%periInput%cnlUtau
  ucenter     = region%periInput%cnlCvel

  tauWall     = global%refDensity*utau**2._RFREAL
  skinFriction= 2._RFREAL*(utau/ucenter)**2._RFREAL

! periodic flow case selection

  WRITE(STDOUT,1005) SOLVER_NAME//' Specific flow case, periodic flow:'
  IF (flowKind==PERI_FLOW_CPR) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//' periodic flow kind   = CPR'
  ELSEIF (flowKind==PERI_FLOW_CHANNEL) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//' periodic flow kind   = CHANNEL'
  ELSEIF (flowKind==PERI_FLOW_BOLA) THEN
    WRITE(STDOUT,1030) SOLVER_NAME//' periodic flow kind   = Boundary Layer'
  ENDIF
  WRITE(STDOUT,1010) SOLVER_NAME//' i,j,k parallel split =', &
                     split(ICOORD),split(JCOORD),split(KCOORD), &
                     ' (1 = split domain, 0 = unsplit)'

! model tools selection

  IF (flowKind==PERI_FLOW_CPR) THEN

! - CPR part 

    WRITE(STDOUT,1030) SOLVER_NAME//' input-data relevant to CPR flow:'

    WRITE(STDOUT,1020) SOLVER_NAME//' derived gas constant =',rgas
    WRITE(STDOUT,1020) SOLVER_NAME//' mass injection rate  =',minjRate
    WRITE(STDOUT,1020) SOLVER_NAME//' ratio inj/bulk mrate =',cprEpsilon
    WRITE(STDOUT,1020) SOLVER_NAME//' head end pressure    =',headPres
    WRITE(STDOUT,1020) SOLVER_NAME//' head end temperature =',headTemp
    WRITE(STDOUT,1020) SOLVER_NAME//' lengthscale delta    =',delta
    WRITE(STDOUT,1020) SOLVER_NAME//' based on above parameters:'
    WRITE(STDOUT,1020) SOLVER_NAME//' bulk mass flux       =',bulkmFlux
    WRITE(STDOUT,1020) SOLVER_NAME//' mean pressure grad.  =',meanPgrad

  ELSEIF (flowKind==PERI_FLOW_CHANNEL) THEN

! - CHANNEL part 

    WRITE(STDOUT,1030) SOLVER_NAME//' input-data relevant to Channel flow:'
    IF (pgradType == CNL_PGRAD_TAUWALL) THEN
      WRITE(STDOUT,1030) SOLVER_NAME//' calc. method dpdx    = tau wall based'
    ELSE
      WRITE(STDOUT,1030) SOLVER_NAME//' calc. method dpdx    = mass flux based'
    ENDIF
    WRITE(STDOUT,1020) SOLVER_NAME//' reTau                =',reTau
    WRITE(STDOUT,1020) SOLVER_NAME//' utau                 =',utau
    WRITE(STDOUT,1020) SOLVER_NAME//' ucenter              =',ucenter
    WRITE(STDOUT,1020) SOLVER_NAME//' tauWall              =',tauWall
    WRITE(STDOUT,1020) SOLVER_NAME//' skinFriction         =',skinFriction
    WRITE(STDOUT,1020) SOLVER_NAME//' based on above parameters:'
    WRITE(STDOUT,1020) SOLVER_NAME//' bulk mass flux       =',bulkmFlux
    WRITE(STDOUT,1020) SOLVER_NAME//' mean pressure grad.  =',meanPgrad

  ELSEIF (flowKind==PERI_FLOW_BOLA) THEN
  ENDIF

! finish ----------------------------------------------------------------------

  CALL DeregisterFunction( global )

1005 FORMAT(/,A)
1010 FORMAT(2X,A,3I4,A)
1020 FORMAT(2X,A,E12.5)
1030 FORMAT(2X,A)

END SUBROUTINE PERI_PrintUserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PERI_PrintUserInput.F90,v $
! Revision 1.4  2008/12/06 08:44:36  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:49  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2004/06/11 21:49:40  wasistho
! prepared for RFLU
!
! Revision 1.1  2004/06/08 23:56:56  wasistho
! changed nomenclature
!
! Revision 1.4  2003/09/22 20:18:37  wasistho
! fix screen output of i,j,ksplit
!
! Revision 1.3  2003/09/18 01:57:23  wasistho
! added ijksplit and pgradType in PERI_PgradUpdate
!
! Revision 1.2  2003/04/02 01:49:13  wasistho
! minimize CPR user input
!
! Revision 1.1.1.1  2003/03/29 03:36:30  wasistho
! install ROCPERI
!
!
!
!******************************************************************************







