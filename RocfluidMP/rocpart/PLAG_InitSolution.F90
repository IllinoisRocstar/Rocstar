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
! Purpose: initialize memory for all variables associated with the Lagrangian 
!          particles (PLAG) for all active regions on current processor.
!
! Description: none.
!
! Input: iReg = current region number
!        region = current region
!
! Output: region%plag = plag variables
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_InitSolution.F90,v 1.3 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_initSolution( iReg, region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModError
  USE ModParameters
  USE ModMPI
  USE PLAG_ModParameters
  USE PLAG_ModInterfaces, ONLY : PLAG_CalcDerivedVariables, &
                                 PLAG_IntrpMixtProperties                           
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

  INTEGER :: iReg
  
! ... loop variables
  INTEGER :: iLev, iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  TYPE(t_level),  POINTER :: pLevel
  TYPE(t_plag),   POINTER :: pPlag
  TYPE(t_global), POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InitSolution.F90,v $ $Revision: 1.3 $'
  
  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_InitSolution',&
  'PLAG_InitSolution.F90' )

  IF ( global%myProcid == MASTERPROC   .AND. & 
       global%verbLevel > VERBOSE_NONE .AND. iReg == 1 ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Initializing Solution for PLAG...'
  END IF ! global%verbLevel

! Loop over all grid levels ---------------------------------------------------

  DO iLev=1,region%nGridLevels

    pLevel    => region%levels(iLev)
    pPlag     => region%levels(iLev)%plag

    IF ( pPlag%nPcls > 0 ) THEN
            
! - Load aivOld variables -----------------------------------------------------
 
      DO iPcls = 1, pPlag%nPcls
        pPlag%aivOld(1:pPlag%nAiv,iPcls) = pPlag%aiv(1:pPlag%nAiv,iPcls)
      END DO ! iPcls 

! - Get derived variables -----------------------------------------------------

      CALL PLAG_calcDerivedVariables( region )
  
! - Invoke interpolation for mixture properties -------------------------------

      CALL PLAG_IntrpMixtProperties( region )
    END IF ! nPcls
       
  ENDDO   ! iLev
  
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_InitSolution

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InitSolution.F90,v $
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:39  fnajjar
! Initial revision after changing case
!
! Revision 1.3  2003/04/17 00:04:30  fnajjar
! Included iReg=1 clause for verbosity
!
! Revision 1.2  2002/12/04 15:35:53  f-najjar
! Recoded to fill appropriate datastructure after PLAG_ReadSolution
!
! Revision 1.1  2002/10/25 14:16:31  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************







