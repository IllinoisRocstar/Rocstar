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
! Purpose: initialize memory for all variables associated 
!          with the Lagrangian particles (PLAG).
!
! Description: none.
!
! Input: region = current region
!
! Output: region%plag = plag variables
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_InitMemory.F90,v 1.5 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_InitMemory( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region, t_level
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModError
  USE ModParameters
  USE ModMPI
  USE PLAG_ModParameters
#ifdef STATS
  USE PLAG_ModEulerian, ONLY: PLAG_InitEulerianField
  USE PLAG_ModStats, ONLY: PLAG_InitStat
#endif

  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region
  
! ... loop variables
  INTEGER :: iCont, iLev

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: nAiv, nArv, nCont, nCv, nDv, nTv
  INTEGER :: iAiv, iArv, iCv, iDv, iPcl, iTv, nPclsMax

  TYPE(t_level),  POINTER :: pLevel
  TYPE(t_plag),   POINTER :: pPlag
  TYPE(t_global), POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_InitMemory.F90,v $ $Revision: 1.5 $'
  
  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_InitMemory',&
  'PLAG_InitMemory.F90' )

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_NONE ) THEN 
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Initializing memory for PLAG...'
  END IF ! global%verbLevel

! Loop over all grid levels ---------------------------------------------------

  DO iLev=1,region%nGridLevels

    pLevel     => region%levels(iLev)
    pPlag      => region%levels(iLev)%plag

    nCont     = region%plagInput%nCont
    nPclsMax  = region%plagInput%nPclsMax
    nAiv      = region%levels(iLev)%plag%nAiv
    nArv      = region%levels(iLev)%plag%nArv
    nCv       = region%levels(iLev)%plag%nCv
    nDv       = region%levels(iLev)%plag%nDv
    nTv       = region%levels(iLev)%plag%nTv
    
    DO iCont = 1, nCont
      WRITE(STDOUT,'(A,3(I2,3X))') 'iCont, cvPlagMass, dvPlagVolu ', &
                                    iCont, pPlag%cvPlagMass(iCont),   &
                                    pPlag%dvPlagVolu(iCont)
    END DO ! iCont   

    pPlag%nPcls       = 0
    pPlag%nextIdNumber= 0
    pPlag%nPclsMax    = nPclsMax

    DO iPcl = 1, nPclsMax 
      DO iAiv = 1, nAiv
!        pPlag%aiv(iAiv,iPcl)    = 0
!        pPlag%aivOld(iAiv,iPcl) = 0
        pPlag%aiv(iAiv,iPcl)    = CRAZY_VALUE_INT 
        pPlag%aivOld(iAiv,iPcl) = CRAZY_VALUE_INT
      END DO ! iAiv

      DO iArv = 1, nArv
!        pPlag%arv(iArv,iPcl)    = 0.0_RFREAL
!        pPlag%arvOld(iArv,iPcl) = 0.0_RFREAL
        pPlag%arv(iArv,iPcl)    = -1.0E+30_RFREAL
        pPlag%arvOld(iArv,iPcl) = -1.0E+30_RFREAL
      END DO ! iArv

      DO iCv  = 1, nCv
!        pPlag%cv(iCv,iPcl)     = 0.0_RFREAL
!        pPlag%rhs(iCv,iPcl)    = 0.0_RFREAL
!        pPlag%rhsSum(iCv,iPcl) = 0.0_RFREAL 
!        pPlag%cvOld(iCv,iPcl)     = 0.0_RFREAL

        pPlag%cv(iCv,iPcl)     = -1.0E+30_RFREAL
        pPlag%rhs(iCv,iPcl)    = -1.0E+30_RFREAL
        pPlag%rhsSum(iCv,iPcl) = -1.0E+30_RFREAL 
        pPlag%cvOld(iCv,iPcl)  = -1.0E+30_RFREAL
      END DO ! iCv

      DO iDv  = 1, nDv
!        pPlag%dv(iDv,iPcl)     = 0.0_RFREAL
        pPlag%dv(iDv,iPcl)     = -1.0E+30_RFREAL
      END DO ! iDv

      DO iTv  = 1, nTv
!        pPlag%tv(iTv,iPcl)     = 0.0_RFREAL
        pPlag%tv(iTv,iPcl)     = -1.0E+30_RFREAL
      END DO ! iTv
    ENDDO ! iPcl
    
! - Initialize memory for statistics -------------------------------------------
    
#ifdef STATS  
    IF ( ( global%flowType == FLOW_UNSTEADY ) .AND. &
         ( global%doStat == ACTIVE )                ) THEN
      CALL PLAG_InitEulerianField( region, pPlag )
      CALL PLAG_InitStat( region, pPlag )
    END IF ! global%flowType        
#endif
 
  ENDDO   ! iLev
  
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_InitMemory

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_InitMemory.F90,v $
! Revision 1.5  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/03/06 23:13:13  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.2  2005/01/08 20:42:13  fnajjar
! Included initialization calls for PLAG statistics
!
! Revision 1.1  2004/12/01 20:57:38  fnajjar
! Initial revision after changing case
!
! Revision 1.5  2004/07/28 18:58:13  fnajjar
! Included new definition for nPclsTot from dynamic memory reallocation
!
! Revision 1.4  2004/04/09 23:11:23  fnajjar
! Made initialization within DO loop construct and added crazy values for proper trapping
!
! Revision 1.3  2003/11/21 22:43:18  fnajjar
! Removed nPclsTot and added nextIdNumber
!
! Revision 1.2  2003/05/27 19:19:57  fnajjar
! Removed distPartBurning and all pertinent LOGICAL datastructure
!
! Revision 1.1  2002/10/25 14:16:31  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************







