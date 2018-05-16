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
! Purpose: check user input for PLAG.
!
! Description: none.
!
! Input: regions = pointer to all regions.
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_CheckUserInput.F90,v 1.7 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CheckUserInput( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY      : t_region
  USE ModGlobal, ONLY          : t_global
  USE ModMixture, ONLY         : t_mixt_input
  USE ModPartLag, ONLY         : t_plag_input
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iCont, iReg

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

  INTEGER :: nRegions

  TYPE(t_global), POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CheckUserInput.F90,v $ $Revision: 1.7 $'

#ifdef RFLO
  global => regions(1)%global
#endif

#ifdef RFLU
  global => regions(0)%global
#endif

  CALL RegisterFunction( global, 'PLAG_CheckUserInput',&
  'PLAG_CheckUserInput.F90' )

! Check for consistency with flow solver -------------------------------------

#ifdef RFLO
  nRegions = global%nRegions
#endif

#ifdef RFLU
  nRegions = global%nRegionsLocal
#endif

  DO iReg = 1, nRegions
    IF ( global%flowType == FLOW_STEADY ) THEN
       WRITE(STDOUT,1000) iReg,global%flowType
       CALL ErrorStop( global,ERR_PLAG_MODULE,__LINE__ )
    END IF ! flowType

    IF ( .NOT. regions(iReg)%mixtInput%computeTv ) THEN
      WRITE(STDOUT,1010) iReg,regions(iReg)%mixtInput%computeTv
      CALL ErrorStop( global,ERR_PLAG_MODULE,__LINE__ )
    END IF ! computeTv

! - viscModel needs to be checked again because computeTv may have been turned
! - on after the previous check (only necessary for RFLO case: the RFLU case
! - this is caught in RFLU_CheckDerivedUserInput)

    IF (  regions(iReg)%mixtInput%computeTv .AND. &
        ( regions(iReg)%mixtInput%viscModel <  VISC_SUTHR .OR.  &
          regions(iReg)%mixtInput%viscModel >  VISC_ANTIB )     ) THEN
      CALL ErrorStop(global,ERR_UNKNOWN_VISCMODEL,__LINE__)
    END IF ! viscModel

#ifdef RFLO
    IF ( regions(iReg)%plagInput%intrplMixtModel /= ZEROTH_ORDER) THEN
      WRITE(STDOUT,1030) iReg,regions(iReg)%plagInput%intrplMixtModel
      CALL ErrorStop( global,ERR_PLAG_INTRPLMODEL,__LINE__ )
    END IF ! intrplMixtModel
#endif 

#ifdef RFLU
    IF ( regions(iReg)%plagInput%intrplMixtModel > FIRST_ORDER) THEN
      WRITE(STDOUT,1030) iReg,regions(iReg)%plagInput%intrplMixtModel
      CALL ErrorStop( global,ERR_PLAG_INTRPLMODEL,__LINE__ )
    END IF ! intrplMixtModel

    IF ( regions(iReg)%plagInput%intrplMixtModel /=ZEROTH_ORDER .AND. &
         regions(iReg)%mixtInput%spaceOrder == 1                      ) THEN 
      CALL ErrorStop(global,ERR_PLAG_INTRPLMODEL,__LINE__)
    END IF ! intrplMixtModel
#endif

    IF ( regions(iReg)%plagInput%breakupModel < PLAG_BREAKUP_NOMODEL .OR. &
         regions(iReg)%plagInput%breakupModel > PLAG_BREAKUP_MODEL1  ) THEN
      WRITE(STDOUT,1040) iReg,regions(iReg)%plagInput%breakupModel
      CALL ErrorStop( global,ERR_PLAG_BREAKUPMODEL,__LINE__ )
    END IF ! breakupModel

    IF ( regions(iReg)%plagInput%breakupWebSwi < PLAG_BREAKUP_NOWEBSWI .OR. &
         regions(iReg)%plagInput%breakupWebSwi > PLAG_BREAKUP_WEBSWI1 ) THEN
      WRITE(STDOUT,1040) iReg,regions(iReg)%plagInput%breakupWebSwi
      CALL ErrorStop( global,ERR_PLAG_BREAKUPWEBSWI,__LINE__ )
    END IF ! breakupModel

    IF ( regions(iReg)%plagInput%injcDiamDist < PLAG_INJC_LOGNORM .OR. &
         regions(iReg)%plagInput%injcDiamDist > PLAG_INJC_PDF) THEN
      WRITE(STDOUT,1050) iReg,regions(iReg)%plagInput%injcDiamDist
      CALL ErrorStop( global,ERR_PLAG_INJCDIAMDIST,__LINE__ )
    END IF ! injcDiamDist

    IF ( regions(iReg)%plagInput%injcDiamDist == PLAG_INJC_LOGSKWD .AND. &
         regions(iReg)%plagInput%injcDiamMean/&
         regions(iReg)%plagInput%injcDiamMax > 0.8_RFREAL ) THEN
      WRITE(STDOUT,1060) iReg,regions(iReg)%plagInput%injcDiamMean/&
                              regions(iReg)%plagInput%injcDiamMax
      CALL ErrorStop( global,ERR_PLAG_INJCDIAM,__LINE__ )
    END IF ! injcDiamDist

    IF ( regions(iReg)%plagInput%injcDiamDist == PLAG_INJC_LOGSKWD .AND. &
         regions(iReg)%plagInput%injcDiamMin/&
         regions(iReg)%plagInput%injcDiamMean > 0.8_RFREAL ) THEN
      WRITE(STDOUT,1070) iReg,regions(iReg)%plagInput%injcDiamMin/&
                              regions(iReg)%plagInput%injcDiamMean
      CALL ErrorStop( global,ERR_PLAG_INJCDIAM,__LINE__ )
    END IF ! injcDiamDist

    IF ( regions(iReg)%plagInput%ejecModel < PLAG_EJEC_MODEL1 .OR. &
         regions(iReg)%plagInput%ejecModel > PLAG_EJEC_CRE  ) THEN
      WRITE(STDOUT,1050) iReg,regions(iReg)%plagInput%ejecModel
      CALL ErrorStop( global,ERR_PLAG_EJECMODEL,__LINE__ )
    END IF ! ejecModel

    IF ( regions(iReg)%plagInput%findPclMethod < FIND_PCL_METHOD_TRAJ_FAST .OR. &
         regions(iReg)%plagInput%findPclMethod > FIND_PCL_METHOD_LOHNER  ) THEN
      WRITE(STDOUT,1050) iReg,regions(iReg)%plagInput%findPclMethod
      CALL ErrorStop( global,ERR_PLAG_FINDPCL,__LINE__ )
    END IF ! ejecModel

  END DO ! iReg

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1000 FORMAT('Region ',I5,', flowType  ',I1)
1010 FORMAT('Region ',I5,', computeTv ',L1)
1030 FORMAT('Region ',I5,', intrplMixtModel ',I1)
1040 FORMAT('Region ',I5,', breakupModel ',I1)
1050 FORMAT('Region ',I5,', injectionModel ',I1)
1060 FORMAT('Region ',I5,', injcDiamMean/injcDiamMax  = ',F12.5)
1070 FORMAT('Region ',I5,', injcDiamMin /injcDiamMean = ',F12.5)

END SUBROUTINE PLAG_CheckUserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CheckUserInput.F90,v $
! Revision 1.7  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2005/11/30 22:18:04  fnajjar
! Added checks for intrplMixtModel
!
! Revision 1.4  2005/07/18 20:46:29  fnajjar
! Bug fix of nRegions for proper handling with MPI
!
! Revision 1.3  2005/04/27 18:35:57  fnajjar
! Included trap error for findPclMethod
!
! Revision 1.2  2005/04/25 18:39:10  luca1
! Imposed PDF from file option for random particle ejection
!
! Revision 1.1  2004/12/01 20:57:22  fnajjar
! Initial revision after changing case
!
! Revision 1.12  2004/07/30 22:47:36  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.11  2004/06/17 15:19:03  fnajjar
! Added infrastructure for ejection model
!
! Revision 1.10  2004/06/17 14:31:25  fnajjar
! Redefined error parameter from ERR_PLAG_INJCMODEL to ERR_PLAG_INCJDIAMDIST
!
! Revision 1.9  2004/06/16 23:03:49  fnajjar
! Renamed variabled for CRE kernel
!
! Revision 1.8  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.7  2004/03/03 23:55:40  jferry
! Allowed particles to be run with Euler case
!
! Revision 1.6  2003/09/18 20:18:42  fnajjar
! Change msg to STDOUT for WRITE statement
!
! Revision 1.5  2003/09/18 20:16:53  fnajjar
! Bug fix for ratio in injection model
!
! Revision 1.4  2003/09/17 21:05:45  fnajjar
! Added error trap for injection model
!
! Revision 1.3  2003/09/13 20:14:21  fnajjar
! Added infrastructure for Breakup model
!
! Revision 1.2  2003/09/10 23:35:50  fnajjar
! Removed flags that are subsumed with Rocinteract
!
! Revision 1.1  2002/10/25 14:14:44  f-najjar
! Initial Import of Rocpart
!
!******************************************************************************







