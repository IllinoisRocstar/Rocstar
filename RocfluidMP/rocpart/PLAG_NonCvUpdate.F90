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
! Purpose: update step for non conserved variables.
!
! Description: none.
!
! Input: region = current region
!        iReg   = current region number
!
! Output: region%plag = plag variables
!
! Notes: This corresponds to Part IX Step 26 in RocfluidMP framework.
!
!******************************************************************************
!
! $Id: PLAG_NonCvUpdate.F90,v 1.3 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_nonCvUpdate( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModMixture, ONLY    : t_mixt_input
  USE ModPartLag, ONLY    : t_plag
  USE PLAG_ModInterfaces, ONLY: PLAG_CalcDerivedVariables,   &
                                PLAG_IntrpMixtProperties
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters         
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: region
  
! ... loop variables

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

#ifdef RFLO
  INTEGER :: iLev
#endif

  INTEGER :: iPcls

  INTEGER      :: iFile
  REAL(RFREAL) :: diamL, massL, tauLR

  TYPE(t_plag),   POINTER :: pPlag  
  TYPE(t_global), POINTER :: global
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_NonCvUpdate.F90,v $ $Revision: 1.3 $'
  
  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_NonCvUpdate',&
  'PLAG_NonCvUpdate.F90' )
  
! get dimensions and set pointer ----------------------------------------------

#ifdef RFLO
  iLev  =  region%currLevel
  pPlag => region%levels(iLev)%plag
#endif
#ifdef RFLU
  pPlag => region%plag
#endif

! Calculate derived variables and interpolate mixture properties --------------
!  Active for non-null number of particles in region --------------------------

  IF ( pPlag%nPcls > 0 ) THEN
            
! Get derived variables -------------------------------------------------------

!    WRITE(STDOUT,'(A)') '    Entering PLAG_CalcDerivedVariables'
    CALL PLAG_calcDerivedVariables( region )
  
! Invoke interpolation for mixture properties ---------------------------------

!    WRITE(STDOUT,'(A)') '    Entering PLAG_IntrpMixtProperties'
    CALL PLAG_IntrpMixtProperties( region )

  END IF ! nPcls

! TEMPORARY
! 
!    iFile = 700
!    DO iPcls = 1, 4 
!      massL = SUM( pPlag%cv(pPlag%cvPlagMass(:),iPcls) )
!      diamL = pPlag%dv(DV_PLAG_DIAM,iPcls)
!
!      IF ( diamL > 1.0E-14_RFREAL ) THEN
!        tauLR = 3.0_RFREAL*global%pi*pPlag%tv(TV_PLAG_MUELMIXT,iPcls)*diamL/massL
!      ELSE
!        tauLR = 0.0_RFREAL
!      ENDIF
! 
!      iFile = iFile+1
!      WRITE(iFile,'(1PE12.5,2X,I5,20(2X,1PE28.19))') global%currentTime+ global%dtMin, &
!            pPlag%aiv(AIV_PLAG_ICELLS,iPcls), &
!            pPlag%cv(CV_PLAG_XPOS:CV_PLAG_ZPOS,iPcls), &
!            pPlag%cv(pPlag%cvPlagMass(:),iPcls),&
!            massL,diamL,1.0_RFREAL/tauLR,pPlag%tv(TV_PLAG_MUELMIXT,iPcls),&
!            pPlag%cv(CV_PLAG_XMOM:CV_PLAG_ENER,iPcls), &
!            pPlag%dv(DV_PLAG_UVEL:DV_PLAG_WVEL,iPcls), &
!            pPlag%dv(DV_PLAG_UVELMIXT:DV_PLAG_WVELMIXT,iPcls)
!    END DO ! iPcls 
! 
! END TEMPORARY

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_NonCvUpdate

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_NonCvUpdate.F90,v $
! Revision 1.3  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:53  fnajjar
! Initial revision after changing case
!
! Revision 1.9  2004/11/14 19:47:43  haselbac
! Changed interface
!
! Revision 1.8  2004/02/26 21:02:16  haselbac
! Removed iReg and iStage arguments, commented out writing to files
!
! Revision 1.7  2004/02/06 21:18:20  fnajjar
! Initial Integration of Rocpart with Rocflu
!
! Revision 1.6  2003/11/21 22:43:54  fnajjar
! Commented out PLAG output files
!
! Revision 1.5  2003/04/18 16:13:14  fnajjar
! Added iReg=10 and error trap for tauLR in I/O
!
! Revision 1.4  2003/04/14 20:16:51  fnajjar
! Included iReg=1 only to write debug file
!
! Revision 1.3  2003/04/14 19:46:52  fnajjar
! Added particle diameter and response time to file output
!
! Revision 1.2  2003/03/25 22:54:28  fnajjar
! Added dtMin in WRITE stamp to align with proper time stamp
!
! Revision 1.1  2003/02/04 19:10:11  f-najjar
! Initial Import
!
!******************************************************************************







