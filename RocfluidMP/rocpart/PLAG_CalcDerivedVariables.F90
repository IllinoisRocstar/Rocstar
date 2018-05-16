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
! Purpose: update PLAG derived variables (dv) arrays.
!
! Description: none.
!
! Input: region = current region.
!
! Output: region%levels%plag%dv
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_CalcDerivedVariables.F90,v 1.3 2008/12/06 08:44:33 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_CalcDerivedVariables( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! ... loop variables
  INTEGER :: iCont, iPcls

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString

#ifdef RFLO
  INTEGER :: iLev
#endif
  INTEGER :: nCont, nPcls
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass, pDvPlagVolu
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv

  REAL(RFREAL) :: heatCapSum, massSum, massSumR, &
                  oneThird, pi, piR, voluSum
  REAL(RFREAL), POINTER, DIMENSION(:)   :: pDens, pMixtVol, pSpcHeat
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv, pDv
  
  TYPE(t_plag)  , POINTER :: pPlag
  TYPE(t_global), POINTER :: global  
  
!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_CalcDerivedVariables.F90,v $ $Revision: 1.3 $'
  
  global => region%global
  
  CALL RegisterFunction( global, 'PLAG_CalcDerivedVariables',&
  'PLAG_CalcDerivedVariables.F90' )

! Get dimensions --------------------------------------------------------------

  oneThird = 1.0_RFREAL/3.0_RFREAL  
  pi       = global%pi
  piR      = 1.0_RFREAL/pi

#ifdef RFLO  
  iLev  = region%currLevel
  nPcls = region%levels(iLev)%plag%nPcls
#endif
#ifdef RFLU
  nPcls = region%plag%nPcls
#endif
  nCont = region%plagInput%nCont
  
! Set pointers ----------------------------------------------------------------

#ifdef RFLO  
  pMixtVol  => region%levels(iLev)%grid%vol     
  pPlag     => region%levels(iLev)%plag
#endif
#ifdef RFLU
  pMixtVol  => region%grid%vol
  pPlag     => region%plag
#endif

  pAiv   => pPlag%aiv  
  pCv    => pPlag%cv
  pDv    => pPlag%dv
  
  pCvPlagMass => pPlag%cvPlagMass
  pDvPlagVolu => pPlag%dvPlagVolu
  
  pDens     => region%plagInput%dens
  pSpcHeat  => region%plagInput%spht
  
! Calculate derived variables for Lagrangian field ----------------------------

  DO iPcls = 1, nPcls

!- Compute particle mass ------------------------------------------------------

    massSum  = SUM( pCv(pPlag%cvPlagMass(:),iPcls) )    
    massSumR = 1.0_RFREAL/massSum
    
!- Compute particle volume for each constituent -------------------------------  
  
    DO iCont = 1, nCont
      pDv(pDvPlagVolu(iCont),iPcls) = pCv(pCvPlagMass(iCont),iPcls) / &
                                          pDens(iCont)
    END DO ! iCont   
      
!- Compute particle volume ----------------------------------------------------

    voluSum  = SUM( pDv(pDvPlagVolu(:),iPcls) )
      
!- Compute particle heat capacity ---------------------------------------------
     
    heatCapSum  = SUM(pCv(pCvPlagMass(:),iPcls) * pSpcHeat(:) )

!- Extract derived variables -------------------------------------------------- 
          
    pDv(DV_PLAG_UVEL,iPcls) = pCv(CV_PLAG_XMOM,iPcls) * massSumR
    pDv(DV_PLAG_VVEL,iPcls) = pCv(CV_PLAG_YMOM,iPcls) * massSumR
    pDv(DV_PLAG_WVEL,iPcls) = pCv(CV_PLAG_ZMOM,iPcls) * massSumR
    
    pDv(DV_PLAG_DENS,iPcls) = massSum/voluSum
    pDv(DV_PLAG_DIAM,iPcls) = (6.0_RFREAL*piR*voluSum)**oneThird
    pDv(DV_PLAG_AREA,iPcls) = (9.0_RFREAL*pi *voluSum/16.0_RFREAL)**oneThird    
    
    pDv(DV_PLAG_SPHT,iPcls) = heatCapSum * massSumR
    pDv(DV_PLAG_TEMP,iPcls) = (           pCv(CV_PLAG_ENER,iPcls)*massSumR &
                            - 0.5_RFREAL*(pDv(DV_PLAG_UVEL,iPcls)**2       &
                                         +pDv(DV_PLAG_VVEL,iPcls)**2       &
                                         +pDv(DV_PLAG_WVEL,iPcls)**2 ) )/  &
                               pDv(DV_PLAG_SPHT,iPcls)
    
    pDv(DV_PLAG_VOLU,iPcls) = voluSum       
  END DO  ! iPcls 

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_CalcDerivedVariables

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_CalcDerivedVariables.F90,v $
! Revision 1.3  2008/12/06 08:44:33  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:45  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 20:57:01  fnajjar
! Initial revision after changing case
!
! Revision 1.6  2004/03/26 21:27:03  fnajjar
! Set particle volume in dv array
!
! Revision 1.5  2004/03/02 21:50:03  jferry
! Changed name of DV_PLAG_HTCP to DV_PLAG_SPHT
!
! Revision 1.4  2004/02/06 21:18:20  fnajjar
! Initial Integration of Rocpart with Rocflu
!
! Revision 1.3  2003/01/16 20:15:11  f-najjar
! Removed iRegionGlobal
!
! Revision 1.2  2002/12/03 19:55:37  f-najjar
! Fixed subroutine string name in RegisterFunction
!
! Revision 1.1  2002/10/25 14:14:44  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************







