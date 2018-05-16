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
! ******************************************************************************
!
! Purpose: Compute cavitation source terms.
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region data 
!
! Output: None. 
!
! Notes: Calculate vapor saturation pressure based from Goff (1957). 
!
! ******************************************************************************
!
! $Id: SPEC_RFLU_SourceTerms_GL.F90,v 1.5 2008/12/06 08:44:40 mtcampbe Exp $
!
! Copyright: (c) 2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE SPEC_RFLU_SourceTerms_GL(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters

  USE ModInterfaces, ONLY: MixtLiq_D_DoBpPPoBtTTo, &
                           MixtPerf_Cv_CpR, &
                           MixtPerf_D_PRT, &
                           MixtPerf_R_CpG, &
                           MixtPerf_R_M

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion
  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================
  
  CHARACTER(CHRLEN) :: RCSIdentString  
  INTEGER :: icg
  REAL(RFREAL):: Bp,Bt,CavNo,cvg,cvl,cvv,Dinf,fac,fvol,K,kf,kv,Linf,N,P, &
                 PinfPe,po,Pv,Rg,rhov,rhol,ro,Rv,rYl,Sv,T,Tauf,Tauv,to,Ts, &
                 Vinf,VLinf,vfv
  REAL(RFREAL), DIMENSION(:), POINTER :: vol          
  REAL(RFREAL), DIMENSION(:,:), POINTER :: pCvMixt,pCvSpec,pDvMixt,rhsMixt, &
                                           rhsSpec

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: SPEC_RFLU_SourceTerms_GL.F90,v $ $Revision: 1.5 $'
  
  global => pRegion%global

  CALL RegisterFunction(global,'SPEC_RFLU_SourceTerms_GL',&
  'SPEC_RFLU_SourceTerms_GL.F90')

! ******************************************************************************
! Set pointers and variables
! ******************************************************************************

  pCvMixt => pRegion%mixt%cv
  pCvSpec => pRegion%spec%cv
  pDvMixt => pRegion%mixt%dv

  rhsMixt => pRegion%mixt%rhs
  rhsSpec => pRegion%spec%rhs
  
  vol  => pRegion%grid%vol

! ******************************************************************************
! Define constants
! ******************************************************************************
    
  ro  = global%refDensityLiq
  po  = global%refPressLiq
  to  = global%refTempLiq
  Bp  = global%refBetaPLiq
  Bt  = global%refBetaTLiq
  cvl = global%refCvLiq

  Rg  = MixtPerf_R_M(pRegion%specInput%specType(1)%pMaterial%molw)
  cvg = MixtPerf_Cv_CpR(pRegion%specInput%specType(1)%pMaterial%spht,Rg)

  Rv  = MixtPerf_R_M(pRegion%specInput%specType(2)%pMaterial%molw)
  cvv = MixtPerf_Cv_CpR(pRegion%specInput%specType(2)%pMaterial%spht,Rv)

  fac = 1.5_RFREAL/global%dtMin

  DO icg = 1,pRegion%grid%nCellsTot
    P   = pDvMixt(DV_MIXT_PRES,icg)
    T   = pDvMixt(DV_MIXT_TEMP,icg)
    rYl = pCvMixt(CV_MIXT_DENS,icg) - pCvSpec(1,icg) - pCvSpec(2,icg)
 
    Ts = 273.16_RFREAL ! steam point temperature in K
    Pv = 100.0_RFREAL*(10.0_RFREAL**(10.79574_RFREAL*(1.0_RFREAL-Ts/T) &
                     - 5.02800_RFREAL*LOG10(T/Ts) &
                     + 1.50475E-04_RFREAL*(1.0_RFREAL &
                     - 10.0_RFREAL**(-8.2969_RFREAL*(T/Ts-1.0_RFREAL))) &
                     + 0.42873E-03_RFREAL &
                     * (10.0_RFREAL**(4.76955*(1.0_RFREAL-Ts/T))-1.0_RFREAL) &
                     + 0.78614))       

! ******************************************************************************
!   Model 1 
! ******************************************************************************

    Dinf   = 882.655_RFREAL  ! Hard Coded for Now
    Vinf   = 23.728_RFREAL
    Linf   = 2.0E-03_RFREAL
    VLinf  = Vinf/Linf
    PinfPe = 0.5_RFREAL*Dinf*Vinf*Vinf
    CavNo  = 0.4_RFREAL      ! 0.3, 0.2
    Tauv   = 0.00001_RFREAL  ! base line rate (b=.001), 
    Tauf   = 0.00001_RFREAL  ! faster: bx10, slower: bx0.1 

    kf = ((0.5_RFREAL*( SIGN(1.0_RFREAL,P-Pv) & 
       - 1.0_RFREAL)*(P-Pv))/(Tauf*(PinfPe)))*VLinf
    kv = ((0.5_RFREAL*(-SIGN(1.0_RFREAL,P-Pv) &
       - 1.0_RFREAL)*(P-Pv))/(Tauv*(PinfPe)))*VLinf

    Sv = kf*rYl + kv*pCvSpec(2,icg)

    rhsSpec(2,icg) = rhsSpec(2,icg) - vol(icg)*Sv

! ******************************************************************************
!   Model 2
! ******************************************************************************

!   rhol  = MixtLiq_D_DoBpPPoBtTTo(ro,Bp,Bt,P,po,T,to)
!   rhov  = MixtPerf_D_PRT(P,Rv,T)
!   vfv   = rYl/rhol 
!    
!   N  = 1.0_RFREAL   ! bubble number density
!   K  = 3.85_RFREAL*(rhol/SQRT(rhov))*(N**0.333_RFREAL)
!    
!   Sv = SIGN(1.0_RFREAL,P-Pv)*K*(vfv**0.666_RFREAL)*SQRT(ABS(P-Pv))
!   rhsSpec(2,icg) = rhsSpec(2,icg) - vol(icg)*Sv
  END DO ! icg

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE SPEC_RFLU_SourceTerms_GL

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: SPEC_RFLU_SourceTerms_GL.F90,v $
! Revision 1.5  2008/12/06 08:44:40  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:53  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/07 15:19:25  haselbac
! Removed tabs
!
! Revision 1.2  2006/03/30 20:52:40  haselbac
! Cosmetics
!
! Revision 1.1  2006/03/26 20:21:05  haselbac
! Initial revision
!
! ******************************************************************************







