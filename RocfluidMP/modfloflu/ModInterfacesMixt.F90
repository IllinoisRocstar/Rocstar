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
! Purpose: set explicit interfaces to subroutines and functions
!          related to the gas mixture.
!
! Description: none
!
! Notes: none.
!
!******************************************************************************
!
! $Id: ModInterfacesMixt.F90,v 1.9 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModInterfacesMixt

  IMPLICIT NONE

  INTERFACE

  FUNCTION MixtGasLiq_C(Cvm,D,P,Dl,Dv,Dg,VFl,VFv,VFg,Cl2,Cv2,Cg2,Bl2,Bv2,Bg2)
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: Bg2,Bl2,Bv2,Cl2,Cv2,Cg2,Cvm,D,Dg,Dl,Dv,P,VFg,VFl,VFv
    REAL(RFREAL) :: MixtGasLiq_C
    REAL(RFREAL) ::  denom,numer,term1,term2,term3
  END FUNCTION MixtGasLiq_C

  FUNCTION MixtGasLiq_Eo_CvmTVm2(Cvm,T,Vm2)
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: Cvm,T,Vm2
    REAL(RFREAL) :: MixtGasLiq_Eo_CvmTVm2
  END FUNCTION MixtGasLiq_Eo_CvmTVm2

  FUNCTION MixtGasLiq_P(DYl,DYv,DYg,Cl2,Cv2,Cg2,D,Do,Po,To,Bp,Bt,T)
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: Bp,Bt,Cg2,Cl2,Cv2,D,Do,DYg,DYl,DYv,Po,T,To
    REAL(RFREAL) :: MixtGasLiq_P
    REAL(RFREAL) :: term1,term2
  END FUNCTION MixtGasLiq_P

  FUNCTION MixtLiq_C_Bp(Bp)
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: Bp
    REAL(RFREAL) :: MixtLiq_C_Bp
  END FUNCTION MixtLiq_C_Bp

  FUNCTION MixtLiq_C2_Bp(Bp)
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: Bp
    REAL(RFREAL) :: MixtLiq_C2_Bp
  END FUNCTION MixtLiq_C2_Bp

  FUNCTION MixtLiq_D_DoBpPPoBtTTo(Do,Bp,Bt,P,Po,T,To)
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: Bp,Bt,Do,P,Po,T,To
    REAL(RFREAL) :: MixtLiq_D_DoBpPPoBtTTo
  END FUNCTION MixtLiq_D_DoBpPPoBtTTo

  FUNCTION MixtPerf_C_Co2GUVW( Co2,G,U,V,W )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: Co2, G, U, V, W
    REAL(RFREAL) :: MixtPerf_C_Co2GUVW
  END FUNCTION MixtPerf_C_Co2GUVW

  FUNCTION MixtPerf_C_DGP( D,G,P )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: D, G, P
    REAL(RFREAL) :: MixtPerf_C_DGP
  END FUNCTION MixtPerf_C_DGP

  FUNCTION MixtPerf_C_GHoVm2( G,Ho,Vm2 )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: G, Ho, Vm2
    REAL(RFREAL) :: MixtPerf_C_GHoVm2
  END FUNCTION MixtPerf_C_GHoVm2

  FUNCTION MixtPerf_C_GRT( G,R,T )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: G, R, T
    REAL(RFREAL) :: MixtPerf_C_GRT
  END FUNCTION MixtPerf_C_GRT

  FUNCTION MixtPerf_Co2_CGUVW( C,G,U,V,W )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: C, G, U, V, W
    REAL(RFREAL) :: MixtPerf_Co2_CGUVW
  END FUNCTION MixtPerf_Co2_CGUVW

  FUNCTION MixtPerf_C2_GRT( G,R,T )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: G, R, T
    REAL(RFREAL) :: MixtPerf_C2_GRT
  END FUNCTION MixtPerf_C2_GRT

  FUNCTION MixtPerf_Cv_CpR( Cp,R )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: Cp, R 
    REAL(RFREAL) :: MixtPerf_Cv_CpR
  END FUNCTION MixtPerf_Cv_CpR

  FUNCTION MixtPerf_D_CGP( C,G,P )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: C, G, P
    REAL(RFREAL) :: MixtPerf_D_CGP
  END FUNCTION MixtPerf_D_CGP

  FUNCTION MixtPerf_D_DoGMa(Do,G,Ma)
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: Do,G,Ma
    REAL(RFREAL) :: MixtPerf_D_DoGMa
  END FUNCTION MixtPerf_D_DoGMa

  FUNCTION MixtPerf_D_PRT( P,R,T )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: P, R ,T 
    REAL(RFREAL) :: MixtPerf_D_PRT
  END FUNCTION MixtPerf_D_PRT

  FUNCTION MixtPerf_Eo_DGPUVW( D,G,P,U,V,W )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: D, G, P, U, V, W
    REAL(RFREAL) :: MixtPerf_Eo_DGPUVW
  END FUNCTION MixtPerf_Eo_DGPUVW

  FUNCTION MixtPerf_Eo_DGPVm( D,G,P,Vm )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: D, G, P ,Vm
    REAL(RFREAL) :: MixtPerf_Eo_DGPVm
  END FUNCTION MixtPerf_Eo_DGPVm

  FUNCTION MixtPerf_Eo_GRTUVW( G,R,T,U,V,W )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: G, R, T, U, V, W
    REAL(RFREAL) :: MixtPerf_Eo_GRTUVW
  END FUNCTION MixtPerf_Eo_GRTUVW

  FUNCTION MixtPerf_G_CpR( Cp,R )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: Cp, R
    REAL(RFREAL) :: MixtPerf_G_CpR
  END FUNCTION MixtPerf_G_CpR

  FUNCTION MixtPerf_Ho_CpTUVW( Cp,T,U,V,W )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: Cp, T, U, V, W
    REAL(RFREAL) :: MixtPerf_Ho_CpTUVW
  END FUNCTION MixtPerf_Ho_CpTUVW

  FUNCTION MixtPerf_M_R( R )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: R
    REAL(RFREAL) :: MixtPerf_M_R
  END FUNCTION MixtPerf_M_R

  FUNCTION MixtPerf_P_DEoGVm2( D,Eo,G,Vm2 )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: D, Eo, G, Vm2
    REAL(RFREAL) :: MixtPerf_P_DEoGVm2
  END FUNCTION MixtPerf_P_DEoGVm2

  FUNCTION MixtPerf_P_DRT(D,R,T)
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: D,R,T
    REAL(RFREAL) :: MixtPerf_P_DRT
  END FUNCTION MixtPerf_P_DRT

  FUNCTION MixtPerf_P_GMaPo( G,Ma,Po )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: G, Ma, Po
    REAL(RFREAL) :: MixtPerf_P_GMaPo
  END FUNCTION MixtPerf_P_GMaPo

  FUNCTION MixtPerf_P_DDoGPo(G,Po,D,Do)  
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: D,Do,G,Po
    REAL(RFREAL) :: MixtPerf_P_DDoGPo
  END FUNCTION MixtPerf_P_DDoGPo

  FUNCTION MixtPerf_P_GPoTTo( G,Po,T,To )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: G, Po, T, To
    REAL(RFREAL) :: MixtPerf_P_GPoTTo
  END FUNCTION MixtPerf_P_GPoTTo

  FUNCTION MixtPerf_Po_GPTTo( G,P,T,To )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: G, P, T, To
    REAL(RFREAL) :: MixtPerf_Po_GPTTo
  END FUNCTION MixtPerf_Po_GPTTo

  FUNCTION MixtPerf_Po_CGPUVW( C,G,P,U,V,W )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: C, G, P, U, V, W
    REAL(RFREAL) :: MixtPerf_Po_CGPUVW
  END FUNCTION MixtPerf_Po_CGPUVW

  FUNCTION MixtPerf_R_CpG( Cp,G )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: Cp, G
    REAL(RFREAL) :: MixtPerf_R_CpG
  END FUNCTION MixtPerf_R_CpG

  FUNCTION MixtPerf_R_M( M )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: M
    REAL(RFREAL) :: MixtPerf_R_M
  END FUNCTION MixtPerf_R_M
  
  FUNCTION MixtPerf_T_CGR( C,G,R )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: C, G, R
    REAL(RFREAL) :: MixtPerf_T_CGR
  END FUNCTION MixtPerf_T_CGR  
  
  FUNCTION MixtPerf_T_CpHoVm2(Cp,Ho,Vm2)
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: Cp,Ho,Vm2
    REAL(RFREAL) :: MixtPerf_T_CpHoVm2
  END FUNCTION MixtPerf_T_CpHoVm2
    
  FUNCTION MixtPerf_T_CvEoVm2(Cv,Eo,Vm2)
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: Cv,Eo,Vm2
    REAL(RFREAL) :: MixtPerf_T_CvEoVm2
  END FUNCTION MixtPerf_T_CvEoVm2
  
  FUNCTION MixtPerf_T_DPR( D,P,R )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: D, P, R
    REAL(RFREAL) :: MixtPerf_T_DPR
  END FUNCTION MixtPerf_T_DPR  
  
  FUNCTION MixtPerf_T_GMaTo( G,Ma,To )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: G, Ma, To
    REAL(RFREAL) :: MixtPerf_T_GMaTo
  END FUNCTION MixtPerf_T_GMaTo
  
  FUNCTION MixtPerf_To_CpTUVW(Cp,T,U,V,W)
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: Cp,T,U,V,W
    REAL(RFREAL) :: MixtPerf_To_CpTUVW
  END FUNCTION MixtPerf_To_CpTUVW  
  
  FUNCTION MixtPerf_Vm_C2Co2G( C2,Co2,G )
    USE ModDataTypes
    REAL(RFREAL), INTENT(IN) :: C2, Co2, G
    REAL(RFREAL) :: MixtPerf_Vm_C2Co2G
  END FUNCTION MixtPerf_Vm_C2Co2G  

  SUBROUTINE MixtureProperties( region,inBeg,inEnd,gasUpdate )
    USE ModDataStruct, ONLY : t_region
    INTEGER :: inBeg, inEnd
    LOGICAL :: gasUpdate
    TYPE(t_region) :: region
  END SUBROUTINE MixtureProperties
  
  SUBROUTINE PerfgasDependentVars( inBeg,inEnd,indCp,indMol,cv,gv,dv )
    USE ModDataTypes
    INTEGER :: inBeg, inEnd, indCp, indMol
    REAL(RFREAL), POINTER :: cv(:,:), gv(:,:), dv(:,:)
  END SUBROUTINE PerfgasDependentVars

  SUBROUTINE PerfgasGasVars( inBeg,inEnd,indCp,indMol,refCp,refGamma,cv,gv )
    USE ModDataTypes
    INTEGER :: inBeg, inEnd, indCp, indMol
    REAL(RFREAL)            :: refCp, refGamma
    REAL(RFREAL), POINTER   :: cv(:,:), gv(:,:)
  END SUBROUTINE PerfgasGasVars

  SUBROUTINE PerfgasTransportVars( inBeg,inEnd,indCp,indMol,viscModel,prLam, &
                                   refVisc,refTemp,suthCoef,cv,dv,gv,tv )
    USE ModDataTypes
    INTEGER :: inBeg, inEnd, indCp, indMol, viscModel
    REAL(RFREAL)          :: prLam, refTemp, refVisc, suthCoef
    REAL(RFREAL), POINTER :: cv(:,:), dv(:,:), gv(:,:), tv(:,:)
  END SUBROUTINE PerfgasTransportVars

  END INTERFACE

END MODULE ModInterfacesMixt

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModInterfacesMixt.F90,v $
! Revision 1.9  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.8  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.7  2006/05/01 21:03:25  haselbac
! Added if for MixtPerf_T_CpHoVm2
!
! Revision 1.6  2006/03/26 20:21:54  haselbac
! Added ifs for new routines
!
! Revision 1.5  2005/07/14 21:40:52  haselbac
! Added interface for MixtPerf_To_CpTUVW
!
! Revision 1.4  2005/03/15 20:43:45  haselbac
! Added interface for MixtPerf_D_CGP
!
! Revision 1.3  2004/04/01 21:28:41  haselbac
! Added entry for MixtPerf_E_GRTUVW
!
! Revision 1.2  2003/09/16 15:04:35  haselbac
! Added interfaces for new MixtPerf functions
!
! Revision 1.1  2003/08/11 21:50:00  jblazek
! Splitted ModInterfaces into 4 sections.
!
!******************************************************************************






