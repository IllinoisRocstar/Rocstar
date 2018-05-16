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
! Purpose: define various paramters pertinent to TURB
!
! Description: none
!
! Notes: none
!
!******************************************************************************
!
! $Id: TURB_ModParameters.F90,v 1.19 2008/12/06 08:44:42 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE TURB_ModParameters

  USE ModDataTypes
  IMPLICIT NONE

! TURB model and class parameters ---------------------------------------------
  
  INTEGER, PARAMETER :: TURB_MODEL_FIXSMAG = 1, &  ! LES basic Smagorinsky
                        TURB_MODEL_SCALSIM = 2, &  ! LES scale similarity
                        TURB_MODEL_DYNSMAG = 3, &  ! LES dynamic Smagorinsky
                        TURB_MODEL_DYNMIXD = 4, &  ! LES dynamic mixed
                        TURB_MODEL_SA      = 5, &  ! RANS Spalart-Allmaras
                        TURB_MODEL_DESSA   = 6, &  ! DES  Spalart-Allmaras
                        TURB_MODEL_HDESSA  = 7     ! Hybrid-DES SA

  INTEGER, PARAMETER :: MODEL_NONE    = 0, &   ! class of model: NONE
                        MODEL_LES     = 1, &   ! class of model: LES
                        MODEL_RANS    = 2      ! class of model: RANS

! LES parameters
  
  INTEGER, PARAMETER :: MAXOUTFLD_LES = 2      ! maximum number LES outfield
  
  INTEGER, PARAMETER :: FILTYPE_UNIFORM = 0, & 
                        FILTYPE_NONUNIF = 1
  
  INTEGER, PARAMETER :: DELTYPE_CBRT = 0, &    ! cube-root formula
                        DELTYPE_SQRT = 1       ! square-root formula
  
  INTEGER, PARAMETER :: FILWIDTH_ZERO = 0, &   ! no filtering 
                        FILWIDTH_ONE  = 1, &   ! one grid spacing
                        FILWIDTH_TWO  = 2, &   ! twice grid spacing
                        FILWIDTH_FOUR = 4      ! four times grid spacing
  
  INTEGER, PARAMETER :: CALCVORT_NO  = 0, &    ! no vorticities computation 
                        CALCVORT_FDT = 1, &    ! computed per fluid timestep
                        CALCVORT_SDT = 2       ! computed per system timestep
  
  INTEGER, PARAMETER :: CALCWDIST_INI = 0, &   ! comp. only initial wall dist.  
                        CALCWDIST_REM = 1, &   ! computed per remesh
                        CALCWDIST_SDT = 2, &   ! computed per system timestep
                        CALCWDIST_FDT = 3      ! computed per fluid timestep

  INTEGER, PARAMETER :: DV_TURB_CDYN = 1       ! dyn.coef. as turb_dv
   
  INTEGER, PARAMETER :: GR_TURB_UX = 1, &      ! gradients of 
                        GR_TURB_VX = 2, &      ! test filtered velocities
                        GR_TURB_WX = 3, &       
                        GR_TURB_UY = 4, &       
                        GR_TURB_VY = 5, &       
                        GR_TURB_WY = 6, &       
                        GR_TURB_UZ = 7, &       
                        GR_TURB_VZ = 8, &       
                        GR_TURB_WZ = 9       

  INTEGER, PARAMETER :: CV_TURB_DENS = 1, &    ! components of filtered cv
                        CV_TURB_XMOM = 2, &    ! (defined at faces)
                        CV_TURB_YMOM = 3, &    
                        CV_TURB_ZMOM = 4, &
                        CV_TURB_NELM = 4

  INTEGER, PARAMETER :: CV_TURB_UVEL = 2, &    ! components of filtered 
                        CV_TURB_VVEL = 3, &    ! primitive variables 
                        CV_TURB_WVEL = 4       ! (defined at cells)

! RANS parameters
  
  INTEGER, PARAMETER :: MAXOUTFLD_RANS = 3     ! maximum number LES outfield

  INTEGER, PARAMETER :: CV_SA_NUTIL = 1, &     ! components cv of SA model
                        CV_SA_NELM  = 1

  INTEGER, PARAMETER :: MC_SA_CB1   = 1, &     ! Spalart-Allmaras  
                        MC_SA_CB2   = 2, &     ! model constants
                        MC_SA_CW1   = 3, &  
                        MC_SA_CW2   = 4, &  
                        MC_SA_CW3   = 5, &  
                        MC_SA_CV1   = 6, &  
                        MC_SA_RSIG  = 7, &  
                        MC_SA_RKAP  = 8, &
                        MC_SA_NELM  = 8     

  INTEGER, PARAMETER :: GR_SA_NUTILX = 1, &    ! components grad(tilde[nu])
                        GR_SA_NUTILY = 2, &
                        GR_SA_NUTILZ = 3

  INTEGER, PARAMETER :: TVT_RANS_MUE  = 1, &   ! components of RaNS tv
                        TVT_RANS_TCO  = 2, &
                        TVT_RANS_NELM = 2

  INTEGER, PARAMETER :: WDIST_DIRECT  = 0, &   ! wall dist.comp. method
                        WDIST_HIERAR  = 1      ! for relevant RaNS models 

  INTEGER, PARAMETER :: RANS_DISCR_CEN = 0, &  ! central discretization
                        RANS_DISCR_UPW = 1     ! upwind

  INTEGER, PARAMETER :: RANS_DISCR_ORD1 = 1, & ! discretization order
                        RANS_DISCR_ORD2 = 2         

  INTEGER, PARAMETER :: SA_FV1_POW3 = 0, &     ! formula for viscous function
                        SA_FV1_POW2 = 1         

! Wall Layer Model parameters

  INTEGER, PARAMETER :: WLM_NSWITCH        = 2     ! number of wlm input switch

  INTEGER, PARAMETER :: WLM_INPUT_MODEL    = 1, &  ! wlm model parameter
                        WLM_INPUT_REFPOINT = 2     ! wlm reference point

  INTEGER, PARAMETER :: WLM_MODEL_NOMODEL  = 0, &  ! wlm model options
                        WLM_MODEL_LOGLAY   = 1, &
                        WLM_MODEL_BNDLAY   = 2, &
                        WLM_MODEL_EXTERN   = 3

  INTEGER, PARAMETER :: WLM_VALS_XIX       =  1, & ! Sij/|S| 1st cell from wall
                        WLM_VALS_ETX       =  2, &
                        WLM_VALS_ZTX       =  3, &
                        WLM_VALS_XIY       =  4, &
                        WLM_VALS_ETY       =  5, &
                        WLM_VALS_ZTY       =  6, &
                        WLM_VALS_XIZ       =  7, &
                        WLM_VALS_ETZ       =  8, &
                        WLM_VALS_ZTZ       =  9, &

                        WLM_VALS_TAUUX     = 10, & ! Tauij at ns-patches
                        WLM_VALS_TAUUY     = 11, &
                        WLM_VALS_TAUUZ     = 12, &
                        WLM_VALS_TAUVX     = 13, &
                        WLM_VALS_TAUVY     = 14, &
                        WLM_VALS_TAUVZ     = 15, &
                        WLM_VALS_TAUWX     = 16, &
                        WLM_VALS_TAUWY     = 17, &
                        WLM_VALS_TAUWZ     = 18, &

                        WLM_VALS_XIV       = 19, & ! bl vel. 1st cell from wall
                        WLM_VALS_ETV       = 20, &
                        WLM_VALS_ZTV       = 21, &

                        WLM_VALS_ROUGH     = 22, & ! roughness size
                        WLM_VALS_UTAU      = 23, & ! friction velocity
                        WLM_VALS_WDIST     = 24, & ! 1st wall distance
                        WLM_VALS_DXI       = 25, & ! xi grid spacing
                        WLM_VALS_DZT       = 26, & ! zeta grid spacing
                        WLM_VALS_DPDXI     = 27, & ! dp/dxi
                        WLM_VALS_DPDZT     = 28, & ! dp/zeta
                        WLM_VALS_DENS      = 29, & ! density
                        WLM_VALS_HFLUX     = 30    ! wall heat flux

! TURB general parameters

  INTEGER, PARAMETER :: ST_TURB_VAR1  = 1, &       ! quantities to be timeavg
                        ST_TURB_VAR2  = 2, &           
                        ST_TURB_VAR3  = 3, &           
                        ST_TURB_VAR4  = 4, &           
                        ST_TURB_VAR5  = 5, &           
                        ST_TURB_VAR6  = 6, &           
                        ST_TURB_NVAR  = 6           

  INTEGER, PARAMETER :: DIRI = 1, &                ! space vector components
                        DIRJ = 2, & 
                        DIRK = 3, &
                        NDIR = 3

  INTEGER, PARAMETER :: E11 = 1, &                 ! symm tensor components
                        E12 = 2, & 
                        E13 = 3, &
                        E22 = 4, &
                        E23 = 5, &
                        E33 = 6, &
                        TENSOR_SYMM_NELM = 6       ! elm number of symm tensor       
  INTEGER, PARAMETER :: A11 = 1, &                 ! whole tensor components
                        A12 = 2, &
                        A13 = 3, &
                        A21 = 4, &
                        A22 = 5, &
                        A23 = 6, &
                        A31 = 7, &
                        A32 = 8, &
                        A33 = 9, &
                        TENSOR_ALL_NELM  = 9       ! elm number of whole tensor     
  INTEGER, PARAMETER :: ZOF_LES_EDDYVIS = 1, &     ! zero-one switch fields
                        ZOF_NELM        = 1

   
  REAL(RFREAL), PARAMETER :: REAL_SMALL    = 1.E-16_RFREAL ! small real nmbr
 
  REAL(RFREAL), PARAMETER :: STOP_VALUE    = 1.E-12_RFREAL ! stop criterium

END MODULE TURB_ModParameters

!******************************************************************************
!
! RCS Revision history:
!
! $Log: TURB_ModParameters.F90,v $
! Revision 1.19  2008/12/06 08:44:42  mtcampbe
! Updated license.
!
! Revision 1.18  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.17  2006/01/17 17:52:32  wasistho
! ZOF_LES_FIXSMAG to ZOF_LES_EDDYVIS
!
! Revision 1.16  2006/01/12 09:49:33  wasistho
! enabled tripping fixed Smagorinsky
!
! Revision 1.15  2005/03/07 05:03:37  wasistho
! install hybrid DESSA turbulence model
!
! Revision 1.14  2004/04/20 20:49:47  wasistho
! added user option for frequency in computing wall distance
!
! Revision 1.13  2004/04/08 20:22:32  wasistho
! increased REAL_SMALL value
!
! Revision 1.12  2004/03/05 21:08:02  wasistho
! changed nomenclature
!
! Revision 1.11  2004/02/19 04:03:18  wasistho
! added new rans/SA parameter VISCFUNCTION
!
! Revision 1.10  2004/02/14 03:42:40  wasistho
! added new WLM parameter: reference point
!
! Revision 1.9  2004/02/11 03:24:32  wasistho
! added feature: variable number of turbulence output fields
!
! Revision 1.8  2003/10/26 00:18:26  wasistho
! added multiple discr.types and order
!
! Revision 1.7  2003/10/07 02:04:30  wasistho
! initial installation of RaNS-SA and DES
!
! Revision 1.6  2003/08/06 15:56:22  wasistho
! added vorticities computation
!
! Revision 1.5  2003/08/02 00:20:02  wasistho
! added parameters for calcVort
!
! Revision 1.4  2003/06/05 19:18:44  wasistho
! implemented heat transfer model
!
! Revision 1.3  2003/05/31 01:45:55  wasistho
! installed turb. wall layer model
!
! Revision 1.2  2003/05/24 02:08:15  wasistho
! turbulence statistics expanded
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************






