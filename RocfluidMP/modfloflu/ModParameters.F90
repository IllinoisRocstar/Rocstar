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
! Purpose: define various parameters
!
! Description: none
!
! Notes:
!   1. DO NOT CHANGE any of the parameters values without good reason and 
!      really knowing what you are doing... 
!   2. Error codes are defined in the module ModError.F90
!
! ******************************************************************************
!
! $Id: ModParameters.F90,v 1.165 2010/02/18 21:47:40 juzhang Exp $
!
! Copyright: (c) 2001-2006 by the University of Illinois
!
! ******************************************************************************

MODULE ModParameters

  IMPLICIT NONE

! ******************************************************************************
! Geometry 
! ******************************************************************************

  INTEGER, PARAMETER :: XCOORD = 1, &
                        YCOORD = 2, &
                        ZCOORD = 3, & 
                        XYZMAG = 4

  INTEGER, PARAMETER :: ICOORD = 1, &
                        JCOORD = 2, &
                        KCOORD = 3

  INTEGER, PARAMETER :: REGION_SHAPE_NORMAL = 0, &
                        REGION_SHAPE_FUNKY  = 1

! ******************************************************************************
! Data structure 
! ******************************************************************************                                  

  INTEGER, PARAMETER :: CELL_TYPE_EXT = -99, & ! Exterior cell (parallel)
                        CELL_TYPE_BND =   0, & ! Boundary cell
                        CELL_TYPE_TET =   1, & ! Tetrahedron
                        CELL_TYPE_HEX =   2, & ! Hexahedron
                        CELL_TYPE_PRI =   3, & ! Prism
                        CELL_TYPE_PYR =   4    ! Pyramid                    
                  
  INTEGER, PARAMETER :: FACE_TYPE_TRI  = 1, & 
                        FACE_TYPE_QUAD = 2                     
                        
  INTEGER, PARAMETER :: FACE_TYPE_NEW = 1, & 
                        FACE_TYPE_OLD = 2                                              

  INTEGER, PARAMETER :: EDGE_TYPE_NEW = 1, & 
                        EDGE_TYPE_OLD = 2                      

  INTEGER, PARAMETER :: FACE_SPLIT_13  = 1, & 
                        FACE_SPLIT_24  = 2, & 
                        FACE_SPLIT_YES = 1, & 
                        FACE_SPLIT_NO  = 2
                        
  INTEGER, PARAMETER :: FACE_NONE     =  0, & ! Must be zero or less
                        OPP_FACE_NONE = -1    ! Must be negative                         

  INTEGER, PARAMETER :: VERT_NONE         = 0, & ! Must be zero because of Charm
                        VERT_KIND_ACTUAL  = 1, & 
                        VERT_KIND_VIRTUAL = 2, & 
                        VERT_KIND_A       = 1, & 
                        VERT_KIND_V       = 2, & 
                        VERT_KIND_AV      = 3
                        
  INTEGER, PARAMETER :: BVERT_CLASS_OUTSIDER = 0, & ! Must be zero or less
                        BVERT_CLASS_INSIDER  = 1    ! Must be positive                  
                        
  INTEGER, PARAMETER :: CELL_KIND_EXT     = -99, & ! Exterior cell (parallel)
                        CELL_KIND_BND     =   0, &
                        CELL_KIND_ACTUAL  =   1, & 
                        CELL_KIND_VIRTUAL =   3

  INTEGER, PARAMETER :: FACE_KIND_AA = 1, & ! Actual-Actual
                        FACE_KIND_AV = 2, & ! Actual-Virtual
                        FACE_KIND_VV = 3, & ! Virtual-Virtual
                        FACE_KIND_VB = 4, & ! Virtual-Boundary
                        FACE_KIND_VX = 5, & ! Virtual-Exterior
                        FACE_KIND_AB = 6    ! Actual-Boundary
                        
  INTEGER, PARAMETER :: EDGE_KIND_AA = 1, & ! Actual-Actual
                        EDGE_KIND_AV = 2, & ! Actual-Virtual
                        EDGE_KIND_VV = 3    ! Virtual-Virtual

  INTEGER, PARAMETER :: DIFF_LOC_CELL = 1, & 
                        DIFF_LOC_FACE = 2, & 
                        DIFF_LOC_VERT = 3
 
  INTEGER, PARAMETER :: DIFF_DEGREE_GRAD = 1, & 
                        DIFF_DEGREE_HESS = 2
 
  INTEGER, PARAMETER :: GRAD_TYPE_UNCONSTRAINED = 0, & 
                        GRAD_TYPE_CONSTRAINED   = 1
 
  INTEGER, PARAMETER :: INTERP_LOC_CELL = 1, & 
                        INTERP_LOC_FACE = 2, &
                        INTERP_LOC_VERT = 3
 
  INTEGER, PARAMETER :: C2F_INIT = -1, & ! Must be negative
                        C2V_INIT =  0, & ! Must be zero or less
                        F2C_INIT =  0                         
  
  INTEGER, PARAMETER :: C2FN_INIT = 0, & ! Must be zero or less
                        C2CS_INIT = 0, & ! Must be zero or less  
                        F2CS_INIT = 0, & ! Must be zero or less 
                        V2CS_INIT = 0    ! Must be zero or less  
 
  INTEGER, PARAMETER :: V2C_BEG = 1, & 
                        V2C_END = 2
                        
  INTEGER, PARAMETER :: X2CS_LAYER_BEG = 1, & 
                        X2CS_LAYER_END = 2                      
 
  INTEGER, PARAMETER :: DERIV_DEGREE_0 = 0, & 
                        DERIV_DEGREE_1 = 1 
 
  INTEGER, PARAMETER :: INT_LIM_LOW = 1, & 
                        INT_LIM_UPP = 2
 
  INTEGER, PARAMETER :: REGION_INDEX_OFFSET = 100
 
#ifdef RFLO 
  INTEGER, PARAMETER :: DEGENERAT_NONE          =  0, & ! degenerative  
                        DEGENERAT_EDGE_IN_PATCH =  1, & ! edges/corners      
                        DEGENERAT_CORN_IN_EDGE  =  1, &                     
                        DEGENERAT_CORN_IN_PATCH =  2, &                     
                        DEGENERAT_DETACH        = -1                        

  INTEGER, PARAMETER :: EDGE_INTERACT_FULL      =  0, & ! all edgeCells interior
                        EDGE_INTERACT_PART      =  1
#endif 
 
! ******************************************************************************
! Coordinate moments 
! ******************************************************************************

  INTEGER, PARAMETER :: XYZ_MOM_11 =  1, & 
                        XYZ_MOM_12 =  2, & 
                        XYZ_MOM_22 =  3, & 
                        XYZ_MOM_13 =  4, & 
                        XYZ_MOM_23 =  5, & 
                        XYZ_MOM_33 =  6, & 
                        XYZ_MOM_14 =  7, & 
                        XYZ_MOM_24 =  8, & 
                        XYZ_MOM_34 =  9, & 
                        XYZ_MOM_44 = 10
 
! ******************************************************************************
! Mixture 
! ******************************************************************************

  INTEGER, PARAMETER :: CV_MIXT_DENS = 1, &    ! density (rho)
                        CV_MIXT_XMOM = 2, &    ! rho * u
                        CV_MIXT_YMOM = 3, &    ! rho * v
                        CV_MIXT_ZMOM = 4, &    ! rho * w
                        CV_MIXT_ENER = 5, &    ! rho * E
                        CV_MIXT_NEQS = 5       ! total no. of equations

#ifdef RFLU
  INTEGER, PARAMETER :: CV_MIXT_XVEL = 2, &    ! u
                        CV_MIXT_YVEL = 3, &    ! v
                        CV_MIXT_ZVEL = 4, &    ! w
                        CV_MIXT_PRES = 5, &    ! p
                        CV_MIXT_TEMP = 5       ! T
                        
  INTEGER, PARAMETER :: CV_MIXT_STATE_CONS  = 1, &
                        CV_MIXT_STATE_PRIM  = 2, & 
                        CV_MIXT_STATE_DUVWP = 3, & 
                        CV_MIXT_STATE_DUVWT = 4
                        
  INTEGER, PARAMETER :: VAR_INFO_POS    = 1, & 
                        VAR_INFO_POSNEG = 2                      
#endif 

#ifdef RFLO
  INTEGER, PARAMETER :: DV_MIXT_UVEL = 1, &    ! velocity components
                        DV_MIXT_VVEL = 2, &
                        DV_MIXT_WVEL = 3, &
                        DV_MIXT_TEMP = 4, &    ! static temperature
                        DV_MIXT_PRES = 5, &    ! static pressure
                        DV_MIXT_SOUN = 6       ! speed of sound
#endif
#ifdef RFLU
  INTEGER, PARAMETER :: DV_MIXT_PRES = 1, &    ! static pressure
                        DV_MIXT_TEMP = 2, &    ! static temperature
                        DV_MIXT_SOUN = 3, &     ! speed of sound
                        DV_MIXT_NVAR = 3
#endif

  INTEGER, PARAMETER :: TV_MIXT_MUEL = 1, &    ! laminar viscosity
                        TV_MIXT_TCOL = 2, &    ! laminar thermal conductivity
                        TV_MIXT_MUET = 3, &    ! turbulent viscosity
                        TV_MIXT_TCOT = 4       ! turbulent thermal conductivity

  INTEGER, PARAMETER :: GV_MIXT_CP   = 1, &    ! heat coeff. at const. pressure
                        GV_MIXT_MOL  = 2       ! molecular weight

#ifdef RFLU
  INTEGER, PARAMETER :: PV_MIXT_SCHL =  1, & 
                        PV_MIXT_SHAD =  2, & 
                        PV_MIXT_INTF =  3, & 
                        PV_MIXT_XVOR =  4, & 
                        PV_MIXT_YVOR =  5, & 
                        PV_MIXT_ZVOR =  6, &
                        PV_MIXT_VCI2 =  7, &
                        PV_MIXT_VCL2 =  8, & 
                        PV_MIXT_VCLI =  9, &
                        PV_MIXT_VCLR = 10, &
                        PV_MIXT_GREX = 11, & 
                        PV_MIXT_GREY = 12, & 
                        PV_MIXT_GREZ = 13, &  
                        PV_PLAG_DIA3 = 14, & 
                        PV_PLAG_DIA4 = 15, & 
                        PV_PLAG_NDNS = 16, & 
                        PV_PLAG_XVEL = 17, & 
                        PV_PLAG_YVEL = 18, & 
                        PV_PLAG_ZVEL = 19, & 
                        PV_PLAG_TEMP = 20, & 
                        PV_PLAG_MASS = 21, &   
                        PV_XXXX_NVAR = 21
#endif

#ifdef RFLO
  INTEGER, PARAMETER :: GR_MIXT_UX   =  1, &   ! Gradients: du/dx
                        GR_MIXT_VX   =  2, &   ! dv/dx
                        GR_MIXT_WX   =  3, &   ! dw/dx
                        GR_MIXT_TX   =  4, &   ! dT/dx
                        GR_MIXT_UY   =  5, &   ! du/dy
                        GR_MIXT_VY   =  6, &   ! dv/dy
                        GR_MIXT_WY   =  7, &   ! dw/dy
                        GR_MIXT_TY   =  8, &   ! dT/dy
                        GR_MIXT_UZ   =  9, &   ! du/dz
                        GR_MIXT_VZ   = 10, &   ! dv/dz
                        GR_MIXT_WZ   = 11, &   ! dw/dz
                        GR_MIXT_TZ   = 12      ! dT/dz
#endif
#ifdef RFLU
  INTEGER, PARAMETER :: GRC_MIXT_DENS =  1, &   ! dr/dx, dr/dy, or dr/dz
                        GRC_MIXT_XVEL =  2, &   ! du/dx, du/dy, or du/dz  
                        GRC_MIXT_YVEL =  3, &   ! dv/dx, dv/dy, or dv/dz 
                        GRC_MIXT_ZVEL =  4, &   ! dw/dx, dw/dy, or dw/dz 
                        GRC_MIXT_PRES =  5, &   ! dP/dx, dP/dy, or dP/dz
                        GRC_MIXT_TEMP =  5      ! dT/dx, dT/dy, or dT/dz

  INTEGER, PARAMETER :: GRBF_MIXT_DENS =  1, &   ! dr/dx, dr/dy, or dr/dz
                        GRBF_MIXT_XVEL =  2, &   ! du/dx, du/dy, or du/dz  
                        GRBF_MIXT_YVEL =  3, &   ! dv/dx, dv/dy, or dv/dz 
                        GRBF_MIXT_ZVEL =  4, &   ! dw/dx, dw/dy, or dw/dz 
                        GRBF_MIXT_PRES =  5, &   ! dP/dx, dP/dy, or dP/dz
                        GRBF_MIXT_TEMP =  5      ! dT/dx, dT/dy, or dT/dz

  INTEGER, PARAMETER :: GRF_MIXT_XVEL =  1, &   ! du/dx, du/dy, or du/dz  
                        GRF_MIXT_YVEL =  2, &   ! dv/dx, dv/dy, or dv/dz 
                        GRF_MIXT_ZVEL =  3, &   ! dw/dx, dw/dy, or dw/dz 
                        GRF_MIXT_TEMP =  4      ! dT/dx, dT/dy, or dT/dz
                        
  INTEGER, PARAMETER :: BF2BG_BEG = 1, & 
                        BF2BG_END = 2                         
#endif

! ******************************************************************************
! Species
! ******************************************************************************

  INTEGER, PARAMETER :: SPEC_SOURCE_TYPE_NONE = 0, & 
                        SPEC_SOURCE_TYPE_CHEM = 1, & 
                        SPEC_SOURCE_TYPE_CAVI = 2

  INTEGER, PARAMETER :: SPEC_METHV_FLUIDVEL  = 0, & ! species vel = fluid vel
                        SPEC_METHV_EQEUL     = 1    ! Eq Eul method

  INTEGER, PARAMETER :: SD_XMOM = 1, &
                        SD_YMOM = 2, &
                        SD_ZMOM = 3

! ******************************************************************************
! Grid motion 
! ******************************************************************************

  INTEGER, PARAMETER :: MOVEGRID_BLOCKS   = 0, &
                        MOVEGRID_GLOBAL   = 1, &
                        MOVEGRID_FRAME    = 2, &
                        MOVEGRID_FOMS     = 3, &
                        MOVEGRID_ELGLOBAL = 4, &
                        MOVEGRID_ELFRAME  = 5, &
                        MOVEGRID_VMS      = 6

  INTEGER, PARAMETER :: PATCH_MOVEMENT_OFF = 0, & 
                        PATCH_MOVEMENT_ON  = 1

  INTEGER, PARAMETER :: PATCH_SMOOTHING_OFF = 0, & 
                        PATCH_SMOOTHING_ON  = 1

  INTEGER, PARAMETER :: MOVEGRID_TYPE_DISP = 1, &
                        MOVEGRID_TYPE_XYZ  = 2, & 
                        MOVEGRID_TYPE_GENX = 3      

  INTEGER, PARAMETER :: MOVEGRID_CONTEXT_MOVESMOOTH = 1, & 
                        MOVEGRID_CONTEXT_ONLYSMOOTH = 2

  INTEGER, PARAMETER :: MOVEGRID_BCTYPE_NONE    = -1,      & 
                        MOVEGRID_BCTYPE_NEUMANN =  0,      & 
                        MOVEGRID_BCTYPE_DIRICHX =  XCOORD, & 
                        MOVEGRID_BCTYPE_DIRICHY =  YCOORD, & 
                        MOVEGRID_BCTYPE_DIRICHZ =  ZCOORD
                        
  INTEGER, PARAMETER :: MOVEPATCH_DIR_NONE = 0, & 
                        MOVEPATCH_DIR_X    = 1, & 
                        MOVEPATCH_DIR_Y    = 2, & 
                        MOVEPATCH_DIR_Z    = 4, & 
                        MOVEPATCH_DIR_XY   = 3, & 
                        MOVEPATCH_DIR_XZ   = 5, & 
                        MOVEPATCH_DIR_YZ   = 6, & 
                        MOVEPATCH_DIR_XYZ  = 7 

#ifdef STATS
! ******************************************************************************
! Statistics 
! ******************************************************************************

  INTEGER, PARAMETER :: STAT_NONE    =  0, &   ! statistics variable code
                        STAT_CV      =  1, &   
                        STAT_DV      =  2, &   
                        STAT_TV      =  3, &   
                        STAT_GV      =  4, &
                        STAT_SV      =  5, &
                        STAT_ST      =  6, &
                        STAT_PLAGEV  =  7, &
                        STAT_MAXTYPE =  9

  INTEGER, PARAMETER :: NSTATS_TEC_MIXT = 9, &  !  mixture nStats to tecplot
                        NSTATS_TEC_TURB = 2      !  turb. nStats to tecplot

#endif

  INTEGER, PARAMETER :: FTYPE_MIXT  =  1, &    ! module labels
                        FTYPE_TURB  =  2, &   
                        FTYPE_PLAG  =  3, &   
                        FTYPE_PEUL  =  4, &   
                        FTYPE_SPEC  =  5, &
                        FTYPE_RADI  =  6, &
                        FTYPE_MAX   =  6       ! highest module label

! ******************************************************************************
! Boundary conditions
! ******************************************************************************

  INTEGER, PARAMETER :: BC_NONREFLECTING     =   0, &
                        BC_REFLECTING        =   1

  INTEGER, PARAMETER :: BC_KIND_SIMPLE       =   0, &
                        BC_KIND_NSCBC        =   1, &
                        BC_KIND_MIN          =   0, &
                        BC_KIND_MAX          =   1

  INTEGER, PARAMETER :: BC_RANGE             =   9, &
                        BC_INFLOW            =  10, &
                        BC_INFLOW_TOTANG     =  10, & 
                        BC_INFLOW_VELTEMP    =  11, & 
                        BC_INFLOW_VELPRESS   =  12, & 
                        BC_OUTFLOW           =  20, &
                        BC_OUTFLOW_FREE      =  21, &
                        BC_OUTFLOW_FIXED     =  22, &
                        BC_OUTFLOW_XSLIDE    =  23, &
                        BC_OUTFLOW_YSLIDE    =  24, &
                        BC_OUTFLOW_ZSLIDE    =  25, &
                        BC_OUTFLOW_XYSLIDE   =  26, &
                        BC_OUTFLOW_XZSLIDE   =  27, &
                        BC_OUTFLOW_YZSLIDE   =  28, &
                        BC_REGIONCONF        =  30, &
                        BC_REGIONINT         =  40, &
                        BC_REGNONCONF        =  50, &
                        BC_SLIPWALL          =  60, &
                        BC_SLIPWALL_FREE     =  61, &
                        BC_SLIPWALL_FIXED    =  62, &
                        BC_SLIPWALL_XSLIDE   =  63, &
                        BC_SLIPWALL_YSLIDE   =  64, &
                        BC_SLIPWALL_ZSLIDE   =  65, &
                        BC_SLIPWALL_XYSLIDE  =  66, &
                        BC_SLIPWALL_XZSLIDE  =  67, &
                        BC_SLIPWALL_YZSLIDE  =  68, &
                        BC_NOSLIPWALL        =  70, &
                        BC_NOSLIPWALL_HFLUX  =  70, & 
                        BC_NOSLIPWALL_TEMP   =  71, & 
                        BC_NOSLIPWALL_FREE   =  71, &
                        BC_NOSLIPWALL_FIXED  =  72, &
                        BC_NOSLIPWALL_XSLIDE =  73, &
                        BC_NOSLIPWALL_YSLIDE =  74, &
                        BC_NOSLIPWALL_ZSLIDE =  75, &
                        BC_NOSLIPWALL_XYSLIDE=  76, &
                        BC_NOSLIPWALL_XZSLIDE=  77, &
                        BC_NOSLIPWALL_YZSLIDE=  78, &
                        BC_FARFIELD          =  80, &
                        BC_INJECTION         =  90, &
                        BC_INJECTION_MRATE   =  90, &
                        BC_INJECTION_APN     =  91, &
                        BC_INJECTION_HT      =  92, &
                        BC_SYMMETRY          = 100, &
                        BC_SYMMETRY_FREE     = 101, &
                        BC_SYMMETRY_FIXED    = 102, &
                        BC_SYMMETRY_XSLIDE   = 103, &
                        BC_SYMMETRY_YSLIDE   = 104, &
                        BC_SYMMETRY_ZSLIDE   = 105, &
                        BC_SYMMETRY_XYSLIDE  = 106, &
                        BC_SYMMETRY_XZSLIDE  = 107, &
                        BC_SYMMETRY_YZSLIDE  = 108, &
                        BC_PERIODIC          = 110, & 
                        BC_TRA_PERI          = 110, &
                        BC_ROT_PERI          = 120, &
                        BC_VIRTUAL           = 130, & 
                        BC_CODE_MIN          =  10, &  ! min. and max. BC number
                        BC_CODE_MAX          = 129, & 
                        BC_SIMPLE_COPY       = 999, &
                        BC_INTERNAL          =   0, &  ! source of BC values
                        BC_EXTERNAL          =   1

  INTEGER, PARAMETER :: BC_NOT_BURNING = 0, & 
                        BC_BURNING     = 1, & 
                        BC_NOT_COUPLED = 2

  INTEGER, PARAMETER :: BCDAT_INFLOW_PTOT  = 1, &
                        BCDAT_INFLOW_TTOT  = 2, &
                        BCDAT_INFLOW_BETAH = 3, &
                        BCDAT_INFLOW_BETAV = 4, &
                        BCDAT_INFLOW_MACH  = 5, &
                        BCDAT_INFLOW_U     = 1, & 
                        BCDAT_INFLOW_V     = 2, & 
                        BCDAT_INFLOW_W     = 3, &
                        BCDAT_INFLOW_T     = 4, &
                        BCDAT_INFLOW_P     = 5, & 
                        BCDAT_INFLOW_NELM  = 5, & 
                        BCSWI_INFLOW_TYPE  = 1, & 
                        BCSWI_INFLOW_FIXED = 2, &
                        BCSWI_INFLOW_MODEL = 3

  INTEGER, PARAMETER :: BCDAT_OUTFLOW_PRESS  = 1, &
                        BCDAT_OUTFLOW_NRCOEF = 2, &
                        BCSWI_OUTFLOW_TYPE   = 1, &
                        BCSWI_OUTFLOW_MODEL  = 2

  INTEGER, PARAMETER :: BCSWI_SLIPW_EXTRAP   = 1

  INTEGER, PARAMETER :: BCDAT_NOSLIP_Q       = 1, &
                        BCDAT_NOSLIP_T       = 1, &  
                        BCDAT_NOSLIP_TWALL   = 1, &
                        BCSWI_NOSLIP_ADIABAT = 1

  INTEGER, PARAMETER :: BCDAT_FARF_MACH      = 1, &
                        BCDAT_FARF_ATTACK    = 2, &
                        BCDAT_FARF_SLIP      = 3, &
                        BCDAT_FARF_PRESS     = 4, &
                        BCDAT_FARF_TEMP      = 5, & 
                        BCSWI_FARF_CORR      = 1

  INTEGER, PARAMETER :: BCDAT_INJECT_MFRATE  = 1, &
                        BCDAT_INJECT_TEMP    = 2, &
                        BCDAT_INJECT_RFVFU   = 3, &
                        BCDAT_INJECT_RFVFV   = 4, &
                        BCDAT_INJECT_RFVFW   = 5, &
                        BCDAT_INJECT_SDENS   = 6, &
                        BCDAT_INJECT_ACOEFF  = 7, &
                        BCDAT_INJECT_NPOWER  = 8, &
                        BCSWI_INJECT_EXTRAP  = 1

  INTEGER, PARAMETER :: BCDAT_CONSTANT = 0, &   ! Values do matter!!!
                        BCDAT_DISTRIB  = 1      ! Used to access vals(:,i) !!!

  INTEGER, PARAMETER :: BCOPT_ADIABAT     = 1, &
                        BCOPT_NON_ADIABAT = 0, &
                        BCOPT_SUBSONIC    = 1, &
                        BCOPT_SUPERSONIC  = 0, &
                        BCOPT_MIXED       = 2, & 
                        BCOPT_FIXED_NO    = 0, & 
                        BCOPT_FIXED_YES   = 1, & 
                        BCOPT_CORR_NO     = 0, & 
                        BCOPT_CORR_YES    = 1, &
                        BCOPT_DEFAULT     = 0, &
                        BCOPT_MODEL1      = 1, &
                        BCOPT_STEADY      = 0, &
                        BCOPT_UNSTEADY    = 1

  INTEGER, PARAMETER :: EXTRAPOL_CONST  = 0, &   ! extrapolation to dummy cells
                        EXTRAPOL_LINEAR = 1

  INTEGER, PARAMETER :: NIJK_INFLOW_INIT = 10000 ! initial intg in recycturb
                        
! ******************************************************************************
! Time-dependent boundary conditions 
! ******************************************************************************

! ==============================================================================
! Aliases for choices
! ==============================================================================

  INTEGER, PARAMETER :: TBC_NONE       = 0, &
                        TBC_SINUSOIDAL = 1, &
                        TBC_STOCHASTIC = 2, &
                        TBC_WHITENOISE = 3, &
                        TBC_PIECEWISE  = 4

! ==============================================================================
! Aliases for indices: params and switches
! ==============================================================================

  INTEGER, PARAMETER :: TBCDAT_ONTIME     = 1, & ! All TBC types
                        TBCDAT_OFFTIME    = 2, &
                        TBCDAT_AMP        = 3    ! (AMP not used for PIECEWISE)

  INTEGER, PARAMETER :: TBCDAT_FREQ       = 4, & ! TBC_SINUSOIDAL
                        TBCDAT_PHASE      = 5

  INTEGER, PARAMETER :: TBCDAT_TIMECOR    = 4, & ! TBC_STOCHASTIC
                        TBCDAT_SHAPE      = 5, &
                        TBCDAT_MINCUT     = 6, &
                        TBCDAT_MAXCUT     = 7

  INTEGER, PARAMETER :: TBCSWI_SUBSTEP    = 1    ! TBC_WHITENOISE

  INTEGER, PARAMETER :: TBCOPT_STEP       = 0, & ! TBCSWI_SUBSTEP
                        TBCOPT_SUBSTEP    = 1


  INTEGER, PARAMETER :: TBCSWI_ORDER      = 1, & ! TBC_PIECEWISE
                        TBCSWI_NJUMPS     = 2, &
                        TBCDAT_DAT0       = 2

  INTEGER, PARAMETER :: TBCOPT_CONSTANT   = 0, & ! TBCSWI_ORDER
                        TBCOPT_LINEAR     = 1

! ==============================================================================
! Aliases for indices: svals and bvals
! ==============================================================================

  INTEGER, PARAMETER :: TBCSTO_VAL        = 1, & ! All TBC types
                        TBCSTO_DVAL       = 2, & ! TBC_STOCHASTIC
                        TBCSTO_FACTOR     = 3

! ******************************************************************************
! Constraints
! ******************************************************************************

  INTEGER, PARAMETER :: CONSTR_NONE     = 0, & 
                        CONSTR_WEIGHTED = 1

  INTEGER, PARAMETER :: CONSTR_TYPE_NONE       = 0, & 
                        CONSTR_TYPE_DIRICHLET  = 1, & 
                        CONSTR_TYPE_VONNEUMANN = 2, & 
                        CONSTR_TYPE_ROBIN      = 3

  INTEGER, PARAMETER :: V_MIXT_XVEL = 1, & 
                        V_MIXT_YVEL = 2, & 
                        V_MIXT_ZVEL = 3, &
                        V_MIXT_DENS = 4, & 
                        V_MIXT_PRES = 5, &  
                        V_MIXT_TEMP = 6
                        
  INTEGER, PARAMETER :: V_SPEC_VAR1 = 10, &
                        V_SPEC_VAR2 = 11, &
                        V_SPEC_VAR3 = 12, &
                        V_SPEC_VAR4 = 13, &
                        V_SPEC_VAR5 = 14, &
                        V_SPEC_VAR6 = 15, &
                        V_SPEC_VAR7 = 16, &
                        V_SPEC_VAR8 = 17, &
                        V_SPEC_VAR9 = 18 

! ******************************************************************************
! Time-stepping schemes 
! ******************************************************************************

  INTEGER, PARAMETER :: FLOW_STEADY      = 0, &
                        FLOW_UNSTEADY    = 1, &
                        SOLV_EXPLICIT    = 1, &
                        SOLV_IMPLICIT    = 2, &
                        SOLV_IMPLICIT_NK = 2, & 
                        TST_HYB5RK       = 1, & ! 5-stage hybrid RK scheme
                        TST_STD4RK       = 2    ! 4-stage classical RK scheme

  INTEGER, PARAMETER :: RK_SCHEME_4_CLASSICAL = 1, & 
                        RK_SCHEME_3_WRAY      = 2
                        
  INTEGER, PARAMETER :: VAR_TYPE_CELL  = 1, & 
                        VAR_TYPE_POINT = 2                                 

! ******************************************************************************
! Physical models
! ******************************************************************************

  INTEGER, PARAMETER :: FLUID_MODEL_INCOMP = 0, & 
                        FLUID_MODEL_COMP   = 1

  INTEGER, PARAMETER :: FLOW_EULER = 0, &
                        FLOW_NAVST = 1

  INTEGER, PARAMETER :: TURB_MODEL_NONE = 0
  
  INTEGER, PARAMETER :: GAS_MODEL_TCPERF      = 1, & 
                        GAS_MODEL_TPERF       = 2, & 
                        GAS_MODEL_MIXT_TCPERF = 3, & 
                        GAS_MODEL_MIXT_TPERF  = 4, & 
                        GAS_MODEL_MIXT_PSEUDO = 5, &
                        GAS_MODEL_MIXT_GASLIQ = 6  

  INTEGER, PARAMETER :: VISC_SUTHR = 0, &
                        VISC_FIXED = 1, &
                        VISC_ANTIB = 2

! ******************************************************************************
! Numerics 
! ******************************************************************************

  INTEGER, PARAMETER :: DISCR_CEN_SCAL     =  0, &
                        DISCR_UPW_ROE      =  1, &
                        DISCR_UPW_MAPS     =  2, & 
                        DISCR_UPW_HLLC     =  3, &
                        DISCR_UPW_AUSMPLUS =  4, & 
                        DISCR_OPT_LES      = 99

  INTEGER, PARAMETER :: DISCR_ORDER_1 = 1, &
                        DISCR_ORDER_2 = 2, &
                        DISCR_ORDER_4 = 4

  INTEGER, PARAMETER :: FLUX_PART_CENTRAL = 0, & 
                        FLUX_PART_DISSIP  = 1, & 
                        FLUX_PART_BOTH    = 2

  INTEGER, PARAMETER :: RECONST_NONE          =  0, &
                        RECONST_WENO_SIMPLE   =  1, & 
                        RECONST_WENO_XYZ      =  2, & 
                        RECONST_LIM_BARTHJESP = 10, & 
                        RECONST_LIM_VENKAT    = 11

  INTEGER, PARAMETER :: PSWITCH_STD = 0, &
                        PSWITCH_TVD = 1

  INTEGER, PARAMETER :: MGCYCLE_NO = 0, &
                        MGCYCLE_V  = 1, &
                        MGCYCLE_W  = 2
                        
  INTEGER, PARAMETER :: FE_AVG_UNIFORM = 0, &
                        FE_AVG_LINEAR  = 1    

! ******************************************************************************
! Probes
! ******************************************************************************

  INTEGER, PARAMETER :: PROBE_REGION = 1, & 
                        PROBE_CELL   = 2, & 
                        PROBE_ILOC   = 2, & 
                        PROBE_JLOC   = 3, & 
                        PROBE_KLOC   = 4 

! ******************************************************************************
! Forces calculation 
! ******************************************************************************

  INTEGER, PARAMETER :: FORCES_NONE  = 0, &
                        FORCES_PRESS = 1, &
                        FORCES_VISC  = 2

  INTEGER, PARAMETER :: COMP_MOM  = 1, &
                        COMP_PRES = 2, &
                        COMP_VISC = 3

! ******************************************************************************
! Thrust calculation 
! ******************************************************************************

  INTEGER, PARAMETER :: THRUST_NONE = 0, &
                        THRUST_MOM  = 1, &
                        THRUST_MOMP = 2

  INTEGER, PARAMETER :: MASS_IN  = 1, &
                        MASS_OUT = 2 

! ******************************************************************************
! File formats
! ******************************************************************************

  INTEGER, PARAMETER :: FORMAT_ASCII  = 0, &
                        FORMAT_BINARY = 1, &
                        FORMAT_HDF    = 2

! ******************************************************************************
! File destinations 
! ******************************************************************************

  INTEGER, PARAMETER :: FILEDEST_INDIR  = 1, & 
                        FILEDEST_OUTDIR = 2

! ******************************************************************************
! Grid source 
! ******************************************************************************

  INTEGER, PARAMETER :: GRID_SRC_CENTAUR_ASCII  =  0, & 
                        GRID_SRC_VGRIDNS        =  1, & 
                        GRID_SRC_MESH3D         =  2, & 
                        GRID_SRC_TETMESH        =  3, & 
                        GRID_SRC_COBALT         =  4, &
                        GRID_SRC_GAMBIT_NEUTRAL =  5, &
                        GRID_SRC_CENTAUR_BINARY = 10

! ******************************************************************************
! Region activation 
! ******************************************************************************

  INTEGER, PARAMETER :: ACTIVE = 1, &
                        OFF    = 0

! ******************************************************************************
! Verbosity levels 
! ******************************************************************************

  INTEGER, PARAMETER :: VERBOSE_NONE = 0, &
                        VERBOSE_LOW  = 1, &
                        VERBOSE_MED  = 2, &
                        VERBOSE_HIGH = 3

! ******************************************************************************
! Checking level 
! ******************************************************************************
 
  INTEGER, PARAMETER :: CHECK_NONE = 0, & 
                        CHECK_LOW  = 1, & 
                        CHECK_HIGH = 2

! ******************************************************************************
! File IDs
! ******************************************************************************

  INTEGER, PARAMETER :: IF_INPUT      = 10, &
                        IF_GRID       = 11, &
                        IF_TOPOL      = 12, &
                        IF_SOLUT      = 13, &
                        IF_CONVER     = 14, &
                        IF_DISTR      = 15, &
                        IF_PLOT       = 16, &
                        IF_REGMAP     = 17, &
                        IF_CONTROL    = 18, &
                        IF_PTMATRIX   = 19, & 
                        IF_BC         = 20, &  ! BC files
#ifdef STATS
                        IF_STAT       = 29, &
#endif
                        IF_INTEG_OLES = 30, &  ! Optimal LES integrals
                        IF_STATS_OLES = 31, &  ! Optimal LES statistics
                        IF_DIMS       = 32, &  ! Dimensions
                        IF_VRS        = 33, &  ! Version
                        IF_POSTINFO   = 34, &  ! Postprocessor info
                        IF_RESTINFO   = 35, &  ! Restart info
                        IF_PATCH_COEF = 36, &  ! Patch coefficients
                        IF_DEGENRT    = 37, &  ! Degeneracy info 
                        IF_RNMB       = 38, &  ! Renumberings
                        IF_CELL_MAPS  = 39, &  ! Cell mappings 
                        IF_MASS       = 40, &  ! Total mass
                        IF_THRUST     = 41, &  ! thrust data
                        IF_RAND_STATE = 42, &  ! State of Random Number Gen
                        IF_FORMOM     = 43, &  ! Forces and moments
                        IF_COMM_LISTS = 44, &  ! Communication lists
                        IF_COLOR      = 45, &  ! Coloring
                        IF_CTRL_VOL   = 46, & 
                        IF_CTRL_SURF  = 47, &
                        IF_EXTR_DATA1 = 48, & 
                        IF_EXTR_DATA2 = 49, & 
                        IF_PROBE      = 50, &  ! + global%nProbes-1 channels
                        IF_ISP        = 51     ! Specific Impulse
                
#ifdef PEUL
  INTEGER, PARAMETER :: IF_PEUL_SOLUT = 73
#endif
#ifdef TURB
  INTEGER, PARAMETER :: IF_TURB_SOLUT = 81
#endif
#ifdef PLAG
  INTEGER, PARAMETER :: IF_PLAG_SURF_STATS = 90, & ! Surface statistics
                        IF_PLAG_STATS      = 91, & ! Statistics
                        IF_PLAG_INJCPDF    = 92    ! Pdf for injection
#endif

  INTEGER, PARAMETER :: IF_ENS_CASE     = 100, & 
                        IF_ENS_GEOMETRY = 101, & 
                        IF_ENS_SCALAR   = 102, & 
                        IF_ENS_VECTOR   = 200

  INTEGER, PARAMETER :: STDIN  = 5, &
                        STDOUT = 6, &
                        STDERR = 6
                        
  INTEGER, PARAMETER :: FILE_STATUS_OLD     = 1, & 
                        FILE_STATUS_UNKNOWN = 2
                        
  INTEGER, PARAMETER :: FILE_POSITION_START = 1, & 
                        FILE_POSITION_END   = 2
                        
! ******************************************************************************
! Output
! ******************************************************************************

  INTEGER, PARAMETER :: POST_OUTPUT_FORMAT_TECPLOT = 1, & 
                        POST_OUTPUT_FORMAT_ENSIGHT = 2

  INTEGER, PARAMETER :: PLOT_GRID_ONLY = 1, & 
                        PLOT_GRID_FLOW = 2

  INTEGER, PARAMETER :: PLOT_FMT_GENERIC  = 1, &
                        PLOT_FMT_TECPLOT  = 2, &
                        PLOT_FMT_TECASCII = 3

! ******************************************************************************
! PETSc
! ****************************************************************************** 

  INTEGER, PARAMETER :: RFLU_PETSC_POISSON_INFO_BEG = 1, & 
                        RFLU_PETSC_POISSON_INFO_A   = 1, & 
                        RFLU_PETSC_POISSON_INFO_B   = 2, & 
                        RFLU_PETSC_POISSON_INFO_X   = 3, &
                        RFLU_PETSC_POISSON_INFO_KSP = 4, &
                        RFLU_PETSC_POISSON_INFO_PC  = 5, &
                        RFLU_PETSC_POISSON_INFO_NSP = 6, &  
                        RFLU_PETSC_POISSON_INFO_END = 6

! ******************************************************************************
! Miscellaneous 
! ******************************************************************************

  INTEGER, PARAMETER :: CRAZY_VALUE_INT = -987654321 ! MUST be negative

  INTEGER, PARAMETER :: MIN_VAL = 1, & 
                        MAX_VAL = 2               
                        
  INTEGER, PARAMETER :: INITFLOW_FROMSCRATCH        = 1, & 
                        INITFLOW_FROMFILE           = 2, & 
                        INITFLOW_FROMHARDCODE       = 3, &
                        INITFLOW_FROMCOMBO_SERIAL   = 4, & 
                        INITFLOW_FROMCOMBO_PARALLEL = 5

  INTEGER, PARAMETER :: PLAG_INIT_FROMSCRATCH     = 1, & ! Initialization Flags for Plag
                        PLAG_INIT_FROMFILE        = 2, & 
                        PLAG_INIT_FROMHARDCODE    = 3, &
                        PLAG_INIT_FROMRANDOMSTATE = 4

  INTEGER, PARAMETER :: PATCH_DIMENS_NPATCHMAX   = 100, & 
                        PATCH_DIMENS_BEG         = 1, &
                        PATCH_DIMENS_IPGLOBAL    = 1, &  
                        PATCH_DIMENS_NBTRIS      = 2, &
                        PATCH_DIMENS_NBTRISTOT   = 3, &
                        PATCH_DIMENS_NBTRISMAX   = 4, & 
                        PATCH_DIMENS_NBQUADS     = 5, & 
                        PATCH_DIMENS_NBQUADSTOT  = 6, &
                        PATCH_DIMENS_NBQUADSMAX  = 7, &
                        PATCH_DIMENS_NBCELLSVIRT = 8, & 
                        PATCH_DIMENS_END         = 8
                        
  INTEGER, PARAMETER :: PATCH_IBORDER_DEFAULT   = 0 ! Must be zero or negative                        
                        
  INTEGER, PARAMETER :: BORDER_INFO_MAX    = 100, & 
                        BORDER_INFO_BEG    = 1, &
                        BORDER_INFO_IRGLOB = 1, & 
                        BORDER_INFO_IBORD  = 2, & 
                        BORDER_INFO_NCSEND = 3, & 
                        BORDER_INFO_NCRECV = 4, &  
                        BORDER_INFO_NVSEND = 5, & 
                        BORDER_INFO_NVRECV = 6, & 
                        BORDER_INFO_NVSHAR = 7, &
                        BORDER_INFO_END    = 7                         
                        
  INTEGER, PARAMETER :: LOCINFO_MODE_SILENT  = 0, & 
                        LOCINFO_MODE_VERBOSE = 1
                        
  INTEGER, PARAMETER :: OUTPUT_MODE_MASTER_ONLY = 0, & 
                        OUTPUT_MODE_ANYBODY     = 1                                             

  INTEGER, PARAMETER :: INFOFILE_READMODE_FLAG = 0, & 
                        INFOFILE_READMODE_DATA = 1

  INTEGER, PARAMETER :: MAPFILE_READMODE_ALL  = 0, &
                        MAPFILE_READMODE_PEEK = 1                                                
  
  INTEGER, PARAMETER :: ALLOC_MODE_ALL       = 0, & 
                        ALLOC_MODE_PRIMARY   = 1, & 
                        ALLOC_MODE_SECONDARY = 2, & 
                        ALLOC_MODE_GRIDGEOM  = 3, & 
                        ALLOC_MODE_FLOWSOL   = 4                                             

  INTEGER, PARAMETER :: NCELLS_SPECIAL_MAX = 100, & 
                        NFACES_SPECIAL_MAX = 100

  INTEGER, PARAMETER :: COMPWTS_MODE_FIXED = 0, & 
                        COMPWTS_MODE_ADAPT = 1
                         
  INTEGER, PARAMETER :: COMPWTS_SCAL_NONE    = 0, & 
                        COMPWTS_SCAL_INVDIST = 1

  INTEGER, PARAMETER :: MODULE_TYPE_NONE     = 0, & 
                        MODULE_TYPE_PART     = 1, & 
                        MODULE_TYPE_SOLVER   = 2, & 
                        MODULE_TYPE_POSTPROC = 3, & 
                        MODULE_TYPE_INIT     = 4 

#ifdef RFLO
  INTEGER, PARAMETER :: REGOFF = 100 ! offset for region number
#endif

! ******************************************************************************
! Preprocessing 
! ******************************************************************************

  INTEGER, PARAMETER :: PARTITION_MODE_PROPER  = 1, & 
                        PARTITION_MODE_IMPOSED = 2
  INTEGER, PARAMETER :: WRITE_GRID_ON  = 1, & 
                        WRITE_GRID_OFF = 2 

! ******************************************************************************
! Postprocessing 
! ******************************************************************************

  INTEGER, PARAMETER :: INTERP_TYPE_NONE   = 0, &  
                        INTERP_TYPE_SIMPLE = 1, & 
                        INTERP_TYPE_PROPER = 2 

! ******************************************************************************
! Random number generator
! ******************************************************************************

  INTEGER, PARAMETER :: RAND_SEED_TYPE_FIXED = 0, &
                        RAND_SEED_TYPE_CLOCK = 1

! ******************************************************************************
! Output prefix
! ******************************************************************************

#ifdef GENX
#ifdef RFLO
  CHARACTER(7), PARAMETER :: SOLVER_NAME = 'Rocflo:'
#endif 
#ifdef RFLU
  CHARACTER(7), PARAMETER :: SOLVER_NAME = 'Rocflu:'
#endif
#else
#ifdef RFLO
  CHARACTER(0), PARAMETER :: SOLVER_NAME = ''
#endif 
#ifdef RFLU
  CHARACTER(0), PARAMETER :: SOLVER_NAME = ''
#endif
#endif

END MODULE ModParameters

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModParameters.F90,v $
! Revision 1.165  2010/02/18 21:47:40  juzhang
! Heat transfer bc for non-propellant surface documented in Rocburn_PY_HT.pdf in Rocburn_PY directory is implemented within Rocburn_PY. Major changes were made to Rocburn, Rocman3, RocfluidMP/genx, RocfluidMP/modflo directories.
!
! Revision 1.164  2009/08/27 14:04:52  mtcampbe
! Updated to enable burning motion with symmetry boundaries and enhanced
! burnout code.
!
! Revision 1.163  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.162  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.161  2008/10/23 18:20:56  mtcampbe
! Crazy number of changes to track and fix initialization and
! restart bugs.  Many improperly formed logical expressions
! were fixed, and bug in allocation for data associated with
! the BC_INFLOWVELTEMP boundary condition squashed in
! RFLO_ReadBcInflowVelSection.F90.
!
! Revision 1.160  2007/03/19 21:41:11  haselbac
! Renamed PV_MIXT_NVAR to PV_XXXX_NVAR
!
! Revision 1.159  2006/10/20 21:30:47  mparmar
! Added parameters for thrust and Isp computations
!
! Revision 1.158  2006/08/28 11:44:48  rfiedler
! Add grid motion constraint types for outflow BC.
!
! Revision 1.157  2006/08/24 13:15:36  rfiedler
! Rocflo now supports XYSLIDE, XZSLIDE, and YZSLIDE instead of TANGEN constraint.
!
! Revision 1.156  2006/08/19 15:54:22  mparmar
! Added parameters for NSCBC implementation
!
! Revision 1.155  2006/08/18 14:00:31  haselbac
! Changed unit numbers to fit in IF_PTMATRIX
!
! Revision 1.154  2006/05/08 22:31:34  wasistho
! added rocprop BC 71-76
!
! Revision 1.153  2006/05/01 20:59:50  haselbac
! Added FLUX_PART parameters
!
! Revision 1.152  2006/04/15 16:57:43  haselbac
! Changed parameters for RECONST_
!
! Revision 1.151  2006/04/07 15:19:18  haselbac
! Removed tabs
!
! Revision 1.150  2006/03/30 20:49:00  haselbac
! Added parameters for cavitation source term
!
! Revision 1.149  2006/03/26 20:21:59  haselbac
! Removed FLUID_MODEL_GASLIQ, added GAS_MODEL_MIXT_GASLIQ
!
! Revision 1.148  2006/03/25 21:47:22  haselbac
! Added PATCH_DIMENS_ parameters bcos of sype patch changes
!
! Revision 1.147  2006/03/04 04:33:15  wasistho
! added REGION_SHAPE_...
!
! Revision 1.146  2006/03/02 01:26:25  wasistho
! split movegrid_epde to elglobal and elframe
!
! Revision 1.145  2006/01/20 06:16:00  wasistho
! added sdens, acoeff and npower in bcdat_inject
!
! Revision 1.144  2006/01/06 22:09:38  haselbac
! Added PV params for gradient errors
!
! Revision 1.143  2005/12/01 17:12:09  fnajjar
! Added parameter defs for randSeedType
!
! Revision 1.142  2005/12/01 08:57:37  wasistho
! added NIJK_INFLOW_INIT
!
! Revision 1.141  2005/11/30 22:16:35  fnajjar
! Added params for PV_PLAG bcos of Eulerian vars in rflupost
!
! Revision 1.140  2005/11/28 20:03:32  wasistho
! added movegrid_epde
!
! Revision 1.139  2005/11/14 16:57:01  haselbac
! Added param for pseudo-gas model
!
! Revision 1.138  2005/11/10 02:21:47  haselbac
! Changed parameters for gas model
!
! Revision 1.137  2005/10/31 21:09:35  haselbac
! Changed specModel and SPEC_MODEL_NONE
!
! Revision 1.136  2005/10/31 19:27:34  haselbac
! Added gas model parameters
!
! Revision 1.135  2005/10/27 18:59:07  haselbac
! Added parameter for constraints
!
! Revision 1.134  2005/10/27 05:11:36  wasistho
! swab moveGridVms and moveGridFoms
!
! Revision 1.133  2005/10/17 22:32:01  wasistho
! added movegrid_foms
!
! Revision 1.132  2005/10/05 20:03:30  haselbac
! Added params for ENSIGHT
!
! Revision 1.131  2005/10/05 13:53:11  haselbac
! Added new params for bc and constrained reconstruction
!
! Revision 1.130  2005/09/20 23:16:43  wasistho
! added BCSWI_INFLOW_MODEL and BCDAT_INFLOW_NELM
!
! Revision 1.129  2005/09/13 21:37:07  haselbac
! Added new init option
!
! Revision 1.128  2005/08/19 02:32:55  haselbac
! Added IF_COLOR
!
! Revision 1.127  2005/08/10 00:34:42  haselbac
! Added more PV_MIXT parameters
!
! Revision 1.126  2005/08/03 18:21:19  hdewey2
! Added SOLV_IMPLICIT_NK parameter
!
! Revision 1.125  2005/07/25 12:22:21  haselbac
! Added PV_MIXT parameters
!
! Revision 1.124  2005/07/14 21:41:30  haselbac
! Added parameter for AUSM flux function
!
! Revision 1.123  2005/07/11 19:27:38  mparmar
! Added reconst parameters
!
! Revision 1.122  2005/06/29 22:51:38  wasistho
! added EDGE_INTERACT parameters
!
! Revision 1.121  2005/06/19 05:31:30  wasistho
! shift index rocprop slipwalls
!
! Revision 1.120  2005/06/13 02:18:43  wasistho
! added new bc_slipwall types including xyz slidewalls
!
! Revision 1.119  2005/06/09 20:18:46  haselbac
! Added MOVEPATCH_DIR parameters
!
! Revision 1.118  2005/06/02 03:21:03  wasistho
! shuffle MoveGridVms with MoveGridFrame
!
! Revision 1.117  2005/05/28 08:06:14  wasistho
! added parameters for moveGridFrame
!
! Revision 1.116  2005/05/21 00:17:35  wasistho
! added MOVE_GRID_VMS
!
! Revision 1.115  2005/05/13 06:06:37  wasistho
! added parameters for bc outflow
!
! Revision 1.114  2005/04/28 05:46:27  wasistho
! added BC_INFLOW_VELPRESS
!
! Revision 1.113  2005/04/27 02:09:32  haselbac
! Added parameters for inflow bc based on velocities and temperature
!
! Revision 1.112  2005/04/25 18:39:08  luca1
! Imposed PDF from file option for random particle ejection
!
! Revision 1.111  2005/04/15 15:06:32  haselbac
! Added MPI parameters and modified module type parameters
!
! Revision 1.110  2005/03/31 16:58:03  haselbac
! Changed SD access parameters
!
! Revision 1.109  2005/03/29 22:29:37  haselbac
! Added INITFLOW_FROMCOMBO
!
! Revision 1.108  2005/03/09 14:55:00  haselbac
! Added parameter for virtual boundary
!
! Revision 1.107  2005/01/08 20:35:03  fnajjar
! Added IF_PLAG_STATS for PLAG statistics file
!
! Revision 1.106  2004/12/29 23:26:30  wasistho
! prepared statistics for PLAG and PEUL
!
! Revision 1.105  2004/12/29 21:03:33  haselbac
! Added parameters for parallelization
!
! Revision 1.104  2004/12/27 23:27:47  haselbac
! Added parameters for farf bc and interpolation
!
! Revision 1.103  2004/12/21 15:02:09  fnajjar
! Added file definition for PLAG surface statistics
!
! Revision 1.102  2004/12/19 15:46:07  haselbac
! Added PETSC parameters for incompressible solver
!
! Revision 1.101  2004/12/04 03:23:55  haselbac
! Added parameters for vertex kinds and file units
!
! Revision 1.100  2004/11/17 16:29:49  haselbac
! Added parameters for RK scheme and variable type
!
! Revision 1.99  2004/11/09 10:56:05  wasistho
! added statistics parameters due to inclusion statistics in rflopost
!
! Revision 1.98  2004/11/03 14:55:24  haselbac
! Added parameter for GAMBIT grid conversion and added comment
!
! Revision 1.97  2004/11/02 02:29:49  haselbac
! Added FLUID_MODEL_ parameters
!
! Revision 1.96  2004/10/26 15:17:32  haselbac
! Added parameters for data extraction
!
! Revision 1.95  2004/10/19 19:29:03  haselbac
! Clean-up
!
! Revision 1.94  2004/10/10 20:03:38  fnajjar
! Added PLAG initialization flags
!
! Revision 1.93  2004/09/27 01:37:16  haselbac
! Added parameters for special and opposing faces
!
! Revision 1.92  2004/09/02 02:33:46  wasistho
! added face-edge averaging input-option parameter in Rocflo
!
! Revision 1.91  2004/08/21 00:30:31  wasistho
! parameters DEGENERAT_...
!
! Revision 1.90  2004/07/30 22:47:35  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.89  2004/07/21 14:55:34  haselbac
! Added INTER_TYPE_ parameters
!
! Revision 1.88  2004/07/06 15:14:16  haselbac
! Moved WRITE_DIMENS parameters to RFLU_ModDimensions
!
! Revision 1.87  2004/06/30 21:08:58  wasistho
! removed ifdef GENX surrounded REGOFF
!
! Revision 1.86  2004/06/30 04:04:44  wasistho
! moved Genx related parameter REGOFF to ModParameters
!
! Revision 1.85  2004/06/16 20:00:56  haselbac
! Added patch coefficient and moment parameters
!
! Revision 1.84  2004/04/01 21:29:03  haselbac
! Added SPEC_SOURCE_* parameters
!
! Revision 1.83  2004/03/17 04:26:28  haselbac
! Added parameters for writing of Rocflu dimensions file
!
! Revision 1.82  2004/03/06 02:32:57  wasistho
! moved mpi tag shifts from ModParameters to ModMPI
!
! Revision 1.81  2004/03/05 21:10:28  wasistho
! added TURB and PEUL mpi-tag-shift
!
! Revision 1.80  2004/02/26 21:02:00  haselbac
! Added parameters for memory allocation
!
! Revision 1.79  2004/02/07 00:53:23  wasistho
! defined separate file id for turbulence solution
!
! Revision 1.78  2004/01/29 22:57:28  haselbac
! Added various new parameters
!
! Revision 1.77  2003/12/04 03:28:27  haselbac
! Added new parameters for data structures and gradients
!
! Revision 1.76  2003/11/21 22:35:51  fnajjar
! Update Random Number Generator
!
! Revision 1.75  2003/10/03 21:39:03  fnajjar
! Bug fix for IF_PEUL_SOLUT by including INTEGER definition
!
! Revision 1.74  2003/09/29 15:29:51  fnajjar
! Removed ampersand and comma after IF_PROBE since build breaks w/o PEUL
!
! Revision 1.73  2003/09/26 22:51:41  jferry
! added parameter for rocsmoke solution files
!
! Revision 1.72  2003/09/15 00:37:06  haselbac
! Added INITFLOW_FROMHARDCODE
!
! Revision 1.71  2003/08/19 22:46:35  haselbac
! Added parameter for COBALT format and C2V_INIT
!
! Revision 1.70  2003/08/11 21:50:00  jblazek
! Splitted ModInterfaces into 4 sections.
!
! Revision 1.69  2003/07/22 02:01:47  haselbac
! Modified and added diff and interp parameters
!
! Revision 1.68  2003/07/03 21:48:45  jblazek
! Implemented dual-time stepping.
!
! Revision 1.67  2003/06/20 22:33:57  haselbac
! Added three new parameters
!
! Revision 1.66  2003/06/10 22:54:43  jferry
! Added Piecewise TBC
!
! Revision 1.65  2003/06/04 22:07:04  haselbac
! Added and removed some parameters
!
! Revision 1.64  2003/06/02 17:11:32  jblazek
! Added computation of thrust.
!
! Revision 1.63  2003/05/24 02:13:17  wasistho
! turbulence statistics expanded
!
! Revision 1.62  2003/05/16 22:07:08  haselbac
! Added HLLC parameter
!
! Revision 1.61  2003/05/07 00:23:48  haselbac
! Added WRITE_GRID parameters
!
! Revision 1.60  2003/04/28 22:42:48  haselbac
! Added PARTITION_MODE_* parameters
!
! Revision 1.59  2003/04/10 23:22:44  fnajjar
! Added Parameters for viscosity models
!
! Revision 1.58  2003/04/10 18:48:24  haselbac
! Changed grid source parameters
!
! Revision 1.57  2003/04/01 19:38:01  haselbac
! Added file status parameters
!
! Revision 1.56  2003/04/01 16:38:50  haselbac
! Added NCELLS_SPECIAL_MAX
!
! Revision 1.55  2003/03/31 16:14:45  haselbac
! Added IF_VRS and MOVEGRID_TYPEs
!
! Revision 1.54  2003/03/19 16:46:23  haselbac
! Added GRID_SOURCE_TETMESH
!
! Revision 1.53  2003/03/18 21:30:59  haselbac
! Added allocation mode parameters
!
! Revision 1.52  2003/03/17 20:40:44  jblazek
! Changed channel numbers (possible conflict with some GenX crap).
!
! Revision 1.51  2003/03/15 17:50:29  haselbac
! Added new parameters, primarily for || RFLU
!
! Revision 1.50  2003/01/28 16:46:54  haselbac
! Added various new parameters
!
! Revision 1.49  2003/01/23 17:48:53  jblazek
! Changed algorithm to dump convergence, solution and probe data.
!
! Revision 1.48  2002/11/26 15:25:15  haselbac
! Added parameters for boundary condition on grid motion
!
! Revision 1.47  2002/11/08 21:25:36  haselbac
! Added parameter for total mass
!
! Revision 1.46  2002/10/27 19:03:20  haselbac
! Added several parameters for edge list and grid motion
!
! Revision 1.45  2002/10/05 19:00:19  haselbac
! Added probe and GENX parameters (coupling types)
!
! Revision 1.44  2002/09/25 18:29:57  jferry
! simplified TBC parameter lists
!
! Revision 1.43  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.42  2002/09/17 13:43:00  jferry
! Added Time-dependent boundary conditions
!
! Revision 1.41  2002/09/09 14:59:21  haselbac
! Added and changed various parameters
!
! Revision 1.40  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.39  2002/08/30 01:47:58  jblazek
! Added support for moving grids.
!
! Revision 1.38  2002/08/18 02:22:35  wasistho
! Moved TURB parameters into rocturb
!
! Revision 1.37  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
! Revision 1.36  2002/08/01 01:29:31  wasistho
! Added gradient parameters for RFLU
!
! Revision 1.35  2002/07/29 17:10:45  jblazek
! Put TURB stuff into #ifdef.
!
! Revision 1.34  2002/07/27 08:08:24  wasistho
! prepared for rocturb preliminary stage
!
! Revision 1.33  2002/07/25 15:12:30  haselbac
! Added various new parameters for diff and OLES, and others
!
! Revision 1.32  2002/07/25 00:39:01  jblazek
! Option for TVD type pressure switch.
!
! Revision 1.31  2002/07/20 00:43:16  jblazek
! Added ASCII Tecplot format.
!
! Revision 1.30  2002/06/30 00:01:44  jblazek
! Removed TAB characters. Grrrrr ...
!
! Revision 1.29  2002/06/27 15:55:47  haselbac
! Added parameters for parallelization
!
! Revision 1.28  2002/06/22 01:13:37  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.27  2002/06/17 13:40:44  haselbac
! Added SOLVER_NAME character parameter
!
! Revision 1.26  2002/06/14 21:34:32  wasistho
! Added time avg statistics
!
! Revision 1.25  2002/06/14 20:15:28  haselbac
! Added parameters for vertex and cell types (for parallel)
!
! Revision 1.24  2002/06/10 21:27:14  haselbac
! Deleted parameters for inflow
!
! Revision 1.23  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.22  2002/05/21 01:48:22  wasistho
! add viscous terms
!
! Revision 1.21  2002/05/04 16:57:28  haselbac
! Added IF_REGMAP AND BC_SIMPLE_COPY
!
! Revision 1.20  2002/04/11 18:54:24  haselbac
! Added new parameters for primitive cv and miscellaneous
!
! Revision 1.19  2002/03/29 23:15:22  jblazek
! Corrected bug in MPI send.
!
! Revision 1.18  2002/03/27 15:52:39  haselbac
! Changed values of GRID_SOURCE_XXX paramaters for ROCFLU
!
! Revision 1.17  2002/03/26 19:19:55  haselbac
! Added parameters for ROCFLU - grid source
!
! Revision 1.16  2002/03/14 19:07:57  haselbac
! Added parameter FACE_SPLIT_NO
!
! Revision 1.15  2002/03/01 16:49:07  haselbac
! Added parameters for data structure and its generation
!
! Revision 1.14  2002/02/27 18:38:20  jblazek
! Changed extrapol. to dummy cells at injection boundaries and slip walls.
!
! Revision 1.13  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.12  2002/02/16 07:16:00  jblazek
! Added implicit residual smoothing.
!
! Revision 1.11  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.10  2002/02/08 15:07:42  haselbac
! Added parameters for plotting of grid and or flow
!
! Revision 1.9  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.8  2002/01/10 00:02:07  jblazek
! Added calculation of mixture properties.
!
! Revision 1.7  2002/01/02 16:20:19  jblazek
! Added flow initialization and dummy cell geometry.
!
! Revision 1.6  2001/12/22 00:09:38  jblazek
! Added routines to store grid and solution.
!
! Revision 1.5  2001/12/19 23:09:21  jblazek
! Added routines to read grid and solution.
!
! Revision 1.4  2001/12/08 00:18:41  jblazek
! Added routines to read BC input file.
!
! Revision 1.3  2001/12/07 18:36:42  jblazek
! Update of ModError and ModParameters.
!
! Revision 1.2  2001/12/07 16:47:44  jblazek
! ModError and ModParameters updated.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
! ******************************************************************************






