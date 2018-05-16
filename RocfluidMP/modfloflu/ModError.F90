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
! Purpose: define error messages and provide an error handling function.
!
! Description: none
!
! Notes: none
!
! ******************************************************************************
!
! $Id: ModError.F90,v 1.151 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2001-2006 by the University of Illinois
!
! ******************************************************************************

MODULE ModError

  IMPLICIT NONE

! ******************************************************************************
! Error codes
! ******************************************************************************

! ==============================================================================
! Common
! ==============================================================================

  INTEGER, PARAMETER :: ERR_NONE             =   0, &  ! basics
                        ERR_REGISTER_FUN     =   1, &
                        ERR_REACHED_DEFAULT  =   2, &
                        ERR_PREVIOUS_ERRORS  =   3, &
                        ERR_EXTERNAL_FUNCT   =   4, &
                        ERR_EXCEEDS_DECL_MEM =   5, &
                        ERR_UNKNOWN_OPTION   =   6, &
                        ERR_ILLEGAL_VALUE    =   7, &
                        ERR_SYSTEM_COMMAND   =   8, &
                        ERR_COMPILE_OPTION   =   9

  INTEGER, PARAMETER :: ERR_ALLOCATE         =  10, &  ! memory allocation
                        ERR_DEALLOCATE       =  11, &
                        ERR_ASSOCIATED       =  12, &
                        ERR_ALLOCATED        =  13

  INTEGER, PARAMETER :: ERR_DEPENDENT_INPUT  =  15     ! input parameters

  INTEGER, PARAMETER :: ERR_FILE_OPEN        =  20, &  ! I/O
                        ERR_FILE_CLOSE       =  21, &
                        ERR_FILE_READ        =  22, &
                        ERR_FILE_WRITE       =  23, &
                        ERR_FILE_EXIST       =  24, &
                        ERR_UNKNOWN_FORMAT   =  25, &
                        ERR_PROBE_LOCATION   =  26, &
                        ERR_PROBE_SPECIFIED  =  27, &
                        ERR_VAL_UNDEFINED    =  28, &
                        ERR_STOPFILE_FOUND   =  29

  INTEGER, PARAMETER :: ERR_UNKNOWN_BC       =  30, &  ! boundary conditions
                        ERR_NO_BCSPECIFIED   =  31, &
                        ERR_BCVAL_MISSING    =  32, &
                        ERR_NO_BCSWITCH      =  33, &
                        ERR_VAL_BCSWITCH     =  34, &
                        ERR_PATCH_RANGE      =  35, &
                        ERR_BC_VARNAME       =  36, &
                        ERR_VAL_BCVAL        =  37, &
                        ERR_TBC_NONINIT      =  38, &
                        ERR_BCVAL_EXTRA      =  39

  INTEGER, PARAMETER :: ERR_TIME_GRID        =  40, &  ! time stepping
                        ERR_TIME_SOLUTION    =  41, &
                        ERR_TIME_GLOBAL      =  42, &
                        ERR_DTIME_NEGATIVE   =  43, &
                        ERR_DITER_NEGATIVE   =  44, &
                        ERR_STEADY           =  45

  INTEGER, PARAMETER :: ERR_MPI_TROUBLE      =  50, &  ! MPI
                        ERR_NO_PROCMAP       =  51, &
                        ERR_NO_PROCMATCH     =  52, &
                        ERR_PATCH_OFFSET     =  53, &
                        ERR_PROC_MISMATCH    =  54

  INTEGER, PARAMETER :: ERR_GRID_LEVEL       =  60     ! multigrid

  INTEGER, PARAMETER :: ERR_GRAD_INDEX       =  65     ! Vel. and T gradients

  INTEGER, PARAMETER :: ERR_FACEVEC_SUM      =  70, &  ! CV geometry
                        ERR_VOLUME_SIZE      =  71, &
                        ERR_FACE_SPLIT       =  72, &
                        ERR_FACE_INVERTED    =  73

#ifdef STATS
  INTEGER, PARAMETER :: ERR_STATS_INPUT      =  75, &  ! time avg statistics
                        ERR_STATS_RESTART    =  76, &
                        ERR_STATS_INDEXING   =  77, &
                        ERR_STATS_TECPLOT    =  78     ! std tecplot stats
#endif

  INTEGER, PARAMETER :: ERR_OPTION_TYPE      =  80     ! multi-physics

  INTEGER, PARAMETER :: ERR_MERGE_SORTED     =  90, &  ! sorting and searching
                        ERR_BINARY_SEARCH    =  91

  INTEGER, PARAMETER :: ERR_GRAD_MISMATCH    = 100     ! solution algorithm

  INTEGER, PARAMETER :: ERR_NEGATIVE_POSDEF  = 200, &
                        ERR_INVALID_VALUE    = 201

  INTEGER, PARAMETER :: ERR_SEC_READ_TWICE   = 240     ! reading sections

  INTEGER, PARAMETER :: ERR_MP_ALLORNONE     = 250     ! MP module use

  INTEGER, PARAMETER :: ERR_MISSING_VALUE    = 260     ! IO value missing 

  INTEGER, PARAMETER :: ERR_UNKNOWN_VISCMODEL = 300, & ! viscosity model
                        ERR_RK_SCHEME_INVALID = 301

  INTEGER, PARAMETER :: ERR_RAND_SEED_TYPE_INVALID = 400 ! random number generator

! ==============================================================================
! ROCFLO specific errors
! ==============================================================================

  INTEGER, PARAMETER :: ERR_REGION_RANGE     = 500, &  ! regions & patches
                        ERR_PATCH_2ALIGN     = 501, &
                        ERR_PATCH_NOALIGN    = 502, &
                        ERR_PATCH_NOSOURCE   = 503, &
                        ERR_PATCH_DIMENS     = 504, &
                        ERR_PATCH_NOTCOVERED = 505, &
                        ERR_WRONG_REGIONFACE = 506, &
                        ERR_PATCH_OVERLAP    = 507, &
                        ERR_PATCH_OVERSPEC   = 508, &
                        ERR_GRID_DIMENSIONS  = 509, &
                        ERR_GRID_DUMCELLS    = 510, &
                        ERR_SRCREGION_OFF    = 511, &
                        ERR_NUMBER_CELLS     = 512, &
                        ERR_REGION_NUMBER    = 513, &
                        ERR_PATCH_NUMBER     = 514

  INTEGER, PARAMETER :: ERR_VOLUME_EDGES     = 900, &  ! CV topology
                        ERR_VOLUME_CORNERS   = 901

! ==============================================================================
! ROCFLU specific errors
! ==============================================================================

  INTEGER, PARAMETER :: LIMIT_INFINITE_LOOP = 1E6

  INTEGER, PARAMETER :: ERR_VERTEX_NUMBER           = 1000, & ! data structure
                        ERR_NBFACES_WRONG           = 1001, &
                        ERR_PATCH_NUMBERING         = 1002, &
                        ERR_CELL_TYPE               = 1003, &
                        ERR_HASHTABLE               = 1004, &
                        ERR_NFACES_WRONG            = 1005, &
                        ERR_NEDGES_ESTIMATE         = 1006, &
                        ERR_NFACES_ESTIMATE         = 1007, &
                        ERR_NCELLS_WRONG            = 1008, &
                        ERR_VOLUME_DIFF             = 1009, &
                        ERR_VOLUME_NEGATIVE         = 1010, &
                        ERR_FACESUM                 = 1011, &
                        ERR_NBVERT_ESTIMATE         = 1012, &
                        ERR_NBVERT_EXTREMA          = 1013, &
                        ERR_BFACE_LIST_EXTREMA      = 1014, &
                        ERR_FACELIST_INVALID        = 1015, &
                        ERR_DIMENS_INVALID          = 1016, &
                        ERR_MOVEPATCH_BC_INVALID    = 1017, &
                        ERR_NEDGES_WRONG            = 1018, &
                        ERR_CELL_KIND_CHECK         = 1019, &
                        ERR_FACE_NVERT_INVALID      = 1020, &
                        ERR_PATCH_RENUMFLAG         = 1021, &
                        ERR_BVERT_LIST_INVALID      = 1022, &
                        ERR_EDGELIST_INVALID        = 1023, &
                        ERR_SVERT_LIST_INVALID      = 1024, &
                        ERR_MOVEPATCH_BC_NOTSET     = 1025, &
                        ERR_NCELLS_SPECIAL_MAX      = 1026, &
                        ERR_NFACES_SPECIAL_MAX      = 1027, & 
                        ERR_C2FLIST_INVALID         = 1028, &
                        ERR_FACE_KIND               = 1029, &
                        ERR_BFACEMEMBS_INVALID      = 1030, &
                        ERR_CELLGRAD_UNAVAILABLE    = 1031, &
                        ERR_BF2BG_INCONSISTENT      = 1032, &
                        ERR_C2CSKEY_INCONSISTENT    = 1033, &
                        ERR_DISTRIB_INVALID         = 1034, &
                        ERR_BCDATA_VALUE_INVALID    = 1035, &
                        ERR_INDMFMIXT_INVALID       = 1036, & 
                        ERR_DENUMBER_LIST           = 1037, & 
                        ERR_PROC2REG_MAPPING        = 1038, & 
                        ERR_OPTIONAL_MISSING        = 1039, & 
                        ERR_NVERT_ESTIMATE          = 1040, & 
                        ERR_NFACESCUT_INVALID       = 1041, & 
                        ERR_CELL_NOT_FOUND          = 1042, & 
                        ERR_VERTEX_NOT_FOUND        = 1043, & 
                        ERR_NBORDERS_INVALID        = 1044, & 
                        ERR_REGION_IDS_INVALID      = 1045, & 
                        ERR_CELL_KIND_INVALID       = 1046, & 
                        ERR_CELL_TYPE_INVALID       = 1047, & 
                        ERR_VERTEX_KIND_INVALID     = 1048, &
                        ERR_PARTITION_INVALID       = 1049, & 
                        ERR_DATADIM_MISMATCH        = 1050, & 
                        ERR_BUFFERDIM_MISMATCH      = 1051, & 
                        ERR_NREGIONS_MISMATCH       = 1052, &                 
                        ERR_NPROCS_MISMATCH         = 1053, & 
                        ERR_NUM_BC_VIRTUAL          = 1054, & 
                        ERR_NVERTSHARED_MISMATCH    = 1055, & 
                        ERR_LUBOUND_MISMATCH        = 1056, & 
                        ERR_ALLOCATE_ADAPTIVE       = 1057, & 
                        ERR_RECONST_INVALID         = 1058, & 
                        ERR_DISCR_INVALID           = 1059, & 
                        ERR_ORDER_INVALID           = 1060, & 
                        ERR_SOLVER_TYPE_INVALID     = 1061, & 
                        ERR_BF2CSORTED_INVALID      = 1062, & 
                        ERR_CONSTR_INVALID          = 1063, & 
                        ERR_GASMODEL_INVALID        = 1064, & 
                        ERR_GASMODEL_DISCR_MISMATCH = 1065, &
                        ERR_FACE_NORMAL_INVALID     = 1066, & 
                        ERR_TOLERICT_INVALID        = 1067, & 
                        ERR_STENCILDIMENS_INVALID   = 1068, &
                        ERR_STENCILMEMBER_INVALID   = 1069, & 
                        ERR_BC_INVALID              = 1070, & 
                        ERR_PATCH_BC_INCONSISTENT   = 1071, & 
                        ERR_REGION_ID_INVALID       = 1072, & 
                        ERR_FLATFLAG_INCONSISTENT   = 1073, & 
                        ERR_VERTEX_MATCH_FAILED     = 1074, & 
                        ERR_BORDER_INDEX_INVALID    = 1075, & 
                        ERR_REGION_ID_NOT_FOUND     = 1076, & 
                        ERR_PATCH_NOT_FLAT          = 1077, & 
                        ERR_PATCH_NOT_ALIGNED       = 1078, &
                        ERR_VIRTUALCELLS_NOTDB2     = 1079, &
                        ERR_BCVAR_VALUE_INVALID     = 1080
  
  INTEGER, PARAMETER :: ERR_INVALID_MARKER   = 1100, & ! I/O
                        ERR_INVALID_NCELLS   = 1101, &
                        ERR_INVALID_NVARS    = 1102

  INTEGER, PARAMETER :: ERR_OLES_STENCIL     = 1200, & ! Optimal LES
                        ERR_OLES_FLOWMODEL   = 1201

#ifdef GENX
  INTEGER, PARAMETER :: ERR_BCCOUPLED_NONE   = 1300, & ! GENx
                        ERR_MDOT_NEGATIVE    = 1301, &
                        ERR_TFLM_NEGATIVE    = 1302, & 
                        ERR_TB_NEGATIVE      = 1303
#endif

  INTEGER, PARAMETER :: ERR_NDIMENS_INVALID  = 1400, & ! Grid conversion
                        ERR_NZONES_INVALID   = 1401, &
                        ERR_FACETYPE_INVALID = 1402, &
                        ERR_C2VLIST_INVALID  = 1403, &
                        ERR_FACE_ORIENT      = 1404, & 
                        ERR_STRING_INVALID   = 1405, & 
                        ERR_NTYPE_INVALID    = 1406, & 
                        ERR_NDP_INVALID      = 1407

  INTEGER, PARAMETER :: ERR_LAPACK_OUTPUT    = 1500, & ! LAPACK
                        ERR_DCUHRE_OUTPUT    = 1510, & ! DCUHRE
                        ERR_TECPLOT_OUTPUT   = 1520, & ! TECPLOT
                        ERR_TECPLOT_FILECNTR = 1521, & 
                        ERR_PETSC_OUTPUT     = 1530, & ! PETSC
                        ERR_MPI_OUTPUT       = 1540, & 
                        ERR_MPI_TAGMAX       = 1541

  INTEGER, PARAMETER :: ERR_INFINITE_LOOP              = 1900, & ! Miscellaneous
                        ERR_PREC_RANGE                 = 1901, &
                        ERR_CV_STATE_INVALID           = 1902, &
                        ERR_FILEDEST_INVALID           = 1903, &
                        ERR_POST_OUTPUT_FORMAT_INVALID = 1904, &
                        ERR_POST_NSERVERS_INVALID      = 1905, &  
                        ERR_EXCEED_DIMENS              = 1906, & 
                        ERR_STRING_READ                = 1907

! ==============================================================================
! ROCTURB specific errors
! ==============================================================================

  INTEGER, PARAMETER :: ERR_TURB_MODULE      = 2001    ! TURB-module activation

#ifdef TURB
  INTEGER, PARAMETER :: ERR_TURB_MODEL       = 2005, & ! turbulence model
                        ERR_TURB_FIXPARAM    = 2006, & ! fixed parameters
                        ERR_TURB_INPUT       = 2007, & ! general turb.input
                        ERR_TURB_STATSINPUT  = 2008, & ! turb statistics input
                        ERR_TURB_LESINPUT    = 2009, & ! LES  input-options
                        ERR_TURB_RANSINPUT   = 2010, & ! RANS input-options
                        ERR_TURB_DESINPUT    = 2011, & ! DES  input-options
                        ERR_TURB_WLMINPUT    = 2012, & ! WLM  input-options
                        ERR_TURB_WLMMETRIC   = 2013    ! WLM metric

  INTEGER, PARAMETER :: ERR_TURB_REGION      = 2020    ! TURB-active regions

  INTEGER, PARAMETER :: ERR_TURB_WORKSPACE   = 2100    ! rocturb workspace
#endif

! ==============================================================================
! ROCPERI specific errors
! ==============================================================================

#ifdef PERI
  INTEGER, PARAMETER :: ERR_PERI_INPUT       = 2501, & ! periodic flow input
                        ERR_PERI_CPRBC       = 2502, & ! cpr boundary condition
                        ERR_PERI_GEO         = 2503, & ! geometry
                        ERR_PERI_FIXPARAM    = 2504, & ! fixed parameters
                        ERR_PERI_PHYSPARAM   = 2505, & ! physical parameters
                        ERR_PERI_MPI         = 2506    ! MPI related errors

  INTEGER, PARAMETER :: ERR_PERI_LAST        = 2600    ! upperbound error id
#endif

! ==============================================================================
! ROCRAD specific errors
! ==============================================================================

#ifdef RADI
  INTEGER, PARAMETER :: ERR_RADI_INPUT       = 2601, & ! radiation user input
                        ERR_RADI_BC          = 2602, & ! boundary conditions
                        ERR_RADI_FIXPARAM    = 2603, & ! fixed parameters
                        ERR_RADI_METRIC      = 2604, & ! radiation metric
                        ERR_RADI_MULPHASE    = 2605, & ! multiphase
                        ERR_RADI_FLDINPUT    = 2606, & ! FLD input option
                        ERR_RADI_REGION      = 2610    ! RADI-active regions

  INTEGER, PARAMETER :: ERR_RADI_LAST        = 2700    ! upperbound error id
#endif

! ==============================================================================
! ROCPART specific errors
! ==============================================================================

  INTEGER, PARAMETER :: ERR_PLAG_MODULE        = 3001    ! PLAG-module activation

#ifdef PLAG
  INTEGER, PARAMETER :: ERR_PLAG_INTRPLMODEL       = 3007, & ! interpolation model
                        ERR_PLAG_CELLINDEX         = 3008, & ! cell indices
                        ERR_PLAG_TILESIZE          = 3009, & ! tile size
                        ERR_PLAG_BREAKUPMODEL      = 3010, & ! breakup model
                        ERR_PLAG_BREAKUPFAC        = 3011, & ! breakup factor
                        ERR_PLAG_BREAKUPWEBSWI     = 3012, & ! breakup Weber switch
                        ERR_PLAG_INJCDIAMDIST      = 3013, & ! injection diameter dist model
                        ERR_PLAG_INJCDIAM          = 3014, & ! injection diameter
                        ERR_PLAG_INVALID_NPCLS     = 3015, &
                        ERR_PLAG_PCL_NOT_FOUND     = 3016, & 
                        ERR_PLAG_NCONT_INVALID     = 3017, &
                        ERR_PLAG_EJECMODEL         = 3018, & ! ejection model
                        ERR_PLAG_MEMOVERFLOW       = 3019, & ! memory overflow
                        ERR_PLAG_FINDPCL           = 3020, & ! find cell method 
                        ERR_PLAG_IPCLSEND_OVERFLOW = 3021, & ! iPclSend memory overflow
                        ERR_PLAG_DSTR_INVALID      = 3022
#endif

! ==============================================================================
! ROCSMOKE specific errors
! ==============================================================================

  INTEGER, PARAMETER :: ERR_PEUL_MODULE      = 4001    ! PEUL-module activation

#ifdef PEUL
  INTEGER, PARAMETER :: ERR_PEUL_NPTYPES     = 4010, & ! input Errors
                        ERR_PEUL_PTYPE       = 4011, &
                        ERR_PEUL_BADVAL      = 4013, &
                        ERR_PEUL_NPMISMATCH  = 4020, &
                        ERR_PEUL_BCVAL_EXTRA = 4030, &
                        ERR_PEUL_EXTERNAL    = 4900, & ! non-implementations
                        ERR_PEUL_DISTRIB     = 4901, &
                        ERR_PEUL_MOVEGRID    = 4902
#endif

! ==============================================================================
! ROCINTERACT specific errors
! ==============================================================================

  INTEGER, PARAMETER :: ERR_INRT_MODULE      = 5001    ! INRT-module activation

#ifdef INRT
  INTEGER, PARAMETER :: ERR_INRT_DEFREAD     = 5010, &
                        ERR_INRT_DEFUNREAD   = 5011, &
                        ERR_INRT_READ        = 5015, &
                        ERR_INRT_MULTPLAGMAT = 5020, &
                        ERR_INRT_MISSPLAGMAT = 5021, &
                        ERR_INRT_MISSINGMAT  = 5025, &
                        ERR_INRT_ALLOCRANGE  = 5030, &
                        ERR_INRT_INDEXRANGE  = 5031, &
                        ERR_INRT_BADSWITCH   = 5040, &
                        ERR_INRT_BADACTV     = 5041, &
                        ERR_INRT_BADPERM     = 5042, &
                        ERR_INRT_BADVAL      = 5043, &
                        ERR_INRT_BADMAT      = 5044, &
                        ERR_INRT_MISSINGVAL  = 5045, &
                        ERR_INRT_ACTVPLAG    = 5050, &
                        ERR_INRT_ACTVSMOKE   = 5055, &
                        ERR_INRT_ACTVDSMOKE  = 5056, &
                        ERR_INRT_PARAMETER   = 5060, &
                        ERR_INRT_NINTL       = 5070, &
                        ERR_INRT_NINPUTEDGES = 5071, &
                        ERR_INRT_CONNECTINTL = 5072, &
                        ERR_INRT_PERMINTL    = 5073, &
                        ERR_INRT_PERMLEVINTL = 5074, &
                        ERR_INRT_BURNING1    = 5080, &
                        ERR_INRT_ONLY1       = 5081, &
                        ERR_INRT_OX_ACTV     = 5082, &
                        ERR_INRT_BOIL_ACTV   = 5083, &
                        ERR_INRT_BOIL_SAME   = 5084, &
                        ERR_INRT_NPCLS       = 5090, &
                        ERR_INRT_ENERVAPOR   = 5100, &
                        ERR_INRT_NOINRT      = 5900    ! non-implementations
#endif


! ==============================================================================
! ROCSPECIES specific errors
! ==============================================================================

  INTEGER, PARAMETER :: ERR_SPEC_MODULE              = 6001, &
                        ERR_SPEC_NTYPES              = 6002, & 
                        ERR_SPEC_MAXEQN              = 6003, & 
                        ERR_SPEC_PROPS_INVALID       = 6004, & 
                        ERR_SPEC_NSPEC_INVALID       = 6005, & 
                        ERR_SPEC_SOURCE_TYPE_INVALID = 6006

! ******************************************************************************
! Error & debug functions
! ******************************************************************************

  CONTAINS

! ==============================================================================
!   Register new function call
! ==============================================================================

    SUBROUTINE RegisterFunction( global,funName,fileName )
      USE ModGlobal, ONLY : t_global
      TYPE(t_global), POINTER :: global
      CHARACTER(*) :: funName, fileName

      IF (global%nFunTree<0 .OR. &
          global%nFunTree>=UBOUND(global%functionTree,2)) THEN ! wrong dimension
        CALL ErrorStop( global,ERR_REGISTER_FUN,__LINE__ )
      ENDIF

      global%nFunTree = global%nFunTree + 1

      global%functionTree(1,global%nFunTree) = funName
      global%functionTree(2,global%nFunTree) = fileName
    END SUBROUTINE RegisterFunction

! ==============================================================================
!   Deregister function call
! ==============================================================================

    SUBROUTINE DeregisterFunction( global )
      USE ModGlobal, ONLY : t_global
      TYPE(t_global), POINTER :: global

      global%functionTree(1,global%nFunTree) = ''
      global%functionTree(2,global%nFunTree) = ''

      global%nFunTree = global%nFunTree - 1
    END SUBROUTINE DeregisterFunction

! ==============================================================================
!   Error function
! ==============================================================================

    SUBROUTINE ErrorStop( global,errorCode,errorLine,addMessage )
      USE ModDataTypes
      USE ModGlobal, ONLY : t_global
      USE ModMPI
      USE ModParameters

      LOGICAL :: flag
      INTEGER :: errorCode,errorCode2,errorFlag,errorLine
      INTEGER :: i, error
      CHARACTER(*), OPTIONAL  :: addMessage
      CHARACTER(2*CHRLEN)     :: message
      TYPE(t_global), POINTER :: global

      errorCode2 = 1

      SELECT CASE (errorCode)

! ------------------------------------------------------------------------------
!       Basics
! ------------------------------------------------------------------------------

        CASE (ERR_REGISTER_FUN)
          message = 'dimension of <global%functionTree> out of bounds.'
        CASE (ERR_REACHED_DEFAULT)
          message = 'reached default statement in select construct.'
        CASE (ERR_PREVIOUS_ERRORS)
          message = 'aborting due to previous errors.'
        CASE (ERR_EXTERNAL_FUNCT)
          message = 'you have to link to an external function.'
        CASE (ERR_EXCEEDS_DECL_MEM)
          message = 'array index or size exceeds declared memory.'
        CASE (ERR_UNKNOWN_OPTION)
          message = 'option unknown or not yet implemented.'
        CASE (ERR_ILLEGAL_VALUE)
          message = 'variable has illegal value.'
        CASE (ERR_SYSTEM_COMMAND)
          message = 'fail to execute system command:'
        CASE (ERR_COMPILE_OPTION)
          message = 'incorrect compilation option.'

! ------------------------------------------------------------------------------
!       Memory allocation
! ------------------------------------------------------------------------------

        CASE (ERR_ALLOCATE)
          message = 'cannot allocate memory.'
        CASE (ERR_DEALLOCATE)
          message = 'cannot deallocate memory.'
        CASE (ERR_ASSOCIATED)
          message = 'pointer not associated.'
        CASE (ERR_ALLOCATED)
          message = 'variable already allocated.'

! ------------------------------------------------------------------------------
!       Input parameters
! ------------------------------------------------------------------------------

        CASE (ERR_DEPENDENT_INPUT)
          message = 'inconsistent dependent input.'

! ------------------------------------------------------------------------------
!       I/O
! ------------------------------------------------------------------------------

        CASE (ERR_FILE_OPEN)
          message = 'cannot open file:'
        CASE (ERR_FILE_CLOSE)
          message = 'cannot close file:'
        CASE (ERR_FILE_READ)
          message = 'cannot read from file:'
        CASE (ERR_FILE_WRITE)
          message = 'cannot write to file:'
        CASE (ERR_FILE_EXIST)
          message = 'file does not exist:'
        CASE (ERR_UNKNOWN_FORMAT)
          message = 'unknown file format.'
        CASE (ERR_PROBE_LOCATION)
          message = 'probe located outside the flow domain.'
        CASE (ERR_PROBE_SPECIFIED)
          message = 'probe(s) already specified.'
        CASE (ERR_VAL_UNDEFINED)
          message = 'value undefined:'
        CASE (ERR_STOPFILE_FOUND)
          message = 'File STOP found. Delete and restart run.'

! ------------------------------------------------------------------------------
!       Boundary conditions
! ------------------------------------------------------------------------------

        CASE (ERR_UNKNOWN_BC)
          message = 'unknown type of boundary condition.'
        CASE (ERR_NO_BCSPECIFIED)
          message = 'boundary condition not specified.'
        CASE (ERR_BCVAL_MISSING)
          message = 'boundary value(s) missing.'
        CASE (ERR_NO_BCSWITCH)
          message = 'BC switch not specified.'
        CASE (ERR_VAL_BCSWITCH)
          message = 'BC switch value is not valid.'
        CASE (ERR_PATCH_RANGE)
          message = 'patch indices outside of region`s dimensions.'
        CASE (ERR_BC_VARNAME)
          message = 'invalid variable name for BC type:'
        CASE (ERR_VAL_BCVAL)
          message = 'BC value is not valid.'
        CASE (ERR_TBC_NONINIT)
          message = 'TBC not initialized.'

! ------------------------------------------------------------------------------
!       Time stepping
! ------------------------------------------------------------------------------

        CASE (ERR_TIME_GRID)
          message = 'wrong physical time in grid file.'
        CASE (ERR_TIME_SOLUTION)
          message = 'wrong physical time in solution file.'
        CASE (ERR_TIME_GLOBAL)
          message = 'physical time differs between processors.'
        CASE (ERR_DTIME_NEGATIVE)
          message = 'current time is later than the max. simulation time.'
        CASE (ERR_DITER_NEGATIVE)
          message = 'iteration number is higher than the max. allowed one.'
        CASE (ERR_STEADY)
          message = 'steady flow not allowed.'

! ------------------------------------------------------------------------------
!       MPI
! ------------------------------------------------------------------------------

        CASE (ERR_MPI_TROUBLE)
          message = 'MPI does not work.'
        CASE (ERR_NO_PROCMAP)
          message = 'no mapping to processors.'
        CASE (ERR_NO_PROCMATCH)
          message = 'no. of regions does not match no. of processors.'
        CASE (ERR_PATCH_OFFSET)
          message = 'no. of regions on processor > MPI_PATCHOFF.'
        CASE (ERR_PROC_MISMATCH)
          message = 'no. of procs does not match specified no. of procs.'

! ------------------------------------------------------------------------------
!       Multigrid
! ------------------------------------------------------------------------------

        CASE (ERR_GRID_LEVEL)
          message = 'no such grid level possible.'

! ------------------------------------------------------------------------------
!       Velocity/temperature gradients
! ------------------------------------------------------------------------------

        CASE (ERR_GRAD_INDEX)
          message = 'inconsistent velocity/temperature gradient indexing.'

! ------------------------------------------------------------------------------
!       Geometry of control volume
! ------------------------------------------------------------------------------

        CASE (ERR_FACEVEC_SUM)
          message = 'sum of face vectors > 0.'
        CASE (ERR_VOLUME_SIZE)
          message = 'volume size below threshold or negative.'
        CASE (ERR_FACE_SPLIT)
          message = 'face splitting detected negative dotproduct of split face vectors.'
        CASE (ERR_FACE_INVERTED)
          message = 'inverted cell-face detected.'

#ifdef STATS
! ------------------------------------------------------------------------------
!       Time averaged statistics
! ------------------------------------------------------------------------------

        CASE (ERR_STATS_INPUT)
          message = 'required statistics input parameter not found:'
        CASE (ERR_STATS_RESTART)
          message = 'inconsistent statistics restart parameters.'
        CASE (ERR_STATS_INDEXING)
          message = 'inconsistent statistics indexing.'
        CASE (ERR_STATS_TECPLOT)
          message = 'violation of fixed nStat and/or statId for tecplot.'
#endif

! ------------------------------------------------------------------------------
!       Multi-physics
! ------------------------------------------------------------------------------
        CASE (ERR_OPTION_TYPE)
          message = 'unknown or unimplemented option.'

! ------------------------------------------------------------------------------
!       Sorting and searching
! ------------------------------------------------------------------------------

        CASE (ERR_MERGE_SORTED)
          message = 'Error in merging sorted lists.'
        CASE (ERR_BINARY_SEARCH)
          message = 'Error in binary search.'

! ------------------------------------------------------------------------------
!       Solution algorithm
! ------------------------------------------------------------------------------

        CASE (ERR_GRAD_MISMATCH)
          message = 'mismatch of no of variables and gradients'

! ------------------------------------------------------------------------------
!       Posivity/validity checking
! ------------------------------------------------------------------------------

        CASE (ERR_NEGATIVE_POSDEF)
          message = 'Negative positive-definite quantity detected.'
        CASE (ERR_INVALID_VALUE)
          message = 'Invalid quantity detected.'

! ------------------------------------------------------------------------------
!       Reading sections
! ------------------------------------------------------------------------------

        CASE (ERR_SEC_READ_TWICE)
          message = 'Attempted two different input sections for same region.'
        CASE (ERR_MISSING_VALUE)
          message = 'Value expected in input but found missing.'

! ------------------------------------------------------------------------------
!       MP module usage
! ------------------------------------------------------------------------------

        CASE (ERR_MP_ALLORNONE)
          message = 'MP modules must be used in all regions or in none.'

! ------------------------------------------------------------------------------
!       Miscellaneous
! ------------------------------------------------------------------------------

        CASE (ERR_UNKNOWN_VISCMODEL)
          message = 'Unknown or unimplemented viscosity model.'
        CASE (ERR_RK_SCHEME_INVALID)
          message = 'Chosen RK scheme invalid.'   
        CASE (ERR_RAND_SEED_TYPE_INVALID)
          message = 'Chosen random seed type invalid.'   

! ------------------------------------------------------------------------------
!       ROCFLO
! ------------------------------------------------------------------------------

! ----- Regions & patches ------------------------------------------------------

        CASE (ERR_REGION_RANGE)
          message = 'number of source region out of range.'
        CASE (ERR_PATCH_2ALIGN)
          message = 'both patch coordinates aligned.'
        CASE (ERR_PATCH_NOALIGN)
          message = 'none of patch coordinates aligned.'
        CASE (ERR_PATCH_NOSOURCE)
          message = 'cannot find matching source patch.'
        CASE (ERR_PATCH_DIMENS)
          message = 'patch dimension does not match distribution.'
        CASE (ERR_PATCH_NOTCOVERED)
          message = 'no boundary condition for some faces of a patch.'
        CASE (ERR_WRONG_REGIONFACE)
          message = 'face number for a region outside 1-6'
        CASE (ERR_PATCH_OVERLAP)
          message = 'boundary patches do overlap.'
        CASE (ERR_PATCH_OVERSPEC)
          message = 'boundary conditions already defined for this patch.'
        CASE (ERR_GRID_DIMENSIONS)
          message = 'grid dimensions not the same as in topology file.'
        CASE (ERR_GRID_DUMCELLS)
          message = 'solution file contains different number of dummy cells.'
        CASE (ERR_SRCREGION_OFF)
          message = 'source region is inactive.'
        CASE (ERR_NUMBER_CELLS)
          message = 'not enough cells to contain dummy cells of adjacent region.'
        CASE (ERR_REGION_NUMBER)
          message = 'got different region number from file.'
        CASE (ERR_PATCH_NUMBER)
          message = 'got different patch number from file.'

! ----- Topology of control volume ----------------------------------------------

        CASE (ERR_VOLUME_EDGES)
          message = 'number of edges out of range (1-12).'
        CASE (ERR_VOLUME_CORNERS)
          message = 'number of corners out of range (1-8).'

! ------------------------------------------------------------------------------
!       ROCFLU
! ------------------------------------------------------------------------------

        CASE (ERR_VERTEX_NUMBER)
          message = 'invalid vertex number in array.'
        CASE (ERR_NBFACES_WRONG)
          message = 'computed number of boundary faces inconsistent.'
        CASE (ERR_PATCH_NUMBERING)
          message = 'patch numbering inconsistent.'
        CASE (ERR_CELL_TYPE)
          message = 'invalid cell type:'
        CASE (ERR_HASHTABLE)
          message = 'internal inconsistency in hash table.'
        CASE (ERR_NFACES_WRONG)
          message = 'computed number of internal faces inconsistent.'
        CASE (ERR_NEDGES_ESTIMATE)
          message = 'Estimate of number of edges too low.'
        CASE (ERR_NFACES_ESTIMATE)
          message = 'Estimate of number of faces too low.'
        CASE (ERR_NCELLS_WRONG)
          message = 'computed number of cells inconsistent.'
        CASE (ERR_VOLUME_DIFF)
          message = 'absolute difference in volumes larger than specified '// &
                    'limit.'
        CASE (ERR_VOLUME_NEGATIVE)
          message = 'Negative volume(s) detected.'
        CASE (ERR_FACESUM)
          message = 'face sum greater than minimum face area.'
        CASE (ERR_NBVERT_ESTIMATE)
          message = 'Estimate of number of boundary vertices too low.'
        CASE (ERR_NBVERT_EXTREMA)
          message = 'Boundary vertex list has invalid extrema.'
        CASE (ERR_BFACE_LIST_EXTREMA)
          message = 'Locally-numbered boundary-face list has invalid extrema.'
        CASE (ERR_FACELIST_INVALID)
          message = 'Face list is invalid.'
        CASE (ERR_DIMENS_INVALID)
          message = 'Dimension invalid'
        CASE (ERR_MOVEPATCH_BC_INVALID)
          message = 'Invalid patch-motion boundary condition.'
        CASE (ERR_NEDGES_WRONG)
          message = 'Computed number of edges inconsistent.'
        CASE (ERR_CELL_KIND_CHECK)
          message = 'Actual cell touches at least one virtual vertex.'
        CASE (ERR_FACE_NVERT_INVALID)
          message = 'Number of vertices in face inconsistent.'
        CASE (ERR_PATCH_RENUMFLAG)
          message = 'Patch vertex renumbering flag is invalid.'
        CASE (ERR_BVERT_LIST_INVALID)
          message = 'Patch vertex list is invalid.'
        CASE (ERR_EDGELIST_INVALID)
          message = 'Edge list is invalid.'
        CASE (ERR_SVERT_LIST_INVALID)
          message = 'Special vertex list is invalid.'
        CASE (ERR_MOVEPATCH_BC_NOTSET)
          message = 'One or more patches without boundary grid-motion bc.'
        CASE (ERR_NCELLS_SPECIAL_MAX)
          message = 'Exceeded maximum allowed number of special cells.'
        CASE (ERR_NFACES_SPECIAL_MAX)
          message = 'Exceeded maximum allowed number of special faces.'          
        CASE (ERR_C2FLIST_INVALID)
          message = 'Cell-to-face list is invalid.'
        CASE (ERR_FACE_KIND)
          message = 'Invalid face kind:'
        CASE (ERR_BFACEMEMBS_INVALID)
          message = 'List of boundary faces for boundary gradients invalid.'
        CASE (ERR_CELLGRAD_UNAVAILABLE)
          message = 'Cell-gradients not available.'
        CASE (ERR_BF2BG_INCONSISTENT)
          message = 'Boundary-face gradient-access list inconsistent.'
        CASE (ERR_C2CSKEY_INCONSISTENT)
          message = 'Cell-to-cell stencil access list inconsistent.'
        CASE (ERR_DISTRIB_INVALID)
          message = 'Distribution parameter invalid:'
        CASE (ERR_BCDATA_VALUE_INVALID)
          message = 'Invalid value read:'
        CASE (ERR_BCVAR_VALUE_INVALID)
          message = 'Invalid value read for boundary variables:'
        CASE (ERR_INDMFMIXT_INVALID)
          message = 'Invalid value for indMfMixt.'
        CASE (ERR_DENUMBER_LIST) 
          message = 'Denumbering index invalid.'
        CASE (ERR_PROC2REG_MAPPING)
          message = 'Process-to-region mapping invalid.'
        CASE (ERR_OPTIONAL_MISSING) 
          message = 'Optional argument missing.'
        CASE (ERR_NVERT_ESTIMATE)
          message = 'Estimate of number of vertices too low.'
        CASE (ERR_NFACESCUT_INVALID)
          message = 'Computed number of cut faces invalid:'
        CASE (ERR_CELL_NOT_FOUND) 
          message = 'Cell not found during search.'
        CASE (ERR_VERTEX_NOT_FOUND) 
          message = 'Vertex not found during search.'
        CASE (ERR_NBORDERS_INVALID)
          message = 'Computed number of borders invalid:'
        CASE (ERR_REGION_IDS_INVALID)
          message = 'Region indices invalid:'
        CASE (ERR_CELL_KIND_INVALID)
          message = 'Cell kind invalid.'
        CASE (ERR_CELL_TYPE_INVALID)
          message = 'Cell type invalid:'
        CASE (ERR_VERTEX_KIND_INVALID)
          message = 'Vertex kind invalid.'
        CASE (ERR_PARTITION_INVALID)
          message = 'Partitioning invalid.'
        CASE (ERR_DATADIM_MISMATCH)
          message = 'Data dimensions mismatch.'
        CASE (ERR_BUFFERDIM_MISMATCH)
          message = 'Data dimensions mismatch.'
        CASE (ERR_NREGIONS_MISMATCH)
          message = 'Numbers of regions do not match.'
        CASE (ERR_NPROCS_MISMATCH)
          message = 'Numbers of processors do not match.' 
        CASE (ERR_NUM_BC_VIRTUAL) 
          message = 'Number of virtual boundaries incorrect.'
        CASE (ERR_NVERTSHARED_MISMATCH)
          message = 'Number of shared vertices does not match.'
        CASE (ERR_LUBOUND_MISMATCH)
          message = 'Lower and/or upper bounds do not match.'
        CASE (ERR_ALLOCATE_ADAPTIVE)
          message = 'Adaptive memory allocation failed.'
        CASE (ERR_RECONST_INVALID)
          message = 'Reconstruction method invalid.'
        CASE (ERR_DISCR_INVALID)
          message = 'Discretization method invalid.'
        CASE (ERR_ORDER_INVALID)
          message = 'Order invalid.'
        CASE (ERR_SOLVER_TYPE_INVALID)
          message = 'Invalid solver type.'      
        CASE (ERR_BF2CSORTED_INVALID)
          message = 'Boundary-face-to-cell list invalid.'
        CASE (ERR_CONSTR_INVALID)
          message = 'Constraint method invalid.'
        CASE (ERR_GASMODEL_INVALID)
          message = 'Gas model invalid:'          
        CASE (ERR_GASMODEL_DISCR_MISMATCH)
          message = 'Gas model and discretization method do not match.' 
        CASE (ERR_FACE_NORMAL_INVALID)
          message = 'Face normal invalid.'
        CASE (ERR_TOLERICT_INVALID)
          message = 'In-cell test tolerance invalid.'
        CASE (ERR_STENCILDIMENS_INVALID)
          message = 'Stencil dimensionality invalid.'
        CASE (ERR_STENCILMEMBER_INVALID)
          message = 'Stencil member invalid.'             
        CASE (ERR_BC_INVALID)
          message = 'Invalid boundary condition.'
        CASE (ERR_PATCH_BC_INCONSISTENT)
          message = 'Boundary condition inconsistent with patch geometry.'
        CASE (ERR_REGION_ID_INVALID)
          message = 'Region index invalid.'
        CASE (ERR_FLATFLAG_INCONSISTENT)
          message = 'Patch flatness flags inconsistent.'
        CASE (ERR_VERTEX_MATCH_FAILED)
          message = 'Vertex matching failed.'
        CASE (ERR_BORDER_INDEX_INVALID)
          message = 'Border index invalid.'
        CASE (ERR_REGION_ID_NOT_FOUND)
          message = 'Region index not found.'   
        CASE (ERR_PATCH_NOT_FLAT)
          message = 'Patch not flat.'
        CASE (ERR_PATCH_NOT_ALIGNED)
          message = 'Patch not aligned with coordinate axes.'       
        CASE (ERR_VIRTUALCELLS_NOTDB2)
          message = 'Virtual cells not divisible by 2.'
          
        CASE (ERR_INVALID_MARKER)
          message = 'invalid section marker:'
        CASE (ERR_INVALID_NCELLS)
          message = 'number of cells invalid.'
        CASE (ERR_INVALID_NVARS)
          message = 'number of variables invalid.'

        CASE (ERR_OLES_STENCIL)
          message = 'inconsistency in stencil construction.'
        CASE (ERR_OLES_FLOWMODEL)
          message = 'Optimal LES approach needs Navier-Stokes flow model.'

#ifdef GENX
        CASE (ERR_BCCOUPLED_NONE)
          message = 'No coupled boundaries defined.'
        CASE (ERR_MDOT_NEGATIVE)
          message = 'Received negative mass flux on burning patch.'
        CASE (ERR_TFLM_NEGATIVE)
          message = 'Received negative flame temperature on burning patch.'
        CASE (ERR_TB_NEGATIVE)
          message = 'Received negative temperature on non-burning patch.'         
#endif

        CASE (ERR_NDIMENS_INVALID)
          message = 'Number of dimensions invalid.'
        CASE (ERR_NZONES_INVALID)
          message = 'Number of zones invalid.'
        CASE (ERR_FACETYPE_INVALID)
          message = 'Face type invalid.'
        CASE (ERR_C2VLIST_INVALID)
          message = 'Cell-to-vertex list is invalid.'
        CASE (ERR_FACE_ORIENT)
          message = 'Face-orientation check failed.'
        CASE (ERR_STRING_INVALID)
          message = 'Section string invalid.'
        CASE (ERR_NTYPE_INVALID)
          message = 'Element type invalid.'
        CASE (ERR_NDP_INVALID)
          message = 'Number of nodes invalid for given element type.'          

        CASE (ERR_LAPACK_OUTPUT)
          message = 'LAPACK routine returned non-zero info variable.'
        CASE (ERR_DCUHRE_OUTPUT)
          message = 'DCUHRE routine returned non-zero error variable.'
        CASE (ERR_TECPLOT_OUTPUT)
          message = 'TECPLOT routine returned non-zero error variable.'
        CASE (ERR_TECPLOT_FILECNTR)
          message = 'TECPLOT file counter exceeds maximum.'
        CASE (ERR_PETSC_OUTPUT) 
          message = 'PETSc routine returned non-zero error variable.'
        CASE (ERR_MPI_OUTPUT)
          message = 'MPI routine returned non-zero error variable.'
        CASE (ERR_MPI_TAGMAX)
          message = 'Exceeded maximum tag value allowed by MPI.'

        CASE (ERR_INFINITE_LOOP)
          message = 'detected what appears to be an infinite loop.'
        CASE (ERR_PREC_RANGE)
          message = 'incompatible precision and range of file data.'
        CASE (ERR_CV_STATE_INVALID)
          message = 'State of solution vector invalid.'
        CASE (ERR_FILEDEST_INVALID)
          message = 'File destination invalid.'
        CASE (ERR_POST_OUTPUT_FORMAT_INVALID)
          message = 'Postprocessing output format invalid.'
        CASE (ERR_POST_NSERVERS_INVALID)
          message = 'Number of servers invalid.'          
        CASE (ERR_EXCEED_DIMENS)
          message = 'Exceeding dimensions of array:'
        CASE (ERR_STRING_READ)
          message = 'Cannot read string.'

! ------------------------------------------------------------------------------
!       ROCTURB specific errors
! ------------------------------------------------------------------------------

        CASE (ERR_TURB_MODULE)   ! reserved for RFLO/U_CheckUserInput
          message = 'inconsistency of flow model and Rocturb module.'
#ifdef TURB
        CASE (ERR_TURB_MODEL)
          message = 'undefined turbulence model selected.'
        CASE (ERR_TURB_FIXPARAM)
          message = 'incorrect fixed parameters setting.'
        CASE (ERR_TURB_INPUT)
          message = 'violation in general turbulence input:'
        CASE (ERR_TURB_STATSINPUT)
          message = 'violation in turbulence statistics input:'
        CASE (ERR_TURB_LESINPUT)
          message = 'violation of LES input-options:'
        CASE (ERR_TURB_RANSINPUT)
          message = 'violation of RANS input-options:'
        CASE (ERR_TURB_DESINPUT)
          message = 'violation of DES input-options:'
        CASE (ERR_TURB_WLMINPUT)
          message = 'violation of WLM input-options:'
        CASE (ERR_TURB_WLMMETRIC)
          message = 'violation of WLM metric values:'
        CASE (ERR_TURB_REGION)
          message = 'turbulence model is not active in this region.'
        CASE (ERR_TURB_WORKSPACE)
          message = 'reserved workspace being used.'
#endif

! ------------------------------------------------------------------------------
!      ROCPERI specific errors
! ------------------------------------------------------------------------------

#ifdef PERI
        CASE (ERR_PERI_INPUT)
          message = 'inconsistency in PERI input parameters.'
        CASE (ERR_PERI_CPRBC)
          message = 'inconsistency in CPR boundary conditions.'
        CASE (ERR_PERI_GEO)
          message = 'inconsistency in PERI geometry.'
        CASE (ERR_PERI_FIXPARAM)
          message = 'inconsistency in PERI fixed parameters.'
        CASE (ERR_PERI_PHYSPARAM)
          message = 'inconsistency in PERI physical parameters.'
        CASE (ERR_PERI_MPI)
          message = 'MPI related error.'
#endif

! ------------------------------------------------------------------------------
!      ROCRAD specific errors
! ------------------------------------------------------------------------------

#ifdef RADI
        CASE (ERR_RADI_INPUT)
          message = 'violation of radiation input-options:'
        CASE (ERR_RADI_BC)
          message = 'inconsistency in radiation bc:'
        CASE (ERR_RADI_FIXPARAM)
          message = 'incorrect fixed parameters setting.'
        CASE (ERR_RADI_METRIC)
          message = 'inconsistency in radiation metric:'
        CASE (ERR_RADI_MULPHASE)
          message = 'multiphase inconsistency from other physics modules:'
        CASE (ERR_RADI_FLDINPUT)
          message = 'violation of FLD radiation input-options:'
        CASE (ERR_RADI_REGION)
          message = 'radiation model is not active in this region.'
#endif

! -------------------------------------------------------------------------------
!       ROCPART specific errors
! ------------------------------------------------------------------------------

        CASE (ERR_PLAG_MODULE)   ! used in PLAG_CheckUserInput
          message = 'inconsistency of flow model and RocPart module.'
#ifdef PLAG
        CASE (ERR_PLAG_INTRPLMODEL)
          message = 'undefined interpolation mixture model selected.'
        CASE (ERR_PLAG_CELLINDEX)
          message = 'unable to locate particle in cell.'
        CASE (ERR_PLAG_TILESIZE)
          message = 'got different tile size from file.'
        CASE (ERR_PLAG_BREAKUPMODEL)
          message = 'undefined breakup model selected.'
        CASE (ERR_PLAG_BREAKUPFAC)
          message = 'inconsistent breakup factor selected.'
        CASE (ERR_PLAG_BREAKUPWEBSWI)
          message = 'inconsistent breakup Weber switch selected.'
        CASE (ERR_PLAG_INJCDIAMDIST)
          message = 'undefined injection diameter distribution model selected.'
        CASE (ERR_PLAG_INJCDIAM)
          message = 'inconsistent injection diameters selected.'
        CASE (ERR_PLAG_INVALID_NPCLS)
          message = 'number of particles invalid:'
        CASE (ERR_PLAG_PCL_NOT_FOUND)
          message = 'Could not find particle:'
        CASE (ERR_PLAG_NCONT_INVALID)
          message = 'Number of constituents invalid.'
        CASE (ERR_PLAG_EJECMODEL)
          message = 'undefined ejection model selected.'
        CASE (ERR_PLAG_MEMOVERFLOW)
          message = 'Maximum dimension exceeded.'
        CASE (ERR_PLAG_FINDPCL)
          message = 'undefined find particle cell method selected.'
        CASE (ERR_PLAG_IPCLSEND_OVERFLOW)
          message = 'Maximum dimension of pBorder%iPclSend exceeded.'
        CASE (ERR_PLAG_DSTR_INVALID)
          message = 'Data structure invalid.'
#endif

! ------------------------------------------------------------------------------
!       ROCSMOKE specific errors
! ------------------------------------------------------------------------------

        CASE (ERR_PEUL_MODULE)
          message = 'inconsistency of flow model and Rocsmoke module.'
#ifdef PEUL
        CASE (ERR_PEUL_NPTYPES)
          message = 'number of Eulerian particle types not constant.'
        CASE (ERR_PEUL_PTYPE)
          message = 'a CONPART_PTYPE section occurs before any CONPART.'
        CASE (ERR_PEUL_BADVAL)
          message = 'invalid input value for Eulerian particles.'
        CASE (ERR_PEUL_NPMISMATCH)
          message = 'number of particle types not equal to nCv.'
        CASE (ERR_PEUL_BCVAL_EXTRA)
          message = 'boundary value specified for non-existent smoke type.'
        CASE (ERR_PEUL_EXTERNAL)
          message = 'external boundaries not yet implemented for PEUL.'
        CASE (ERR_PEUL_DISTRIB)
          message = 'distribution boundaries not yet implemented for PEUL.'
        CASE (ERR_PEUL_MOVEGRID)
          message = 'moving grids not yet implemented for PEUL.'
#endif

! ------------------------------------------------------------------------------
!       ROCINTERACT specific errors
! ------------------------------------------------------------------------------

        CASE (ERR_INRT_MODULE)
          message = 'inconsistency of flow model and Rocinteract module.'
#ifdef INRT
        CASE (ERR_INRT_DEFREAD)
          message = 'cannot read INRT_DEFAULT section twice for a region.'
        CASE (ERR_INRT_DEFUNREAD)
          message = 'cannot read an interaction without an INRT_DEFAULT.'
        CASE (ERR_INRT_READ)
          message = 'cannot read an interaction section twice for a region.'
        CASE (ERR_INRT_MULTPLAGMAT)
          message = 'two distinct PLAG constituents are the same material.'
        CASE (ERR_INRT_MISSPLAGMAT)
          message = 'material name read does not match any PLAG material.'
        CASE (ERR_INRT_MISSINGMAT)
          message = 'missing material name in input deck.'
        CASE (ERR_INRT_ALLOCRANGE)
          message = 'array allocated to the wrong size.'
        CASE (ERR_INRT_INDEXRANGE)
          message = 'index out of range.'
        CASE (ERR_INRT_BADSWITCH)
          message = 'invalid input value for an interaction switch.'
        CASE (ERR_INRT_BADACTV)
          message = 'invalid value for Activeness.'
        CASE (ERR_INRT_BADPERM)
          message = 'invalid value for Permission level.'
        CASE (ERR_INRT_BADVAL)
          message = 'invalid input value for an interaction.'
        CASE (ERR_INRT_BADMAT)
          message = 'invalid input name for a material.'
        CASE (ERR_INRT_MISSINGVAL)
          message = 'missing input data for an interaction.'
        CASE (ERR_INRT_ACTVPLAG)
          message = 'Lagrangian particle constituents differ in Activeness.'
        CASE (ERR_INRT_ACTVSMOKE)
          message = 'Smoke type cannot be more active than Gas.'
        CASE (ERR_INRT_ACTVDSMOKE)
          message = 'whether is as active as Gas was altered for interaction.'
        CASE (ERR_INRT_PARAMETER)
          message = 'inconsistent values assigned to parameters.'
        CASE (ERR_INRT_NINTL)
          message = 'number of Internal Nodes must be 0 or 1.'
        CASE (ERR_INRT_NINPUTEDGES)
          message = 'Internal Node needs least one input and one output Edge.'
        CASE (ERR_INRT_CONNECTINTL)
          message = 'Node is misconnected to Internal Node.'
        CASE (ERR_INRT_PERMINTL)
          message = 'invalid Permission token on Internal Node.'
        CASE (ERR_INRT_PERMLEVINTL)
          message = 'invalid Permission Level for Internal Node.'
        CASE (ERR_INRT_BURNING1)
          message = 'inconsistency in Activeness for Burning interaction.'
        CASE (ERR_INRT_ONLY1)
          message = 'only one oxidizer or boiling product smoke type allowed.'
        CASE (ERR_INRT_OX_ACTV)
          message = 'oxidizer smoke type cannot be active.'
        CASE (ERR_INRT_BOIL_ACTV)
          message = 'boiling product smoke type must be active.'
        CASE (ERR_INRT_BOIL_SAME)
          message = 'boiling input and output need same physical properties.'
        CASE (ERR_INRT_NPCLS)
          message = 'interactions with particles exist, but not particles.'
        CASE (ERR_INRT_ENERVAPOR)
          message = 'vapor energy positive, but nothing should be creating it.'
        CASE (ERR_INRT_NOINRT)
          message = 'interaction not implemented.'
#endif

! ------------------------------------------------------------------------------
!       ROCSPECIES specific errors
! ------------------------------------------------------------------------------

        CASE (ERR_SPEC_MODULE)
          message = 'Inconsistency of flow model and Rocspecies module.'
        CASE (ERR_SPEC_NTYPES)
          message = 'Number of SPECIES_TYPE sections unequal to NSPECIES'
        CASE (ERR_SPEC_MAXEQN) 
          message = 'NSPECIES inconsistent with MAXEQN:'
        CASE (ERR_SPEC_PROPS_INVALID) 
          message = 'MOLW and SPHT invalid: gamma out of bounds' 
        CASE (ERR_SPEC_NSPEC_INVALID)
          message = 'Number of species invalid:' 
        CASE (ERR_SPEC_SOURCE_TYPE_INVALID)
          message = 'Source type invalid.'     

! ------------------------------------------------------------------------------
!       If everything fails...
! ------------------------------------------------------------------------------

        CASE DEFAULT
          message = 'reason unknown.'
      END SELECT

! ------------------------------------------------------------------------------
!     Write error message
! ------------------------------------------------------------------------------

      IF (PRESENT(addMessage)) THEN
        message = TRIM(message)//' '//TRIM(addMessage)
      ENDIF ! PRESENT

      WRITE(STDERR,'(A)') SOLVER_NAME
      WRITE(STDERR,'(A,1X,A,I5.5,A,A)') SOLVER_NAME,'ERROR (proc. ', &
                                        global%myProcid,') - ',TRIM(message)

      WRITE(STDERR,'(A,1X,5(A),I4)') &
        SOLVER_NAME,'Function: ',TRIM(global%functionTree(1,global%nFunTree)), &
        ', file: ',TRIM(global%functionTree(2,global%nFunTree)), &
        ', line: ',errorLine

      DO i=global%nFunTree-1,1,-1    ! write out call tree
        WRITE(STDERR,'(A,1X,4(A))') &
          SOLVER_NAME,'Called from: ',TRIM(global%functionTree(1,i)), &
          ', file: ',TRIM(global%functionTree(2,i))
      ENDDO

      WRITE(STDERR,'(A)') SOLVER_NAME

! ------------------------------------------------------------------------------
!     Stop the run
! ------------------------------------------------------------------------------

#ifdef RFLU
      CALL MPI_Initialized(flag,errorFlag)

      IF ( flag .EQV. .TRUE. ) THEN
        IF ( global%nProcAlloc == 1 ) THEN
          CALL MPI_Finalize(errorFlag)
        ELSE
          CALL MPI_Abort(global%mpiComm,errorCode2,errorFlag)
        END IF ! global%nProcAlloc
      END IF ! flag

      STOP 1
#endif

#ifdef RFLO
#ifdef MPI
      IF (global%nProcAlloc == 1) THEN
        CALL MPI_Finalize( error )
      ELSE
        CALL MPI_Abort( error )
      ENDIF
#endif
      STOP 1
#endif
    END SUBROUTINE ErrorStop

END MODULE ModError

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ModError.F90,v $
! Revision 1.151  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.150  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.149  2007/03/27 00:18:07  haselbac
! Added new error condition for PLAG
!
! Revision 1.148  2006/10/20 21:29:54  mparmar
! Added ERR_BCVAR_VALUE_INVALID for NSCBC
!
! Revision 1.147  2006/08/21 16:46:01  haselbac
! Added error condition
!
! Revision 1.146  2006/04/07 15:19:18  haselbac
! Removed tabs
!
! Revision 1.145  2006/04/07 14:45:03  haselbac
! Added new error cases for patch flatness and orientation
!
! Revision 1.144  2006/03/30 20:48:39  haselbac
! Added error treatment of spec src terms
!
! Revision 1.143  2006/03/26 20:21:49  haselbac
! Removed ifdefs on SPEC, required for error checking
!
! Revision 1.142  2006/03/25 21:45:24  haselbac
! Added error conditions for sype patches
!
! Revision 1.141  2006/03/24 23:34:10  wasistho
! added ERR_DEPENDENT_INPUT
!
! Revision 1.140  2006/03/22 03:04:08  wasistho
! added ERR_PATCH_NUMBER
!
! Revision 1.139  2006/01/06 22:07:49  haselbac
! Added stencil error conditions, changed ERR_ASSOCIATED message
!
! Revision 1.138  2005/12/24 21:27:00  haselbac
! Added error treatment for ICT tolerance
!
! Revision 1.137  2005/12/01 18:37:54  fnajjar
! Added error trap for missing value
!
! Revision 1.136  2005/12/01 17:11:37  fnajjar
! Added error trap for randSeedType
!
! Revision 1.135  2005/11/14 16:56:35  haselbac
! Modified message for invalid gas model
!
! Revision 1.134  2005/11/10 02:20:54  haselbac
! Added error treatment for gas model and species
!
! Revision 1.133  2005/11/04 14:08:35  haselbac
! Added face normal error treatment
!
! Revision 1.132  2005/10/31 19:26:31  haselbac
! Added error treatment for gas model
!
! Revision 1.131  2005/10/27 18:57:18  haselbac
! Added err treatment of invalid constraints
!
! Revision 1.130  2005/10/14 14:03:31  haselbac
! Added ERR_TB_NEGATIVE
!
! Revision 1.129  2005/10/05 20:04:02  haselbac
! Added error treatment for ENSIGHT filter
!
! Revision 1.128  2005/09/20 15:47:53  fnajjar
! Added error definition for iPclSend memory overflow
!
! Revision 1.127  2005/08/05 15:27:04  haselbac
! Added new condition, cleaned up
!
! Revision 1.126  2005/08/03 18:54:12  hdewey2
! Added parameter for invalid solver type
!
! Revision 1.125  2005/07/14 21:40:26  haselbac
! Added new error conditions for invalid DISCR and ORDER
!
! Revision 1.124  2005/07/11 19:24:54  mparmar
! Aded error treatment for invalid reconst option
!
! Revision 1.123  2005/07/04 17:20:46  haselbac
! Bug fix: Proper error treatment depending on MPI
!
! Revision 1.122  2005/06/14 17:46:14  haselbac
! Added ERR_ALLOCATE_ADAPTIVE parameter and treatment
!
! Revision 1.121  2005/05/26 22:01:29  haselbac
! Fixed bug: MPI_Abort expects three arguments
!
! Revision 1.120  2005/05/16 20:41:35  haselbac
! Changed calling of MPI_Finalize and MPI_Abort
!
! Revision 1.119  2005/04/27 18:36:30  fnajjar
! Added trap error for findPclMethod
!
! Revision 1.118  2005/04/25 04:58:27  wasistho
! added ERR_FACE_INVERTED
!
! Revision 1.117  2005/04/20 02:49:25  wasistho
! added ERR_COMPILE_OPTION
!
! Revision 1.116  2005/04/15 15:06:27  haselbac
! Removed Charm/FEM error parameters, added MPI error parameters
!
! Revision 1.115  2005/03/09 23:16:20  gzheng
! when compiled under charm (CHARM=1), ErrorStop should also call MPI_Abort instead of calling STOP.
!
! Revision 1.114  2005/03/09 14:54:18  haselbac
! Added error treatment for virtual boundaries
!
! Revision 1.113  2005/01/17 19:55:56  haselbac
! Added error condition and treatment
!
! Revision 1.112  2005/01/14 21:12:49  haselbac
! Added error condition for MPI
!
! Revision 1.111  2004/12/19 15:44:15  haselbac
! Added PETSC error condition
!
! Revision 1.110  2004/12/04 03:22:57  haselbac
! Added error condition for estimate of number of vertices
!
! Revision 1.109  2004/11/30 20:10:50  fnajjar
! Added error definition for RK schemes
!
! Revision 1.108  2004/11/11 14:50:45  haselbac
! Removed CHARM section for writing error message, broken on popovich in serial
!
! Revision 1.107  2004/11/09 10:55:28  wasistho
! added error option due to inclusion statistics in rflopost
!
! Revision 1.106  2004/11/03 16:59:23  haselbac
! Removed error treatment related to HACK_PERIODIC
!
! Revision 1.105  2004/11/03 14:54:53  haselbac
! Added error conditions for GAMBIT grid conversion
!
! Revision 1.104  2004/10/19 19:28:40  haselbac
! Added new error conditions, cosmetics
!
! Revision 1.103  2004/09/29 00:52:29  wasistho
! added Radiation error-msg: flux limited diffusion input
!
! Revision 1.102  2004/09/27 01:36:13  haselbac
! Added error message for special faces
!
! Revision 1.101  2004/07/28 18:54:57  fnajjar
! Added overflow memory error trap for PLAG
!
! Revision 1.100  2004/07/23 22:43:43  wasistho
! added ERR_SYSTEM_COMMAND
!
! Revision 1.99  2004/06/17 15:18:14  fnajjar
! Included proper error trapping for ejection model
!
! Revision 1.98  2004/06/17 14:30:56  fnajjar
! Redefined error parameter from ERR_PLAG_INJCMODEL to ERR_PLAG_INCJDIAMDIST
!
! Revision 1.97  2004/06/16 20:00:49  haselbac
! Added Tecplot error condition
!
! Revision 1.96  2004/04/01 21:27:22  haselbac
! Added error condition ERR_SPEC_MAXEQN
!
! Revision 1.95  2004/03/05 23:21:27  haselbac
! Added two new PLAG error conditions
!
! Revision 1.94  2004/03/05 22:09:01  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.93  2004/03/02 21:47:28  jferry
! Added After Update interactions
!
! Revision 1.92  2004/02/26 21:01:56  haselbac
! Improved readability and added PLAG error conditions
!
! Revision 1.91  2004/01/29 22:57:20  haselbac
! Added three new error conditions
!
! Revision 1.90  2003/12/05 16:53:54  haselbac
! Added and changed error parameters
!
! Revision 1.89  2003/12/04 03:28:23  haselbac
! Added various error conditions
!
! Revision 1.88  2003/11/25 21:03:09  haselbac
! Added error support for rocspecies
!
! Revision 1.87  2003/11/21 22:38:57  fnajjar
! Added generic error messages
!
! Revision 1.86  2003/09/25 15:48:43  jferry
! implemented Boiling Regulation interaction
!
! Revision 1.85  2003/09/19 20:35:25  jferry
! Implemented oxidizer species for burning interaction
!
! Revision 1.84  2003/09/17 21:06:25  fnajjar
! Included error traps for injection model
!
! Revision 1.83  2003/09/13 20:16:50  fnajjar
! Added error traps for Breakup model
!
! Revision 1.82  2003/09/10 23:36:24  fnajjar
! Removed flags that are subsumed with Rocinteract
!
! Revision 1.81  2003/08/19 22:45:46  haselbac
! Added code for COBALT conversion errors
!
! Revision 1.80  2003/08/06 15:50:36  wasistho
! added turb. input error code
!
! Revision 1.79  2003/07/30 22:19:54  wasistho
! enter part and smoke data into radiation
!
! Revision 1.78  2003/07/17 01:00:18  wasistho
! initial activation rocrad
!
! Revision 1.77  2003/07/08 21:21:37  jblazek
! Modified start up procedure for dual-time stepping.
!
! Revision 1.76  2003/05/31 01:42:57  wasistho
! installed turb. wall layer model
!
! Revision 1.75  2003/05/24 02:13:51  wasistho
! turbulence statistics expanded
!
! Revision 1.74  2003/05/13 23:48:14  haselbac
! Added error treatment for negative flame temperature
!
! Revision 1.73  2003/05/01 20:44:55  haselbac
! Added ERR_MDOT_NEGATIVE and corresponding CASE
!
! Revision 1.72  2003/04/10 23:18:07  fnajjar
! Included error trap for unknown viscosity model
!
! Revision 1.71  2003/04/09 22:51:44  jferry
! removed peul_save and peul_verify structures
!
! Revision 1.70  2003/04/07 14:21:53  haselbac
! Added param and message for c2f list
!
! Revision 1.69  2003/04/03 21:10:17  jferry
! implemented additional safety checks for rocinteract
!
! Revision 1.68  2003/04/02 22:32:03  jferry
! codified Activeness and Permission structures for rocinteract
!
! Revision 1.67  2003/04/01 17:03:24  haselbac
! Added error condition for special cells
!
! Revision 1.66  2003/03/29 03:27:18  wasistho
! install ROCPERI
!
! Revision 1.65  2003/03/24 23:30:52  jferry
! overhauled rocinteract to allow interaction design to use user input
!
! Revision 1.64  2003/03/15 18:48:04  haselbac
! Added error for gm
!
! Revision 1.63  2003/03/15 17:44:10  haselbac
! Added several new error conditions
!
! Revision 1.62  2003/03/04 22:12:34  jferry
! Initial import of Rocinteract
!
! Revision 1.61  2003/02/25 21:11:11  fnajjar
! Added Error for PLAG Tile size
!
! Revision 1.60  2003/02/20 19:48:32  haselbac
! Added error conditions
!
! Revision 1.59  2003/02/12 20:49:51  jferry
! Moved Rocsmoke range to 4000-4999
!
! Revision 1.58  2003/02/11 22:52:50  jferry
! Initial import of Rocsmoke
!
! Revision 1.57  2003/02/06 19:30:24  haselbac
! Added ERR_MOVEPATCH_BC_INVALID
!
! Revision 1.56  2003/01/28 16:40:12  haselbac
! Added three new error conditions
!
! Revision 1.55  2002/10/27 19:01:01  haselbac
! Added several error conditions
!
! Revision 1.54  2002/10/25 14:03:51  f-najjar
! Define PLAG Error Message
!
! Revision 1.53  2002/10/17 14:12:10  haselbac
! Added error condition for number of coupled boundaries
!
! Revision 1.52  2002/10/12 19:11:20  haselbac
! Added new message and fixed bug
!
! Revision 1.51  2002/10/07 14:10:14  haselbac
! Removed tab
!
! Revision 1.50  2002/10/05 18:58:03  haselbac
! Added error condition for boundary vertex list
!
! Revision 1.49  2002/10/04 20:36:05  jblazek
! Extended range check of nFunTree.
!
! Revision 1.48  2002/09/25 18:29:57  jferry
! simplified TBC parameter lists
!
! Revision 1.47  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.46  2002/09/17 22:49:44  jferry
! Deleted tabs
!
! Revision 1.45  2002/09/17 13:43:00  jferry
! Added Time-dependent boundary conditions
!
! Revision 1.44  2002/09/13 14:54:09  haselbac
! Whoops. Deleted a few lines too many last time...
!
! Revision 1.43  2002/09/09 14:52:42  haselbac
! Added several error flags and proper output for parallel runs with FEM FW
!
! Revision 1.42  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.41  2002/08/30 01:47:58  jblazek
! Added support for moving grids.
!
! Revision 1.40  2002/08/23 03:16:29  wasistho
! modify ERR_TURB_MODULE
!
! Revision 1.39  2002/08/18 02:23:23  wasistho
! Added some error msg pertinent to TURB
!
! Revision 1.38  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
! Revision 1.37  2002/07/29 17:10:45  jblazek
! Put TURB stuff into #ifdef.
!
! Revision 1.36  2002/07/27 08:08:42  wasistho
! prepared for rocturb preliminary stage
!
! Revision 1.35  2002/07/25 15:13:11  haselbac
! Added various new error conditions for OLES, DCUHRE, and gradients
!
! Revision 1.34  2002/06/30 00:01:44  jblazek
! Removed TAB characters. Grrrrr ...
!
! Revision 1.33  2002/06/27 15:54:48  haselbac
! Added ERR_NCELLS_WRONG
!
! Revision 1.32  2002/06/22 01:13:37  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.31  2002/06/17 15:39:09  haselbac
! Prefixed SOLVER_NAME to all screen output
!
! Revision 1.30  2002/06/17 15:20:30  jblazek
! Added ERR_PREVIOUS_ERRORS flag.
!
! Revision 1.29  2002/06/14 21:34:32  wasistho
! Added time avg statistics
!
! Revision 1.28  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.27  2002/06/05 18:50:22  haselbac
! Added treatment for CHARM errors
!
! Revision 1.26  2002/05/28 13:52:00  haselbac
! Added FEM framework error message
!
! Revision 1.25  2002/05/21 01:48:22  wasistho
! add viscous terms
!
! Revision 1.24  2002/05/04 16:58:32  haselbac
! Added ERR_PROC_MISMATCH and file name for file errors
!
! Revision 1.23  2002/04/11 18:55:42  haselbac
! Added various new error codes
!
! Revision 1.22  2002/03/26 19:16:01  haselbac
! Added ROCFLU error conditions
!
! Revision 1.21  2002/03/21 18:07:15  jblazek
! Added check of MPI_PATCHOFF (for tags).
!
! Revision 1.20  2002/03/18 23:07:19  jblazek
! Finished multiblock and MPI.
!
! Revision 1.19  2002/03/01 16:42:53  haselbac
! Added some more error conditions
!
! Revision 1.18  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.17  2002/02/09 01:47:01  jblazek
! Added multi-probe option, residual smoothing, physical time step.
!
! Revision 1.16  2002/02/08 15:06:26  haselbac
! Added data structure errors
!
! Revision 1.15  2002/01/31 20:23:59  jblazek
! Added treatment of edge & corner cells.
!
! Revision 1.14  2002/01/12 00:02:48  jblazek
! Added postprocessor.
!
! Revision 1.13  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.12  2002/01/10 18:21:29  jblazek
! Added iteration number and initial residual to solution file.
!
! Revision 1.11  2002/01/10 00:02:07  jblazek
! Added calculation of mixture properties.
!
! Revision 1.10  2002/01/08 22:09:16  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.9  2002/01/02 16:20:19  jblazek
! Added flow initialization and dummy cell geometry.
!
! Revision 1.8  2001/12/22 00:09:38  jblazek
! Added routines to store grid and solution.
!
! Revision 1.7  2001/12/21 23:04:54  haselbac
! Added ROCFLU error parameters
!
! Revision 1.6  2001/12/19 23:09:21  jblazek
! Added routines to read grid and solution.
!
! Revision 1.5  2001/12/10 15:28:26  jblazek
! Fix to output format.
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






