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
! Purpose: define the overall data structure related to regions
!          and to grid levels.
!
! Description: none
!
! Notes: none
!
!******************************************************************************
!
! $Id: ModDataStruct.F90,v 1.61 2008/12/06 08:44:18 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModDataStruct

  USE ModDataTypes
  USE ModParameters
  USE ModRandom, ONLY     : t_rand_data
  USE ModBndPatch, ONLY   : t_patch
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY       : t_grid
  USE ModInteract, ONLY   : t_inrt_input
  USE ModMixture, ONLY    : t_mixt, t_mixt_input

#ifdef RFLO
#ifdef RADI
  USE ModRadiation, ONLY  : t_radi, t_radi_input
#endif
#ifdef PLAG
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModPartLag, ONLY    : t_buffer_plag
#endif
#ifdef SPEC
!  USE ModSpecies, ONLY    : t_spec, t_spec_input
#endif
#ifdef PEUL
  USE ModPartEul, ONLY    : t_peul, t_peul_input
#endif
#ifdef TURB
  USE ModTurbulence, ONLY : t_turb, t_turb_input
#endif
#ifdef PERI
  USE ModPeriodic, ONLY   : t_peri, t_peri_input
#endif
#endif

#ifdef RFLU
  USE ModRadiation, ONLY: t_radi,t_radi_input
  USE ModPartLag, ONLY: t_plag,t_plag_input
  USE ModPartLag, ONLY: t_buffer_plag
  USE ModPlotting, ONLY: t_plot
  USE ModSpecies, ONLY: t_spec,t_spec_input
  USE ModTurbulence, ONLY: t_turb,t_turb_input
  USE ModPeriodic, ONLY: t_peri,t_peri_input
#endif

#ifdef PETSC

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscsnes.h"
#include "include/finclude/petscis.h"
#include "include/finclude/petscao.h"

#endif

  IMPLICIT NONE

#ifdef RFLO

! edge & corner cells

  TYPE t_dCellSrc
    LOGICAL :: rotate
    INTEGER :: srcRegion, srcCell
    INTEGER :: srcIndexMapMat(3,4)
#ifdef PLAG
    TYPE(t_buffer_plag) :: bufferExchPlag
#endif
  END TYPE t_dCellSrc

  TYPE t_dCell
    LOGICAL :: interact
    INTEGER :: degenrt, interType
    TYPE(t_dCellSrc), POINTER :: cells(:)
  END TYPE t_dCell

  TYPE t_dCellTransf
    INTEGER :: nCells, iRequest
    REAL(RFREAL), POINTER :: buff(:)
    INTEGER :: iRequestMetrics
    REAL(RFREAL), POINTER :: buffMetrics(:)
#ifdef PLAG
    INTEGER :: nBuffSizePlag, iRequestPlag
    INTEGER, POINTER :: buffPlagI(:)
    REAL(RFREAL), POINTER :: buffPlagR(:)
#endif
  END TYPE t_dCellTransf

! grid level related data

  TYPE t_level
    REAL(RFREAL), POINTER :: dt(:)

    TYPE(t_grid) :: grid, gridOld, gridOld2
    TYPE(t_mixt) :: mixt
#ifdef TURB
    TYPE(t_turb) :: turb
#endif
#ifdef SPEC
!    TYPE(t_spec) :: spec
#endif
#ifdef RADI
    TYPE(t_radi) :: radi
#endif
#ifdef PEUL
    TYPE(t_peul) :: peul
#endif
#ifdef PLAG
    TYPE(t_plag) :: plag, plagTemp
#endif
#ifdef PERI
    TYPE(t_peri) :: peri
#endif

    TYPE(t_patch), POINTER :: patches(:)

    TYPE(t_dCell)                :: edgeCells(12), cornerCells(8)
    TYPE(t_dCellTransf), POINTER :: sendEcCells(:), recvEcCells(:)
#ifdef RADI
    TYPE(t_dCellTransf), POINTER :: sndRadiEcCells(:), rcvRadiEcCells(:)
#endif
#ifdef TURB
    TYPE(t_dCellTransf), POINTER :: sndTurbEcCells(:), rcvTurbEcCells(:)
#endif
#ifdef PEUL
    TYPE(t_dCellTransf), POINTER :: sndPeulEcCells(:), rcvPeulEcCells(:)
#endif
  END TYPE t_level

! region related data

  TYPE t_region
    INTEGER :: localNumber, active, procid, iRegionGlobal
    INTEGER :: nDumCells, nEdgeCells, nPatches, nGridLevels
    INTEGER :: startLevel, currLevel, irkStep
    INTEGER :: blockShape
    INTEGER :: dimWork1D, dimWork2D(2)

    REAL(RFREAL), POINTER :: work1D(:), work2D(:,:)

    TYPE(t_rand_data) :: randData

    TYPE(t_mixt_input)          :: mixtInput
#ifdef TURB
    TYPE(t_turb_input)          :: turbInput
#endif
#ifdef SPEC
!    TYPE(t_spec_input)          :: specInput
#endif
#ifdef PEUL
    TYPE(t_peul_input)          :: peulInput
#ifndef PLAG
    TYPE(t_inrt_input), POINTER :: inrtInput
#endif
#endif
#ifdef PLAG
    TYPE(t_plag_input)          :: plagInput
    TYPE(t_inrt_input), POINTER :: inrtInput
#endif
#ifdef RADI
    TYPE(t_radi_input)          :: radiInput
#endif
#ifdef PERI
    TYPE(t_peri_input)          :: periInput
#endif
    TYPE(t_level),      POINTER :: levels(:)
    TYPE(t_global),     POINTER :: global

  END TYPE t_region
#endif

! *****************************************************************************

#ifdef RFLU

! region-related data

  TYPE t_region
    LOGICAL :: postActiveFlag
    LOGICAL, DIMENSION(:), POINTER :: thrustFlagsGlobal
    INTEGER :: dtMinLoc,iRegionGlobal,irkStep
    INTEGER :: fieldFlagCoord,fieldFlagGmDisp,fieldFlagGmRhs,fieldFlagMixt, &
               fieldFlagSpec
    INTEGER :: dimWork1D, dimWork2D(2)

    REAL(RFREAL):: dtMin
    REAL(RFREAL), DIMENSION(:), POINTER :: dt
    REAL(RFREAL), DIMENSION(:,:), POINTER :: massCoeffsGlobal, &
                                             specImpulseGlobal, &
                                             specImpulseVacGlobal, &
                                             thrustGlobal, &
                                             thrustVacGlobal
    REAL(RFREAL), DIMENSION(:,:,:), POINTER :: forceCoeffsGlobal, &
                                               forceVacCoeffsGlobal, &
                                               momentCoeffsGlobal

    TYPE(t_rand_data) :: randData

    TYPE(t_grid) :: grid,gridOld,gridOld2
    TYPE(t_mixt) :: mixt
    TYPE(t_turb) :: turb
    TYPE(t_spec) :: spec
    TYPE(t_radi) :: radi
!    TYPE(t_peul) :: peul
    TYPE(t_peri) :: peri
    TYPE(t_plag) :: plag,plagTemp
    TYPE(t_plot) :: plot

    TYPE(t_patch), DIMENSION(:), POINTER :: patches
    TYPE(t_global), POINTER :: global
    TYPE(t_region), POINTER :: pRegion

    TYPE(t_mixt_input) :: mixtInput
    TYPE(t_turb_input) :: turbInput
    TYPE(t_spec_input) :: specInput
!    TYPE(t_peul_input) :: peulInput
    TYPE(t_plag_input) :: plagInput
    TYPE(t_radi_input) :: radiInput
    TYPE(t_peri_input) :: periInput
    TYPE(t_inrt_input), POINTER :: inrtInput

#ifndef NO_TECPLOT
    REAL(RFREAL), DIMENSION(:,:), POINTER :: varCellTEC,varVertTEC
#endif

! TEMPORARY
!    PetscFortranAddr :: poissonInfoPETSc(RFLU_PETSC_POISSON_INFO_BEG: & 
!                                RFLU_PETSC_POISSON_INFO_END)
! END TEMPORARY

#ifdef PETSC
    Vec :: x, r
    Mat :: A, preA
    SNES :: snes
    MatFDColoring :: fdcolor
    AO :: ao
#endif

  END TYPE t_region

! level-related data

  TYPE t_level
    TYPE(t_region), DIMENSION(:), POINTER :: regions
  END TYPE t_level

  TYPE(t_level), DIMENSION(:), ALLOCATABLE :: levels
#endif

END MODULE ModDataStruct

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModDataStruct.F90,v $
! Revision 1.61  2008/12/06 08:44:18  mtcampbe
! Updated license.
!
! Revision 1.60  2008/11/19 22:17:29  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.59  2007/03/19 21:46:49  haselbac
! Added t_plot, cosmetics
!
! Revision 1.58  2006/10/21 04:15:21  mparmar
! Rearranged variables
!
! Revision 1.57  2006/10/20 21:29:04  mparmar
! Added Global thrustFlags, forceVac, massCoeffs, specImpulse, specImpulseVac, thrust, thrustVac
!
! Revision 1.56  2006/03/04 04:32:52  wasistho
! added blockShape in rocflo data structure
!
! Revision 1.55  2006/02/15 21:06:39  wasistho
! defined inrt_input in peul if not plag
!
! Revision 1.54  2006/02/13 21:50:48  wasistho
! added ifdefs for phys modules in RFLO
!
! Revision 1.53  2006/02/08 21:01:58  hdewey2
! Added gridOld2
!
! Revision 1.52  2005/08/17 20:15:26  hdewey2
! Added PETSc Application Ordering context
!
! Revision 1.51  2005/08/02 18:16:07  hdewey2
! Added PETSc variables, uncommented PETSC include statements
!
! Revision 1.50  2005/06/29 22:49:33  wasistho
! added interType in dCell type
!
! Revision 1.49  2005/01/13 21:47:19  haselbac
! Comment out changes for now - lead to compilation error
!
! Revision 1.48  2005/01/13 21:41:36  haselbac
! Bug fix: Different declaration for PETSc info variables
!
! Revision 1.47  2004/12/19 15:43:44  haselbac
! Added poissonInfoPETSc
!
! Revision 1.46  2004/09/27 22:29:37  wasistho
! added sndRadiEcCells and rcvRadiEcCells
!
! Revision 1.45  2004/08/21 00:28:11  wasistho
! added degenrt in edge/corner data structure
!
! Revision 1.44  2004/07/28 18:54:00  fnajjar
! Added plagTemp for dynamic memory reallocation
!
! Revision 1.43  2004/07/23 22:43:16  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.42  2004/06/16 20:00:48  haselbac
! Added force coefficients and Tecplot variables
!
! Revision 1.41  2004/03/10 23:08:34  fnajjar
! Added infrastructure for MPI corner-edge cells
!
! Revision 1.40  2004/03/02 21:44:05  jferry
! Added corner and edge cell data structures and routines
!
! Revision 1.39  2004/02/10 21:21:16  fnajjar
! Added index mapping matrix for proper geometry transform
!
! Revision 1.38  2004/01/31 03:57:21  haselbac
! Changed instance of t_plag in RFLU, clean-up
!
! Revision 1.37  2004/01/23 00:30:54  wasistho
! added turb. edge/corners data structures under type t_level
!
! Revision 1.36  2004/01/15 21:08:42  fnajjar
! Included datastructure for corner-edge cells metrics
!
! Revision 1.35  2003/11/25 21:03:08  haselbac
! Added fieldFlagSpec for rocspecies support
!
! Revision 1.34  2003/11/12 21:17:49  fnajjar
! Included Corner-Edge Cells Infrastructure for PLAG
!
! Revision 1.33  2003/10/15 02:40:49  haselbac
! Added dtMinLoc and dtMin as members of region
!
! Revision 1.32  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.31  2003/08/07 15:31:40  haselbac
! Changed var name
!
! Revision 1.30  2003/07/22 02:00:05  haselbac
! Added pRegion
!
! Revision 1.29  2003/06/04 22:05:17  haselbac
! Added activeFlag
!
! Revision 1.28  2003/03/31 16:13:17  haselbac
! Added fieldFlagGmDisp
!
! Revision 1.27  2003/03/29 03:27:10  wasistho
! install ROCPERI
!
! Revision 1.26  2003/03/15 17:41:40  haselbac
! Added field flags for gm communication
!
! Revision 1.25  2003/03/04 22:12:34  jferry
! Initial import of Rocinteract
!
! Revision 1.24  2003/02/18 14:57:31  jferry
! Implemented portable random number generator ModRandom
!
! Revision 1.23  2003/02/14 22:32:36  jblazek
! Finished implementation of corener and edge cells.
!
! Revision 1.22  2003/02/05 21:07:30  jblazek
! Coordinated stop of a run works now for MPI.
!
! Revision 1.21  2003/02/03 19:20:47  jblazek
! Added treatment of edge and corner cells for one processor.
!
! Revision 1.20  2002/10/25 14:03:17  f-najjar
! Redefine plag without POINTER description as single-class Lagrangian
! particles are considered
!
! Revision 1.19  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.18  2002/09/17 22:51:23  jferry
! Removed Fast Eulerian particle type
!
! Revision 1.17  2002/09/09 14:51:42  haselbac
! Input types now under regions
!
! Revision 1.16  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.15  2002/08/30 19:08:58  jblazek
! Dimensions of work arrays now set in derivedInputValues.
!
! Revision 1.14  2002/06/27 15:53:25  haselbac
! Added fieldFlagMixt for communication
!
! Revision 1.13  2002/05/04 16:56:21  haselbac
! Added dt and irkStep for RFLU
!
! Revision 1.12  2002/03/01 16:39:45  haselbac
! Deleted old regions, changed domains to regions
!
! Revision 1.11  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.10  2002/02/08 15:05:29  haselbac
! Added iDomainGlobal variable
!
! Revision 1.9  2002/02/01 00:00:24  jblazek
! Edge and corner cells defined for each level.
!
! Revision 1.8  2002/01/31 20:23:59  jblazek
! Added treatment of edge & corner cells.
!
! Revision 1.7  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.6  2002/01/23 03:51:24  jblazek
! Added low-level time-stepping routines.
!
! Revision 1.5  2002/01/11 17:20:19  jblazek
! Added time stamp or iteration number to file names.
!
! Revision 1.4  2002/01/08 22:09:16  jblazek
! Added calculation of face vectors and volumes.
!
! Revision 1.3  2002/01/02 16:20:19  jblazek
! Added flow initialization and dummy cell geometry.
!
! Revision 1.2  2001/12/21 23:01:38  haselbac
! Added ROCFLU data types
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
!******************************************************************************






