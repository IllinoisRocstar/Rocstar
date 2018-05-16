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
! Purpose: define data types related to turbulence model.
!
! Description: -distinction is made between LES, RaNS/DES specific, or
!               common LES and RaNS/DES input and data
!              -in RaNS/DES equations are solved for turbulence quantities
!              -in LES subgrid viscous terms are computed according to 
!               selected SGS model
!
! Notes:  leswv: LES work variables
!         altwv: altered wv: allocated/used differently in LES, RaNS regions 
!         comwv: common wv for LES and RaNS
!         altv : altered permanent variables: allocated/used differently
!         comv : common permanent variables for LES and RaNS
!
!******************************************************************************
!
! $Id: ModTurbulence.F90,v 1.29 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModTurbulence

  USE ModDataTypes
  IMPLICIT NONE

! input -----------------------------------------------------------------------

  TYPE t_turb_input

! - general input

    INTEGER             :: modelClass  ! class of model, LES or RANS
    INTEGER             :: nOutField   ! number of solution fields output
    INTEGER             :: nZof        ! number zero-one switch fields

! - RaNS/DES input

    INTEGER             :: nTurbEqs    ! number of turb.equations
    INTEGER             :: wDistMethod ! wall dist.: 0 direct, 1 hierarchical
    INTEGER             :: functV1     ! visc.function, SA: 0 pow3(std), 1 pow2
    INTEGER             :: spaceDiscr  ! spatial discretization type
    INTEGER             :: spaceOrder  ! spatial discretization order
    REAL(RFREAL)        :: vis2        ! 2nd order dissipation coefficient
    REAL(RFREAL)        :: vis4        ! 4th order dissipation coefficient
    REAL(RFREAL)        :: smoocf      ! residual smoothing coefficient    
    REAL(RFREAL)        :: cDes        ! DES lengthscale coefficient    
    REAL(RFREAL), POINTER :: const(:)  ! RaNS model constants

! - LES input

    REAL(RFREAL) :: cSmag             ! Smagorinsky constant, typically 0.1-0.2
    REAL(RFREAL) :: xyzSmag(3)        ! tripping location of fixed Smagorinsky 
    REAL(RFREAL) :: wallRough         ! max. wall roughness size in curr. region 

    REAL(RFREAL) :: delFac2           ! filter width scaling factor
    INTEGER      :: filterType        ! filter type, uniform or non-uniform
    INTEGER      :: deltaType         ! model length scale formula, sqrt or cbrt
    INTEGER      :: filterWidth(3)    ! base filter width, 0,1, or 2 grd-spacing
    INTEGER      :: homDir(3)         ! identifier if i,j, or k dir. homogeneous

    INTEGER      :: engModel          ! switch for subgrid energy models
    INTEGER      :: calcVort          ! compute vortic.0=no 1=per-fdt 2=per-sdt
    INTEGER      :: nCv,nDv,nSv,nGrad ! number of cv, dv, sv and gradients
    INTEGER      :: nSt,nFixSt        ! number of st and permanent turb.nstats
    INTEGER      :: wallModel         ! most general WLM in current region 
    INTEGER      :: wlmRefPoint       ! WLM reference point above the wall

! - POST input

    INTEGER      :: nPostv            ! number of postprocessed variables
  END TYPE t_turb_input

! data ------------------------------------------------------------------------

  TYPE t_turb

! - RaNS/DES data

    REAL(RFREAL), POINTER :: cv(:,:), cvOld(:,:)
    REAL(RFREAL), POINTER :: rhs(:,:), rhsSum(:,:), diss(:,:), dsterm(:,:)
    REAL(RFREAL), POINTER :: lens(:)
#ifdef RFLO
    REAL(RFREAL), POINTER :: srad(:,:), epsIrs(:,:)
#endif
    REAL(RFREAL), POINTER :: tv(:,:)                             ! ranswv

! - additional Rans/DES data with dual-timestepping
    REAL(RFREAL), POINTER :: cvn(:,:), cvn1(:,:), cvn2(:,:), sDual(:,:)

! - LES data

#ifdef RFLO
    REAL(RFREAL), POINTER :: ccCofi1(:,:) ,ccCofi2(:,:) ,ccCofi4(:,:)
    REAL(RFREAL), POINTER :: ccCofj1(:,:) ,ccCofj2(:,:) ,ccCofj4(:,:)
    REAL(RFREAL), POINTER :: ccCofk1(:,:) ,ccCofk2(:,:) ,ccCofk4(:,:)
    REAL(RFREAL), POINTER :: ffCofi1I(:,:),ffCofi2I(:,:),ffCofi4I(:,:)
    REAL(RFREAL), POINTER :: ffCofi1J(:,:),ffCofi2J(:,:),ffCofi4J(:,:)
    REAL(RFREAL), POINTER :: ffCofi1K(:,:),ffCofi2K(:,:),ffCofi4K(:,:)
    REAL(RFREAL), POINTER :: ffCofj1I(:,:),ffCofj2I(:,:),ffCofj4I(:,:)
    REAL(RFREAL), POINTER :: ffCofj1J(:,:),ffCofj2J(:,:),ffCofj4J(:,:)
    REAL(RFREAL), POINTER :: ffCofj1K(:,:),ffCofj2K(:,:),ffCofj4K(:,:)
    REAL(RFREAL), POINTER :: ffCofk1I(:,:),ffCofk2I(:,:),ffCofk4I(:,:)
    REAL(RFREAL), POINTER :: ffCofk1J(:,:),ffCofk2J(:,:),ffCofk4J(:,:)
    REAL(RFREAL), POINTER :: ffCofk1K(:,:),ffCofk2K(:,:),ffCofk4K(:,:)
    REAL(RFREAL), POINTER :: fvolI(:)  ,fvolJ(:)  ,fvolK(:)
    REAL(RFREAL), POINTER :: fISij(:,:),fJSij(:,:),fKSij(:,:)    !leswv
#endif
#ifdef RFLU
    REAL(RFREAL), POINTER :: avgCoI(:,:), bAvgCoI(:,:)
    REAL(RFREAL), POINTER :: ccCofi1(:,:) ,ccCofi2(:,:) ,ccCofi4(:,:)
    REAL(RFREAL), POINTER :: ffCofi1I(:,:),ffCofi2I(:,:),ffCofi4I(:,:)
    REAL(RFREAL), POINTER :: bffCofi1I(:,:),bffCofi2I(:,:),bffCofi4I(:,:)
    REAL(RFREAL), POINTER :: fvolI(:)  ,bfVolI(:)
    REAL(RFREAL), POINTER :: fISij(:,:),bfISij(:,:)              !leswv
#endif
    REAL(RFREAL), POINTER :: lij(:,:),mij(:,:)                   !leswv
    REAL(RFREAL), POINTER :: fVar(:,:),ffVar(:,:),ccVar(:,:)     !leswv
    REAL(RFREAL), POINTER :: trace(:),coef(:,:),mueT(:,:)        !leswv
#ifdef RFLU
    REAL(RFREAL), POINTER :: bLij(:,:),bMij(:,:)                 !leswv
    REAL(RFREAL), POINTER :: bfVar(:,:),bffVar(:,:)              !leswv
    REAL(RFREAL), POINTER :: bCoef(:,:),bMueT(:,:)               !leswv
#endif

! - common LES and RaNS/DES data

#ifdef RFLO
    REAL(RFREAL), POINTER :: gradi(:,:),gradj(:,:),gradk(:,:)    !altwv
    REAL(RFREAL), POINTER :: mISij(:,:),mJSij(:,:),mKSij(:,:)    !comwv
    REAL(RFREAL), POINTER :: workI(:,:),workJ(:,:),workK(:,:)    !comwv
#endif
#ifdef RFLU
    REAL(RFREAL), POINTER :: gradi(:,:,:),bGradi(:,:,:)          !altwv
    REAL(RFREAL), POINTER :: mISij(:,:),bmISij(:,:)              !comwv
#endif
    REAL(RFREAL), POINTER :: dv(:,:)                             !altv
    REAL(RFREAL), POINTER :: sv(:,:),vort(:,:)                   !comv
#ifdef RFLO
    REAL(RFREAL), POINTER :: zofi(:,:,:),zofj(:,:,:),zofk(:,:,:) !comwv
#endif
#ifdef RFLU
    REAL(RFREAL), POINTER :: zofi(:,:,:),bZofi(:,:,:)            !comwv
#endif

! - statistics data

#ifdef STATS
    REAL(RFREAL), POINTER :: tav(:,:),st(:,:)                    !comv
    REAL(RFREAL), POINTER :: stwork(:,:)                         !wv
#endif

! - postprocessing data

    REAL(RFREAL), POINTER :: postv(:,:), postvVert(:,:)          !comv
#ifdef STATS
    REAL(RFREAL), POINTER :: tavVert(:,:)                        !comv
#endif

! - in summary, RFLU boundary data: 
!   bAvgCoI, bfVolI, bGradi                  (allocMemory)
!   bMuet, bfISij, bmISij                    (coViscFluxesFlu)
!   bLij, bMij, bfVar, bffVar, bCoef         (lesCalcEddyVis)

  END TYPE t_turb

END MODULE ModTurbulence

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModTurbulence.F90,v $
! Revision 1.29  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.28  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.27  2006/01/12 09:42:19  wasistho
! added post proc variables
!
! Revision 1.26  2004/12/09 22:16:06  wasistho
! added data structures
!
! Revision 1.25  2004/10/22 23:14:25  wasistho
! added work variable for statistics
!
! Revision 1.24  2004/08/04 02:48:51  wasistho
! removed turb%avgCoI,J,K as it is defined as grid%c2fCoI,J,K
!
! Revision 1.23  2004/05/28 02:02:23  wasistho
! update unstructured grid LES
!
! Revision 1.22  2004/03/25 04:40:03  wasistho
! added boundary data needed for RFLU
!
! Revision 1.21  2004/03/19 02:41:12  wasistho
! prepared for RFLU
!
! Revision 1.20  2004/03/13 03:06:18  wasistho
! prepared for RFLU
!
! Revision 1.19  2004/02/26 21:15:24  wasistho
! delete esg1Sum
!
! Revision 1.18  2004/02/19 04:02:42  wasistho
! added new rans/SA parameter VISCFUNCTION
!
! Revision 1.17  2004/02/14 03:42:06  wasistho
! added new WLM parameter: reference point
!
! Revision 1.16  2004/02/11 03:23:31  wasistho
! added feature: variable number of turbulence output fields
!
! Revision 1.15  2003/10/26 00:07:44  wasistho
! added multiple discr.types and order
!
! Revision 1.14  2003/10/21 20:32:55  wasistho
! added dt relaxation in steady flow due to RANS source term
!
! Revision 1.13  2003/10/15 03:40:33  wasistho
! added 2nd order dissipation coeff. k2
!
! Revision 1.12  2003/10/09 20:49:20  wasistho
! added DES lengthscale coefficient CDES
!
! Revision 1.11  2003/10/03 20:14:29  wasistho
! initial installation of turbModel SA and DES
!
! Revision 1.10  2003/08/09 02:21:13  wasistho
! restore version 1.8
!
! Revision 1.9  2003/08/08 02:57:31  wasistho
! added Genx buffers
!
! Revision 1.8  2003/08/06 15:51:15  wasistho
! added vorticities computation
!
! Revision 1.7  2003/08/01 22:12:56  wasistho
! changed soluFile to calcVort
!
! Revision 1.6  2003/07/22 02:55:07  wasistho
! prepare accurate rocturb restart
!
! Revision 1.5  2003/05/31 01:40:21  wasistho
! installed wall layer model
!
! Revision 1.4  2003/05/24 02:06:34  wasistho
! turbulence statistics expanded
!
! Revision 1.3  2002/11/02 01:54:45  wasistho
! Added TURB statistics
!
! Revision 1.2  2002/10/14 23:51:12  wasistho
! Install Rocturb
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
!******************************************************************************






