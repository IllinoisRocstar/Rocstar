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
! Purpose: define data types related to radiation.
!
! Description: Input data are grouped in t_radi_input data type, while
!              computed radiation field quantities in t_radi data type.
!
! Notes: none
!
!******************************************************************************
!
! $Id: ModRadiation.F90,v 1.13 2008/12/06 08:44:19 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE ModRadiation

  USE ModDataTypes
  IMPLICIT NONE

! input -----------------------------------------------------------------------

  TYPE t_radi_input

! - general input

    INTEGER               :: radiModel      ! radiation model
    INTEGER               :: media          ! type of participating media
    INTEGER               :: nDv            ! number of derived vars. comp. 
    INTEGER               :: nAng           ! # intensity angles to be output
    REAL(RFREAL)          :: stBoltz        ! Stefan-Boltzman constant
    REAL(RFREAL), POINTER :: angles(:,:)    ! selected angle of rad.intensity
    REAL(RFREAL), POINTER :: optConst(:,:)  ! optical constants
    CHARACTER(256)        :: line(2)        ! angle readin line-paramters

! - diffusion model input

    INTEGER      :: fluxLim       ! type of flux limiter
    INTEGER      :: nCv, nGrad    ! number of cv and gradient components
    INTEGER      :: spaceDiscr    ! spatial discretization type
    INTEGER      :: spaceOrder    ! spatial discretization order
    REAL(RFREAL) :: vis2          ! 2nd order dissipation coefficient
    REAL(RFREAL) :: vis4          ! 4th order dissipation coefficient
    REAL(RFREAL) :: smoocf        ! residual smoothing coefficient    

! - RTE solver input

    INTEGER :: solMethod     ! numerical method to solve RTE
    INTEGER :: nOrdin        ! number of ordinate directions 
    INTEGER :: nPol          ! number of polar angles
    INTEGER :: nAzi          ! number of azimuthal angles
  END TYPE t_radi_input

! data ------------------------------------------------------------------------

  TYPE t_radi
    REAL(RFREAL), POINTER :: dv(:,:)                            ! all
    REAL(RFREAL), POINTER :: cv(:,:), cvOld(:,:)                ! FLDTRAN+RTE
    REAL(RFREAL), POINTER :: rhs(:,:), rhsSum(:,:), diss(:,:)   ! FLDTRAN+RTE
    REAL(RFREAL), POINTER :: gradi(:,:), gradj(:,:), gradk(:,:) ! FLDTRAN+RTE(?)
#ifdef RFLO
    REAL(RFREAL), POINTER :: srad(:,:), epsIrs(:,:)             ! FLDTRAN+RTE
#endif
    REAL(RFREAL), POINTER :: fluxlim(:)                         ! FLD+ROSS
    REAL(RFREAL), POINTER :: ptens(:,:), eddFact(:)             ! FLD
    REAL(RFREAL), POINTER :: qri(:), qrj(:), qrk(:), goFact(:)     ! all
    REAL(RFREAL), POINTER :: radInt(:,:), radCoef(:,:)             ! all
    REAL(RFREAL), POINTER :: wvInt(:,:)                            ! allwv 
    REAL(RFREAL), POINTER :: dWghti(:,:), dWghtj(:,:), dWghtk(:,:) ! RTE

  END TYPE t_radi

END MODULE ModRadiation

!******************************************************************************
!
! RCS Revision history:
!
! $Log: ModRadiation.F90,v $
! Revision 1.13  2008/12/06 08:44:19  mtcampbe
! Updated license.
!
! Revision 1.12  2008/11/19 22:17:30  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.11  2004/09/30 17:06:47  wasistho
! prepared for full FLD radiation model
!
! Revision 1.10  2004/09/22 01:30:06  wasistho
! switch LFD to FLD for flux limited diffusion
!
! Revision 1.9  2004/09/18 17:40:06  wasistho
! install Limited Flux Diffusion radiation
!
! Revision 1.8  2004/09/14 00:51:09  wasistho
! prepared for LimitedFluxDiffusion (LFD)
!
! Revision 1.7  2003/08/09 02:20:33  wasistho
! restore version 1.5
!
! Revision 1.6  2003/08/08 02:57:54  wasistho
! added Genx buffers
!
! Revision 1.5  2003/07/30 22:21:30  wasistho
! enter part and smoke data into radiation
!
! Revision 1.4  2003/07/18 01:40:11  wasistho
! removed bcModel from input data structure
!
! Revision 1.3  2003/07/17 01:01:13  wasistho
! initial activation rocrad
!
! Revision 1.2  2001/12/11 21:59:29  jblazek
! memory allocation added.
!
! Revision 1.1.1.1  2001/12/03 21:44:05  jblazek
! Import of RocfluidMP
!
!******************************************************************************






