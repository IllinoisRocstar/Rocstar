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
! Purpose: define various parameters pertinent to RADI
!
! Description: none
!
! Notes: none
!
!******************************************************************************
!
! $Id: RADI_ModParameters.F90,v 1.7 2008/12/06 08:44:37 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

MODULE RADI_ModParameters

  USE ModDataTypes
  IMPLICIT NONE

! RADI integer parameters -----------------------------------------------------

! simplified radiation models:  

  INTEGER, PARAMETER :: RADI_MODEL_NONE     = 0, &   ! non active
                        RADI_MODEL_ROSS     = 1, &   ! Rosseland diff appr.
                        RADI_MODEL_FLDSRC   = 2, &   ! FLD by just E-src term
                        RADI_MODEL_FLDTRAN  = 3, &   ! FLD with solving Er eq.
                        RADI_MODEL_RTEGRAY  = 4, &   ! RTE gray
                        RADI_MODEL_RTEBAND  = 5      ! RTE non-gray

! radiation cv, dv, pressure-tensor, flux-limiter, grads, etc:

  INTEGER, PARAMETER :: CV_RADI_ENER        = 1      ! radiation energy

  INTEGER, PARAMETER :: DV_RADI_TEFF        = 1      ! effective temperature

  INTEGER, PARAMETER :: E11 = 1, &                 ! symm tensor components
                        E12 = 2, & 
                        E13 = 3, &
                        E22 = 4, &
                        E23 = 5, &
                        E33 = 6, &
                        TENSOR_SYMM_NELM = 6       ! elm number of symm tensor

  INTEGER, PARAMETER :: GR_RADI_EX        = 1, &   ! gradients of
                        GR_RADI_EY        = 2, &   ! radiation energy
                        GR_RADI_EZ        = 3      

  INTEGER, PARAMETER :: FLD_LIM_NONE      = 0, &   ! no limiter (pure diff =1/3)
                        FLD_LIM_LP        = 1      ! number of cv elements

! FLD numeric:

  INTEGER, PARAMETER :: FLD_DISCR_CEN  = 0, &  ! central discretization
                        FLD_DISCR_UPW  = 1     ! upwind

  INTEGER, PARAMETER :: FLD_DISCR_ORD1 = 1, &  ! discretization order
                        FLD_DISCR_ORD2 = 2         

! radiation bc model:

  INTEGER, PARAMETER :: RADI_BC_NONE        = 0,  &   ! no specific bc
                        RADI_BC_RANGE       = 9,  &   ! range interval
                        RADI_BC_DIFFUS      = 10, &   ! diffuse boundary
                        RADI_BC_REFRAC      = 20, &   ! specular refraction
                        RADI_BC_CYCLIC      = 30      ! cyclic bc

! radiation source media to determine extinction coeffs.:

  INTEGER, PARAMETER :: RADI_MEDIA_ARTIF    = 1, &   ! artificial ext. coef
                        RADI_MEDIA_REAL     = 2      ! real media

! media phase of optical constants for extinct. coeffs. calculation:

  INTEGER, PARAMETER :: RADI_PHASE_GAS      = 1, &   ! pure air phase
                        RADI_PHASE_DISPART  = 2, &   ! Al2O3 phase
                        RADI_PHASE_CONPART  = 3, &   ! Al particles (no smoke)
                        NPHASE              = 3      ! number of media phase

! phase properties of optical constants for extinct. coeffs. calculation:

  INTEGER, PARAMETER :: PHASE_PROP_V        = 1, &   ! cell/field avg.vol.frac.
                        PHASE_PROP_D        = 2, &   ! cell/field avg.diameter
                        PHASE_PROP_Q        = 3, &   ! ext. efficiency
                        NPROPERTY           = 3      ! number of properties

! numerical method:

  INTEGER, PARAMETER :: RADI_NUM_NONE       = 0, &   ! no solver required 
                        RADI_NUM_DOM4       = 1, &   ! discrete ordinate S4
                        RADI_NUM_DOM8       = 2, &   ! S8
                        RADI_NUM_DOM16      = 3, &   ! S16
                        RADI_NUM_FVM        = 4      ! finite volume

! radiation transfer coefficients:

  INTEGER, PARAMETER :: RADI_COEFF_EXTINCT  = 1, &   ! extinction coeff.
                        RADI_COEFF_SCATTER  = 2, &   ! scattering coeff.
                        RADI_COEFF_PLANCK   = 3, &   ! scattering coeff.
                        RADI_COEFF_NCOMP    = 3      ! no of components

! angular direction:

  INTEGER, PARAMETER :: RADI_ANGLE_POLAR    = 1, &   ! polar direction
                        RADI_ANGLE_AZIMU    = 2, &   ! azimuthal direction
                        RADI_ANGLE_NCOMP    = 2      ! no of components

! PERI real parameters --------------------------------------------------------

  REAL(RFREAL), PARAMETER :: RADI_REAL_SMALL  = 1.E-16_RFREAL ! small real nmbr
  REAL(RFREAL), PARAMETER :: RADI_REAL_STOP   = 1.E-8_RFREAL  ! stop criterium 
  REAL(RFREAL), PARAMETER :: RADI_REAL_MICRON = 1.E-6_RFREAL  ! micron/meter 
  REAL(RFREAL), PARAMETER :: RADI_REAL_ECMIN  = 100._RFREAL   ! trsh xtinc.coef 

END MODULE RADI_ModParameters

!******************************************************************************
!
! RCS Revision history:
!
! $Log: RADI_ModParameters.F90,v $
! Revision 1.7  2008/12/06 08:44:37  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:50  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2004/09/30 17:10:12  wasistho
! prepared for full FLD radiation model
!
! Revision 1.4  2004/09/22 01:30:24  wasistho
! switch LFD to FLD for flux limited diffusion
!
! Revision 1.3  2004/09/18 17:40:36  wasistho
! install Limited Flux Diffusion radiation
!
! Revision 1.2  2003/07/30 22:22:20  wasistho
! enter part and smoke data into radiation
!
! Revision 1.1  2003/07/17 01:16:59  wasistho
! initial activation rocrad
!
!
!
!******************************************************************************






