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
! Purpose: set far field boundary condition for one cell.
!
! Description: the subsonic boundary condition is based on non-reflecting,
!              characteristics method of Whitfield and Janus: Three-Dimensional
!              Unsteady Euler Equations Solution Using Flux Vector Splitting.
!              AIAA Paper 84-1552, 1984. The supersonic boundary condition
!              consists of simple extrapolation.
!
! Input: machInf  = given Mach number at "infinity"
!        alphaInf = angle of attack
!        betaInf  = slip angle
!        pInf     = given static pressure
!        tInf     = given static temperature
!        rho      = density at boundary cell
!        rhou/v/w = density * velocity components at boundary cell
!        rhoe     = density * total energy at boundary cell
!        press    = static pressure at boundary cell
!        sx/y/zn  = components of ortho-normalized face vector (outward facing)
!        csound   = speed of sound at boundary cell
!        cpgas    = specific heat at constant pressure (boundary cell)
!        mol      = molecular mass at boundary cell
!
! Output: rhob      = density at boundary
!         rhou/v/wb = density * velocity components at boundary
!         rhoeb     = density * total energy at boundary
!         pb        = pressure at boundary
!
! Notes: this condition is valid only for a thermally and calorically
!        perfect gas (supersonic in/outflow valid for all gases).
!
!******************************************************************************
!
! $Id: BcondFarfieldPerf.F90,v 1.3 2008/12/06 08:44:08 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE BcondFarfieldPerf( machInf,alphaInf,betaInf,pInf,tInf, &
                              sxn,syn,szn,cpgas,mol, &
                              rho,rhou,rhov,rhow,rhoe,press, &
                              rhob,rhoub,rhovb,rhowb,rhoeb,pb )

  USE ModDataTypes
  USE ModParameters
  USE ModInterfaces, ONLY: MixtPerf_C_DGP, MixtPerf_C_GRT, MixtPerf_D_PRT, &
        MixtPerf_Eo_DGPUVW, MixtPerf_Eo_DGPVm, MixtPerf_G_CpR, MixtPerf_R_M, &
        MixtPerf_P_DEoGVm2

  IMPLICIT NONE

! ... parameters
  REAL(RFREAL) :: machInf, alphaInf, betaInf, pInf, tInf
  REAL(RFREAL) :: rho, rhou, rhov, rhow, rhoe, press
  REAL(RFREAL) :: sxn, syn, szn, cpgas, mol
  REAL(RFREAL) :: rhob, rhoub, rhovb, rhowb, rhoeb, pb

! ... local variables
  REAL(RFREAL) :: rgas, gamma, rhoInf, qInf, uInf, vInf, wInf
  REAL(RFREAL) :: re, ue, ve, we, pe, qn, crho0
  REAL(RFREAL) :: ra, ua, va, wa, pa, sgn, csound

!******************************************************************************
! gas properties

  rgas  = MixtPerf_R_M( mol )
  gamma = MixtPerf_G_CpR( cpgas,rgas )

! flow values at "infinity"

  rhoInf = MixtPerf_D_PRT( pInf,rgas,tInf )
  qInf   = machInf * MixtPerf_C_GRT( gamma,rgas,tInf )
  uInf   = qInf * COS(alphaInf) * COS(betaInf)
  vInf   = qInf * SIN(alphaInf) * COS(betaInf)
  wInf   = qInf *                 SIN(betaInf)

! flow values at a reference location (= interior cell)

  re  = rho
  ue  = rhou/rho
  ve  = rhov/rho
  we  = rhow/rho
  pe  = press
  qn  = sxn*ue + syn*ve + szn*we

! subsonic flow (qn<0: inflow / qn>0: outflow)

  IF (machInf < 1._RFREAL) THEN
    csound = MixtPerf_C_DGP( re,gamma,pe )
    crho0  = csound*re

    IF (qn < 0._RFREAL) THEN
      ra  = rhoInf
      ua  = uInf
      va  = vInf
      wa  = wInf
      pa  = pInf
      sgn = -1._RFREAL
      pb  = 0.5_RFREAL*(pa+pe-crho0*(sxn*(ua-ue)+syn*(va-ve)+szn*(wa-we)))
    ELSE
      ra  = re
      ua  = ue
      va  = ve
      wa  = we
      pa  = pe
      sgn = +1._RFREAL
      pb  = pInf
    ENDIF
    rhob  = ra + (pb-pa)/(csound*csound)
    rhoub = rhob*(ua+sgn*sxn*(pa-pb)/crho0)
    rhovb = rhob*(va+sgn*syn*(pa-pb)/crho0)
    rhowb = rhob*(wa+sgn*szn*(pa-pb)/crho0)
    rhoeb = rhob*MixtPerf_Eo_DGPUVW( rhob,gamma,pb,rhoub/rhob,rhovb/rhob, & 
                                     rhowb/rhob )

! supersonic flow (qn<0: inflow / qn>0: outflow)

  ELSE
    IF (qn < 0._RFREAL) THEN
      rhob  = rhoInf
      rhoub = rhoInf*uInf
      rhovb = rhoInf*vInf
      rhowb = rhoInf*wInf
      rhoeb = rhob*MixtPerf_Eo_DGPVm( rhoInf,gamma,pInf,qInf )
      pb    = pInf
    ELSE
      rhob  = rho
      rhoub = rhou
      rhovb = rhov
      rhowb = rhow
      rhoeb = rhoe
      pb    = press
    ENDIF
  ENDIF

END SUBROUTINE BcondFarfieldPerf

!******************************************************************************
!
! RCS Revision history:
!
! $Log: BcondFarfieldPerf.F90,v $
! Revision 1.3  2008/12/06 08:44:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:47:53  haselbac
! Initial revision after changing case
!
! Revision 1.6  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.3  2002/06/29 23:05:19  jblazek
! Removed gam1 because not used.
!
! Revision 1.2  2002/06/22 00:49:50  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.1  2002/06/10 21:19:34  haselbac
! Initial revision
!
!******************************************************************************






