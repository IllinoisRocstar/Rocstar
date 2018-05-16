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
! Purpose: set inflow boundary condition for one cell.
!
! Description: The subsonic boundary condition is based on the extrapolation of 
!   the Riemann invariant from the interior field (Holmes, D.G.: Inviscid 2D 
!   Solutions on Unstructured, Adaptive Grids. VKI Lecture Series 1989-06, 1989). 
!   The supersonic inflow boundary condition computes the conservative variables
!   from given velocity components, density and pressure.
!
! Input: bcOptType  = boundary treatment: subsonic, supersonic, or mixed
!        bcOptFixed = whether _computed_ inflow angle should be fixed or not
!        ptot       = given total pressure
!        ttot       = given total temperature
!        betah      = given inlet angle wrp. to y-axis
!        betav      = given inlet angle wrp. to z-axis
!        sx/y/zn    = components of normalized face vector (outward facing)
!        cpgas      = specific heat at constant pressure (boundary cell)
!        mm         = molecular mass at boundary cell
!        rl         = given density
!        ru/v/wl    = given velocity components
!
! Output: rr      = density at boundary
!         ru/v/wr = density * velocity components at boundary
!         rer     = density * total internal energy at boundary
!         pr      = pressure at boundary
!
! Notes: 
!   1. This condition is valid only for thermally and calorically perfect  
!      gas.
!   2. Important to avoid division by MakeNonZero(sl) when computing eta 
!      because that computation can become undefined for reservoir inflow
!      conditions, i.e., for a vanishing velocity vector. 
!
!******************************************************************************
!
! $Id: BcondInflowPerf.F90,v 1.3 2008/12/06 08:44:08 mtcampbe Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
!******************************************************************************

SUBROUTINE BcondInflowPerf(bcOptType,bcOptFixed,ptot,ttot,betah,betav,mach, & 
                           sxn,syn,szn,cpgas,mm,rl,rul,rvl,rwl,rr,rur,rvr, & 
                           rwr,rer,pr)

  USE ModDataTypes
  USE ModParameters
  USE ModTools, ONLY     : MakeNonZero
  USE ModInterfaces, ONLY: MixtPerf_C_Co2GUVW, MixtPerf_C_DGP, MixtPerf_C_GRT, &
        MixtPerf_Co2_CGUVW, MixtPerf_C2_GRT,MixtPerf_D_PRT, & 
        MixtPerf_Eo_DGPUVW, MixtPerf_Eo_DGPVm, MixtPerf_G_CpR, &
        MixtPerf_P_GMaPo, MixtPerf_P_GPoTTo, MixtPerf_Po_GPTTo, & 
        MixtPerf_Po_CGPUVW, MixtPerf_R_M, MixtPerf_T_CGR, MixtPerf_T_GMaTo, & 
        MixtPerf_Vm_C2Co2G

  IMPLICIT NONE

! ... parameters
  INTEGER, INTENT(IN) :: bcOptFixed,bcOptType

  REAL(RFREAL), INTENT(IN) :: betah, betav, cpgas, mach, mm, sxn, syn, szn, &
                              ptot, rl, rul, rvl, rwl, ttot
  REAL(RFREAL), INTENT(OUT) :: rer, rr, rur, rvr, rwr, pr

! ... local variables
  REAL(RFREAL) :: al, ar, a02, cp, disc, eta, g, gm1, igm1, ql, rgas, Rm, &
                  sl, sr , tr, ul, ur, vl, vr, wl, wr

!******************************************************************************
! gas properties

  rgas = MixtPerf_R_M(mm)
  g    = MixtPerf_G_CpR(cpgas,rgas)

! subsonic or mixed -----------------------------------------------------------   

  IF ( bcOptType == BCOPT_SUBSONIC .OR. bcOptType == BCOPT_MIXED ) THEN     
    gm1  = g - 1.0_RFREAL
    igm1 = 1.0_RFREAL/gm1

    ul = rul/rl
    vl = rvl/rl
    wl = rwl/rl

    a02 = MixtPerf_C2_GRT(g,rgas,ttot)

    al = MixtPerf_C_Co2GUVW(a02,g,ul,vl,wl) ! make al consistent with a02
    ql = ul*sxn + vl*syn + wl*szn

! - subsonic

    IF ( ABS(ql) < al ) THEN
      sl = SQRT(ul*ul + vl*vl + wl*wl)
      
      IF ( bcOptFixed == BCOPT_FIXED_NO ) THEN 
        IF ( sl > 1.0E-6_RFREAL ) THEN ! Avoid ill-defined angle computation
          eta = ql/sl        
        ELSE 
          eta = -1.0_RFREAL
        END IF ! sl
      ELSE 
        eta = -1.0_RFREAL
      END IF ! bcOptFixed

      Rm   = al - 0.5_RFREAL*gm1*ql
      disc = 0.5_RFREAL*gm1*eta**2* & 
             (a02/(Rm*Rm)*(1.0_RFREAL + 0.5_RFREAL*gm1*eta**2) - 1.0_RFREAL)     
        
      IF ( disc < 0.0_RFREAL ) THEN ! discriminant cannot be negative
        ar = SQRT(a02)
        tr = ttot
        pr = ptot
        sr = 0.0_RFREAL
      ELSE
        ar = Rm/(1.0_RFREAL + 0.5_RFREAL*gm1*eta*eta)*(1.0_RFREAL+SQRT(disc))                   
                      
        tr = MixtPerf_T_CGR(ar,g,rgas)
        pr = MixtPerf_P_GPoTTo(g,ptot,tr,ttot)
        sr = MixtPerf_Vm_C2Co2G(ar*ar,a02,g)
      END IF ! disc    
                  
      rr = MixtPerf_D_PRT( pr,rgas,tr )
      ur = sr*COS(betah)*COS(betav)
      vr = sr*SIN(betah)
      wr = sr*COS(betah)*SIN(betav)

      rer = rr*MixtPerf_Eo_DGPVm(rr,g,pr,sr)
      rur = rr*ur
      rvr = rr*vr
      rwr = rr*wr

! - supersonic

    ELSE
      IF ( bcOptType == BCOPT_MIXED ) THEN
        pr = mixtPerf_P_GMaPo(g,mach,ptot)
        tr = mixtPerf_T_GMaTo(g,mach,ttot)
        rr = mixtPerf_D_PRT(pr,rgas,tr)
        ar = mixtPerf_C_GRT(g,rgas,tr)

        ur = mach*ar*COS(betah)*COS(betav)
        vr = mach*ar*SIN(betah)
        wr = mach*ar*COS(betah)*SIN(betav)

        rer = rr*MixtPerf_Eo_DGPUVW(rr,g,pr,ur,vr,wr)
        rur = rr*ur
        rvr = rr*vr
        rwr = rr*wr
      END IF ! bcOptType
    END IF ! ql < al

! supersonic ------------------------------------------------------------------

  ELSE ! bcOptType == BCOPT_SUPERSONIC
    pr = mixtPerf_P_GMaPo(g,mach,ptot)
    tr = mixtPerf_T_GMaTo(g,mach,ttot)
    rr = mixtPerf_D_PRT(pr,rgas,tr)
    ar = mixtPerf_C_GRT(g,rgas,tr)

    ur = mach*ar*COS(betah)*COS(betav)
    vr = mach*ar*SIN(betah)
    wr = mach*ar*COS(betah)*SIN(betav) 

    rer = rr*MixtPerf_Eo_DGPUVW(rr,g,pr,ur,vr,wr)
    rur = rr*ur
    rvr = rr*vr
    rwr = rr*wr
  END IF ! bcOptType

END SUBROUTINE BcondInflowPerf

!******************************************************************************
!
! RCS Revision history:
!
! $Log: BcondInflowPerf.F90,v $
! Revision 1.3  2008/12/06 08:44:08  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:17:22  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2004/12/01 16:47:56  haselbac
! Initial revision after changing case
!
! Revision 1.7  2004/01/29 22:52:36  haselbac
! Added bcOptFixed, fixed bug, clean-up
!
! Revision 1.6  2003/12/04 03:22:56  haselbac
! Fixed bug in formulation, added partial fix for eta
!
! Revision 1.5  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.2  2002/06/22 00:49:50  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.1  2002/06/10 21:19:34  haselbac
! Initial revision
!
!******************************************************************************






